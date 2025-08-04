#!/usr/bin/env nextflow

// Enable DSL2 syntax (introduces modularity and improved data flow manipulation)
nextflow.enable.dsl=2

/*
 * Split the big Lotta et al 2021 file by metabolite
 */
 process split_Lotta {
    publishDir "${projectDir}/data/Lotta_et_al_2021", mode: "copy"

    input:
        path Lotta_txt_gz

    output:
        path "*.txt.gz"

    script:
        """
        python3 ${projectDir}/bin/split_Lotta.py \
            --filename="${Lotta_txt_gz}"
        """
 }

 /*
 * The mediation analysis: take harmonized data, find intersection of variants, filter by p-value and Steiger, clump by instrument strength, and perform the analysis
 */
process prepare_mediation_data {
    input:
        tuple val(exposure_ID), val(exposure_type), path(in_exposure_gz), val(outcome_ID), val(outcome_type), path(in_outcome_gz), val(mediator_IDs), val(mediator_type), path(in_mediators_gz), val(instrument_pval_threshold)

    output:
        tuple val(exposure_ID), val(exposure_type), val(outcome_ID), val(outcome_type), val(mediator_type), val(instrument_pval_threshold), path("mediation_data.RDS")

    script:
        mediator_dir = "mediators"

        """
        mkdir ${mediator_dir}
        mv ${in_mediators_gz} ${mediator_dir}/

        Rscript ${projectDir}/bin/prepare_mediation_data.R \
            --in_exposure_gz=${in_exposure_gz} \
            --in_outcome_gz=${in_outcome_gz} \
            --mediator_dir=${mediator_dir} \
            --instrument_pval_threshold=${instrument_pval_threshold} \
            --out_mediation_data_RDS="mediation_data.RDS"
        """
}

/*
 * Concatenate MR-IVW results by exposure
 */
 process concat_mr_ivw {
    input:
        tuple val(exposure_ID), val(outcome_ID), path("mr_ivw*.txt")

    output:
        tuple val(exposure_ID), path("mr_ivw_${exposure_ID}.txt")

    script:
        """
        awk '{if (FNR!=1 || NR==1) print \$0}' mr_ivw*.txt > "mr_ivw_${exposure_ID}.txt"
        """
 }

include { mr_ivw; inner_join } from "./mr_ivw.nf"

 workflow multi_file_mediator_pipeline {
    take:
        UKBB_pairs    // exposure_ID, exposure_type, exposure_filename, outcome_ID, outcome_type, outcome_filename
        MR_UKBB_pairs // exposure_ID, exposure_type, exposure_filename (harmonized), outcome_ID, outcome_type, outcome_filename (harmonized)
        reference_dir
        instrument_selection_tools_R
        instrument_pval_threshold

    main:
        if (params.mediators == "Lotta_et_al_2021") {
            if (params.split_Lotta) {
                Lotta_txt_gz = Channel.fromPath("data/raw/metabolites/Lotta_et_al_2021/all.grch37.tabix.gz").first()

                split_Lotta(Lotta_txt_gz) // Split by metabolite
                mediators = split_Lotta.out.flatten()
            }
            else {
                mediators = Channel.fromPath("data/Lotta_et_al_2021/*.txt.gz")
            }
            mediators = mediators.map{file -> [file.simpleName, params.mediators, file]} // mediator_ID, mediator_type, mediator_filename
        }

        // Construct exposure-mediator pairs for MR (note: out of all the UKBB exposures, not just MR-significant)
        exposure_mediator_pairs = UKBB_pairs
            .map{[it.exposure_ID, it.exposure_type, it.exposure_filename]}     // exposure_ID, exposure_type, exposure_filename
            .unique()                                                          // select all exposures just once
            .combine(mediators)
            .map{[exposure_ID:it[0], exposure_type:it[1], exposure_filename:it[2], outcome_ID:it[3], outcome_type:it[4], outcome_filename:it[5]]}

        // Construct mediator-outcome pairs for MR (note: out of all the UKBB outcomes, not just MR-significant)
        mediator_outcome_pairs = UKBB_pairs
            .map{[it.outcome_ID, it.outcome_type, it.outcome_filename]}        // outcome_ID, outcome_type, outcome_filename
            .unique()                                                          // select all outcomes just once
            .combine(mediators)
            .map{[exposure_ID:it[3], exposure_type:it[4], exposure_filename:it[5], outcome_ID:it[0], outcome_type:it[1], outcome_filename:it[2]]}

        // Perform MR IVW between all exposure-mediator pairs (note: harmonization will be done again for the UKBB traits -- this is inefficient! -- but the outcome will be the same)
        mr_ivw(exposure_mediator_pairs.concat(mediator_outcome_pairs), reference_dir, instrument_selection_tools_R, instrument_pval_threshold)

        MR_UKBB_mediator_pairs = mr_ivw.out.mr_ivw_results.filter{it[3] == params.mediators}
        if (params.mediator_pval_threshold_from_exposure == "Bonferroni") {
            MR_UKBB_mediator_pairs = MR_UKBB_mediator_pairs.groupTuple(by: [0, 2, 3, 4])    // exposure_ID, [mediator IDs], exposure_type, mediator_type, exposure_filename, [mediator filenames], [exposure-mediator MR P-values]
                .map{[*it, it[1].size()]}                                                   // N_tests = number of mediators for each exposure
                .transpose()                                                                // remove the grouping
                .filter{Double.valueOf(it[6]) <= 0.05 / it[7]}
        }
        else {
            MR_UKBB_mediator_pairs = MR_UKBB_mediator_pairs.filter{Double.valueOf(it[6]) <= Double.valueOf(params.mediator_pval_threshold_from_exposure)}
        }

        MR_mediator_UKBB_pairs = mr_ivw.out.mr_ivw_results.filter{it[2] == params.mediators}
        if (params.mediator_pval_threshold_to_outcome == "Bonferroni") {
            MR_mediator_UKBB_pairs = MR_mediator_UKBB_pairs.groupTuple(by: [1, 2, 3, 5])    // [mediator IDs], outcome_ID, mediator_type, outcome_type, [mediator filenames], outcome_filename, [mediator-outcome MR P-values]
                .map{[*it, it[0].size()]}                                                   // N_tests = number of mediators for each outcome
                .transpose()                                                                // remove the grouping
                .filter{Double.valueOf(it[6]) <= 0.05 / it[7]}
        }
        else {
            MR_mediator_UKBB_pairs = MR_mediator_UKBB_pairs.filter{Double.valueOf(it[6]) <= Double.valueOf(params.mediator_pval_threshold_to_outcome)}
        }
        
        exposure_groups = MR_UKBB_mediator_pairs.map{[it[0], it[2], it[1], it[3], it[5]]}.groupTuple(by: [0, 1])    // exposure_ID, exposure_type, [mediator IDs], [mediator types], [mediator filenames]
        exposure_groups = inner_join(MR_UKBB_pairs, exposure_groups, 2)                                             // exposure_ID, exposure_type, exposure_filename, outcome_ID, outcome_type, outcome_filename, [mediator IDs], [mediator types], [mediatior filenames]

        outcome_groups = MR_mediator_UKBB_pairs.map{[it[1], it[3], it[0], it[2], it[4]]}.groupTuple(by: [0, 1])     // outcome_ID, outcome_type, [mediator IDs], [mediator types], [mediator filenames]
        outcome_groups = inner_join(MR_UKBB_pairs.map{[*it[3..5], *it[0..2]]}, outcome_groups, 2)
            .map{[*it[3..5], *it[0..2], *it[6..8]]}                                                                 // exposure_ID, exposure_type, exposure_filename, outcome_ID, outcome_type, outcome_filename, [mediator IDs], mediator_type, [mediatior filenames]

        // Prepare mediation data (I/O can take considerable time and we don't want to repeat it)
        mediation_triplets = exposure_groups.concat(outcome_groups)
            .transpose()
            .groupTuple(by: [0..8])                             // Group by everything
            .filter{it[0].size() == 2}                          // Both exposure->mediator and mediator->outcome associations are required
            .transpose()
            .unique()
            .groupTuple(by: [0, 1, 2, 3, 4, 5, 7])              // exposure_ID, exposure_type, exposure_filename, outcome_ID, outcome_type, outcome_filename, [mediator IDs], mediator_type, [mediatior filenames]
            .combine(instrument_pval_threshold)
        prepare_mediation_data(mediation_triplets)

        // Add MR-IVW results by exposure for consistency with the single_file_mediator_pipeline
        concat_mr_ivw(mr_ivw.out.mr_ivw_files.groupTuple())
        mediation_data = inner_join(prepare_mediation_data.out, concat_mr_ivw.out, 1)
            .combine(Channel.value(true))   // Indicate "clump = true" to clump mediator data later down the line

    emit:
        mediation_data = mediation_data
 }
