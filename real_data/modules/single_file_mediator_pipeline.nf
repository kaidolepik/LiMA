#!/usr/bin/env nextflow

// Enable DSL2 syntax (introduces modularity and improved data flow manipulation)
nextflow.enable.dsl=2

/*
 * Find common SNPs from groups of individual harmonized gz-files (e.g. UKBB trait files)
 */
process common_trait_SNPs {
    input:
        path files

    output:
        path "common_trait_SNPs.txt"
    
    script:
        """
        Rscript ${projectDir}/scripts/common_SNPs.R \
            --in_files_dir=. \
            --type=harmonized_trait \
            --out_SNPs="common_trait_SNPs.txt"
        """
}

/*
 * Find common mediator and trait SNPs for the analyses
 */
process common_analysis_SNPs {
    input:
        tuple val(prefix), path(besd), path(epi), path(esi)
        path common_trait_SNPs

    output:
        path "common_analysis_SNPs.txt"

    script:
        """
        Rscript ${projectDir}/scripts/common_SNPs.R \
            --in_files_dir=. \
            --type=${prefix} \
            --out_SNPs="common_analysis_SNPs.txt"
        """
}

/*
 * Find significant eQTL instruments for all mediators (proteins) and save to a single file
 */
process eQTL_analysis {
    input:
        tuple val(prefix), path(besd), path(epi), path(esi)
        path common_analysis_SNPs
        val instrument_pval_threshold
        val cis

    output:
        tuple val(prefix), path("significant_eQTLs.txt"), emit: eQTLs
        path "instrumented_mediators.txt", emit: eMediators_file


    script:
        """
        awk '{print \$2}' ${common_analysis_SNPs} > common_${prefix}_SNPs.txt

        if [[ ${cis} = true ]]
        then
            smr --beqtl-summary ${prefix} --extract-snp common_${prefix}_SNPs.txt --query ${instrument_pval_threshold} --cis-wind 2000 --out significant_eQTLs
        else
            smr --beqtl-summary ${prefix} --extract-snp common_${prefix}_SNPs.txt --query ${instrument_pval_threshold} --out significant_eQTLs
        fi

        awk 'NR>1 {print \$7}' "significant_eQTLs.txt" | sort | uniq > instrumented_mediators.txt
        """
}

/*
 * Harmonize and clump significant mediator eQTLs (we do it now for computational efficiency)
 */
process clump_mediators {
    input:
        each mediator_ID
        tuple val(type), path(significant_eQTLs)
        path reference_dir
        path instrument_selection_tools_R
        val instrument_pval_threshold

    output:
        tuple val(mediator_ID), val(type), path("${mediator_ID}_clumped.txt")

    script:
        """
        Rscript ${projectDir}/scripts/harmonize_data.R \
            --trait_ID=${mediator_ID} \
            --trait_type=${type} \
            --in_trait_gz=${significant_eQTLs} \
            --reference_dir=${reference_dir} \
            --instrument_selection_tools_R=${instrument_selection_tools_R} \
            --out_trait_gz=${mediator_ID}_harmonized.txt.gz

        Rscript ${projectDir}/scripts/clump_harmonized_data.R \
            --in_harmonized_gz=${mediator_ID}_harmonized.txt.gz \
            --reference_dir=${reference_dir} \
            --instrument_selection_tools_R=${instrument_selection_tools_R} \
            --instrument_pval_threshold=${instrument_pval_threshold} \
            --out_clumped_txt=${mediator_ID}_clumped.txt
        """
}

/*
 * Clump exposure instruments (we do it now for computational efficiency)
 */
process clump_exposures {
    input:
        tuple val(exposure_ID), val(exposure_type), path(harmonized_gz)
        path common_analysis_SNPs
        path reference_dir
        path instrument_selection_tools_R
        val instrument_pval_threshold

    output:
        tuple val(exposure_ID), val(exposure_type), path("${exposure_ID}_clumped.txt")

    script:
        """
        gzip -cd ${harmonized_gz} | \
            awk 'FNR==NR { ID_trait[\$3]; next } (FNR==1 || \$3 in ID_trait)' ${common_analysis_SNPs} - | \
            gzip \
        > filtered_${harmonized_gz}

        Rscript ${projectDir}/scripts/clump_harmonized_data.R \
            --in_harmonized_gz=filtered_${harmonized_gz} \
            --reference_dir=${reference_dir} \
            --instrument_selection_tools_R=${instrument_selection_tools_R} \
            --instrument_pval_threshold=${instrument_pval_threshold} \
            --out_clumped_txt=${exposure_ID}_clumped.txt
        """
}

/*
 * Concatenate clumped mediator exposure/instrument data
 */
process concat_mediator_exposure_data {
    publishDir "${projectDir}/data/${params.mediators}", mode: "copy", name: "mediators_as_exposure.txt.gz"

    input:
        path mediator_clump

    output:
        path "mediators_as_exposure.txt.gz"

    script:
        """
        awk '{ if (NR == 1 || FNR != 1) print \$0 }' *.txt | \
            gzip --stdout > mediators_as_exposure.txt.gz
        """
}

/*
 * Get mediator data for trait clumps (traits could be both exposures and other mediators)
 */
process mediator_output_data {
    input:
        tuple val(trait_ID), val(trait_type), path(trait_clump)
        tuple val(prefix), path(besd), path(epi), path(esi)
        path eMediators_file
        path common_analysis_SNPs
        path reference_dir
        path instrument_selection_tools_R

    output:
        tuple val(trait_ID), path("mediator_data_of_${trait_ID}_clumps.txt.gz")

    script:
        """
        awk 'FNR==NR { ID_trait[\$3]; next } { if (\$3 in ID_trait) print \$2 }' ${trait_clump} ${common_analysis_SNPs} > SNPs_esi.txt

        smr --beqtl-summary ${prefix} --extract-snp SNPs_esi.txt --extract-probe ${eMediators_file} --query 1 --out mediator_data_of_${trait_ID}_clumps
        
        Rscript ${projectDir}/scripts/harmonize_data.R \
            --trait_ID=${prefix} \
            --trait_type=${prefix} \
            --in_trait_gz=mediator_data_of_${trait_ID}_clumps.txt \
            --reference_dir=${reference_dir} \
            --instrument_selection_tools_R=${instrument_selection_tools_R} \
            --out_trait_gz=mediator_data_of_${trait_ID}_clumps.txt.gz
        """
}

/*
 * Concatenate harmonized mediator output data of instrument clumps
 * Saving without compressing because this file is read into many following R processes which uncompress it in tmp/;
 * for many concurrent processes, there is not enough storage to do so.
 */
process concat_mediator_output_data {
    publishDir "${projectDir}/data/${params.mediators}", mode: "copy", name: "mediators_as_output.txt"

    input:
        path harmonized_clumps_txt_gz

    output:
        path "mediators_as_output.txt"

    script:
        """
        cat *.txt.gz | \
            gzip -cd | \
            sort -u | \
            sort -n -k1,1 -k2,2 \
            > mediators_as_output.txt
        """
}

/*
 * Prepare mediation data
 */
process prepare_mediation_data_single {
    input:
        tuple val(exposure_ID), val(exposure_type), path(exposure_gz), val(outcome_ID), val(outcome_type), path(outcome_gz),
            path(mediator_exposure_data), path(mediator_output_data), val(mediator_type), path(mr_ivw_results_exposure), path(mr_ivw_results_outcome)
        val instrument_pval_threshold
        val mediator_pval_threshold_from_exposure
        val mediator_pval_threshold_to_outcome

    output:
        tuple val(exposure_ID), val(exposure_type), val(outcome_ID), val(outcome_type), val(mediator_type), val(instrument_pval_threshold), 
            path("mediation_data.RDS"), path(mr_ivw_results_exposure), val(false) // Indicate "clump = false" so as to not clump again later down the line

    script:
        """
        Rscript ${projectDir}/scripts/prepare_mediation_data_single.R \
            --exposure_ID=${exposure_ID} \
            --outcome_ID=${outcome_ID} \
            --in_exposure_gz=${exposure_gz} \
            --in_outcome_gz=${outcome_gz} \
            --in_mediators_as_exposure=${mediator_exposure_data} \
            --in_mediators_as_outcome=${mediator_output_data} \
            --in_mr_ivw_results_exposure=${mr_ivw_results_exposure} \
            --in_mr_ivw_results_outcome=${mr_ivw_results_outcome} \
            --mediator_pval_threshold_from_exposure=${mediator_pval_threshold_from_exposure} \
            --mediator_pval_threshold_to_outcome=${mediator_pval_threshold_to_outcome} \
            --out_mediation_data_RDS="mediation_data.RDS"
        """
}

def files_exist(files) {
    all_files_exist = true

    for (file in files) {
        all_files_exist = all_files_exist && (new File("${projectDir}/" + file)).exists()
    }

    return(all_files_exist)
}

include { mr_ivw_analysis; inner_join } from "./mr_ivw.nf"

 workflow single_file_mediator_pipeline {
    take:
        MR_UKBB_pairs // exposure ID, exposure type, exposure filename (harmonized), outcome ID, outcome type, outcome filename (harmonized)
        reference_dir
        instrument_selection_tools_R
        instrument_pval_threshold

    main:
        // Find common SNPs of harmonized traits
        trait_harmonized_files = MR_UKBB_pairs.map{[it[2]]}
            .mix(MR_UKBB_pairs.map{[it[5]]})
            .unique()
        common_trait_SNPs(trait_harmonized_files.collect())

        // Read mediator data into a channel
        if (params.mediators == "INTERVAL_plasma_proteins") {
            mediators = Channel.fromFilePairs("data/raw/proteins/INTERVAL_plasma/INTERVAL_plasma_proteins.{besd,epi,esi}", size: 3, flat: true)
        }

        // Define waypoint files
        file_mediators_exposure = "data/${params.mediators}/mediators_as_exposure.txt.gz"
        file_mediators_output = "data/${params.mediators}/mediators_as_output.txt"
        file_mr_ivw_exposure = "results/mr_ivw/UKBB_quantitative_${params.mediators}_5E-8.txt"
        file_mr_ivw_outcome = "results/mr_ivw/${params.mediators}_5E-8.txt" // This would have to be created manually from ${params.mediators}_UKBB_binary_5E-8.txt" and ${params.mediators}_UKBB_quantitative_5E-8.txt"

        // If using waypoints is not desired or waypoint files are not available then create the waypoint files
        if (!params.use_waypoints || !files_exist([file_mediators_exposure, file_mediators_output, file_mr_ivw_exposure, file_mr_ivw_outcome])) {
            // Find common mediator and trait SNPs
            common_analysis_SNPs(mediators, common_trait_SNPs.out)

            // Find all significant eQTLs and corresponding eMediators
            eQTL_analysis(mediators, common_analysis_SNPs.out, instrument_pval_threshold, params.cis)
            eMediator_IDs = eQTL_analysis.out.eMediators_file
                .splitCsv(header: false, strip: true)
                .collect()

            // Clump the harmonized trait data
            exposures = MR_UKBB_pairs.map{[it[0], it[1], it[2]]}.unique()
            clump_exposures(exposures, common_analysis_SNPs.out.first(), reference_dir, instrument_selection_tools_R, instrument_pval_threshold)

            // Clump the mediator exposure/instrument data (harmonize in the process) and concatenate
            clump_mediators(eMediator_IDs, eQTL_analysis.out.eQTLs, reference_dir, instrument_selection_tools_R, instrument_pval_threshold)
            concat_mediator_exposure_data(clump_mediators.out.map{it[2]}.collect())

            // Get mediator output data for clumped instruments and concatenate
            clumped_instruments = clump_mediators.out.concat(clump_exposures.out)
            mediator_output_data(clumped_instruments, mediators.first(), eQTL_analysis.out.eMediators_file.first(), common_analysis_SNPs.out.first(), reference_dir, instrument_selection_tools_R)
            concat_mediator_output_data(mediator_output_data.out.map{it[1]}.collect())

            // Construct exposure-mediator pairs for MR
            exposure_mediator_pairs = clump_exposures.out
                .join(mediator_output_data.out)
                .map{[it[0], params.mediators, it[1], params.mediators, it[2], it[3]]} // // exposure_ID, outcome_ID (mediator type here is a placeholder), exposure_type, outcome_type, exposure_filename, outcome_filename

            // Construct mediator-outcome pairs for MR
            mediator_outcome_pairs = MR_UKBB_pairs
                .map{[it[3], it[4], it[5]]}        // outcome_ID, outcome_type, outcome_filename
                .unique()                          // select all outcomes just once
                .combine(concat_mediator_exposure_data.out)
                .map{[params.mediators, it[0], params.mediators, it[1], it[3], it[2]]} // exposure_ID (mediator type here is a placeholder), outcome_ID, exposure_type, outcome_type, exposure_filename, outcome_filename

            // Perform MR IVW analysis between exposures-mediators (multi-mediator file for each exposure), mediators-outcomes (multi-mediator file for each outcome), and save the results
            mr_ivw_data = exposure_mediator_pairs
                .concat(mediator_outcome_pairs)
                .combine(instrument_pval_threshold)

            mr_ivw_analysis(mr_ivw_data, reference_dir, instrument_selection_tools_R)
            mr_ivw_analysis.out.mr_ivw_results_txt.map{it[2]}
                .collectFile(storeDir: "results/mr_ivw", keepHeader: true) { file -> ["${file.baseName.tokenize(':').get(0)}_${file.baseName.tokenize(':').get(1)}_${file.baseName.tokenize(':').get(2)}.txt", file] }

            // Create exposure-mediators-outcome triplets and prepare mediation data
            mediation_triplets = MR_UKBB_pairs
                .combine(concat_mediator_exposure_data.out)
                .combine(concat_mediator_output_data.out)

            mediation_triplets = inner_join(mediation_triplets, mr_ivw_analysis.out.mr_ivw_results_txt, 1) // Add mediator type and MR results from exposure to mediators
                .map{ [it[3], *it] } // Put outcome ID first for the following join
            mediation_triplets = inner_join(mediation_triplets, mr_ivw_analysis.out.mr_ivw_results_txt.map{[it[1], it[2]]}, 1) // Add MR results from mediators to outcome
                .map{ [*it[1..11]] } // Remove the leading outcome ID
        }
        else {
            mediation_triplets = MR_UKBB_pairs
                .combine(Channel.fromPath(file_mediators_exposure))
                .combine(Channel.fromPath(file_mediators_output))
                .combine(Channel.value(params.mediators))
                .combine(Channel.fromPath(file_mr_ivw_exposure))
                .combine(Channel.fromPath(file_mr_ivw_outcome))
        }

        prepare_mediation_data_single(mediation_triplets, instrument_pval_threshold, params.mediator_pval_threshold_from_exposure, params.mediator_pval_threshold_to_outcome)

    emit:
        mediation_data = prepare_mediation_data_single.out
 }
