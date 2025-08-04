#!/usr/bin/env nextflow

// Enable DSL2 syntax (introduces modularity and improved data flow manipulation)
nextflow.enable.dsl=2

 /*
 * Harmonize the data in the beginning as it's needed by all downstream analyses
 * (effect allele is the A1 allele – usually minor allele – of the reference UK10K)
 */
process harmonize_data {
    input:
        tuple val(trait_ID), val(trait_type), path(in_trait_gz)
        path reference_dir
        path instrument_selection_tools_R

    output:
        tuple val("${trait_ID}"), path("harmonized_${trait_ID}.txt.gz"), emit: trait

    script:
        """
        Rscript ${projectDir}/scripts/harmonize_data.R \
            --trait_ID=${trait_ID} \
            --trait_type=${trait_type} \
            --in_trait_gz=${in_trait_gz} \
            --reference_dir=${reference_dir} \
            --instrument_selection_tools_R=${instrument_selection_tools_R} \
            --out_trait_gz="harmonized_${trait_ID}.txt.gz"
        """
}

/*
 * The MR IVW analysis: take harmonized data, find intersection of variants, filter by exposure p-value and Steiger, clump, and perform the analysis
 */
process mr_ivw_analysis {
    input:
        tuple val(exposure_ID), val(outcome_ID), val(exposure_type), val(outcome_type), path(in_exposure_gz), path(in_outcome_gz), val(instrument_pval_threshold)
        path reference_dir
        path instrument_selection_tools_R

    output:
        tuple val(exposure_ID), val(outcome_ID), path("${exposure_type}:${outcome_type}:${instrument_pval_threshold}.txt"), emit: mr_ivw_results_txt

    script:
        """
        Rscript ${projectDir}/scripts/mr_ivw_analysis.R \
            --exposure_ID=${exposure_ID} \
            --in_exposure_gz=${in_exposure_gz} \
            --outcome_ID=${outcome_ID} \
            --in_outcome_gz=${in_outcome_gz} \
            --reference_dir=${reference_dir} \
            --instrument_selection_tools_R=${instrument_selection_tools_R} \
            --instrument_pval_threshold=${instrument_pval_threshold} \
            --out_mr_ivw_results_txt="${exposure_type}:${outcome_type}:${instrument_pval_threshold}.txt"
        """
}

def inner_join(duplicated_keys_ch, unique_keys_ch, N_keys) {
    unique_keys_ch.cross(duplicated_keys_ch)
        .map{ [*it[0][0..N_keys-1], *it[1][N_keys..-1], *it[0][N_keys..-1]] } // [keys, all other items of duplicated_keys_ch, all other items of unique_keys_ch]
}

workflow mr_ivw {
    take:
        trait_pairs_ch
        reference_dir
        instrument_selection_tools_R
        instrument_pval_threshold

    main:
        // Harmonize the trait data (both exposures and outcomes)
        traits_ch = trait_pairs_ch.map{[it.exposure_ID, it.exposure_type, it.exposure_filename]}
            .mix(trait_pairs_ch.map{[it.outcome_ID, it.outcome_type, it.outcome_filename]})
            .unique()

        harmonize_data(traits_ch, reference_dir, instrument_selection_tools_R)

        // To perform MR, add harmonized data to trait pairs together with the instrument P-value threshold
        MR_trait_pairs_ch = inner_join(trait_pairs_ch.map{[it.exposure_ID, it.outcome_ID, it.exposure_type, it.outcome_type]}, harmonize_data.out.trait, 1) // Add exposure harmonized file by exposure ID
        MR_trait_pairs_ch = inner_join(MR_trait_pairs_ch.map{[it[1], it[0], it[2], it[3], it[4]]}, harmonize_data.out.trait, 1) // Add outcome harmonized file by outcome ID
            .map{[it[1], it[0], it[2], it[3], it[4], it[5]]} // exposure_ID, outcome_ID, exposure_type, outcome_type, exposure_filename, outcome_filename
            .combine(instrument_pval_threshold)

        mr_ivw_analysis(MR_trait_pairs_ch, reference_dir, instrument_selection_tools_R)

        // Save and collect the MR IVW results together
        mr_ivw_results = mr_ivw_analysis.out.mr_ivw_results_txt.map{it[2]}
            .collectFile(storeDir: "results/" + (params.test ? "test/" : "") + "mr_ivw", keepHeader: true) { file -> ["${file.baseName.tokenize(':').get(0)}_${file.baseName.tokenize(':').get(1)}_${file.baseName.tokenize(':').get(2)}.txt", file] }
            .splitCsv(header: true, sep: " ", strip: true)
            .map{[it.exposure, it.outcome, it.pval]}
            .filter({it[2] != "NA"})
            .join(MR_trait_pairs_ch, by: [0, 1])
            .map{[it[0], it[1], it[3], it[4], it[5], it[6], it[2]]} // exposure_ID, outcome_ID, exposure_type, outcome_type, exposure_filename, outcome_filename, pval

    emit:
        mr_ivw_results = mr_ivw_results
        mr_ivw_files = mr_ivw_analysis.out.mr_ivw_results_txt
}
