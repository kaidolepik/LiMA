#!/usr/bin/env nextflow

// Enable DSL2 syntax (introduces modularity and improved data flow manipulation)
nextflow.enable.dsl=2

/*
 * Convert xlsx to csv
 */
process xlsx_to_csv {
    input:
        path xlsx

    output:
        path "${xlsx.baseName}.csv"

    script:
        """
        #!/usr/bin/env Rscript

        readr::write_delim(openxlsx::read.xlsx("${xlsx}"), "${xlsx.baseName}.csv", delim = "\t")
        """
}

/*
 * The mediation analysis: take harmonized data, find intersection of variants, filter by p-value and Steiger, clump by instrument strength, and perform the analysis
 */
process prepare_summary_statistics {
    input:
        tuple val(exposure_ID), val(exposure_type), val(outcome_ID), val(outcome_type), val(mediator_type), 
            val(instrument_pval_threshold), path(mediation_data_RDS), path(mr_results_txt), val(clump)
        val(mediator_cor_threshold)
        path reference_dir
        path instrument_selection_tools_R
        path ML_methods_R

    output:
        tuple val(exposure_ID), val(exposure_type), val(outcome_ID), val(outcome_type), val(mediator_type), 
            val(instrument_pval_threshold), val(mediator_cor_threshold), path("summary_statistics_withIDs.RDS")

    script:
        """
        Rscript ${projectDir}/bin/prepare_summary_statistics.R \
            --exposure_ID=${exposure_ID} \
            --outcome_ID=${outcome_ID} \
            --instrument_pval_threshold=${instrument_pval_threshold} \
            --mediator_cor_threshold=${mediator_cor_threshold} \
            --clump=${clump} \
            --reference_dir=${reference_dir} \
            --instrument_selection_tools_R=${instrument_selection_tools_R} \
            --ML_methods_R=${ML_methods_R} \
            --in_mediation_data_RDS=${mediation_data_RDS} \
            --in_mr_results_txt=${mr_results_txt} \
            --out_summary_statistics_RDS="summary_statistics_withIDs.RDS"
        """
}

/*
 * Perform the mediation analysis
 */
process mediation_analysis {
    input:
        tuple val(exposure_ID), val(exposure_type), val(outcome_ID), val(outcome_type), val(mediator_type), 
            val(instrument_pval_threshold), val(mediator_cor_threshold), path(summary_statistics_RDS), val(mediation_method)
        val pleiotropy
        path ML_methods_R

    output:
        tuple val(exposure_ID), val(exposure_type), val(outcome_ID), val(outcome_type), val(mediator_type), 
            val(instrument_pval_threshold), val(mediator_cor_threshold), path(summary_statistics_RDS), val(mediation_method),
            path("${exposure_ID}_${outcome_ID}_${mediation_method}_${mediator_type}_${instrument_pval_threshold}_${mediator_cor_threshold}_${pleiotropy}.txt")

    script:
        """
        Rscript ${projectDir}/bin/mediation_analysis.R \
            --exposure_ID=${exposure_ID} \
            --exposure_type=${exposure_type} \
            --outcome_ID=${outcome_ID} \
            --outcome_type=${outcome_type} \
            --mediator_type=${mediator_type} \
            --instrument_pval_threshold=${instrument_pval_threshold} \
            --mediator_cor_threshold=${mediator_cor_threshold} \
            --mediation_method=${mediation_method} \
            --pleiotropy=${pleiotropy} \
            --ML_methods_R=${ML_methods_R} \
            --in_summary_statistics_RDS=${summary_statistics_RDS} \
            --out_mediation_results_txt="${exposure_ID}_${outcome_ID}_${mediation_method}_${mediator_type}_${instrument_pval_threshold}_${mediator_cor_threshold}_${pleiotropy}.txt"
        """
}

/*
 * Concatenate the results of mediation analyses
 */
process concat_mediation_analyses {
    publishDir "results/" + (params.test ? "test/" : "") + "mediation_by_method", mode: "copy"

    input:
        tuple val(out_filename), path(mediation_analyses)

    output:
        path "${out_filename}"

    script:
        """
        #!/usr/bin/env Rscript

        library(magrittr)

        list.files() %>%
            lapply(data.table::fread) %>%
            data.table::rbindlist(fill = TRUE) %>%
            readr::write_delim("${out_filename}", delim = "\t")
        """
}

// Include workflows
include { mr_ivw; inner_join } from "./modules/mr_ivw.nf"
include { multi_file_mediator_pipeline } from "./modules/multi_file_mediator_pipeline.nf"
include { single_file_mediator_pipeline } from "./modules/single_file_mediator_pipeline.nf"

workflow {
    // Construct data channels
    reference_dir = Channel.fromPath("data/raw/reference/UK10K").first()
    instrument_selection_tools_R = Channel.fromPath("../scripts/instrument_selection_tools.R").first()
    ML_methods_R = Channel.fromPath("../scripts/ML_methods.R").first()
    instrument_pval_threshold = Channel.value(params.instrument_pval_threshold) // For backward compatibility
    mediation_method = Channel.of("naive", "original", "integrated_fixed_contained")

    UKBB_annotations_ch = Channel.fromPath("data/raw/conf/UKBB_Neale_traits_" + (params.test ? "test" : "ALL") + ".xlsx").first()
    UKBB_files_ch = Channel.fromPath(["data/raw/UKBB_Neale/**gwas.imputed_v3.both_sexes*"])
        .map{file -> [file.baseName.split("(_irnt|\\.)")[0], file.name.split("\\.")[-3], file]} // trait_ID, trait_sex, filename

    // Construct a channel of UKBB traits
    xlsx_to_csv(UKBB_annotations_ch)
    UKBB_traits = xlsx_to_csv.out
        .splitCsv(header: true, sep: "\t", strip: true)
        .map{row -> [row.Field_ID, row.Sex, "UKBB_" + row.Type]}
        .join(UKBB_files_ch, by: [0, 1])
        .map{[it[0], it[2], it[3]]} // trait_ID, trait_type, filename

    // Construct exposure-outcome UKBB pairs for MR
    UKBB_pairs = UKBB_traits.combine(UKBB_traits) // exposure_ID, exposure_type, exposure_filename, outcome_ID, outcome_type, outcome_filename
        .filter{it[0] != it[3] && it[1] != "UKBB_binary"} // Do not allow exposure and outcome to be the same, nor a binary exposure
        .map{[exposure_ID:it[0], exposure_type:it[1], exposure_filename:it[2], outcome_ID:it[3], outcome_type:it[4], outcome_filename:it[5]]}

    // Perform MR analysis
    mr_ivw(UKBB_pairs, reference_dir, instrument_selection_tools_R, instrument_pval_threshold)

    // Filter UKBB exposure-outcome pairs by genome-wide Bonferroni and allow only one-sided significance
    MR_UKBB_pairs = mr_ivw.out.mr_ivw_results
        //.combine(mr_ivw.out.mr_ivw_results.count())
        .filter{Double.valueOf(it[6]) <= 0.05 / 1000000}
        .map{[*it[0..5], [it[0], it[1]].sort().join("-")]}  // add a key made of trait pairs
        .groupTuple(by: 6)                                  // group trait pairs by the key
        .filter{it[0].size() < 2}                           // remove trait pairs that are significant in both directions
        .transpose()                                        // remove the grouping
        .map{[it[0], it[2], it[4], it[1], it[3], it[5]]}    // exposure_ID, exposure_type, exposure_filename, outcome_ID, outcome_type, outcome_filename

    // Prepare mediation data
    if (params.mediators == "Lotta_et_al_2021") {
        multi_file_mediator_pipeline(UKBB_pairs, MR_UKBB_pairs, reference_dir, instrument_selection_tools_R, instrument_pval_threshold)

        mediation_data = multi_file_mediator_pipeline.out.mediation_data
    }
    else if (params.mediators == "INTERVAL_plasma_proteins") {
        single_file_mediator_pipeline(MR_UKBB_pairs, reference_dir, instrument_selection_tools_R, instrument_pval_threshold)

        mediation_data = single_file_mediator_pipeline.out.mediation_data
    }

    // Prepare summary statistics
    prepare_summary_statistics(mediation_data, params.mediator_cor_threshold, reference_dir, instrument_selection_tools_R, ML_methods_R)

    // Run mediation analyses, then collect the results together and save
    mediation_data = prepare_summary_statistics.out.combine(mediation_method)
    mediation_analysis(mediation_data, params.pleiotropy, ML_methods_R)

    // Concatenate and publish the mediation analyses
    analyses_groups = mediation_analysis.out
        .map{ [params.mediators + "_" + params.instrument_pval_threshold + "_" + params.mediator_pval_threshold_from_exposure + "_" + params.mediator_pval_threshold_to_outcome + "_" + params.mediator_cor_threshold + "_" + params.pleiotropy + ".txt", it[9]] }
        .groupTuple(by: 0)
    concat_mediation_analyses(analyses_groups)
}
