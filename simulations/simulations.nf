#!/usr/bin/env nextflow

// Enable DSL2 syntax (introduces modularity and improved data flow manipulation)
nextflow.enable.dsl=2

/*
 * Describe mediators: distributions for heritabilities and QTL counts from eQTLGen cis-eQTL data,
 *                     correlations (Sigma) from Colaus gene expression RPKM data.
 */
process describe_mediators {
    publishDir "${projectDir}/figures/mediator_distributions", pattern: "*.pdf", mode: "copy"
    publishDir "${projectDir}/data", pattern: "mediator_QTLs.txt", mode: "copy"
    publishDir "${projectDir}/data", pattern: "Sigma.RDS", mode: "copy"

    input:
        path mediator_summary_data
        path mediator_individual_data
        path heritability_tools_R

    output:
        path "mediator_heritabilities.pdf", emit: mediator_heritabilities_pdf
        path "mediator_QTLs.pdf", emit: mediator_QTLs_pdf
        path "mediator_QTLs.txt", emit: QTLs_txt
        path "Sigma.RDS", emit: Sigma_RDS

    script:
        """
        Rscript "${projectDir}/scripts/mediator_description.R" \
            --heritability_tools_R="${heritability_tools_R}" \
            --in_mediator_summary_data="${mediator_summary_data}" \
            --in_mediator_individual_data="${mediator_individual_data}" \
            --out_mediator_heritabilities_pdf="mediator_heritabilities.pdf" \
            --out_mediator_QTLs_pdf="mediator_QTLs.pdf" \
            --out_mediator_QTLs_txt="mediator_QTLs.txt" \
            --out_mediator_Sigma_RDS="Sigma.RDS"
        """
}

/*
 * Read and process the simulations scenarios (set of configurations)
 */
process sim_scenarios {
    input:
        path sims
        val N_sims

    output:
        path "sim_scenarios_${N_sims}.txt", emit: scenarios
        path "sims_${N_sims}.txt", emit: sims

    script:
        """
        Rscript ${projectDir}/scripts/sim_scenarios.R \
            --N_sims=${N_sims} \
            --in_sims_xlsx=${sims} \
            --out_scenarios_txt="sim_scenarios_${N_sims}.txt" \
            --out_sims_txt="sims_${N_sims}.txt"
        """
}

/*
 * Create true data and estimates for simulation scenarios
 */
process sim_generator {
    input:
        tuple val(seed), val(start), val(end), val(sim_name), val(m), val(nX), val(nZ), val(nY), val(h2X), val(total_causal), val(p_direct), val(cor_mediated), 
            val(var_explained_ZtoY), val(k), val(k_sig_perc), val(p_value_mediators), val(Sigma_type), val(sigmaC2), val(sigmac2), val(sigmab2), val(out_simulation_RDS)
        path mediator_QTLs_txt
        path mediator_Sigma_RDS
        path heritability_tools_R
        path ML_methods_R

    output:
        path "${out_simulation_RDS}"

    script:
        """
        Rscript "${projectDir}/scripts/sim_generator.R" \
            --m="${m}" \
            --nX="${nX}" \
            --nZ="${nZ}" \
            --nY="${nY}" \
            --h2X="${h2X}" \
            --total_causal="${total_causal}" \
            --p_direct="${p_direct}" \
            --cor_mediated="${cor_mediated}" \
            --var_explained_ZtoY="${var_explained_ZtoY}" \
            --k="${k}" \
            --k_sig_perc="${k_sig_perc}" \
            --p_value_mediators="${p_value_mediators}" \
            --Sigma_type="${Sigma_type}" \
            --sigmaC2="${sigmaC2}" \
            --sigmac2="${sigmac2}" \
            --sigmab2="${sigmab2}" \
            --seed="${seed}" \
            --start="${start}" \
            --end="${end}" \
            --heritability_tools_R="${heritability_tools_R}" \
            --ML_methods_R="${ML_methods_R}" \
            --in_mediator_QTLs_txt="${mediator_QTLs_txt}" \
            --in_mediator_Sigma_RDS="${mediator_Sigma_RDS}" \
            --out_simulation_RDS="${out_simulation_RDS}"
        """
}

/*
 * Simulations with the naive, original likelihood and integrated likelihood methods
 */
process simulate {
    publishDir { "${projectDir}/results/" + params.type + "/" + "${simulation_RDS}".split("_")[3] }, mode: "copy"

    input:
        tuple path(simulation_RDS), val(method)
        path ML_methods_R

    output:
        path "${method}_${simulation_RDS}"

    script:
        """
        Rscript "${projectDir}/scripts/simulate.R" \
            --method="${method}" \
            --ML_methods_R="${ML_methods_R}" \
            --in_simulation_RDS=${simulation_RDS} \
            --out_simulation_RDS="${method}_${simulation_RDS}"
        """
}

workflow {
    // Set the input channels
    heritability_tools_R = Channel.fromPath("../scripts/heritability_tools.R").first()
    ML_methods_R = Channel.fromPath("../scripts/ML_methods.R").first()
    methods_ch = Channel.of("naive", "original", "integrated")

    // Read and process the simulations scenarios from xlsx to channel
    sims_ch = Channel.fromPath("data/raw/conf/" + params.type + "/[!~]*.xlsx")

    sim_scenarios(sims_ch, params.N_sims)
    scenarios_ch = sim_scenarios.out.scenarios
        .splitCsv(header: true, sep: "\t", strip: true)
        .map{row -> [row.seed, row.start, row.end, row.sim_name, row.m, row.nX, row.nZ, row.nY, row.h2X, row.total_causal, row.p_direct, row.cor_mediated, 
                     row.var_explained_ZtoY, row.k, row.k_sig_perc, row.p_value_mediators, row.Sigma_type, row.sigmaC2, row.sigmac2, row.sigmab2, row.out_simulation_RDS]}

    // Use eQTLGen cis-eQTL summary data and Colaus gene expression RPKM data to describe mediators approximately
    mediator_summary_data = Channel.fromPath("data/raw/mediators/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt")
    mediator_individual_data = Channel.fromPath("data/raw/mediators/all.colaus.rpkms.csv")
    describe_mediators(mediator_summary_data, mediator_individual_data, heritability_tools_R)

    // Generate the simulation scenario (corresponding true data and effect estimates)
    sim_generator(scenarios_ch, describe_mediators.out.QTLs_txt.first(), describe_mediators.out.Sigma_RDS.first(), heritability_tools_R, ML_methods_R)

    // Run the simulations
    simulate(sim_generator.out.combine(methods_ch), ML_methods_R)
}
