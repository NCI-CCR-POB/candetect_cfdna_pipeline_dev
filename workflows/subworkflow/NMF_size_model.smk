#########1#########2#########3#########4#########5#########6#########7#########8
#                                                                              #
#     NMF size Fragmentomic Analysis of Cell-free DNA Whole Genome Sequencing  #
#                                                                              #
#########1#########2#########3#########4#########5#########6#########7#########8

rule nmf:
    input:
        freq_csv = nmf_agg_hist_dir + "/freq_histogram.30_700.csv",
        counts = nmf_agg_hist_dir + "/count_histogram.csv",

        raw_metrics = expand(nmf_raw_hist_dir + "/{lib}_metrics.txt", lib=FRAG_LIBS),
        raw_plots = expand(nmf_raw_hist_dir + "/{lib}_histogram.pdf", lib=FRAG_LIBS)
    output:
        probs = nmf_dir + "/probs/{fnames}.n={nc}.csv",
        coeffs = nmf_dir + "/coeffs/{fnames}.n={nc}.csv",
        model = nmf_dir + "/model/{fnames}.n={nc}.pkl",
    conda:
        ENV_FRAGENTOMICS,
    log:
        logdir_nmf + "/nmf/model/{fnames}.n={nc}.log",
    params:
        out_dir = nmf_dir,
        log_dir = logdir_nmf + "/nmf/model",
        script = frag_scriptdir + "/nmf_model.py",
        samplelist = samples,
        existing_model = lambda wildcards: (
            "no"
            if config.get("existing_nmf_model", "no") == "no"
            else os.path.join(
                config["existing_nmf_model"],
                f"{wildcards.fnames}.n={wildcards.nc}.pkl"
            )
        ),
    shell:
        """
        mkdir -p {params.out_dir}/{{model,probs,coeffs}} {params.log_dir}
        #module load python/3.9
        python3 {params.script} \
            --freq_csv {input.freq_csv} \
            --samplelist {params.samplelist} \
            --n_components {wildcards.nc} \
            --cohort1 $(echo "{wildcards.fnames}" | sed 's/_v_.*//') \
            --cohort2 $(echo "{wildcards.fnames}" | sed 's/.*_v_//') \
            --out_dir {params.out_dir} \
            --existing_model {params.existing_model} \
            --name {wildcards.fnames}.n={wildcards.nc} \
            &> {log}
        """
rule nmf_aucs:
    input:
        expand(nmf_probs,
                nc = range(min_components,max_components+1)
        ),
        samplelist = samples,
    output:
        nmf_dir + "/aucs/aucs.csv",
        nmf_dir + "/aucs/best_med_aucs.csv",
        expand(nmf_scores,
                nc = range(min_components,max_components+1)
        ),
    conda:
        ENV_FRAGENTOMICS,
    log:
        logdir_nmf + "/nmf/auc/aucs.log"
    params:
        out_dir = nmf_dir,
        log_dir = logdir_nmf + "/nmf/auc",
        script  = frag_scriptdir + "/nmf_aucs.py",
        probs_dir = nmf_dir + "/probs",
        min_components = min_components,
        max_components = max_components,
        cohorts = ",".join(config["cohorts"]),
    shell:
        """
        mkdir -p {params.out_dir}
        mkdir -p {params.log_dir}

        python3 {params.script} \
        --probs_dir {params.probs_dir} \
        --min_components {params.min_components} \
        --max_components {params.max_components} \
        --samplelist {input.samplelist} \
        --out_dir {params.out_dir} \
        --cohorts {params.cohorts} \
        &> {log}
        """

rule plot_nmf_components:
    input:
        model = nmf_dir + "/model/{fnames}.n={nc}.pkl",
        freq_csv = nmf_agg_hist_dir + "/freq_histogram.30_700.csv",
    output:
        plot = nmf_dir + "/plots/components/{fnames}.n={nc}.pdf",
    conda:
        ENV_FRAGENTOMICS,
    log:
        logdir_nmf + "/nmf/plots/components/{fnames}.n={nc}.log"
    params:
        out_dir = nmf_dir + "/plots/components",
        log_dir = logdir_nmf + "/nmf/plots/components",
        script = frag_scriptdir + "/plot_nmf_components.py",
    shell:
        """
        mkdir -p {params.out_dir}
        mkdir -p {params.log_dir}

        python3 {params.script} \
            --model {input.model} \
            --freq_csv {input.freq_csv} \
            --output {output.plot} \
            --cohort1 $(echo "{wildcards.fnames}" | sed 's/_v_.*//') \
            --cohort2 $(echo "{wildcards.fnames}" | sed 's/.*_v_//') \
            --n_components {wildcards.nc} \
            &> {log}
        """
