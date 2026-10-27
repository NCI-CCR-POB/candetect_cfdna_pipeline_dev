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
        #coeffs = nmf_dir + "/coeffs/{fnames}.n={nc}.csv",
        model = nmf_dir + "/model/{fnames}.n={nc}.pkl",
    conda:
        ENV_FRAGENTOMICS,
    log:
        logdir_nmf + "/nmf/model/{fnames}.n={nc}.log",
    params:
        out_dir = nmf_dir,
        log_dir = logdir_nmf + "/nmf/model",
        script = frag_scriptdir + "/nmf_model_mod3.py",
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
rule nmf_apply_cutpoints:
    input:
        probs = nmf_dir + "/probs/{fnames}.n={nc}.csv",
    output:
        scores = nmf_dir + "/scores/{fnames}.n={nc}.csv",
    conda:
        ENV_FRAGENTOMICS,
    log:
        logdir_nmf + "/nmf/apply_cutpoints/{fnames}.n={nc}.log"
    params:
        existing_nmf_dir = config.get("existing_nmf_scores_dir", "no"),
        script = frag_scriptdir + "/apply_cutpoints.py",
        samplelist = samples,
        log_dir = logdir_nmf + "/nmf/apply_cutpoints",
    shell:
        """
        mkdir -p {params.log_dir}
        mkdir -p $(dirname {output.scores})

        python3 {params.script} \
        --probs_csv {input.probs} \
        --scores_csv {params.existing_nmf_dir}/scores/$(echo "{wildcards.fnames}" | sed 's/health/health/g').n={wildcards.nc}.csv \
        --samplelist {params.samplelist} \
        --cohort1 $(echo "{wildcards.fnames}" | sed 's/_v_.*//') \
        --cohort2 $(echo "{wildcards.fnames}" | sed 's/.*_v_//') \
        --out_csv {output.scores} \
        &> {log}
        """
