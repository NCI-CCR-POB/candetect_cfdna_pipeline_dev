############################
###   binwise_frags modeling Snakefile   ###
############################

#########1#########2#########3#########4#########5#########6#########7#########8
#
# A snakefile for binwise prediction

rule binwise_model:
    input:
        ratio = bin_wise_frag + "/binwise_ratio.norm_rows.clean.csv",
    output:
        probs = bin_wise_frag + "/model/probs/{c1}_v_{c2}.csv",
        #coeffs = bin_wise_frag + "/model/coeffs/{c1}_v_{c2}.csv",
    conda:
        ENV_FRAGENTOMICS
    log:
        logdir_binwise + "/model/binwise_model_{c1}_v_{c2}.log"
    params:
        arm_z = bin_wise_frag + "/binwise_armz.clean.csv",
        sample_list = samples,
        model = lambda wildcards: (
            "no"
            if config.get("existing_delfi_model", "no") == "no"
            else os.path.join(
                config["existing_delfi_model"],f"{wildcards.c1}_v_{wildcards.c2}.pkl"
            )),
        out_dir = bin_wise_frag + "/model",
        out_dirprob = bin_wise_frag + "/model/prob",
        out_dircoeffs = bin_wise_frag + "/model/coeffs",
        script = frag_scriptdir + "/binwise_model_mod3.py",
        c1 = lambda wildcards: wildcards.c1,
        c2 = lambda wildcards: wildcards.c2,
        name = "{c1}_v_{c2}",
        log_dir = logdir_binwise + "/model",
    shell:
        """
        mkdir -p {params.out_dir} {params.out_dirprob} {params.out_dircoeffs} {params.log_dir}
        echo "Using existing model from config: {params.model}"
            python3 {params.script} \
            --armz_csv {params.arm_z} \
            --ratio_csv {input.ratio} \
            --samplelist {params.sample_list} \
            --cohort1 {params.c1} \
            --cohort2 {params.c2} \
            --name {params.name} \
            --out_dir {params.out_dir} \
            --existing_model {params.model} \
            &> {log}
        """
rule binwise_apply_cutpoints:
    input:
        probs = bin_wise_frag + "/model/probs/{c1}_v_{c2}.csv",
    output:
        scores = bin_wise_frag + "/scores/{c1}_v_{c2}.csv",
    conda:
        ENV_FRAGENTOMICS
    log:
        logdir_binwise + "/apply_cutpoints/binwise_apply_cutpoints_{c1}_v_{c2}.log"
    params:
        existing_scores_dir = config.get("existing_binwise_scores_dir", "no"),
        script = frag_scriptdir + "/apply_cutpoints.py",
        samplelist = samples,
        log_dir = logdir_binwise + "/apply_cutpoints",
    shell:
        """
        mkdir -p {params.log_dir}
        mkdir -p $(dirname {output.scores})

        python3 {params.script} \
        --probs_csv {input.probs} \
        --scores_csv {params.existing_scores_dir}/scores/$(echo "{wildcards.c1}" | sed 's/health/health/g')_v_{wildcards.c2}.csv \
        --samplelist {params.samplelist} \
        --cohort1 {wildcards.c1} \
        --cohort2 {wildcards.c2} \
        --out_csv {output.scores} \
        &> {log}
        """