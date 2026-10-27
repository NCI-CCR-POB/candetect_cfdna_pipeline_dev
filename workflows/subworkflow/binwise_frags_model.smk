############################
###   binwise_frags modeling Snakefile   ###
############################

#########1#########2#########3#########4#########5#########6#########7#########8
#
# A snakefile for binwise fragment feature modeling

rule binwise_model:
    input:
        ratio = bin_wise_frag + "/binwise_ratio.norm_rows.clean.csv",
    output:
        probs = bin_wise_frag + "/model/probs/{c1}_v_{c2}.csv",
        coeffs = bin_wise_frag + "/model/coeffs/{c1}_v_{c2}.csv",
    conda:
        ENV_FRAGENTOMICS
    log:
        logdir_binwise + "/model/binwise_model_{c1}_v_{c2}.log"
    params:
        arm_z = bin_wise_frag + "/binwise_armz.clean.csv",
        sample_list = samples,
        model = lambda wildcards: (
            "None"
            if config.get("existing_delfi_model", "no") == "no"
            else os.path.join(
                config["existing_delfi_model"],f"{wildcards.c1}_v_{wildcards.c2}.pkl"
            )),
        out_dir = bin_wise_frag + "/model",
        out_dirprob = bin_wise_frag + "/model/prob",
        out_dircoeffs = bin_wise_frag + "/model/coeffs",
        script = frag_scriptdir + "/binwise_model.py",
        c1 = lambda wildcards: wildcards.c1,
        c2 = lambda wildcards: wildcards.c2,
        name = "{c1}_v_{c2}",
        log_dir = logdir_binwise + "/model",
    shell:
        """
        mkdir -p {params.out_dir} {params.out_dirprob} {params.out_dircoeffs} {params.log_dir}
        if [ "{params.model}" != "None" ]; then
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
        else
            echo "Training new model"
            python3 {params.script} \
            --armz_csv {params.arm_z} \
            --ratio_csv {input.ratio} \
            --samplelist {params.sample_list} \
            --cohort1 {params.c1} \
            --cohort2 {params.c2} \
            --name {params.name} \
            --out_dir {params.out_dir} \
            &> {log}
        fi
        """

rule binwise_aucs:
    input:
        expand(bin_wise_frag + "/model/probs/{c1}_v_{c2}.csv",
        zip,
            c1 = COHORT_1,
            c2 = COHORT_2
        ),
        sample_list = samples,
    output:
        bin_wise_frag + "/aucs/aucs/aucs.csv",
        expand(bin_wise_frag + "/scores/{c1}_v_{c2}.csv",
                zip,
                c1 = COHORT_1,
                c2 = COHORT_2
        ),
    conda:
        ENV_FRAGENTOMICS
    log:
        logdir_binwise + "/auc/binwise_aucs.log"
    params:
        out_dir = bin_wise_frag,
        out_dir_auc = bin_wise_frag + "/aucs",
        out_dir_scores = bin_wise_frag + "/scores",
        script = frag_scriptdir + "/binwise_auc.py",
        probs_dir = bin_wise_frag + "/model/probs",
        log_dir = logdir_binwise + "/auc",
        cohorts = ",".join(config["cohorts"]),
    shell:
        """
        mkdir -p {params.out_dir_auc}
        mkdir -p {params.out_dir_scores}
        mkdir -p {params.log_dir}

        python3 {params.script} \
        --probs_dir {params.probs_dir} \
        --samplelist {input.sample_list} \
        --out_dir {params.out_dir} \
        --cohorts {params.cohorts} \
        &> {log}
        """

rule plot_probs:
    input:
        expand(bin_wise_frag + "/model/probs/{c1}_v_{c2}.csv",
        zip,
            c1 = COHORT_1,
            c2 = COHORT_2
        ),
        sample_list = samples,
    output:
        expand(bin_wise_frag + "/plots/prob_per_repeat/{c1}_v_{c2}.pdf",
        zip,
            c1 = COHORT_1,
            c2 = COHORT_2
        ),
        expand(bin_wise_frag + "/plots/prob_per_sample/{c1}_v_{c2}.pdf",
        zip,
            c1 = COHORT_1,
            c2 = COHORT_2
        ),
    conda:
        ENV_FRAGENTOMICS
    log:
        logdir_binwise + "/plot_probs/plot_probs.log"
    params:
        out_dir = bin_wise_frag,
        log_dir = logdir_binwise + "/plot_probs",
        script  = frag_scriptdir + "/plot_probs.py",
        probs_dir = bin_wise_frag + "/model/probs",
    shell:
        """
        mkdir -p {params.out_dir}
        mkdir -p {params.log_dir}

        python3 {params.script} \
        --probs_dir {params.probs_dir} \
        --samplelist {input.sample_list} \
        --out_dir {params.out_dir} \
        &> {log}
        """

rule plot_coeff_heatmaps:
    input:
        arm_z = bin_wise_frag + "/binwise_armz.clean.csv",
        coeffs = expand(bin_wise_frag + "/model/coeffs/{c1}_v_{c2}.csv", zip, c1=COHORT_1, c2=COHORT_2),
    output:
        bin_wise_frag + "/binwise/plots/coeff_heatmap/binwise_armz.clean.pdf",
    conda:
        ENV_FRAGENTOMICS
    log:
        logdir_binwise + "/coeff_heatmaps/binwise_plot_coeff_heatmaps.log",
    params:
        log_dir = logdir_binwise + "/coeff_heatmaps",
        script  = frag_scriptdir + "/binwise_plot_coefficient_heatmaps.py",
        data_dir = bin_wise_frag,
    shell:
        """
        mkdir -p {params.log_dir}

        python3 {params.script} \
        --data_dir {params.data_dir} \
        --armz_csv {input.arm_z} \
        &> {log}
        """