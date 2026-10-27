#########1#########2#########3#########4#########5#########6#########7#########8
#                                                                              #
#     Create a master dataframe with fragmentomic outputs (R version)          #
#                                                                              #
#########1#########2#########3#########4#########5#########6#########7#########8

def nmf_score_path(c1, c2):
    nc = BEST_NC.get((c1, c2), min_components)
    return f"{nmf_dir}/scores/{c1}_v_{c2}.n={nc}.csv"

rule combine_datasets:
    input:
        tf = lambda wildcards: data_dir + "/analysis/cfdna_wgs/ichorCNA/tf.csv" if os.path.exists(data_dir + "/analysis/cfdna_wgs/ichorCNA/tf.csv") else [],
        fragle = lambda wildcards: data_dir + "/FRAGLE/master_compiled_csv_2025-04-22.csv" if os.path.exists(data_dir + "/FRAGLE/master_compiled_csv_2025-04-22.csv") else [],
        #fragle = lambda wildcards: data_dir + "/FRAGLE/master_compiled_csv_2025-04-13.csv" if os.path.exists(data_dir + "/FRAGLE/master_compiled_csv_2025-04-13.csv") else [],
        mds = lambda wildcards: data_dir + "/analysis/cfdna_wgs/features/motif_diversity_scores.csv" if os.path.exists(data_dir + "/analysis/cfdna_wgs/features/motif_diversity_scores.csv") else [],
        fprofile_auc = lambda wildcards: (
            f"{data_dir}/analysis/cfdna_wgs/features/fprofile_AUCROC/"
            f"fprofile_AUC_{wildcards.c1}_v_{wildcards.c2}.csv"
            if os.path.exists(
                f"{data_dir}/analysis/cfdna_wgs/features/fprofile_AUCROC/"
                f"fprofile_AUC_{wildcards.c1}_v_{wildcards.c2}.csv"
            )
            else []
        ),
        fprofiles = lambda wildcards: data_dir + "/analysis/cfdna_wgs/features/fprofiles.tsv" if os.path.exists(data_dir + "/analysis/cfdna_wgs/features/fprofiles.tsv") else [],
        #nmf = lambda wildcards: (nmf_dir + f"/scores/{wildcards.c1}_v_{wildcards.c2}.n={BEST_NMF_NC}.csv" if os.path.exists(nmf_dir + f"/scores/{wildcards.c1}_v_{wildcards.c2}.n={BEST_NMF_NC}.csv") else []),
        nmf = lambda wildcards: nmf_score_path(wildcards.c1, wildcards.c2),
        binwise = data_dir + "/analysis/cfdna_wgs/bin_wise_fragmentomics/scores/{c1}_v_{c2}.csv"
    output:
        combined_dataset = data_dir + "/analysis/cfdna_wgs/combined_dataset_{c1}_v_{c2}.csv",
        status_file = data_dir + "/analysis/cfdna_wgs/combined_dataset_status_{c1}_v_{c2}.txt"
    conda:
        ENV_FRAGENTOMICS,
    log:
        log_dir_combine + "/combine_datasets_{c1}_v_{c2}.log"
    params:
        tf_arg = lambda w, input: f"--tf='{input.tf}'" if input.tf else "",
        fragle_arg = lambda w, input: f"--fragle='{input.fragle}'" if input.fragle else "",
        mds_arg = lambda w, input: f"--mds='{input.mds}'" if input.mds else "",
        fprofile_auc_arg = lambda w, input: f"--fprofile_auc='{input.fprofile_auc}'" if input.fprofile_auc else "",
        fprofiles_arg = lambda w, input: f"--fprofiles='{input.fprofiles}'" if input.fprofiles else "",
        nmf_arg = lambda w, input: f"--nmf='{input.nmf}'" if input.nmf else "",
        binwise_arg = lambda w, input: f"--binwise='{input.binwise}'" if input.binwise else ""
    shell:
        """
        Rscript scripts/combine_datasets_pred.R \
            --config='{config[samples]}' \
            --data_dir='{data_dir}' \
            --cohort1='{wildcards.c1}' \
            --cohort2='{wildcards.c2}' \
            --normal_control='{normal_control}' \
            --output='{output.combined_dataset}' \
            --status='{output.status_file}' \
            --log='{log}' \
            {params.tf_arg} \
            {params.fragle_arg} \
            {params.mds_arg} \
            {params.fprofile_auc_arg} \
            {params.fprofiles_arg} \
            {params.nmf_arg} \
            {params.binwise_arg} 2>&1 | tee -a {log}
        """

rule plot_subject_trajectories:
    input:
        combined_dataset = data_dir + "/analysis/cfdna_wgs/combined_dataset_{c1}_v_{c2}.csv"
    output:
        plot = plot_dir_combine + "/subject_trajectories_{c1}_v_{c2}.png"
    conda:
        ENV_FRAGENTOMICS,
    log:
        log_dir_combine + "/plot_subject_trajectories_{c1}_v_{c2}.log"
    shell:
        """
        mkdir -p {plot_dir_combine}
        Rscript scripts/plot_subject_trajectories.R \
            --input {input.combined_dataset} \
            --output {output.plot} \
            --cohort1 {wildcards.c1} \
            --cohort2 {wildcards.c2} \
            --normal_control {normal_control} \
            --log {log} 2>&1 | tee -a {log}
        """

rule plot_auc_comparison:
    input:
        combined_dataset = data_dir + "/analysis/cfdna_wgs/combined_dataset_{c1}_v_{c2}.csv"
    output:
        plot = plot_dir_combine + "/auc_comparison_{c1}_v_{c2}.png"
    conda:
        ENV_FRAGENTOMICS,
    log:
        log_dir_combine + "/plot_auc_comparison_{c1}_v_{c2}.log"
    shell:
        """
        mkdir -p {plot_dir_combine}
        Rscript scripts/plot_auc.R \
            --input {input.combined_dataset} \
            --output {output.plot} \
            --cohort1 {wildcards.c1} \
            --cohort2 {wildcards.c2} \
            --normal_control {normal_control} \
            --log {log} 2>&1 | tee -a {log}
        """

rule plot_metric_distributions:
    input:
        combined_dataset = data_dir + "/analysis/cfdna_wgs/combined_dataset_{c1}_v_{c2}.csv"
    output:
        directory(plot_dir_combine + "/distributions_{c1}_v_{c2}")
    conda:
        ENV_FRAGENTOMICS,
    log:
        log_dir_combine + "/plot_metric_distributions_{c1}_v_{c2}.log"
    shell:
        """
        mkdir -p {output}
        Rscript scripts/plot_metric_distributions.R \
            --input {input.combined_dataset} \
            --output_dir {output} \
            --cohort1 {wildcards.c1} \
            --cohort2 {wildcards.c2} \
            --normal_control {normal_control} \
            --log {log} 2>&1 | tee -a {log}
        """

rule plot_voting:
    input:
        combined_dataset = data_dir + "/analysis/cfdna_wgs/combined_dataset_{c1}_v_{c2}.csv",
        tf = data_dir + "/analysis/cfdna_wgs/ichorCNA/tf.csv" ,
        binwise_score = data_dir + "/analysis/cfdna_wgs/bin_wise_fragmentomics/scores/{c1}_v_{c2}.csv",
        #nmf_score = lambda wildcards: (nmf_dir + f"/scores/{wildcards.c1}_v_{wildcards.c2}.n={BEST_NMF_NC}.csv"),
        nmf_score = lambda wildcards: nmf_score_path(wildcards.c1, wildcards.c2),
    output:
        data_dir + "/analysis/cfdna_wgs/plots/voting/roc_calls_heatmap_{c1}_v_{c2}.pdf",
        data_dir + "/analysis/cfdna_wgs/plots/voting/roc_calls_heatmap_{c1}_v_{c2}.png",
        data_dir + "/analysis/cfdna_wgs/plots/voting/roc_calls_heatmap_{c1}_v_{c2}_majority_vote.csv",
    conda:
        ENV_FRAGENTOMICS,
    log:
        log_dir_combine + "/plot_voting_{c1}_v_{c2}.log"
    params:
        output_dir = data_dir + "/analysis/cfdna_wgs/plots/voting"
    shell:
        """
        mkdir -p {params.output_dir}
        Rscript scripts/create_voting_heatmap_pred.R \
            --config='{config[samples]}' \
            --input='{input.combined_dataset}' \
            --data_dir='{data_dir}' \
            --cohort1='{wildcards.c1}' \
            --cohort2='{wildcards.c2}' \
            --normal_control='{normal_control}' \
            --output_dir='{params.output_dir}' \
            --log='{log}' \
            --nmf_score='{input.nmf_score}' \
            --binwise_score='{input.binwise_score}' \
            --ichor_tf='{input.tf}' 2>&1 | tee -a {log}
        """

rule html_report:
    input:
        tf = data_dir + "/analysis/cfdna_wgs/ichorCNA/tf.csv",
        voting_png = expand(
            data_dir + "/analysis/cfdna_wgs/plots/voting/roc_calls_heatmap_{c1}_v_{c2}.png",
            zip,
            c1 = COHORT_1,
            c2 = COHORT_2
        ),
        voting_csv = expand(
            data_dir + "/analysis/cfdna_wgs/plots/voting/roc_calls_heatmap_{c1}_v_{c2}_majority_vote.csv",
            zip,
            c1 = COHORT_1,
            c2 = COHORT_2
        )
    output:
        html = data_dir + "/analysis/cfdna_wgs/report/cfdna_wgs_report.html"
    conda:
        ENV_FRAGENTOMICS,
    log:
        log_dir_combine + "/html_report.log"
    params:
        script = frag_repo + "/scripts/create_html_report_pred.py",
        voting_dir = data_dir + "/analysis/cfdna_wgs/plots/voting",
        ichor_dir = data_dir + "/analysis/cfdna_wgs/ichorCNA",
        report_dir = data_dir + "/analysis/cfdna_wgs/report"
    shell:
        """
        mkdir -p {params.report_dir}
        python {params.script} \
            --data-dir '{data_dir}' \
            --output '{output.html}' \
            --tf '{input.tf}' \
            --voting-dir '{params.voting_dir}' \
            --ichor-dir '{params.ichor_dir}' \
            > {log} 2>&1
        """
