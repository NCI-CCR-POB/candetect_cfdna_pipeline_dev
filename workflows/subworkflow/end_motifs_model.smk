#########1#########2#########3#########4#########5#########6#########7#########8
#                                                                              #
#     Fragmentomic Analysis of Cell-free DNA Whole Genome Sequencing           #
#                                                                              #
#########1#########2#########3#########4#########5#########6#########7#########8

rule f_profile:
    input: data_dir + "/analysis/cfdna_wgs/features/combined_motif_fractions_matrix.tsv",
    log: logdir_motif + "/motifs/f_profile/f-profile.log",
    output: 
        fprofile_per_lib = data_dir + "/analysis/cfdna_wgs/features/fprofiles.tsv",
        motif_per_fprofile = data_dir + "/analysis/cfdna_wgs/features/motif_per_fprofile.tsv",
        model_dir = directory(data_dir + "/analysis/cfdna_wgs/features/fprofile_model/"),
    conda:
        ENV_FRAGENTOMICS,
    params:
        script = frag_scriptdir + "/fprofiles_pkl.py",
        log_dir = logdir_motif + "/motifs/f_profile",
        existing_model = config["existing_motif_model"],
    shell:
        r"""
            set -euo pipefail
            mkdir -p {params.log_dir} {output.model_dir}

            # If loading an existing model
            if [[ "{params.existing_model}" != "no" ]]; then
                model_file="{params.existing_model}"
                model_out="{output.model_dir}/fprofile_model.pkl"

                if [[ -f "$model_file" ]]; then
                    cp "$model_file" "$model_out"
                    python {params.script} \
                        --data_file {input} \
                        --output_fprof_per_lib {output.fprofile_per_lib} \
                        --output_motif_per_fprof {output.motif_per_fprofile} \
                        --load_model "$model_out" \
                        > {log} 2>&1
                else
                    echo "Specified model file not found: $model_file" >&2
                    exit 1
                fi
            else
                # Train a new model
                python {params.script} \
                    --data_file {input} \
                    --output_fprof_per_lib {output.fprofile_per_lib} \
                    --output_motif_per_fprof {output.motif_per_fprofile} \
                    --save_model {output.model_dir}/fprofile_model.pkl \
                    > {log} 2>&1
            fi
            """


rule f_profile_plot:
    input: data_dir + "/analysis/cfdna_wgs/features/fprofiles.tsv",
    log: logdir_motif + "/f-profile_plot.log",
    output: 
        data_dir + "/analysis/cfdna_wgs/features/fprofiles.pdf",
    conda:
        ENV_FRAGENTOMICS,
    params:
        script = frag_scriptdir + "/fprofile_stackedbarplot.py",
        plot_title = config["fprofile_stacked_title"],
        metadata_file = metadata_file, #Note, if not using e.g. different cycles, set as "none" in config
        colorscheme = "viridis",
    shell:
        """
        #module load python/3.9
        python {params.script} --input_file {input} --metadata_file {params.metadata_file} \
        --output_file {output} --title "{params.plot_title}" --colormap {params.colorscheme} \
        > {log} 2>&1
	    """      

rule fprofile_roc_plots:
    input: data_dir + "/analysis/cfdna_wgs/features/fprofiles.tsv",
    log: logdir_motif + "/motifs/roc_plots/roc_f-profiles_{c1}_v_{c2}.log",
    output:
        plot_out = data_dir + "/analysis/cfdna_wgs/features/fprofile_AUCROC/fprofile_ROC_{c1}_v_{c2}.pdf",
        csv_out  = data_dir + "/analysis/cfdna_wgs/features/fprofile_AUCROC/fprofile_AUC_{c1}_v_{c2}.csv",
    conda:
        ENV_FRAGENTOMICS,
    params:
        libraries_tsv = config["samples"],
        script = frag_scriptdir + "/motif_ROC_plot.R",
        out_dir = data_dir + "/analysis/cfdna_wgs/features",
        log_dir = logdir_motif + "/motifs/roc_plots",
    shell:
        """
        mkdir -p {params.log_dir}
        mkdir -p {params.out_dir}
        Rscript {params.script} \
        {input} \
        {params.libraries_tsv} \
        {output.plot_out} \
        {output.csv_out} \
        &> {log}
        """
