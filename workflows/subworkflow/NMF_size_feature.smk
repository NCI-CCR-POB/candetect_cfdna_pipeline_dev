#########1#########2#########3#########4#########5#########6#########7#########8
#                                                                              #
#     NMF size Fragmentomic Feature Analysis of Cell-free DNA Whole Genome Sequencing  #
#                                                                              #
#########1#########2#########3#########4#########5#########6#########7#########8

rule raw_histograms:
    input:
        bam = frag_bams + "/{lib}_filt.bam",
    output:
        metrics   = nmf_raw_hist_dir + "/{lib}_metrics.txt",
        histogram = nmf_raw_hist_dir + "/{lib}_histogram.pdf",
    conda:
        ENV_FRAGENTOMICS,
    log:
        logdir_nmf + "/nmf/raw_histograms/{lib}.log",
    params:
        out_dir = nmf_raw_hist_dir,
        log_dir = logdir_nmf + "/nmf/raw_histograms",
        width   = config["histogram_width"],
        #picardpath = config["picardpath"],
    shell:
        """
        mkdir -p {params.out_dir}
        mkdir -p {params.log_dir}

        #module load picard
        #module load R

        picard CollectInsertSizeMetrics \
        -I {input.bam} \
        -O {output.metrics} \
        -H {output.histogram} \
        -W {params.width} \
        &> {log}
        """
rule agg_histograms:
    input:
        metrics = expand(nmf_raw_hist_dir + "/{lib}_metrics.txt", lib=FRAG_LIBS),
        plots = expand(nmf_raw_hist_dir + "/{lib}_histogram.pdf", lib=FRAG_LIBS)
    output:
        counts = nmf_agg_hist_dir + "/count_histogram.csv",
        freqs = nmf_agg_hist_dir + "/freq_histogram.30_700.csv",
    conda:
        ENV_FRAGENTOMICS,
    log:
        logdir_nmf + "/nmf/agg_histograms.log",
    params:
        hist_dir = nmf_raw_hist_dir,
        out_dir = nmf_agg_hist_dir,
        logs_dir = logdir_nmf + "/nmf",
        script = frag_scriptdir + "/agg_histograms.py"
    shell:
        """
        mkdir -p {params.out_dir}
        mkdir -p {params.logs_dir}
        #module load python/3.9
        python3 {params.script} \
        --hist_dir {params.hist_dir} \
        --output {params.out_dir} \
        --verbose \
        &> {log}
        """
