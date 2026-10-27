#########1#########2#########3#########4#########5#########6#########7#########8
#                                                                              #
#    End_motif Features for Fragmentomic Analysis of Cell-free DNA Whole Genome Sequencing  #
#                                                                              #
#########1#########2#########3#########4#########5#########6#########7#########8

rule bam_to_motif:
    input: frag_bams + "/{library}_filt.bam",
    log: logdir_motif + "/motifs/bam_to_motif/{library}_bam_to_motif.log",
    output: motif_dir + "/{library}_motifs.tsv",
    conda:
        ENV_FRAGENTOMICS,
    params:
        fasta = genome_fasta,
        script = frag_scriptdir + "/extract_5p_motifs.sh",
        threads = getthreads("bam_to_motif"),
        out_dir = motif_dir,
        log_dir = logdir_motif + "/motifs/bam_to_motif",
    shell:
        """
        mkdir -p {params.log_dir}
        mkdir -p {params.out_dir} 
        #module load python/3.9 samtools bedtools
        unset PYTHONPATH
        bash {params.script} \
        {input} \
        {params.fasta} \
        {output} \
        {params.threads} \
        > {log} 2>&1
        """
# Merge short and long fragment counts by genomic poistion for all libraries
rule motif_merge:
    input: expand(motif_dir + "/{library}_motifs.tsv",  library = FRAG_LIBS),
    log: logdir_motif + "/motifs/motif_merge/motif_merge.log",
    output:  
        combined_counts = data_dir + "/analysis/cfdna_wgs/features/combined_motif_counts_matrix.tsv",
        combined_fractions = data_dir + "/analysis/cfdna_wgs/features/combined_motif_fractions_matrix.tsv",
    conda:
        ENV_FRAGENTOMICS,
    params:
        motif_tsv_dir = motif_dir,
        script = frag_scriptdir + "/aggregate_motif_counts.py",
        out_dir = data_dir + "/analysis/cfdna_wgs/features",
        log_dir = logdir_motif + "/motifs/motif_merge",
    shell:
        """
        mkdir -p {params.log_dir}
        mkdir -p {params.out_dir} 
        python {params.script} \
        {params.motif_tsv_dir} \
        {output.combined_counts} \
        {output.combined_fractions} \
        &> {log}
        """
rule mds:
    input: data_dir + "/analysis/cfdna_wgs/features/combined_motif_fractions_matrix.tsv",
    log: logdir_motif + "/motifs/mds/motif_diversity_score.log",
    output:  
        mdsRDS = data_dir + "/analysis/cfdna_wgs/features/mds.RData",
        CSVout = data_dir + "/analysis/cfdna_wgs/features/motif_diversity_scores.csv",
    conda:
        ENV_FRAGENTOMICS,
    params:
        annotation_table = config["samples"],
        script = frag_scriptdir + "/motif_diversity_score.R",
        out_dir = data_dir + "/analysis/cfdna_wgs/features",
        log_dir = logdir_motif + "/motifs/mds",
    shell:
        """
        mkdir -p {params.log_dir}
        mkdir -p {params.out_dir} 
        Rscript  {params.script} \
        {input} \
        {params.annotation_table} \
        {output.mdsRDS} \
        {output.CSVout} \
        &> {log}
        """

