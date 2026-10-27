##################################################################
###   Integration testing snakefile for WGS cfDNA Processing   ###
##################################################################


rule symlink_inputs:
    input:
        lambda wildcards: frag_libs[frag_libs["library"] == wildcards.library]["file"],
    output:
        read1 = frag_fastqs + "/{library}_raw_"+ PE1 + ".fastq.gz",
        read2 = frag_fastqs + "/{library}_raw_"+ PE2 + ".fastq.gz",
    params:
        outdir = frag_fastqs,
        script = frag_scriptdir+ "/detect_read_pairs.sh",
        log_dir = logdir + "/symlink",
        pair_detection_dir = logdir + "/pair_detection",
    log:
        cmd = logdir + "/symlink/{library}_symlink.log",
        debug = logdir + "/symlink/{library}_symlink_debug.log",
        pair_detection = logdir + "/pair_detection/{library}_read_pair_detection.log",
    run:
        # Create log directories
        os.makedirs(params.log_dir, exist_ok=True)
        os.makedirs(params.pair_detection_dir, exist_ok=True)
        
        # Get source file path
        source_file = frag_libs[frag_libs["library"] == wildcards.library]["file"].iloc[0]
        
        # Log debug information
        with open(log.debug, 'w') as debug_log:
            debug_log.write(f"=== Symlink Debug Info ===\n")
            debug_log.write(f"Library: {wildcards.library}\n")
            debug_log.write(f"Source file: {source_file}\n")
            debug_log.write(f"Output PE1: {output.read1}\n")
            debug_log.write(f"Output PE2: {output.read2}\n")
            debug_log.write(f"Output directory: {params.outdir}\n")
            debug_log.write(f"File exists: {os.path.exists(source_file)}\n")
            debug_log.write(f"File readable: {os.access(source_file, os.R_OK)}\n")
            debug_log.write(f"All libraries in FRAG_LIBS: {FRAG_LIBS}\n")
            debug_log.write(f"Current library in wildcards: {wildcards.library}\n")
            debug_log.write(f"Library exists in FRAG_LIBS: {wildcards.library in FRAG_LIBS}\n")
            debug_log.write(f"=== End Debug Info ===\n")
        
        # Check if source file exists
        if not os.path.exists(source_file):
            error_msg = f"Source file not found: {source_file}"
            with open(log.debug, 'a') as debug_log:
                debug_log.write(f"ERROR: {error_msg}\n")
            raise FileNotFoundError(error_msg)
            
        # Create output directory
        os.makedirs(params.outdir, exist_ok=True)
        
        # Create symlinks
        shell("""
            echo "Creating symlinks for {wildcards.library}" > {log.cmd}
            echo "Source: {input}" >> {log.cmd}
            echo "Destination R1: {output.read1}" >> {log.cmd}
            echo "Destination R2: {output.read2}" >> {log.cmd}
            
            {params.script} \
            {input} \
            {output.read1} \
            {output.read2} \
            {params.outdir} \
            >> {log.pair_detection} 2>&1
            
            echo "Symlink creation completed" >> {log.cmd}
            ls -l {output.read1} {output.read2} >> {log.cmd} 2>&1
        """)
#rule frag_index:
#    benchmark: bench_dir + "/frag_index_"+ fasta_base +".benchmark.txt",
#    input: genome_fasta,
#    log: logdir + "/reference_index/frag_index_"+ fasta_base +".log",
#    output: refdir + "/"+ fasta_base +".sa",
#    params:
#        out_prefix = refdir + "/" + fasta_base,
#        script = frag_scriptdir+ "/frag_index.sh",
#        log_dir = logdir + "/reference_index",
#    shell:
#        """
#        mkdir -p {params.log_dir}
#        {params.script} {input} {params.out_prefix} &> {log}
#        """

# Adapter-trim and QC reads with fastp
rule frag_fastp:
    benchmark: bench_dir + "/{library}_frag_fastp.benchmark.txt",
    input:
        read1 = frag_fastqs + "/{library}_raw_"+ PE1 + ".fastq.gz",
        read2 = frag_fastqs + "/{library}_raw_"+ PE2 + ".fastq.gz", 
    log:
        cmd = logdir + "/fastp/{library}_frag_fastp.log",
        html = logdir + "/fastp/{library}_frag_fastp.html",
        json = logdir + "/fastp/{library}_frag_fastp.json",
    output:
        read1 = frag_fastqs + "/{library}_processed_"+ PE1 + ".fastq.gz",
        read2 = frag_fastqs + "/{library}_processed_"+ PE2 + ".fastq.gz",
        failed = frag_fastqs + "/{library}_failed_fastp.fastq.gz",
        unpaired1 = frag_fastqs + "/{library}_unpaired_"+ PE1 + ".fastq.gz",
        unpaired2 = frag_fastqs + "/{library}_unpaired_"+ PE2 + ".fastq.gz",
    params:
        script = frag_scriptdir+ "/fastp.sh",
        threads = getthreads("frag_fastp"),
        tmp = config.get("tmp", "/tmp"),
        log_dir = logdir + "/fastp",
    resources:
        tmpdir = config.get("tmp", "/tmp"),
    conda:
        ENV_FRAG_READS
    shell:
        """
        mkdir -p {params.log_dir}
        fastp --detect_adapter_for_pe \
            --disable_length_filtering \
            --failed_out {output.failed} \
            --in1 {input.read1} \
            --in2 {input.read2} \
            --html {log.html} \
            --json {log.json} \
            --out1 {output.read1} \
            --out2 {output.read2} \
            --unpaired1 {output.unpaired1} \
            --unpaired2 {output.unpaired2} \
            --thread {params.threads}\
        > {log.cmd} 2>&1
        """
# Align reads with BWA
rule frag_align:
    benchmark: bench_dir + "/{library}_frag_align.benchmark.txt",
    input:
        ref = genome_fasta,
        read1 = frag_fastqs + "/{library}_processed_"+ PE1 + ".fastq.gz",
        read2 = frag_fastqs + "/{library}_processed_"+ PE2 + ".fastq.gz",
    log: logdir + "/align/{library}_frag_align.log",
    output:
        sort = frag_bams + "/{library}_raw.bam",
        index = frag_bams + "/{library}_raw.bam.bai",
    params:
        script = frag_scriptdir+ "/align.sh",
        threads = getthreads("frag_align"),
        log_dir = logdir + "/align",
    conda:
        ENV_FRAG_READS
    shell:
        """
        mkdir -p {params.log_dir}
        {params.script} \
        {input.ref} \
        {input.read1} \
        {input.read2} \
        {params.threads} \
        {output.sort} &> {log}
        """
# Remove PCR duplicates from aligned reads
rule frag_dedup:
    benchmark: bench_dir + "/{library}_frag_dedup.benchmark.txt",
    input: frag_bams + "/{library}_raw.bam",
    log: logdir + "/dedup/{library}_frag_dedup.log",
    output: frag_bams + "/{library}_dedup.bam",
    params:
        script = frag_scriptdir + "/dedup.sh",
        threads = getthreads("frag_dedup"),
        log_dir = logdir + "/dedup",
    conda:
        ENV_FRAG_READS
    shell:
        """
        mkdir -p {params.log_dir}
        {params.script} \
        {input} \
        {output} \
        {params.threads} &> {log}
        """

# Filter de-duplicated alignments.
# Remove unmapped, not primary, and duplicate reads. Additional location filter by config bedfile variable.
checkpoint frag_filter_alignment:
    benchmark: bench_dir + "/{library}_frag_filter_alignment.benchmark.txt",
    input: frag_bams + "/{library}_dedup.bam",
    log: logdir + "/filter/{library}_frag_filter_alignment.log",
    output: frag_bams + "/{library}_filt.bam",
    params:
        script = frag_scriptdir+ "/filter_alignment.sh",
        threads = getthreads("frag_filter_alignment"),
        log_dir = logdir + "/filter",
    conda:
        ENV_FRAG_READS
    shell:
        """
        mkdir -p {params.log_dir}
        {params.script} \
        {input} \
        {params.threads} \
        {output} &> {log}
        """
# Get read quality by FASTQC
rule frag_fastqc:
    benchmark: bench_dir + "/{library}_{processing}_{read}_frag_fastqc.benchmark.txt",
    input: frag_fastqs + "/{library}_{processing}_{read}.fastq.gz",
    log: logdir + "/fastqc/{library}_{processing}_{read}_frag_fastqc.log",
    output:
        qcdir + "/{library}_{processing}_{read}_fastqc.html",
        qcdir + "/{library}_{processing}_{read}_fastqc.zip",
    params:
        outdir = qcdir,
        script = frag_scriptdir+ "/fastqc.sh",
        threads = getthreads("frag_fastqc"),
        log_dir = logdir + "/fastqc"
    conda:
        ENV_FRAG_READS
    shell:
        """
        mkdir -p {params.log_dir}
        {params.script} \
        {input} \
        {params.outdir} \
        {params.threads} &> {log}
        """
# Get alignment QC using samtools
rule frag_alignment_qc:
    input: frag_bams + "/{library}_{processing}.bam",
    log:
        flagstat = logdir + "/alignment_qc/{library}_{processing}_flagstat.log",
        samstat = logdir + "/alignment_qc/{library}_{processing}_samstats.log",
    output:
        flagstat = qcdir + "/{library}_{processing}_flagstat.txt",
        samstat = qcdir + "/{library}_{processing}_samstats.txt",
    params:
        script = frag_scriptdir+ "/alignment_qc.sh",
        threads = getthreads("frag_alignment_qc"),
        log_dir = logdir + "/alignment_qc",
    conda:
        ENV_FRAG_READS
    shell:
        """
        mkdir -p {params.log_dir}
        {params.script} \
        {input} \
        {log.flagstat} \
        {log.samstat} \
        {output.flagstat} \
        {output.samstat} \
        {params.threads} &> {log}
        """
# Sequencing depth metrics via Picard
rule frag_picard_depth:
    benchmark: bench_dir + "/{library}_frag_picard_depth.benchmark.txt",
    input: frag_bams + "/{library}_filt.bam",
    log: logdir + "/picard_depth/{library}_picard_depth.log",
    output: qcdir + "/{library}_picard_depth.txt",
    params:
        script = frag_scriptdir+ "/picard_depth.sh",
        threads =  getthreads("frag_picard_depth"),
        log_dir = logdir + "/picard_depth",
    conda:
        ENV_FRAG_READS
    shell:
        """
        mkdir -p {params.log_dir}
        {params.script} \
        {input} \
        {config[genome_fasta]} \
        {output} \
        {params.threads} &> {log}
        """
# Get fragment sizes using deepTools
rule frag_bampefragsize:
    benchmark: bench_dir + "/frag_bampefragsize.benchmark.txt",
    input: expand(frag_bams + "/{library}_filt.bam", library = FRAG_LIBS),
    log: logdir + "/bampefragsize/bampefragsize.log",
    output:
        raw = qcdir + "/deeptools_frag_lengths.txt",
        hist = qcdir + "/deeptools_frag_lengths.png",
    params:
        blacklist = config["blklist"],
        script = frag_scriptdir+ "/bampefragsize.sh",
        threads = getthreads("frag_bampefragsize"),
        log_dir = logdir + "/bampefragsize",
    conda:
        ENV_FRAG_READS
    shell:
        """
        mkdir -p {params.log_dir}
        {params.script} \
        "{input}" \
        {log} \
        {output.hist} \
        {output.raw} \
        {params.blacklist} \
        {params.threads} &> {log}
        """
# Make deeptools bamCoverage bedfile
rule frag_bamcoverage:
    benchmark: bench_dir + "/{library}_frag_bamcoverage.benchmark.txt",
    input: frag_bams + "/{library}_filt.bam",
    log: logdir + "/bamcoverage/{library}_bamcoverage.log",
    output: qcdir + "/{library}_bamcoverage.bg",
    params:
        bin = "10000",
        blacklist = config["blklist"],
        script = frag_scriptdir+ "/bamcoverage.sh",
        threads = getthreads("frag_bamcoverage"),
        read_type = config.get("seq_params", {}).get("read_type", "PE"),
        read_length = config.get("seq_params", {}).get("read_length", 150),
        log_dir = logdir + "/bamcoverage",
    conda:
        ENV_FRAG_READS
    shell:
        """
        mkdir -p {params.log_dir}
        {params.script} \
        {input} \
        {output} \
        {params.bin} \
        {params.blacklist} \
        {params.threads} \
        {params.read_type} \
        {params.read_length} &> {log}
        """
# Make deepTools plotCoverage coverage maps for all filtered bams
rule frag_plotcoverage:
    benchmark: bench_dir + "/{library}_frag_plotcoverage.benchmark.txt",
    input: 
        frag_bams + "/{library}_filt.bam"
    log: logdir + "/plotcoverage/{library}_plotcoverage.log",
    output:
        raw = qcdir + "/{library}_coverage.tsv",
        plot = qcdir + "/{library}_coverage.pdf",
    params:
        blacklist = config["blklist"],
        script = frag_scriptdir+ "/plotcoverage.sh",
        threads = getthreads("frag_plotcoverage"),
        log_dir = logdir + "/plotcoverage",
    conda:
        ENV_FRAG_READS
    shell:
        """
        mkdir -p {params.log_dir}
        {params.script} \
        "{input}" \
        {params.blacklist} \
        {params.threads} \
        {output.raw} \
        {output.plot} &> {log}
        """
# Aggregate QC files using MultiQC
rule frag_multiqc:
    benchmark: bench_dir + "/frag_multiqc.benchmark.txt"
    input:
        fastp = expand(logdir + "/fastp/{library}_frag_fastp.json", library=FRAG_LIBS),  
        fastqc = expand(qcdir + "/{library}_{processing}_{read}_fastqc.zip", 
                       library=FRAG_LIBS, 
                       processing=["raw", "processed", "unpaired"], 
                       read=[PE1,PE2]),
        samstats = expand(qcdir + "/{library}_{processing}_samstats.txt", 
                         library=FRAG_LIBS, 
                         processing=["raw","filt"]),
        flagstat = expand(qcdir + "/{library}_{processing}_flagstat.txt", 
                         library=FRAG_LIBS, 
                         processing=["raw","filt"]),
        picard = expand(qcdir + "/{library}_picard_depth.txt", library=FRAG_LIBS),
        frag_len = qcdir + "/deeptools_frag_lengths.txt",
        coverage = expand(qcdir + "/{library}_coverage.tsv", library=FRAG_LIBS)
    output:
        html = qcdir + "/frag_multiqc.html",
        fastqc = qcdir + "/frag_multiqc_data/multiqc_fastqc.txt",
        samstats = qcdir + "/frag_multiqc_data/multiqc_samtools_stats.txt",
        picard = qcdir + "/frag_multiqc_data/multiqc_picard_wgsmetrics.txt",
        flagstat = qcdir + "/frag_multiqc_data/multiqc_samtools_flagstat.txt"
    log: logdir + "/multiqc/multiqc.log"
    params:
        script = frag_scriptdir + "/multiqc.sh",
        input_files = lambda w, input: " ".join([
            *input.fastp,
            *input.fastqc,
            *input.samstats,
            *input.flagstat,
            *input.picard,
            input.frag_len,
            *input.coverage
        ]),
        log_dir = logdir + "/multiqc",
    conda:
        ENV_FRAG_READS
    shell:
        """
        mkdir -p {qcdir}/frag_multiqc_data
        mkdir -p {params.log_dir}
        
        {params.script} \
            '{params.input_files}' \
            frag_multiqc \
            {qcdir} \
            &> {log}
        """
# Make a tab-separated aggregate QC table
checkpoint frag_make_qc_tsv:
    benchmark: bench_dir + "/frag_make_qc_tsv.benchmark.txt",
    input:
        fq = qcdir + "/frag_multiqc_data/multiqc_fastqc.txt",
        mqsam = qcdir + "/frag_multiqc_data/multiqc_samtools_stats.txt",
        mqflag = qcdir + "/frag_multiqc_data/multiqc_samtools_flagstat.txt",
        picard = qcdir + "/frag_multiqc_data/multiqc_picard_wgsmetrics.txt",
        deeptools_frag = qcdir + "/deeptools_frag_lengths.txt",
        deeptools_cov = expand(qcdir + "/{library}_coverage.tsv", library = FRAG_LIBS)
    log: logdir + "/make_qc_tsv/make_qc_tsv.log",
    output:
        readqc = qcdir + "/summary/frag_read_qc.tsv",
        fraglen = qcdir + "/summary/frag_frag_len.tsv",
        deeptools_cov = qcdir +"/frag_coverage.tsv",
    params:
        script = frag_scriptdir+ "/make_qc_tsv.R",
        script_cov = frag_scriptdir + "/combine_deeptool_coverage.R",
        log_dir = logdir + "/make_qc_tsv",
    conda:
        ENV_FRAG_READS
    shell:
        """
        mkdir -p {params.log_dir}
        mkdir -p {qcdir}/summary

        Rscript {params.script_cov} --out {output.deeptools_cov} {input.deeptools_cov}

        Rscript {params.script} \
        {input.fq} \
        {input.mqsam} \
        {input.mqflag} \
        {input.picard} \
        {input.deeptools_frag} \
        {output.deeptools_cov} \
        {output.readqc} \
        {output.fraglen} >& {log}
        """
