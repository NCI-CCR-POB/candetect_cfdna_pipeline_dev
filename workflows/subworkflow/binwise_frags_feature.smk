############################
###   binwise_frags features Snakefile   ###
############################

#########1#########2#########3#########4#########5#########6#########7#########8
#
# A snakefile for binwise fragment feature profiling

rule filt_bam_to_frag_bed:
    benchmark: benchdir + "/{library}_filt_bam_to_frag_bed.benchmark.txt",
    input: frag_bams + "/{library}_filt.bam",
    log: logdir_binwise + "/bam_to_bed/{library}_filt_bam_to_frag_bed.log",
    output: bin_wise_beds + "/{library}_filt.bed",
    conda:
        ENV_FRAGENTOMICS
    params:
        fragdir = frag,
        bin_wise_fragdir = bin_wise_frag,
        out_dir = bin_wise_beds,
        fasta = genome_fasta,
        script = frag_scriptdir + "/filt_bam_to_frag_bed.sh",
        threads = getthreads("filt_bam_to_frag_bed"),
        log_dir = logdir_binwise + "/bam_to_bed",
    shell:
        """
        mkdir -p {params.log_dir} 
        mkdir -p {params.fragdir} 
        mkdir -p {params.bin_wise_fragdir} 
        mkdir -p {params.out_dir} 

        {params.script} \
        {input} \
        {params.fasta} \
        {params.threads} \
        {output}
        """

# Make GC distributions
rule gc_distro:
    benchmark: benchdir + "/{library}_frag_gc_distro.benchmark.txt",
    input: bin_wise_beds + "/{library}_filt.bed",
    log: logdir_binwise + "/gc_distro/{library}_frag_gc_distro.log",
    output: bin_wise_gc_distros + "/{library}_gc_distro.csv",
    conda:
        ENV_FRAGENTOMICS
    params:
        script = frag_scriptdir + "/gc_distro.R",
        out_dir = bin_wise_gc_distros,
        log_dir = logdir_binwise + "/gc_distro",
    shell:
        """
        mkdir -p {params.log_dir}
        mkdir -p {params.out_dir}

        Rscript {params.script} \
        {input} \
        {output} > {log} 2>&1
        """

# Make healthy GC distributions summary file
rule healthy_gc:
    benchmark: benchdir + "/frag_healthy_gc.benchmark.txt",
    input: expand(bin_wise_gc_distros + "/{library}_gc_distro.csv", library = FRAG_HEALTHY_LIBRARIES),
    log: logdir_binwise + "/healthy_gc/frag_healthy_gc.log",
    output: bin_wise_gc_distros + "/healthy_med.rds",
    conda:
        ENV_FRAGENTOMICS
    params:
        distro_dir = bin_wise_gc_distros,
        log_dir = logdir_binwise + "/healthy_gc",
        script = frag_scriptdir + "/healthy_gc.R",
        existing_rds = healthy_med_rds
    shell:
        """
        mkdir -p {params.log_dir} 
        if [ -n "{params.existing_rds}" ] && [ -s "{params.existing_rds}" ]; then
            echo "Using existing healthy GC RDS: {params.existing_rds}" > {log}
            cp "{params.existing_rds}" "{output}"
        else
            echo "No existing healthy GC RDS found. Running healthy_gc.R..." > {log}
            Rscript {params.script} \
                {params.distro_dir} \
                "{input}" \
                {output} > {log} 2>&1
        fi

        """

# Sample fragments by healthy GC proportions
rule frag_gc_sample:
    benchmark: benchdir + "/{library}_frag_gc_sample.benchmark.txt",
    input:
        frag_bed = bin_wise_beds+ "/{library}_filt.bed",
        healthy_med = bin_wise_gc_distros + "/healthy_med.rds",
    log: logdir_binwise + "/gc_sample/{library}_frag_gc_sample.log",
    output: bin_wise_beds+ "/{library}_sampled_frag.bed",
    conda:
        ENV_FRAGENTOMICS
    params:
        log_dir = logdir_binwise + "/gc_sample",
        out_dir = bin_wise_beds,
        script = frag_scriptdir + "/gc_sample.R",
    shell:
        """
        mkdir -p {params.log_dir} 
        mkdir -p {params.out_dir} 

        Rscript {params.script} \
        {input.healthy_med} \
        {input.frag_bed} \
        {output} > {log} 2>&1
        """
# Sum fragments in short and long length groups
rule frag_sum:
    benchmark: benchdir + "/{library}_frag_sum.benchmark.txt",
    input:
        bin_wise_beds+ "/{library}_sampled_frag.bed",
    log:
        logdir_binwise + "/frag_sum/{library}_bin_wise_frag_window_sum.log",
    output:
        short = bin_wise_beds+ "/{library}_norm_short.bed",
        long =  bin_wise_beds+ "/{library}_norm_long.bed",
    conda:
        ENV_FRAGENTOMICS
    params:
        log_dir = logdir_binwise + "/frag_sum",
        script = frag_scriptdir + "/frag_window_sum.sh",
        threads = getthreads("frag_sum"),
    shell:
        """
        mkdir -p {params.log_dir}
        {params.script} \
        {input} \
        {output.short} {output.long} &> {log}
        """

# Count short and long fragments intersecting kept genomic windows
rule frag_window_count:
    benchmark: benchdir + "/{library}_bin_wise_frag_window_int.benchmark.txt",
    input:
        short = bin_wise_beds+ "/{library}_norm_short.bed",
        long = bin_wise_beds+ "/{library}_norm_long.bed",
        matbed = refdir + "/keep_5mb.bed",
        genome = genome_fasta + ".fai",
    log: logdir_binwise + "/frag_window_count/{library}_bin_wise_frag_window_int.log",
    output:
        short = bin_wise_counts + "/{library}_cnt_short.tmp",
        long = bin_wise_counts + "/{library}_cnt_long.tmp",
    conda:
        ENV_FRAGENTOMICS
    params:
        script = frag_scriptdir + "/frag_window_int.sh",
        threads = getthreads("frag_window_count"),
        out_dir = bin_wise_counts,
        log_dir = logdir_binwise + "/frag_window_count",
    shell:
        """
        # Create output directories
        mkdir -p {params.out_dir} 
        mkdir -p {params.log_dir}

        # Log start of processing
        echo "Starting processing for library {wildcards.library}" >> {log} 2>&1
        echo "Input files:" >> {log} 2>&1
        echo "Short fragments: {input.short}" >> {log} 2>&1
        echo "Long fragments: {input.long}" >> {log} 2>&1
        echo "Matrix bed: {input.matbed}" >> {log} 2>&1
        echo "Genome file: {input.genome}" >> {log} 2>&1

        # Process short fragments
        echo "Processing short fragments..." >> {log} 2>&1
        {params.script} \
            {input.short} \
            {input.matbed} \
            {output.short} \
            {input.genome} >> {log} 2>&1 || {{
            echo "ERROR: Failed to process short fragments for {wildcards.library}" >> {log} 2>&1
            echo "Check the log file for details: {log}" >> {log} 2>&1
            exit 1
        }}

        # Process long fragments
        echo "Processing long fragments..." >> {log} 2>&1
        {params.script} \
            {input.long} \
            {input.matbed} \
            {output.long} \
            {input.genome} >> {log} 2>&1 || {{
            echo "ERROR: Failed to process long fragments for {wildcards.library}" >> {log} 2>&1
            echo "Check the log file for details: {log}" >> {log} 2>&1
            exit 1
        }}

        # Verify outputs
        echo "Verifying output files..." >> {log} 2>&1
        for out_file in {output.short} {output.long}; do
            if [[ ! -s $out_file ]]; then
                echo "ERROR: Output file is empty: $out_file" >> {log} 2>&1
                exit 1
            fi
            echo "Successfully created: $out_file" >> {log} 2>&1
            echo "First few lines of output:" >> {log} 2>&1
            head -n 3 $out_file >> {log} 2>&1
        done

        echo "Successfully completed processing for {wildcards.library}" >> {log} 2>&1
        """

# Merge short and long fragment counts by genomic poistion for all libraries
rule frag_count_merge:
    benchmark:
        benchdir + "/frag_count_merge.benchmark.txt",
    input:
        expand(bin_wise_counts + "/{library}_cnt_{length}.tmp",
            library = FRAG_LIBS,
            length = ["short","long"],
        ),
    log:
        logdir_binwise + "/frag_count_merge/frag_count_merge.log",
    output:
        bin_wise_frag + "/frag_counts.tsv",
    conda:
        ENV_FRAGENTOMICS
    params:
        counts_dir = bin_wise_counts,
        script = frag_scriptdir + "/count_merge.sh",
        log_dir = logdir_binwise + "/frag_count_merge",
    shell:
        """
        mkdir -p {params.log_dir}
        {params.script} \
        {params.counts_dir} \
        {output} \
        &> {log}
        """
rule unit_cent_sd:
    benchmark: benchdir + "/unit_cent_sd.benchmark.txt",
    input: bin_wise_frag + "/frag_counts.tsv",
    log: logdir_binwise + "/unit_cent_sd/unit_cent_sd.log",
    output: bin_wise_frag + "/ratios.tsv",
    conda:
        ENV_FRAGENTOMICS
    params:
        script = frag_scriptdir + "/make_ratios.R",
        log_dir = logdir_binwise + "/unit_cent_sd",
    shell:
        """
        mkdir -p {params.log_dir}
        Rscript {params.script} \
        {input} {output} > {log} 2>&1
        """
        
rule binwise_armz_raw:
    input:
        bin_wise_frag + "/frag_counts.tsv",
    output:
        bin_wise_frag + "/binwise_armz.raw.csv",
    conda:
        ENV_FRAGENTOMICS
    log:
        logdir_binwise + "/armz_raw/binwise_armz_raw.log",
    params:
        log_dir = logdir_binwise + "/armz_raw",
        script = frag_scriptdir + "/binwise_armz.R",
        cytoband = config["cytoband"],
        sample_list = samples,
        control_term = normal_control,
        existing_healthy = config.get("existing_binwise_healthy_mean_sd", ""),
    shell:
        """
        mkdir -p {params.log_dir}
        if [ -n "{params.existing_healthy}" ] && [ -s "{params.existing_healthy}" ]; then
            echo "Using existing healthy mean and sd: {params.existing_healthy}" > {log}
            cp "{params.existing_healthy}" "{bin_wise_frag}"
            Rscript {params.script} \
                {params.cytoband} \
                {params.sample_list} \
                {params.control_term} \
                {input} \
                {output} \
                {params.existing_healthy} \
                &> {log}
        else
            echo "No existing healthy mean and sd in yaml, using the one generated in this pipeline." > {log}
            Rscript {params.script} \
                {params.cytoband} \
                {params.sample_list} \
                {params.control_term} \
                {input} \
                {output} \
                &> {log}
        fi
        """
rule binwise_features:
    input:
        arm_z = bin_wise_frag +  "/binwise_armz.raw.csv",
        ratio = bin_wise_frag + "/ratios.tsv",
    output:
        arm_z = bin_wise_frag + "/binwise_armz.clean.csv",
        ratio = bin_wise_frag + "/binwise_ratio.norm_rows.clean.csv",
    conda:
        ENV_FRAGENTOMICS
    log:
        logdir_binwise + "/features/binwise_features.log",
    params:
        script  = frag_scriptdir + "/binwise_features.py",
        out_dir = bin_wise_frag,
        log_dir = logdir_binwise + "/features",
    shell:
        """
        mkdir -p {params.log_dir}

        python3 {params.script} \
        --armz_csv {input.arm_z} \
        --ratio_tsv {input.ratio} \
        --out_dir {params.out_dir} \
        &> {log}
        """
