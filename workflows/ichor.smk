###################################################################################################
######rules for ichor cna analysis
###################################################################################################

rule size_select:
    input:
        bam_file = data_dir + "/bams/{samples}_filt.bam"
    output:
        ss_bam = size_selected_dir + "/{samples}.ss.bam",
    log:
        logdir_ichor + "/size_select/{samples}.log",
    conda:
        ENV_FRAGENTOMICS
    params:
        out_dir = size_selected_dir,
        logdir_ichor = logdir_ichor + "/size_select",
        min_size = config["min_fragment_size"],
        max_size = config["max_fragment_size"],
        threads = getthreads("size_select")
    shell:
        """
        mkdir -p {params.out_dir}
        mkdir -p {params.logdir_ichor}
        
        #module load samtools
        samtools view -h {input.bam_file} -@ {params.threads} | \
        awk 'substr($0,1,1)=="@" || ($9>= {params.min_size} && $9<={params.max_size}) || ($9<=-{params.min_size} && $9>=-{params.max_size})' | \
        samtools view -b -@ {params.threads} > {output.ss_bam}
        samtools index -@ {params.threads} {output.ss_bam}
        """

rule read_counter:
    input:
        bam_file = size_selected_dir + "/{samples}.ss.bam",
    output:
        wig_file = wig_dir + "/{samples}.wig",
    log:
        logdir_ichor + "/wig_files/{samples}.log",
    conda:
        ENV_FRAGENTOMICS
    params:
        window		= config["window"],
        quality		= config["quality"],
        chrs		= config["chromosomes"],
        out_dir = wig_dir,
        logdir_ichor = logdir_ichor + "/wig_files",
    shell:
        """
        mkdir -p {params.out_dir}
        mkdir -p {params.logdir_ichor}

        readCounter \
        -w {params.window} \
        -q {params.quality} \
        -c {params.chrs} \
        {input.bam_file} \
        > {output.wig_file} \
        2> {log}
        """
        
rule ichorCNA:
    input:
        tum_wig = wig_dir + "/{tumor}.wig",
    output:
        cna = ichor_dir + "/{tumor}/{tumor}.cna.seg",
        #corrDepth="results/ichorCNA/{tumor}/{tumor}.correctedDepth.txt",
        #param="results/ichorCNA/{tumor}/{tumor}.params.txt",
        #segTxt="results/ichorCNA/{tumor}/{tumor}.seg.txt",
        #seg="results/ichorCNA/{tumor}/{tumor}.seg",
        #rdata="results/ichorCNA/{tumor}/{tumor}.RData"
    log:
        logdir_ichor + "/ichorCNA/{tumor}.log",
    conda:
        ENV_FRAGENTOMICS
    params:
        rscript           = config["ichorCNA_rscript"],
        libdir            = config["ichorCNA_libdir"], # added
        id                = "{tumor}",
        gcwig             = config["ichorCNA_gcWig"],
        mapwig            = config["ichorCNA_mapWig"],
        centromere        = config["ichorCNA_centromere"],
        normalpanel       = config["ichorCNA_normalPanel"],
        genomeBuild       = config["ichorCNA_genomeBuild"], # added
        genomeStyle       = config["ichorCNA_genomeStyle"], # added
        chrs              = config["ichorCNA_chrs"],
        chrTrain          = config["ichorCNA_chrTrain"],
        estimateNormal    = config["ichorCNA_estimateNormal"],
        estimatePloidy    = config["ichorCNA_estimatePloidy"],
        estimateClonality = config["ichorCNA_estimateClonality"],
        normal            = config["ichorCNA_normal"],
        ploidy            = config["ichorCNA_ploidy"],
        maxCN             = config["ichorCNA_maxCN"],
        scStates          = config["ichorCNA_scStates"],
        txnE              = config["ichorCNA_txnE"],
        txnStrength       = config["ichorCNA_txnStrength"],
        out_dir           = ichor_dir + "/{tumor}",
        logdir_ichor           = logdir_ichor + "/ichorCNA",

        #### UNUSED OPTIONS
        # includeHOMD=config["ichorCNA_includeHOMD"],
        # fracReadsChrYMale=config["ichorCNA_fracReadsInChrYForMale"],
        # minMapScore=config["ichorCNA_minMapScore"],
        # maxFracGenomeSubclone=config["ichorCNA_maxFracGenomeSubclone"],
        # maxFracCNASubclone=config["ichorCNA_maxFracCNASubclone"],
        # exons=config["ichorCNA_exons"],
        # plotFileType=config["ichorCNA_plotFileType"],
        # plotYlim=config["ichorCNA_plotYlim"],
        # 

    shell:
        """
        mkdir -p {params.out_dir}
        mkdir -p {params.logdir_ichor}
        
        Rscript {params.rscript} \
        --libdir {params.libdir} \
        --WIG {input.tum_wig} \
        --id {params.id} \
        --gcWig {params.gcwig} \
        --mapWig {params.mapwig} \
        --centromere {params.centromere} \
        --normalPanel {params.normalpanel} \
        --chrs \"{params.chrs}\" \
        --chrTrain \"{params.chrTrain}\" \
        --estimateNormal {params.estimateNormal} \
        --estimatePloidy {params.estimatePloidy} \
        --estimateScPrevalence {params.estimateClonality} \
        --normal \"{params.normal}\" \
        --ploidy \"{params.ploidy}\" \
        --maxCN {params.maxCN} \
        --scStates \"{params.scStates}\" \
        --txnE {params.txnE} \
        --txnStrength {params.txnStrength} \
        --outDir {params.out_dir} \
        --genomeStyle {params.genomeStyle} \
        --genomeBuild {params.genomeBuild} \
        &> {log}
        """

        ### UNUSED OPTIONS:
        # """
        # --includeHOMD {params.includeHOMD} \
        # --exons.bed {params.exons} \
        # --minMapScore {params.minMapScore} \
        # --fracReadsInChrYForMale {params.fracReadsChrYMale} \
        # --maxFracGenomeSubclone {params.maxFracGenomeSubclone} \
        # --maxFracCNASubclone {params.maxFracCNASubclone} \
        # --plotFileType {params.plotFileType} \
        # --plotYLim \"{params.plotYlim}\" \
        # """


rule features: 
    input:
        ichor_dir = directory(ichor_dir),
        cna_files = expand(ichor_dir + "/{samples}/{samples}.cna.seg", samples=FRAG_LIBS)
    output:
        tf_csv = ichor_dir + "/tf.csv"
    log:
        logdir_ichor + "/features/features.log"
    conda:
        ENV_FRAGENTOMICS
    params:
        script = "scripts/features.py",
        logdir_ichor = logdir_ichor + "/features"
    shell:
        """
        # Create log directory
        mkdir -p {params.logdir_ichor}

        # Check if Python script exists
        if [ ! -f {params.script} ]; then
            echo "ERROR: Python script not found at {params.script}" >> {log} 2>&1
            exit 1
        fi

        # Check input directory exists and has content
        if [ ! -d {input.ichor_dir} ]; then
            echo "ERROR: Input directory {input.ichor_dir} does not exist" >> {log} 2>&1
            exit 1
        fi

        # Create output directory if it doesn't exist
        mkdir -p $(dirname {output.tf_csv})

        echo "Starting feature extraction..." >> {log} 2>&1
        echo "Input directory: {input.ichor_dir}" >> {log} 2>&1
        echo "Output file: {output.tf_csv}" >> {log} 2>&1
        
        # Load Python module if needed
        #module load python/3.9 || echo "Warning: Failed to load Python module" >> {log} 2>&1

        # Run the Python script with error handling
        python {params.script} \
            --ichor_dir {input.ichor_dir} \
            --output {output.tf_csv} \
            2>> {log} || {{
            echo "ERROR: Feature extraction failed" >> {log} 2>&1
            exit 1
        }}

        # Verify output was created
        if [ ! -f {output.tf_csv} ]; then
            echo "ERROR: Output file was not created" >> {log} 2>&1
            exit 1
        fi

        echo "Feature extraction completed successfully" >> {log} 2>&1
        """
