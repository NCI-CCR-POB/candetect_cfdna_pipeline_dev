
# Introduction
This repository hosts a snakemake workflow for basic processing and fragmentomic analysis of whole-genome sequencing reads from cell-free DNA. These workflows, written by Jeff Szymanski, R. Taylor Sundby, and Alex C. Pan, were used for the following publications:

1.	Szymanski JJ, Sundby RT, Jones PA, Srihari D, Earland N, Harris PK, et al. Cell-free DNA ultra-low-pass whole genome sequencing to distinguish malignant peripheral nerve sheath tumor (MPNST) from its benign precursor lesion: A cross-sectional study. PLOS Medicine. 2021;18(8):e1003734-e. doi: 10.1371/JOURNAL.PMED.1003734.

2.	Sundby RT, Szymanski JJ, Pan A, Jones PA, Mahmood SZ, Reid OH, et al. Early detection of malignant and pre-malignant peripheral nerve tumors using cell-free DNA fragmentomics. Clin Cancer Res. 2024. Epub 20240802. doi: 10.1158/1078-0432.CCR-24-0797. PubMed PMID: 39093127.

# Use
-   See the example configuration yaml `config/default.yaml` and wrapper workflow `inputs/library_list_format.tsv` for necessary run conditions.
-   Workflow can be run individually (Default_launch_SMK.sh) or as an integrated workflow (Default_launch_master_workflow.sh) which will run pre-processing first and then, if completed without errors, all downstream analsysis (end motifs, cfDNA size NMF, bin-wise fragmentomics) workflows in parallel.  If running individually, frag_reads.smk (pre-processing) MUST be run first.
-   In its current itteration, software is loaded from the host HPC is there is not a saved image. Snakemakes and configs can be changed to use a docker container (see `docker` and <https://hub.docker.com/repository/docker/jeszyman/cfdna_wgs/general> for software versions used in publicatons above)
-   Several reference files are required and must be downloaded seperately:
        * Reference genome e.g. GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
        * Example cytoBand.txt, gc5mb.bed, keep_5mb.bed. hg38_chr_arms.tsv, hg38-blasklist.v2.bed, and hg39.chrom.sizes are included in `ref/` and are generated to hg38 using the healthy samples from the publications above. These may need to be re-generated for project specific needs.  
