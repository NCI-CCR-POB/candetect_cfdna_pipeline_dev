# CanDetect cfDNA fragmentomics

This repository contains Snakemake workflows for paired-end cell-free DNA
whole-genome sequencing preprocessing, fragmentomics feature extraction,
modeling, prediction, and reporting.

The current workflow produces four feature streams:

- **Preprocessing and QC**: fastp, BWA alignment, duplicate removal, BAM
  filtering, FastQC, Picard/samtools QC, bamPEFragmentSize, deepTools coverage,
  and MultiQC.
- **End-motif profiles**: strand-aware 5-prime 4-mer extraction from filtered
  BAMs, aggregation to motif count and fraction matrices, motif-diversity
  scoring, and F-profile modeling.
- **Bin-wise fragmentomics**: filtered BAM to fragment BED conversion, GC
  normalization with healthy controls, short and long fragment counts across
  kept 5 Mb windows, chromosome-arm z scores, ratio features, and logistic
  modeling.
- **Fragment-size NMF**: Picard insert-size histograms, aggregate
  30-700 bp histograms, NMF component models
- **CNA and tumor fraction with ichorCNA**:  90-150 bp size-selected ichorCNA,
  and tumor-fraction output.

The final report is written to:

```text
<data_dir>/analysis/cfdna_wgs/report/cfdna_wgs_report.html
```

## Key updates

- Added the integrated modeling workflow `workflows/fragmentomics_model.smk`
  for feature generation, model training or validation, AUC generation,
  ichorCNA, combined tables, voting plots, and HTML reporting.
- Added the integrated prediction workflow `workflows/fragmentomics_pred.smk`
  for applying existing bin-wise, NMF, and F-profile models plus training-set
  cutpoints to new samples.
- Added Mode 3 cutpoint application through `scripts/apply_cutpoints.py`; this
  mode applies existing cutpoints and writes predicted class calls without
  recomputing ROC curves or AUC-based thresholds.
- Updated end-motif extraction to call `scripts/extract_5p_motifs.sh` for
  strand-aware 5-prime cfDNA 4-mers and `scripts/aggregate_motif_counts.py`
  for motif count and fraction matrices with motif rows by library columns.
- Added support for reusing training-set healthy-control quantities during
  prediction with `existing_health_gc_rds` and
  `existing_binwise_healthy_mean_sd`.
- Added prediction-specific combined-output and report rules in
  `workflows/combine_datasets_pred.smk`, `scripts/create_voting_heatmap_pred.R`,
  and `scripts/create_html_report_pred.py`.

## Repository layout

```text
config/
  default.yaml              # copy and edit for each run
  cluster.yaml              # per-rule resources used inside Snakemake files
envs/
  frag_reads.yaml           # conda environment for preprocessing rules
  fragmentomics_py3.9.15.yaml
inputs/
  library_list_format.tsv   # sample manifest template
profile/
  config.yaml               # Snakemake cluster profile
ref/
  cytoBand.txt, keep_5mb.bed, gc5mb.bed, hg38-blacklist.v2.bed, ...
scripts/
  shell, Python, and R scripts called by the workflow rules
workflows/
  fragmentomics_model.smk   # train or validate models and make a full report
  fragmentomics_pred.smk    # apply existing models and cutpoints
  frag_reads.smk            # preprocessing-only workflow
  subworkflow/*.smk         # feature, model, and prediction rule groups
```

## Inputs

### FASTQs

The preprocessing workflow expects paired-end FASTQs. Put the R1 FASTQ path for
each library in the manifest `file` column; `scripts/detect_read_pairs.sh`
locates the matching R2 file and creates standardized symlinks in:

```text
<data_dir>/fastq/
```

The read-pair tokens are controlled by `PE1` and `PE2` in the run YAML:

```yaml
PE1: "R1"
PE2: "R2"
```

### Sample manifest

Copy `inputs/library_list_format.tsv` and provide at least these columns:

| column | description |
| --- | --- |
| `library` | Unique library ID. This becomes the sample prefix for BAMs, feature files, and scores. |
| `cohort` | Cohort label used for pairwise model comparisons. Values must match `cohorts` in the YAML exactly for training set. |
| `subject` | Subject or patient ID used by longitudinal plotting. |
| `file` | Path to the input R1 FASTQ. |
| `library_type` | Assay label, for example `wgs`. |
| `isolation_type` | Input molecule label, for example `cfdna`. |

Additional columns, such as `specimen`, may be present and are carried through
where downstream scripts use them.

### Reference files

Set the hg38 reference paths in the YAML. The current rules require:

- `genome_fasta`: reference FASTA used by BWA, GC extraction, and motif
  extraction.
- `ref_dir`: directory containing FASTA/BWA index files and `keep_5mb.bed`.
- `cytoband`: chromosome-arm annotations for bin-wise arm z scores.
- `blklist`: blacklist BED used by deepTools coverage rules.
- `gc5mb`: 5 Mb GC-content BED.
- ichorCNA support files:
  - `ichorCNA_gcWig`
  - `ichorCNA_mapWig`
  - `ichorCNA_centromere`
  - `ichorCNA_normalPanel`

Example hg38 helper files are included under `ref/`, but project-specific runs
should point to the vetted reference directory on the analysis system.

## Configure a run

Start from `config/default.yaml`:

```bash
cp config/default.yaml config/my_run.yaml
```

Update these required paths first:

```yaml
data_dir: "/path/to/output/run"
frag_repo: "/path/to/candetect_cfdna_pipeline_dev"
samples: "/path/to/input_libraries.tsv"

genome_fasta: "/path/to/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
ref_dir: "/path/to/ref"
cytoband: "/path/to/ref/cytoBand.txt"
blklist: "/path/to/ref/hg38-blacklist.v2.bed"
gc5mb: "/path/to/ref/gc5mb.bed"

env_frag_reads: "/path/to/candetect_cfdna_pipeline_dev/envs/frag_reads.yaml"
env_fragentomics: "/path/to/candetect_cfdna_pipeline_dev/envs/fragmentomics_py3.9.15.yaml"
tmp: "/path/to/output/run/tmp"
```

Then set the cohort labels. The workflow runs every pairwise comparison in the
order returned by `itertools.combinations(cohorts, 2)`.

```yaml
normal_control: "healthy"
cohorts: ["healthy", "PN", "AN", "MPNST"]
```

`normal_control` must match the healthy-control label in the manifest
`cohort` column because those libraries define the bin-wise GC summary and
chromosome-arm mean/SD used for standardization.

## Choose a model mode

The bin-wise, NMF, and F-profile branches can run in three modes.

### Mode 1: train new models

Use this mode for a training set. The workflow trains models, chooses
cutpoints, computes AUCs, and creates the HTML report.

```yaml
existing_nmf_model: "no"
existing_delfi_model: "no"
existing_motif_model: "no"

existing_binwise_scores_dir: "no"
existing_nmf_scores_dir: "no"
existing_motif_scores_dir: "no"
```

Run with `-m model` after setting these values.

### Mode 2: validate existing models and compute new AUCs

Use this mode when model files already exist but the new cohort still has
ground-truth labels and should receive a new validation AUC.

Set the model locations and leave the existing score directories as `"no"`:

```yaml
existing_nmf_model: "/path/to/training/analysis/cfdna_wgs/nmf/model"
existing_delfi_model: "/path/to/training/analysis/cfdna_wgs/bin_wise_fragmentomics/model/model"
existing_motif_model: "/path/to/training/analysis/cfdna_wgs/features/fprofile_model/fprofile_model.pkl"

existing_binwise_scores_dir: "no"
existing_nmf_scores_dir: "no"
existing_motif_scores_dir: "no"
```

Run with `-m model` after setting these values.

### Mode 3: predict with existing models and cutpoints

Use this mode for new samples where the training-set thresholds should be
frozen. The workflow applies:

- bin-wise model pickle files from `existing_delfi_model`
- NMF model pickle files from `existing_nmf_model`
- the F-profile pickle from `existing_motif_model`
- bin-wise and NMF score cutpoints from the original `scores/` directories
- the best NMF component count from
  `<existing_nmf_scores_dir>/aucs/best_med_aucs.csv`
- original F-profile AUC/threshold files from `existing_motif_scores_dir`

Example:

```yaml
existing_nmf_model: "/path/to/training/analysis/cfdna_wgs/nmf/model"
existing_delfi_model: "/path/to/training/analysis/cfdna_wgs/bin_wise_fragmentomics/model/model"
existing_motif_model: "/path/to/training/analysis/cfdna_wgs/features/fprofile_model/fprofile_model.pkl"

existing_binwise_scores_dir: "/path/to/training/analysis/cfdna_wgs/bin_wise_fragmentomics"
existing_nmf_scores_dir: "/path/to/training/analysis/cfdna_wgs/nmf"
existing_motif_scores_dir: "/path/to/training/analysis/cfdna_wgs/features"

existing_health_gc_rds: "/path/to/training/bin_wise_fragmentomics/distros/healthy_med.rds"
existing_binwise_healthy_mean_sd: "/path/to/training/bin_wise_fragmentomics/binwise_armz_healthy_mean_sd.csv"
```

When applying an NMF model with a fixed component count, set both NMF bounds to
that count:

```yaml
min_nmf_components: 18
max_nmf_components: 18
```

## Run the wrapper

Use `-m model` to dry-run the integrated modeling workflow:

```bash
./run_snakemake.py -m model -cf config/my_run.yaml
```

Add `--rerun` to execute the model run:

```bash
./run_snakemake.py -m model -cf config/my_run.yaml --rerun
```

Use `-m predict` to dry-run the integrated prediction workflow:

```bash
./run_snakemake.py -m predict -cf config/my_run.yaml
```

Add `--rerun` to execute the prediction run:

```bash
./run_snakemake.py -m predict -cf config/my_run.yaml --rerun
```

## Run with sbatch

Use `submit_model.sh` for a model-training or validation run:

```bash
sbatch submit_model.sh
```

Use `submit_predict.sh` for an existing-model prediction run:

```bash
sbatch submit_predict.sh
```

Before submitting, update the `--configfile` argument in the selected submit
script so it points to the run YAML you prepared:

```bash
--configfile config/my_run.yaml
```

`submit_model.sh` already calls `workflows/fragmentomics_model.smk`, and
`submit_predict.sh` already calls `workflows/fragmentomics_pred.smk`.

## Main outputs

All primary outputs are written under:

```text
<data_dir>/analysis/cfdna_wgs/
```

Important subdirectories and files include:

| path | contents |
| --- | --- |
| `bin_wise_fragmentomics/` | Bin-wise features, probabilities, scores, and plots. |
| `features/combined_motif_counts_matrix.tsv` | End-motif count matrix, with motif rows and library columns. |
| `features/combined_motif_fractions_matrix.tsv` | End-motif fraction matrix, with motif rows and library columns. |
| `features/fprofiles.tsv` | F-profile contribution matrix per library. |
| `nmf/` | Insert-size histograms, NMF models, NMF probabilities, best component counts, and NMF scores. |
| `ichorCNA/tf.csv` | ichorCNA tumor-fraction summary used by the combined report. |
| `combined_dataset_<cohort1>_v_<cohort2>.csv` | Per-comparison merged fragmentomics and ichorCNA metrics. |
| `plots/voting/` | Majority-vote heatmaps and majority-vote CSVs. |
| `report/cfdna_wgs_report.html` | Final interactive HTML report. |

Logs are written under:

```text
<data_dir>/logs/
```

## Standalone workflows

The integrated workflows above are the preferred entry points for full model
or prediction runs. Individual workflows are still available when debugging one
branch at a time:

```bash
snakemake -s workflows/frag_reads.smk --configfile config/my_run.yaml --profile profile/
snakemake -s workflows/end_motifs.smk --configfile config/my_run.yaml --profile profile/
snakemake -s workflows/binwise_frags.smk --configfile config/my_run.yaml --profile profile/
snakemake -s workflows/NMF_size.smk --configfile config/my_run.yaml --profile profile/
```

Run `workflows/frag_reads.smk` first because the end-motif, bin-wise, NMF, and
ichorCNA branches consume filtered BAMs from:

```text
<data_dir>/bams/<library>_filt.bam
```

## Citation

If you use this workflow, please cite:

1. Szymanski JJ, Sundby RT, Jones PA, Srihari D, Earland N, Harris PK, et al.
   Cell-free DNA ultra-low-pass whole genome sequencing to distinguish
   malignant peripheral nerve sheath tumor from its benign precursor lesion: a
   cross-sectional study. *PLOS Medicine*. 2021;18(8):e1003734.
   doi:10.1371/JOURNAL.PMED.1003734

2. Sundby RT, Szymanski JJ, Pan A, Jones PA, Mahmood SZ, Reid OH, et al. Early
   detection of malignant and pre-malignant peripheral nerve tumors using
   cell-free DNA fragmentomics. *Clinical Cancer Research*. 2024.
   doi:10.1158/1078-0432.CCR-24-0797
