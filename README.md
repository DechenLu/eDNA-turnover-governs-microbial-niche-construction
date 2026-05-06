# Metagenome analysis pipeline: assembly, binning, refinement, annotation, statistical analysis, and community metabolic modeling

This repository contains the bioinformatics, ecological/statistical, and metabolic modeling workflow used for metagenome analysis in our study. The pipeline processes paired-end Illumina reads through quality control, assembly, binning, bin refinement, genome quality assessment, taxonomic assignment, functional annotation, downstream ecological/statistical analysis, and genome-scale metabolic modeling.

## Repository structure

```text
.
├── README.md
├── Rawdata/                  # Raw paired-end sequencing reads; not included in the repository
├── scripts/                  # Shell, Python, and R scripts used in the analysis
├── Genome/                   # Refined genome bins or MAGs used for downstream analysis
├── annotation/               # Prokka, RAST, GTDB-Tk, eggNOG, and KEGG annotation outputs
├── statistics/               # R scripts and statistical analysis outputs
├── metabolic_modeling/       # CarveMe and COBRApy scripts or model files
└── results/                  # Processed results and summary tables
```

Large raw sequencing files, intermediate assembly files, and external databases are not included in this repository. Users should modify input paths, database paths, thread numbers, and memory settings according to their own computing environment.

### Online websites and web-based analysis platforms

The following online platforms were used for genome annotation, genome comparison, differential abundance analysis, and visualization:

- RAST, for genome annotation
- ChunLab's ANI Calculator, for average nucleotide identity analysis
- GGDC 3.0, for digital DNA-DNA hybridization and genome relatedness estimation
- LEfSe, for biomarker and differential abundance analysis
- Chiplot, for data visualization
- 

### R packages

The following R packages were used for microbial community analysis, ecological statistics, mixed-effect modeling, post hoc comparisons, community assembly analysis, and microbial network analysis:

- `microeco`
- `vegan`
- `lme4`
- `lmerTest`
- `emmeans`
- `glmmTMB`
- `iCAMP`
- `SPIEC-EASI`

### Conda environments used

The following Conda environments were used in this workflow. Environment names may need to be adjusted depending on the local installation.

- `BBTools`
- `metawrap-env2`
- `CheckM`
- `Prokka_1.4`
- `gtdbtk-2.3.2`
- `eggnog`
- `carveme`
- `cobra`

> Note: Custom Python scripts were used for KEGG pathway completeness analysis. Users should modify input paths, database locations, and variable names according to their local environment.

## Ecological and statistical analysis in R

R was used for downstream microbial community analysis, ecological statistics, differential abundance analysis, mixed-effect modeling, and visualization. Packages used in the analysis included `microeco`, `vegan`, `lme4`, `lmerTest`, `emmeans`, `glmmTMB`, `iCAMP`, and `SPIEC-EASI`.

Typical analyses included:

- alpha and beta diversity analysis
- ordination analysis
- permutational multivariate analysis of variance
- differential abundance analysis
- linear mixed-effect models and generalized linear mixed models
- estimated marginal means and post hoc comparisons
- community assembly process analysis
- microbial association network inference

Example R package loading code:

```r
library(microeco)
library(vegan)
library(lme4)
library(lmerTest)
library(emmeans)
library(glmmTMB)
library(iCAMP)
library(SpiecEasi)
```

LEfSe and Chiplot were used as additional online platforms for biomarker discovery, differential abundance analysis, and visualization.

## Genome-scale metabolic modeling

Genome-scale metabolic models were reconstructed using CarveMe. The reconstructed GSMs were assembled into a shared-environment community model, with the sum of member biomass fluxes set as the objective function. Steady-state fluxes and flux variability were computed using COBRApy.

Flux balance analysis, parsimonious flux balance analysis, and flux variability analysis were performed under specified exchange constraints:

- FBA, for estimating feasible steady-state flux distributions
- pFBA, for estimating parsimonious flux distributions
- FVA, for estimating the variability range of individual reaction fluxes

Example CarveMe command:

```bash
conda activate carveme

carve input_genome.faa \
-o metabolic_modeling/input_genome.xml
```

Example COBRApy analysis outline:

```python
import cobra
from cobra.flux_analysis import pfba, flux_variability_analysis

model = cobra.io.read_sbml_model("metabolic_modeling/community_model.xml")

# Set exchange constraints according to the experimental or environmental condition.
# Example:
# model.reactions.EX_glc__D_e.lower_bound = -10

solution = model.optimize()
pfba_solution = pfba(model)
fva_result = flux_variability_analysis(model)
```

## Software, online platforms, and R packages

### Command-line software and bioinformatics tools
- MEGA 11, for phylogenetic or evolutionary analysis where applicable
- GNU Parallel, for parallel execution of batch jobs
- CarveMe, for genome-scale metabolic model reconstruction
- COBRApy, for constraint-based metabolic modeling, FBA, pFBA, and FVA
- Conda, for environment management
- BBTools / BBduk, for read quality control
- MetaWRAP, for metagenomic assembly, binning, and bin refinement
- MEGAHIT, used through MetaWRAP for metagenomic assembly
- SOAPdenovo v2.04, for genome assembly or comparative assembly analysis where applicable
- MetaBAT2, MaxBin2, and CONCOCT, for genome binning through MetaWRAP
- CheckM, for genome quality assessment
- Prokka, for genome annotation
- GTDB-Tk, for taxonomic classification
- eggNOG-mapper, for functional annotation


## Metagenomic and MAGs analysis Overview

The workflow includes the following major steps:

1. Quality control of raw reads using BBduk
2. Metagenomic assembly using MEGAHIT via MetaWRAP
3. Additional genome assembly or comparative analysis using SOAPdenovo v2.04 and MEGA 11, where applicable
4. Genome binning using MetaBAT2, MaxBin2, and CONCOCT via MetaWRAP
5. Bin refinement using `metawrap bin_refinement`
6. Genome quality assessment using CheckM
7. Genome annotation using Prokka and RAST
8. Taxonomic classification using GTDB-Tk
9. ANI and genome relatedness analysis using ChunLab's ANI Calculator and GGDC 3.0
10. Functional annotation using eggNOG-mapper and KEGG pathway completeness analysis
11. Microbial community, diversity, differential abundance, and ecological statistical analyses using R packages and online visualization/statistical platforms
12. Genome-scale metabolic model reconstruction using CarveMe, followed by community metabolic modeling using COBRApy

## Input data

Place raw paired-end reads in the `Rawdata/` folder using the following naming convention:

```text
sample1_1.fq.gz
sample1_2.fq.gz
sample2_1.fq.gz
sample2_2.fq.gz
```

The sample name is defined as the prefix before `_1.fq.gz` or `_2.fq.gz`.

## Workflow

The entire workflow can be executed step by step as described below.

### 1. Prepare sample list

```bash
cd Rawdata
for i in *_1.fq.gz; do echo "${i%_1.fq.gz}"; done > sample_names.txt
```

### 2. Quality trimming using BBduk

```bash
conda activate BBTools

cat sample_names.txt | parallel --link -j2 --delay 120 \
"bbduk.sh in1={}_1.fq.gz in2={}_2.fq.gz \
out1={}_1.fastq.gz out2={}_2.fastq.gz \
ktrim=r k=28 mink=12 hdist=1 tbo=t tpe=t \
qtrim=rl trimq=20 minlen=100 stats={}_detailed_stats.txt"
```

### 3. Assembly with MEGAHIT through MetaWRAP

```bash
conda activate metawrap-env2

cat sample_names.txt | parallel -j2 --delay 20 \
"metawrap assembly -1 {}_1.fastq.gz -2 {}_2.fastq.gz \
--megahit -o {}_new -m 700 -t 50"
```

### 4. Binning using MetaBAT2, MaxBin2, and CONCOCT

```bash
conda activate metawrap-env2

cat sample_names.txt | parallel -j2 --delay 20 \
"metawrap binning -a {}_new/final_assembly.fasta -o {} \
-t 20 -l 500 --metabat2 --maxbin2 --concoct \
{}_1.fastq.gz {}_2.fastq.gz"
```

### 5. Bin refinement

Bins were refined using a minimum completeness threshold of 50% and a maximum contamination threshold of 10%.

```bash
conda activate metawrap-env2

cat sample_names.txt | parallel -j2 --delay 20 \
"metawrap bin_refinement -o {}_refinement -t 20 \
-A {}/metabat2_bins/ -B {}/maxbin2_bins/ -C {}/concoct_bins/ \
-c 50 -x 10"
```

### 6. Genome quality assessment using CheckM

```bash
conda activate CheckM

checkm lineage_wf \
-f checkm_output/checkm_output.txt \
--tab_table \
-x fasta \
-t 20 \
--pplacer_threads 20 \
Genome/ checkm_output
```

### 7. Genome annotation using Prokka

```bash
conda activate Prokka_1.4

cd Genome
for i in *.fasta; do echo "${i%.fasta}"; done > sample_names.txt

cat sample_names.txt | parallel -j5 --delay 12 \
"prokka --outdir {}_prokka --prefix {} {}.fasta"
```

Genome annotation was also performed using the RAST online server where applicable.

### 8. Taxonomic classification using GTDB-Tk

```bash
conda activate gtdbtk-2.3.2

gtdbtk classify_wf \
-x fasta \
--genome_dir Genome/ \
--out_dir classify_wf \
--cpus 40 \
--pplacer_cpus 40
```

### 9. Functional annotation using eggNOG-mapper

```bash
conda activate eggnog

for i in *.faa; do echo "${i%.faa}"; done > sample_names.txt

cat sample_names.txt | parallel -j5 --delay 20 \
"emapper.py -i {}.faa --itype proteins \
--data_dir /path/to/eggnog-mapper-data \
-o {}.output --output_dir ./ \
-m diamond --cpu 20 --seed_ortholog_evalue 1e-5"
```

### 10. KEGG pathway completeness analysis

Custom Python scripts were used to estimate KEGG pathway completeness from functional annotation outputs. Before running these scripts, users should check and modify:

- input annotation file paths
- KEGG module or pathway database paths
- output directory paths
- sample names and variable names

Example command:

```bash
python scripts/kegg_pathway_completeness.py \
--input annotation/eggnog_outputs/ \
--database /path/to/kegg/database \
--output results/kegg_pathway_completeness.tsv
```

### 11. ANI and genome relatedness analysis

Average nucleotide identity and genome relatedness analyses were conducted using online tools, including ChunLab's ANI Calculator and GGDC 3.0. Input genome FASTA files were selected from the refined genome bins or MAGs after quality assessment.

## Notes

- Raw sequencing reads and large intermediate files are not included in this repository.
- External databases, including GTDB-Tk, eggNOG-mapper, KEGG, and other annotation databases, should be downloaded and configured separately.
- Paths in the scripts are examples and should be modified before use.
- Thread numbers and memory settings should be adjusted according to the local computing environment.
- Online tools may change over time; users should follow the current instructions provided by each website.

## Citation

If you use this workflow, please cite the associated study and the software tools used in the analysis.