# SEA_project
All scripts for data and IBD analysis

The analysis pipeline consists of three main stages, ranging from raw sequencing data processing to graphical visualization.

### 1. Variant Discovery Pipeline (GATK)
Raw sequencing data were processed using a custom pipeline based on GATK best practices. The workflow is automated through four shell scripts located in `Scripts/GATK/`:

1.  **Alignment & Pre-processing** (`1_fastq_to_bam.sh`): Converts raw FASTQ files into aligned BAM format.
2.  **Variant Calling** (`2_variants_calling.sh`): Performs initial variant calling per sample to generate GVCFs.
3.  **Joint Genotyping** (`3_combine_gvcf.sh`): Merges individual GVCF files into a comprehensive multi-sample VCF.
4.  **Variant Filtering** (`4_filtering_variants.sh`): Applies quality filters to retain high-confidence variants for downstream analysis.

### 2. IBD Inference
Pairwise Identity by Descent (IBD) was inferred using **hmmIBD**. The filtered multi-sample VCFs were first converted using `vcf2hmm.py` before running the `hmmIBD` executable.

> **Tool Reference:**
> Schaffner, S.F., Taylor, A.R., Wong, W. *et al.* hmmIBD: software to infer pairwise identity by descent between haploid genotypes. *Malar J* **17**, 196 (2018).
> 🔗 [Source Code & Documentation](https://github.com/glipsnort/hmmIBD) | [Paper](https://doi.org/10.1186/s12936-018-2349-7)

### 3. Visualization
The resulting IBD outputs were combined with metadata files (stored in `Data/Metadata/`). These datasets were processed using the R script `Scripts/IBD/plot_IBD_graph.R` to generate the final graphical visualizations.
