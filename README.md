# SEA_project
All scripts for data and IBD analysis

The analysis pipeline consists of three main stages, ranging from raw sequencing data processing to graphical visualization.

### 1. Variant Calling
Post-sequencing data were processed using four dedicated scripts located in the `Scripts/GATK/` directory. These scripts generate **multi-sample VCF files** required for downstream analysis.

### 2. IBD Inference
Pairwise Identity by Descent (IBD) was inferred using **hmmIBD**. The multi-sample VCFs were first converted using `vcf2hmm.py` before running the `hmmIBD` executable.

> **Tool Reference:** > Schaffner, S.F., Taylor, A.R., Wong, W. *et al.* hmmIBD: software to infer pairwise identity by descent between haploid genotypes. *Malar J* **17**, 196 (2018).  
> 🔗 [Source Code & Documentation](https://github.com/glipsnort/hmmIBD) | [Paper](https://doi.org/10.1186/s12936-018-2349-7)

### 3. Visualization
The resulting IBD data were combined with metadata files (stored in `Data/Metadata/`). These datasets were processed using the R script `Scripts/IBD/plot_IBD_graph.R` to generate the final graphical visualizations.
