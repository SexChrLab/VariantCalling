# VariantCalling

Project to assess best practices for variant calling and filtering on human samples of African origin and across autosomes, the sex chromosome, and mtDNA.

## Install Conda and create conda environmnet

To run these scripts, you will need to install some programs

```
conda create --name varCalling
source activate varCalling
conda install -c conda-forge r-base
conda install -c r r-ggplot2
conda install -c bioconda gatk4=4.1.0.0
conda install -c bioconda snakemake
conda install -c conda-forge r-ggpubr
```

Note: If using gatk3 `conda install -c bioconda gatk`.

See: https://anaconda.org/bioconda/gatk
