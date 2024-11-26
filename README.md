---
title: 'Effect of Antibiotic Gut Microbiome Depletion on TBI'
disqus: hackmd
---

Effect of Antibiotic Gut Microbiome Depletion on TBI
===

## Table of Contents

[TOC]

## Introduction

In this paper we used Illumina short-read 16S rRNA sequencing and ONT nanopore long-read metagenomic sequencing to characterize differences in the gut microbiome on TBI mice, all the while depleting one groups microbiome with an **intense** cocktail of antibiotics. 

>A summary of this dataset:
>2. The 16S data visualization was performed in Rstudio and the code is available on our [GitHub pages site]().
>3. Nanopore metagenomic data was sequenced with R9.4.1 chemistry but basecalled with dorado v0.7.0 using the super accuracy model.
>4. A custom analysis pipeline was used on the nanopore metagenomic data and the nitty gritty is detailed below!

Dependencies
---
This was run within two mamba environments the first containing these dependencies
```
# 2024_11_26 metassem env creation
mamba create -n metassem

mamba activate metassem

mamba install -c bioconda nanoplot hostile flye chopper medaka

# exact environment dependencies can be installed using this text file on our github
metassem.txt
```
```
# 2024_11_26 aviary env creation
mamba create -n aviary -c bioconda aviary

mamba activate aviary

aviary --version
version 0.10.0
```

Analysis Workflow
---
```sequence
Note left of Input: Raw pod5 data
Input->Output: Basecalling w/ dorado
Note right of Output: Fastq data
Output->Input: QC & filt w/ NanoPack2 and Hostile
Note left of Input: Clean fastqs
Input->Output: Metagenomic assembly w/ flye
Note right of Output: Assembled reads
Output->Input: Bin and refine w/ aviary recover
Note left of Input: Rough MAGs
Input->Output: MAG stats w/ CheckM2
Note right of Output: Final MAGs
Output->Output: Gene annotation
```

Example code
---
```
# I use a special env to basecall on my HPC
mamba activate dorado

# To process the raw pod5 data (gotta use a GPU)
dorado basecaller ...

# now further process with other envs
```

Once fastqs are obtained use this!
```
# Activate metassem env made earlier
mamba activate metassem

# Let's check our data visually (with host still present)
NanoPlot ...

# Now we can filter out the host DNA
hostile ...

# And filter by Qscore / read lengths
chopper

# Then we can assemble with metaflye

flye --meta ...

# After assembly is complete, polish with medaka (idk if --bacteria works but its a new argument I wanna try)
medaka --bacteria
```

Once a polished assembly has been obtained we now use the aviary recover pipeline for the rest of the analysis.
```
# Activate the env
mamba activate aviary

# Now run the pipeline
aviary recover ...
```

:::info
This should give you final MAGs with stats and annotations!
:::

> Read more about mermaid here: http://mermaid-js.github.io/mermaid/
> 
###### tags: `Templates` `Documentation`


