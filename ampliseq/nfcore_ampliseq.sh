#! /bin/bash

nextflow run nf-core/ampliseq -r 2.8.0 --input Samplesheet.tsv --metadata Depletion_Metadata.tsv --outdir DEPL_ampliseq -profile singularity --picrust --metadata_category Treatment,Timepoint --FW_primer AGAGTTTGATYMTGGCTCAG --RV_primer ATTACCGCGGCKGCTGG --trunclenf 275 --trunclenr 265 --skip_cutadapt --max_cpus 12 -resume