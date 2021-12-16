# PCAWG-to-MutSigCV

This project is for processing PCAWG somatic mutation and coverage files. The processed files are used as input for MutSigCV to identify significantly mutated gene in cancer.

### Some raw inputs that are public(processed and ready-to-use files available upon request):
1. Coverage: Retrieve .wig file for processing coverage files at [PCAWG portal](https://dcc.icgc.org/releases/PCAWG/consensus_snv_indel/wig_files)
2. Mutation: Retrieve PCAWG SNV MAF files(The ICGA part could be retrieved at [PCAWG portal](https://dcc.icgc.org/releases/PCAWG/consensus_snv_indel/), the TCGA part is controlled access so make sure you have the access, then go to GDC website to retrieve the data)
3. Covariate and mutation dictionary input files. All of them are available to download at [MutSigCV website](https://software.broadinstitute.org/cancer/cga/mutsig). **Updated covariate file, mutation type dictionary file available upon request**
4. Reference files(Genome sequence, annotation and annotated intermediate files), processed files availble upon requests. Some raw files are available to found at:
    1. [Gencode v19 annotation and GRCh37 sequence file](https://www.gencodegenes.org/human/release_19.html)
    2. [PCAWG RNA-seq files](https://dcc.icgc.org/releases/PCAWG/transcriptome/gene_expression)
    3. [Cancer Gene Census gene list v95](https://cancer.sanger.ac.uk/census)
    4. [SynMICdb whole database](http://synmicdb.dkfz.de/rsynmicdb/)
    5. [PCAWG patient fusion file](https://dcc.icgc.org/releases/PCAWG/transcriptome/fusion)
    6. [PCAWG patient sample info- supplementary table1 of the PCAWG paper](https://www.nature.com/articles/s41586-020-1969-6)
    7. [PCAWG cancer drivers](https://dcc.icgc.org/releases/PCAWG/driver_mutations)



### Directory layout
    .
    ├── data                            # Input and reference files(All content available upon request, some are available to download in official website).  
    │   ├── anno_refs                   # Including unprocessed annotation and anlysis inputs.  
    │   │   ├── gencode_v19             # Gencode v19 annotation and GRCh37 sequence file.  
    │   │   ├── pcawg_rnaseq            # PCAWG RNA-seq files.  
    │   │   ├── other scatter files     # CGC gene list, PCAWG RNA-seq/fusion/cancer driver/sample info, CGC genes, SynMICdb database.
    │   ├── cov                         # Coverage and related files
    │   │   ├── wig_input               # Raw decompressed wig track file downloaded from PCAWG(see above for download link, due to the extremely large datasize, this is not available upon request)
    │   │   ├── intermediate            # Processed positional wig track file and calculated individual patient coverage file (available upon request)
    │   │   ├── histology               # Finished patient coverage files by PCAWG histologies (ready-to-use, available upon request)
    │   │   ├── example                 # Example coverage files for MutSigCVsyn test run(ready-to-use, available upon requestion)
    │   ├── maf                         # Mutation annotation and related files
    │   │   ├── maf_raw                 # Raw maf files from ICGC and TCGA(see above for download link, require TCGA access)
    │   │   ├── intermediate            # Processed individual patient mutation annotation files, with and without mutation categ information (available upon request)
    │   │   ├── histology               # Finished patient mutation annotation files by PCAWG histology (ready-to-use, available upon request)
    │   │   ├── example                 # Example mutation annotation files for MutSigCVsyn test run (ready-to-use, available upon request)
    │   ├── MutSigCVsyn_inputs          # Mutation 

















    ├── out                  # Documentation files (alternatively `doc`).  
    ├── notbook              # Source files (alternatively `lib` or `app`).  
    ├── source               # Automated tests (alternatively `spec` or `tests`).  
    ├── supp                 # Tools and utilities. 
    └── README.md
