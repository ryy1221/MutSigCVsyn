# PCAWG-to-MutSigCV

This project is for processing PCAWG somatic mutation and coverage files. The processed files are used as input for MutSigCV to identify significantly mutated gene in cancer.

Inputs requires:
1. Coverage: Retrieve .wig file for processing coverage files at [PCAWG portal](https://dcc.icgc.org/releases/PCAWG/consensus_snv_indel/wig_files)
2. Mutation: Retrieve PCAWG SNV MAF files(The ICGA part could be retrieved at [PCAWG portal](https://dcc.icgc.org/releases/PCAWG/consensus_snv_indel/), the TCGA part is controlled access so make sure you have the access, then go to GDC website to retrieve the data)
3. Reference files(Genome sequence, annotation and annotated intermediate files) availble upon requests
### Directory layout
    .
    ├── data                        # Input and reference files(All content available upon request, some are available to download in official website).  
    │   ├── anno_refs               # Including unprocessed annotation and anlysis inputs.  
    │   │   ├── gencode_v19         # Gencode v19 annotation and GRCh37 sequence file(downloaded).  
    │   │   ├── pcawg_rnaseq        # PCAWG RNA-seq files (tophat processed file: fpkm).  
    │   │   ├── other scatter files # CGC gene list(from COSMIC.  
    ├── out                  # Documentation files (alternatively `doc`).  
    ├── notbook              # Source files (alternatively `lib` or `app`).  
    ├── source               # Automated tests (alternatively `spec` or `tests`).  
    ├── supp                 # Tools and utilities. 
    └── README.md
