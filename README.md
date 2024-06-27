# MutSigCVsyn
MutSigCVsyn is an algorithm identified from [MutSigCV, Lawrence et al 2013](https://www.nature.com/articles/nature12213) for identifying synonymous cancer drivers. Running MutSigCVsyn using whole genome cancer mutation and coverage file. The current version is majorly designed for utilizing [PCAWG](https://www.nature.com/articles/s41586-020-1969-6) database. Other data can be used by making modification. 

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
    │   ├── MutSigCVsyn_inputs          # MutSigCVsyn other inputs: covariate file, mutational dictionary file, chr19 chain files etc.
    │   ├── proc_refs                   # All processed annotation files, final patient,gene and histology set(available upon request)
    ├── out                             # MutSigCVsyn output folder
    │   ├── syn                         # MutSigCVsyn synonymous output
    │   ├── nsyn                        # MutSigCVsyn non-synonymous output
    │   ├── anlyze                      # Pickled result folder paths for quickly reading result
    │   ├── example                     # Output folder for example MutSigCVsyn test run
    ├── source                          # MutSigCVsyn source code and input processing scripts
    │   ├── MutSigCVsyn                 # MutSigCVsyn source code
    │   ├── MutSigCVsyn_nsyn            # MutSigCVsyn source code for identifying non-synonymous candidates
    │   ├── run_MutSigCVsyn*            # Scripts for running MutSigCVsyn
    │   │   ├── histology_syn           # Identify synonymous candidates in PCAWG histologies, include a source python code for parsing inputs, a pbs file and a jupyternotebook for final submitting multiple histologies in parallel. An example notebook can show how to only run for 1 histology. 
    │   │   ├── histology_syn           # Same as the histology_syn but the nonsynonymous version
    │   ├── other scatter files         # Other scatter files for processing annotation, coverage and mutation files. Processed files are stored in ./data/proc_refs/
    ├── environment.yml    # conda environment
    ├── requirements.txt   # package requirements
    └── README.md

### How to set up the environment
```
conda env create -f environment.yml
conda activate Msigsyn
```
### How to run MutSigCVsyn(Example for single cohort, running in Jupyter notebook)
1. Obtain input files(available upon request)
    1. Example coverage file should be stored at `./data/cov/example`
    2. Example mutation annotation file should be stored at `./data/maf/example`
    3. Other files are stored at `./data/MutSigCVsyn_inputs`
2. Go to `./source/run_MutSigCVsyn/histology_syn/Example-single_MutSigCvsyn.ipynb` and run the notebook.
3. Expect less than 15min to finish.(Options to specify output file name can be find in the notebook)
4. Output files are stored at `../out/example/`

### Alternative option to run MutSigCVsyn by parallel for multiple histology types
1. Source code are stored at `./source/run_MutSigCVsyn/histology_syn` (For non-synonymous result, there is a parallel folder `histology_nsyn`). Jobs can be submitted in Jupyter notebook(`submit_MutSigCVsyn.ipynb`, where you can input parameters and test), but bash command can also be used. 
2. The main bash command is `python run_histology_syn-new.py ${HISTOLOGY_NAME} --c ${COVARIATE} --hyperm ${HYPERMUTATOR}`, which is already included if you run the Jupyter notebook. 

| Parameter | Value | Annotation/Description | 
| ------------- | ------------- | -------------------- | 
| HISTOLOGY_NAME  | The PCAWG histology names  | Also linked to the name of maf and coverage files |
| COVARIATE  | y/n;Yes/No | Whether to use the new covariate file for Gencode v10 |
| HYPERMUTATOR |y/n;Yes/No | Whether to exclude hypermutators in the cohort.Yes means excluding |

### Processed and ready-to-use files available upon request(except for controlled access data), but here are some public available resources for reference:
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
