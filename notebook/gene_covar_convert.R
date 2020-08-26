# this code is for 

#get gene list
gene.list = read.csv('~/Downloads/pcawg_genelist.txt', header = F)

### BioMart get paralogs
library(biomaRt)
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host= "grch37.ensembl.org")
searchAttributes(mart = human, pattern = "syn")

res = getBM(attributes = c("external_gene_name","external_synonym"),
            filters = "external_gene_name",
            values = gene.list,
            mart = human)

### Format result
agg.res = res[!duplicated(res$external_gene_name),]
agg.res[,"external_synonym"] = aggregate(external_synonym~external_gene_name, data=res, toString)[,2]
agg.res[,"external_synonym"] = gsub(",", "\t", agg.res$external_synonym)
agg.res[,"external_synonym"] = gsub("^\t ", "", agg.res$external_synonym)

### Write output
write.table(res, 'gene.covariates.external.txt', append = FALSE, sep = "\t", dec = "\t",
            row.names = FALSE, col.names = c("org","syn"), quote = FALSE)
