### This script is for adjusting p-values using control (non-expressed gene result)
# This script fit kernel density curve to cdf of  p-val distribution
library(sROC)
library(ggplot2)

### Sig genes input folder------------------------------------
# dir.sig = '/gpfs/group/epo2/default/yur97/gitlab/pcawg-to-mutsigcv/notebook/figure4/p-val_distribution'
dir.sig = './figure4/p-val_distribution'

# Read file
feature_type = 'histology'
df.all.feat = read.csv(file.path(dir.sig,paste0(feature_type,'.syn.df_all_forFDR.062121.csv')),sep = ',')
df.all.feat = df.all.feat[order(df.all.feat$p),]
df.exp.feat = df.all.feat[df.all.feat$exp.nonexp =='exp',] ; df.ne.feat = df.all.feat[df.all.feat$exp.nonexp =='nonexp',]
df.sig= na.omit(df.exp.feat[df.exp.feat$p<0.05,])
print('Finish Reading files')

x.all = na.omit(df.all.feat$p)
x.ne = na.omit(df.ne.feat$p)

print('Start Calculating CDF')
fn.all= kCDF(x.all,xgrid =x.all,from = 0, to = 1)
print('Finish Calculating CDF')
saveRDS(fn.all,file.path('./figure4/FDR_calc/default_fn.all.RDS'))
print('Finish saving function')
