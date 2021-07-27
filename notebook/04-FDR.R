### This script is for adjusting p-values using control (non-expressed gene result)
# This script fit kernel density curve to cdf of  p-val distribution
library(sROC)
library(ggplot2)
require(gridExtra)
library(dplyr)

#-----Functions-----
kcdf.estimate.all = function(df.all, df.ne,df.sig.exp, bw.all, bw.ne, proba.h0){
  x.all = na.omit(df.all$p)
  x.ne = na.omit(df.ne$p)
  x.sig.exp = na.omit(df.sig.exp$p)
  
  print('Start Calculating CDF')
  fn.all= kCDF(x.all,bw = bw.all, xgrid =x.all,from = 0, to = 1)
  # For the nonexp distribution, we estimate the value of pvalues in all gene run
  fn.ne = kCDF(x.ne, bw = bw.ne, xgrid = x.all, from = 0, to = 1) 
  print('Finish Caluclating CDF')
  
  for(i in 1:(nrow(df.sig.exp))){
    x = df.sig.exp$p[i] # observed p-value
    denom = fn.all[["Fhat"]][i] # Denominator = p(x)
    proba.FP.x = fn.ne[["Fhat"]][i]
    numer = proba.h0*proba.FP.x # Numerator = P(H0,x)
    proba.x.FP = numer/denom
    df.all$FDR[i] = proba.x.FP
  }
  print('Finish Calulating FDR')
  return(list(f.all = fn.all, f.ne = fn.ne ,df.res = df.all,bandwidth.all = bw.all, bandwidth.ne = bw.ne))
}

kcdf.estimate.sig = function(df.all, df.ne,df.sig.exp, bw.all, bw.ne, proba.h0){
  x.all = na.omit(df.all$p)
  x.ne = na.omit(df.ne$p)
  x.sig.exp = na.omit(df.sig.exp$p)
  
  print('Start Calculating CDF')
  fn.all= kCDF(x.all,bw = bw.all, xgrid =x.sig.exp,from = 0, to = 1)
  # For the nonexp distribution, we estimate the value of pvalues in all gene run
  fn.ne = kCDF(x.ne, bw = bw.ne, xgrid = x.sig.exp, from = 0, to = 1) 
  print('Finish Caluclating CDF')
  
  for(i in 1:(nrow(df.sig.exp))){
    x = df.sig.exp$p[i] # observed p-value
    denom = fn.all[["Fhat"]][i] # Denominator = p(x)
    proba.FP.x = fn.ne[["Fhat"]][i]
    numer = proba.h0*proba.FP.x # Numerator = P(H0,x)
    proba.x.FP = numer/denom
    df.sig.exp$FDR[i] = proba.x.FP
  }
  print('Finish Calulating FDR')
  return(list(f.all = fn.all, f.ne = fn.ne ,df.res = df.sig.exp,bandwidth.all = bw.all, bandwidth.ne = bw.ne))
}

#-----Read file-----
dir.sig = './figure4/p-val_distribution'
dir.out.fig = './figure4/fdr_figs'
dir.out.func = './figure4/fdr_funcs'
dir.out.res = './figure4/fdr_res'
feature_type = 'histology'; run = 'cohort_072221'; threshold = 1;
df.all.feat = read.csv(file.path(dir.sig,paste0(feature_type,'.syn.df_all_forFDR.',
                                                run, '.',threshold,'.csv')),sep = ',')
df.all.feat = df.all.feat[order(df.all.feat$p),]
df.exp.feat = df.all.feat[df.all.feat$exp.nonexp =='exp',] ; df.ne.feat = df.all.feat[df.all.feat$exp.nonexp =='nonexp',]

###-----Global FDR-----
df.sig = na.omit(df.all.feat[df.all.feat$p<0.05,])
res = kcdf.estimate.sig(df.all.feat,df.ne.feat,df.sig, 0.00001,0.00001,1)

###-----Result part-----
df.res.sig = res$df.res
df.res.sig <- dplyr::filter(df.res.sig, !grepl("^PCDH", gene))

x= res$f.all$x; y.all = res$f.all$Fhat; y.ne = res$f.ne$Fhat
df = data.frame(x, y.all, y.ne)
par(mfrow=c(1,2))
p1 = ggplot(df, aes(x))+
  geom_point(aes(y = y.all, colour = 'all genes'))+
  geom_line(y = y.all, colour = 'grey')+
  geom_point(aes(y = y.ne, colour = 'non exp genes'))+
  geom_line(y = y.ne, colour = 'grey')+
  coord_cartesian(xlim=c(0,0.05), ylim = c(0,0.05))+
  labs(colour = '', x ='x(p-val)', y = 'F(x)', title = 'Kernel Estimated Cdf(zoom in )')+
  theme_classic()
p2 = ggplot(df, aes(x))+
  geom_point(aes(y = y.all, colour = 'all genes'))+
  geom_line(y = y.all, colour = 'grey')+
  geom_point(aes(y = y.ne, colour = 'non exp genes'))+
  geom_line(y = y.ne, colour = 'grey')+
  coord_cartesian(xlim=c(0,0.005), ylim = c(0,0.005))+
  labs(colour = '', x ='x(p-val)', y = 'F(x)', title = 'Kernel Estimated Cdf')+
  theme_classic()
grid.arrange(p2, p1, ncol=2)
fig = arrangeGrob(p2, p1, ncol=2)
saveRDS(res,file.path(dir.out.func,'fdr.063021.rds'))
write.csv(df.res.sig,file.path(dir.out.res, 'fdr.063021.csv'))
ggsave(file=file.path(dir.out.fig,'fdr.063021.png'), fig)

### The CDF function of all values---------
set.seed(500)
kcdf.estimate.plt = function(df.all, df.ne,df.sig.exp, bw.all, bw.ne, proba.h0){
  x.all = na.omit(df.all$p)
  x.ne = na.omit(df.ne$p)
  x.plt = sample(x.all,10000)
  x.plt = x.plt[order(x.plt)]
  
  print('Start Calculating CDF')
  fn.all= kCDF(x.all,bw = bw.all, xgrid =x.plt)
  # For the nonexp distribution, we estimate the value of pvalues in all gene run
  fn.ne = kCDF(x.ne, bw = bw.ne, xgrid = x.plt) 
  print('Finish Caluclating CDF')
  
  for(i in 1:(nrow(df.sig.exp))){
    x = df.sig.exp$p[i] # observed p-value
    denom = fn.all[["Fhat"]][i] # Denominator = p(x)
    proba.FP.x = fn.ne[["Fhat"]][i]
    numer = proba.h0*proba.FP.x # Numerator = P(H0,x)
    proba.x.FP = numer/denom
    df.sig.exp$FDR[i] = proba.x.FP
  }
  print('Finish Calulating FDR')
  # df.pass = df.sig.exp[df.sig.exp$FDR<0.1,]
  #  return(nrow(df.pass))
  return(list(f.all = fn.all, f.ne = fn.ne ,df.res = df.sig.exp,bandwidth.all = bw.all, bandwidth.ne = bw.ne))
}
res = kcdf.estimate.plt(df.all.feat,df.ne.feat,df.sig, 0.00001,0.00001,0.99)

###-----Plot-----
x= res$df.all$x; y.all = res$f.all$Fhat; y.ne = res$f.ne$Fhat
df = data.frame(x, y.all, y.ne)
p2 = ggplot(df, aes(x))+
  geom_point(aes(y = y.all, colour = 'all genes'))+
  geom_line(y = y.all, colour = 'grey')+
  geom_point(aes(y = y.ne, colour = 'non exp genes'))+
  geom_line(y = y.ne, colour = 'grey')+
  coord_cartesian(xlim=c(0,0.05), ylim = c(0,0.05))+
  labs(colour = '', x ='x(P-value)', y = 'F(x)', title = 'Zoom in')+
  theme_classic()+ theme(legend.position = "none",text=element_text(size=8)) 
p1 = ggplot(df, aes(x))+
  geom_point(aes(y = y.all, colour = 'all genes'))+
  geom_line(y = y.all, colour = 'grey')+
  geom_point(aes(y = y.ne, colour = 'non exp genes'))+
  geom_line(y = y.ne, colour = 'grey')+
  coord_cartesian(xlim=c(0,1), ylim = c(0,1))+
  labs(colour = '', x ='x(P-value)', y = 'F(x)', title = 'Kernel Estimated CDF')+
  theme_classic()+theme(legend.position = "none",text=element_text(size=15)) +
  geom_rect(aes(
      xmin = -0.005,xmax = 0.01,ymin = -0.005,ymax = 0.01),fill = NA,color = "blue",size = 0.25
  )+
  annotation_custom(ggplotGrob(p2), xmin = 0.05, xmax = 0.5, 
                    ymin = 0.1, ymax = 0.6)
print(p1)
ggsave(file=file.path(dir.out.fig,'fdr.all.sample.070721.png'))
