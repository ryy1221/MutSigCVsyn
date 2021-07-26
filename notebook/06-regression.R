# load package
library(ggplot2)
require(foreign)
require(ggplot2)
require(MASS)
require(Hmisc)
require(reshape2)
library(dplyr)


# set working directory
setwd('/gpfs/group/epo2/default/yur97/gitlab/pcawg-to-mutsigcv/notebook')

# read file
fpath = './figure6/data'
fname = 'VCAN.csv'
df.vcan = na.omit(read.csv(file.path(fpath,fname)))

# describe exp_uq and transform into categorical
summ.exp = summary(df.vcan$exp_uq) 
boxplot(df.vcan$exp_uq,data= df.vcan)
cut1 = quantile(df.vcan$exp_uq,.25, na.rm = T)
cut2 = quantile(df.vcan$exp_uq,.5, na.rm = T)
cut3 = quantile(df.vcan$exp_uq,.75, na.rm = T)
df.vcan$exp_categ <- cut(df.vcan$exp_uq,
                       breaks = c(0, cut1,cut2,cut3,max(df.vcan$exp_uq, na.rm = T)), # cut points
                       right = TRUE # closed on the righ, open on the left
)
df.vcan$exp_ncateg = ntile(df.vcan$exp_uq, 4)

# Convert copy number variation into binary
df.vcan$cnv_binary[df.vcan$cnv>2] <- 1
df.vcan$cnv_binary[df.vcan$cnv==2] <- 0

# See the data
ggplot(df.vcan, aes(x = exp_categ, y = as.factor(has_syn), fill = exp_categ)) +   geom_boxplot(size = .75) +   
  facet_grid(has_nsyn ~ cnv_binary, margins = FALSE) +   
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

# Ordinal logisitic regression
ord.logit <- polr(exp_categ ~ as.factor(has_syn) + as.factor(has_nsyn) + as.factor(cnv_binary), data = df.vcan, Hess=TRUE)
summary(ord.logit)
ctable <- coef(summary(ord.logit))
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
ctable <- cbind(ctable, "p value" = p)
ctable
# odds ratio
exp(coef(ord.logit))
new_data <- data.frame("has_syn"= 0,"has_nsyn"="0","cnv_binary"="1")
round(predict(ord.logit,new_data,type = "p"), 3)

# evaluate logit model
sf <- function(y) {
  c('exp in (42.2,111]' = qlogis(mean(y == '(42.2,111]')),
    'exp in (0,3.73]' = qlogis(mean(y == '(0,3.73]')),
    'exp in (3.73,42.2]' = qlogis(mean(y == '(3.73,42.2]')))
}
na.omit(df.vcan)
s <- with(na.omit(df.vcan), summary(exp_categ ~ has_syn + has_nsyn + cnv_binary, fun=sf))
s


lm.co <- lm(exp_ncateg~has_syn + has_nsyn + cnv_binary, data=na.omit(df.vcan))
summary(lm.co)

lm.syn <- lm(exp_uq~has_syn, data=df.vcan)
summary(lm.syn)
ggplot(data = df.vcan, aes(y = exp_uq,x = cnv, color = has_nsyn)) +geom_point() +theme_minimal()

lm.cnv <- lm(exp_uq~cnv, data=df.vcan)
summary(lm.cnv)

glm.fit <- glm(exp_uq~ has_syn + has_nsyn + cnv, data=df.vcan, family = binomial)


ggplot(df.vcan, aes(x = has_syn, y = exp_uq)) +
  geom_boxplot(aes(group = 2),size = .75) +
  geom_jitter(alpha = .5) +
  facet_grid(has_nsyn ~ cnv, margins = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
