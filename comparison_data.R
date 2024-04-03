# model application - RISCC Deliverable Jan 2024
source("functions_PICmodel.R")
library(dplyr)
library(survminer)

#---------------------------
# POBASCAM:
pob_int <- readRDS("pob_int.RData")
pob_int$age <- ifelse(pob_int$age > 39, 1, 0)
pob_int$agegrp <- NULL
pob_int$cyt <- pob_int$cyt - 1 # makes more sense to have cyt = 0/1/2 rather than 1/2/3
pob_int$cyt_abnormal <- ifelse(pob_int$cyt > 0, 1, 0 )
colnames(pob_int) <- c("left", "right", "z", "age", "hpv16", "cyt", "cyt_abnormal") # rename hpv16 column name
pob.fit.h <- model.fit(c("hpv16"), c(), c("age", "hpv16", "cyt_abnormal"), pob_int, include.h = T, silent = F, two.step.h = T, short.runs=30, short.iter=20)
pob.fit.h$summary

#---------------------------
# VUSA-Screen:
vusa_int <- readRDS("vusa_int.RData")
vusa_int$age <- ifelse(vusa_int$age > 39, 1, 0)
vusa_int$cyt <- case_when(
  vusa_int$cyt =='pap1' ~1,
  vusa_int$cyt =='BMD'~2,
  vusa_int$cyt =='>BMD'~3
)
vusa_int$cyt <- vusa_int$cyt - 1
vusa_int$cyt_abnormal <- ifelse(vusa_int$cyt > 0, 1, 0 )

vusa.fit.h <- model.fit(c("hpv"), c(), c("age", "hpv", "cyt_abnormal"), vusa_int, include.h = T, silent = F,
                        two.step.h = T, short.runs=20, short.iter=10)
vusa.fit.h$summary

#---------------------------
# slovenia:
slov_int <- readRDS("slov_int.RData")
slov_int$cyt_abnormal <- ifelse(slov_int$cyt_bmd | slov_int$cyt_severe, 1, 0 )

slov.fit.h <-  model.fit(c("hpv16"), c(), c("age", "hpv16", "cyt_abnormal"), slov_int, include.h = F,
                         silent = F, two.step.h = F, short.runs=30, short.iter=20)
slov.fit.h

#---------------------------
# Italy:
italy_int <- readRDS("italy_int.RData")
italy_int$cyt_abnormal <- ifelse(italy_int$cyt_bmd | italy_int$cyt_severe, 1, 0 )

italy.fit.h <-  model.fit(c("hpv16"), c(), c("age", "hpv16", "cyt_abnormal"), italy_int, include.h = T,
                          silent = F, two.step.h = T, short.runs=30, short.iter=10)
italy.fit.h

#---------------------------
# Sweden:
swed_int <- readRDS("swed_int.RData")
swed_int <- swed_int[swed_int$left!=0 | swed_int$right!=Inf,]
swed_int <- swed_int[(swed_int$right - swed_int$left > 0.01) | swed_int$right ==0,]
swed_int$cyt_abnormal <- ifelse(swed_int$cyt_bmd | swed_int$cyt_severe, 1, 0 )

swed.fit.h <-  model.fit(c("hpv16"), c(), c("hpv16", "cyt_abnormal"), swed_int, include.h = F,
                         silent = F, two.step.h = F, short.runs=20, short.iter=5)

swed.fit.h
