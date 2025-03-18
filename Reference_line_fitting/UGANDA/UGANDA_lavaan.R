library(lavaan)
library(blavaan)
library(MplusAutomation)
library(here)
library(dplyr)

data <- read.table(here("IRT_CFA - free all","Mplus", "Fitting", "Baseline", "UGANDA", "data.dat"), header = F, sep = "\t")
colnames(data) <- c("UAVLTOT", "UAVLTIR", "UAVLTDELR", "UAVLTREC", "GPEGDOM", "GPEGNDOM", "COLOR1", "COLOR2", "SYDM")


data_cohort <- data

model <- '
  MEM =~ UAVLTOT + UAVLTIR + UAVLTDELR + UAVLTREC
  SPEED =~ COLOR1 + COLOR2 + SYDM
  MOT =~ GPEGDOM + GPEGNDOM
  G =~ MEM + SPEED + MOT
'

mydp <- dpriors(nu="normal(0,100)", alpha="normal(0,100)", lambda="normal(0,100)", beta="normal(0,100)")
bfit <- bcfa(model, data = data_cohort, burnin = 3000, sample = 3000, dp = mydp, meanstructure = TRUE, seed = 0)
summary(bfit)


null_model <- '
  MEM =~ 1*UAVLTOT + 1*UAVLTIR + 1*UAVLTDELR + 1*UAVLTREC
  MEM ~~ 0*MEM

  SPEED =~ 1*COLOR1 + 1*COLOR2 + 1*SYDM
  SPEED ~~ 0*SPEED

  MOT =~ 1*GPEGDOM + 1*GPEGNDOM
  MOT ~~ 0*MOT
'
bfit_null <- bcfa(null_model, data = data_cohort, burnin = 3000, sample = 3000, dp = mydp, meanstructure = TRUE, seed = 0)
summary(bfit_null)

bfit_all <- blavFitIndices(bfit, baseline.model = bfit_null, pD = c("loo"), rescale = "devM", fit.measures = c("BCFI", "BGammaHat"))
summary(bfit_all, central.tendency = c("mean","median","mode"), prob = .95)
