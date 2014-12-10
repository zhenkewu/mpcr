rm(list=ls())

library(sandwich)
library(mpcr)

## read in data set:
#dat = read.csv('/Users/zhenkewu/Dropbox/ZW/working_projects/ZW_R_packages/mpcr_test/data/data_5.csv',header=TRUE)

res = mpcr(dat,
           arm = "tx",
           cluster = "team",
           pair    = "sitenew",
           outcome = "sf36pcs32",
           X_nm_all = c("race","ageatint","hcc","livesalone","education","s10",
                        "sf36mcs","sf36pcs","gender","h1"),
           X_nm_binary = c("livesalone","education","gender"),
           X_nm_cat    = c("race","s10","h1"),
           X_nm_cont   = c("ageatint","hcc","sf36mcs","sf36pcs"),
           BAYES_STATUS = FALSE,
           digit.round = 2)

cat("============See figures and tables.===========","\n")

print(res)


