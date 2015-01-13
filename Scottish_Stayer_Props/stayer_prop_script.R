# Script for handing stayer prop data
# 13 / 1/ 2015

rm(list=ls())


raw_simplified <- read.csv(
    file="Scottish_Stayer_Props/excel/stayers_2001_simplified.csv",
    stringsAsFactors=FALSE
)

require(stringr)
raw_simplified[,2]  <- str_replace_all(raw_simplified[,2], "-", "0")
raw_simplified[,3]  <- str_replace_all(raw_simplified[,3], "-", "0")
raw_simplified[,2] <- as.numeric(raw_simplified[,2])
raw_simplified[,3] <- as.numeric(raw_simplified[,3])

head(raw_simplified)


