# Script for handing stayer prop data
# 13 / 1/ 2015

rm(list=ls())

require(stringr)
require(dplyr)


raw_simplified <- read.csv(
    file="Scottish_Stayer_Props/excel/stayers_2001_less_simplified.csv",
    stringsAsFactors=FALSE
) 

dim(raw_simplified)

for (i in 2:ncol(raw_simplified)){
    raw_simplified[,i] <- raw_simplified[,i] %>% 
        str_replace_all("-", "0") %>%
        as.numeric()
}

derived <- raw_simplified %>% mutate(
    datazone=datazone,
    stayer= Inflow.All.people.in.households.Lived.at.same.address + Inflow.All.people.in.households.Moved.within.Datazone.area +
        Inflow.All.people.in.communal.establishments.Lived.at.same.address + Inflow.All.people.in.communal.establishments.Moved.within.Datazone.area
    )

derived$total= apply(raw_simplified[,-c(1, 16)], 1, sum)

derived <- derived %>% select(datazone, stayer,total)

derived <- derived %>% mutate(
    mover = total - stayer, 
    one_year_mover_prop = mover / total,
    ten_year_stayer_prop = (1 - one_year_mover_prop)^10
    )

write.csv(derived, file="Scottish_Stayer_Props/derived/ten_year_stayer_prop_all_ages.csv", row.names=F)


# NOTE : Qualifications figures are for people aged 25 to 74
# And so are employment status figures

# To do: the above, but for different age groups and class groups

# Combinations

# Age
# 0 - 18
# 19 - 29
# 30 - 44
# 45 - 64
# 65 plus


# let's just produce a simple adjustment curve





