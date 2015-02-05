# tidier script for running contraction mapping algorithm as part of Timmins' structural sorting models

# Jon Minton
# 5/2/2015


rm(list=ls())

require(compiler)
require(Rcpp)
require(plyr)
require(dplyr)
require(stringr)
require(tidyr)
require(ggplot2)
require(quantreg)
require(repmis)
# functions
sourceCpp("contraction_mapping_algorithm/contraction_mapping_algorithm_in_cpp.cpp")
source("contraction_mapping_algorithm/moving_cost_function.r")



# Management of Scottish data ---------------------------------------------


# Files needed: tenure type by datazone


tenure_households <- source_DropboxData(
    file="tenure_households.csv",
    key="kng5wc40le9kapj"    
) %>% tbl_df() %>% select(
    datazone, year, 
    all_households=HO.allhouseholds,
    council_houses=HO.council,
    rented_from_employer=HO.employ,
    owned_with_mortgage=HO.ownmortloan,
    owned_outright=HO.ownoutright,
    private_rented=HO.privlet,
    rented_from_relative=HO.relative,
    shared_ownership=HO.sharedown,
    other_social_rented=HO.social
) %>% 
    mutate(
        social=council_houses + other_social_rented,
        rented=rented_from_employer + private_rented+ rented_from_relative,
        owned=owned_with_mortgage + owned_outright + shared_ownership
    ) %>%
    mutate(
        council_houses=council_houses/all_households,
        rented_from_employer=rented_from_employer/all_households,
        owned_with_mortgage=owned_with_mortgage/all_households,
        owned_outright=owned_outright/all_households,
        private_rented=private_rented/all_households,
        rented_from_relative=rented_from_relative/all_households,
        shared_ownership=shared_ownership/all_households,
        other_social_rented=other_social_rented/all_households,
        social = social/all_households,
        rented = rented/all_households,
        owned=owned/all_households
    )



## An issue with adapting the model to the Scottish context becomes with classifying the social group: neither 
# 'renters' nor 'owners' in the groupings above. 

# Median values 


house_prices <- source_DropboxData(
    file="house_sales_and_prices.csv",
    key="pypifzbldk1wa6a"
    ) %>% 
    tbl_df() %>%
    select(datazone, year, 
           lower_quartile=HO.hpricelquartile,
           median=HO.hpricemedian,
           upper_quartile=HO.hpriceuquartile,
           mean=HO.hpricemean,
           n_sales=HO.hsalesno
           )

#second problem: what about areas with no house prices recorded?

 # One option
#  find contiguous regions, and use average of these regions
#  however, what if there really are very big disparities 
# between neighbouring regions (as in many parts of Glasgow)?
 

# Population counts


# Derived data: see script 
# generate_persons_by_gender_age_and_year.r
# to reproduce from SNS files

persons <- source_DropboxData(
    file="persons_by_gender_year_and_age.csv",
    key="p134pw625yd4f80"
    ) %>% tbl_df()



################

# Number of households only available for 2001
# Want to know it for 2011 as well

# Answer: use 2011 census table

households_2011 <- source_DropboxData(
    file = "QS116SC.csv",
    key="hb0phfnf38c5si3"
    ) %>% tbl_df() %>%
    select(datazone=X, all_households=All.households) %>%
    filter(grepl("^S01", datazone)) 

households_2011 <- households_2011 %>% 
    mutate(year=2011) %>%
    select(datazone, year, all_households)

households_2001 <- tenure_households %>% 
    select(datazone, year, all_households)

households <- households_2001 %>% bind_rows(households_2011)

rm(households_2001, households_2011)

# Indirect appraoch would be to estimate it based on number of persons
# Assuming number of persons stays constant over time. 


# Stayer probabilities: 
# First pass: 
#  - the median age in Scotland in 2001 was 38
#  - the estimated 10 year stayer proportion for this group is 0.82
prent <- tenure_households %>%
    select(owned) %>%
    transmute(rented = 1 - owned)
prent <- prent$rented
    
medval <- house_prices %>%
    filter(year==2001) %>%
    select(median) %>% 
    mutate(median=ifelse(median==0,68000, median))  # replace 0 values from no sales with national median value of 68000
medval <- medval$median


moving_costs_matrix <- generate_moving_costs(
    proportion_renting=prent,
    median_value=medval
)

###################################################################################################################################
# Input structure for Rcpp
###################################################################################################################################
oldcnts <- households %>% filter(year==2001) %>% select(all_households)
oldcnts <- oldcnts$all_households

newcnts <- households %>% filter(year==2011) %>% select(all_households)
newcnts <- newcnts$all_households

## Issue _ different number of datazones : 6505 for old and 6500 for new

In <- list(
    data = list(
        stayerprop= 0.82,
        counts = list(
            oldcounts=oldcnts ,
            newcounts=newcnts      
        )
    ),
    movingcostmatrix = moving_costs_matrix,
    utils = list(
        deltas = rep(0, length(oldcnts) + 1) # Added one for outside region
    ),
    params = list(
        tol_delta = 10^-2,
        tol_mu = 10^-1,
        mu_upper = 0.50,
        mu_lower = 0,
        maxit_outer = 1000,
        maxit_inner = 100000,
        use_logmu=F
    ),
    dbg = list(
        dbgflag= T,
        verboseflag=T
        
    )
)

rm(moving_costs_matrix)

x <- new(Contraction, In)


x$run()

x$extractests() -> outputs

rm(x)



####################################################################################################################################
##
###################################################################################################################################
# Input Data
###################################################################################################################################










# Example using US Data ---------------------------------------------------



# The code below has been run to separate the data in inputs_combined into a consistent format:

# data_messy <- read.csv("Data/inputs/us/inputs_combined.csv") %>% tbl_df()
# data_counts <- data_messy %>% select(-percent.rents, -medval)
# data_other <- data_messy %>% select(area, prop_rents=percent.rents, median_value=medval)
# 
# data_counts <- data_counts %>% 
#     gather(group_name, count, -area) 
# tmp <- data_counts$group_name %>% str_split("\\.") %>% ldply() 
# names(tmp) <- c("race", "extra")
# data_counts <- data_counts %>% bind_cols(tmp)
# rm(tmp)
# data_counts <- data_counts %>% 
#     tbl_df() %>%
#     mutate(income=str_sub(extra, 1L, 1L),
#            year=str_sub(extra, 2L)
#     ) %>% 
#     select(-group_name, -extra) %>% 
#     select(area, year, race, income, count)
# rm(data_messy)
# 
# write.csv(data_counts, file="Data/inputs/us/tidied_counts.csv")
# write.csv(data_other, file="Data/inputs/us/tidied_other_data.csv")

data_counts <- read.csv("Data/inputs/us/tidied_counts.csv") %>% tbl_df()
data_other <- read.csv("Data/inputs/us/tidied_other_data.csv") %>% tbl_df()

data_counts$year <- revalue(as.factor(data_counts$year), c("0"="t1", "9"="t2"))

# start with white, medium income
data_counts_by_year <- data_counts %>% 
    group_by(area, year) %>%
    filter(race=="white" & income=="m") %>%
    summarise(count=sum(count)) %>%
    spread(year, count)

correction  <- 0.0001

data_counts_by_year <- data_counts_by_year %>% 
    mutate(
        t1=t1 + correction,
        t2=t2 + correction
        )

stayer_proportions <- read.csv("Data/inputs/us/stayer_proportions.csv") %>% tbl_df()

stayer_proportions %>% 
    filter(race=="white" & income_level=="medium") %>%
    select(stayer_proportion) %>% as.numeric()
############################################################################################################################

# 
# moving_costs_matrix <- generate_moving_costs(
#     proportion_renting=data_other$prop_rents,
#     median_value=data_other$median_value
# )
# 
# ###################################################################################################################################
# # Input structure for Rcpp
# ###################################################################################################################################
# 
# 
# In <- list(
#     data = list(
#         stayerprop= 0.4174,
#         counts = list(
#             oldcounts=data_counts_by_year$t1 ,
#             newcounts=data_counts_by_year$t2      
#         )
#     ),
#     movingcostmatrix = moving_costs_matrix,
#     utils = list(
#         deltas = rep(0, length(data_counts_by_year$t1) + 1) # Added one for outside region
#     ),
#     params = list(
#         tol_delta = 10^-2,
#         tol_mu = 10^-1,
#         mu_upper = 0.50,
#         mu_lower = 0,
#         maxit_outer = 1000,
#         maxit_inner = 100000,
#         use_logmu=F
#     ),
#     dbg = list(
#         dbgflag= T,
#         verboseflag=T
#         
#     )
# )
# 
# rm(moving_costs_matrix)
# 
# x <- new(Contraction, In)
# 
# 
# x$run()
# 
# x$extractests() -> outputs
# 
# rm(x)
# 
# 
# 
# 
# 
# 
# setwd("~/Google Drive/PROJECTS/Sorting Models/SortingModelRcpp/")
# #load(file="Contraction_Mapping_Three_Stage__2014_04_16.rData")
# load("Results_White_lower.Rdata")
# 
# # Quantile regression in R
# 
# 
# 
# deltas_to_use <- outputs_final$deltas
# attach(Data012)
# 
# qreg <- rq(outputs_final$deltas ~ c(vcrime, 0), tau=0.5)
# 
# detach(Data012)
# 
# 
