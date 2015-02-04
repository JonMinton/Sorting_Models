# tidier script for running contraction mapping algorithm as part of Timmins' structural sorting models

# Jon Minton
# 4/2/2015


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

# by individual or by household?

persons <- source_DropboxData(
    file="persons.csv",
    key="vcz7qngb44vbynq"
    ) %>% 
    tbl_df() %>%
    select(datazone, year,
           contains("hspeop")
    )

names(persons) <- names(persons) %>%str_replace_all( "GR.hspeop", "count_both_")
persons <- persons %>% gather(age_group, count, -datazone, -year)

persons <- persons %>% mutate(gender="both", age_group=str_replace_all(age_group, "count_both_", ""))
persons$age_group <- persons$age_group %>% revalue(
    c(
        "1619" = "16_19",
        "2024" = "20_24",
        "2529" = "25_29", 
        "3034" = "30_34",
        "3539" = "35_39",
        "4044" = "40_44",
        "4549" = "45_49",
        "5054" = "50_54",
        "5559" = "55_59",
        "6064" = "60_64",
        "6569" = "65_69",
        "7074" = "70_74",
        "7579" = "75_79",
        "8084" = "80_84",
        "85over" = "85_100"
    )
)

# deal with "" separately as revalue can't cope
persons$age_group[nchar(persons$age_group)==0] <- "all"

fn <- function(x){
    tmp <- str_split(x, "_") %>% ldply()
    if (length(tmp)==2){
        lower <- tmp[1] %>% as.numeric()
        upper <- tmp[2] %>% as.numeric()
        if (!is.na(lower) & !is.na(upper)){
            out <- data.frame(lower=lower, upper=upper)
        } else {out <- data.frame(lower=NA, upper=NA)}
    } else {out <- data.frame(lower=NA, upper=NA)}
    return(out)
}

lower_upper <- ldply(persons$age_group, fn)

# Stayer probabilities: 
# First pass: 
#  - the median age in Scotland in 2001 was 38
#  - the estimated 10 year stayer proportion for this group is 0.82





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
