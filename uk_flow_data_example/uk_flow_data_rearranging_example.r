# Script for working with 2011 UK census flow data table

rm(list=ls())
gc()

require(tidyr)
require(dplyr)
require(stringr)

input_table <- read.csv("uk_flow_data_example/input/bulk.csv", check.names=F) %>% tbl_df()

#check.names means that characters are preserved in column names; the ';' character will be needed later

input_table <- input_table %>% 
    gather(key="category", value="count", -c(1:3)) %>%
    rename(year=date)

input_table$num_semicolons <- 
    str_count(input_table$category, ";")

input_table %>% select(num_semicolons) %>% table()
# so, most have three semi-colons; a smaller proportion have four


# split input_table into three and four semicolon tables

input_table_3 <- input_table %>% filter(num_semicolons==3)
    
input_table_4 <- input_table %>% filter(num_semicolons==4)

input_table_3 <- input_table_3 %>% separate(
    category,
    sep=";",
    into=c(
        "ns_sec",
        "ethnicity",
        "migration",
        "measure"
    )
)

    
input_table_4 <- input_table_4 %>% separate(
    category,
    sep=";",
    into=c(
        "ns_sec",
        "ethnicity", 
        "migration",
        "migration_2",
        "measure"
        )
    )

# now need to merge back migration and migration_2

input_table_4 <- input_table_4 %>% unite(
    migration_united,
    migration, migration_2,
    sep=""
    ) %>% rename(migration=migration_united)
names(input_table_3)
names(input_table_4)

input_table <- input_table_3 %>% bind_rows(input_table_4)
rm(input_table_3, input_table_4)

input_table <- input_table %>% select(-measure, -num_semicolons)

# now want to neaten category names in ns_sec, ethnicity and migration columns

input_table$ns_sec <- input_table$ns_sec %>% 
    str_replace_all("NS-SeC|[\\:]", "") %>% 
    str_trim()

input_table$ethnicity <- input_table$ethnicity %>% 
    str_replace_all("Ethnic [gG]roup|[\\:]", "") %>% 
    str_trim()

input_table$migration <- input_table$migration %>% 
    str_replace_all("Migration|[\\:]", "") %>% 
    str_trim() 


names(input_table)[3] <- "areal_unit"

# Now to save

write.csv(input_table, file="uk_flow_data_example/output/tidied.csv", row.names=FALSE)
