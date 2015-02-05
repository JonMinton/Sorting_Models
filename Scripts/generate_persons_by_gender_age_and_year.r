
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

persons_by_age <- persons %>% filter(grepl("_", age_group)) # this is how to filter by contents of age_group

persons_by_age <- persons_by_age %>% 
    group_by(age_group) %>% 
    mutate(
        lower_age = str_split(age_group, "_")[[1]][1] %>% as.numeric(),
        upper_age = str_split(age_group, "_")[[1]][2] %>% as.numeric()
    ) # slow version

# Faster version, suggested by
# http://stackoverflow.com/questions/28342012/how-to-make-creating-upper-and-lower-age-variables-from-existing-variable-much-f#28343290
   # replace the mutate function with:
# do(cbind(.,matrix(rep(unlist(strsplit(as.character(.[1,3]), "_")),nrow(.)),ncol=2,byrow=TRUE)))
 #  However additional variables are not named. 

# add age 0 

persons_by_age <- persons %>% 
    filter(age_group==0) %>%
    mutate(
        lower_age=0,
        upper_age=0) %>% 
    bind_rows(persons_by_age)


# To do the same with males and females


males <- source_DropboxData(
    file="males.csv",
    key="n77e4r376gz368e"
) %>% 
    tbl_df() %>%
    select(datazone, year,
           contains("hsmal")
    )

names(males) <- names(males) %>%str_replace_all( "GR.hsmal", "count_males_")
males <- males %>% gather(age_group, count, -datazone, -year)

males <- males %>% mutate(gender="male", age_group=str_replace_all(age_group, "count_males_", ""))
males$age_group <- males$age_group %>% revalue(
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
males$age_group[nchar(males$age_group)==0] <- "all"

males_by_age <- males %>% filter(grepl("_", age_group)) # this is how to filter by contents of age_group

males_by_age <- males_by_age %>% 
    group_by(age_group) %>% 
    mutate(
        lower_age = str_split(age_group, "_")[[1]][1] %>% as.numeric(),
        upper_age = str_split(age_group, "_")[[1]][2] %>% as.numeric()
    )

# add age 0 

males_by_age <- persons %>% 
    filter(age_group==0) %>%
    mutate(
        lower_age=0,
        upper_age=0) %>% 
    bind_rows(males_by_age)

persons_by_age <- persons_by_age %>% bind_rows(males_by_age)



females <- source_DropboxData(
    file="females.csv",
    key="qq083qk9iah5txz"
) %>% 
    tbl_df() %>%
    select(datazone, year,
           contains("hsfem")
    )

names(females) <- names(females) %>%str_replace_all( "GR.hsfem", "count_females_")
females <- females %>% gather(age_group, count, -datazone, -year)

females <- females %>% mutate(gender="female", age_group=str_replace_all(age_group, "count_females_", ""))
females$age_group <- females$age_group %>% revalue(
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
females$age_group[nchar(females$age_group)==0] <- "all"

females_by_age <- females %>% filter(grepl("_", age_group)) # this is how to filter by contents of age_group

females_by_age <- females_by_age %>% 
    group_by(age_group) %>% 
    mutate(
        lower_age = str_split(age_group, "_")[[1]][1] %>% as.numeric(),
        upper_age = str_split(age_group, "_")[[1]][2] %>% as.numeric()
    )

# add age 0 

females_by_age <- persons %>% 
    filter(age_group==0) %>%
    mutate(
        lower_age=0,
        upper_age=0) %>% 
    bind_rows(females_by_age)

persons_by_age <- persons_by_age %>% bind_rows(females_by_age)


# now to write this as a csv file so I don't have to do the coding each time
#as the file is too large for github I'm saving to dropbox instead

write.csv(persons_by_age, "E:/Dropbox/Data/SNS/derived_data/persons_by_gender_year_and_age.csv", row.names=F)
