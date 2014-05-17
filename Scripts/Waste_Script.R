#Code to to do some exploratory analysis of waste siting by type

# Prerequisites
require("xlsx")
require("ggplot2")
require("repmis")
require("devtools")
require("RCurl")
# Downloading data from Dropbox location

# Whitespace trimmer
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

FlagFlipper <- function(
    tokens, 
    aliases
    ){
    out <- if (any(aliases %in% tokens)) {T} else {F}
    return(out)
}

data_key <- "ka8d94o91sekwrf/"
data_filename <- "List%20of%20waste%20sites%20%26%20capacity%20-%202012.xls"
    
link <- paste0(
        "https://www.dropbox.com/s/",
        "ka8d94o91sekwrf/",
        "List%20of%20waste%20sites%20%26%20capacity%20-%202012.xls"
    )

url <- getURL(link)


temp <- tempfile()
download.file(url, temp)

file.create("temp.xlsx")
download.file(textConnection(tmp), "temp.xlsx")
Data_from_Excel <- read.xlsx(
    tmp, 
    sheetName="Sheet1")



#install.packages("repmis")


# Giving up on this for now, just looking at the data from a local location

Data_from_Excel <- read.xlsx(
    "Data/Excel_Files//List of waste sites & capacity - 2012.xls",
    sheetName="All Sites - 2012", stringsAsFactors=F
    )

# Some data cleaning needed, but not too much
tmp <- Data_from_Excel

tmp <- tmp[,1:36]
N.vars <- length(names(tmp))

# want to identify categories and supercategories
# supercategories change whenever they are not NA

row1 <- names(tmp)
row2 <- tmp[1, , drop=T]
row1.blanks <- grep("^NA", row1)
row2.blanks <- which(is.na(row2))
row3 <- rep("", N.vars)

for (i in 1:N.vars){
    if (i %in% row1.blanks){ row1[i] <- row1[i-1]}
    if (i %in% row2.blanks){ 
        row3[i] <- row1[i]
    } else {
        row3[i] <- paste(row1[i], row2[i], sep="_")
    }   
}

Data_Cleaned <- tmp[-1,]
names(Data_Cleaned) <- row3


waste.types <- Data_Cleaned[,"Licensed.Waste.Types_Waste Type"]
#Want to separate by "/" and remove leading and trailing white spaces

tokens <- strsplit(waste.types, "/")
tokens <- lapply(tokens, trim)
# go through each list element and look for unique tokens,
# if a unique token exists then add it to the lists of tokens to search for

unique.tokens <- c()

N.obs <- length(tokens)

for (i in 1:N.obs){
    this.row <- tokens[[i]]
    N.tokens <- length(this.row)
    
    for (j in 1:N.tokens){
        this.token <- this.row[j]
        if (!(this.token %in% unique.tokens)){
            unique.tokens <- c(unique.tokens, this.token)
        }
    }
}

is.Household <- rep(NA, N.obs)
is.Commercial <- rep(NA, N.obs)
is.Industrial <- rep(NA, N.obs)
is.Inert <- rep(NA, N.obs)
is.Asbestos <- rep(NA, N.obs)
is.Special <- rep(NA, N.obs)

Household.aliases <- c(
    "Household",
    "Househol"
    )
Commercial.aliases <- c(
    "Commercial", 
    "Commercal"
    )
Industrial.aliases <- c(
    "Industrial"
    )
Inert.aliases <- c(
    "Inert"
    )
Asbestos.aliases <- c(
    "Special asbestos",
    "Special Asbestos"
    ) # Should 'Special' be part of this or is it a distinct category?
Special.aliases <- c(
    "Special"
    )


for (i in 1:N.obs){
    this.tokens <- tokens[[i]]
    
    is.Household[i] <- FlagFlipper(this.tokens, Household.aliases)
    is.Commercial[i] <- FlagFlipper(this.tokens, Commercial.aliases)
    is.Industrial[i] <- FlagFlipper(this.tokens, Industrial.aliases)
    is.Inert[i] <- FlagFlipper(this.tokens, Inert.aliases)
    is.Asbestos[i] <- FlagFlipper(this.tokens, Asbestos.aliases)
    is.Special[i] <- FlagFlipper(this.tokens, Special.aliases)    
    
}

table(is.Household)
table(is.Commercial)
table(is.Industrial)
table(is.Inert)
table(is.Asbestos)
table(is.Special)

# Capacity per annum

tmp1 <- as.numeric(Data_Cleaned[,23])
tmp2 <- as.character(Data_Cleaned[,24])
tmp3 <- ifelse(tmp2=="annual" | tmp2=="Annual", 
               tmp1,
               ifelse(
                   tmp2=="Daily",
                   tmp1*365,
                   NA
                   )
               )


DF_tmp <- data.frame(
    size=tmp3,
    x=Data_Cleaned[,6],
    y=Data_Cleaned[,7]
    )



qplot(x=x,y=y, size=size, data=DF_tmp)
p1 <- ggplot(x=x, y=y, size=size, data=DF_tmp)
p1 + geom_point()
?geom_point
p1 + geom_point(aes(x=x, y=y, size=size))
p1 + geom_point(aes(x=x, y=y, size=size, colour=size))
