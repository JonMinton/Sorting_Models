#Requirements 

require("foreign")
require("RCurl")
require("gdata")

# RoS GSPC Merge

https://www.dropbox.com/s/qycuvbout0mder3/GUdata2014q1_UoGprocessed_a_I.dta

Dropbox_Prefix <- "https://www.dropbox.com/s"

Dropbox_Data <- list(
    GSPC=list(
        GUData2014q1_UoGprocessed_a_1=list(
            key="qycuvbout0mder3",
            file="GUdata2014q1_UoGprocessed_a_I.dta"
        )
    ),
    Nationwide=list(
        F_2005_v1a=list(
            key="fb3ho467fywjxb2",
            file="2005_v1a.dta"
        ),
        F_2006_v1a=list(
            key="6x2w1a9rgzg3sos",
            file="2006_v1a.dta"
        ),
        F_2007_v1a=list(
            key="jrhb7p3to2rf7rs",
            file="2007_v1a.dta"
        ),
        F_2007=list(
            key="nyhpr6x08ngwu31",
            file="2007.dta"
        ),
        F_2008=list(
            key="j6nv1c7h9h07lhh",
            file="2008.dta"
        ),
        F_2009=list(
            key="kf1hwu2oup6qyhx",
            file="2009.dta"
        ),
        F_2010=list(
            key="00sa1z3tq7fb9wv",
            file="2010.dta"
        ),
        F_2011=list(
            key="vylisf5s2x58151",
            file="2011.dta"
        ),
        F_2012=list(
            key="0vy1jp97s309zwm",
            file="2012.dta"
        ),
        HPI_Data_Dictionary=list(
            key="1y4wqcjfue8zqvl",
            file="HPI%20data%20dictionary.xls"
        ),
        NBSDS_2005=list(
            key="jwqo0gwpkljmp4u",
            file="NBSDS%202005.sav"
        ),
        NBSDS_2006=list(
            key="iv4795udhxyf4rj",
            file="NBSDS%202006.sav"
        )
    ),
    RoS_1990_2010=list(
            key="v3kdcguqbza5ggv",
            file="RoS19902010_with_OS_geocodes.dta"
    ),
    RoS_Match=list(
        Data_Sales_Geocoded=list(
            key="mcia30b91nlpr0c",
            file="data_sales_geocoded.dta"
        ),
        Post_Merged=list(
            key="wjrrxd569q7ey5s",
            file="postIDTOrosID_merged_1to21.zip"
        ),
        Readme_Match=list(
            key="9p3yse4g3egyl7w",
            file="readme-match_ROS_detail.txt"
        ),
        Second_Address_Merge=list(
            key="fg3cbz7pvfkcnfs",
            file="SecondAddressMerge.zip"
        )
    )
)


link <- paste("https://dl.dropboxusercontent.com/s",
              "{THE SHA-1 KEY}",
              "{THE FILE NAME}",
              sep="/")


tmp <- tempfile()
download.file(link, 
              destfile="ExcelFile.xlsx", method="curl",
              userpwd = "jonathan.minton@glasgow.ac.uk:Foxhead00"
              )

#read.xls(tmp, skip=2)
unlink(tmp)


download.file(url=getURL(url),destfile=Dropbox_Data[[2]][[1]]$file)
closeAllConnections()


# Given up (again...)

# Let's just look at one of the files
require("foreign")

GSPC_Data <- read.dta(
    file="/Users/JonMinton/Dropbox/Data/GSPC/GUdata2014q1_UoGprocessed_a_I.dta"
    )
