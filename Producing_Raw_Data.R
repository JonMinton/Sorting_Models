# Converting all files to plain text
require(foreign)


basedir <- "/Users/JonMinton/Dropbox/Data/"

# Files are : 

#GSPC
"GSPC/GUdata2014q1_UoGprocessed_a_I.dta"

#Nationwide
"2005_v1a.dta"            
"2006_v1a.dta"            
"2007_v1a.dta"            
"2007.dta"               
"2008.dta"                
"2009.dta"                
"2010.dta"                
"2011.dta"               
"2012.dta"                

"HPI data dictionary.xls" 

"NBSDS 2005.sav"          
"NBSDS 2006.sav"

#RoS 1990to2009
"RoS19902010_with_OS_geocodes.dta"

#RoS-Match
"data_sales_geocoded.dta"        
"postIDTOrosID_merged_1to21.zip" 
"readme-match_ROS_detail.txt"   
"SecondAddressMerge/ postIDTOrosID_merged.dta"             
"SecondAddressMerge.zip"    

###################################################
###################################################

#GSPC
Data <- read.dta(
    paste0(
        basedir, "GSPC/GUdata2014q1_UoGprocessed_a_I.dta"
        )
    )

write.table(
    Data,
    file=paste0(
        basedir,
        "raw/GSPC_Data_2014Q1.txt"
        )
    )

#Nationwide
"2005_v1a.dta"    
Data <- read.dta(
    paste0(
        basedir,
        "Nationwide/2005_v1a.dta"
        )
    )
write.table(
    Data,
    file=paste0(
        basedir,
        "raw/Nationwide_2005.txt"
        )
    )


"2006_v1a.dta"            
Data <- read.dta(
    paste0(
        basedir,
        "Nationwide/2006_v1a.dta"
    )
)
write.table(
    Data,
    file=paste0(
        basedir,
        "raw/Nationwide_2006.txt"
    )
)

"2007_v1a.dta"            
Data <- read.dta(
    paste0(
        basedir,
        "Nationwide/2007_v1a.dta"
    )
)
write.table(
    Data,
    file=paste0(
        basedir,
        "raw/Nationwide_2007_v1a.txt"
    )
)

"2007.dta"               
Data <- read.dta(
    paste0(
        basedir,
        "Nationwide/2007.dta"
    )
)
write.table(
    Data,
    file=paste0(
        basedir,
        "raw/Nationwide_2007.txt"
    )
)

"2008.dta"                
Data <- read.dta(
    paste0(
        basedir,
        "Nationwide/2008.dta"
    )
)
write.table(
    Data,
    file=paste0(
        basedir,
        "raw/Nationwide_2008.txt"
    )
)

"2009.dta"                
Data <- read.dta(
    paste0(
        basedir,
        "Nationwide/2009.dta"
    )
)
write.table(
    Data,
    file=paste0(
        basedir,
        "raw/Nationwide_2009.txt"
    )
)


"2010.dta"                
Data <- read.dta(
    paste0(
        basedir,
        "Nationwide/2010.dta"
    )
)
write.table(
    Data,
    file=paste0(
        basedir,
        "raw/Nationwide_2010.txt"
    )
)

"2011.dta"               
Data <- read.dta(
    paste0(
        basedir,
        "Nationwide/2011.dta"
    )
)
write.table(
    Data,
    file=paste0(
        basedir,
        "raw/Nationwide_2011.txt"
    )
)

"2012.dta"      

Data <- read.dta(
    paste0(
        basedir,
        "Nationwide/2012.dta"
    )
)
write.table(
    Data,
    file=paste0(
        basedir,
        "raw/Nationwide_2012.txt"
    )
)

"HPI data dictionary.xls" 

"NBSDS 2005.sav"     
Data <- read.spss(
    paste0(
        basedir,
        "Nationwide/NBSDS 2005.sav"
    ),
    to.data.frame=T
)
write.table(
    Data,
    file=paste0(
        basedir,
        "raw/Nationwide_NBSDS_2005.txt"
        )
    )

"NBSDS 2006.sav"
Data <- read.spss(
    paste0(
        basedir,
        "Nationwide/NBSDS 2006.sav"
    ),
    to.data.frame=T
)
write.table(
    Data,
    file=paste0(
        basedir,
        "raw/Nationwide_NBSDS_2006.txt"
    )
)

#RoS 1990to2009
"RoS19902010_with_OS_geocodes.dta"
Data <- read.dta(
    paste0(
        basedir,
        "RoS 1990to2009/RoS19902010_with_OS_geocodes.dta"
        )
    )
write.table(
    Data,
    file=paste0(
        basedir,
        "raw/RoS_with_OS_geocodes_1990_to_2010.txt"
    )
)


#RoS-Match
"data_sales_geocoded.dta"        
Data <- read.dta(
    paste0(
        basedir,
        "RoS-Match/data_sales_geocoded.dta"
        )
    )
# This will not load on the mac:
# Error in read.dta(paste0(basedir, "RoS-Match/data_sales_geocoded.dta")) : 
#     cannot have attributes on a CHARSXP


write.table(
    Data,
    file=paste0(
        basedir,
        "raw/RoS_Data_Sales_Geocoded.txt"
        )
    )

"postIDTOrosID_merged_1to21.zip" 
"readme-match_ROS_detail.txt"   
"SecondAddressMerge/ postIDTOrosID_merged.dta"             
"SecondAddressMerge.zip"    

Data <- read.dta(
    paste0(
        basedir,
        "RoS-Match/SecondAddressMerge/postIDTOrosID_merged.dta"
        )
    )

write.table(
    Data,
    
    file=paste0(
        basedir,
        "raw/RoS_Match__Second_Address_Merge__Post_IDTO_RoS_ID_Merged.txt"
        )
    )

# After some unpacking
Data <- read.dta(
    paste0(
        basedir,
        "RoS-Match/PostIDTOrosID_merged_1to21.dta"
        )
    )

write.table(
    Data,
    file=paste0(
        basedir,
        "raw/RoS_Match__Post_IDTO_RoS_ID_Merged__01_to_21.txt"
        )
    )
