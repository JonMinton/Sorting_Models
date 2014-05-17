# Pollution data
# Source 
# http://uk-air.defra.gov.uk/data/pcm-data#ozone


require("xlsx")
require("ggplot2")

source("Scripts/Secure_File_Downloader.R")
# link to Excel workbook providing links:

Links_Workbook_Link <- paste0(
    "https://www.dropbox.com/s/",
    "l1q8jpzclfm414l/",
    "Links_To_Modelled_Background_Pollution_Data.xlsx"
    )

Download_File_url(
    url=Links_Workbook_Link,
    outfile="Links_To_Modelled_Background_Pollution_Data.xlsx"
    )

Links_DF <- read.xlsx(
    file="Links_To_Modelled_Background_Pollution_Data.xlsx",
    sheetName="Main"
                      )
Links_DF$Link <- as.character(Links_DF$Link)
Links_DF$Comments <- as.character(Links_DF$Comments)

#Want to be able to read in a csv file directly


N.files <- dim(Links_DF)[1]


for (i in 1:N.files){
    
    
    DF_from_csv <- source_data(
        url=Links_DF$Link[i],
        skip=5,
        na.strings="MISSING"
    )
    filename <- paste0(
        "Figures/Pollution/",
        DF_from_csv$Header.Label[i],
        ".png"
    )
    
    p <- qplot(x=x, y=y, 
               colour=eval(parse(text=(names(DF_from_csv)[4]))), 
               data=DF_from_csv
    ) 
    p <- p + labs(x="Easting", y="Northing", 
                  title=names(DF_from_csv[4])
    ) 
    p <- p + scale_colour_gradient(
        limits=c(0,30), low="white", high="blue",
        guide=guide_legend(
            title=NULL
        )
    )
    p <- p + theme_minimal()
    
    
    print(p)
    browser()
#     ggsave(file=filename,
#            width=900,
#            height=1000
#            )
    
    
}

