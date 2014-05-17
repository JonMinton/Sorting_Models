Download_File_url <- function (
    url, 
    outfile,
    ..., sha1 = NULL) 
    {
        require(RCurl)
        require(devtools)
        require(repmis)
        require(httr)
        require(digest)
        
        stopifnot(is.character(url), length(url) == 1)
        filetag <- file(outfile, "wb")
        request <- GET(url)
        stop_for_status(request)
        writeBin(content(request, type = "raw"), filetag)
        close(filetag)
    }


