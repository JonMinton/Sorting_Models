### calculate a house moving cost matrix. THERE are a lot of factor to consider: physical cost mainly the house removal
### cost and financial costs including mainly tax and fees. 1. Financial costs mainly include (from rightmove.com) Stamp duty
### fee; Estate agent fee and legal fees; 2. Physical cost is mainly influenced by the weight of shipment and the distance moved

### This matrix is annualised cost of house moving, see the white flight paper.

### First calculate the financial cost adjusted by house ownership rate in each datazone.



### Read in the datazone shapefile data
require(spdep)
require(maptools)
require(plyr)
require(dplyr)
require(tidyr)
require(sp)
require(rgeos)

# Import the Scottish datazone shapefile 
datazone_shp <- readShapePoly("moving_cost_estimation/shapefiles/datazone_correct.shp")
row.names(datazone_shp) <- as.character(datazone_shp$datazone)
# assign the correct projection
proj4string(datazone_shp) <- CRS("+init=epsg:27700")

# import the house price data
houseprice <- read.csv("moving_cost_estimation/csv_files/house_sales_and_prices.csv",
                       header=TRUE) %>% tbl_df()
# for estimating a robust median house price for year 2001, we use averages of house prices transacted between 1999 and 2003
houseprice <- houseprice %>%
  select(datazone,HO.hpricemedian,year) %>%
  filter(year >= 1999 & year <= 2003)
houseprice <- houseprice %>% group_by(datazone) %>% summarise(meanprice=mean(HO.hpricemedian))
# there are 21 datazones with house price equal 0
houseprice %>% filter(meanprice==0)

# where are they?
row.names(houseprice) <- as.character(houseprice$datazone)
houseprice <- data.frame(houseprice)
# Simply assme the mean house price of the whole scotland to these areas. Could always change this.
houseprice$meanprice[houseprice$meanprice == 0] <- mean(houseprice$meanprice) 
summary(houseprice$meanprice)
plot(density(houseprice$meanprice))

# home ownership data
tenure <- read.csv("moving_cost_estimation/csv_files/tenure_households.csv",
                       header=TRUE) %>% tbl_df()
tenure <- as.numeric(tenure$HO.percown/100)

# spaitalpolygondataframe
datazone.df <- SpatialPolygonsDataFrame(datazone_shp,houseprice)
# extract the coordinates of each datazone polygon
datazones_centroids <- gCentroid(datazone.df,byid=TRUE)
coords <- datazones_centroids@coords/1000
#### A house moving cost matrix
# stamp duty fee in 2001; 
# <= 60,000   0%
# 60,000 ~ 250,000 1%
# 250,000 ~ 500,000 3%
# >= 500,000 4%
# this is different from the current rates of stamp duty fee

housemovecostMatrix <- function(ownershipVector,
                                coords,# coordinates of each datazone used for weights of physical costs
                                w.min=0.1,
                                w.max=3.6, # these two speficy the range of weights adjusted for removal cost
                                agentfeeRate=0.018,#agent fee rate from rightnove
                                othercostRate=0.005,# rightmove uses a fixed rate of house value being sold
                                legalfeeRate=0.005,
                                interestRate=0.0232,# set real interest rate, default mean of real interest rate 2000-2004
                                discountYear=25,# set longth of discount
                                housepriceDF# the house price data
                                ){
  # the output is a house move cost matrix, 6505*6505
  n <- dim(coords)[1]#6505
  costMatrix <- matrix(nrow=n,ncol=n)
  # note, it is not symmetric
  # 
  houseprice <- housepriceDF
  # link with the stamp duty fee rate
  houseprice$stamp.rate[houseprice$meanprice <= 60000] <- 0
  houseprice$stamp.rate[houseprice$meanprice > 60000 & houseprice$meanprice <= 250000] <- 0.01
  houseprice$stamp.rate[houseprice$meanprice > 250000 & houseprice$meanprice <= 500000] <- 0.03
  houseprice$stamp.rate[houseprice$meanprice > 500000] <- 0.04
  # proportion of house ownership
  ownership <- ownershipVector
  
  meanprice <- houseprice$meanprice
  stamp.rate <- houseprice$stamp.rate
  
  # set up the distance weighted physical moving costs, note later in the loop there is further adjustment
  # by the proportion of home ownership
  # In the function, the othercostRate is an average rate of house removal cost (see rightmove house moving cost
  # caculator). we consider this rate as the removal cost when the moving distance is the overall average of  
  # the distance matrix among Scottish datazones. so we need to adjust this average rate by the ratio of real distance
  # between two datazones to the average distance
  
  # distance matrix
  dist.M <- spDists(coords)
  # weights associated with distance matrix 
  mean.distance <- median(c(dist.M))
  # weights matrix
  weights.M <- dist.M/mean.distance
  min.temp <- min(c(weights.M))
  max.temp <- max(c(weights.M))
  # the elememts of this matrix ranges from 0 to about 9.
  # Note the removal cost at the average distance is about 0.5% and the agent fee rate is about 1.8%.
  # we donot want the romoval cost larger than the agent fee rate. so we linearly rescale the weights matrix
  # to a range of [0.1,3.6]. note 3 is about 1.8% / 0.5%. therefore, we assume the largest removal cost is as high as
  # the estate agent fee. However, we can change the ratio when conducting sensitive analyese.
  
  # the linearly rescale the weights matrix to a range of [0.1,3.6]
  aa <- (w.max - w.min)/(max.temp - min.temp)
  bb <- w.max - aa*max.temp
  weights.M <- aa*weights.M + bb
  
  ### start the calculation
  
  for(i in 1:n) {
    for(j in 1:n) {
      # estate agent fee only apply to house sellers. NOTE it is adjusted by house ownership rate in area i
      agentfee.temp <- meanprice[i]*agentfeeRate*ownershipVector[i]
      # removal cost, see calculation method in rightmove. we also consider the distance
      # note even for renters who want to move, they also incur removal cost although the cost should be smaller than
      # house owners. We hope by adjusting the home ownership, this bias could become smaller.
      removal.temp <- meanprice[i]*othercostRate*ownershipVector[i]*weights.M[i,j]
      #legal fee applies to the value of house to buy
      legal.temp <- meanprice[j]*legalfeeRate*ownershipVector[j]
      # stamp duty fee applies to the value of house to buy
      stamp.temp <- meanprice[j]*stamp.rate[j]*ownershipVector[j]
      
      ### the sum
      sum.temp <- agentfee.temp + removal.temp + legal.temp + stamp.temp
      ### calculate the annualised moving cost using UK real interest rates, see the white flight paper, 
      inv.rates <- 1/(interestRate + 1)
      temp <- (1-inv.rates^discountYear)/(1-inv.rates)
      annual.temp <- sum.temp/temp
      
      costMatrix[i,j] <- annual.temp
    }
    if(i%%50 == 0) cat(".")
  }
  return(costMatrix)
} 


 costMatrix <-  housemovecostMatrix(ownershipVector=tenure,
                            coords=coords,housepriceDF=houseprice)
  
# The moving cost of moving to or from a single "catch-all" location N+1
# following the whiteflight paper, we assume the cost of moving outside of Scotland is the larget element of costMatrix
# The dimension of the house moving cost matrix is therefore 6506 (6505+1) * 6506
summary(c(costMatrix))
rm(list=setdiff(ls(),"costMatrix"))
gc()
max.cost <- max(c(costMatrix))
moving.cost <- matrix(nrow=6506,ncol=6506) 
moving.cost[1:6505,1:6505] <- costMatrix
moving.cost[6506,] <- max.cost
moving.cost[,6506] <- max.cost

save(list="moving.cost",file="housemovingcost.RData")
### just load("housemovingcost.RData") 