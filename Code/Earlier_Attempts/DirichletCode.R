# For producing simulated area counts given the available count data, but also incorporating
# a noninformative prior of 0.5 to all cell counts

# Function Arguments:
# ===================
# InputData : Data frame containing cell counts only. N observations deep
# M : number of hypothetical datasets to generate
# prior: continuity correction to add to all cell counts
# as.freqs: flags whether the output should give frequencies (T) or row-wise proportions (F)


# Function Output:
# A list containing M dataframes.
# Each dataframe should contain a hypothetical dataset of the same dimensions
# as the real dataset. OLS and other processes could then be applied to each 
# of these datasets, to produce a bootstrapped family of quantities of interest 
# (e.g. distribution of medians, means, upper quartiles etc)

MakeDirichletFreqs <- function(InputData, M, prior=0.5, as.freqs=F){
  # the following line checks whether a package called MCMCPack has already been loaded, 
  # and if it has not then loads it. (MCMCpack contains the rdirichlet function used later)
  require(MCMCpack)
  
  # The next command calculates the number of rows the dataframe object InputData has.
  # It passes InputData to the dim() function. The output of the dim function is a 
  # vector whose length is equal to the number of dimensions of the object passed to it
  # A standard dataframe has two dimensions. The first dimension is the number of rows
  # The second dimension is the number of columns
  # The [1] suffix extracts the first element from the output of dim(InputData), namely the number 
  # of rows! 
  N <- dim(InputData)[1]  
  
  # This creates a list object of length N
  tmp <- vector("list", N)

  # If the output from the function should be frequencies rather than proportions,
  # then calculate the sum of each row, add the prior to each row,
  # then save the output to a vector object called area.totals
  # This uses the apply function
  #   the arguments of an apply function are:
  #    1) InputData: Object to be acted upon
  #    2) 1: Dimension of object to be acted over (1 for 'by row', 2 for 'by column')
  #    3) sum: function to perform on each row (1) or column (2) of the object acted upon
  # So, the line says: "For each row, calculate the sum"
  if (as.freqs==T){
    area.totals <- apply(InputData,1, sum) + prior
  }
  
  # by default, the input to rdirichlet is a vector of the frequencies in each area
  # the output from rdirichlet is a dataframe where each row is a different estimate 
  # for that area. 
  # This will require rearranging later
  
  # n.b. the double square brackets [[i]] are used to access positions in list objects
  # as opposed to [i,] : this returns the whole of the ith row from the dataframe InputData
  # [i,j] would return the element in the ith row and jth column within InputData
  # [,j] would return the element in the jth column from InputData
  # The t() function means transpose, which is needed to 'flip' InputData[i,] 
  # into a K x 1 vector which rdirichlet expects as its second argument
  # rather than the 1 x K matrix which InputData[i,]  returns by default
  
  for (i in 1:N){
    tmp[[i]] <- rdirichlet(M, t(InputData[i,] + prior) )
  }
  
  # This code rearranges the output tmp from above
  # so that the final output from the function is
  # M sets of estimated frequencies, where each row 
  # indicates a different area
  
  # The command Output[[i]] <- t(sapply(tmp, function(x) x[i,] ))
  # may need some unpacking. It's equivalent to:
  
  # Output[[i]] <-  t(                        
  #                   sapply(                 
  #                     tmp,                  
  #                     function(x){          
  #                                           
  #                       x[i,]             
  #                     }
  #                   )
  #                 )

  # However I'll skip over trying to describe it for now...
  
  Output <- vector("list", M)
  for (i in 1:M){
    Output[[i]] <- t(sapply(tmp, function(x) x[i,] ))
    if (as.freqs==T){
      Output[[i]] <- Output[[i]] * area.totals
    }
  }
  
  return(Output)
}
