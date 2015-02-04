###################################################################################################################################
# R Functions
###################################################################################################################################

generate_moving_costs <- function(
    proportion_renting,
    median_value,
    vend_cost=0.03,
    movecost_numerator=2910, # moving cost numerator
    movecost_denominator=24.556, # moving cost denominator
    movecost_outside=528.71 # cost of moving out of country      
){
    
    # Error checking
    if (length(proportion_renting)!=length(medval)) stop("proprents and medval of different lengths")
    N <- length(proportion_renting)
    
    out <- matrix(nrow=N+1, ncol=N+1)
    
    for (i in 1:(N+1)){
        for (j in 1:(N+1)){      
            # If not moving
            if (i==j){
                out[i,j] <- 0
            } 
            # If not moving out of the county
            if ((i != j) & (i != (N+1)) & (j != (N+1))){
                out[i,j] <- movecost_numerator + vend_cost *(
                    (1 - proportion_renting[i]) * median_value[i] + (1 - proportion_renting[j]) * median_value[j]
                ) / movecost_denominator
            }
            # If moving out of the county
            if ((i != j) & ((i==(N+1)) | (j==(N+1)))){
                out[i,j] <- movecost_outside    
            }  
            
        }
    }
    return (out) 
}



generate_moving_costs <- cmpfun(generate_moving_costs)
