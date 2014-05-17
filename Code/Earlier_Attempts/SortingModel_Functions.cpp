#include <Rcpp.h>
#include <iostream>
#include <vector>
#include <cstdio> 

using namespace Rcpp;

// STRUCTURE OF INPUT OBJECT TO DeltaContraction Cl
//In <- list(
//  counts = list(
//    oldcounts=countsold,
//    newcounts=countsnew
//  ),
//  movingcostmatrix = movingcostmatrix,
//  utils = list(
//    mu = 0.5,
//    deltas = rep(0, length(countsold))
//  ),
//  dbg = list(
//    dbgflag= T,
//    verboseflag=T
//    )  
//)



// What do I want 

// 1) A function which returns an  n by n double matrix
// n is the number of areal units
// the row is the from 
// the col is the to

// [[Rcpp::export]] 
Class DeltaContraction{
  
public:
  DeltaContraction (SEXP Input_);  // Constructor function

  SEXP ReturnDeltas(); 
  
private:

  Rcpp::IntegerVector oldcounts, newcounts; 
  Rcpp::NumericMatrix movingcostmatrix; // Unsure whether this can be dynamically allocated
  double mu; 
  bool dbgflag, verboseflag;
  vector<double> deltas; // Ideally, this should be a double vector as it is where the most loading and unloading will take place
  
}


// MakeMovingCostMatrix needs to resize the matrix movingcostmatrix according to the 
// size of the data passed to it
// Use the following as an examplar: 
// http://stackoverflow.com/questions/1403150/how-do-you-dynamically-allocate-a-matrix


DeltaContraction::DeltaContraction (SEXP Input_)
{
    Rcpp::create::List                     Input(Input_)                  ; 
    Rcpp::List            utils           (Input["utils"])                ;
    Rcpp::List            counts          (Input["counts"])               ;      
    Rcpp::List            dbg             (Input["dbg"])                  ;

    // Utils list subset
    mu  =   Rcpp::as<double>              (utils["mu"])                   ;
    Rcpp::NumericVector   deltas_nv       (utils["deltas"] )              ;
    
    deltas.resize(deltas_nv.size())                                       ;
    for (int i = 0; i < deltas.size(); i++)
      deltas.push_back() = deltas_nv(i)                                   ;
    
    // counts list subset
    oldcounts = as<IntegerVector>         (counts["oldcounts"])           ;
    newcounts = as<IntegerVector>         (counts["newcounts"])           ; 
    
    // dbg list subset
    dbgflag     = as<bool>                (dbg["dbgflag"])                ;
    verboseflag = as<bool>                (dbg["verboseflag"])            ;
    
    // movingcostmatrix
    
    movingcostmatrix = as<NumericVector>  (Input["movingcostmatrix"])     ;
    
}


// double check how to return list as SEXP
SEXP DeltaContraction::ReturnDeltas(){
  
  // Returns whatever the output is
  
  return Rcpp::List::Create(
      Rcpp::named("newdeltas")  = deltas;
      // other objects to return when needed
    ); 
}


  double CalcDenom(
    int areak,
    double mu,
    NumericVector deltas,
    NumericVector percentRents,
    NumericVector medval
  ){    
    double output = 0;
    areak = areak - 1;
    int N = deltas.size() ;

    for (int i = 0; i < N; ++i){
      double tmp2 = CalcMovingCost(
        areak + 1,
        i + 1,
        percentRents, 
        medval);
      output += exp(deltas(i) - deltas(areak) - mu * tmp2);
    }

    return output;
  } 
  
  
double CalcForceAttract(
  int fromArea,
  int toArea,
  double mu,
  NumericVector deltas,
  NumericVector percentRents,
  NumericVector medval
  ){
    toArea = toArea - 1 ;
    fromArea = fromArea - 1 ;
    
    double output ;
    
    double tmp = CalcMovingCost(
        fromArea + 1,
        toArea + 1,
        percentRents, 
        medval);
    output = exp(deltas(toArea) - deltas(fromArea) - mu * tmp) ; 
    
    return output ;
  } 
  
  double CalcPredictedShare(
    int fromArea,
    int toArea,
    double mu,
    NumericVector deltas,
    NumericVector percentRents,
    NumericVector medval
    ){
      
      
      double denom = CalcDenom(
          fromArea, 
          mu, 
          deltas, 
          percentRents, 
          medval
          ) ;

      double forceAttract = CalcForceAttract(
          fromArea,
          toArea,
          mu,
          deltas,
          percentRents,
          medval
          );
      
    
      double output = forceAttract / denom ; 
      return output ; 
      
    }

Rcpp::NumericVector DeltaContractionMapping(
    Rcpp::NumericVector initialDeltas,
    Rcpp::NumericVector realOldCounts,
    Rcpp::NumericVector realNewCounts,
    Rcpp::NumericVector percentRents,
    Rcpp::NumericVector medval,
    double currentMu,
    double tolerance,
    bool verbose=false
  )
  {
    Rcpp::NumericVector currentDeltas, previousDeltas, predictedOldShares, realOldShares, realnewShares;
    
    int inner, i, j, N;
    
    if (verbose)
      Rcpp::Rcout << "Entered DeltaContractionMapping" << std::endl;

    currentDeltas = initialDeltas;
    // how do I get sum of counts?
    
    previousDeltas = currentDeltas + tolerance * 100;
    
    realOldShares = realOldCounts/sum(realOldCounts);
    realNewShares = realNewCounts/sum(realNewCounts);
    
    
    
    inner = 0; 
    N = initialDeltas.size();
    
    for (i=0, i < N, i++){
      
      
    }
    
  }

//  

//  

//  
//  
//   # to make sure this object exists and fails first condition
//  
//  # This may need changing again later, because there is not a closed system
//  realOldShares <- realOldCounts/sum(realOldCounts)
//  realNewShares <- realNewCounts/sum(realNewCounts)
//  
//  inner <- 0 
//  
//  repeat
//  {
//    # BUG CATCHING
//    if (any(is.na(previousDeltas))){
//      stop("Some elements in previousDeltas are NA")
//    }
//    
//    if (any(is.na(currentDeltas))){
//      stop("Some elements in currentDeltas are NA")
//    }
//    
//    tmp <- abs(previousDeltas-currentDeltas) 
//    tmp <- tmp < tolerance
//        
//    if (all(tmp))
//    {
//      
//      if (verbose){
//        print("exiting DeltaContractionMapping loop: condition met")
//      }
//      
//      break
//    } 
//    else 
//    {
//      
//      if (verbose){
//        tmp <- abs(previousDeltas - currentDeltas)
//        print("repeating DeltaContractionMapping loop: condition not met")
//        cat("Discrepancy: mean\t", mean(tmp), "\tmax:\t", max(tmp),"\n")
//        cat("inner increased: ", inner, "\n")
//      }
//      
//      inner <- inner + 1   
//      
//      
//      predictedNewShares <- PredictNewShares(
//        deltas=currentDeltas, 
//        mu=currentMu, 
//        oldCounts=realOldCounts, 
//        percentRents=percentRents, 
//        medval=medval,
//        verbose=verbose)     
//      
//      updatedDeltas <- currentDeltas + log(realNewShares) - log(predictedNewShares)
//      updatedDeltas <- updatedDeltas - mean(updatedDeltas)
//      
//      tmp <- which(abs(updatedDeltas - currentDeltas) >= tolerance)
//      
//      if (verbose){
//        cat("Changing ", length(tmp), " elements of delta vector\n")  
//      }
//      
//      previousDeltas <- currentDeltas
//      currentDeltas[tmp] <- updatedDeltas[tmp] # change the deltas where the discrepency is greater than the tolerance
//    }
//  }
//  
//  return(currentDeltas)  
//}
//
//
//DeltaContractionMapping <- cmpfun(DeltaContractionMapping)



// [[Rcpp::export]]     
void StageOne(
  Rcpp::NumericVector oldCounts,
  Rcpp::NumericVector newCounts,  
  Rcpp::NumericVector percentRents,
  Rcpp::NumericVector medval,
  double realStayerProp,
  bool  verbose = false, 
  double initLowerLimit=-0.15,
  double initUpperLimit=0.10,
  double muTol=0.001,
  double deltaTol=0.001,
  double maxIt=1000

  ){
    int N, iteration;
    double upperLim, lowerLim, curMu, finalMu;
    Rcpp::NumericVector newDeltas, curDeltas, initDeltas, finalDeltas;
    
    
    Rcpp::Rcout << "Length of oldCounts is " << oldCounts.size() << std::endl;
    
    if (verbose) 
      Rcpp::Rcout << "Starting FindMu" << std::endl;

//  # ERROR CHECKING
//  if (initLowerLimit > initUpperLimit)
//    stop("initLowerLimit greater than initUpperLimit");
//  
//  
//  if (any(percentRents > 1) | any(percentRents < 0)){
//    stop("percentRents not a valid proportion")
//  }
//  
    
    
    N = medval.size();
    iteration = 0;
    upperLim = initUpperLimit;
    lowerLim = initLowerLimit;
    curMu = (upperLim + lowerLim) / 2;
    
    if (verbose){
      Rcpp::Rcout << "UpperLim:\t" << upperLim << "\tlowerLim:\t" << lowerLim << "\tcurMu:\t" << curMu << std::endl;
      Rcpp::Rcout << "Starting first contraction mapping" << std::endl;
    }
    
    newDeltas = DeltaContractionMapping(
      initialDeltas=initDeltas,
      currentMu=curMu,
      tolerance=deltaTol,
      realOldCounts=oldCounts,
      realNewCounts=newCounts,
      percentRents=percentRents,
      medval=medval,
      verbose=verbose
    )


    return Rcpp::List::create(
        Rcpp::Named("mu", finalMu),
        Rcpp::Named("deltas", finalDeltas)
      );
  }
  
//  
//  
//  

//  
//  
//  
//  

//  
//  

//  
//  

//

//  
//  
//  if (verbose){
//    print("Finished first contraction mapping")
//    print("About to enter CalcPredictedStayerProp")
//  }
//  
//  predictedStayerProp <- CalcPredictedStayerProp(popT2=newCounts, deltas=newDeltas, mu=curMu, percentRents=percentRents, medval=medval)
//  difStayerProp <- realStayerProp - predictedStayerProp
//  
//  if (verbose){
//    cat("StayerProp. Real:\t", realStayerProp, "\tPredicted:\t",predictedStayerProp, "\tDif:\t", difStayerProp, "\n")
//    print("About to enter while loop")
//
//  }
//  while( (iteration < maxIt) & (abs(difStayerProp) > muTol)){
//
//    if (verbose){
//      print("In while loop within FindMu")
//      cat("FindMu iteration:\t", iteration, "\n")
//    }    
//    
//    if (difStayerProp > 0){
//      lowerLim <- curMu
//      if (verbose){
//        print("difStayerProp is positive")
//      }
//      
//    } else if (difStayerProp < 0){
//      upperLim <- curMu
//      if (verbose){
//        print("difStayerProp is negative")
//      }
//    } else {
//      stop("difStayerProp neither positive nor negative")
//    }
//    
//    if (verbose){
//      cat("curMu is ", curMu, "\n")
//    }
//    
//    curMu <- (lowerLim + upperLim)/2
//    
//    if (verbose){
//      cat("curMu is now ", curMu, "\n")
//    }
//    
//    prevDeltas <- newDeltas
//    
//    if (verbose){
//      print("Entering DeltaContractionMapping for a subsequent time")
//    }
//    
//    newDeltas <- DeltaContractionMapping(
//      initialDeltas=prevDeltas,
//      currentMu=curMu,
//      tolerance=deltaTol,
//      realOldCounts=oldCounts,
//      realNewCounts=newCounts,
//      percentRents=percentRents,
//      medval=medval,
//      verbose=verbose
//    )
//    
//    if (verbose){
//      print("Exited DeltaContractionMapping a subsequent time")
//      print("Entering CalcPredictedStayerProp")
//    }
//    
//    predictedStayerProp <- CalcPredictedStayerProp(
//      popT2=newCounts, 
//      deltas=newDeltas, 
//      mu=curMu, 
//      percentRents=percentRents, 
//      medval=medval)
//    
//    difStayerProp <- realStayerProp - predictedStayerProp
//    if (verbose){
//      cat("difStayerProp now:\t", difStayerProp, "\n")
//    }
//    iteration <- iteration + 1
//    
//    if (verbose){
//      cat("iteration now ", iteration, "\n")
//      print("Reached end of while loop")
//    }
//    
//  }
//  # Now need to calculate different between predicted and actual new shares
//  
//  output <- list(
//    mu=curMu,
//    deltas=newDeltas
//    )
//  
//  return(output)
//}
//
//      
//      
//  
