#include<cstdio>
#include<Rcpp.h>
#include<vector>

using namespace Rcpp;
using namespace std; 


class Contraction{

// Private fields 

// DATA  

  NumericMatrix movingcost;
  NumericMatrix s; //transition matrix
  
  std::vector<int> oldcounts;
  std::vector<int> newcounts;
  
  
  double stayerprop_real;
  double predictednewprop;

// PARAMETERS
  double mu_upper, mu_lower;
  double tol_mu;
  double tol_delta;
  double correction; // the amount to add to each cellcount
  
  int maxit_outer, maxit_inner;

// DERIVED QUANTITIES  

  int totpop;
  std::vector<double> oldshares;
  std::vector<double> newshares;


// ESTIMATED QUANTITIES
  std::vector<double> deltas;
  double mu;
  
  

// INTERMEDIATE QUANTITIES  

  std::vector<double> newshares_predicted; 
  std::vector<double> newshares_real_predicted_dif;
  
  double stayerprop_predicted;
  
  double stayerprop_dif;
  
  std::vector<double> newcounts_predicted;   

  
// CONTROL QUANTITIES
  
  bool dbgflag;
  bool verboseflag; 
  
  bool deltas_OK;
  
  
  
  

  
  
  
// Private methods 



  double CalcDenom(int area);
  
  
  void ClosedSystemCounts(); //Closed system count estimator
  void ClosedSystemShares();  

  
  
public:


  //Constructors
  Contraction();
  Contraction(List In_);
  
  //Destructor
  
  ~Contraction() {Rcout<< "Destructor called\n";};


  void CalcS();
  //RunContractionMapping
  void RunContraction(); // Run contraction
  // Cascades to RunOuter (adjusts mu)
  void RunOuter();
  // Cascades to RunInner (adjusts deltas)
  void RunInner();  
  
//  NumericVector CalcPredictedShares();
  void PredictNewCounts();
  void PredictNewShares();
  // Return output
//  NumericVector GetOutputs();

  void CalcRealPredictedShareDif();
  
  void PredictStayerProp();
  
  List ExtractEstimates();
  
};


//Constructors functions

Contraction::Contraction(List In_)
{
//  Rcout << "Input List constructor\n";
  // Tracking vals
  
  
  // Main object 
  List In = In_;
  
  // One layer down
  List data = In["data"];
  List utils = In["utils"];
  List params = In["params"];
  List dbg = In["dbg"];
  

  List counts = data["counts"];
  
  stayerprop_real = as<double> (data["stayerprop"]);
  
  movingcost = as<NumericMatrix> (In["movingcostmatrix"]);
  
  // Counts objects
  oldcounts = Rcpp::as<std::vector<int> > (counts["oldcounts"]);
  newcounts = Rcpp::as<std::vector<int> > (counts["newcounts"]);
  
  // Calculate total population, redefine oldcounts and newcounts 
  // as closed system
  ClosedSystemCounts();
  ClosedSystemShares();
  
  mu_upper = as<double>(params["mu_upper"]);
  mu_lower = as<double>(params["mu_lower"]);
  
  maxit_outer = as<int> (params["maxit_outer"]);
  maxit_inner = as<int> (params["maxit_inner"]);

  // utils objects
    
  deltas = Rcpp::as<std::vector<double> > (utils["deltas"]);

 // params objects  

  tol_delta = as<double> (params["tol_delta"]);
  tol_mu = as<double> (params["tol_mu"]);
  
  // dbg objects  
  
  dbgflag = as<bool> (dbg["dbgflag"]);
  verboseflag = as<bool> (dbg["verboseflag"]);
  
  if (verboseflag){
      Rcout << "Constructor completed" << endl;
  }

}  

Contraction::Contraction()
{
  Rcout << "Empty constructor\n";
  // Main object 
  
}  

double Contraction::CalcDenom(int area)
{

  double tmp = 0;
  
  for (unsigned int l = 0; l < deltas.size(); ++l)
    tmp += exp(deltas[l] - deltas[area] - mu * movingcost(l, area));
    
  return tmp;  
}



// Closed system: this function adds a fictitous areal unit N+1
// to account for moves out of the geography of interest.
// It also produces an estimate of totpop 
void Contraction::ClosedSystemCounts()
{
    if (verboseflag){
        Rcout << "ClosedSystemCounts" << endl; 
    }
//  Rcout << "Entered ClosedSystemCounts\n";
  int newcounts_sum = 0;
  int oldcounts_sum = 0;
  
  for (unsigned int i = 0; i < oldcounts.size(); ++i)
    oldcounts_sum += oldcounts[i];
  
//  Rcout << "Sum of oldcounts is " << oldcounts_sum << endl;
  
  for (unsigned int i = 0; i < newcounts.size(); ++i)
    newcounts_sum += newcounts[i];
//  Rcout << "Sum of newcounts is " << newcounts_sum << endl;  
    
  int dpop = newcounts_sum - oldcounts_sum;
//  Rcout << "dpop " << dpop << endl;  
  
  int tmp = 2 * abs(dpop);
  oldcounts.push_back(tmp);
  
  int tmp2 = tmp - dpop;
  newcounts.push_back(tmp2);
    
  totpop = newcounts_sum + tmp2;
//  Rcout << "totpop " << totpop << endl;
}

// This populates oldshares and newshares based on oldcounts and newcounts 
// respectively. Totpop is needed, which is set in ClosedSystemCounts,
// so ClosedSystemCounts must be run first.
void Contraction::ClosedSystemShares()
{
    if (verboseflag){
        Rcout << "ClosedSystemShares\n";
    }
  
  oldshares.resize(oldcounts.size());
  newshares.resize(newcounts.size());
  
  for (unsigned int i = 0; i < oldshares.size(); ++i)
    oldshares[i] = (double) oldcounts[i] / (double) totpop;
    
  for (unsigned int i = 0; i < newshares.size(); ++i)
    newshares[i] = (double) newcounts[i] / (double) totpop;
  
}

// This sets the algorithm rolling.
// It starts RunOuter, 
// Which cascades to RunInner,
// And uses the method of bisection to find Mu
void Contraction::RunContraction()
{
    if (verboseflag){
        Rcout << "RunContraction\n";
    }
    
  // Empty for now
  
//  Rcout << "RunContraction run\n";
  int outercount = 0;
  bool mu_OK = false;
  
  // Keep running outer until mu within acceptable tolerance
  while (!mu_OK){
    RunOuter();// Affected by current value of Mu
    outercount++;
    if (outercount > maxit_outer) mu_OK = true;
    Rcout << "|";
    PredictStayerProp();
    stayerprop_dif = stayerprop_predicted - stayerprop_real;
    if (verboseflag){
        Rcout << "stayerprop_predicted: " << stayerprop_predicted << endl;
        Rcout << "stayerprop_real: " << stayerprop_real << endl;
        Rcout << "stayerprop_dif: " << stayerprop_dif << endl;    
    }
    if (stayerprop_dif > 0) {
        
      // More people are predicted to stay than actually stay
      // So.. predicted effect of moving costs is too high
      // So.. lower mu by moving upper
      mu = mu_upper; // THIS IS THE RIGHT WAY AROUND!!!!
    } else {
      // More people stay than are predicted to stay
      // So.. predicted effect of moving costs is too low
      // So.. raise my by moving lower
      mu = mu_lower;
    }
    Rcout << endl;
    if (abs(stayerprop_dif) < tol_mu) mu_OK = true;
  }
}

// This recalculates mu, and runs RunInner
// while the condition for all deltas has not been met
void Contraction::RunOuter()
{
  
  mu = (mu_upper + mu_lower) / 2;

    if (verboseflag){
        Rcout << "RunOuter\n";
        Rcout << "mu_upper: " << mu_upper << endl;
        Rcout << "mu_lower: " << mu_lower << endl;
        Rcout << "mu: " << mu << endl;
    }  
  
  
  int innercount = 0;
  deltas_OK = false;
  
  while(!deltas_OK){
    Rcout << ".";
    RunInner();
    innercount++;
    if (innercount > maxit_inner) deltas_OK = true;
  
  }
  if (verboseflag) Rcout << "innercount: " << innercount << endl;
    
}

void Contraction::RunInner()
{

  
  // Calc predicted proportions given the current mu and delta estimates
  
  PredictNewShares();
  // compare predicted with real proportions 
  int numchanged = 0;
  CalcRealPredictedShareDif();
  
  deltas_OK = true;
  
  for (unsigned int i = 0; i < (deltas.size() - 1) ; ++i){ // deltas.size() - 1 because not interested
  // in last areal unit (catch-all for outside areas)
    if (abs(newshares_real_predicted_dif[i]) > tol_delta){
      deltas_OK = false; // if any of the differences between predicted and real area shares
      // is greater than the tolerance on delta, then deltas_OK is set to false,
      // which means RunOuter() will call RunInner() again
      numchanged++; 
      // Contraction mapping of delta
      deltas[i] = deltas[i] + (log(newshares[i]) - log(newshares_predicted[i]));
    }
  }
  if (verboseflag){
      Rcout << "RunInner\n";
      Rcout << "Deltas changed: " << numchanged << endl;
  }
  
}


void Contraction::CalcRealPredictedShareDif()
{
    if (verboseflag){
        Rcout << "CalcRealPredictedShareDif\n";
    }
  
  newshares_real_predicted_dif.resize(newshares.size());
  
  for (unsigned int i = 0; i < newshares_real_predicted_dif.size(); ++i){
//    Rcout << i << ": " << newshares_predicted[i] << " " << newshares[i] << " " << endl;
    double tmp = newshares_predicted[i] - newshares[i];
//    Rcout << tmp << " ";
   newshares_real_predicted_dif[i] = tmp;
  }
//    Rcout << endl;
}

  


void Contraction::PredictStayerProp()
{
  double tmp = 0 ;
  
  for (unsigned int i = 0; i < (deltas.size() - 1); ++i )
  {
    tmp += s(i,i) * (double) oldcounts[i];
  }
  
  stayerprop_predicted = tmp / (double) totpop;
  Rcout << "PredictStayerProp debugging\n";
  Rcout << "tmp: " << tmp << endl;
  Rcout << "totpop: " << totpop << endl;
  Rcout << "Predicted stayerprop: " << stayerprop_predicted << endl;
}




void Contraction::CalcS()
{
//  Rcout << "Entered CalcS\n";
  s = Rcpp::clone(movingcost);
  
  for (int j = 0; j < movingcost.nrow(); ++j){
    for (int k = 0; k < movingcost.ncol(); ++k){
//      Rcout << movingcost(j, k) << " ";
    }
//    Rcout << endl;
  }
  
  for (int j = 0; j < s.nrow(); ++j){
    for (int k = 0; k < s.ncol(); ++k){
      double tmp = exp(deltas[j] - deltas[k] - mu * movingcost(j, k));
      tmp = tmp / CalcDenom(k);
      s(j,k) = tmp;     
//      Rcout << tmp << " ";
    }
//    Rcout << endl;
  }
}



void Contraction::PredictNewCounts()
{
//  Rcout << "Entered PredictNewCounts\n";
  newcounts_predicted.resize(deltas.size());
  
  for (unsigned int j = 0; j < newcounts_predicted.size(); ++j){
    double tmp = 0;
    for (unsigned int k = 0; k < deltas.size(); ++k){
      double tmp2 = exp(deltas[j] - deltas[k] - mu * movingcost(j, k));
      tmp2 = tmp2 * (double) oldcounts[k] / CalcDenom(k);
      tmp += tmp2;
    }
//    Rcout << tmp << " ";
    newcounts_predicted[j] = tmp;
  }
}



void Contraction::PredictNewShares()
{
//  Rcout << "Entered PredictNewShares\n";
  newshares_predicted.resize(deltas.size());
  
  for (unsigned int j = 0; j < newshares_predicted.size(); ++j){
    double tmp = 0;
    for (unsigned int k = 0; k < deltas.size(); ++k){
      double tmp2 = exp(deltas[j] - deltas[k] - mu * movingcost(j, k));
      tmp2 = tmp2 * (double) oldcounts[k] / CalcDenom(k);
      tmp += tmp2;
    }
    tmp = tmp / totpop;
//    Rcout << tmp << " ";
    newshares_predicted[j] = tmp;
  }
}

List Contraction::ExtractEstimates()
{
  return(Rcpp::List::create(
     Rcpp::Named("deltas") = deltas,
     Rcpp::Named("mu") = mu
    )
  );
}



RCPP_MODULE(ContractionClass){
  using namespace Rcpp;
  
  class_<Contraction> ("Contraction")
  
  .constructor<List>()
  .constructor()
  .method("run", &Contraction::RunContraction)
//  .method("getoutputs", &DeltaContraction::GetOutputs)
  .method("CalcS", &Contraction::CalcS)
  .method("PredictCounts", &Contraction::PredictNewCounts)
  .method("PredictShares", &Contraction::PredictNewShares)
  .method("CompareRealPredicted", &Contraction::CalcRealPredictedShareDif)
  .method("extractests", &Contraction::ExtractEstimates)
  ;
}
