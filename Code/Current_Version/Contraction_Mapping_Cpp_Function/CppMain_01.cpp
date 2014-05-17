#include<cstdio>
#include<Rcpp.h>
#include<vector>

using namespace Rcpp;
using namespace std; 


class Contraction{

// Private fields 

// DATA  

  NumericMatrix movingcost;
//  NumericMatrix s; //transition matrix
  NumericMatrix util; 
  
  std::vector<long double> oldcounts;
  std::vector<long double> newcounts;
  
  std::vector<long double> xf;
  std::vector<long double> x1;
  std::vector<long double> x;
  std::vector<long double> xguess;   
//  std::vector<double> oldcounts_added;
//  std::vector<double> newcounts_added;
  
  long double stayerprop_real;
  long double stayerprop_predicted;
  long double stayerprop_dif;
  
  long double stayp; 
  long double stay;
  long double pstay; 
  long double mu_diff;
  long double mu_guess; 

// PARAMETERS
  long double mu_upper, mu_lower;
  long double tol_mu;
  long double tol_delta;
  long double correction; // the amount to add to each cellcount
  
  int maxit_outer, maxit_inner;

// DERIVED QUANTITIES  

  long double totpop;
  long double totpop_t1;
  long double totpop_t2;
  long double dpop; 
  
  int n;
  int n1;
  
  std::vector<long double> oldshares;
  std::vector<long double> newshares;
    

    

// ESTIMATED QUANTITIES
  std::vector<long double> deltas;
  std::vector<long double> deltas_new; 
  

  

// INTERMEDIATE QUANTITIES  

  std::vector<long double> newshares_predicted; 
  std::vector<long double> newshares_real_predicted_dif;
  

  
  std::vector<long double> newcounts_predicted;   

  
// CONTROL QUANTITIES
  
  bool dbgflag;
  bool verboseflag; 
  
  bool deltas_OK;
  
  
// METHOD TRACKERS  
  int timescalled_CalcDenom;
  int timescalled_ClosedSystemCounts;
  int timescalled_ClosedSystemShares;
  int timescalled_CalcS;
  int timescalled_RunContraction;
  int timescalled_RunOuter;
  int timescalled_RunInner; 
  int timescalled_PredictNewCounts;
  int timescalled_PredictNewShares;
  int timescalled_CalcRealPredictedShareDif;
  int timescalled_PredictStayerProp;
  int timescalled_TimminsInner;
  int timescalled_TimminsOuter;
  int timescalled_TimminsOuterTest;

  
  
  
// Private methods 



//  double CalcDenom(int area);
  
  
//  void ClosedSystemCounts(); //Closed system count estimator
//  void ClosedSystemShares();  

//  void CalcUtil();
  
  
public:


  //Constructors
  Contraction();
  Contraction(List In_);
  
  //Destructor
  
  ~Contraction() {Rcout<< "Destructor called\n";};


//  void CalcS();
  //RunContractionMapping
//  void RunContraction(); // Run contraction
  // Cascades to RunOuter (adjusts mu)

//  void RunOuter();
  // Cascades to RunInner (adjusts deltas)
//  void RunInner();  
  
//  NumericVector CalcPredictedShares();
//  void PredictNewCounts();
//  void PredictNewShares();
  // Return output
//  NumericVector GetOutputs();

//  void CalcRealPredictedShareDif();
  
//  void PredictStayerProp();
  void InitialiseTrackers();
  
  void TimminsInner(); // Attempt to replicate Timmins code blindly from 300 onwards
  void TimminsOuterTest(); 
  void TimminsOuter();
  void TimminsRunContraction();
  void TimminsPop();
  
  List ExtractEstimates();
  
};

void Contraction::TimminsRunContraction()
{
    Rcout << "TRC\n";
    TimminsPop();
    mu_diff = 99999.99; 
    stayp = stayerprop_real; 
    
    TimminsOuter();
}


void Contraction::TimminsOuter()
{
    Rcout << "TO\n";
    
    timescalled_TimminsOuter++;
   
    xguess.resize(n1, 0.0);
    Rcout << n1 << endl;
    
//   deltas_new.resize(n1, 0.0);

   
    bool repeat_loop = true; 
    unsigned int repeat_nums = 0; 
    
    while (repeat_loop){
        repeat_nums++;

        mu_guess = (mu_upper + mu_lower) / 2.0;
        
        TimminsOuterTest();
         
        if (mu_diff > 0.0){
            mu_lower = mu_guess;
        } else {
            mu_upper = mu_guess;
        }
        Rcout << "200\n";

        for (int j = 0; j < n; j++){
            xguess[j] = x[j];
        }
        

        if (repeat_nums > maxit_outer){
            repeat_loop = false;
        } 
        
        if (abs(mu_diff) < tol_mu){
            repeat_loop = false;
        }
        
    }
    
    for (unsigned int j = 0; j < n; j++){
        deltas[j] = x[j];
    }
    
}


void Contraction::TimminsOuterTest()
{
    Rcout << "TOT\n";           
    timescalled_TimminsOuterTest++;
    Rcout << n1 << endl;
    x1.resize(n1);
    xf.resize(n1);
    x.resize(n1);
    
    TimminsInner();    
    
    for (int j = 0; j < n1; j++){

        long double test = abs(xf[j] - x1[j]);
        if (test > tol_delta){
            for (int k = 0; k < n1; k++){
                x1[k] = xf[k];
            }
            TimminsInner();
        }
        
    }
   
    mu_diff = stayp - pstay;
    
}

void Contraction::TimminsInner()
{
    Rcout << "TI\n";
    timescalled_TimminsInner++; Rcout << "240\n";

    
    Rcout << "243\n";
    for (int i = 0; i < n1; i++){
        x1[i] = xguess[i];
        xf[i] = x[i];
    }
    Rcout << "248\n";
   
    util = Rcpp::clone(movingcost);

    Rcout << "249\n";
    for (int j =0; j < n; j++){
        for (int k = 0; k < n; k++){
//            Rcout << x1[j] << "\t" << x1[j] << "\t" << mu_guess << "\t" << movingcost(j,k) << endl;
//                Rcout << "[ " << j << ", " << k << "]\t" << x1[k] - x1[j] - mu_guess * movingcost(j, k) << endl;

            util(k, j) = x1[k] - x1[j] - mu_guess * movingcost(j, k);
        }
    }
    Rcout << "273\n";
     
    Rcout << "256\n";
    std::vector <long double> denom;
    denom.resize(n1, 0.0);
    
    Rcout << "260\n";
    for (unsigned int j = 0; j < n; j++){
        for (unsigned int k = 0; k < n; k++){
            denom[j] += exp(util(k, j));
        }
    }
   
    Rcout << "267\n";
    std::vector <long double> rhs;
    rhs.resize(n1, 0.0);
    
    Rcout << "271\n";
     
    for (unsigned int i = 0; i < n; i++){
        for (unsigned int j = 0; j < n; j++){
            double numer = exp(util(i, j));
            rhs[i] += (numer/denom[j])* oldcounts[i];
        }
    }

    Rcout << "279\n";
    
    stay = 0;
    for (unsigned int i = 0; i < (n -1); i++){
        stay += exp(util(i,i))/oldcounts[i];
    }
    stay = stay / (totpop - newcounts[newcounts.size()]);
    
    pstay = 0;
    Rcout << "308\n";
    
    for (unsigned int i = 0; i < (n-1); i++){
        pstay+= (exp(util(i,i)) / denom[i])/(n1 - 1);
    }
    Rcout << "313\n";
    
    newshares.resize(n1);
    newshares_predicted.resize(n1);
    
     for (unsigned int i=0; i < n; i++){
         newshares[i] = newcounts[i] / totpop;
         newshares_predicted[i] = rhs[i] / totpop;         
     }
     Rcout << "319\n";

     long double xavg = 0; 
     deltas_new.resize(n1);
     for (unsigned int j =0; j < n1; ++j){
         xf[j] = x1[j] + (log(newshares[j]) - log(newshares_predicted[j]));
         xavg += xf[j] / (long double) n1;
     }
     Rcout << "330\n";
         
     for (unsigned int i = 0; i < n1; ++i){
         xf[i] = xf[i] - xavg;
     }
    
    Rcout << "336\n";     
}


void Contraction::TimminsPop()
{
    Rcout << "TP\n";
    totpop_t1 = 0;
    totpop_t2 = 0;
    
    for (int i = 0; i < n1; ++i){
        totpop_t1 += oldcounts[i];
        totpop_t2 += newcounts[i];
    }

    dpop = totpop_t2 - totpop_t1;
    
    long double tmp = 2.0 * abs(dpop);
    oldcounts.push_back(tmp);
    newcounts.push_back(tmp - dpop);
    

}

/*
void Contraction::CalcUtil()
{
    util = Rcpp::clone(movingcost);
    
    
    for (unsigned int j = 0; j < deltas.size(); ++j){
        for (unsigned int k = 0; k < deltas.size(); ++k){
            util(k,j) = deltas[k]-deltas[j]- mu_guess*movingcost(k,j);
        }
    }
}
*/

//Constructors functions

Contraction::Contraction(List In_)
{
    Rcout << "Con\n";
    InitialiseTrackers();
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
  
  stayerprop_real = as<long double> (data["stayerprop"]);
  
  movingcost = as<NumericMatrix> (In["movingcostmatrix"]);
  
  // Counts objects
  oldcounts = Rcpp::as<std::vector<long double> > (counts["oldcounts"]);
  newcounts = Rcpp::as<std::vector<long double> > (counts["newcounts"]);
  
  // Calculate total population, redefine oldcounts and newcounts 
  // as closed system
  //ClosedSystemCounts();
  //ClosedSystemShares();
  
  mu_upper = as<long double>(params["mu_upper"]);
  mu_lower = as<long double>(params["mu_lower"]);
  
  maxit_outer = as<int> (params["maxit_outer"]);
  maxit_inner = as<int> (params["maxit_inner"]);

  // utils objects
    
  deltas = Rcpp::as<std::vector<long double> > (utils["deltas"]);
//  Rcout << deltas.size() << endl; 

    n = (int) deltas.size();
    n1 = n + 1;
    
 // params objects  

  tol_delta = as<long double> (params["tol_delta"]);
  tol_mu = as<long double> (params["tol_mu"]);
  
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

void Contraction::InitialiseTrackers()
{
  timescalled_CalcDenom = 0;
  timescalled_ClosedSystemCounts = 0;
  timescalled_ClosedSystemShares = 0;
  timescalled_CalcS = 0;
  timescalled_RunContraction = 0;
  timescalled_RunOuter = 0;
  timescalled_RunInner = 0; 
  timescalled_PredictNewCounts = 0 ;
  timescalled_PredictNewShares = 0;
  timescalled_CalcRealPredictedShareDif = 0;
  timescalled_PredictStayerProp = 0;
  timescalled_TimminsInner = 0;
  timescalled_TimminsOuter = 0;
  timescalled_TimminsOuterTest = 0;
}





List Contraction::ExtractEstimates()
{
  return(Rcpp::List::create(
     Rcpp::Named("deltas") = deltas,
     Rcpp::Named("mu") = mu_guess,
     Rcpp::Named("newshares_predicted") = newshares_predicted,
     Rcpp::Named("newshares") = newshares
    )
  );
}



RCPP_MODULE(ContractionClass){
  using namespace Rcpp;
  
  class_<Contraction> ("Contraction")
  
  .constructor<List>()
  .constructor()
  .method("run", &Contraction::TimminsRunContraction)
//  .method("run", &Contraction::RunContraction)
//  .method("getoutputs", &DeltaContraction::GetOutputs)
//  .method("CalcS", &Contraction::CalcS)
//  .method("PredictCounts", &Contraction::PredictNewCounts)
//  .method("PredictShares", &Contraction::PredictNewShares)
//  .method("CompareRealPredicted", &Contraction::CalcRealPredictedShareDif)
  .method("extractests", &Contraction::ExtractEstimates)
  ;
}


/*
double Contraction::CalcDenom(int area)
{
  timescalled_CalcDenom++;

  
  double tmp = 0;
  
  for (unsigned int l = 0; l < deltas.size(); ++l)
    tmp += exp(deltas[l] - deltas[area] - mu * movingcost(l, area));
  
//    tmp += exp(deltas[l] - deltas[area] - mu * movingcost(area, l));

    
  return tmp;  
}
*/


// Closed system: this function adds a fictitous areal unit N+1
// to account for moves out of the geography of interest.
// It also produces an estimate of totpop 
/*
void Contraction::ClosedSystemCounts()
{
    timescalled_ClosedSystemCounts++;
    
    if (verboseflag){
        Rcout << "ClosedSystemCounts" << endl; 
    }
//  Rcout << "Entered ClosedSystemCounts\n";
  double newcounts_sum = 0;
  double oldcounts_sum = 0;
  
  for (unsigned int i = 0; i < oldcounts.size(); ++i)
    oldcounts_sum += oldcounts[i];
  
//  Rcout << "Sum of oldcounts is " << oldcounts_sum << endl;
  
  for (unsigned int i = 0; i < newcounts.size(); ++i)
    newcounts_sum += newcounts[i];
//  Rcout << "Sum of newcounts is " << newcounts_sum << endl;  
    
  double dpop = newcounts_sum - oldcounts_sum;
//  Rcout << "dpop " << dpop << endl;  
  
  double tmp = 2 * abs(dpop);
  oldcounts.push_back(tmp); // this adds tmp to the back of the vector oldcounts
  // i.e. the vector length changes from N to N+1
  // and the value of the last vector is tmp
  
  double tmp2 = tmp - dpop;
  newcounts.push_back(tmp2);
    
  totpop = newcounts_sum + tmp2;
//  Rcout << "totpop " << totpop << endl;
}
*/
// This populates oldshares and newshares based on oldcounts and newcounts 
// respectively. Totpop is needed, which is set in ClosedSystemCounts,
// so ClosedSystemCounts must be run first.
/*
void Contraction::ClosedSystemShares()
{
    timescalled_ClosedSystemShares++;
    
    if (verboseflag){
        Rcout << "ClosedSystemShares\n";
    }
  
  oldshares.resize(oldcounts.size());
  newshares.resize(newcounts.size());
  
    // (double) in the lines below coerces oldcounts from type int to type double
  for (unsigned int i = 0; i < oldshares.size(); ++i)
    oldshares[i] = (double) oldcounts[i] / (double) totpop;
    
  for (unsigned int i = 0; i < newshares.size(); ++i)
    newshares[i] = (double) newcounts[i] / (double) totpop;
  
}
*/
// This sets the algorithm rolling.
// It starts RunOuter, 
// Which cascades to RunInner,
// And uses the method of bisection to find Mu
/*
void Contraction::RunContraction()
{
    timescalled_RunContraction++;
    
    if (verboseflag){
        Rcout << "RunContraction\n";
    }
    
  // Empty for now
  
//  Rcout << "RunContraction run\n";
  int outercount = 0;
  bool mu_OK = false;
  
  // Keep running outer until mu within acceptable tolerance
  while (!mu_OK){
    CalcS();
    PredictStayerProp();
    RunOuter();// Affected by current value of Mu
    outercount++;
    if (outercount > maxit_outer) mu_OK = true;
    Rcout << "|";
    stayerprop_dif = stayerprop_predicted - stayerprop_real;
    if (verboseflag){
        Rcout << "stayerprop_predicted: " << stayerprop_predicted << endl;
        Rcout << "stayerprop_real: " << stayerprop_real << endl;
        Rcout << "stayerprop_dif: " << stayerprop_dif << endl;    
    }
    if (stayerprop_dif > 0) {
      if (verboseflag){
          Rcout << "stayerprop_dif is > 0 so changing mu_upper from ";
          Rcout << mu_upper << " to " << mu << endl;
      }  
      // More people are predicted to stay than actually stay
      // So.. predicted effect of moving costs is too high
      // So.. lower mu by moving upper
      
      mu_upper = mu; // THIS IS THE RIGHT WAY AROUND!!!!
    } else {
        if (verboseflag){
          Rcout << "stayerprop_dif is < 0 so changing mu_lower from ";
          Rcout << mu_lower << " to " << mu << endl;
        }  

      // More people stay than are predicted to stay
      // So.. predicted effect of moving costs is too low
      // So.. raise my by moving lower
      mu_lower = mu;
    }
    Rcout << endl;
    if (abs(stayerprop_dif) < tol_mu) mu_OK = true;
    if (((mu_upper - mu_lower)/2) < tol_mu) mu_OK = true;
  }
  
  Rcout << "\nTimes each method called:\n";
  Rcout << "CalcDenom: " << timescalled_CalcDenom << endl;
  Rcout << "ClosedSystemCounts: " << timescalled_ClosedSystemCounts << endl;
  Rcout << "ClosedSystemShares: " << timescalled_ClosedSystemShares << endl;    
  Rcout << "CalcS: " << timescalled_CalcS << endl;
  Rcout << "RunContraction: " << timescalled_RunContraction << endl;
  Rcout << "RunOuter: " << timescalled_RunOuter << endl;
  Rcout << "RunInner: " << timescalled_RunInner << endl;
  Rcout << "PredictNewCounts: " << timescalled_PredictNewCounts << endl;
  Rcout << "PredictNewShares: " << timescalled_PredictNewShares << endl;
  Rcout << "CalcRealPredictedShareDif: " << timescalled_CalcRealPredictedShareDif << endl;
  Rcout << "PredictStayerProp: " << timescalled_PredictStayerProp << endl;
  
}
*/


/*
// This recalculates mu, and runs RunInner
// while the condition for all deltas has not been met
void Contraction::RunOuter()
{
    timescalled_RunOuter++;
    
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
*/

/*
void Contraction::RunInner()
{
    timescalled_RunInner++;
  
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
  
  // Normalisation
  double tmp = 0;
  for (unsigned int i =0; i < (deltas.size() - 1) ; ++i){
      tmp +=deltas[i]; // look for a smarter way of doing this
  }
  
  tmp = tmp / (deltas.size() - 1);
  
  for (unsigned int i = 0; i < (deltas.size() - 1) ; ++i){
      deltas[i] = deltas[i] - tmp;
  }
  
  if (verboseflag){
      Rcout << "RunInner\n";
      Rcout << "Deltas changed: " << numchanged << endl;
  }
  
}
*/

/*
void Contraction::CalcRealPredictedShareDif()
{
    timescalled_CalcRealPredictedShareDif++;
    
    if (verboseflag){
        Rcout << "CalcRealPredictedShareDif\n";
    }
  
  newshares_real_predicted_dif.resize(newshares.size());
  
  for (unsigned int i = 0; i < newshares_real_predicted_dif.size(); ++i){
    double tmp = newshares_predicted[i] - newshares[i];
    newshares_real_predicted_dif[i] = tmp;
  }
//    Rcout << endl;
}
*/

/*
void Contraction::PredictStayerProp()
{
    timescalled_PredictStayerProp++;
    
  double tmp = 0 ;
  
    // tmp += something
    // is equivalent to
    // tmp = tmp + something
    // using (deltas.size() - 1) because the last areal unit is the catch-all region
  for (unsigned int i = 0; i < (deltas.size() - 1); ++i )
  {
    tmp += s(i,i) * (double) oldcounts[i];
  }
  
  stayerprop_predicted = tmp / (double) totpop;
//  Rcout << "PredictStayerProp debugging\n";
//  Rcout << "tmp: " << tmp << endl;
//  Rcout << "totpop: " << totpop << endl;
//  Rcout << "Predicted stayerprop: " << stayerprop_predicted << endl;
}
*/



/*
void Contraction::CalcS()
{
    timescalled_CalcS++;
    
//  Rcout << "Entered CalcS\n";
  s = Rcpp::clone(movingcost);
  
  
  for (int j = 0; j < s.nrow(); ++j){
    for (int k = 0; k < s.ncol(); ++k){
        
      double tmp = exp(deltas[j] - deltas[k] - mu * movingcost(j, k));
//      double tmp = exp(deltas[j] - deltas[k] - mu * movingcost(k, j));
      
      tmp = tmp / CalcDenom(k);
      s(j,k) = tmp;     
      
//      s(k,j) = tmp;     
//      Rcout << tmp << " ";
    }
//    Rcout << endl;
  }
}
*/



/*
void Contraction::PredictNewCounts()
{
    timescalled_PredictNewCounts++;
    
//  Rcout << "Entered PredictNewCounts\n";
  newcounts_predicted.resize(deltas.size());
  
  for (unsigned int j = 0; j < newcounts_predicted.size(); ++j){
    double tmp = 0;
    for (unsigned int k = 0; k < deltas.size(); ++k){

      double tmp2 = exp(deltas[j] - deltas[k] - mu * movingcost(j, k));
//      double tmp2 = exp(deltas[j] - deltas[k] - mu * movingcost(k, j));

      tmp2 = tmp2 * (double) oldcounts[k] / CalcDenom(k);
      tmp += tmp2;
    }
//    Rcout << tmp << " ";
    newcounts_predicted[j] = tmp;
  }
}
*/



/*
void Contraction::PredictNewShares()
{
    timescalled_PredictNewShares++;
//  Rcout << "Entered PredictNewShares\n";
  newshares_predicted.resize(deltas.size());
  
  for (unsigned int j = 0; j < newshares_predicted.size(); ++j){
    double tmp = 0;
    for (unsigned int k = 0; k < deltas.size(); ++k){
        
      double tmp2 = exp(deltas[j] - deltas[k] - mu * movingcost(j, k));
//      double tmp2 = exp(deltas[j] - deltas[k] - mu * movingcost(k, j));

      tmp2 = tmp2 * (double) oldcounts[k] / CalcDenom(k);
      tmp += tmp2;
    }
    tmp = tmp / totpop;
//    Rcout << tmp << " ";
    newshares_predicted[j] = tmp;
  }
}
*/