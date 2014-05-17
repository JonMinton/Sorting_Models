#include<cstdio>
#include<Rcpp.h>
#include<vector>

using namespace Rcpp;
using namespace std; 


class Contraction{

// Private fields 

// DATA  

  NumericMatrix movingcost;
  NumericMatrix util; 
  
  std::vector<long double> oldcounts;
  std::vector<long double> newcounts;
  
  std::vector<long double> newshares;
  std::vector<long double> newshares_predicted;
  
  std::vector<long double> xf;
  std::vector<long double> x1;
  std::vector<long double> x;
  std::vector<long double> xguess;   
  
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
  
    

    

// ESTIMATED QUANTITIES
  std::vector<long double> deltas;
  std::vector<long double> deltas_new; 
  

  

// INTERMEDIATE QUANTITIES  

  

  
  std::vector<long double> newcounts_predicted;   

  
// CONTROL QUANTITIES
  
  bool dbgflag;
  bool verboseflag; 
  
  bool deltas_OK;
    bool mu_ok;  
  
// METHOD TRACKERS  
  int timescalled_TimminsInner;
  int timescalled_TimminsOuter;
  int timescalled_TimminsOuterTest;
    int timescalled_TimminsRunContraction; 
  
  
  
// Private methods 



  
  
  
  
public:


  //Constructors
  Contraction();
  Contraction(List In_);
  
  //Destructor
  
  ~Contraction() {Rcout<< "Destructor called\n";};


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
    timescalled_TimminsRunContraction++;
    Rcout << "TRC\n";
    TimminsPop();
    mu_diff = tol_mu * 99999.9; 
    stayp = stayerprop_real; 
    
    TimminsOuter();
}


void Contraction::TimminsOuter()
{
    Rcout << "TO\n";
    
    timescalled_TimminsOuter++;
   
    xguess.resize(n1, 0.0);
    Rcout << n1 << endl;
    
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
    
    TimminsInner();    
    for (int k = 0; k < n1; k++){
        x1[k] = xf[k];
    }    
    for (int j = 0; j < n1; j++){

        long double test = abs(xf[j] - x1[j]);
        // xf is nan, so test is nan
        Rcout << j << "\t" << xf[j] << "\t" << x1[j] << "\t" << test << endl;
        if (test > tol_delta){
            for (int k = 0; k < n1; k++){
                x1[k] = xf[k];
            }
            TimminsInner();
        }
        
    }
   
    Rcout << "stayp:\t" << stayp << "\tpstay:\t" << pstay << endl;
    mu_diff = stayp - pstay;
    
}

void Contraction::TimminsInner()
{
    Rcout << "TI\n";
    timescalled_TimminsInner++; Rcout << "240\n";

    for (int i = 0; i < n1; i++){
//        Rcout << "x1: " << x1[i] << "\txguess: " << xguess[i] << "\txf: " << xf[i] << "\tx: " << x[i] << endl; 
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
    std::vector <long double> rhs(n1, 0.0);
    
    Rcout << "271\n";
     
    for (unsigned int i = 0; i < n; i++){
        for (unsigned int j = 0; j < n; j++){
            double numer = exp(util(i, j));
            rhs[i] += (numer/denom[j])* oldcounts[i];
//            Rcout << rhs[i] << endl;
        }
    }

    Rcout << "279\n";
    
    stay = 0;
    for (unsigned int i = 0; i < (n -1); i++){
        stay += exp(util(i,i))/oldcounts[i];
    }
    stay = stay / (totpop_t2 - dpop);
    
    pstay = 0;
    Rcout << "308\n";
    
    Rcout << "About to calc pstay: ";
    for (unsigned int i = 0; i < (n-1); i++){
/*        
        Rcout << pstay << " ";
        if ((i % 10) == 0) {
            Rcout << endl;
        }
*/ // pstay appears to be working correctly        
        pstay+= (exp(util(i,i)) / denom[i])/(n);
    }
    
    Rcout << "313\n";
    

    
     for (unsigned int i=0; i < n; i++){
         newshares[i] = newcounts[i] / totpop_t2;
         newshares_predicted[i] = rhs[i] / totpop_t2;         
//         Rcout << "totpop_t2: " << totpop_t2 << "\tnewcounts: " << newcounts[i] << "\trhs: " << rhs[i] << "\tnewshares: " << newshares[i] << "\tnewshares_predicted: " << newshares_predicted[i] <<endl;
     }
//     Rcout << "319\n";

     long double xavg = 0;    
     for (unsigned int j =0; j < n1; ++j){
         xf[j] = x1[j] + (log(newshares[j]) - log(newshares_predicted[j]));
         xavg += xf[j]; 
//         Rcout << "xf: " << xf[j] << "\tx1: " << x1[j] << "\txavg: " << xavg << endl;
     }
     xavg = xavg / (long double) n1;
//     Rcout << "330\n";
         
     for (unsigned int i = 0; i < n1; ++i){
         xf[i] = xf[i] - xavg;
     }
    
//    Rcout << "336\n";     
}


void Contraction::TimminsPop()
{
    Rcout << "TP\n";
    totpop_t1 = 0;
    totpop_t2 = 0;
    
    for (int i = 0; i < (n-1); i++){
        totpop_t1 += oldcounts[i];
        totpop_t2 += newcounts[i];
        Rcout << "[" << i << "] " << oldcounts[i] << " " << newcounts[i] << "\ttotpop_t1: " << totpop_t1 << "\ttotpop_t2: " << totpop_t2 << endl;
    }

    dpop = totpop_t2 - totpop_t1;

    Rcout << "totpop_t1: " << totpop_t1 << "\ttotpop_t2: " << totpop_t2 << "\tdpop: " << dpop<< endl;
    
    long double tmp = 2.0 * abs(dpop);
    oldcounts.push_back(tmp);
    newcounts.push_back(tmp - dpop);
    

}


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
  xguess = Rcpp::as<std::vector<long double> > (utils["deltas"]);
//  Rcout << deltas.size() << endl; 

    n = (int) deltas.size();
    n1 = n + 1;
    
 // params objects  

  tol_delta = as<long double> (params["tol_delta"]);
  tol_mu = as<long double> (params["tol_mu"]);
  
  // dbg objects  
  
  dbgflag = as<bool> (dbg["dbgflag"]);
  verboseflag = as<bool> (dbg["verboseflag"]);
  
    x1.resize(n1);
    xf.resize(n1);
    x.resize(n1);
    

    newshares.resize(n1);
    newshares_predicted.resize(n1);
     deltas_new.resize(n1);
  
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

  timescalled_TimminsInner = 0;
  timescalled_TimminsOuter = 0;
  timescalled_TimminsOuterTest = 0;
  timescalled_TimminsRunContraction = 0;  
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
  .method("extractests", &Contraction::ExtractEstimates)
  ;
}


