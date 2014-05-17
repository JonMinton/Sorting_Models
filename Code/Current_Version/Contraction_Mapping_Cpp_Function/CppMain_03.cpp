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

  
  std::vector<long double> denom;
  std::vector<long double> rhs; 
  
  long double stayerprop_real;
  long double stayerprop_predicted;
  long double stayerprop_dif;
  
  long double stayp; 
  long double stay;
  long double pstay; 
  long double mu_diff;
  long double mu_guess; 
  long double mu; 

// PARAMETERS
  long double mu_upper, mu_lower;
  long double tol_mu;
  long double tol_delta;
  
  int maxit_outer, maxit_inner;

// DERIVED QUANTITIES  

  long double totpop_t1;
  long double totpop_t2;
  long double dpop; 
  
  int n;
  int n1;
  
    

    

// ESTIMATED QUANTITIES
  std::vector<long double> deltas;

  

  

// INTERMEDIATE QUANTITIES  

  

  

  
// CONTROL QUANTITIES
  
  bool dbgflag;
  bool verboseflag; 
  
  bool deltas_OK;
    bool mu_ok;  
    bool repeat_inner; 
  
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

//Constructors functions

Contraction::Contraction(List In_)
{
    Function browser("browser");    
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
  
  stayp = as<long double> (data["stayerprop"]);
  
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

    n1 = (int) deltas.size();
    n = n1 - 1;
    
 // params objects  

  tol_delta = as<long double> (params["tol_delta"]);
  tol_mu = as<long double> (params["tol_mu"]);
  
  // dbg objects  
  
  dbgflag = as<bool> (dbg["dbgflag"]);
  verboseflag = as<bool> (dbg["verboseflag"]);
  
    x1.resize(n1, 0.0);
    xf.resize(n1, 0.0);
    
    denom.resize(n1, 0.0);    

    newshares.resize(n1);
    newshares_predicted.resize(n1);

    
    util = Rcpp::clone(movingcost);
  
  if (verboseflag){
      Rcout << "Constructor completed" << endl;
  }

}  

void Contraction::TimminsRunContraction()
{
    Function browser("browser");
    timescalled_TimminsRunContraction++;
    Rcout << "TRC\n";
    TimminsPop();
    mu_diff = tol_mu * 99999.9;  
    Rcout << "mu_diff: " << mu_diff << "\ttol_mu: " << tol_mu << "\tstayp: " << stayp << endl;
    browser();
    TimminsOuter();
}


void Contraction::TimminsOuter()
{
    Function browser("browser");    
    Rcout << "TO\n";
    
    timescalled_TimminsOuter++;
   
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
//        Rcout << "200\n";

        for (int j = 0; j < n; j++){
            x1[j] = xf[j];
        }
        

        if (repeat_nums > maxit_outer){
            Rcout << " Repeat_loop false as repeat_nums > maxit_outer\n";
            repeat_loop = false;
        } 
        
        if (abs(mu_diff) < tol_mu){
            Rcout << "Repeat_loop false as abs(mu_diff) < tol_mu\n";
            Rcout << "mu_diff: " << mu_diff << "\ttol_mu: " << tol_mu << endl;
            repeat_loop = false;
        }
        
    }
    
    for (unsigned int j = 0; j < n1; j++){
        deltas[j] = xf[j];
    }
    
}


void Contraction::TimminsOuterTest()
{
    Function browser("browser");
    Rcout << "TOT\n";           
    timescalled_TimminsOuterTest++;
   
/* /// FORTRAN CODE
      DO j = 1,n1
         test(j) = dabs(xf(j)-x1(j))
         IF (test(j).gt.tolc) then
            DO k = 1,n1
               x1(k) = xf(k)
            ENDDO
            GO TO 300
         ENDIF
      ENDDO
      mudiff = stayp - pstay
*/
    while(repeat_inner){
        TimminsInner();
        repeat_inner = false;
        
        for (int j = 0; j < n1; j++){
            long double test = abs(xf[j] - x1[j]);
            if (test > tol_delta){
                for (int k = 0; k < n1; k++){
                    x1[k] = xf[k];
                }
                repeat_inner= true;
            }
        }
        
    }

   
    mu_diff = stayp - pstay;
//    Rcout << "stayp: " << stayp << "\tpstay: " << pstay << "\tmu_diff: " << mu_diff << endl;
//    browser();
}

void Contraction::TimminsInner()
{
    // TO DO
    // 1 ) Look for operations that access and modify xf 
    
    Function browser("browser");     
    Rcout << "TI\n";
    timescalled_TimminsInner++; 
    
    mu = mu_guess;
    

    // At this stage all values are 0.0

//    browser(); 
    
   



    
/*    // FORTRAN CODE
      do j = 1,n1
         do k = 1,n1
            mc = 0.d0
            if ((j.ne.k).and.(j.ne.n1).and.
     &         (k.ne.n1)) mc = (2910.d0+(0.03d0*(1.d0-
     &         rperc(j))*medval(j))+(0.03d0*(1.d0-
     &         rperc(k))*medval(k)))/24.556d0
            if ((j.ne.k).and.((j.eq.n1).or.
     &         (k.eq.n1))) mc = 528.71d0
            util(k,j) = x1(k)-x1(j)-mu*mc
         enddo
      enddo    
*/
// C++ Code    
    for (int j =0; j < n1; j++){
        for (int k = 0; k < n1; k++){
            
//            Rcout << x1[j] << "\t" << x1[j] << "\t" << mu_guess << "\t" << movingcost(j,k) << endl;
//                Rcout << "[ " << j << ", " << k << "]\t" << x1[k] - x1[j] - mu_guess * movingcost(j, k) << endl;

            util(k, j) = x1[k] - x1[j] - mu * movingcost(k, j); // utility values seem OK
        }
    }
    
////////////////////////////////
/* /// FORTRAN CODE

      do j = 1,n1
         denom(j) = 0.d0
         do k = 1,n1
            denom(j) = denom(j)+dexp(util(k,j))
         enddo
      enddo
*/

// C++ code
      
//    Rcout << "260\n";
    for (unsigned int j = 0; j < n1; j++){
        for (unsigned int k = 0; k < n1; k++){
            denom[j] += exp(util(k, j));
        } // denom values seem OK
//        Rcout << "denom[" << j << "]\t" << denom[j] << endl;
    }
/////////////////////////////

/*  /// FORTRAN CODE
      do i = 1,n1
         rhs(i) = 0.d0
         do j = 1,n1
            numer = dexp(util(i,j))
            rhs(i) = rhs(i)+(numer/denom(j))*
     &               pop1(j)
         enddo
      enddo
*/

// C++ CODE
//    Rcout << "267\n";
    std::vector <long double> rhs(n1, 0.0);
    
//    Rcout << "271\n";
     
    for (unsigned int i = 0; i < n1; i++){
        for (unsigned int j = 0; j < n; j++){
            long double numer = exp(util(i, j));
            rhs[i] += (numer/denom[j])* oldcounts[i];
        }
//        Rcout << "rhs[" << i << "]\t" << rhs[i] << endl; 

    }
//////////////////////////////

/*    /// FORTRAN CODE
      stay = 0.d0
      do i = 1,n
         stay = stay + (dexp(util(i,i))/denom(i))*pop1(i)
      enddo
      stay = stay/(totpop-pop2(n1))
*/
// C++ Code
    stay = 0;
    for (unsigned int i = 0; i < n; i++){
        stay += (exp(util(i,i))/denom[i]) * oldcounts[i];
//        Rcout << "stay: " << stay << "\tutil(i,i): " << util(i,i)<< endl;
    }

    stay = stay / (totpop_t2 - dpop);
//    Rcout << "stay after dividing: " << stay << endl;
 
////////////
/*
      pstay = 0.d0
      do i = 1,n
         pstay = pstay + (dexp(util(i,i))/denom(i))/(1.d0*n)
      enddo
*/      
    pstay = 0;

    for (unsigned int i = 0; i < (n-1); i++){
 // pstay appears to be working correctly        
        pstay+= (exp(util(i,i)) / denom[i])/((long double) n);
    }
/////////////////////

/* /// FORTRAN CODE
      do i = 1,n1
         share(i) = pop2(i)/totpop
         pshare(i) = rhs(i)/totpop
      enddo
*/
// C++ Code    
     for (unsigned int i=0; i < n1; i++){
         newshares[i] = newcounts[i] / totpop_t2;
         newshares_predicted[i] = rhs[i] / totpop_t2;         
//         Rcout << "totpop_t2: " << totpop_t2 << "\tnewcounts: " << newcounts[i] << "\trhs: " << rhs[i] << "\tnewshares: " << newshares[i] << "\tnewshares_predicted: " << newshares_predicted[i] <<endl;
     }
//    browser(); 

///////////////
/* // FOTRAN CODE
      xavg = 0.d0
      DO j = 1,n1
         xf(j) = x1(j) + (dlog(share(j)) - dlog(pshare(j)))
c         write(*,*) j,share(j),pshare(j)
         xavg = xavg+(xf(j)/(1.d0*n1))
      ENDDO
 225  format(i10,f15.6,f15.6)
 */

/// C++ Code 
     long double xavg = 0;    
     for (unsigned int j =0; j < n1; ++j){
         xf[j] = x1[j] + (log(newshares[j]) - log(newshares_predicted[j]));
         xavg += xf[j]; 
//         Rcout << "xf: " << xf[j] << "\tx1: " << x1[j] << "\txavg: " << xavg << "\tnewshares: ";
//         Rcout << log(newshares[j]) << "\tnewshares_predicted: " << log(newshares_predicted[j]) << endl;
     }
     xavg = xavg / (long double) n1;
//     Rcout << "final xavg: " << xavg << endl;;
//    browser();
/////////////////////////
/* /// FORTRAN CODE
      temp = xavg
      DO j = 1,n1
         xf(j) = xf(j)-temp
      enddo
*/
/// C++ Code
     for (unsigned int i = 0; i < n1; ++i){
         xf[i] = xf[i] - xavg;
//         Rcout << "final xf[" << i << "] " << xf[i] << endl;
     }
//    browser();
//    Rcout << "336\n";     
}


void Contraction::TimminsPop()
{
    Function browser("browser"); 
    
    Rcout << "TP\n";
    totpop_t1 = 0;
    totpop_t2 = 0;
    
    for (int i = 0; i < n; i++){
        totpop_t1 += oldcounts[i];
        totpop_t2 += newcounts[i];
//      Rcout << "[" << i << "] " << oldcounts[i] << " " << newcounts[i] << "\ttotpop_t1: " << totpop_t1 << "\ttotpop_t2: " << totpop_t2 << endl;
    }

    dpop = totpop_t2 - totpop_t1;

    Rcout << "totpop_t1: " << totpop_t1 << "\ttotpop_t2: " << totpop_t2 << "\tdpop: " << dpop<< endl;
    browser();
    long double tmp = 2.0 * abs(dpop);
    oldcounts.push_back(tmp);
    newcounts.push_back(tmp - dpop);
    

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
     Rcpp::Named("newshares") = newshares,
     Rcpp::Named("trackers") = Rcpp::List::create(
         Rcpp::Named("inner") = timescalled_TimminsInner,
         Rcpp::Named("outer") = timescalled_TimminsOuter,
         Rcpp::Named("outertest") = timescalled_TimminsOuterTest,
         Rcpp::Named("runcontraction") = timescalled_TimminsRunContraction
         )
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


