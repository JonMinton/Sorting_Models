#include<cstdio>
#include<Rcpp.h>
#include<vector>


using namespace Rcpp;
using namespace std; 


class Contraction{

// Private fields 

// ESTIMATED QUANTITIES
  std::vector<long double> deltas;
  std::vector<long double> xf;
  std::vector<long double> x1;
  std::vector<long double> denom;
  std::vector<long double> rhs; 

// DATA  

  NumericMatrix movingcost;
  NumericMatrix util; 
  
  std::vector<long double> oldcounts;
  std::vector<long double> newcounts;  
  std::vector<long double> newshares;
  std::vector<long double> newshares_predicted;
  

  
  
  long double stayp; 
  long double stay;
  long double pstay; 
  long double mu_diff;
  long double mu_guess; 
  long double mu; 
  long double xavg; 

  long double mu_upper, mu_lower;
  long double tol_mu;
  long double tol_delta;
  
  int maxit_outer, maxit_inner;

  long double totpop_t1;
  long double totpop_t2;
  long double dpop; 
  
  int n;
  int n1;
  
  
// CONTROL QUANTITIES
  
  bool dbgflag;
  bool verboseflag; 
  
  bool deltas_OK;
    bool mu_ok;  
    bool repeat_inner; 
    int calls_to_inner;
  
// METHOD TRACKERS  
  int timescalled_TimminsInner;
  int timescalled_TimminsOuter;
  int timescalled_TimminsInnerTest;
    int timescalled_TimminsRunContraction; 
  
public:


  //Constructors
  Contraction();
  Contraction(List In_);
  
  //Destructor
  
  ~Contraction() {Rcout<< "Destructor called\n";};

  void InitialiseTrackers();
  
  void TimminsInner(); // Attempt to replicate Timmins code blindly from 300 onwards
  void TimminsInnerTest(); 
  void TimminsOuter();
  void TimminsRunContraction();
  void TimminsPop();
  
  List ExtractEstimates();
  
};

//Constructors functions

Contraction::Contraction(List In_)
{
/* 
This is the constructor function. It loads in inputs passed to it from R
*/
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
  
  Rcout << "Lengths of vectors oldcounts and newcounts are " << oldcounts.size() << " and " << newcounts.size();
  Rcout << " respectively. Is this correct?" << endl;
  browser();

  
  mu_upper = as<long double>(params["mu_upper"]);
  mu_lower = as<long double>(params["mu_lower"]);
  
  maxit_outer = as<int> (params["maxit_outer"]);
  maxit_inner = as<int> (params["maxit_inner"]);

  // utils objects
    
  deltas = Rcpp::as<std::vector<long double> > (utils["deltas"]);
//  Rcout << deltas.size() << endl; 

    n1 = (int) deltas.size(); // At present the size of the deltas vector passed externally includes the open region
    n = n1 - 1;
    
 // params objects  

  
  tol_delta= as<long double> (params["tol_delta"]);
  tol_mu = as<long double> (params["tol_mu"]);
  
  // dbg objects  
  
  dbgflag = as<bool> (dbg["dbgflag"]);
  verboseflag = as<bool> (dbg["verboseflag"]);
  
    x1.resize(n1, 0.0); // resizes the vector x1 to be n1 elements long, with each element initialised to 0.0
    xf.resize(n1, 0.0);
    
    denom.resize(n1, 0.0);    
    rhs.resize(n1, 0.0);

    newshares.resize(n1);  // Here the vectors are resized, but the elements are not initialised to a specific value
    newshares_predicted.resize(n1);

    
    util = Rcpp::clone(movingcost); // A quick (but maybe computationally expensive) way of creating a second 
    // matrix with the same dimensions as a first matrix. However, care has to be taken to avoid leaving 
    // any of the elements of util unadjusted.
  
  if (verboseflag){
      Rcout << "Constructor completed" << endl;
  }

}  

void Contraction::TimminsRunContraction()
{
/* 
    High level function for controlling main stages of the contraction mapping
*/    
    Function browser("browser");
    timescalled_TimminsRunContraction++; // tracks number of times this function has been called
    Rcout << "TRC\n";
    TimminsPop(); //Run routines which calculate counts for open region (n1st element)
    mu_diff = tol_mu * 99999.9;  // set mu_diff to a value big enough that it exceeds mu_tol and the outer loop will 
    // not complete on first iteration
    Rcout << "mu_diff: " << mu_diff << "\ttol_mu: " << tol_mu << "\tstayp: " << stayp << endl;
    browser();// Allows code execution to be stopped until user presses return
    TimminsOuter(); // Runs main part of contraction mapping algorithm
}


void Contraction::TimminsOuter()
{
/* 
Cascades from TimminsRunContraction
*/
    Function browser("browser");    
    timescalled_TimminsOuter++;
    Rcout << "TO " << timescalled_TimminsOuter << endl;
    
    bool repeat_loop = true; 
    unsigned int repeat_nums = 0; 
    // Using a while loop instead of a go to statement
    while (repeat_loop){
        repeat_nums++;

        mu_guess = (mu_upper + mu_lower) / 2.0;
        Rcout << "Outer: " << repeat_nums << "\tmu: " << mu_guess << endl;
        TimminsInnerTest(); // Similar to but not exactly the same as fcn1
        // in particular the variables x1 and xf are global and so changing them
        // anywhere will change them in the object
         
         // bifurcation
        if (mu_diff > 0.0){
            mu_lower = mu_guess;
        } else {
            mu_upper = mu_guess;
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
    
    // Load the vector deltas with xf once loop has completed as they are the best guess
    for (unsigned int j = 0; j < n1; j++){
        deltas[j] = xf[j];
    }
    
}


void Contraction::TimminsInnerTest()
{
    Function browser("browser");
    Rcout << "TOT\n";           
    timescalled_TimminsInnerTest++;
   
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
*/ // I NOW THINK THE COMPARISON SHOULD BE BETWEEN STAY AND STAYP AND NOT STAYP AND PSTAY
    calls_to_inner = 0; 
    repeat_inner = true;
    while(repeat_inner){
        TimminsInner();
        calls_to_inner++;
        repeat_inner = false;
        bool test_x = true;
        int j = 0;
        while (test_x){
            long double test = abs(xf[j] - x1[j]);
            j++;
            if (test > tol_delta){
                repeat_inner = true;
                test_x = false;
            }
        }
        if (calls_to_inner > maxit_inner){ 
            repeat_inner = false;
        }
    }
/*
 To be discussed: calculation and interpretion of pstay and stay
 */
 
//    mu_diff = stayp - pstay;
    mu_diff = stayp - stay; 
//    Rcout << "stayp: " << stayp << "\tstay: " << stay << "\tmu_diff: " << mu_diff << endl;
//    browser();
}

void Contraction::TimminsInner()
{

    
    Function browser("browser");     
//    Rcout << "TI\n";
    timescalled_TimminsInner++; 
    
    mu = mu_guess;
    
                
    for (int k = 0; k < n1; k++){
        x1[k] = xf[k];
    }
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
/* 
NOTE : THis is different as moving costs are calculated in R then passed to C++ as a matrix
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
    for (int j = 0; j < n1; j++){
        denom[j] = 0.0;
        for (int k = 0; k < n1; k++){  
            denom[j] += exp(util(k, j));
        } 
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

/* 
Note: the var oldcounts is used instead of pop1
oldcounts is the vector of counts (with correction) in the first period
*/
// C++ CODE
     
    for (int i = 0; i < n1; i++){
        rhs[i] = 0.0;
        for (int j = 0; j < n1; j++){
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

/* totpop here is assumed to be the total number of counts observed in time period 2
Is this correct?
*/
// CHECK WITH CHRIS THAT MY INTERPRETATION OF TOTPOP IS CORRECT
// C++ Code
    stay = 0;
    for (int i = 0; i < n; i++){
        stay += (exp(util(i,i))/denom[i]) * oldcounts[i];
//        Rcout << "stay: " << stay << "\tutil(i,i): " << util(i,i)<< endl;
    }


    stay = stay / (totpop_t2 - dpop);
//    Rcout << "stay after dividing: " << stay << endl;
//    browser(); 


////////////
/* // FORTRAN CODE
      pstay = 0.d0
      do i = 1,n
         pstay = pstay + (dexp(util(i,i))/denom(i))/(1.d0*n)
      enddo
*/      
    pstay = 0;
    for (int i = 0; i < n; i++){
        pstay+= (exp(util(i,i)) / denom[i])/((long double) n);
    }

/* /// FORTRAN CODE
      do i = 1,n1
         share(i) = pop2(i)/totpop
         pshare(i) = rhs(i)/totpop
      enddo
*/
long double tmp = 0.0;
// C++ Code    
     for (int i=0; i < n1; i++){
         newshares[i] = newcounts[i] / totpop_t2;
         newshares_predicted[i] = rhs[i] / totpop_t2;  
         tmp += pow(newshares_predicted[i] - newshares[i], 2);
//         Rcout << "totpop_t2: " << totpop_t2 << "\tnewcounts: " << newcounts[i] << "\trhs: " << rhs[i] << "\tnewshares: " << newshares[i] << "\tnewshares_predicted: " << newshares_predicted[i] <<endl;
     }
    tmp = tmp / (long double) newshares.size();
    tmp = sqrt(tmp);
    if ((calls_to_inner % 100) == 0){
        Rcout << tmp << endl;        
    }


///////////////

/* 
The contraction mapping itself
*/

/* // FORTRAN CODE
      xavg = 0.d0
      DO j = 1,n1
         xf(j) = x1(j) + (dlog(share(j)) - dlog(pshare(j)))
c         write(*,*) j,share(j),pshare(j)
         xavg = xavg+(xf(j)/(1.d0*n1))
      ENDDO
 225  format(i10,f15.6,f15.6)
 */

/// C++ Code 
     xavg = 0;    
     for (int j =0; j < n1; ++j){
         xf[j] = x1[j] + (log(newshares[j]) - log(newshares_predicted[j]));
         xavg += xf[j] /((long double) n1);  
     }
 
/////////////////////////
/* 
Subtration of averages before passing back to mapping procedure again
*/ 

/* /// FORTRAN CODE
      temp = xavg
      DO j = 1,n1
         xf(j) = xf(j)-temp
      enddo
*/
/// C++ Code
     for (int i = 0; i < n1; i++){
         xf[i] = xf[i] - xavg;
//         Rcout << "final xf[" << i << "] " << xf[i] << endl;
     }
     
//     Rcout << xavg << endl;
//    browser();
//    Rcout << "336\n";     
}


void Contraction::TimminsPop()
{
    Function browser("browser"); 
    /* // FORTRAN CODE
    
          do i = 1,n
         pop1t(i) = white_l00(i)+0.001d0
         pop2t(i) = white_l09(i)+0.001d0
         pop00 = pop00+pop1t(i)
         pop09 = pop09+pop2t(i)
      enddo
 555  format(i9,3f21.7)
      dpop = pop09-pop00
      write(*,*) 'dpop =',dpop
      do i = 1,n
         pop1(i) = pop1t(i)
         pop2(i) = pop2t(i)
      enddo
      pop1(n1) = 2.d0*dabs(dpop)
      pop2(n1) = pop1(n1)-dpop
     */
      
    Rcout << "TP\n";
    // C++ CODE
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
    
    // the method Y.push_back(x) adds an element with the value x to the end of the vector Y
    // As oldcounts and newcounts are the data themselves, they are vectors of length n
    // adding an element to the end therefore changes them to vectors of length n1
    oldcounts.push_back(tmp);
    newcounts.push_back(tmp - dpop);
    
    Rcout << "The total lengths of oldcounts and newcounts are " << oldcounts.size() << " and " <<  newcounts.size() << endl;
    Rcout << "They should be 1922 in the current example. Are they?" << endl;
    browser();
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
  timescalled_TimminsInnerTest = 0;
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
         Rcpp::Named("innertest") = timescalled_TimminsInnerTest,
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


