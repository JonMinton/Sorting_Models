#include<cstdio>
#include<Rcpp.h>
#include<vector>

using namespace Rcpp;
using namespace std; 


// This version uses the Rcpp classes only


// [[Rcpp::export]]
void Load(List input_)
{
  // Loading from list
  NumericVector   A = as<NumericVector>   (input_["A"]);
  CharacterVector B = as<CharacterVector> (input_["B"]);
  List            C = as<GenericVector>   (input_["C"]);

  NumericVector   D = as<NumericVector>   (C["D"]);
  NumericVector   E = as<NumericVector>   (C["E"]);
  
  // displaying contents
  
  Rcout << "Contents of A are: ";
  for (int i=0; i < A.size(); ++i)
    Rcout << A[i] << ' ';
  Rcout << endl;
  
  Rcout << "Contents of B are: ";
  for (int i=0; i < B.size(); ++i)
    Rcout << B[i] << ' ';
  Rcout << endl;
  

  
  Rcout << "Contents of D are: ";
  for (int i=0; i < D.size(); ++i)
    Rcout << D[i] << ' ';
  Rcout << endl;

  Rcout << "Contents of E are: ";
  for (int i=0; i < E.size(); ++i)
    Rcout << E[i] << ' ';
  Rcout << endl;


}



class Base{
  int secretnumber;
  std::vector<int> A;
  CharacterVector B;
  List C;
  NumericVector D;
  std::vector<double> E;
  Rcpp::NumericMatrix M;
  
public:
  Base();  // Constructor function  
  void setSecret(int inval) {secretnumber = inval;};
  int getSecret() {return secretnumber;};
  ~Base(){Rcout<< "Destructor called\n";};
  void LoadVec(List input_);
  Rcpp::List ReturnList();
  void LoadMatrix(NumericMatrix M_);
  void DisplayMatrix();
};


Base::Base()
{
  Rcout << "Base Constructed\n";
  secretnumber = 0;
  Rcout << "The secret number is " << secretnumber << std::endl;
};

void Base::LoadMatrix(NumericMatrix M_)
{
  M = M_;
  
  Rcout << "LoadMatrix completed\n";
}

Rcpp::List Base::ReturnList()
{
  return(

    Rcpp::List::create(
      Rcpp::Named("NewB") = B,
      Rcpp::Named("NewD") = D,
      Rcpp::Named("NewM") = M

    )
  );
}
void Base::LoadVec(List input_)
{
  // Loading from list
  A     = Rcpp::as<std::vector<int> > (input_["A"]);
  B     = as<CharacterVector> (        input_["B"]);
  C     = as<GenericVector>   (        input_["C"]);

  D     = as<NumericVector>   (             C["D"]);
  
  E     = Rcpp::as<std::vector<double> > (  C["E"]);
  

  Rcout << "Contents of A are: ";
  for (unsigned int i=0; i < A.size(); ++i)
    Rcout << A[i] << ' ';
  Rcout << endl;
  
  Rcout << "Contents of B are: ";
  for (int i=0; i < B.size(); ++i)
    Rcout << B[i] << ' ';
  Rcout << endl;
  

  
  Rcout << "Contents of D are: ";
  for (int i=0; i < D.size(); ++i)
    Rcout << D[i] << ' ';
  Rcout << endl;

  Rcout << "Contents of E are: ";
  for (unsigned int i=0; i < E.size(); ++i)
    Rcout << E[i] << ' ';
  Rcout << endl;

}


void Base::DisplayMatrix()
{
  int nrow, ncol, i, j;
  
  nrow = M.nrow();
  ncol = M.ncol();
  
  
  Rcout << "The matrix has " << nrow << " rows\n and " << ncol << " columns\n";
  
  for (i = 0; i < nrow; ++i){
    for (j = 0; j < ncol; ++j){
      Rcout << "[" << i << ", "<< j << "]: "<< M(i,j) << endl;    
    }
  }
  
}


RCPP_MODULE(BaseClass){
using namespace Rcpp;
  
  class_<Base> ("Base")
  .constructor()
  .method("get", &Base::getSecret)
  .method("set", &Base::setSecret)
  .method("loadvec", &Base::LoadVec)
  .method("returnlist", &Base::ReturnList)
  .method("loadmatrix", &Base::LoadMatrix)
  .method("displaymatrix", &Base::DisplayMatrix)
  ;
}


