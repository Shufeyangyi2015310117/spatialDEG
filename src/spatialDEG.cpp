// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rmath.h>
#include "mt_paral_testing_job.hpp"

// [[Rcpp::plugins(cpp11)]]

using namespace arma;  
using namespace Rcpp;



// get the length-scale/periodicity parameter
void get_kscale_vec6( mat& spa, vec& kscale_vec, mat& sqdist, int num_ls) {
    int N=spa.n_rows; 
    int n_pair = N*(N-1)/2;
    vec pair_dist(n_pair, fill::zeros); 
    int counter=0;
    double m1=0, m2=0; 

    for (int i=0; i<N-1; i++) {
        for (int j=i+1; j<N; j++) {
            sqdist(i,j) = pow( spa(i,0)-spa(j,0), 2 ) + pow( spa(i,1)-spa(j,1), 2 );
            //dist(i,j) = sqrt( sqdist(i,j) );

            pair_dist(counter) = sqrt( sqdist(i,j) );
            counter++; 
        }
    }
    
    vec p = {0, 1}; 
    vec mvec = quantile(pair_dist, p);
    m1 = mvec(0);
    m2 = mvec(1);
    //m1 = min( pair_dist );
    //m2 = max( pair_dist );
    Rprintf("m1 %f m2 %f \n", m1, m2);
    
    double lm1 = log10(m1); //m1 = min( pair_dist ); lm1=log10(0.5*m1)
    double lm2 = log10(m2); //m2 = max( pair_dist ); lm2=log10(2.0*m2)
    double ldiff = lm2 - lm1; 
    double diff = m2 - m1;
    double temp = 0.0;
    for (int i=0; i< (num_ls/2); i++) {
        //temp = lm1 + i/9.0*ldiff; 
        //kscale_vec(i) = pow(10.0, temp);
        kscale_vec(i) = m1 + i/((double)(num_ls/2) - 1.0)*diff;
        kscale_vec(i+(num_ls/2)) = kscale_vec(i);
    }

}



// get gaussian kernel matrix K
void get_gauss_kern( mat& spa, double kscale, mat& sqdist, mat& K ) {
    int N=spa.n_rows; 
    double denom = 2.0*pow(kscale,2); 

    K.diag().ones();

    for (int i=0; i<N-1; i++) {
        for (int j=i+1; j<N; j++) {
            K(i,j) = exp( - sqdist(i,j)/denom );
            K(j,i) = K(i,j);
        }
    }
}

// get cosine kernel matrix K
void get_cosine_kern( mat& spa, double kscale, mat& sqdist, mat& K ) {
    int N=spa.n_rows; 
    double pi = 3.141592653589793238463; 
    double coef = 2.0*pi/kscale;
    K.diag().ones();
    
    for (int i=0; i<N-1; i++) {
        for (int j=i+1; j<N; j++) {
            K(i,j) = cos( sqrt(sqdist(i,j)) * coef );
            K(j,i) = K(i,j);
        }
    }    
}

void get_cosine_kern2( mat& spa, double kscale, mat& sqdist, mat& K ) {
    int N=spa.n_rows; 
    double pi = 3.141592653589793238463; 
    double coef = 2.0*pi/kscale;
    K.diag().ones();
    
    double temp=0.0;
    for (int i=0; i<N-1; i++) {
        for (int j=i+1; j<N; j++) {
            temp = pow( cos( sqrt(sqdist(i,j)) * coef ), 2.0 ); 
            K(i,j) = exp(- temp );
            K(j,i) = K(i,j);
        }
    }    
}




//' @author Yi Yang
//' @title
//' SysMatEigen2
//' @description
//' SysMatEigen2 to transform the matrix to a positive definite matrix
//' @param M  a matrix
//' @return a List
//'
//'
//' @details
//' \code{spatialDEG_test} run spatialDEG model.
//' @export
// [[Rcpp::export]]
List SysMatEigen2(arma::mat M) {
try {
		arma::vec eigval = zeros<vec>( M.n_rows );
		arma::mat eigvec = zeros<mat>( size(M) );
		eig_sym(eigval, eigvec, M, "dc");
		const uvec idx = find(eigval < 1e-8 );
		arma::vec tmp_value = ones<arma::vec>(idx.n_elem);
		eigval.elem( idx ) = tmp_value * 1e-8;
		arma::mat M1 = eigvec.each_row() % eigval.t();
		M = M1 * eigvec.t();
		// return values
		return List::create(Named("eigval") = eigval, Named("eigvec") = eigvec, Named("kernel_mat") = M);
		//return List::create(Named("eigval") = eigval, Named("eigvec") = eigvec);
	}
	catch (std::exception &ex)
	{
		forward_exception_to_r(ex);
	}
	catch (...)
	{
		::Rf_error("C++ exception (unknown reason)...");
	}
	return R_NilValue;
}




//' @author Yi Yang
//' @title
//' spatialDEG
//' @description
//' spatialDEG to perform differentially expressed gene analysis by leveraging spatial information with true kernels.
//' @param W  covariates for the two datasets
//' @param Y  stacked gene expression
//' @param K1  kernel for the first dataset.
//' @param K2  kernel for the second dataset.
//' @param Initial_theta  initial value of theta
//' @param mu_fixed logical value specifying whether mu is fixed
//' @param verbose logical value specifying whether to print the parameters during the spatialDEG analysis or not
//' @param check_positive a logical value specifying whether to check the positive definiteness of kernels
//'
//' @return List of model parameters
//'
//'
//' @details
//' \code{spatialDEG_test} run spatialDEG model.
//' @export
// [[Rcpp::export]]
List spatialDEG_true_kernel_test(arma::mat& W, arma::mat& Y, arma::mat& K1, arma::mat& K2, arma::vec& Initial_theta, bool mu_fixed = true, bool verbose = true, bool check_positive = false) {
  // input: gene expression Y = [y1 y2 ... y_nphe] at N pixels
  // spa is an N by 2 matrix containing the (x,y) coordinates of the N pixels
  // if number of pixels N is too large, do not store sqdist and dist -- recalculate dist for every length-scale parameter
  // Initial_theta 
    
  // eigen decomposion of Kernel
    
  
  int n1 = K1.n_rows; 
  int n2 = K2.n_rows;
    
  int p = W.n_cols;
     
    vec D1;
    mat K1_vector;
    
    if (check_positive == true){
        List res = SysMatEigen2(K1);
        D1 = as<vec>(res["eigval"]);
        K1_vector = as<mat>(res["eigvec"]);
    }else{
        eig_sym(D1, K1_vector, K1);
    }
    
    cout << "D1 " << D1 << endl;
    
  vec D2;
  mat K2_vector;
    
    if (check_positive == true){
        List res2 = SysMatEigen2(K2);
        D2 = as<vec>(res2["eigval"]);
        K2_vector = as<mat>(res2["eigvec"]);
    }else{
        eig_sym(D2, K2_vector, K2);
    }
    
        cout << "D2 " << D2 << endl;
    
  int ngene = Y.n_cols;
  int Ncell = Y.n_rows;
    
  mat I1 = eye(n1, n1);   // identity matrix
  mat I2 = eye(n2, n2);
    
  mat output = zeros(ngene, 13);
  int count = 0;
  double tol = 0.0001;
    
  vec z = join_cols(zeros(n1, 1), ones(n2, 1));
      
  mat K1_vector_zeros = join_rows(K1_vector, zeros(n1, n2));
  mat zeros_K2_vector = join_rows(zeros(n2, n1), K2_vector);
  mat U = join_cols(K1_vector_zeros, zeros_K2_vector);
    
  mat Ytilde = U.t()*Y;
  mat Wtilde = U.t()*W;
  vec ztilde = U.t()*z;
    
  vec H, S1, S2;
  vec beta = zeros(p, 1);
  mat beta_all = zeros(p, ngene);
  double mu = 0;
  List ret;
  double diff;
  vec etilde;
    
  for (int i = 0; i<ngene; i++){

      colvec theta = zeros(3, 1);
      theta = Initial_theta; // sigma_e, sigma_1, sigma_2
      colvec theta_old = theta;
      
      
      double sigma_e = 0;
      sigma_e = theta(0,0);
      
      vec sigma = zeros(2, 1); 
      sigma(0,0) = theta(1,0);
      sigma(1,0) = theta(2,0);
      
      vec Omega_11, Omega_22;
      vec H_11, H_22, H;
      vec S1, S2;
      vec score = zeros(3, 1);
      mat AImat = zeros(3, 3);  //average information matrix
        

      if (verbose){
        cout <<  "Run Newton-Raphson algorithm" << endl;
      }

      int iter = 0;
      diff = datum::inf;
      while(diff > tol && iter < 100){

        Omega_11 = sigma(0,0)*D1 + sigma_e*diagvec(I1);  // N-by-1 vector
        Omega_22 = sigma(1,0)*D2 + sigma_e*diagvec(I2);  // N-by-1 vector

        if (verbose){
            cout <<  "calculate inversion of Omega" << endl;
        }

        H_11 = 1.0/Omega_11;
        H_22 = 1.0/Omega_22;
          
        if (verbose){
            cout <<  "calculate H" << endl;
        }
          
        H = join_cols(H_11, H_22);
        S1 = join_cols(D1, zeros(n2, 1));
        S2 = join_cols(zeros(n1, 1), D2);
          
        if (verbose){
            cout <<  "calculate etilde" << endl;
        }
          
        etilde = Ytilde.col(i) - Wtilde*beta - ztilde*mu; 
          
        if (verbose){
            cout <<  "calculate score" << endl;
        }
            
        score(0,0) = 0.5*sum(etilde%H%H%etilde) - 0.5*sum(H);
        score(1,0) = 0.5*sum(etilde%H%S1%H%etilde) - 0.5*sum(H%S1);
        score(2,0) = 0.5*sum(etilde%H%S2%H%etilde) - 0.5*sum(H%S2);
                   
        if (verbose){
            cout << "score is " << score << endl;
            cout <<  "calculate AI matrix" << endl;
        }
   
        AImat(0,0) = sum(etilde%H%H%H%etilde);
        AImat(1,1) = sum(etilde%H%S1%H%S1%H%etilde);
        AImat(2,2) = sum(etilde%H%S2%H%S2%H%etilde);

        AImat(0,1) = sum(etilde%H%H%S1%H%etilde);
        AImat(1,0) = AImat(0,1);
        AImat(0,2) = sum(etilde%H%H%S2%H%etilde);
        AImat(2,0) = AImat(0,2);

        AImat(1,2) = sum(etilde%H%S1%H%S2%H%etilde);
        AImat(2,1) = AImat(1,2);
          
          
        vec Dtheta =  solve(0.5*AImat, score);
        if (verbose){
            cout << "Dtheta is " << Dtheta << endl;
        }  
        theta = theta + Dtheta;
          
        beta = inv_sympd(Wtilde.t()%repmat(H, 1, p).t()*Wtilde)*(Wtilde.t()%repmat(H, 1, p).t())*(Ytilde.col(i) - ztilde*mu);    
          
        if (mu_fixed == false){
           mu = 1.0/sum(ztilde%H%ztilde)*sum(ztilde.t()%H.t()*(Ytilde.col(i) - Wtilde*beta)); 
        }
  
        if (verbose){
            cout << "AImat is " << AImat << endl;
            cout <<  "update theta" << endl;
        }
          
        if (theta(0)<0.001){
            theta(0) = 0.001;
        }
        if (theta(1)<0.001){
            theta(1) = 0.001;
        }
        if (theta(2)<0.001){
            theta(2) = 0.001;
        }

        if (verbose){
            cout << "theta is " << theta << endl;
        }         
          
        sigma_e = theta(0);        
        sigma(0,0) = theta(1);
        sigma(1,0) = theta(2);
          
        diff = norm(theta - theta_old);
        theta_old = theta;
        if (verbose){
            cout <<  "diff is " << diff << endl;
        }    
        iter++;
    
      }

        double loglik = - 0.5*sum(log(1.0/H)) - 0.5*sum(etilde%H%etilde);
  
        output(count, 0) = loglik;
        output(count, 1) = theta(0);
        output(count, 2) = theta(1);
        output(count, 3) = theta(2);
        output(count, 4) = mu;
        beta_all.col(count) = beta;        

        count++;
        if (verbose){
            cout <<  "count is " <<  count << endl;
        }

  }
    
    
  ret["output"] = output;
  ret["beta"] = beta_all;

  // construct return object
  
  
  return ret;
  
}




//' @author Yi Yang
//' @title
//' spatialDEG_paral_test
//' @description
//' spatialDEG to perform differentially expressed gene analysis by leveraging spatial information in parallel.
//'
//' @param spa1  spatial location for the first dataset.
//' @param spa2  spatial location for the second dataset.
//' @param W  covariates for the two datasets
//' @param Y  stacked gene expression
//' @param Initial_theta  initial value of theta
//' @param num_ls total number of length scale and periodicity
//' @param iteration The maximum iteration
//' @param Kernel_fixd a logical value specifying whether the kernel is fixed
//' @param verbose a logical value specifying whether to print the parameters during the spatialDEG analysis or not
//' @param kernel_matched a logical value specifying whether the the kernels under both null hypothesis and alternative hypothesis are restricted to be the same type of kernel
//' @param coreNum Number of cores used for parallel computation
//' @param tol tolence to stop the algorithm
//'
//' @return List of model parameters
//'
//'
//' @details
//' \code{spatialDEG_paral_test} run spatialDEG model.
//' @export
// [[Rcpp::export]]
List spatialDEG_paral_test(arma::mat& spa1, arma::mat& spa2, arma::mat& W, arma::mat& Y, arma::vec& Initial_theta, int num_ls = 20, int iteration = 100, bool Kernel_fixd = false, bool verbose = false, bool kernel_matched = false, const int coreNum = 1, double tol = 0.0001) {
  // input: gene expression Y = [y1 y2 ... y_nphe] at N pixels
  // spa is an N by 2 matrix containing the (x,y) coordinates of the N pixels
  // if number of pixels N is too large, do not store sqdist and dist -- recalculate dist for every length-scale parameter
  // Initial_theta 


  // get the number of cells for each dataset, genes and covariates
  int n1 = spa1.n_rows; 
  int n2 = spa2.n_rows;
  int p = W.n_cols;
  int ngene = Y.n_cols;
  int Ncell = Y.n_rows;
    
  // define 
  mat sqdist1(n1,n1, fill::zeros), K1(n1,n1, fill::zeros); 
  mat sqdist2(n2,n2, fill::zeros), K2(n2,n2, fill::zeros);
  vec kscale1_vec(num_ls, fill::zeros); 
  vec kscale2_vec(num_ls, fill::zeros);  

    
  // calculate all the distances of pair cells i and j and and define 10 lengthscale
  get_kscale_vec6(spa1, kscale1_vec, sqdist1, num_ls);
  get_kscale_vec6(spa2, kscale2_vec, sqdist2, num_ls);
    
    
  cout << "The " << num_ls << " lengthscales for the first dataset are " << kscale1_vec << endl;
  cout << "The " << num_ls << " lengthscales for the second dataset are " << kscale2_vec << endl;
    
  // define the iterative variable for lengthscale and initialize it
  double kscale1 = kscale1_vec(0);
  double kscale2 = kscale2_vec(0);
    
  // define the matrix for storing all the eigen values, rows for cells and columns for kernels    
  mat D1_all(n1, num_ls);
  mat D2_all(n2, num_ls);
    
  // define cubes for storing the pseudo gene expression and covariates and matrix for storing the pseudo indicator
  cube Ytilde_all(Ncell, ngene, num_ls, fill::zeros);
  cube Wtilde_all(Ncell, p, num_ls, fill::zeros );
  mat ztilde_all(Ncell, num_ls);
   
  // 
  cout << "Get the pseudo gene expression, covariates and indicators" << endl;
    

    for (int k = 0; k < num_ls/2; k++){

        kscale1 = kscale1_vec(k);
        kscale2 = kscale2_vec(k);

        get_gauss_kern(spa1, kscale1, sqdist1, K1); 
        get_gauss_kern(spa2, kscale2, sqdist2, K2); 

        vec D1_tmp = zeros(n1, 1);
        vec D2_tmp = zeros(n2, 1);
        mat K1_vector_tmp = zeros(n1, n1);
        mat K2_vector_tmp = zeros(n2, n2);

        eig_sym(D1_tmp, K1_vector_tmp, K1);
        eig_sym(D2_tmp, K2_vector_tmp, K2);

        D1_all.col(k) = D1_tmp;
        D2_all.col(k) = D2_tmp;

        vec z = join_cols(zeros(n1, 1), ones(n2, 1));
        mat K1_vector_zeros = join_rows(K1_vector_tmp, zeros(n1, n2));
        mat zeros_K2_vector = join_rows(zeros(n2, n1), K2_vector_tmp);
        mat U = join_cols(K1_vector_zeros, zeros_K2_vector);

        Ytilde_all.slice(k) = U.t()*Y;
        Wtilde_all.slice(k) = U.t()*W;
        ztilde_all.col(k) = U.t()*z;
    }
    
    
    for (int k = num_ls/2; k < num_ls; k++){

        kscale1 = kscale1_vec(k);
        kscale2 = kscale2_vec(k);

        get_cosine_kern(spa1, kscale1, sqdist1, K1); 
        get_cosine_kern(spa2, kscale2, sqdist2, K2);         
        
        List res = SysMatEigen2(K1);
        mat K1_vector_tmp = as<mat>(res["eigvec"]);
        vec D1_tmp = as<vec>(res["eigval"]);
        
        List res2 = SysMatEigen2(K2);
        mat K2_vector_tmp = as<mat>(res2["eigvec"]);
        vec D2_tmp = as<vec>(res2["eigval"]);

        //eig_sym(D1_tmp, K1_vector_tmp, K1);
        //eig_sym(D2_tmp, K2_vector_tmp, K2);    

        D1_all.col(k) = D1_tmp;
        D2_all.col(k) = D2_tmp;

        vec z = join_cols(zeros(n1, 1), ones(n2, 1));
        mat K1_vector_zeros = join_rows(K1_vector_tmp, zeros(n1, n2));
        mat zeros_K2_vector = join_rows(zeros(n2, n1), K2_vector_tmp);
        mat U = join_cols(K1_vector_zeros, zeros_K2_vector);

        Ytilde_all.slice(k) = U.t()*Y;
        Wtilde_all.slice(k) = U.t()*W;
        ztilde_all.col(k) = U.t()*z;
    }
    
            
   cout << "### Start running SpDEG on genes ... " << endl;;

   mat out_param = -99 * ones<mat>(ngene, 20);

   //set parallel structure object
   parGene_SpDEG parObj(D1_all, D2_all, Ytilde_all, Wtilde_all, ztilde_all, 
                        Initial_theta, kscale1_vec, kscale2_vec, 
                        num_ls, n1, n2, p, ngene, iteration, 
                        Kernel_fixd, verbose, kernel_matched, out_param, tol);

   const int n_thread = coreNum;
   std::vector<std::thread> threads(n_thread);
   for (int i_thread = 0; i_thread < n_thread; i_thread++){
      threads[i_thread] = std::thread(&parGene_SpDEG::update_by_thread_SpDEG, &parObj, i_thread);
   }

   for (int i = 0; i < n_thread; i++){
      threads[i].join();
   }
   out_param = parObj.out_param;
    
  List ret;
    
  ret["out_param"] = out_param;
  ret["kscale1"] = kscale1_vec;
  ret["kscale2"] = kscale2_vec;

  // construct return object
  
  
  return ret;
  
}








//' @author Yi Yang
//' @title
//' spatialDEG
//' @description
//' spatialDEG to perform differentially expressed gene analysis by leveraging spatial information.
//'
//' @param spa1  spatial location for the first dataset.
//' @param spa2  spatial location for the second dataset.
//' @param W  covariates for the two datasets
//' @param Y  stacked gene expression
//' @param Initial_theta  initial value of theta
//' @param num_ls total number of length scale and periodicity
//' @param mu_fixed logical value specifying whether mu is fixed
//' @param verbose logical value specifying whether to print the parameters during the spatialDEG analysis or not
//'
//' @return List of model parameters
//'
//'
//' @details
//' \code{spatialDEG_test} run spatialDEG model.
//' @export
// [[Rcpp::export]]
List spatialDEG_test(arma::mat& spa1, arma::mat& spa2, arma::mat& W, arma::mat& Y, arma::vec& Initial_theta, int num_ls = 20, bool mu_fixed = true, bool verbose = true) {
  // input: gene expression Y = [y1 y2 ... y_nphe] at N pixels
  // spa is an N by 2 matrix containing the (x,y) coordinates of the N pixels
  // if number of pixels N is too large, do not store sqdist and dist -- recalculate dist for every length-scale parameter
  // Initial_theta 


  // get the number of cells for each dataset, genes and covariates
  int n1 = spa1.n_rows; 
  int n2 = spa2.n_rows;
  int p = W.n_cols;
  int ngene = Y.n_cols;
  int Ncell = Y.n_rows;
    
  // define 
  mat sqdist1(n1,n1, fill::zeros), K1(n1,n1, fill::zeros); 
  mat sqdist2(n2,n2, fill::zeros), K2(n2,n2, fill::zeros);
  vec kscale1_vec(num_ls, fill::zeros); 
  vec kscale2_vec(num_ls, fill::zeros);  

    
  // calculate all the distances of pair cells i and j and and define 10 lengthscale
  get_kscale_vec6(spa1, kscale1_vec, sqdist1, num_ls);
  get_kscale_vec6(spa2, kscale2_vec, sqdist2, num_ls);
  cout << "The 10 lengthscale for the first dataset are " << kscale1_vec << endl;
  cout << "The 10 lengthscale for the second dataset are " << kscale2_vec << endl;
    
  // define the iterative variable for lengthscale and initialize it
  double kscale1 = kscale1_vec(0);
  double kscale2 = kscale2_vec(0);
    
  // define the matrix for storing all the eigen values, rows for cells and columns for kernels    
  mat D1_all(n1, num_ls);
  mat D2_all(n2, num_ls);
    
  // define cubes for storing the pseudo gene expression and covariates and matrix for storing the pseudo indicator
  cube Ytilde_all(Ncell, ngene, num_ls, fill::zeros);
  cube Wtilde_all(Ncell, p, num_ls, fill::zeros );
  mat ztilde_all(Ncell, num_ls);
   
  // 
  cout << "Get the pseudo gene expression, covariates and indicators" << endl;
    for (int k = 0; k < num_ls/2; k++){
        
        kscale1 = kscale1_vec(k);
        kscale2 = kscale2_vec(k);
        
        get_gauss_kern(spa1, kscale1, sqdist1, K1); 
        get_gauss_kern(spa2, kscale2, sqdist2, K2); 
        
        vec D1_tmp = zeros(n1, 1);
        vec D2_tmp = zeros(n2, 1);
        mat K1_vector_tmp = zeros(n1, n1);
        mat K2_vector_tmp = zeros(n2, n2);
        
        eig_sym(D1_tmp, K1_vector_tmp, K1);
        eig_sym(D2_tmp, K2_vector_tmp, K2);
        
        D1_all.col(k) = D1_tmp;
        D2_all.col(k) = D2_tmp;
        
        vec z = join_cols(zeros(n1, 1), ones(n2, 1));
        mat K1_vector_zeros = join_rows(K1_vector_tmp, zeros(n1, n2));
        mat zeros_K2_vector = join_rows(zeros(n2, n1), K2_vector_tmp);
        mat U = join_cols(K1_vector_zeros, zeros_K2_vector);
        
        Ytilde_all.slice(k) = U.t()*Y;
        Wtilde_all.slice(k) = U.t()*W;
        ztilde_all.col(k) = U.t()*z;
    }
    
    
    
    for (int k = num_ls/2; k < num_ls; k++){
        
        kscale1 = kscale1_vec(k);
        kscale2 = kscale2_vec(k);
        
        get_cosine_kern(spa1, kscale1, sqdist1, K1); 
        get_cosine_kern(spa2, kscale2, sqdist2, K2); 
        
        vec D1_tmp = zeros(n1, 1);
        vec D2_tmp = zeros(n2, 1);
        mat K1_vector_tmp = zeros(n1, n1);
        mat K2_vector_tmp = zeros(n2, n2);
        
        eig_sym(D1_tmp, K1_vector_tmp, K1);
        eig_sym(D2_tmp, K2_vector_tmp, K2);
        
        D1_all.col(k) = D1_tmp;
        D2_all.col(k) = D2_tmp;
        
        vec z = join_cols(zeros(n1, 1), ones(n2, 1));
        mat K1_vector_zeros = join_rows(K1_vector_tmp, zeros(n1, n2));
        mat zeros_K2_vector = join_rows(zeros(n2, n1), K2_vector_tmp);
        mat U = join_cols(K1_vector_zeros, zeros_K2_vector);
        
        Ytilde_all.slice(k) = U.t()*Y;
        Wtilde_all.slice(k) = U.t()*W;
        ztilde_all.col(k) = U.t()*z;
    }

    
  // define the iterative vector of eigen values and constant vector of ones
  vec D1 = D1_all.col(0);
  vec D2 = D2_all.col(0);
  vec I1 = ones(n1, 1);  
  vec I2 = ones(n2, 1);

  // define the iterative matrices for pseudo gene expression and pseudo covariates and vector for pseudo indicators and initilize them
  mat Ytilde = Ytilde_all.slice(0);
  mat Wtilde = Wtilde_all.slice(0);
  vec ztilde = ztilde_all.col(0);
    
  // define the output list 
  List ret;
  mat output = zeros(ngene, 13);
  mat beta_all = zeros(p, ngene);
    
  int count = 0;
  double tol = 0.0001;
  double diff;

  vec H, S1, S2, H_square;
  vec etilde, etilde_square;
    
  vec etilde1, etilde2, etilde1_square, etilde2_square;
  uword idx1, idx2;
    
  double loglik; 
  vec loglik1_all, loglik2_all;
    
  cout <<  "Perform the spatial DEGs algorithm for all genes" <<  endl;
  for (int i = 0; i<ngene; i++){

      // define the model parameters
      colvec theta = zeros(3, 1);
      theta = Initial_theta; // sigma_e, sigma_1, sigma_2
      colvec theta_old = theta;
      vec beta = zeros(p, 1);
      double mu = 0;      
      
      vec Omega_11, Omega_22;
      vec H_11, H_22;
      
      // define the score and AI matrix
      vec score = zeros(3, 1);
      mat AImat = zeros(3, 3);  
        

      if (verbose){
        cout <<  "Run Newton-Raphson algorithm" << endl;
      }

      int iter = 0;
      diff = datum::inf;
      while(diff > tol && iter < 100){

        // get the diagonal elements of covariance matrix
        Omega_11 = theta(1)*D1 + theta(0)*I1;  // N-by-1 vector
        Omega_22 = theta(2)*D2 + theta(0)*I2;  // N-by-1 vector

        if (verbose){
            cout <<  "calculate inversion of Omega" << endl;
        }
          
        // get the diagonal elements of precision matrix
        H_11 = 1.0/Omega_11;
        H_22 = 1.0/Omega_22;
          
        if (verbose){
            cout <<  "calculate H" << endl;
        }
        
        // calculate the predefined quantities used in the AI algorithm
        H = join_cols(H_11, H_22);
        S1 = join_cols(D1, zeros(n2, 1));
        S2 = join_cols(zeros(n1, 1), D2);
          
        if (verbose){
            cout <<  "calculate etilde" << endl;
        }
          
        // estimate the error
        etilde = Ytilde.col(i) - Wtilde*beta - ztilde*mu; 
          
        if (verbose){
            cout <<  "calculate score" << endl;
        }
            
          
        etilde_square = etilde%etilde;
        H_square = H%H;
          
        // calculate the score
        score(0,0) = 0.5*sum(etilde_square%H_square) - 0.5*sum(H);
        score(1,0) = 0.5*sum(etilde_square%H_square%S1) - 0.5*sum(H%S1);
        score(2,0) = 0.5*sum(etilde_square%H_square%S2) - 0.5*sum(H%S2);
                   
        if (verbose){
            cout << "score is " << score << endl;
            cout <<  "calculate AI matrix" << endl;
        }
   
        // calculate the AI matrix
        AImat(0,0) = sum(etilde_square%H_square%H);
        AImat(1,1) = sum(etilde_square%H_square%S1%H%S1);
        AImat(2,2) = sum(etilde_square%H_square%S2%H%S2);

        AImat(0,1) = sum(etilde_square%H_square%S1%H);
        AImat(1,0) = AImat(0,1);
        AImat(0,2) = sum(etilde_square%H_square%S2%H);
        AImat(2,0) = AImat(0,2);

        AImat(1,2) = sum(etilde_square%H_square%S1%H%S2);
        AImat(2,1) = AImat(1,2);
          
          
        vec Dtheta =  solve(0.5*AImat, score);
        if (verbose){
            cout << "Dtheta is " << Dtheta << endl;
        }  
        theta = theta + Dtheta;
          
        // estimate the model parmater beta and mu  
        beta = inv_sympd(Wtilde.t()%repmat(H, 1, p).t()*Wtilde)*(Wtilde.t()%repmat(H, 1, p).t())*(Ytilde.col(i) - ztilde*mu);    
          
        if (mu_fixed == false){
           mu = 1.0/sum(ztilde%H%ztilde)*sum(ztilde.t()%H.t()*(Ytilde.col(i) - Wtilde*beta)); 
        }
  
        if (verbose){
            cout << "AImat is " << AImat << endl;
            cout <<  "update theta" << endl;
        }
          
        if (theta(0)<0.001){
            theta(0) = 0.001;
        }
        if (theta(1)<0.001){
            theta(1) = 0.001;
        }
        if (theta(2)<0.001){
            theta(2) = 0.001;
        }

        if (verbose){
            cout << "theta is " << theta << endl;
        }         
        
           
        // update the hyperparameter lengthscale

        cout << "calculate the loglikelihood for the first dataset and second dataset respectively, on a sequence of 10 candidate lengthscale" << endl;
        loglik1_all = zeros(num_ls, 1);
        loglik2_all = zeros(num_ls, 1);
        for (int k1 = 0; k1 < kscale1_vec.n_elem; k1++){ 
            Omega_11 = theta(1)*D1_all.col(k1) + theta(0)*I1;  // N-by-1 vector   
            Omega_22 = theta(2)*D2_all.col(k1) + theta(0)*I2;  // N-by-1 vector
            H_11 = 1.0/Omega_11;  
            H_22 = 1.0/Omega_22;
            etilde = Ytilde_all.slice(k1).col(i) - Wtilde_all.slice(k1)*beta - ztilde_all.col(k1)*mu;
            etilde1 = etilde.rows(0, n1-1);
            etilde2 = etilde.rows(n1, n1+n2-1);
            etilde1_square = etilde1%etilde1;
            etilde2_square = etilde2%etilde2; 
            loglik1_all(k1) = - 0.5*sum(log(1.0/H_11)) - 0.5*sum(etilde1_square%H_11); 
            loglik2_all(k1) = - 0.5*sum(log(1.0/H_22)) - 0.5*sum(etilde2_square%H_22);
        }
        if (verbose){
          cout << "search for the optimal lengthscale by maximizing the marginal loglikelihood on a sequence of user-defined lengthscales" << endl;     }
          idx1 = loglik1_all.index_max();      
          D1 = D1_all.col(idx1);
          if (verbose){
            cout << "search for the optimal lengthscale by maximizing the marginal loglikelihood on a sequence of user-defined lengthscale" << endl;  
          }
          idx2 = loglik2_all.index_max();      
          D2 = D2_all.col(idx2);

               
          if (verbose){
          cout << "lengthscale for the first dataset is " << kscale1 << endl;
          cout << "The logliklihood for the first dataset is " << loglik1_all << endl;
          }

          
        if (verbose){
            cout << "update Ytilde, Wtilde and ztilde"<< endl;
        }
          Ytilde = join_cols(Ytilde_all.slice(idx1).rows(0, n1-1), Ytilde_all.slice(idx2).rows(n1, n1+n2-1));
          Wtilde = join_cols(Wtilde_all.slice(idx1).rows(0, n1-1), Wtilde_all.slice(idx2).rows(n1, n1+n2-1));
          ztilde = join_cols(ztilde_all.col(idx1).rows(0, n1-1), ztilde_all.col(idx2).rows(n1, n1+n2-1));

          if (verbose){
          cout << "lengthscale for the second dataset is " << kscale2 << endl;
          cout << "The logliklihood for the second dataset is " << loglik2_all << endl;
          }
          

        diff = norm(theta - theta_old);
        theta_old = theta;
        if (verbose){
            cout <<  "diff is " << diff << endl;
        }
        iter++;
      }
      
        double loglik = loglik1_all(idx1) + loglik2_all(idx2);

        output(count, 0) = loglik;
        output(count, 1) = theta(0);
        output(count, 2) = theta(1);
        output(count, 3) = theta(2);
        output(count, 4) = mu;
        output(count, 5) = kscale1_vec(idx1);
        output(count, 6) = kscale2_vec(idx2); 
        output(count, 7) = idx1;
        output(count, 8) = idx2;  
        beta_all.col(count) = beta;        

        count++;
        if (verbose){
            cout <<  "count is " <<  count << endl;
        }

  }
    
    
  ret["output"] = output;
  ret["beta"] = beta_all;
  ret["kscale1"] = kscale1_vec;
  ret["kscale2"] = kscale2_vec;

  // construct return object
  
  
  return ret;
  
}





