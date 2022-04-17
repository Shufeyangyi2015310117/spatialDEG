#include "mt_paral_testing_job.hpp"


using namespace std;
using namespace arma;
using namespace Rcpp;


struct Obj_Single_gene_test{
	colvec output;
    colvec beta;
	uvec idx_K;
    
};


//' Do inverse of sysmetric matrix 
//' @param Min A sysmetric matrix
//' 
//' @return A list
//' 
//' @export
// [[Rcpp::export]]
arma::mat SysMatEigen3(arma::mat M) {
		arma::vec eigval = zeros<vec>( M.n_rows );
		arma::mat eigvec = zeros<mat>( size(M) );
		eig_sym(eigval, eigvec, M, "dc");
		const uvec idx = find(eigval < 1e-8 );
		arma::vec tmp_value = ones<arma::vec>(idx.n_elem);
		eigval.elem( idx ) = tmp_value * 1e-8;
		arma::mat M1 = eigvec.each_row() % eigval.t();
		M = M1 * eigvec.t();
		// return values
		return M;
		//return List::create(Named("eigval") = eigval, Named("eigvec") = eigvec);
}// end func


Obj_Single_gene_test Single_gene_test(mat& D1_all, mat& D2_all, arma::cube& Wtilde_all, arma::cube& Ytilde_all, arma::mat& ztilde_all, vec& kscale1_vec, vec& kscale2_vec, uvec& idx_K, arma::vec& Initial_theta, int num_ls = 40, int n1 = 10, int n2 = 10, int p = 1, int ngene = 10, int g = 0, int iteration = 100, bool mu_fixed = false, bool Kernel_fixd = false, bool verbose = false, string kernel1 = "all", string kernel2 = "all", double tol = 0.0001){
     
      colvec output = zeros(13, 1);
     
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
      mat beta_all = zeros(p, ngene);

      int count = 0;
      double diff;

      vec H, S1, S2, H_square;
      vec etilde, etilde_square;

      vec etilde1, etilde2, etilde1_square, etilde2_square;
      uword idx1, idx2;

      double loglik, loglik_old; 
      vec loglik1_all, loglik2_all;
      mat loglik_all = zeros(ngene, iteration);
    
      double kscale1, kscale2;
       
     
    // define the model parameters          
      colvec theta = zeros(3, 1);
      theta = Initial_theta; // sigma_e, sigma_1, sigma_2
      loglik_old = -datum::inf;
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
      while(diff > tol && iter < iteration){

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
        etilde = Ytilde.col(g) - Wtilde*beta - ztilde*mu; 
          
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
          
          
          
        mat AImat2 = 0.5*AImat;
          /*
        while (cond(AImat2) < 1e-10){
            AImat2 = AImat2 + eye(3,3)*1e-10;
        }*/
          
          if (verbose){
              cout << "rcond AImat2 is " << rcond(AImat2) << endl;
              cout << "cond AImat2 is " << rcond(AImat2) << endl;
          }
          
          
          
        if (verbose){
            cout << AImat2 << endl;
        }
          
          vec eigval;
          eig_sym(eigval, AImat2);
          if (min(eigval) < 0){
            AImat2 = SysMatEigen3(AImat2); 
          }
          
        if (verbose){
            cout << eigval << endl;
        }
          
        if (verbose){
            cout << AImat2 << endl;
        }
              
        vec Dtheta =  solve(AImat2, score);
        
          
          
          
        if (verbose){
            cout << "Dtheta is " << Dtheta << endl;
        }  
        vec theta0 = theta;
          
          double stepsize = 1;
          if (iter<100){
             theta = theta + stepsize*Dtheta; 
          }else if (iter > 99 & iter <200){
              stepsize = 0.1;
              theta = theta + stepsize*Dtheta;
          } else if (iter > 199 & iter <300){
               stepsize = 0.01;
              theta = theta + stepsize*Dtheta;
          }else if (iter > 299 & iter < 400){
              stepsize = 0.001;
              theta = theta + stepsize*Dtheta;
          }else if (iter > 399 & iter < 500){
                 stepsize = 0.0001;
              theta = theta + stepsize*Dtheta;
          }else if (iter > 499 & iter < 600){
              stepsize = 0.00001;
               theta = theta + stepsize*Dtheta;
          }else{
              stepsize = 0.000001;
              theta = theta + stepsize*Dtheta;
          }

        if (verbose){
          cout << "stepsize is " << stepsize << endl;
            }
        /*  
        double step = 1;
        while (any(theta < 0.0)) {
				step *= 0.5;
				theta = theta0 + step * Dtheta;
			} // end while
            */
          //theta.elem(find(theta < 0.000001)).zeros();
          
        // estimate the model parmater beta and mu  
        beta = inv_sympd(Wtilde.t()%repmat(H, 1, p).t()*Wtilde)*(Wtilde.t()%repmat(H, 1, p).t())*(Ytilde.col(g) - ztilde*mu);    
          
        if (mu_fixed == false){
           mu = 1.0/sum(ztilde%H%ztilde)*sum(ztilde.t()%H.t()*(Ytilde.col(g) - Wtilde*beta)); 
        }
  
        if (verbose){
            cout << "AImat is " << AImat << endl;
            cout <<  "update theta" << endl;
        }
          
          
        if (theta(0)<0.000001){
            theta(0) = 0.000001;
        }
        if (theta(1)<0.000001){
            theta(1) = 0.000001;
        }
        if (theta(2)<0.000001){
            theta(2) = 0.000001;
        }

        if (verbose){
            cout << "theta is " << theta << endl;
        }         
        
           
        // update the hyperparameter lengthscale
        if (verbose){
        cout << "calculate the loglikelihood for the first dataset and second dataset respectively, on a sequence of 10 candidate lengthscale" << endl;}
        loglik1_all = zeros(num_ls, 1);
        loglik2_all = zeros(num_ls, 1);       
        for (int k1 = 0; k1 < kscale1_vec.n_elem; k1++){ 
            Omega_11 = theta(1)*D1_all.col(k1) + theta(0)*I1;  // N-by-1 vector   
            Omega_22 = theta(2)*D2_all.col(k1) + theta(0)*I2;  // N-by-1 vector
            H_11 = 1.0/Omega_11;  
            H_22 = 1.0/Omega_22;
            etilde = Ytilde_all.slice(k1).col(g) - Wtilde_all.slice(k1)*beta - ztilde_all.col(k1)*mu;
            etilde1 = etilde.rows(0, n1-1);
            etilde2 = etilde.rows(n1, n1+n2-1);
            etilde1_square = etilde1%etilde1;
            etilde2_square = etilde2%etilde2; 
            loglik1_all(k1) = - 0.5*sum(log(1.0/H_11)) - 0.5*sum(etilde1_square%H_11); 
            loglik2_all(k1) = - 0.5*sum(log(1.0/H_22)) - 0.5*sum(etilde2_square%H_22);
        }
        if (verbose){
          cout << "search for the optimal lengthscale by maximizing the marginal loglikelihood on a sequence of user-defined lengthscales" << endl;     }
          
          if (Kernel_fixd == true){
              idx1 = idx_K(0);
          }else{
              if (kernel1 == "all"){
                idx1 = loglik1_all.index_max(); 
              }else if(kernel1 == "Gaussian"){
                idx1 = loglik1_all(span(0, num_ls/2-1)).index_max(); 
              }else if (kernel1 == "cosine"){
                idx1 = loglik1_all(span(num_ls/2, num_ls-1)).index_max(); 
                idx1 = idx1 + num_ls/2;
              }else{
                  cout << "please specify the correct kernel: all, Gaussian, cosine" << endl;
              }

          }
          
          /*
          if (theta(1)==0.000001){
              idx1 = 0;
          }*/
          
          kscale1 = kscale1_vec(idx1); 
          D1 = D1_all.col(idx1);
          if (verbose){
            cout << "search for the optimal lengthscale by maximizing the marginal loglikelihood on a sequence of user-defined lengthscale" << endl;  
          }
          if (Kernel_fixd == true){
              idx2 = idx_K(1);
          }else{
              if (kernel1 == "all"){
                idx2 = loglik2_all.index_max(); 
              }else if(kernel2 == "Gaussian"){
                idx2 = loglik2_all(span(0, num_ls/2-1)).index_max(); 
              }else if (kernel2 == "cosine"){
                idx2 = loglik2_all(span(num_ls/2, num_ls-1)).index_max(); 
                idx2 = idx2 + num_ls/2;
              }else{
                cout << "please specify the correct kernel: all, Gaussian, cosine" << endl;
              }
          }
          /*
          if (theta(2)==0.000001){
              idx2 = 0;
          }*/

          kscale2 = kscale2_vec(idx2);
          D2 = D2_all.col(idx2);

          
          loglik_all(g, iter) = loglik1_all(idx1) + loglik2_all(idx2);
          
          
            if (verbose){
                cout << "loglik_all is " << loglik_all(g, iter) << endl;
            }
          
          
          if (verbose){
          cout << "lengthscale for the first dataset is " << kscale1 << endl;      
          cout << "The logliklihood for the first dataset is " << fixed << setprecision(9) << loglik1_all << endl;
          loglik1_all.raw_print(std::cout);
          loglik1_all.save("loglik1_all.txt", raw_ascii);
          }

          
        if (verbose){
            cout << "update Ytilde, Wtilde and ztilde"<< endl;
        }
          Ytilde = join_cols(Ytilde_all.slice(idx1).rows(0, n1-1), Ytilde_all.slice(idx2).rows(n1, n1+n2-1));
          Wtilde = join_cols(Wtilde_all.slice(idx1).rows(0, n1-1), Wtilde_all.slice(idx2).rows(n1, n1+n2-1));
          ztilde = join_cols(ztilde_all.col(idx1).rows(0, n1-1), ztilde_all.col(idx2).rows(n1, n1+n2-1));

          if (verbose){
          cout << "lengthscale for the second dataset is " << kscale2 << endl;
          cout << "The logliklihood for the second dataset is " << fixed << setprecision(9) << loglik2_all << endl;
          loglik2_all.raw_print(std::cout);
          }
          if (verbose){
            cout << "diff is " << loglik_all(g, iter) - loglik_old << endl;
          }

 
        diff = abs(loglik_all(g, iter) - loglik_old);
        loglik_old = loglik_all(g, iter);
          
        if (verbose){
        if (iter%1 == 0){
            cout << "iteration is " << iter << endl;
            cout << "theta is " << theta.t() << endl;
            cout <<  "absdiff is " << diff << endl; 
        }
        }
          
         
        iter++;
      }
     
        loglik = loglik1_all(idx1) + loglik2_all(idx2);
    
        output(0) = loglik;
        output(1) = theta(0);
        output(2) = theta(1);
        output(3) = theta(2);
        output(4) = mu;
        output(5) = kscale1_vec(idx1);
        output(6) = kscale2_vec(idx2);   
        output(7) = iter;
        output(8) = idx1+1;
        output(9) = idx2+1;
    
        Obj_Single_gene_test out;
        out.output = output;
        out.beta = beta;
        out.idx_K = {idx1, idx2};
     
     return out;
}






void parGene_SpDEG::loop_by_gene_SpDEG(int g){


	// set the parameter
    uvec idx_K = {0,0};
	Obj_Single_gene_test out_alter = Single_gene_test(D1_all, D2_all, Wtilde_all, Ytilde_all, ztilde_all, kscale1_vec, kscale2_vec, idx_K, Initial_theta, num_ls, n1, n2, p, ngene, g, iteration, false, false, verbose, "all", "all", tol);
    
    
    idx_K = {out_alter.idx_K(0), out_alter.idx_K(1)};
    string kernel1, kernel2;
    if (kernel_matched == true){
        
        if (out_alter.output(8) <= num_ls/2){
            kernel1 = "Gaussian";
        }else{
            kernel1 = "cosine";
        }
        
        if (out_alter.output(9) <= num_ls/2){
            kernel2 = "Gaussian";
        }else{
            kernel2 = "cosine";
        }
        
    }else{
        kernel1 = "all";
        kernel2 = "all";
    }

	Obj_Single_gene_test out_null = Single_gene_test(D1_all, D2_all, Wtilde_all, Ytilde_all, ztilde_all, kscale1_vec, kscale2_vec, idx_K, Initial_theta, num_ls, n1, n2, p, ngene, g, iteration, true, Kernel_fixed, verbose, kernel1, kernel2, tol);


    for (int i = 0; i < 10; i++){
        out_param(g, i) = out_alter.output(i);
        out_param(g, i+10) = out_null.output(i);
    }

	// reset


	if ((g + 1) % 100 == 0 && (g + 1) != 0){
		cout << g + 1 << "-th Gene starts working ..." << endl;
	}

}

std::mutex _mtx22;
int parGene_SpDEG::next_SpDEG(){
	std::lock_guard<std::mutex> lockGuard(_mtx22);
	if (current_idx >= ngene){
		return -1;
	}
	current_idx++;
	return current_idx - 1;
}

void parGene_SpDEG::update_by_thread_SpDEG(int thread_id){
	while (true){
		int idx = next_SpDEG();
		if (idx == -1){
			break;
		}
		loop_by_gene_SpDEG(idx);
	}
}



