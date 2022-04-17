#ifndef mt_paral_testing_job_hpp
#define mt_paral_testing_job_hpp

#include <stdio.h>
#include <math.h>
#include <RcppArmadillo.h>
#include <thread>
#include <mutex>


using namespace std;
using namespace arma;
using namespace Rcpp;


class parGene_SpDEG{
public:
	int current_idx = 0;
	int Ngene;

    cube Ytilde_all, Wtilde_all;
	mat D1_all, D2_all, ztilde_all;
    vec Initial_theta, kscale1_vec, kscale2_vec;
    int num_ls, n1, n2, p, ngene, iteration;
    bool Kernel_fixed, verbose, kernel_matched;
	mat out_param;
    double tol;

	parGene_SpDEG(mat& D1_all, mat& D2_all, cube& Ytilde_all, cube& Wtilde_all, mat& ztilde_all, 
                  vec& Initial_theta, vec& kscale1_vec, vec& kscale2_vec, 
                  int num_ls, int n1, int n2, int p, int ngene, int iteration, 
                  bool Kernel_fixed, bool verbose, bool kernel_matched, mat& out_param, double tol){

        this->D1_all = D1_all;
        this->D2_all = D2_all;
		this->Ytilde_all = Ytilde_all;
		this->Wtilde_all = Wtilde_all;
        this->ztilde_all = ztilde_all;
        this->Initial_theta = Initial_theta;
        this->kscale1_vec = kscale1_vec;
        this->kscale2_vec = kscale2_vec;
        this->num_ls = num_ls;
        this->n1 = n1;
        this->n2 = n2;
        this->p = p;
        this->ngene = ngene;
        this->iteration = iteration;
        this->Kernel_fixed = Kernel_fixed;
        this->verbose = verbose;
        this->kernel_matched = kernel_matched;
        this->out_param = out_param;
        this->tol =tol;
        
	}

	void loop_by_gene_SpDEG(int g);
	void update_by_thread_SpDEG(int thread_id);
	int  next_SpDEG();

};


#endif 

