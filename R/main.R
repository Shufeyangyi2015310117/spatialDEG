#' spatialDEG.
#' 
#' @description
#' spatialDEG used to perform differentially expressed gene analysis by leveragin spatial informaiton. 
#'
#' @details spatialDEG used to perform differentially expressed gene analysis by leveragin spatial informaiton. 
#' @param Y1 is a vector of gene expression for the first dataset.
#' @param Y2 is a vector of gene expression for the second dataset.
#' @param spa1 is a matrix of spatial location for the first dataset.
#' @param spa2 is a matrix of spatial location for the second dataset.
#' @param W1 is a matrix of covariates for the first dataset.
#' @param W2 is a matrix of covariates for the second dataset.
#' @param Initial_theta is a vector of initial value for theta. The default vector is ones
#' @param num_ls is a integer specifying the number of length scale or periodicity parameters.
#' @param Kernel1 is a matrix of Kernel for the first dataset.
#' @param Kernel2 is a matrix of Kernel for the second dataset.
#' @param parallel is a logical value specifying whether to run spatialDEG in parallel.
#' @param coreNum is an integer specifying the number of cores used for parallel computation.
#' @param max_Iter is an integer specifying the maximum iteration of the AI algorithm.
#' @param eps a value specifying when the AI average convergence.
#' @param kernel_fixed Kernel_fixd a logical value specifying whether the kernel is fixed
#' @param kernel_matched a logical value specifying whether the the kernels under both null hypothesis and alternative hypothesis are restricted to be the same type of kernel
#' @param check_positive a logical value specifying whether to check the positive definiteness of kernels
#' @param verbose is logical value specifying whether to print the parameters during spatialDEG analysis
#' @return a list
#' @export
spatialDEG <- function(Y1, Y2, spa1=NULL, spa2=NULL, W1, W2, Initial_theta = matrix(1, 3, 1), num_ls = 10, Kernel1=NULL, Kernel2=NULL, parallel = TRUE, kernel_fixed = FALSE, kernel_matched = FALSE, coreNum = 1, max_Iter = 1000, check_positive = TRUE, eps = 0.00001, verbose = FALSE) {
  Y = rbind(Y1, Y2)
  n1 = dim(W1)[1]
  n2 = dim(W2)[1]
  q1 = dim(W1)[2]
  q2 = dim(W2)[2]
  W11 = W1
  W12 = matrix(0, n1, q2)
  W21 = matrix(0, n2, q1)
  W22 = W2
  W = cbind(rbind(W11, W21), rbind(W12, W22))
  tol = eps
  iteration = max_Iter
  verbose = verbose
  Initial_theta = Initial_theta
  Kernel_fixd = kernel_fixed
  kernel_matched = kernel_matched
  K1 = Kernel1
  K2 = Kernel2
  check_positive = check_positive
  if (!is.null(spa1) & !is.null(spa2) & is.null(Kernel1) & is.null(Kernel2)){
    if (parallel == TRUE){
      out = spatialDEG_paral_test(spa1, spa2, W, Y, Initial_theta, num_ls, iteration, Kernel_fixd, verbose, kernel_matched, coreNum, tol) 
      chi_stat = 2*(out$out_param[,1] - out$out_param[,11])
      pvalue = 1 - pchisq(chi_stat, 1)
      fit = list(chi_stat = chi_stat, 
                 pvalue = pvalue,
                 sigma_e_alternative = out$out_param[,2],
                 sigma_kernel1_alternative=out$out_param[,3],
                 sigma_kernel2_alternative=out$out_param[,4],
                 mu=out$out_param[,5],
                 kernel1_length_scale_periodicity_alternative = out$out_param[,6],
                 kernel2_length_scale_periodicity_alternative = out$out_param[,7],
                 kernel1_length_scale_periodicity_idx_alternative = out$out_param[,9],
                 kernel2_length_scale_periodicity_idx_alternative = out$out_param[,10],
                 kernel1_all_vec = out$kscale1,
                 kernel2_all_vec = out$kscale2,
                 sigma_e_null=out$out_param[,12],
                 sigma_kernel1_null=out$out_param[,13],
                 sigma_kernel2_null=out$out_param[,14],
                 kernel1_length_scale_periodicity_null = out$out_param[,16],
                 kernel2_length_scale_periodicity_null = out$out_param[,17],
                 kernel1_length_scale_periodicity_idx_null = out$out_param[,19],
                 kernel2_length_scale_periodicity_idx_null = out$out_param[,20])
    }else{
      out = spatialDEG_test(spa1, spa2, W, Y, Initial_theta, num_ls, mu_fixed = FALSE, verbose)
      out2 = spatialDEG_test(spa1, spa2, W, Y, Initial_theta, num_ls, mu_fixed = TRUE, verbose) 
      chi_stat = 2*(out$output[,1] - out2$output[,1])
      pvalue = 1 - pchisq(chi_stat, 1)
      fit = list(chi_stat = chi_stat, 
                 pvalue = pvalue,
                 sigma_e_alternative = out$output[,2],
                 sigma_kernel1_alternative=out$output[,3],
                 sigma_kernel2_alternative=out$output[,4],
                 mu=out$output[,5],
                 kernel1_length_scale_periodicity_alternative = out$output[,6],
                 kernel2_length_scale_periodicity_alternative = out$output[,7],
                 kernel1_length_scale_periodicity_idx_alternative = out$output[,8],
                 kernel2_length_scale_periodicity_idx_alternative = out$output[,9],
                 beta_alternative = out$beta,
                 sigma_e_null=out2$output[,2],
                 sigma_kernel1_null=out2$output[,3],
                 sigma_kernel2_null=out2$output[,4],
                 beta_null = out2$beta,
                 kernel1_length_scale_periodicity_null = out$output[,6],
                 kernel2_length_scale_periodicity_null = out$output[,7],
                 kernel1_length_scale_periodicity_idx_null = out$output[,8],
                 kernel2_length_scale_periodicity_idx_null = out$output[,9])
    }
  }else{
    if (is.null(spa1) & is.null(spa2) & !is.null(Kernel1) & !is.null(Kernel2)){
      out = spatialDEG_true_kernel_test(W, Y, K1, K2, Initial_theta, mu_fixed = FALSE, verbose, check_positive)
      out2 = spatialDEG_true_kernel_test(W, Y, K1, K2, Initial_theta, mu_fixed = TRUE, verbose, check_positive)
      chi_stat = 2*(out$output[,1] - out2$output[,1])
      pvalue = 1 - pchisq(chi_stat, 1)
      fit = list(chi_stat = chi_stat, 
                 pvalue = pvalue,
                 sigma_e_alternative = out$output[,2],
                 sigma_kernel1_alternative=out$output[,3],
                 sigma_kernel2_alternative=out$output[,4],
                 beta_alternative = out$beta,
                 mu=out$output[,5],
                 sigma_e_null=out2$output[,2],
                 sigma_kernel1_null=out2$output[,3],
                 sigma_kernel2_null=out2$output[,4],
                 beta_null = out2$beta)
    }else{
      print("error")
    }

  }
  fit
}


