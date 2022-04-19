library(MASS);
#' get_sq_dist.
#' 
#' @description
#' get_sq_dist used to get Euclidean distance of all spots. 
#'
#' @details get_sq_dist used to get Euclidean distance of all spots. 
#' @param spa_mat is a matrix of Euclidean distance of all spots
#' @return a matrix of Euclidean distance
#' @export
get_sq_dist <- function(spa_mat) {
  n <- dim(spa_mat)[1]
  sq_dist <- matrix(0, n, n)
  diag(sq_dist) <- 0
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      sq_dist[i,j] <- sum( (spa_mat[i,]-spa_mat[j,])^2 )
      sq_dist[j,i] <- sq_dist[i,j]
    }
  }
  return(sq_dist)
}


#' get_K_gauss.
#' 
#' @description
#' get_K_gauss used to get the Gaussian kernel. 
#'
#' @details get_K_cos used to get the Gaussian kernel. 
#' @param sq_dist is a matrix of Euclidean distance of all spots
#' @param length_scale is a length scale parameter of the Gaussian kernel
#' @return a Gaussian kernel
#' @export
get_K_gauss <- function(sq_dist, length_scale) {
  K_gauss <- exp(-sq_dist / length_scale^2 / 2)
  return(K_gauss)
} 




#' get_K_cos.
#' 
#' @description
#' get_K_cos used to get the Cosine kernel. 
#'
#' @details get_K_cos used to get the Cosine kernel. 
#' @param sq_dist is a matrix of Euclidean distance of all spots
#' @param phi is a periodicity parameter of the Cosine kernel
#' @return a cosine kernel
#' @export
get_K_cos <-function(sq_dist, phi) {
  K_cos <- cos(2*pi*sqrt(sq_dist)/phi)
  return(K_cos)
}


#' select_gene.
#' 
#' @description
#' select_gene used to select the genes that have the same type of kernels under both null and alternative hypothesis. 
#'
#' @details select_gene used to select the genes that have the same type of kernels under both null and alternative hypothesis. 
#' @param fit a list
#' @return a list 
#' @export
select_gene <-function(fit) {
  n1 = length(fit$kernel1_all_vec)
  n2 = length(fit$kernel2_all_vec)
  
  idx = which(abs(fit$kernel1_length_scale_periodicity_idx_alternative - fit$kernel1_length_scale_periodicity_idx_null)<n1 &
  abs(fit$kernel2_length_scale_periodicity_idx_alternative - fit$kernel2_length_scale_periodicity_idx_null)<n2)

  result = list(chi_stat = fit$chi_stat[idx], 
             pvalue = fit$pvalue[idx],
             sigma_e_alternative = fit$sigma_e_alternative[idx],
             sigma_kernel1_alternative=fit$sigma_kernel1_alternative[idx],
             sigma_kernel2_alternative=fit$sigma_kernel2_alternative[idx],
             mu=fit$mu[idx],
             kernel1_length_scale_periodicity_alternative = fit$kernel1_length_scale_periodicity_alternative[idx],
             kernel2_length_scale_periodicity_alternative = fit$kernel2_length_scale_periodicity_alternative[idx],
             kernel1_length_scale_periodicity_idx_alternative = fit$kernel1_length_scale_periodicity_idx_alternative[idx],
             kernel2_length_scale_periodicity_idx_alternative = fit$kernel2_length_scale_periodicity_idx_alternative[idx],
             kernel1_all_vec = fit$kernel1_all_vec,
             kernel2_all_vec = fit$kernel1_all_vec,
             sigma_e_null=fit$sigma_e_null[idx],
             sigma_kernel1_null=fit$sigma_kernel1_null[idx],
             sigma_kernel2_null=fit$sigma_kernel2_null[idx],
             kernel1_length_scale_periodicity_null = fit$kernel1_length_scale_periodicity_null[idx],
             kernel2_length_scale_periodicity_null = fit$kernel2_length_scale_periodicity_null[idx],
             kernel1_length_scale_periodicity_idx_null = fit$kernel1_length_scale_periodicity_idx_null[idx],
             kernel2_length_scale_periodicity_idx_null = fit$kernel2_length_scale_periodicity_idx_null[idx])
  
  result
}



