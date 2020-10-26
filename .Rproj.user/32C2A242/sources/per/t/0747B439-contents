#' Generate stats from the rewiring test.
#'
#' @description
#'
#' When performing the rewiring test, some stats must be generated in order to proceed with the
#' rewiring method. `dave_test()` generates the correlations to work with. `dave_test_pair_detail()`
#' can produce the correlation plots.
#'
#' @param x data matrix containing subjects as rows and genes as columns.
#' @param grp indicator of group membership
#' @param perm number of permutations for the test.
#'
#' @examples
#' \dontrun{
#'   example here
#' }
#'
#' @export dave_test
dave_test = function(x, grp, perm = 500) {
  p = ncol(x)
  numgrp = length(unique(grp))

  ## keep only rows with variance within groups
  bools = rep(TRUE, p)
  for (g in unique(grp)){
    bools = bools & (apply(x[grp ==g, ], 2, stats::sd) != 0)
  }
  if(sum(bools)<p){
    methods::show(paste(collapse = " ", c("...Rewiring Test - Dropping 0 Variance Genes:", paste(collapse=",",names(which(bools==FALSE))))))
#    show(c(p,sum(bools),which(bools==FALSE)))
    x = x[,bools]
    p = ncol(x)
  }


  ## test stat
  T = sum(apply(utils::combn(numgrp, 2), 2, function(pair) {
    as.vector((stats::cor(x[grp == pair[1], 1:p]) -
                 stats::cor(x[grp == pair[2], 1:p]))^2)
  }))

  ## p
  T_star = rep(NA, perm)
  for (j in 1:perm) {
    grp_perm = sample(grp)
    T_star[j] = sum(apply(utils::combn(numgrp, 2), 2, function(pair) {
      as.vector((stats::cor(x[grp_perm == pair[1],]) -
                   stats::cor(x[grp_perm == pair[2],]))^2)
    }))
  }
  #show(c(T, T_star))
  return(mean(c(T, T_star) >= T, na.rm=TRUE))
}
#' @export
#' @rdname dave_test
dave_test_pair_detail = function(x, grp, perm = 500) {
  p = ncol(x)
  numgrp = length(unique(grp))

  ## keep only rows with variance within groups
  bools = rep(TRUE, p)
  for (g in unique(grp)){
    bools = bools & (apply(x[grp ==g, ], 2, stats::sd) != 0)
  }
  if(sum(bools)<p){
    methods::show(paste(collapse = " ", c("...Detailed Rewiring Test - Dropping 0 Variance Genes:", paste(collapse=",",names(which(bools==FALSE))))))
    #show(c(p,sum(bools),which(bools==FALSE)))
    x = x[,bools]
    p = ncol(x)
  }

  ## test stat
  T = sum(apply(utils::combn(numgrp, 2), 2, function(pair) {
    as.vector((stats::cor(x[grp == pair[1], 1:p]) -
                 stats::cor(x[grp == pair[2], 1:p]))^2)
  }))

  ## perm
  T_star = rep(NA, perm)
  for (j in 1:perm) {
    grp_perm = sample(grp)
    T_star[j] = sum(apply(utils::combn(numgrp, 2), 2, function(pair) {
      as.vector((stats::cor(x[grp_perm == pair[1],]) -
                   stats::cor(x[grp_perm == pair[2],]))^2)
    }))
  }

  pval = mean(c(T, T_star) >= T, na.rm=TRUE)
  return(list(pval=pval, T_star=T_star, T=T))
}
