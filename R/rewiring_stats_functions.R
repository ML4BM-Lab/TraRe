#' Generate stats for the rewiring method.
#'
#' @description
#'
#' When performing the rewiring test, some stats must be generated in order to proceed with the
#' rewiring method. `dave_test()` performs a permutation test from a data matrix and a group membership.
#' This function is used in the `runrewiring()` method.
#'
#' @param x data matrix containing subjects as rows and genes as columns.
#' @param grp array indicating the subject group membership.
#' @param perm number of permutations for the test.
#'
#' @return list containing the pvalue associated to the rewiring test, the initial t-stat
#' and the t-stat after the permutation test.
#'
#' @examples
#'
#'   ## We are going to generate a random matrix `x` and `grp` for the example.
#'   ## We will use 40 samples and 100 genes. First 25 samples will belong to one group
#'   ## and the rest (15) to the other.
#'
#'   mat <- matrix(stats::rnorm(40*100),40,100)
#'   group_m <- c(rep(1,25),rep(0,15))
#'
#'   ## Note: the `rewiring_test()` works with group membership (1,2) instead of (0,1)
#'
#'   results <- rewiring_test(x=mat,grp=1+group_m)
#'
#' @export
#'
#'
#'
rewiring_test <- function(x, grp, perm = 500) {

    p <- ncol(x)
    numgrp <- length(unique(grp))

    ## keep only rows with variance within groups
    bools <- rep(TRUE, p)
    names(bools) = colnames(x)
    for (g in unique(grp)) {
        bools <- bools & (apply(x[grp ==g, ], 2, stats::sd) != 0)
    }
    if (sum(bools) < p) {
        methods::show(paste(collapse = " ", c("...Rewiring Test - Dropping 0 Variance Genes:", paste(collapse = ",", names(which(bools ==
            FALSE))))))
        # show(c(p,sum(bools),which(bools==FALSE)))
        x <- x[, bools]
        p <- ncol(x)
    }


    ## test stat
    TS <- sum(apply(utils::combn(numgrp, 2), 2, function(pair) {
        as.vector((stats::cor(x[grp == pair[1], seq_len(p)]) - stats::cor(x[grp == pair[2], seq_len(p)]))^2)
    }))

    ## p
    T_star <- rep(NA, perm)
    for (j in seq_len(perm)) {
        grp_perm <- sample(grp)
        T_star[j] <- sum(apply(utils::combn(numgrp, 2), 2, function(pair) {
            as.vector((stats::cor(x[grp_perm == pair[1], ]) - stats::cor(x[grp_perm == pair[2], ]))^2)
        }))
    }
    # show(c(T, T_star))
    return(mean(c(TS, T_star) >= TS, na.rm = TRUE))

}
#' @export
#' @rdname rewiring_test
rewiring_test_pair_detail <- function(x, grp, perm = 500) {
    p <- ncol(x)
    numgrp <- length(unique(grp))

    ## keep only rows with variance within groups
    bools <- rep(TRUE, p)
    names(bools) = colnames(x)
    for (g in unique(grp)) {
        bools <- bools & (apply(x[grp ==g, ], 2, stats::sd) != 0)
    }
    if (sum(bools) < p) {
        methods::show(paste(collapse = " ", c("...Detailed Rewiring Test - Dropping 0 Variance Genes:", paste(collapse = ",", names(which(bools ==
            FALSE))))))
        # show(c(p,sum(bools),which(bools==FALSE)))
        x <- x[, bools]
        p <- ncol(x)
    }

    ## test stat
    TS <- sum(apply(utils::combn(numgrp, 2), 2, function(pair) {
        as.vector((stats::cor(x[grp == pair[1], seq_len(p)]) - stats::cor(x[grp == pair[2], seq_len(p)]))^2)
    }))

    ## perm
    T_star <- rep(NA, perm)
    for (j in seq_len(perm)) {
        grp_perm <- sample(grp)
        T_star[j] <- sum(apply(utils::combn(numgrp, 2), 2, function(pair) {
            as.vector((stats::cor(x[grp_perm == pair[1], ]) - stats::cor(x[grp_perm == pair[2], ]))^2)
        }))
    }

    pval <- mean(c(TS, T_star) >= TS, na.rm = TRUE)
    return(list(pval = pval, T_star = T_star, TS = TS))
}

