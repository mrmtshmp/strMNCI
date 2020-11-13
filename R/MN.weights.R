#' computes the MN wights necessary for CI for stratified data
#'
#' @references  Cochran-Mantel-Haenszel Weighted Miettinen & Nurminen Method for Confidence Intervals of the Difference in Binomial Proportions from Stratified 2x2 Samples, Kaifeng Lu (2008)
#' @param  delta
#' @param  dataset dataset needs to have centers as rows and number of successes in group1, number of trials in group 1, number of successes in group2, number of trials in group 2 as columns.
#' @param  bias.correct Bessel's correction for a finite sample.
#'
#' @export
#'
MN.weights <-
  function(delta, dataset, bias.correct=FALSE) {
    w =  # weight_CMH = n1*n2/(n1+n2) for each strata
      apply(
        dataset[,c(2,4)],1,
        prod) /
      rowSums(dataset[,c(2,4)])

    repeat {
          w.old = w
          w = w/sum(w)
          p=t( #' 2. Compute the weighted averages of constrained MLEs, Ri , using (6).
            apply(
              dataset,1,
              function(ro)
                resMLE.MN( # Restricted MLE (Miettinen, Nurminen 1985)
                  delta,
                  y=ro[c(1,3)],
                  n=ro[c(2,4)]
                  )
              )
            )
    #' contains the restricted MLE's for p1 and p2 (as a vector of length 2 ([1]:p1 [2]:p2))
    #' for each strata for given delta

    R = colSums(w*p)  # the weighted averages of constrained MLEs, Ri (i = 1,2)
    RR = R*(1-R)

    #' Update W_j by *Eq.16* in Miettinen, Nurminen 1985
    #' with or without bias correction for a finite sample.
    if (bias.correct)
      w = # updated weights with bias correction (for a finite sample)
      apply(
        dataset[,c(2,4)], 1,
        function(ro)   # ro: a vector with length 2 of the number of subjects in each groups (by strata-wise)
          (1 / (
              RR[1]/RR[2]/ro[1] + # RR = R*(1-R)
                1/ro[2]
              )) *
          (sum(ro)-1)/ sum(ro) # (not in Eq.16) Bias correction for a finite sample.
        )

    else w = apply( # updated weights without bias correction
      dataset[,c(2,4)],1,
      function(ro)
        (1 / (
          RR[1]/RR[2]/ro[1] +
            1/ro[2]
          ))
      )

    # Check the stop condition. (epsil <10^(-4))
    sum(abs(w-w.old))
    if (sum(abs(w-w.old))<10^(-4))
      break
    }
  return (w)
  }
# End
