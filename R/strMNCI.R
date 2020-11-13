#' R program code for the computation of the Miettinen-Nurminen confidence limits for stratified 2x2 data 
#' by Bernhard Klingenberg, Department of Mathematics & Statistics, Williams College,Williamtown, MA, USA, 31.8.2012.
#' For the code to work, each row has to contain y1, n1, y2 and n2, in that order,
#' where y1 is the binomial response out of n1 trials for the treatment group and y2 is
#' the binomial response out of n2 trials for the control group. Each row presents a
#' stratum. The procedure then computes the CI for the common difference in the
#' proportions of the treated group minus the control group.
#' 
#' Below are the necessary functions, the very last one is the main function that you
#' call to get the CI. E.g., if we had three strata

dataset = 
  matrix(
    c(23, 50, 32, 53, 25, 51, 35, 49, 25, 53, 31, 48),
    ncol=4,byrow=TRUE
    )

scoreint.strat.MN(dataset)

# shows the 95% CI for the common difference to be [-0.29, -0.07].
# One can choose different weights through e.g.,

scoreint.strat.MN(dataset, weight="Cochran")

# One can also use the bias correction to be applied to the weights, i.e., the weights
# are now unbiased estimators under the null of no differences of the variance of the
# difference of proportion.

scoreint.strat.MN(dataset, bias.correct=TRUE)

# End