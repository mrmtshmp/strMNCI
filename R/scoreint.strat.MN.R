#' computes the CI for a common risk difference from several stratified 2x2 tables
#'
#' @importFrom  stats qchisq
#' @importFrom  stats uniroot
#' @param  dataset a matrix-class object with each row has to contain y1, n1, y2 and n2, in that order, where y1 is the binomial response out of n1 trials for the treatment group and y2 is the binomial response out of n2 trials for the control group. Each row presents a stratum.
#' @param  weight  "inv.Var", "Add4" or "Cochran". MN: "Miettinen Nurminen (1985)" The solution for "scorestat.MN" via uniroot().
#' @param  conflev
#' @param  bias.correct Bessel's correction for a finite sample.
#'
#' @references  Miettinen And Nurminen 1985, Stat. in Medicine
#'
#' @export

scoreint.strat.MN <-
  function(
    dataset, weight= "MN", conflev=0.95,
    bias.correct=FALSE) {

    if(!(weight %in% c("inv.Var", "Add4", "Cochran","MN"))) stop("wrong specification in the argument 'weight'")

    delta.max = 1
    delta.min = -1 # For uniroot() solution for scorestat.MN.
    eps=10^(-4)    # Stop condition for weight update by resMLE.

    #getting initial estimte of point estimate to start numerical search:
    if (weight == "inv.Var") {
      p1 = dataset[,1]/dataset[,2]
      p2 = dataset[,3]/dataset[,4]
      w = (p1*(1-p1)/dataset[,2] + p2*(1-p2)/dataset[,4])^(-1)
      } else
        if (weight=="Add4") {
          p1 = (dataset[,1]+1)/(dataset[,2]+2)
          p2 = (dataset[,3]+1)/(dataset[,4]+2)
          w = (p1*(1-p1)/dataset[,2] + p2*(1-p2)/dataset[,4])^(-1)
          } else
            w = apply(
              dataset[,c(2,4)],1,
              prod
              ) / rowSums(dataset[,c(2,4)])
            #w =n1*n2/(n1+n2)

    w = w/sum(w) # standardized weights
    d =          # sample risk differences
      dataset[,1]/dataset[,2] -
      dataset[,3]/dataset[,4]

    point.est = sum(w*d)

    if(point.est==delta.min) #to handle extreme cases
      delta.lb=delta.min

    else delta.lb =
      uniroot( # uniroot() function searches the interval from lower to upper for a root (i.e., zero) of the function f with respect to its first argument.. (http://www.endmemo.com/r/uniroot.php)
        f = scorestat.MN,
        interval=
          c(
            delta.min+eps,
            point.est-eps
            ),
        dataset=dataset,
        weight=weight,
        conflev=conflev,
        bias.correct=bias.correct
        )$root

    if(point.est==delta.max) #to handle y=0 or y=n
      delta.ub=delta.max
    else delta.ub =
      uniroot(
        scorestat.MN,
        interval =
          c(point.est+eps,delta.max-eps),
        dataset =
          dataset, weight=weight, conflev=conflev, bias.correct=bias.correct
        )$root
    return(
      c(delta.lb,delta.ub)
      )
  }
