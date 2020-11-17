#' computes restricted MLEs for pi1, pi2 under restriction pi1-pi2=delta for a given 2x2 table and returns components of MN chi-square function
#'
#' @param delta risk difference(Restriction)
#' @param y the number of cases (Restriction)
#' @param n the number of subjects (Restriction)
#'
#' @export
resMLE.MN1 <- function(delta, y, n) {
  N=sum(n)
  C=sum(y)
  L3=N
  L2=(n[1]+2*n[2])*delta-N-C
  L1=(n[2]*delta-N-2*y[2])*delta+C
  L0=y[2]*delta*(1-delta)
  c=L2^3/(3*L3)^3 - L1*L2/(6*L3^2) + L0/(2*L3)
  b=ifelse(c>=0,1,-1)*sqrt(L2^2/(3*L3)^2 - L1/(3*L3))
  a=(pi+acos(c/b^3))/3
  p2=2*b*cos(a)-L2/(3*L3)
  p1=p2+delta
  numer = y[1]/n[1] - y[2]/n[2] - delta
  denom = (p1*(1-p1)/n[1] + p2*(1-p2)/n[2])*sum(n)/sum(n-1)
  return(c(numer,denom))
}
