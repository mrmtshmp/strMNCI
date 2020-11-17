#' gives MLEs for p1,p2 under restriction p1-p2=delta
#'
#' @param delta risk difference(Restriction)
#' @param y the number of cases (Restriction)
#' @param n the number of subjects (Restriction)
#'
#' @references Eqs.27 and 28 in "Comparative analysis of two rates", MIETTINEN AND NURMINEN, Stat. Med.(1985)
#'
#' @export
#'
resMLE.MN <- function(delta, y, n) {
  N=sum(n)
  C=sum(y)
  L3=N #e use   ??

  L2=(
    n[1] + 2*n[2]
    ) *
    delta - N - C

  L1=(
    n[2]*delta - N - 2*y[2]
    ) *
    delta + C

  L0 =
    y[2] * delta * (1-delta)

  # p2 <- "a", "b", "c" ;
  # a unique closed-form solution-following the guidelincs set forth by Bronshtein and Samendyayev
  # Bronshtein, J. N. and §emendyayev, K. A. A Guide Book to Mathematics, Verlag Harri Deutch, Frankfurt/Ivlain, 1973, pp. 161-163.

  c = # q in eq.28
    L2^3/(3*L3)^3 -
    L1*L2/(6*L3^2) +
    L0/(2*L3)

  b = # p in eq.28
    ifelse(c>=0,1,-1) * # "with the sign of p chosen so as to have it coincide with that of q."
    sqrt(
      L2^2/(3*L3)^2 - L1/(3*L3)
      )

  a =
    (pi+acos(c/b^3))/3

  p2 =
    2 * b*cos(a) - L2/(3*L3)

  #Miettinen and Nurminen
  return(c(p2+delta,p2))    #p1=p2+delta
  }
# End

 # Mee
  #z=(y[1]/n[1]-y[2]/n[2]-delta)/sqrt(p1*(1-p1)/n[1] + p2*(1-p2)/n[2])
  #z2=(y[1]/n[1]-y[2]/n[2]-delta)^2 / ((p1*(1-p1)/n[1] + p2*(1-p2)/n[2])*N/(N-1))
  #return(z2-qchisq(conflev,1))
