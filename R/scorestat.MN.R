#' computes the chi-square function
#' @references  Miettinen and Nurminen 1985, Stat. in Medicine, page 218
#' @param delta0   a given null value
#' @param dataset  stratified 2x2 tables, needs to have centers as rows. number of successes in group1, number of trials in group 1,  number of successes in group2, number of trials in group 2
# as columns.
#' @param weight Weight specification.
#' @param conflev Confidence level
#' @param bias.correct Bessel's correction for a finite sample.
#'
#' @export

scorestat.MN <- function(
  delta0, dataset, weight = "MN", conflev,
  bias.correct=FALSE
  ) {
  w = MN.weights(delta0,dataset,bias.correct)
  w = switch(
    weight,
    "sample" = rowSums(dataset[,c(2,4)]),
    "Cochran" = apply(dataset[,c(2,4)],1,prod) / rowSums(dataset[,c(2,4)]),
    #w= n1*n2/(n1+n2)
    "MH" = apply(dataset[,c(2,4)],1,prod) / (rowSums(dataset[,c(2,4)])-1),
    #w =n1*n2/(n1+n2-1)

    "inv.Var" = {},
    "Add4" = {},
    "Add2" = {
      p1=dataset[,1]/dataset[,2]
      p2=dataset[,3]/dataset[,4]
      (p1*(1-p1)/dataset[,2] +
          p2*(1-p2)/dataset[,4])^(-1)
      p1=(dataset[,1]+1)/(dataset[,2]+2)
      p2=(dataset[,3]+1)/(dataset[,4]+2)
      (p1*(1-p1)/dataset[,2] + p2*(1-p2)/dataset[,4])^(-1)
      p1=(dataset[,1]+0.5)/(dataset[,2]+1)
      p2=(dataset[,3]+0.5)/(dataset[,4]+1)
      (p1*(1-p1)/dataset[,2] + p2*(1-p2)/dataset[,4])^(-1)
      },
    MN.weights(delta0,dataset,bias.correct)
    )
  helpi =
    apply(
      dataset, 1,
      function(ro)
        resMLE.MN1(
          delta0, y=ro[c(1,3)], n=ro[c(2,4)]
          )
      )
  chisqf = sum(w*helpi[1,])^2 / sum(w^2*helpi[2,])
  return(chisqf - qchisq(conflev,df=1))
}
