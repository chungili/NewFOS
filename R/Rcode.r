#' Data Preprocessing
#'
#' @param x Input data vector
#' @param trim logical; if TRUE, the fraction (0 to 0.01) of observations to be trimmed
#'
#' @return The preprocessed data
#' @export
#'
prepro = function(x, trim=TRUE){
  #filters out the missing values
  xfull = x[!is.na(x)]
  if (trim) {
    xfinal = DescTools::Trim(xfull, trim = 0.01)
  } else{
    xfinal = xfull
  }
  return(xfinal)
}

#' The classic FOS method
#'
#' @param x Input data vector
#' @param alpha False alarm rate
#'
#' @return the lower and upper confidence bounds ans pn value.
#' \itemize{
#'  \item{ limits:}{ The lower and upper confidence limit}
#'  \item{ pn :}{ The nominal coverage probability}
#'  }
#' @export
#'
#' @examples
#' FOS(x1, alpha = 0.0027)
FOS = function(x, alpha = 0.0027){
  n = length(x)
  Rx <- sort(x)
  Pn_FOS = function(n){
    if (n > 1851) {0.95
    } else if (n > 1481) {0.90
    } else if (n > 1111) {0.80
    } else if (n > 740)  {0.70
    } else if (n > 370)  {0.50
    } else {0.3}
  }
  limit = function(u){
    n1 = 1+n
    if (0<u && u<=1/n1) {
      Rx[1]+(Rx[2]-Rx[1])*log((n1)*u)
    } else if (1/n1<u && u<n/n1) {
      (1-(n1*u-floor(n1*u)))*Rx[floor(n1*u)]+(n1*u-floor(n1*u))*Rx[floor(n1*u)+1]
    } else {
      Rx[n]-(Rx[n]-Rx[n-1])*log(n1*(1-u))
    }
  }
  fun1 = function(r){pbeta(alpha/2, r, n-r+1)-(1+Pn_FOS(n))/2}
  fun2 = function(r){pbeta(1-alpha/2, r, n-r+1)-(1-Pn_FOS(n))/2}
  kl = uniroot(fun1, c(0, n+1))$root
  ku = uniroot(fun2, c(0, n+1))$root
  ur = kl/(n+1)
  us = ku/(n+1)
  pn = Pn_FOS(n)
  LCL = limit(u=ur)
  UCL = limit(u=us)
  limits = cbind(LCL, UCL)
  return(list(limits=limits, pn=pn)
}

#' The adaptive FOS method
#'
#' @param x Input data vector
#' @param alpha False alarm rate
#' @param pn Nominal coverage probability
#'
#' @return the lower and upper confidence bounds
#' \itemize{
#'  \item{ LCL:}{ The lower confidence limit}
#'  \item{ UCL:}{ The upper confidence limit}
#'  }
#' @export
#'
#' @examples
#' FOS_ad(x6, alpha = 0.0027, pn = 0.3)
FOS_ad = function(x, alpha = 0.0027, pn = 0.9){
  n = length(x)
  Rx <- sort(x)
  n1 = n + 1
  m = function(pn){
    if (0.5<pn && pn<=0.7){370}
    else if (0.7<=pn && pn<=0.8){740}
    else if (0.8<pn && pn<=0.95){1111}
  }
  sigma = function(pn){
    if (0.5<pn && pn<=0.7){2}
    else if (0.7<pn && pn<=0.8){1}
    else if (0.8<pn && pn<=0.95){0.5}
  }
  pp0 = function(n){
    if(0<pn && pn<=0.5){1}
    else {pnorm(n/m(pn)-1, mean =0, sd =sigma(pn))}
  }
  limit = function(u){
    if (0<u && u<=1/(pp0(n)*n1)) {
      Rx[1]+(Rx[2]-Rx[1])*log(pp0(n)*n1*u)-(Rx[2]-Rx[1])*(log(pp0(n)*n1*u))^2
    } else if (1/(n1*pp0(n))<u && u<1-1/(n1*pp0(n))) {
      (1-(n1*u-floor(n1*u)))*Rx[floor(n1*u)]+(n1*u-floor(n1*u))*Rx[floor(n1*u)+1]
    } else {
      Rx[n]-(Rx[n]-Rx[n-1])*log(pp0(n)*n1*(1-u))+(Rx[n]-Rx[n-1])*(log(pp0(n)*n1*(1-u)))^2
    }
  }
  fun1 = function(r){pbeta(alpha/2, r, n-r+1)-(1+pn)/2}
  fun2 = function(r){pbeta(1-alpha/2,r, n-r+1)-(1-pn)/2}
  kl <- uniroot(fun1,c(0, n+1))$root
  ku <- uniroot(fun2,c(0, n+1))$root
  ur = kl/(n+1)
  us = ku/(n+1)
  LCL = limit(u=ur)
  UCL = limit(u=us)
  limits = cbind(LCL, UCL)
  return(limits)
}

#' The 3-term FOS method
#'
#' @param x Input data vector
#' @param alpha False alarm rate
#' @param pn Nominal coverage probability
#'
#' @return The lower and upper confidence bounds
#' \itemize{
#'  \item{ LCL:}{ The lower confidence limit}
#'  \item{ UCL:}{ The upper confidence limit}
#'  }
#' @export
#'
#' @examples
#' FOS_3terms(x6, alpha = 0.0027, pn = 0.3)
FOS_3terms <- function(x, alpha = 0.0027, pn = 0.9){
  n = length(x)
  Rx = sort(x)
  n1 = n + 1
  limit = function(u){
    if (0<u && u<=1/(n1)) {
      Rx[1]+(Rx[2]-Rx[1])*log(n1*u)-(Rx[2]-Rx[1])*(log(n1*u))^2+(Rx[2]-Rx[1])*(log(n1*u))^3
    } else if (1/n1<u && u<1-1/n1) {
      (1-(n1*u-floor(n1*u)))*Rx[floor(n1*u)]+(n1*u-floor(n1*u))*Rx[floor(n1*u)+1]
    } else {
      Rx[n]-(Rx[n]-Rx[n-1])*log(n1*(1-u))+(Rx[n]-Rx[n-1])*(log(n1*(1-u)))^2-(Rx[n]-Rx[n-1])*(log(n1*(1-u)))^3
    }
  }
  fun1 = function(r){pbeta(alpha/2, r, n-r+1)-(1+pn)/2}
  fun2 = function(r){pbeta(1-alpha/2, r, n-r+1)-(1-pn)/2}
  kl = uniroot(fun1,c(0, n+1))$root
  ku = uniroot(fun2,c(0, n+1))$root
  ur = kl/(n+1)
  us = ku/(n+1)
  LCL = limit(u=ur)
  UCL = limit(u=us)
  limits = cbind(LCL, UCL)
  return(limits)
}

#' The ad-AF method involves fine-tuning parameters σ and m based on chosen parameter distributions.
#'
#' @param Sd Standard derivation for the CDF of normal distribution: search domains \eqn{\sigma \in \{0.5, 1, 2\}}
#' @param m The inflection points for the CDF of normal distribution: search domains \eqn{ m \in \{370, 740, 1111\}}
#' @param n Sample size
#' @param pn Nominal coverage probability
#' @param dist The trial-run distributions are: Normal(0,1), Chi-square(4), t(4), Beta(9,1), and lognormal(0,1) represented by c("normal", "chisq", "t", "beta", "lnorm"), respectively.
#' @param alpha False alarm rate
#' @param B Number of simulation replications
#'
#' @return A list of the optimal tuning parameters
#' \itemize{
#'  \item{ Opt.m:}{ The optimal inflection point}
#'  \item{Opt.sd:}{ The optimal standard derivation}
#'  }
#' @export
#'
#' @examples
#' FineTune(n = 370, pn = 0.6, dist="t")
FineTune = function(Sd = c(0.5, 1, 2), m = c(370, 740, 1111), n, pn,
                     dist = c("normal", "chisq", "t", "beta", "lnrom"),
                     alpha = 0.0027, B = 10000){
  n.m = length(m)
  n.Sd = length(Sd)
  FOS = function(x, pn, sd, m, alpha = alpha){
    Rx = sort(x)
    n = length(x)
    n1 = n + 1
    p0 = pnorm(n/m-1, mean = 0, sd = sd)
    Xu = function(u){
      logwl = log(p0*n1*u)
      logwu = log(p0*n1*(1-u) )
      if (0<u && u<=1/(p0*n1)) {
        Rx[1]+(Rx[2]-Rx[1])*logwl-(Rx[2]-Rx[1])*(logwl)^2
      } else if (1/(n1*p0)<u && u<1-1/(n1*p0)) {
        (1-(n1*u-floor(n1*u)))*Rx[floor(n1*u)]+(n1*u-floor(n1*u))*Rx[floor(n1*u)+1]
      } else {
        Rx[n]-(Rx[n]-Rx[n-1])*logwu+(Rx[n]-Rx[n-1])*(logwu)^2
      }
    }
    findkl = function(r, n, alpha = alpha, pn = pn ) pbeta(alpha/2, r, n-r+1)-(1+pn)/2
    findku = function(r, n, alpha = alpha, pn = pn ) pbeta(1-alpha/2, r, n-r+1)-(1-pn)/2
    kl <- uniroot(findkl, c(0, n1), n = n, alpha = alpha, pn = pn)$root
    ku <- uniroot(findku, c(0, n1), n = n, alpha = alpha, pn = pn)$root
    ur = kl/n1
    us = ku/n1
    LCL = Xu(ur)
    UCL = Xu(us)
    return(list(UCL=UCL, LCL=LCL))
  }
  library(doParallel)
  cl = makeCluster(detectCores()-2)
  registerDoParallel(cl)
  result = foreach( m = rep(m, n.Sd), Sd = rep(Sd, each=n.m),
                    .combine='cbind') %:%
    foreach (1:B, .combine = c) %dopar% {
      x0 = switch(dist,
                  "normal" = rnorm(n, mean = 0, sd = 1),
                  'chisq' = rchisq(n, df = 4),
                  "t" = rchisq(n, df = 4),
                  "beta" = rbeta(n, shape1 = 9, shape2 = 1),
                  "lnorm" = rlnorm(n, meanlog = 0, sdlog = 1))
      x1 = switch(dist,
                  "normal" = rnorm(n, mean = 0, sd = 1),
                  'chisq' = rchisq(n, df = 4),
                  "t" = rchisq(n, df = 4),
                  "beta" = rbeta(n, shape1 = 9, shape2 = 1),
                  "lnorm" = rlnorm(n, meanlog = 0, sdlog = 1))
      FOS_limits = FOS(x0, pn = pn, sd = Sd, m = m, alpha = alpha)
      mean( x1 > FOS_limits$UCL | x1 < FOS_limits$LCL ) < alpha
    }
  stopCluster(cl)
  pa = matrix( colMeans(result), nrow = n.Sd, ncol = n.m, byrow = T)
  loc = which(pa==min(pa), arr.ind = T)
  return(list(Opt.m = m[loc[2]], Opt.sd = Sd[loc[1]]) )
}

#' The ad-AF method: fine-tune parameters σ and m based on the given training and testing data.
#'
#' @param Sd Standard derivation for the CDF of normal distribution: search domains \eqn{\sigma \in \{0.5, 1, 2\}}
#' @param m The inflection points for the CDF of normal distribution: search domains \eqn{ m \in \{370, 740, 1111\}}
#' @param pn Nominal coverage probability
#' @param Train The provided training data
#' @param Test The provided testing data
#' @param alpha False alarm rate
#' @param B Number of simulation replications
#'
#' @return A list of the optimal tuning parameters
#' \itemize{
#'  \item{ Opt.m:}{ The optimal inflection point}
#'  \item{Opt.sd:}{ The optimal standard derivation}
#'  }
#' @export
#'
#' @examples
#' Train = matrix(rnorm(10000*370), nrow = 10000 )
#' Test = matrix(rnorm(10000*370), nrow = 10000 )
#' FineTuneWithData(Train = Train, Test=Test, pn = 0.7)
FineTuneWithData = function(Sd = c(0.5, 1, 2), m = c(370, 740, 1111), pn, Train, Test,
                            alpha = 0.0027, B = 10000){
  n = dim(Train)[2]
  n.m = length(m)
  n.Sd = length(Sd)
  FOS = function(x, pn, sd, m, alpha = alpha){
    Rx = sort(x)
    n = length(x)
    n1 = n + 1
    p0 = pnorm(n/m-1, mean = 0, sd = sd)
    Xu = function(u){
      logwl = log(p0*n1*u)
      logwu = log(p0*n1*(1-u) )
      if (0<u && u<=1/(p0*n1)) {
        Rx[1]+(Rx[2]-Rx[1])*logwl-(Rx[2]-Rx[1])*(logwl)^2
      } else if (1/(n1*p0)<u && u<1-1/(n1*p0)) {
        (1-(n1*u-floor(n1*u)))*Rx[floor(n1*u)]+(n1*u-floor(n1*u))*Rx[floor(n1*u)+1]
      } else {
        Rx[n]-(Rx[n]-Rx[n-1])*logwu+(Rx[n]-Rx[n-1])*(logwu)^2
      }
    }
    findkl = function(r, n, alpha = alpha, pn = pn ) pbeta(alpha/2, r, n-r+1)-(1+pn)/2
    findku = function(r, n, alpha = alpha, pn = pn ) pbeta(1-alpha/2, r, n-r+1)-(1-pn)/2
    kl <- uniroot(findkl, c(0, n1), n = n, alpha = alpha, pn = pn)$root
    ku <- uniroot(findku, c(0, n1), n = n, alpha = alpha, pn = pn)$root
    ur = kl/n1
    us = ku/n1
    LCL = Xu(ur)
    UCL = Xu(us)
    return(list(UCL=UCL, LCL=LCL))
  }
  library(doParallel)
  cl <- makeCluster(detectCores()-2)
  registerDoParallel(cl)
  result <- foreach( m = rep(m, n.Sd), Sd = rep(Sd, each=n.m),
                     .combine='cbind') %:%
    foreach (1:B, .combine = c) %dopar% {
      locTrain = sample(x = 1:B, size = 1, replace = TRUE)
      locTest = sample(x = 1:B, size = 1, replace = TRUE)
      x0 = Train[locTrain, ]
      x1 = Test[locTest, ]
      FOS_limits <- FOS(x0, pn = pn, sd = Sd, m = m, alpha = alpha)
      mean( x1 > FOS_limits$UCL | x1 < FOS_limits$LCL ) < alpha
    }
  stopCluster(cl)
  pa = matrix( colMeans(result), nrow = n.Sd, ncol = n.m, byrow = T)
  pa[pa-pn < 0]=1
  loc = which(pa==min(pa), arr.ind = T)
  return(list(Opt.m = m[loc[2]], Opt.sd = Sd[loc[1]]) )
}
