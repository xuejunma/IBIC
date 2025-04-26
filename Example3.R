rm(list=ls())
library(glmnet)
library(MASS)
library(doParallel)
library(foreach)
Dat1 <- function(n , p, sig0=0.5){
  sigm <- sig0^abs(outer(1:(p),1:(p),"-"))
  muz <- rep(0,p)
  x <- mvrnorm(n = n, mu=muz, Sigma=sigm)
  ####
  beta <-c(3, 1.5, 0, 0, 2, rep(0, p-5))
  xbeta <-  x %*% beta 
  py <- exp(xbeta) / (1+exp(xbeta))
  y <- rbinom(n, 1, prob=py)
  return(list(y=y, x=x, beta=beta))
}

lamsel <- function(x, y, lam=0.1){
  n <- dim(x)[1]
  p <- dim(x)[2]
  pk <- matrix(0, 1, 4)
  colnames(pk) <- c("EBIC", "IBIC", "ERIC", "GIC")
  ##
  fit  <- glmnet::glmnet(x=x, y=y, family = "binomial", lambda = lam)
  coef <- fit$beta
  k    <- sum(coef!=0)
  xbeta_fit <-  x %*% coef + fit$a0
  p_fit  <- 1 / (1 + exp(-xbeta_fit))
  loglik <- p_fit ^ y * (1 - p_fit) ^ (1 - y)
  ## 
  
  lsigm  <- -2 * sum( log(loglik) )
  pk[1, "EBIC"] <- lsigm +  k * log(n) +  2 * 1 * k * log(p)  #EBIC  Chen and Chen 2012
  pk[1, "IBIC"] <- lsigm +  k * (n + p) / ( p) * log( p ^ 2) #IBIC
  pk[1, "ERIC"] <- lsigm +  2 * 1 * k * log(n/lam)   # ERIC  
  pk[1, "GIC"]  <- lsigm +  k * log(log(n)) * log(p)  #GIC
  return(pk)
}


myf <- function(x){ 
  min(which(x==min(x)))
}

#############################under correct overfitted function
whix <- function(x) {
  zx  <- which(abs(x) > 0) 
  udf <- 1- all(c(1, 2, 5) %in% zx )    # under-fit
  crf <- setequal(c(1, 2, 5), zx) # correct-fit
  if(udf==1){
    orf <- 0
  }else{
    orf <-  1- crf # # over-fit
  }
  return(c(udf, crf, orf))
}


cl <- makeCluster(8)
registerDoParallel(cl)
p <- 200 #100 200
n <- 300 #300 500

sig0 <- 0


m <- 500

labname <- c("EBIC", "IBIC", "ERIC", "GIC")
nlab <- length(labname)

####under correct and over fit
undfit <- matrix(0, m, nlab)
colnames(undfit) <- labname
corfit = ovefit =undfit


for(i in 1:m){
  dat <- Dat1(n = n, p = p,sig0 = sig0)
  x <- dat$x
  y <- dat$y
  ####################################### my begin
  nlam <- 100
  lam.min <- {if (nrow(x) > ncol(x)) 5*1e-4 else .05} ###le-4
  lam.max <- 0.4
  lam.v <- seq(from=lam.min, to=lam.max, length=nlam)
  pk.lam <- foreach(lam=lam.v, .combine="rbind", .packages = "glmnet" ) %dopar% {lamsel(y=y, x=x, lam=lam)}
  ind.lam <- apply(pk.lam, 2, myf)
  coef.full <- matrix(0, nlab, p)
  rownames(coef.full) <- names(ind.lam)
  
  ######################
  print(lam.v[ind.lam])
  par(mfrow=c(2,2))
  for(tt in 1:4){
    plot(pk.lam[, tt])
    title(names(ind.lam)[tt])
  }
  
  ############# check 
  for(k in 1:(nlab)){
    inde <- as.numeric(ind.lam[k])
    fit.my <- glmnet::glmnet(x=x, y=y, family = "binomial", lambda =  lam.v[ inde ])
    coef.full[k, ] <-  coef(fit.my)[-1]
  }
  ###################check 
  ####################################### my end
  
  ####under correct over fitted
  fit.uco <- apply(coef.full, 1, whix)
  undfit[i, ] <- fit.uco[1,]
  corfit[i, ] <- fit.uco[2,]
  ovefit[i, ] <- fit.uco[3,]
  

  print(i)
}

stopCluster(cl)



uco <- rbind(apply(undfit,2, mean),apply(corfit, 2, mean),apply(ovefit, 2, mean))
rownames(uco) <- c("under","correct","over")
uco



write.csv(t(uco), file="D:\\ICuco.csv") 
print(sig0)
print(n)
print(p)
