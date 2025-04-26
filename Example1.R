rm(list=ls())
library(grpreg)
library(MASS)
library(doParallel)
library(foreach)
Dat1 <- function(n , p, sig0=0.5, rho=1){
  sigm <- sig0^abs(outer(1:(p),1:(p),"-"))
  muz <- rep(0,p)
  x <- mvrnorm(n = n, mu=muz, Sigma=sigm)
  e <- rnorm(n,0,1)
  beta <-c(3, 1.5, 0, 0, 2, rep(0, p-5))
  y <-  x %*% beta + rho * e
  return(list(y=y, x=x, beta=beta))
}

lamsel <- function(x, y, lam=0.1){
  n <- dim(x)[1]
  p <- dim(x)[2]
  pk <- matrix(0, 1, 6)
  colnames(pk) <- c("GIC", "EBIC", "HBIC", "IBIC", "MBIC","ERIC")
  ##fitted begin
  fit <- grpreg(X=x, y=y, penalty = "grSCAD", lambda=lam)
  coef <- fit$beta[-1]
  k <- sum(coef!=0)
  res <- y - x %*% coef- fit$beta[1]
  ##fitted end
  
  sigm <- sum(res^2) / n
  lsigm  <- log(sigm)
  
  pk[1, "EBIC"] <- lsigm +  (log(n) +  2 * 1 * k * log(p)) /n  #EBIC  
  pk[1, "HBIC"] <- lsigm +  (2 * 1.25 * k * log(p)) / n   # HBIC 
  pk[1, "IBIC"] <- lsigm +   k * (n + p) / (n * p) * log( p ^ 2) #IBIC
  pk[1, "MBIC"] <- lsigm +   k * log(n) / n * log(log(p)) # MBIC
  pk[1, "ERIC"] <- lsigm +  2 * 1 * k * log(n * sigm/ lam) / n ##ERIC
  pk[1, "GIC"]  <- lsigm +  k * log(log(n)) * log(p) /n #GIC
  
  return(pk)
}


myf <- function(x){ 
  min(which(x==min(x)))
}

#############################under correct overfitted function begin
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
#############################under correct overfitted function end



cl <- makeCluster(8)
registerDoParallel(cl)

p    <- 8 ## 8 60
n    <- 50 ## 50 70
rho  <- 2 ## 1 2
sig0 <- 0 ## 0 0.25


m <- 500
labname <- c("GIC", "EBIC", "HBIC", "IBIC", "MBIC","ERIC","AIC", "BIC","GCV")
nlab <- length(labname)

####under correct and over 
undfit <- matrix(0, m, nlab)
colnames(undfit) <- labname
corfit = ovefit =undfit


for(i in 1:m){
  dat <- Dat1(n = n, p = p, rho = rho,sig0 = sig0)
  x <- dat$x
  y <- dat$y
  beta.ture <- dat$beta
  beta.ture1 <- matrix(rep(beta.ture, nlab), nlab, p, byrow = T)
  ######
  nlam <- 100
  lam.min <- 0.05
  lam.max <- 0.8
  lam.v   <- seq(from=lam.min, to=lam.max, length=nlam)
  pk.lam  <- foreach(lam=lam.v, .combine="rbind", .packages = "grpreg" ) %dopar% {lamsel(y=y, x=x, lam=lam)}
  ind.lam <- apply(pk.lam, 2, myf)
  coef.my <- matrix(0, (nlab-3), (p+1))
  rownames(coef.my) <- names(ind.lam)
  
  #########################check begin
  print( lam.v[ind.lam] )
  par(mfrow=c(3,2))
  for(tt in 1:6){
    plot(pk.lam[, tt])
    title(names(ind.lam)[tt])
  }
  ########################check end
  
  for(k in 1:(nlab-3)){
    inde <- as.numeric(ind.lam[k])
    fit.my <- grpreg(X=x, y=y, penalty = "grSCAD",lambda = lam.v[ inde ] )
    coef.my[k, ] <- fit.my$beta
  }
  #######################################
  #############################  package begin
  coef.it <- matrix(0, 3, (p+1))
  rownames(coef.it) <-c("AIC", "BIC", "GCV")
  fit <- grpreg(X=x, y=y, penalty = "grSCAD")
  coef.it[1, ] <- grpreg::select(fit,"AIC")$beta
  coef.it[2, ] <- grpreg::select(fit,"BIC")$beta
  coef.it[3, ] <- grpreg::select(fit,"GCV")$beta
  #############################  package end
  coef.full  <- rbind(coef.my, coef.it )
  coef.full1 <- coef.full[,-1]
  
  ####under correct over 
  fit.uco <- apply(coef.full1, 1, whix)
  undfit[i, ] <- fit.uco[1,]
  corfit[i, ] <- fit.uco[2,]
  ovefit[i, ] <- fit.uco[3,]
  
  print(i)
}

stopCluster(cl)


uco <- rbind(apply(undfit,2, mean),apply(corfit, 2, mean),apply(ovefit, 2, mean))
rownames(uco) <- c("under","correct","over")
uco


write.csv(t(uco), file="D:/ICuco.csv") 


print(sig0)
print(n)
print(p)
