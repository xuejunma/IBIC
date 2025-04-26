rm(list=ls())
library(grpreg)
library(MASS)
library(doParallel)
library(foreach)



lamsel <- function(x, y, lam=0.1){
  n  <- dim(x)[1]
  p  <- dim(x)[2]
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


#############################under correct overfitted function
whix <- function(x,index.123) {
  zx  <- which(abs(x) > 0) 
  udf <- 1- all(index.123 %in% zx )    # under-fit
  crf <- setequal(index.123, zx) # correct-fit
  if(udf==1){
    orf <- 0
  }else{
    orf <-  1- crf # # over-fit
  }
  
  return(c(udf, crf, orf))
}



########################################screening begin
SIS.my <- function (x, y, d){
  p <- dim(x)[2]
  n <- dim(x)[1]
  a <- apply(x, 2, function(x)(cor(x,y=y)))
  
  w.sis <- abs(a)
  M.sis <- order(w.sis, decreasing = T)[1:d]
  return(list(w = w.sis, M = M.sis))
}

##################################### screening end
myf <- function(x){
  min(which(x==min(x)))
}



Data1<- function(n, pp, rho=0.5)
{
  sigm <- rho ^ abs(outer(1:(pp),1:(pp),"-"))
  muz <- rep(0, pp)
  x <- mvrnorm(n = n, mu=muz, Sigma=sigm)
  e <- rnorm(n, 0, 1)
  ############
  yture <- 5 * x[, 1] + 5 * x[, 2] + 5 * x[, 3]
  y <- yture + e
  dat = list(y=y, x=x)
  return(dat)    
}
###
cl <- makeCluster(8)
registerDoParallel(cl)

n <- 70  ## 50
pp <- 1000 ## 1000 200
rho <- 0   ## 0 0.25

m   <- 500 
d <- n-1


labname <- c("GIC", "EBIC", "HBIC", "IBIC", "MBIC","ERIC", "AIC", "BIC","GCV")
nlab <- length(labname)


####under correct and over 
undfit <- matrix(0, m, nlab)
colnames(undfit) <- labname
corfit = ovefit = undfit


####### B
fgl <- NULL

for(i in 1:m){
  dat <- Data1(n = n, p = pp, rho = rho)
  xx <- dat$x
  yy <- dat$y
  fit.sis <- SIS.my(x=xx, y=yy, d=d)
  inde <- fit.sis$M
  fgl[i] <- sum(inde %in% c(1,2,3))==3

  beta.ture  <- rep(0, n-1)
  index.1 <- which(inde==1)
  beta.ture[index.1] <-5
  index.2 <- which(inde==2)
  beta.ture[index.2] <-5
  index.3 <- which(inde==3)
  beta.ture[index.3] <-5
  index.123 <- c(index.1,index.2,index.3)
  p <- n - 1
  x <- xx[, inde]
  y <- yy
  beta.ture1 <- matrix(rep(beta.ture, nlab), nlab, p, byrow = T)
  ####################################### my begin
  nlam <- 100
  lam.min <- 0.05
  lam.max <- 1.5
  lam.v <- seq(from=lam.min, to=lam.max, length=nlam)
  pk.lam <- foreach(lam=lam.v, .combine="rbind", .packages = "grpreg" ) %dopar% {lamsel(y=y, x=x, lam=lam)}
  ind.lam <- apply(pk.lam, 2, myf)
  coef.my <- matrix(0, (nlab-3), (p+1))
  rownames(coef.my) <- names(ind.lam)
  
  
  #########################check
  print( lam.v[ind.lam] )
  par(mfrow=c(3,2))
  for(tt in 1:6){
    plot(pk.lam[, tt])
    title(names(ind.lam)[tt])
  }
  ########################check 
  
  for(kk in 1:(nlab-3)){
    inde <- as.numeric(  unlist(ind.lam[kk])  ) 
    fit.my <- grpreg(X=x, y=y, penalty = "grSCAD",lambda = lam.v[ inde ] )
    coef.my[kk, ] <- fit.my$beta
  }
  
  ####################################### my end
  
  ############################# package begin
  coef.it <- matrix(0, 3, (p+1))
  rownames(coef.it) <-c("AIC", "BIC", "GCV")
  fit <- grpreg(X=x, y=y, penalty = "grSCAD")
  coef.it[1, ] <- grpreg::select(fit,"AIC")$beta
  coef.it[2, ] <- grpreg::select(fit,"BIC")$beta
  coef.it[3, ] <- grpreg::select(fit,"GCV")$beta
  ############################# package end
  
  coef.full <- rbind(coef.my, coef.it )


  ####under correct over fitted
  fit.uco <- apply(coef.full[, -1], 1, whix, index.123=index.123)
  
  undfit[i, ] <- fit.uco[1, ]
  corfit[i, ] <- fit.uco[2, ]
  ovefit[i, ] <- fit.uco[3, ]
  
  print(i)
}
stopCluster(cl)

###B
mean(fgl)



uco <- rbind(apply(undfit,2, mean),apply(corfit, 2, mean),apply(ovefit, 2, mean))
rownames(uco) <- c("under","correct","over")
uco



write.csv(t(uco), file="D:/ICuco.csv") 

print(rho)
print(n)
print(pp)

