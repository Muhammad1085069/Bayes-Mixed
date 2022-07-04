data_generate <- function(n, d1,lz,ff) {
  x0 <- randomLHS(n,d1)
  z <- as.matrix(gen.factorial(lz,length(lz), center=F))
  x <- matrix(x0,nrow(z)*n,ncol(x0),byrow=T)
  TT <- rep(1:nrow(z),each=n)
  z <- z[TT,]
  
  y <- ff(x,z)
  
  data <- NULL
  data$x0 <- x0
  data$x <- x
  data$z <- z
  data$TT <- TT
  data$y <- y
  data$lf <- matrix(1, nrow(x), 1)
  data$lfm <- model.matrix(y~as.factor(TT))
  data
} 

borehole<-function(x,z) {
  if(is.vector(x)) {
    x=matrix(x,1,length(x))
    z=matrix(z,1,length(z))
    
  }
  rw=x[,1]*(0.15-0.05)+0.05
  Tu=x[,2]*(115600-63070)+63070
  Tl=(x[,3])*(116-63.1)+63.1
  L=(x[,4])*(1680-1120)+1120
  
  r=z[,1]*20000/(max(z[,1])-1)
  H=z[,2]*200/(max(z[,2])-1)
  Kw=(z[,3]-1)*2000/(max(z[,3])-1)+10000
  
  
  
  y=2*pi*Tu*H
  y=y/log(r/rw)
  y1=1+Tu/Tl+(2*L*Tu)/(Kw*rw^2*log(r/rw))
  y=y/y1
  
  log(y)
}
# Exchangeable
maxlik_ec<-function(data, kf) {
  d1 <- ncol(data$x)
  
  #correlation parameter
  theta <- rep(0.5,d1)
  
  #bounds of the correlation parameter
  thetaup <- theta*100
  thetalo <- theta*0.01
  idtr <- data$idtr
  x <- as.matrix(data$x[idtr,])
  y <- data$y[idtr]
  TT <- data$TT[idtr]
  lf <- as.matrix(data$lfm[idtr,])
  lik0 <- function(param) lik_cc(x,y,lf,TT, param)
  opts <- list( "algorithm" = "NLOPT_LN_NELDERMEAD", 
                "xtol_rel" = 1.0e-3, 
                "maxeval" = 5)
  param.opt <- nloptr(c(0.5, theta), lik0, 
                      eval_grad_f = NULL, 
                      lb = c(0, thetalo), 
                      ub =c(0.999, thetaup),
                      opts=opts)$solution
  thetaopt <- param.opt[-1]
  tauopt <- param.opt[1]
  
  
  x <- x%*%diag(x=thetaopt, nrow=length(theta), ncol=length(theta))
  R <- kernelMatrix(kf, x)
  Tau0 <- matrix(tauopt, max(TT), max(TT))
  diag(Tau0) <- 1
  Tau <- Tau0[TT,TT]
  R <- R*Tau
  invR <- solve(R)
  l <- log(det(R))
  FinvR<-t(lf)%*%invR
  sigma2 <- t(y)%*%(invR-t(FinvR)%*%solve(FinvR%*%lf)%*%FinvR)%*%y/length(y)
  maxlikre <- NULL
  maxlikre$thetaopt <- thetaopt
  maxlikre$Sigma0 <- Tau0*as.vector(sigma2)
  maxlikre
}

lik_cc<-function(x, y, lf, TT,  param) {
  theta <- param[-1]
  tau <- param[1]
  x <- x%*%diag(x=theta, nrow=length(theta), ncol=length(theta))
  R <- kernelMatrix(kf, x)
  Tau <- matrix(tau, max(TT), max(TT))
  diag(Tau) <- 1
  Tau <- Tau[TT,TT]
  R <- R*Tau
  invR<-solve(R)
  l<-log(det(R))
  FinvR<-t(lf)%*%invR
  l+log(t(y)%*%(invR-t(FinvR)%*%solve(FinvR%*%lf)%*%FinvR)%*%y/length(y))*length(y)
}
  # Main code
  rm(list=ls(all=TRUE))
  install.packages("SparseM")
library(SparseM)
  install.packages("kernlab")
  library(kernlab)
  install.packages("nloptr")
  library(nloptr)
  install.packages("spcov")
  library(spcov)
  install.packages("lhs")
  library(lhs)
  install.packages("AlgDesign")
  library(AlgDesign)
  install.packages("knitr")
 library(knitr)
  source("data_generate.r")
  source("ic.r")
  source("rc.r")
  source("ec.r")
  source("ec_sep.r")
  source("uc.r")
  source("pred.r")
  
  
  set.seed(10)
  #initial setting
  #number of total sample
  n <- 4
  n.test <- 100
  
  #dimension of continuous variable
  d1 <- 4
  
  #number of levels in the categorical variable
  lz <- c(3, 3, 3)      #case 1: k=27
  # lz <- c(2, 2, 10)   #case 2: k=40
  # lz <- c(2, 2, 15)   #case 3: k=60
  # lz <- c(2, 2, 20)   #case 4: k=80
  
  
  
  #function generate the data
  ff <- borehole
  
  #penalty parameter
  ll <- c(seq(0.001, 0.01, 0.003), 
          seq(0.02, 0.1, 0.02))
  
  #correlation function
  kf <- laplacedot(sigma =1)
  
  
  # test data generation
  data.test <- data_generate(n.test, d1, lz, ff)
  
  #data generation
  data <- data_generate(2*n, d1, lz, ff)
  id<-1:nrow(data$x)
  idtrs<-tapply(id,data$TT,function(u) sample(u,n))
  idtr<-c(do.call('cbind', idtrs))
  data$idtr <- idtr
  
  
  
  
  computing.times <- rep(0, 5)
  mse.results <- matrix(0, prod(lz), 5)
  
  # IC method
  computing.times[2] <- system.time(
    ncfit<-ic(data, kf)
  )[3]
  
  # making prediction
  data$thetaopt <- apply(ncfit$thetaopts, 2, median)
  sigma2s <- ncfit$sigma2s
  data$Sigma0 <- diag(sigma2s)
  haty <- pred(data, data$thetaopt, kf, data.test)
  mse.results[ ,2] <- tapply((haty-data.test$y)^2, data.test$TT, mean)
  
  
  # EC
  computing.times[4] <- system.time(
    ecfit<-maxlik_ec(data, kf)
  )[3]
  
  # making predictions 
  thetaopt <- ecfit$thetaopt
  data$Sigma0 <-  ecfit$Sigma0
  haty <- pred(data, thetaopt, kf, data.test)
  mse.results[ ,4] <- tapply((haty-data.test$y)^2, data.test$TT, mean)
  
  #EC_sep
  computing.times[5] <- system.time(
    ecsepfit<-maxlik_ec_sep(data, kf)
  )[3]
  
  # making predictions 
  thetaopt <- ecsepfit$thetaopt
  data$Sigma0 <-  ecsepfit$Sigma0
  haty <- pred(data, thetaopt, kf, data.test)
  mse.results[ ,5] <- tapply((haty-data.test$y)^2, data.test$TT, mean)
  
  # UC
  computing.times[3] <- system.time(
    ucfit <- maxlik_uc(data, kf)
  )[3]
  
  # making predictions
  thetaopt <- ucfit$thetaopt
  data$Sigma0 <-  ucfit$Sigma0
  haty <- pred(data, thetaopt, kf, data.test)
  mse.results[ ,3] <- tapply((haty-data.test$y)^2, data.test$TT, mean)
  
  # rc method
  
  rc <- function(data, kf, lambda, data.test) {
    d1 <- ncol(data$x)
    #correlation parameter
    theta <- rep(1,d1)
    
    #bounds of the correlation parameter
    thetaup <- theta*100
    thetalo <- theta*0.01
    idtr <- data$idtr
    S <- list()
    A <- diag(rep(1, nrow(data$x)))
    S$A <- A[,idtr]
    S$lf <- data$lfm[idtr,]
    S$y <- data$y[idtr]
    S$Sigma <- data$Sigma
    S$x0 <- data$x0
    S$Phi <- kernelMatrix(kf, 
                          S$x0%*%diag(x=data$thetaopt,
                                      nrow=length(theta),
                                      ncol=length(theta)))
    # Sigma.last <- S$Sigma
    # theta.last <- theta
    spcovre<-iss(S,lambda)
    S$Sigma<-spcovre$Sigma
    S$Phi<-kernelMatrix(kf, S$x0%*%diag(x=data$thetaopt,
                                        nrow=length(theta),
                                        ncol=length(theta)))
    data$Sigma0<-S$Sigma
    lik0 <- function(param) lik.theta(S, param)
    opts <- list( "algorithm" = "NLOPT_LN_NELDERMEAD", 
                  "xtol_rel" = 1.0e-3, 
                  "maxeval" = 5)
    data$thetaopt <- nloptr(data$thetaopt, lik0, 
                            eval_grad_f = NULL, 
                            lb = thetalo, 
                            ub =thetaup,
                            opts=opts)$solution
    
    
    haty <- pred(data, data$thetaopt, kf, data.test)
    
    ## compute 2-fold cross-validation mse
    n <- length(idtr)/max(data$TT) 
    ids <- sapply(unique(data$TT), 
                  function(u) sample(idtr[data$TT[idtr]==u], n/2))
    fold1 <- as.vector(ids)
    fold2 <- setdiff(idtr, fold1)
    data.cvtrain <- data
    data.cvtrain$idtr <- fold1
    data.cvtest <- data
    data.cvtest$x <- data$x[fold2, ]
    data.cvtest$TT <- data$TT[fold2]
    data.cvtest$lfm <- data$lfm[fold2, ]
    cvmse <- mean((pred(data.cvtrain, data$thetaopt, kf, data.cvtest)-data$y[fold2])^2)
    
    data.cvtrain <- data
    data.cvtrain$idtr <- fold2
    data.cvtest <- data
    data.cvtest$x <- data$x[fold1, ]
    data.cvtest$TT <- data$TT[fold1]
    data.cvtest$lfm <- data$lfm[fold1, ]
    cvmse <- cvmse + mean((pred(data.cvtrain, data$thetaopt, kf, data.cvtest)-data$y[fold1])^2)
    cvmse <- cvmse/2
    list(Sigma=S$Sigma, thetaopt=data$thetaopt, haty=haty, cvmse=cvmse)
  }
  
  lik.theta <- function(S, theta0) {
    S0<-S
    S0$Phi<-kernelMatrix(kf, S$x0 %*% diag(x=theta0,
                                           ncol=length(theta0),
                                           nrow=length(theta0)))
    ComputeLikelihood(S0$Sigma, S0)
  }
  
  iss<- function(S, lambda) {
    n.outer.steps <- 1e+2
    n.inner.steps <- 1e+2
    tol.outer <- 1e-3
    thr.inner <- 1e-3
    Sigma<-S$Sigma
    step.size<-1
    beta<-0.5
    mean.abs.S <- mean(abs(Sigma))
    del<-min(eigen(Sigma, symmetric=T, only.values=T)$val)
    likere<-ComputeLikelihood(Sigma,S)
    objective <- likere + sum(abs(lambda*Sigma))
    objectives<-objective
    n.iter <- NULL # number of inner iterations on each step
    for (i in seq(n.outer.steps)) {
      S$Sigma <- Sigma
      gg <- GGDescent(S, lambda,n.inner.steps,step.size,thr.inner * mean.abs.S,del,beta)
      Sigma <- gg$Sigma
      mean.abs.S <- mean(abs(Sigma))
      likere<-ComputeLikelihood(Sigma,S)
      objective <- likere + sum(abs(lambda*Sigma))
      objectives<-c(objectives,objective)
      n.iter <- c(n.iter, gg$niter)
      if((objectives[i + 1] - objectives[i])/abs(objectives[i]) >- tol.outer) {
        cat("-Majorization Maximization converged in", i, "steps!", fill=TRUE)
        break
      }
      if(max(gg$step.sizes)<step.size*beta^2) {
        step.size<-step.size*beta
      }
    }
    list(n.iter=n.iter, Sigma=gg$Sigma, obj=objectives)
  }
  
  GGDescent <- function(S, lambda, nsteps,step.size,tol,del,beta) {
    require(spcov)
    exit<-FALSE
    tt <-step.size
    tts<-c()
    Sigma<-S$Sigma
    gradre<-Computegrad(Sigma,S)
    g<-gradre$g
    tangent<-gradre$tangent
    for (i in seq(nsteps)) {
      Sigma.last<-Sigma
      if(!exit) {
        Sigma.temp<-ProxADMM(Sigma.last-tt*g,del,1,P=lambda*tt)$Z
        gen.grad<-(Sigma.last-Sigma.temp)/tt
        decent.d<-sum(gen.grad*g)
        right<-tangent-tt*decent.d+tt*sum(gen.grad^2)
        gradre.temp<-Computegrad(Sigma.temp,S)
        left<-gradre.temp$tangent
        if(decent.d>1e-3) {
          while((tt>1e-5)&(left>right)) {
            tt<-tt*beta
            Sigma.temp<-Sigma.last-gen.grad*tt
            gradre.temp<-Computegrad(Sigma.temp,S)
            left<-gradre.temp$tangent
            right<-tangent-tt*decent.d+tt*sum(gen.grad^2)
          }
        }
        tts<-c(tts,tt)
      }
      Sigma.temp<-ProxADMM(Sigma.last-tt*g,del,1,P=lambda*tt)$Z
      g<-gradre.temp$g
      tangent<-gradre.temp$tangent
      Sigma<-Sigma.temp
      if (mean(abs(Sigma - Sigma.last)) < tol) {
        cat(" --- Gradient Descent converged in", i, "steps!",fill=TRUE)
        break
      }
      
      
    }
    list(Sigma=Sigma, niter=i,step.sizes=tts)
  }
  
  
  ComputeLikelihood<-function(Sigma, S) {
    R<-Sigma%x%S$Phi
    R<-t(S$A)%*%R%*%S$A
    log.det.R<-sum(log(eigen(R+diag(rep(0.000001, nrow(R))))$values))
    invR<-solve(R+diag(rep(0.000001, nrow(R))))
    FinvR<-t(S$lf)%*%invR
    H<-invR-t(FinvR)%*%solve(FinvR%*%S$lf)%*%FinvR
    l1<-t(S$y)%*%H%*%S$y
    l2<-log.det.R
    l<-l1+l2
    l
  }
  
  Computegrad<-function(Sigma,S) {
    K<-nrow(Sigma)
    n<-nrow(S$Phi)
    R0<-S$Sigma%x%S$Phi
    R0<-t(S$A)%*%R0%*%S$A
    invR0<-solve(R0+diag(rep(0.000001, nrow(R0))))
    B<-S$A%*%invR0%*%t(S$A)
    R<-Sigma%x%S$Phi
    tangent<-sum(diag(B%*%R))
    tB<-matrix(0,K,K)
    id<-matrix(1:(n*K),n,K)
    for(i in 1:K) {
      for(j in 1:K) {
        tB[i,j]<-sum(diag(B[id[,i],id[,j]]%*%S$Phi))
      }
    }
    R<-t(S$A)%*%R%*%S$A
    invR<-solve(R+diag(rep(0.000001, nrow(R))))
    FinvR<-t(S$lf)%*%invR
    H<-invR-t(FinvR)%*%solve(FinvR%*%S$lf)%*%FinvR
    tangent<-tangent+t(S$y)%*%H%*%S$y
    ty<-matrix(S$A%*%H%*%S$y,n,K)
    g<-tB-t(ty)%*%S$Phi%*%ty
    return(list(g=g,tangent=tangent))
  }
  
  compute.theta<-function(S,kf,theta,thetaup,thetalo) {
    lik.theta <- function(theta0) {
      S0<-S
      S0$Phi<-kernelMatrix(kf, S$x0 %*% diag(x=theta0,
                                             ncol=length(theta0),
                                             nrow=length(theta0)))
      ComputeLikelihood(S0$Sigma, S0)
    }
    opts <- list( "algorithm" = "NLOPT_LN_NELDERMEAD", 
                  "xtol_rel" = 1.0e-3, 
                  "maxeval" = 5)
    thetaopt<-nloptr(theta,lik.theta, 
                     eval_grad_f = NULL, 
                     lb = thetalo, 
                     ub =thetaup,
                     opts=opts)$solution
    thetaopt
  }
  
  # RC
  Sigmas <- list()
  hatys <- NULL
  cvmses <- c()
  times <- c()
  data$Sigma0 <- diag(rep(1,nrow(data$Sigma0)))
  #scal.mat<- diag(1/sqrt(diag(ucfit$Sigma0)))
  #data$Sigma0 <- scal.mat %*% ucfit$Sigma0 %*% scal.mat
  
  for(i in 1:length(ll)) {
    cat("lambda equals", ll[i], fill=T)
    lambda <- ll[i]*matrix(1,nrow(data$Sigma0),nrow(data$Sigma0))
    lambda <- lambda-diag(diag(lambda))
    if(i>1) {
      data$Sigma0 <- Sigmas[[i-1]]
    }
    times[i] <- system.time(
      rcfit<-rc(data, kf, lambda, data.test)
    )[3]
    cvmses <- c(cvmses, rcfit$cvmse)
    Sigmas[[i]] <- rcfit$Sigma
    hatys <- cbind(hatys, rcfit$haty)
  }
  
  computing.times[1] <- sum(times)
  mse.results[,1] <- tapply((hatys[ , which.min(cvmses)]-data.test$y)^2, 
                            data.test$TT, mean)
  
  results <- NULL
  results <- rbind(results, apply(mse.results, 2, mean))
  results <- rbind(results, apply(mse.results, 2, sd))
  results <- rbind(results, computing.times/60)
  results <- data.frame(results)
  names(results) <- c("RC", "IC", "UC", "EC", "EC_sep")
  row.names(results) <-c("mean", "sd", "times")
  print(kable(results))
  
  #store the results
  save(computing.times, mse.results, cvmses, Sigmas,
       file=paste("/results/re",prod(lz),".RData",sep=""))
  
  
  install.packages("corrplot") 
  
  
  
  maxlik_ec_sep<-function(data, kf) {
    d1 <- ncol(data$x)
    
    #correlation parameter
    theta <- rep(1,d1)*0.5
    
    #bounds of the correlation parameter
    thetaup <- theta*100
    thetalo <- theta*0.01
    idtr <- data$idtr
    x <- as.matrix(data$x[idtr,])
    z <- as.matrix(data$z[idtr,])
    nf <- ncol(z)
    y <- data$y[idtr]
    TT <- data$TT[idtr]
    lf <- as.matrix(data$lfm[idtr,])
    lik0 <- function(param) lik_cc_sep(x, z, y, lf, TT,  param)
    opts <- list( "algorithm" = "NLOPT_LN_NELDERMEAD", 
                  "xtol_rel" = 1.0e-3, 
                  "maxeval" = 5)
    param.opt <- nloptr(c(rep(0.5, ncol(z)), theta), lik0, 
                        eval_grad_f = NULL, 
                        lb = c(rep(0, ncol(z)), thetalo), 
                        ub =c(rep(0.999, ncol(z)), thetaup),
                        opts=opts)$solution
    thetaopt <- param.opt[-c(1:nf)]
    tauopt <- param.opt[1:nf]
    
    
    x <- x%*%diag(x=thetaopt, nrow=length(thetaopt), ncol=length(thetaopt))
    R <- kernelMatrix(kf, x)
    Taus <- list()
    Taus[[1]] <- matrix(tauopt[nf], max(z[,1]), max(z[,1])) 
    diag(Taus[[1]]) <- 1
    Tau <- Taus[[1]]
    for(i in 2:nf) {
      Taus[[i]] <- matrix(tauopt[nf-i+1], max(z[,i]), max(z[,i])) 
      diag(Taus[[i]]) <- 1
      Tau <- kronecker(Tau, Taus[[i]])
    }
    Sigma0 <- Tau
    Tau <- Tau[TT,TT]
    R <- R*Tau
    invR<-solve(R)
    l<-log(det(R))
    FinvR<-t(lf)%*%invR
    sigma2 <- t(y)%*%(invR-t(FinvR)%*%solve(FinvR%*%lf)%*%FinvR)%*%y/length(y)
    maxlikre <- NULL
    maxlikre$thetaopt <- thetaopt
    maxlikre$Sigma0 <- Sigma0*as.vector(sigma2)
    maxlikre
  }
  
  lik_cc_sep<-function(x, z, y, lf, TT, param) {
    nf <- ncol(z)
    theta <- param[-c(1:nf)]
    tau <- param[1:nf]
    x <- x%*%diag(x=theta, nrow=length(theta), ncol=length(theta))
    R <- kernelMatrix(kf, x)
    Taus <- list()
    Taus[[1]] <- matrix(tau[nf], max(z[,1]), max(z[,1])) 
    diag(Taus[[1]]) <- 1
    Tau <- Taus[[1]]
    for(i in 2:nf) {
      Taus[[i]] <- matrix(tau[nf-i+1], max(z[,i]), max(z[,i])) 
      diag(Taus[[i]]) <- 1
      Tau <- kronecker(Tau, Taus[[i]])
    }
    Tau <- Tau[TT,TT]
    R <- R*Tau
    invR<-solve(R)
    l<-log(det(R))
    FinvR<-t(lf)%*%invR
    l+log(t(y)%*%(invR-t(FinvR)%*%solve(FinvR%*%lf)%*%FinvR)%*%y/length(y))*length(y)
  }
  x=as.matrix(xn,zp)
 data1=data.frame(x,nwx,y2,wa) 
  maxlik_ec_sep<-function(data1, kf) {
    d1 <- ncol(x)
    
    #correlation parameter
    theta <- rep(1,d1)*0.5
    
    #bounds of the correlation parameter
    thetaup <- theta*100
    thetalo <- theta*0.01
    idtr <- data1$x
    x <- as.matrix(x)
    z <- as.matrix(data$z[idtr,])
    nf <- ncol(z)
    y <- data1$y2
    TT <- data1$wx
    lf <- as.matrix(data1$wa)
    lik0 <- function(param) lik_cc_sep(x, z, y, lf, TT,  param)
    opts <- list( "algorithm" = "NLOPT_LN_NELDERMEAD", 
                  "xtol_rel" = 1.0e-3, 
                  "maxeval" = 5)
    param.opt <- nloptr(c(rep(0.5, ncol(z)), theta), lik0, 
                        eval_grad_f = NULL, 
                        lb = c(rep(0, ncol(z)), thetalo), 
                        ub =c(rep(0.999, ncol(z)), thetaup),
                        opts=opts)$solution
    thetaopt <- param.opt[-c(1:nf)]
    tauopt <- param.opt[1:nf]
    
    
    x <- x%*%diag(x=thetaopt, nrow=length(thetaopt), ncol=length(thetaopt))
    R <- kernelMatrix(kf, x)
    Taus <- list()
    Taus[[1]] <- matrix(tauopt[nf], max(z[,1]), max(z[,1])) 
    diag(Taus[[1]]) <- 1
    Tau <- Taus[[1]]
    for(i in 2:nf) {
      Taus[[i]] <- matrix(tauopt[nf-i+1], max(z[,i]), max(z[,i])) 
      diag(Taus[[i]]) <- 1
      Tau <- kronecker(Tau, Taus[[i]])
    }
    Sigma0 <- Tau
    Tau <- Tau[TT,TT]
    R <- R*Tau
    invR<-solve(R)
    l<-log(det(R))
    FinvR<-t(lf)%*%invR
    sigma2 <- t(y)%*%(invR-t(FinvR)%*%solve(FinvR%*%lf)%*%FinvR)%*%y/length(y)
    maxlikre <- NULL
    maxlikre$thetaopt <- thetaopt
    maxlikre$Sigma0 <- Sigma0*as.vector(sigma2)
    maxlikre
  }
  
  lik_cc_sep<-function(x, z, y, lf, TT, param) {
    nf <- ncol(z)
    theta <- param[-c(1:nf)]
    tau <- param[1:nf]
    x <- x%*%diag(x=theta, nrow=length(theta), ncol=length(theta))
    R <- kernelMatrix(kf, x)
    Taus <- list()
    Taus[[1]] <- matrix(tau[nf], max(z[,1]), max(z[,1])) 
    diag(Taus[[1]]) <- 1
    Tau <- Taus[[1]]
    for(i in 2:nf) {
      Taus[[i]] <- matrix(tau[nf-i+1], max(z[,i]), max(z[,i])) 
      diag(Taus[[i]]) <- 1
      Tau <- kronecker(Tau, Taus[[i]])
    }
    Tau <- Tau[TT,TT]
    R <- R*Tau
    invR<-solve(R)
    l<-log(det(R))
    FinvR<-t(lf)%*%invR
    l+log(t(y)%*%(invR-t(FinvR)%*%solve(FinvR%*%lf)%*%FinvR)%*%y/length(y))*length(y)
  }
  K <- 8
  n <- length(y)
  n.taus <- (K*K-K)/2
  taus <- param[1:n.taus]
  theta <- param[-c(1:n.taus)]
  taus.mat <- matrix(0, K, K)
  for(k in 2: K) {
    taus.mat[k, 1:(k-1)] <- taus[1:(k-1)]
    taus <- taus[-(1:(k-1))]
  }
  LL <- matrix(0, K, K)
  LL[1, 1] <- 1
  for(r in 2:nrow(LL)) {
    LL[r, 1] <- cos(taus.mat[r, 1])
    LL[r, r] <- prod(sin(taus.mat[r,1:(r-1)]))
    if(r>2) {
      for(s in 2:(r-1)) {
        LL[r, s] <- prod(sin(taus.mat[r,1:(s-1)]))*cos(taus.mat[r, s])
      }
    }
  }
  Tau <- LL %*% t(LL)
  
  
  
  param.opt <- nloptr(c(rep(0.5, ncol(z)), theta), lik0, 
                      eval_grad_f = NULL, 
                      lb = c(rep(0, ncol(z)), thetalo), 
                      ub =c(rep(0.999, ncol(z)), thetaup),
                      opts=opts)$solution
  
  
  maxlik_uc<-function(data, kf) {
    d1 <- ncol(data$x)
    #correlation parameter
    theta <- rep(1,d1)
    
    #bounds of the correlation parameter
    thetaup <- theta*100
    thetalo <- theta*0.01
    idtr <- data$idtr
    x <- as.matrix(data$x[idtr,])
    y <- data$y[idtr]
    TT <- data$TT[idtr]
    lf <- as.matrix(data$lfm[idtr,])
    K <- max(TT)
    n.taus <- (K*K-K)/2
    lik0 <- function(param) lik_uc(x,y,lf,TT, param)
    opts <- list( "algorithm" = "NLOPT_LN_NELDERMEAD", 
                  "xtol_rel" = 1.0e-3, 
                  "maxeval" = 5)
    param.opt <- nloptr(c(rep(pi/4, n.taus), data$thetaopt), lik0, 
                        eval_grad_f = NULL, 
                        lb = c(rep(0, n.taus), thetalo), 
                        ub =c(rep(pi, n.taus), thetaup),
                        opts=opts)$solution
    thetaopt <- param.opt[-c(1:n.taus)]
    taus <- param.opt[1:n.taus]    
    taus.mat <- matrix(0, K, K)
    for(k in 2: K) {
      taus.mat[k, 1:(k-1)] <- taus[1:(k-1)]
      taus <- taus[-(1:(k-1))]
    }
    LL <- matrix(0, K, K)
    LL[1, 1] <- 1
    for(r in 2:nrow(LL)) {
      LL[r, 1] <- cos(taus.mat[r, 1])
      LL[r, r] <- prod(sin(taus.mat[r,1:(r-1)]))
      if(r>2) {
        for(s in 2:(r-1)) {
          LL[r, s] <- prod(sin(taus.mat[r,1:(s-1)]))*cos(taus.mat[r, s])
        }
      }
    }
    Tauopt <- LL %*% t(LL)
    
    x <- x%*%diag(x=thetaopt, nrow=length(theta), ncol=length(theta))
    R <- kernelMatrix(kf, x)
    Tau <- Tauopt[TT,TT]
    R <- R*Tau
    invR<-solve(R+diag(rep(0.000001, nrow(R))))
    l<-log(det(R))
    FinvR<-t(lf)%*%invR
    sigma2 <- t(y)%*%(invR-t(FinvR)%*%solve(FinvR%*%lf)%*%FinvR)%*%y/length(y)
    maxlikre <- NULL
    maxlikre$thetaopt <- thetaopt
    maxlikre$Sigma0 <- Tauopt*as.vector(sigma2)
    maxlikre
  }
  
  lik_uc<-function(x, y, lf, TT,  param) {
    K <- max(TT)
    n <- length(y)
    n.taus <- (K*K-K)/2
    taus <- param[1:n.taus]
    theta <- param[-c(1:n.taus)]
    taus.mat <- matrix(0, K, K)
    for(k in 2: K) {
      taus.mat[k, 1:(k-1)] <- taus[1:(k-1)]
      taus <- taus[-(1:(k-1))]
    }
    LL <- matrix(0, K, K)
    LL[1, 1] <- 1
    for(r in 2:nrow(LL)) {
      LL[r, 1] <- cos(taus.mat[r, 1])
      LL[r, r] <- prod(sin(taus.mat[r,1:(r-1)]))
      if(r>2) {
        for(s in 2:(r-1)) {
          LL[r, s] <- prod(sin(taus.mat[r,1:(s-1)]))*cos(taus.mat[r, s])
        }
      }
    }
    Tau <- LL %*% t(LL)
    x <- x%*%diag(x=theta, nrow=length(theta), ncol=length(theta))
    R <- kernelMatrix(kf, x)
    Tau <- Tau[TT,TT]
    R <- R*Tau
    invR<-solve(R+diag(rep(0.000001, nrow(R))))
    l<-log(det(R))
    FinvR<-t(lf)%*%invR
    l+log(t(y)%*%(invR-t(FinvR)%*%solve(FinvR%*%lf, FinvR))%*%y/n)*n
  }