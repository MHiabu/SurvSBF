


mhcovariate <- function(n,d=5,rho=0,seed){
    if (!missing(seed)) set.seed(seed)
    stddev<-rep(1,d-1)
    corMat<-matrix(rho, nrow=d-1,ncol=d-1)
    corMat[col(corMat)==row(corMat)]<-1
    covMat<-stddev %*% t(stddev) * corMat
    Z<-mvrnorm(n=n,  mu=rep(0,d-1), Sigma = covMat, empirical = FALSE)
    Z<-2.5*atan(Z)/pi
    Z
}



# data generating mechanism
mhrate <- function(Z,model=1,violate.cox=TRUE){
     Z<-as.matrix(Z)
    d <- NCOL(Z)
    phi <- vector(length=d,mode="list")
    if (violate.cox==TRUE){
        # Cox violated
        for(k in 1:d){
            if ((k%%2)==1)  {phi[[k]]<- function(z) 1*sin(pi*z)} else {phi[[k]]<-function(z) -1*sin(pi*z)}
        }
    }else{
        # Cox satisfied
        for(k in 1:d){
            if ((k%%2)==1)  {phi[[k]]<- function(z) -2*z} else {phi[[k]]<-function(z) 2*z}
        }
    }
    
    
    
    
    top <- rep(0,NROW(Z))
    for (k in 1:d)
    {
        top<-  top+phi[[k]](Z[,k])
    }
    
  
    
    if (model==1)  surv.function <- function(t){ pexp(t,rate=exp(top),lower.tail = FALSE)}
    if (model==2)  surv.function<-function(t){ pmakeham(t,scale=1, shape=exp(top),lower.tail =                                                  FALSE)}
    
    #### true rate parameter for the hazard function
    # true_par<-0
    # for (k in 1:d){
    #     true_par<-  true_par+phi[[k]](Z[[k]])
    # }
    # true_par<-exp(true_par)
     return(list(top=top, surv.function=surv.function))
}

mhdata <- function(n=200,d=5,rho=0,model=1,violate.cox=TRUE,seed){
    if (!missing(seed)) set.seed(seed)
    Z <- mhcovariate(n=n,d=d,rho=rho)
    # regression coefficients
    top <- mhrate(Z,model,violate.cox=violate.cox)$top
    if (model==1){
    Time<-rexp(n,exp(top))
    C<-rexp(n,exp(top)/1.75)
 
    }

  if (model==2)
  #survivaldistr=='makeham' 
  { 
     beta<-1
     alpha<- 1
     alpha <- exp(top)*alpha
    Time<-rmakeham(n,beta, alpha)
    C<-rmakeham(n,beta, alpha/1.75)
  }
    
  TT<-Time*(Time<=C)+C*(Time>C)
  status<-(Time<=C)*1  ## censoring indicator
  data<-data.frame(time=TT,status=status,as.data.frame(Z))
  data

}



######################################################################
### mhdata.R ends here
