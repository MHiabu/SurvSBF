
SBF.MH.CLL<-function(formula,data,bandwidth,weight='sw',x.grid=NULL,n.grid.additional=0, x.min=NULL, x.max=NULL, integral.approx='right',it=100,kern=function(u){return(0.75*(1-u^2                         )*(abs(u)<1))},initial=NULL,kcorr=kcorr,LC)
{   data<-data[,c(1,3:ncol(data),2)]

# 
# formula <- formula(formula)
# if (class(formula) != "formula") {
#   stop("Error: Invalid formula.")
# }
# data.selected <- as.list(attr(terms(frmla), "variables"))[-1]
#   

Terms <- terms(x=formula,data=data)
mm <- na.omit(get_all_vars(formula(Terms),data=data))
if (NROW(mm) == 0) stop("No (non-missing) observations")

response <- model.response(model.frame(update(formula,".~1"),data=mm))
X       <- prodlim::model.design(Terms,
                                 data=mm,
                                 maxOrder=1,
                                 dropIntercept=TRUE)[[1]]


time <- as.vector(response[,"time"])
status <- as.vector(response[,"status"])
X <- cbind( time,X)

# }}}

smooth.alpha<-function(alpha,K.X.b,k.X.b,K.b,k.b,d,n)
{
  alpha.smooth.i<-array(dim=c(n,d))
  
  alpha.smooth.i[,1] <-rep( 1,n)
  
  alpha.smooth.i.0<- (alpha[[1]]*dx[[1]]) %*% K.b
  
  
  for (k in 2:d){
      alpha.smooth.i[,k] <- (alpha[[k]]*dx[[k]]) %*% K.X.b[[k]]
  }
  
  
  
  
  
  return(list(alpha.smooth.i=alpha.smooth.i,alpha.smooth.i.0=alpha.smooth.i.0))
}

taylor.alpha<-function(alpha,alpha.1,dx,K.X.b,K.b,dX.b,dX0.b,x.grid,n.grid,d,n)
{
  
  taylor.alpha.i <- array(dim=c(n,d))
  
  taylor.alpha.i[ ,1] <- rep(1,n)
  
  
  #taylor.alpha.i.0<- ((dX0.b*K.b/k.b)%*%(alpha.1[[1]]*dx[[1]])) + ((K.b/k.b)%*%(alpha[[1]]*dx[[1]]))
 # temp <-taylor.alpha.i.0
  
  taylor.alpha.i.0 <-   as.numeric( t(dx[[1]])   %*%   (( alpha[[1]]+(dX0.b*alpha.1[[1]] ))*K.b) )
  
  
   # plot(alpha.1[[1]])
   # plot(taylor.alpha.i.0)
   # lines(temp,col='red')
   # 
  # 
  # taylor.alpha.i.0.square <- t(alpha[[1]]-(dX0.b)*alpha.1[[1]])
  # taylor.alpha.i.0.square <- as.numeric((taylor.alpha.i.0.square^2*K.b/k.b)%*%dx[[1]])
  
  for (k in 2:d){
    
  taylor.alpha.i[,k] <-   dx[[k]] %*% ((alpha[[k]]+alpha.1[[k]]*dX.b[[k]])* K.X.b[[k]])  
    
  }
  
  
  # for (k in 2:d){
  #   for (i in 1:n)
  #   {
  #     taylor.alpha.i.square[k,i] <-(  (K.X.b[[k]][i,]/k.X.b[[k]][i])* (alpha[[k]]+alpha.1[[k]]*dX.b[[k]][i,])^2 ) %*% dx[[k]]  
  #   }
  # }
  
  
  return(list( taylor.alpha.i= taylor.alpha.i,taylor.alpha.i.0=taylor.alpha.i.0))
}

get.alpha.new<-function(O1,time,status,alpha,alpha.1,dx,K.X.b,K.b,dX.b,dX0.b,Y,k,x.grid,n.grid,d,n,l)
{
  
    taylor.alpha <- taylor.alpha(alpha,alpha.1,dx,K.X.b,K.b,dX.b,dX0.b,x.grid,n.grid,d,n)
    
    if (k==1)  taylor.alpha.minusk<-numeric(n) else   taylor.alpha.minusk<-array(dim=c(n,n.grid[1])) 
   
   
    if (k==1) taylor.alpha.minusk <-  apply(taylor.alpha$taylor.alpha.i[,-1],1,prod) else  
        {taylor.alpha.minusk <- apply(taylor.alpha$taylor.alpha.i[,-k],1,prod) * matrix(taylor.alpha$taylor.alpha.i.0, ncol=n, nrow=n.grid[[1]], byrow = TRUE)
    }
  
    if (k==1)
    {
      D <-     colSums( (taylor.alpha.minusk*Y)    %*%  (t(K.b) *dx[[1]])  ) 
    }else
      D <-   K.X.b[[k]]  %*%  (taylor.alpha.minusk*Y)    %*%  dx[[1]] 
    
      
    
    if (k==1)
    { 
      O2 <- colSums( (taylor.alpha.minusk*Y)    %*%  (t(dX0.b*K.b) *dx[[1]])  )  
    }else 
      O2 <- ( dX.b[[k]]*K.X.b[[k]] ) %*%  (taylor.alpha.minusk*Y)    %*%  dx[[1]] 

    
    
    O<-O1[[k]] - alpha.1[[k]]*O2
  
  

  
  return(as.numeric(O/D))
} 

get.alpha.1.new<-function(O1.1,time,status,alpha,alpha.1,dx,K.X.b,K.b,dX.b,dX0.b,Y,k,x.grid,n.grid,d,n)
{
  
    taylor.alpha <- taylor.alpha(alpha,alpha.1,dx,K.X.b,K.b,dX.b,dX0.b,x.grid,n.grid,d,n)
    
    
    if (k==1)  taylor.alpha.minusk<-numeric(n) else   taylor.alpha.minusk<-array(dim=c(n,n.grid[1])) 
    
    if (k==1) taylor.alpha.minusk <-  apply(taylor.alpha$taylor.alpha.i[,-1],1,prod) else  
    {taylor.alpha.minusk <- apply(taylor.alpha$taylor.alpha.i[,-k],1,prod) * matrix(taylor.alpha$taylor.alpha.i.0, ncol=n, nrow=n.grid[[1]], byrow = TRUE)
    }
    
    
    
    
    if (k==1)
    {
      D.1 <-     colSums( (taylor.alpha.minusk*Y)    %*%  (t((dX0.b)^2*K.b) *dx[[1]])  ) 
    }else
      D.1 <-   ((dX.b[[k]])^2*K.X.b[[k]])  %*%  (taylor.alpha.minusk*Y)    %*%  dx[[1]] 
    
    
    
    
    
    
    if (k==1)
    { 
      O2.1 <- colSums( (taylor.alpha.minusk*Y)    %*%  (t(dX0.b*K.b) *dx[[1]])  )  
    }else 
      O2.1 <- ( dX.b[[k]]*K.X.b[[k]] ) %*%  (taylor.alpha.minusk*Y)    %*%  dx[[1]] 
    

   
    O.1<-O1.1[[k]] - alpha[[k]]*O2.1
    
    

  
  
  
  
  
  
  
  
  
  return(as.numeric(O.1/D.1))
} 

get.O1<-function(k,status,K.X.b){
  
  O1 <- as.numeric(  K.X.b[[k]] %*% status )
  
  return(O1)
}


get.O1.1<-function(k,status,K.X.b,dX.b){
  O1<- as.numeric((dX.b[[k]]*K.X.b[[k]]) %*% status)
  return(O1)
}






d <- ncol(X) 
n <- nrow(X)

if(is.null(x.grid)) x.grid<-lapply(1:d,function(k) X[order(X[,k]),k])
if(is.null(x.min)) x.min<-sapply(x.grid,head,1)
if(is.null(x.max)) x.max<-sapply(x.grid,tail,1)



x.grid.additional <- lapply(1:d, function(k) seq(x.min[k],x.max[k], length=n.grid.additional))
x.grid <- lapply(1:d, function(k) sort(c(x.grid[[k]], x.grid.additional[[k]])))
n.grid <- sapply(x.grid, length)

if (integral.approx=='midd'){
  dx<-lapply(1:d,function(k){  ddx<-diff(x.grid[[k]])
  c(  (ddx[1]/2)+(x.grid[[k]][1]-x.min[k])  ,(ddx[-(n.grid[k]-1)]+ddx[-1])/2,  
      (ddx[n.grid[k]-1]/2)+(x.max[k]-x.grid[[k]][n.grid[k]])  
  )
  }
  )
}

if (integral.approx=='left'){
  dx<-lapply(1:d,function(k){  ddx<-diff(x.grid[[k]])
  c(  ddx,  x.max[k]-x.grid[[k]][n.grid[k]])  
  }
  )
}

if (integral.approx=='right'){
  dx<-lapply(1:d,function(k){  ddx<-diff(x.grid[[k]])
  c(  (x.grid[[k]][1]-x.min[k])  ,ddx)
  }
  )
}




K.b<-array(0,dim=c(n.grid[1],n.grid[1]))
k.b<-array(0,dim=c(n.grid[1]))



x.grid.array<-matrix(rep(x.grid[[1]], times=n.grid[1]),nrow = n.grid[1], ncol = n.grid[1],byrow=FALSE)   # 


u<- x.grid.array-t(x.grid.array) # u[t,s]= t-s for all t,s on the grid considered
K.b<-apply(u/bandwidth[1],1:2,kern)/(bandwidth[1])
k.b<-colSums(dx[[1]]*K.b)                      # small k for normailzing kernel, since grid points are symmetric normalization can be used for both row and column
if (kcorr==FALSE) {k.b<-rep(1,n.grid[1])}                     


dX0.b<- u#/bandwidth[1]     # dX0.b[t,s]= t-s for all t,s on the grid considered
K.b <- K.b %*% diag(1/k.b)  ### row-wise division --> colSums(dx[[1]]*K.b)=1


### define exposure=Y[i,s]
Y<-t(sapply(1:n,function(i) { temp<-numeric(n.grid[1])
for (l in 1:n.grid[1]) {temp[l]<-as.numeric((x.grid[[1]][l]<=time[i]))
} 
return(temp) 
}
))

dX.b<-K.X.b<-k.X.b<-list()
for( k in 1:d){
  K.X.b[[k]]<-array(0,dim=c(n,n.grid[k]))
  k.X.b[[k]]<-numeric(n)
  
  x.grid.array<-matrix(-X[,k],nrow =n.grid[k], ncol =n ,byrow=TRUE)   # 
  u<- x.grid.array+x.grid[[k]]  ####  u[,i]=x_k - X_{ik}
  K.X.b[[k]]<-apply(u/bandwidth[k],1:2,kern)/(bandwidth[k])
  k.X.b[[k]]<-colSums(dx[[k]]*K.X.b[[k]])   
  if (kcorr==FALSE) k.X.b[[k]] <- rep(1,n)
  dX.b[[k]]<- u
  K.X.b[[k]]<- K.X.b[[k]] %*%  diag(1/k.X.b[[k]])  ### row-wise division ---> col integration=1
}


O1<-O1.1<-list()
O1 <- lapply(1:d, get.O1, status,K.X.b)
O1.1 <- lapply(1:d, get.O1.1,status,K.X.b,dX.b)


alpha_backfit<-list()
alpha.1_backfit<-list()

if (is.null(initial)){
  for(k in 1:d){
    alpha_backfit[[k]]<-rep(1, n.grid[k])
  }
} else  alpha_backfit<-initial

for(k in 1:d){
  alpha.1_backfit[[k]]<-rep(0, n.grid[k])
  
}


count<-rep(1,d)
for (l in 2:it)
{  
  
  
  
  alpha_backfit_old<-alpha_backfit
  alpha.1_backfit_old<-alpha.1_backfit
  for(k in 1:d)
  {
    # if(k==1) alpha.temp<-lapply(1:d,function(z){alpha_backfit[[z]][b,]}) else 
    #   { alpha.temp[1:(k-1)]<-lapply(1:(k-1),function(z){alpha_backfit[[z]][b,]})
    #     #alpha.temp[k:d]<-lapply(k:d,function(z){alpha_backfit_old[[z]][b,]})
    #     }
    
    
    
    
    #for(inner in 1:5){
    
    alpha_backfit[[k]]<- get.alpha.new(O1,time,status,alpha_backfit,alpha.1_backfit,dx,K.X.b,K.b,dX.b,dX0.b,Y,k,x.grid,n.grid,d,n,l)
    
    
    #for (j in 2:d) {  alpha_backfit[[j]]<- alpha_backfit[[j]]/ alpha_backfit[[j]][1]
    # alpha.1_backfit[[k]]<- alpha.1_backfit[[k]]/ alpha_backfit[[k]][1]
    # }
    
    alpha_backfit[[k]][1]<- alpha_backfit[[k]][2]
    alpha_backfit[[k]][length(x.grid[[k]])]<- alpha_backfit[[k]][length(x.grid[[k]])-1]
    
    
    alpha_backfit[[k]][is.nan(alpha_backfit[[k]])]<-1
    alpha_backfit[[k]][is.na(alpha_backfit[[k]])]<-1
    alpha_backfit[[k]][alpha_backfit[[k]]==Inf]<-1
    alpha_backfit[[k]][alpha_backfit[[k]]<=0]<-0.0001
    #  }
    #  for(k in 1:d)
    #   {
    #    if (k!=1) {
    
    #     }
    if (LC==TRUE) alpha.1_backfit[[k]]<-rep(0,n.grid[k]) else{
      alpha.1_backfit[[k]]<- get.alpha.1.new(O1.1,time,status,alpha_backfit,alpha.1_backfit,dx,K.X.b,K.b,dX.b,dX0.b,Y,k,x.grid,n.grid,d,n)}
    alpha.1_backfit[[k]][is.nan(alpha.1_backfit[[k]])]<-0
    alpha.1_backfit[[k]][is.na(alpha.1_backfit[[k]])]<-0
    # alpha.1_backfit[[k]][alpha.1_backfit[[k]]>30]<-30
    #  alpha.1_backfit[[k]][alpha.1_backfit[[k]]<=(-30)]<--30
    #   print(!(prod(alpha.1_backfit_old[[k]] == 0) ))
    
    
    if ((prod(alpha.1_backfit_old[[k]] == 0) ))  {### if old vaues all zero
      if (max(abs(alpha.1_backfit_old[[k]]- alpha.1_backfit[[k]])) >(100))    {alpha.1_backfit[[k]] <-rep(0, n.grid[k])
      #  alpha.1_backfit[[k]]<-rep(0, n.grid[k])
      count[k]<-1
      }}
    
    if (!(prod(alpha.1_backfit_old[[k]] == 0) )){ ### if old vaues NOT all zero
      if (max(abs(alpha.1_backfit_old[[k]]- alpha.1_backfit[[k]])) >(100/log(count[k])))    alpha.1_backfit[[k]] <-rep(0, n.grid[k])
      #  alpha.1_backfit[[k]]<-rep(0, n.grid[k])
    }
    
    count[k]<-count[k]+1
    # plot(x.grid[[k]],alpha.1_backfit[[k]],lty=3,col=1,lwd=2) 
    
  }
  #}
  
  
  
  #plot(x.grid[[2]],phi[[1]](x.grid[[2]]),lty=3,col=1,lwd=2) 
  #  if (l==2) plot(x.grid[[2]],log(alpha_backfit[[2]]),ylim=c(-2,2) ) else lines(x.grid[[2]],log(alpha_backfit[[2]] ) )
  #plot(x.grid[[1]], log(alpha_backfit[[1]]),col='blue',lwd=2) 
  
  
  if (max(max(abs(unlist(alpha_backfit_old)-unlist(alpha_backfit)),na.rm=TRUE),max(abs(unlist(alpha.1_backfit_old)-unlist(alpha.1_backfit)),na.rm=TRUE))<= 0.001) break
  for(k in 1:d){
    print(c(l,max(abs(unlist(log(alpha_backfit_old[[k]]))-unlist(log(alpha_backfit[[k]]))),na.rm=TRUE),max(abs(unlist((alpha.1_backfit_old[[k]]))-unlist((alpha.1_backfit[[k]]
    ))),na.rm=TRUE)))
  }
}

for(k in 2:d){
  if (alpha_backfit[[k]][1]!=0){
    alpha_backfit[[1]]<-alpha_backfit[[1]]*alpha_backfit[[k]][1]
    alpha_backfit[[k]]<-alpha_backfit[[k]]/ alpha_backfit[[k]][1]
    alpha.1_backfit[[k]]<-alpha.1_backfit[[k]]/ alpha_backfit[[k]][1]}
  
  # if (max(abs((alpha_backfit_old[[k]])-(alpha_backfit[[k]])),na.rm=TRUE) >5) alpha_backfit[[k]]<-alpha_backfit_old[[k]]
  #if (max(abs((alpha.1_backfit_old[[k]])-(alpha.1_backfit[[k]])),na.rm=TRUE) >5) alpha.1_backfit[[k]]<-alpha.1_backfit_old[[k]]
}
return(list(alpha_backfit=alpha_backfit,alpha.1_backfit=alpha.1_backfit,l=l,x.grid=x.grid))
}
