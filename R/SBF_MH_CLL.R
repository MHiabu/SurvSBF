
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
  alpha.smooth.i<-array(dim=c(d,n))
  
  alpha.smooth.i[1,] <-rep( 1,n)
  alpha.smooth.i.0<-(K.b/k.b)%*%(alpha[[1]]*dx[[1]])
  
  for (k in 2:d){
    for (i in 1:n)
    {
      alpha.smooth.i[k,i] <- (K.X.b[[k]][i,]/k.X.b[[k]][i])%*%(alpha[[k]]*dx[[k]])
    }}
  return(list(alpha.smooth.i=alpha.smooth.i,alpha.smooth.i.0=alpha.smooth.i.0))
}

taylor.alpha<-function(alpha,alpha.1,dx,K.X.b,k.X.b,K.b,k.b,dX.b,dX0.b,x.grid,n.grid,d,n)
{
  
  taylor.alpha.i<-taylor.alpha.i.square<-array(dim=c(d,n))
  taylor.alpha.i[1,] <- taylor.alpha.i.square[1,]<- rep(1,n)
  
  
  taylor.alpha.i.0<-((dX0.b*K.b/k.b)%*%(alpha.1[[1]]*dx[[1]])) + ((K.b/k.b)%*%(alpha[[1]]*dx[[1]]))
  
  
  temp <-taylor.alpha.i.0
  taylor.alpha.i.0<-t(alpha[[1]]-(dX0.b)*alpha.1[[1]])
  taylor.alpha.i.0<- as.numeric((taylor.alpha.i.0*K.b/k.b)%*%dx[[1]])
  
  # plot(alpha.1[[1]])
  # plot(taylor.alpha.i.0)
  # lines(temp,col='red')
  
  
  taylor.alpha.i.0.square <- t(alpha[[1]]-(dX0.b)*alpha.1[[1]])
  taylor.alpha.i.0.square <- as.numeric((taylor.alpha.i.0.square^2*K.b/k.b)%*%dx[[1]])
  
  for (k in 2:d){
    for (i in 1:n)
    {
      taylor.alpha.i[k,i] <- ( (K.X.b[[k]][i,]/k.X.b[[k]][i]) * (alpha[[k]]+alpha.1[[k]]*dX.b[[k]][i,]) )    %*% dx[[k]]  
    }
  }
  
  
  for (k in 2:d){
    for (i in 1:n)
    {
      taylor.alpha.i.square[k,i] <-(  (K.X.b[[k]][i,]/k.X.b[[k]][i])* (alpha[[k]]+alpha.1[[k]]*dX.b[[k]][i,])^2 ) %*% dx[[k]]  
    }
  }
  
  
  return(list( taylor.alpha.i= taylor.alpha.i,taylor.alpha.i.square= taylor.alpha.i.square,taylor.alpha.i.0=taylor.alpha.i.0,taylor.alpha.i.0.square=taylor.alpha.i.0.square             ))
}

get.alpha.new<-function(time,status,alpha,alpha.1,dx,K.X.b,k.X.b,K.b,k.b,dX.b,dX0.b,Y,k,x.grid,n.grid,d,n,l)
{
  
  if (weight=='sw'){
    taylor.alpha<-taylor.alpha(alpha,alpha.1,dx,K.X.b,k.X.b,K.b,k.b,dX.b,dX0.b,x.grid,n.grid,d,n)
    
    if (k==1)  taylor.alpha.minusk<-numeric(n) else   taylor.alpha.minusk<-array(dim=c(n,n.grid[1])) 
    for (i in 1:n)
    {
      if (k==1) taylor.alpha.minusk[i] <-  prod(taylor.alpha$taylor.alpha.i[-1,i]) else
      {taylor.alpha.minusk[i,] <-  prod(taylor.alpha$taylor.alpha.i[-k,i])*taylor.alpha$taylor.alpha.i.0
      }
    } 
    
    if (k==1)
    {
      D<- rowSums(sapply(1:n, function (i) { return( (K.b/k.b) %*%  (dx[[1]]*(taylor.alpha.minusk[i]*Y[i,]))  ) }))
    }else D<- rowSums(sapply(1:n, function (i) { return(as.numeric(dx[[1]]%*%(Y[i,]*(taylor.alpha.minusk[i,])))*(K.X.b[[k]][i,]/k.X.b[[k]][i]))}))
    
    
    if (k==1)
    {
      O2<- rowSums(sapply(1:n, function (i) { return(   (t(dX0.b)*K.b/k.b)  %*% (dx[[1]]*(taylor.alpha.minusk[i]*Y[i,]))   )}))
    }else O2<- rowSums(sapply(1:n, function (i) { return(as.numeric(dx[[1]]%*%(Y[i,]*(taylor.alpha.minusk[i,])))*(dX.b[[k]][i,]*K.X.b[[k]][i,]/k.X.b[[k]][i]))}))
    
    
    
    #  O<- rowSums(sapply(1:n, function (i) { return((alpha.minusk.smooth[kk,bb,i]*data$status[i])*(K.X.b[[kk]][bb,i,]/k.X.b[[kk]][bb,i]))})) ### for different waiting
    O1<- rowSums(sapply(1:n, function (i) { return((status[i])*(K.X.b[[k]][i,]/k.X.b[[k]][i]))}))
    O<-O1-alpha.1[[k]]*O2
  }
  
  
  
  if (weight=='no'){
    taylor.alpha<-taylor.alpha(alpha,alpha.1,dx,K.X.b,k.X.b,K.b,k.b,dX.b,dX0.b,x.grid,n.grid,d,n)
    # plot(taylor.alpha$taylor.alpha.i.0,ylab='nosquare')
    # plot(taylor.alpha$taylor.alpha.i.0.square,ylab='square')
    
    
    if (k==1)  taylor.alpha.minusk<-numeric(n) else   taylor.alpha.minusk<-array(dim=c(n,n.grid[1])) 
    for (i in 1:n)
    {
      if (k==1) taylor.alpha.minusk[i] <-  prod(taylor.alpha$taylor.alpha.i[-1,i]) else
      {taylor.alpha.minusk[i,] <-  prod(taylor.alpha$taylor.alpha.i[-k,i])*taylor.alpha$taylor.alpha.i.0
      }
    } 
    
    
    
    if (k==1)  taylor.alpha.minusk.square<-numeric(n) else   taylor.alpha.minusk.square<-array(dim=c(n,n.grid[1])) 
    for (i in 1:n)
    {
      if (k==1) taylor.alpha.minusk.square[i] <-  prod(taylor.alpha$taylor.alpha.i.square[-1,i]) else
      {taylor.alpha.minusk.square[i,] <-  prod(taylor.alpha$taylor.alpha.i.square[-k,i])*taylor.alpha$taylor.alpha.i.0.square
      }
    }
    
    
    
    if (k==1)
    {
      D<- rowSums(sapply(1:n, function (i) { return(  (K.b/k.b)  %*%  ((dx[[1]]*taylor.alpha.minusk.square[i]*Y[i,]))   )}))
    }else D<- rowSums(sapply(1:n, function (i) { return(as.numeric(dx[[1]]%*%(Y[i,]*(taylor.alpha.minusk.square[i,])))*(K.X.b[[k]][i,]/k.X.b[[k]][i]))}))
    
    
    
    if (k==1)
    {
      O2<- rowSums(sapply(1:n, function (i) { return((t(dX0.b)*K.b/k.b)  %*% (dx[[1]]*(taylor.alpha.minusk.square[i]*Y[i,]))  )}))
    }else O2<- rowSums(sapply(1:n, function (i) { return(as.numeric(dx[[1]]%*%(Y[i,]*(taylor.alpha.minusk.square[i,])))*(dX.b[[k]][i,]*K.X.b[[k]][i,]/k.X.b[[k]][i]))}))
    
    
    
    if (k==1)
    {
      O1<- rowSums(sapply(1:n, function (i) { return((status[i])*taylor.alpha.minusk[i]*(K.X.b[[k]][i,]/k.X.b[[k]][i]))}))
    }else O1<- rowSums(sapply(1:n, function (i) { return((status[i])*taylor.alpha.minusk[i,which(x.grid[[1]]==time[i])]*(K.X.b[[k]][i,]/k.X.b[[k]][i]))}))
    
    
    O<-O1-alpha.1[[k]]*O2
  }
  
  
  
  
  return(O/D)
} 

get.alpha.1.new<-function(time,status,alpha,alpha.1,dx,K.X.b,k.X.b,K.b,k.b,dX.b,dX0.b,Y,k,x.grid,n.grid,d,n)
{
  
  if (weight=='sw'){
    taylor.alpha<-taylor.alpha(alpha,alpha.1,dx,K.X.b,k.X.b,K.b,k.b,dX.b,dX0.b,x.grid,n.grid,d,n)
    if (k==1)  taylor.alpha.minusk<-numeric(n) else   taylor.alpha.minusk<-array(dim=c(n,n.grid[1])) 
    for (i in 1:n)
    {
      if (k==1) taylor.alpha.minusk[i] <-  prod(taylor.alpha$taylor.alpha.i[-1,i]) else
      {taylor.alpha.minusk[i,] <-  prod(taylor.alpha$taylor.alpha.i[-k,i])*taylor.alpha$taylor.alpha.i.0
      }
    } 
    
    
    if (k==1)
    {
      D<- rowSums(sapply(1:n, function (i) { return(((t(dX0.b))^2*K.b/k.b)%*%((dx[[1]]*taylor.alpha.minusk[i]*Y[i,])))}))
    }else D<- rowSums(sapply(1:n, function (i) { return(as.numeric(dx[[1]]%*%(Y[i,]*(taylor.alpha.minusk[i,])))*((dX.b[[k]][i,])^2*K.X.b[[k]][i,]/k.X.b[[k]][i]))}))
    
    
    if (k==1)
    {
      O2<- rowSums(sapply(1:n, function (i) { return((t(dX0.b)*K.b/k.b)%*%(dx[[1]]*(taylor.alpha.minusk[i]*Y[i,])))}))
    }else O2<- rowSums(sapply(1:n, function (i) { return(as.numeric(dx[[1]]%*%(Y[i,]*(taylor.alpha.minusk[i,])))*(dX.b[[k]][i,]*K.X.b[[k]][i,]/k.X.b[[k]][i]))}))
    
    
    #  O<- rowSums(sapply(1:n, function (i) { return((alpha.minusk.smooth[kk,bb,i]*data$status[i])*(K.X.b[[kk]][bb,i,]/k.X.b[[kk]][bb,i]))})) ### for different waiting
    O1<- rowSums(sapply(1:n, function (i) { return((status[i])*(dX.b[[k]][i,]*K.X.b[[k]][i,]/k.X.b[[k]][i]))}))
    O<-O1-alpha[[k]]*O2
  }
  
  
  
  if (weight=='no'){
    taylor.alpha<-taylor.alpha(alpha,alpha.1,dx,K.X.b,k.X.b,K.b,k.b,dX.b,dX0.b,x.grid,n.grid,d,n)
    
    if (k==1)  taylor.alpha.minusk<-numeric(n) else   taylor.alpha.minusk<-array(dim=c(n,n.grid[1])) 
    for (i in 1:n)
    {
      if (k==1) taylor.alpha.minusk[i] <-  prod(taylor.alpha$taylor.alpha.i[-1,i]) else
      {taylor.alpha.minusk[i,] <-  prod(taylor.alpha$taylor.alpha.i[-k,i])*taylor.alpha$taylor.alpha.i.0
      }
    } 
    
    if (k==1)  taylor.alpha.minusk.square<-numeric(n) else   taylor.alpha.minusk.square<-array(dim=c(n,n.grid[1])) 
    for (i in 1:n)
    {
      if (k==1) taylor.alpha.minusk.square[i] <-  prod(taylor.alpha$taylor.alpha.i.square[-1,i]) else
      {taylor.alpha.minusk.square[i,] <-  prod(taylor.alpha$taylor.alpha.i.square[-k,i])*taylor.alpha$taylor.alpha.i.0.square
      }
    }
    
    
    if (k==1) { D<- rowSums(sapply(1:n, function (i) { return(  ((t(dX0.b))^2*K.b/k.b) %*% ((dx[[1]]*taylor.alpha.minusk.square[i]*Y[i,])) ) }))
    }else D<- rowSums(sapply(1:n, function (i) { return(as.numeric(dx[[1]]%*%(Y[i,]*(taylor.alpha.minusk.square[i,])))*((dX.b[[k]][i,])^2*K.X.b[[k]][i,]/k.X.b[[k]][i]))}))
    
    
    if (k==1)
    {
      O2<- rowSums(sapply(1:n, function (i) { return((t(dX0.b)*(K.b/k.b))%*%(dx[[1]]*(taylor.alpha.minusk.square[i]*Y[i,])))}))
    }else O2<- rowSums(sapply(1:n, function (i) { return(as.numeric(dx[[1]]%*%(Y[i,]*(taylor.alpha.minusk.square[i,])))*((dX.b[[k]][i,])*K.X.b[[k]][i,]/k.X.b[[k]][i]))}))
    
    
    #  O<- rowSums(sapply(1:n, function (i) { return((alpha.minusk.smooth[kk,bb,i]*data$status[i])*(K.X.b[[kk]][bb,i,]/k.X.b[[kk]][bb,i]))})) ### for different waiting
    if (k==1)
    {
      O1<- rowSums(sapply(1:n, function (i) { return((status[i])*taylor.alpha.minusk[i]*((dX.b[[k]][i,])*K.X.b[[k]][i,]/k.X.b[[k]][i]))}))
    }else O1<- rowSums(sapply(1:n, function (i) { return((status[i])*taylor.alpha.minusk[i,which(x.grid[[1]]==time[i])]*((dX.b[[k]][i,])*K.X.b[[k]][i,]/k.X.b[[k]][i]))}))
    
    O<-O1-alpha[[k]]*O2
    
  }
  
  
  
  
  
  return(O/D)
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


u<-x.grid.array-t(x.grid.array) # u is t-s for all t,s on the grid considered
K.b<-apply(u/bandwidth[1],1:2,kern)/(bandwidth[1])
k.b<-colSums(dx[[1]]*K.b)                      # small k for normailzing kernel, sincs grid points are symmetric normalization can be used for both row and column
if (kcorr==FALSE) {k.b<-rep(1,n.grid[1])}                     


dX0.b<-u#/bandwidth[1]

X<-t(X)  # status is in the last column of data


Y<-t(sapply(1:n,function(i) { temp<-numeric(n.grid[1])
for (l in 1:n.grid[1]) {temp[l]<-as.numeric((x.grid[[1]][l]<=time[i]))
} 
return(temp) 
}
))

K.X.b<-k.X.b<-list()
for( k in 1:d){
  K.X.b[[k]]<-array(0,dim=c(n,n.grid[k]))
  k.X.b[[k]]<-numeric(n)
  for (i in 1:n)
  {
    u<-x.grid[[k]]
    u<-(X[k,i]-u)
    K.X.b[[k]][i,]<-sapply(u/bandwidth[k],kern)/(bandwidth[k])
    k.X.b[[k]][i]<-sum(dx[[k]]*K.X.b[[k]][i,])
    if (kcorr==FALSE) {k.X.b[[k]][i]<-1}
    
  } }


dX.b<-list()
for( k in 1:d){
  dX.b[[k]]<-array(0,dim=c(n,n.grid[k]))
  u<-x.grid[[k]]
  u<-X[k,]-matrix(u,nrow=n,ncol=length(u),byrow=TRUE)
  dX.b[[k]]<- (u)#/bandwidth[k])
}







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
    
    alpha_backfit[[k]]<- get.alpha.new(time,status,alpha_backfit,alpha.1_backfit,dx,K.X.b,k.X.b,K.b,k.b,dX.b,dX0.b,Y,k,x.grid,n.grid,d,n,l)
    
    
    #for (j in 2:d) {  alpha_backfit[[j]]<- alpha_backfit[[j]]/ alpha_backfit[[j]][1]
    # alpha.1_backfit[[k]]<- alpha.1_backfit[[k]]/ alpha_backfit[[k]][1]
    # }
    
    alpha_backfit[[k]][1]<- alpha_backfit[[k]][2]
    alpha_backfit[[k]][length(x.grid[[k]])]<- alpha_backfit[[k]][length(x.grid[[k]])-1]
    
    
    alpha_backfit[[k]][is.nan(alpha_backfit[[k]])]<-1
    alpha_backfit[[k]][is.na(alpha_backfit[[k]])]<-1
    alpha_backfit[[k]][alpha_backfit[[k]]==Inf]<-1
    alpha_backfit[[k]][alpha_backfit[[k]]<0]<-0.0001
    #  }
    #  for(k in 1:d)
    #   {
    #    if (k!=1) {
    
    #     }
    if (LC==TRUE) alpha.1_backfit[[k]]<-rep(0,n.grid[k]) else{
      alpha.1_backfit[[k]]<- get.alpha.1.new(time,status,alpha_backfit,alpha.1_backfit,dx,K.X.b,k.X.b,K.b,k.b,dX.b,dX0.b,Y,k,x.grid,n.grid,d,n)}
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
