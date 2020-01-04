SBF.MH.LC<-function(formula,data,bandwidth,x.grid=NULL,n.grid.additional=0, x.min=NULL, x.max=NULL, integral.approx='midd',it=100,kern=function(u){return(0.75*(1-u^2)*(abs(u)<1))},initial=NULL)
{       
  
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
  
  
  smooth.alpha<-function(alpha,K.X.b,k.X.b,K.b,k.b,x.grid,dx,n.grid,d,n)
  {
    alpha.smooth.i<-array(dim=c(d,n))
    
    alpha.smooth.i[1,] <-rep( 1,n)
    alpha.smooth.i.0<-(K.b/k.b)%*%(alpha[[1]]*dx[[1]])
    
    for (k in 2:d){
      for (i in 1:n)
      {
        alpha.smooth.i[k,i] <- as.numeric((K.X.b[[k]][i,]/k.X.b[[k]][i])%*%(alpha[[k]]*dx[[k]]))
      }}
    return(list(alpha.smooth.i=alpha.smooth.i,alpha.smooth.i.0=alpha.smooth.i.0))
  }
  get.new.alpha<-function(status,alpha,K.X.b,k.X.b,K.b,k.b,Y,k,x.grid,dx,n.grid,d,n)
  {
    
    alpha.smooth.i<-smooth.alpha(alpha,K.X.b,k.X.b,K.b,k.b,x.grid,dx,n.grid,d,n)
    
    
    if (k==1) alpha.minusk.smooth<-numeric(n) else  alpha.minusk.smooth<-array(dim=c(n,n.grid[1])) 
    
    for (i in 1:n)
    {
      if (k==1) alpha.minusk.smooth[i] <-  prod(alpha.smooth.i$alpha.smooth.i[-1,i]) else
      {alpha.minusk.smooth[i,] <-  prod(alpha.smooth.i$alpha.smooth.i[-k,i])*alpha.smooth.i$alpha.smooth.i.0
      }
    } 
    
    if (k==1)
    {
      D<- rowSums(sapply(1:n, function (i) { return((dx[[1]]*(alpha.minusk.smooth[i]*Y[i,]))%*%(K.b/k.b))}))
    }else D<- rowSums(sapply(1:n, function (i) { return(as.numeric(dx[[1]]%*%(Y[i,]*(alpha.minusk.smooth[i,])))*(K.X.b[[k]][i,]/k.X.b[[k]][i]))}))
    
    O<- rowSums(sapply(1:n, function (i) { return((status[i])*(K.X.b[[k]][i,]/k.X.b[[k]][i]))}))
    
    
    return(O/D)
  }
  ##### Note: The adjusted kernel,  K.b[k,b,,]/k.b[k,b,], has rowSums equal one.
  
  # K.b<-array(0,dim=c(n.grid[1],n.grid[1]))
  # k.b<-array(0,dim=c(n.grid[1]))
  d <- ncol(X) 
  n <- nrow(X)
  
  
  
  if(is.null(x.grid)) x.grid<-lapply(1:d,function(k) X[order(X[,k]),k])
  
  
  if(is.null(x.min)) x.min<-sapply(x.grid,head,1)
  if(is.null(x.max)) x.max<-sapply(x.grid,tail,1)
  
  x.grid2<-lapply(1:d, function(k) seq(x.min[k],x.max[k],length=n.grid.additional))
  x.grid<-lapply(1:d,function(k) sort(c(x.grid[[k]],x.grid2[[k]])))
  
  
  n.grid<-sapply(x.grid, length)
  
  
  
  
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
  
  x.grid.array<-matrix(rep(x.grid[[1]], times=n.grid[1]),nrow = n.grid[1], ncol = n.grid[1],byrow=FALSE)   # 
  u<-x.grid.array-t(x.grid.array) # u is x_j-u_j for all x,u on the grid considered
  
  
  K.b<-apply(u/bandwidth[1],1:2,kern)/(bandwidth[1])
  k.b<-colSums(dx[[1]]*apply(u/bandwidth[1],1:2,kern)/(bandwidth[1]))
  
  
  X<-t(X)
  
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
      k.X.b[[k]][i]<-sum(dx[[k]]*sapply(u/bandwidth[k],kern)/(bandwidth[k]))
    } }
  
  
  alpha_backfit<-list()
  
  if (is.null(initial)){
    for(k in 1:d){
      alpha_backfit[[k]]<-rep(1, n.grid[k])
    }
  } else  alpha_backfit<-initial
  
  for (l in 2:it)
  {  
    alpha_backfit_old<-alpha_backfit
    for(k in 1:d)
    {
      
      
      
      alpha_backfit[[k]]<- get.new.alpha(status,alpha_backfit,K.X.b,k.X.b,K.b,k.b,Y,k,x.grid,dx,n.grid,d,n)
      
      
      # for (j in 1:(d-1))
      #{   alpha_backfit[[j]]<- 1*(alpha_backfit[[j]]/ alpha_backfit[[j]][1]) }
      
      
      alpha_backfit[[k]][is.nan(alpha_backfit[[k]])]<-0
      alpha_backfit[[k]][alpha_backfit[[k]]==Inf]<-0
      
    }
    
    if (max(abs(unlist(alpha_backfit_old)-unlist(alpha_backfit)),na.rm=TRUE)<=0.001) break
    print(c(l,max(abs(unlist(alpha_backfit_old)-unlist(alpha_backfit)),na.rm=TRUE)))
  }
  return(list(alpha_backfit=alpha_backfit,l=l,x.grid=x.grid))
}
