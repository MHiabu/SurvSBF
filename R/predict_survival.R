predict.survival<-function(result, data, times=NULL){
  data<-as.matrix(data)
  x.grid<-result$x.grid
  alpha_sbf<-result$alpha_backfit
  if (is.null(times)) surv.times<-x.grid[[1]][-1] else surv.times<-times
  surv.prob <- matrix(nrow=nrow(data), ncol=length(surv.times))
  
  hazard.values<-predict.hazard(result,cbind(rep(1,nrow(data)),data))
  
  for(i in 1:nrow(data)){
    # find.index<- error<-numeric(length(alpha_sbf))
    # for (j in 2:length(alpha_sbf)){
    #   error[j] <- min(abs(as.numeric(as.numeric(data[i,j-1]))-x.grid[[j]]))
    #   find.index[j]<- which.min(abs(as.numeric(data[i,j-1])-x.grid[[j]]))
    # }
    
    
    
    par<-1
    for (j in 2:length(alpha_sbf)){
      par <- par * hazard.values[i,j]
    }
    
    if (is.null(times)){
      surv.prob[i,]<-exp(-par*cumsum(alpha_sbf[[1]][-1]*diff(x.grid[[1]])))
    } else{ 
      for (j in 1:length(times)) {
        index <- which.min(abs(result$x.grid[[1]]-times[j]))
        error <-  result$x.grid[[1]][index] - times[j]
        
        if (error==0)  surv.prob[i,j] <- exp(-par*sum(alpha_sbf[[1]][2:index]*diff(x.grid[[1]][1:index]))) else{   if (error>0&index==2) surv.prob[i,j] <- exp(-par*sum(alpha_sbf[[1]][2:index]*diff(x.grid[[1]][1:index])))
        if ((error>0&index>2) | (error<0& index>1)) {
          if (error>0) {index2<-index-1} 
          if (error<0)  {index2<-index+1}
          error2<-   result$x.grid[[1]][index2] - times[j]          
          temp1<-  exp(-par*sum(alpha_sbf[[1]][2:index]*diff(x.grid[[1]][1:index])))
          temp2<-  exp(-par*sum(alpha_sbf[[1]][2:index2]*diff(x.grid[[1]][1:index2])))
          surv.prob[i,j] <- (abs(error)*temp1 + abs(error2)*temp2)/(abs(error)+abs(error2))}
        }
        if (index==1) surv.prob[i,j]  <- 1
      }}
    
  }
  return(list(surv.prob=surv.prob, surv.times= surv.times))
}
