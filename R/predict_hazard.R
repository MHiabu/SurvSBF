predict.hazard<-function(result, data){
  data<-as.matrix(data)
  if (ncol(data)==1) data<-t(data)
  x.grid<-result$x.grid
  alpha_sbf<-result$alpha_backfit
  hazard <- matrix(nrow=nrow(data), ncol=(length( alpha_sbf)+1))
  
  for(i in 1:nrow(data)){
    for (j in 1:length(alpha_sbf)){
      index<- which.min(abs(as.numeric(data[i,j])-x.grid[[j]]))
      error<-  x.grid[[j]][index] - as.numeric(as.numeric(data[i,j]))
      if (error>0&index>=2) {index2<-index-1} 
      if (error<0) {index2<-index+1} 
      error2<-  result$x.grid[[j]][index2] -   as.numeric(data[i,j])
      
      if (error!=0&index2>=1&index2<=length(x.grid[[j]]))  {
        hazard[i,j] <- (abs(error)*alpha_sbf[[j]][index]+abs(error2)*alpha_sbf[[j]][index2])/(abs(error)+abs(error2))} else{
          hazard[i,j] <- alpha_sbf[[j]][index]
        }
    } 
    hazard[i,length( alpha_sbf)+1]<-prod(hazard[i,-(length( alpha_sbf)+1)])
  }
  
  
  return(hazard)
}
