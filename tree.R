require(MASS)
require(ranger)
require(survival)
library(randomForestSRC)
library(riskRegression)
library(prodlim)
# load some functions 
source("~/tmp/mhdata.R")
set <- list(d=5,rho=0,model=1,violate.cox=1,seed=10)
large.data <- do.call("mhdata",c(list(n=20000),set))
train.data <- do.call("mhdata",c(list(n=500),set))
pred <- paste0("X", 1:(set$d-1))
frmla<-reformulate(pred,"Surv(time, status)")
# plot(prodlim(Hist(time,status)~1,data=large.data),xlim=c(0,100))


# ranger
forest1 <-ranger(frmla, data=train.data,importance = "none",num.trees=1000,mtry=set$d,verbose = TRUE)
# randomForestSRC
forest2 <-rfsrc(frmla,data=train.data,importance = "none",num.trees=1000,mtry=set$d)
# cox
cox.fit<-coxph(frmla , data=train.data)
# munir's stuff
it=30
b.grid<-numeric(set$d)
b.grid[1]<-0.000001
b.grid[2:set$d]<-rep(0.2,set$d-1)
alpha_backfit<-SBF.MH.LC(train.data,b.grid,it=it)
alpha_backfit2<-SBF.MH.CLL(train.data,b.grid,weight='sw',it=it,x.grid=NULL,integral.approx='midd',kcorr='TRUE',LC='FALSE')


#####  now fit trained models
set.seed(82423)   #363, 30462
seed<-sample(1:99999,1) 
##### predicted survival for a single new subject
test.data <- do.call("mhdata",c(list(n=2),replace(set,"seed",seed)))[,-c(1,2)]
forest1.fit<-1-predictRisk(forest1, newdata=test.data[1:2,,drop=FALSE],times=sort(unique(train.data$time)))[1,]
forest2.fit<-predict(forest2, test.data[1,,drop=FALSE])
cox.pred<-survfit(cox.fit, test.data[1,,drop=FALSE])
predict.sbf.LC<- predict.sbf(alpha_backfit,test.data,1)
predict.sbf.LL<- predict.sbf(alpha_backfit2,test.data,1)

# plot true survival curve
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#D55E00", "#0072B2", "#CC79A7", "#F0E442")
true_par <- exp(mhrate(Z=test.data[1,],set$violate.cox))
plot.points.index<-pexp(forest1$unique.death.times,rate=true_par,lower.tail = FALSE)>0.001
plot.points<-forest1$unique.death.times[plot.points.index]
plot(plot.points,pexp(plot.points,rate=true_par,lower.tail = FALSE),type='l',lwd=3)
lines(sort(unique(train.data$time)),forest1.fit,col=cbbPalette[3],lwd=3)
lines(cox.pred$time,cox.pred$surv,col=cbbPalette[2], lwd=3)
lines(forest2.fit$time,as.numeric(forest2.fit$survival),col=cbbPalette[6], lwd=3)
lines(predict.sbf.LC[[2]],predict.sbf.LC[[1]][1,],col=cbbPalette[4],lwd=3)
lines(predict.sbf.LL[[2]],predict.sbf.LL[[1]][1,],col=cbbPalette[5],lwd=3)
legend('topright', c('true','cox', 'random forest', 'sbf.LC','sbf.CLL'),col=cbbPalette, lty=1,lwd=3)





