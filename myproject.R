library(MASS)
library(randomForestSRC)

ro<-0.5
n=100
p=20  

design.cov <- function(p, cov.x, equalCorr=F){
  if (!equalCorr){
    x <- matrix(0,p,p)            
    for(i in 1:p){
      x[i,(i:p)] <- cov.x^(0:(p - i))
      x[i,(1:(i - 1))] <- x[(1:(i - 1)),i]
    }
  }
  else{
    x <- matrix(cov.x,p,p)
    diag(x) <- 1
  }               
  return(x)
}

mvngeneration <- function(n, p, varcov){
  choleski <- chol(varcov)
  tcholeski <- t(choleski)
  z <- matrix(rnorm(p * n, 0, 1), p, n)
  tranx <- tcholeski %*% z 
  x <- t(tranx)
  x
}	

simdata<- mvngeneration(n, p, design.cov(p, ro))

simdata<-t(simdata)
rownames(simdata)<-paste("x",1:20,sep="")
colnames(simdata)<-paste("id",1:100,sep="")

result<-rep(0,n)
for(i in 1:p){result<-result+simdata[i,]}
result<-result+10*simdata[4,]*simdata[17,]

training<-t(simdata)

data<-data.frame(result,training)
names(data)[1]<-c("y")

formula<-y~.

###try setting 1
sim.out<-rfsrc(formula,data)
sim.int<-find.interaction(sim.out,method="maxsubtree")

sim.diag<-diag(sim.int)
sim.colsum<-colSums(sim.int)-sim.diag
sim.colmean<-sim.colsum/(length(sim.diag)-1)
threshold<-mean(sim.colmean)
sim.colmean[which(sim.colmean<threshold)]<-threshold

sim.out<-rfsrc(formula,data,predictorWt=
                 (1-sim.colmean[match(names(data),names(sim.colmean))][-c(1,2)]))
sim.int<-find.interaction(sim.out,method="maxsubtree")


###try setting 2
fNames <- all.vars(formula, max.names=1e7)
var.columns <- (1:ncol(data))[-c(which(names(data) == fNames[1]))]

sim.out<-rfsrc(formula,data)
sim.int<-find.interaction(sim.out,method="maxsubtree")

sim.diag<-diag(sim.int)
sim.colsum<-colSums(sim.int)-sim.diag
sim.colmean<-sim.colsum/(length(sim.diag)-1)
threshold<-mean(sim.colmean)
sim.colmean[which(sim.colmean>threshold)]<-max(sim.colmean)

sim.out<-rfsrc(formula,data,predictorWt=
                 (1-sim.colmean[match(names(data),names(sim.colmean))][-c(1,2)]))
sim.int<-find.interaction(sim.out,method="maxsubtree")

sim.diag<-diag(sim.int)
threshold1<-mean(sim.diag)
var.pt1<-which(sim.diag<threshold1)
idx<-match(colnames(sim.int),colnames(data))
var.pt<-idx[var.pt1]
drop.var.pt <- setdiff(var.columns, var.pt)
data<-data[,-drop.var.pt]








library(iRF)
x<-data[,-1]
y<-data[,1]
iRF(x,y)


