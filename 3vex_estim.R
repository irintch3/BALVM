
D.ex<-read.table(file="DavidEx.txt", header=FALSE)
p=ncol(D.ex)
n=nrow(D.ex)

ord=4
m=10
S <- grid.size <- 300
grid.eta<-seq(0.00000001,0.999999,length.out=S)
eta.seq<-seq(0.00000001,0.999999,length.out=n)


Bj1<-bs(u1,knots=AknotsI[(ord+1):(length(AknotsI)-ord)], intercept = TRUE, df=m,
        Boundary.knots = c(1e-10, (1-1e-10)))
#bs(eta.seq, intercept = TRUE, df=m)
Bjk<-bs(u2, knots=Aknots[(ord+1):(length(Aknots)-ord)], intercept = FALSE, df=m,
        Boundary.knots = c(1e-10, (1-1e-10))) #
#bs(eta.seq, intercept = FALSE, df=m)


B.ar <- array(0,c(n, m*p, p))
for (j in 1:p){
  B.ar[,1:m,j] <- Bj1
  cat("j=",j)
  if (j>=2){ for (k in 2:j){
    cat(" k=",k,"\n")
    B.ar[,(1+m*(k-1)):(m*k),j] <- Bjk}
  }}


data=D.ex
sp1<-lm(data[,1]~0+Bj1)
beta11<-as.vector(sp1$coefficients)
beta21<-c(rep(0.001,m))
sp2<-lm(data[,2]~0+Bjk, offset=Bj1%*%beta21)
beta22<-sp2$coefficients
sp3<-lm(data[,3]~0+Bjk, offset=Bj1%*%beta11+Bjk%*%beta22)
beta31<-beta11
beta32<-beta22
beta33<-sp3$coefficients


Rallbeta<-beta <- c(beta11, beta21, beta22, beta31,beta32, beta33)



eta <- matrix(eta.seq, ncol=p, nrow=n)
Rsigma.vec <- c(rep(1,p))


intercept=FALSE
nIknots <-m-ord+(1-intercept)
knots=seq.int(from=0, to=1,length.out=nIknots+2)[-c(1,nIknots+2)]
knots<-stats::quantile(eta.seq, knots)
Boundary.knots = c(1e-8, (1-1e-8))
Aknots=sort(c(rep(Boundary.knots, ord), knots))
nknots=length(Aknots)

intercept=TRUE
nIknots <-m-ord+(1-intercept)
knotsI=seq.int(from=0, to=1,length.out=nIknots+2)[-c(1,nIknots+2)]
knotsI<-stats::quantile(eta.seq, knotsI)
Boundary.knots = c(1e-8, (1-1e-8))
AknotsI=sort(c(rep(Boundary.knots, ord), knotsI))

dyn.load("Naive.so")
Rtot_draws <- as.integer(2000)#0000
Rdens=0
Rdrws.step <- 1
RBurn<-1900
setwd("/home/staff/ii9/Duke/R/Dex123")

mu.all.my <- .Call("gibbs_BALVM_naive", Rtot_draws, Aknots, 
                   as.integer(nknots), AknotsI, as.integer(ord),
                   as.integer(n), as.integer(m), as.integer(p), 
                   as.integer(S), grid.eta, as.matrix(D.ex), c(1,0.05),
                   c(1,0.05), eta, Rallbeta, as.integer(Rdrws.step), 
                   as.integer(RBurn) )
setwd("/home/staff/ii9/Duke/R")
write.table(file="Sig123.txt", t(mu.all.my), col.names=FALSE, row.names=FALSE)
dyn.unload("Naive.so")



newd<-cbind(D.ex[,3], D.ex[,1],D.ex[,2])
dyn.load("Naive.so")
Rtot_draws <- as.integer(2000)#0000
Rdens=0
Rdrws.step <- 1
RBurn<-1900
setwd("/home/staff/ii9/Duke/R/Dex312")
mu.all.my <- .Call("gibbs_BALVM_naive", Rtot_draws, Aknots, 
                   as.integer(nknots), AknotsI, as.integer(ord),
                   as.integer(n), as.integer(m), as.integer(p), 
                   as.integer(S), grid.eta, as.matrix(newd), c(1,0.05),
                   c(1,0.05), eta, Rallbeta, as.integer(Rdrws.step), 
                   as.integer(RBurn) )
setwd("/home/staff/ii9/Duke/R")
write.table(file="Sig312.txt", t(mu.all.my), col.names=FALSE, row.names=FALSE)
dyn.unload("Naive.so")
