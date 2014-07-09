n <- sample.size <- 20
p <- 3


# DATA SIMULATION
uniqueness=1
mix.w1=0.7
mix.w2=0.3

mix.mu1=c(2,1,1)
mix.mu2=c(-1,-1, -2)


w1<-mix.w1
w2<-mix.w2
mu1=mix.mu1
mu2=mix.mu2

  U<-runif(n=sample.size, min=0, max=1)
  U[U[]<=w1]<-1
  U[U[]!=1]<-0
  library(MASS)
  F1<-mvrnorm(sample.size, mu=mu1, Sigma=matrix(c(1,.8,0.5,.8, 1, 0.8, 0.5,0.8,1), ncol=3, nrow=3, byrow=TRUE))
  F2<-mvrnorm(sample.size, mu=mu2, Sigma=diag(c(1,1,1)))
  data<-U*F1+(1-U)*F2+mvrnorm(sample.size, mu=c(0,0,0), Sigma=diag(c(1,1,1)))
  ii<-order(y1<-data[,1], y2<-data[,2], y3<-data[,3])
  data1<-cbind(y1[ii],y2[ii], y3[ii])
  y1<-as.vector(data1[,1])
  y2<-as.vector(data1[,2])
  y3<-as.vector(data1[,3])
plot(as.data.frame(data1))
plot(y1); plot(y2);plot(y3);

#checking what is this mixture of normals
  mixnorm<-function(y1,y2)
  {
  rep<-c()
    for (i in 1:length(y1)){

   rep<-c(rep, w1*dmvnorm(c(y1[i], y2[i]), mean=mu1, sigma=matrix(c(1,.8,.8, 1), ncol=2, nrow=2, byrow=TRUE), log=FALSE)+
   w2*dmvnorm(c(y1[i], y2[i]), mean=mu2, sigma=diag(c(1,1)), log=FALSE))

   }
   rep
  }
  

  yy1<-seq(-5,5,by=0.1)
  yy2<-seq(-5,5,by=0.1)
  mu1=c(2,1);mu2=c(-1,-1);
  yy3<-outer(yy1,yy2,mixnorm)

  persp(x=yy1,y=yy2,yy3, theta =20, phi = 20, expand = 0.5, col = "white",
  ticktype="detailed",nticks=5,cex.main=1, main=NULL)

  mu1=c(1,1);mu2=c(-1,-2);
  yy3<-outer(yy1,yy2,mixnorm)

  persp(x=yy1,y=yy2,yy3, theta =20, phi = 20, expand = 0.5, col = "white",
  ticktype="detailed",nticks=5,cex.main=1, main=NULL)

  mu1=c(2,1);mu2=c(-1,-2);
  yy3<-outer(yy1,yy2,mixnorm)

  persp(x=yy1,y=yy2,yy3, theta =20, phi = 20, expand = 0.5, col = "white",
  ticktype="detailed",nticks=5,cex.main=1, main=NULL)

data=data1
# end of DATA SIMULATION


# COMPUTE INITIAL VALUES

S <- grid.size <- 60
grid.eta<-seq(0.001,0.99,length.out=S)
eta.seq<-seq(0.001,0.99,length.out=100)

library(splines)
K <- n.knots <- 10
M <- 4 # order (=degree+1) of B-spline polynomial
m<-M+K
Im<-diag(rep(1,m))
Bj1<-bs(eta.seq, intercept = TRUE, df=m)
#this B is also B11(eta), B21, B31 etc in the initial computation
Bjk<-bs(eta.seq, intercept = FALSE, df=m)
#this B is also B22(eta), B32,B33 in the initial computation

 B.ar <- array(0,c(n, (M+K)*p, p))
# n rows, (M+K)*p columns, p tables
# B.ar[,1:(j*(M+K)),j] is Bj(eta)
# B.ar[i,,] is B(etai)

for (j in 1:p){
  B.ar[,1:m,j] <- Bj1
  B.ar[,(1+m*(j-1)):(m*j),j] <- Bjk
 }

# checking:
# B.ar[,1:(M+K),1]
# B.ar[,1:2*(M+K),2]
# B.ar[5,,2]   #B.ar[i,,j]
# B.ar[50,,]


sp1<-lm(data[,1]~0+Bj1)
summary(sp1)
plot(eta.seq, data[,1], xlab = "x", ylab = "y")
lines(eta.seq, predict(sp1, as.data.frame(eta.seq)), col="red", lwd=2, lty="dotted")
beta11<-as.vector(sp1$coefficients)


sp2<-lm(data[,2]~0+Bjk, offset=Bj1%*%beta11)
plot(eta.seq, data[,2], xlab = "x", ylab = "y")
lines(eta.seq, Bj1%*%beta11, col="blue", lwd=2, lty="dotted")
beta21<-beta11
beta22<-sp2$coefficients
lines(eta.seq, Bj1%*%beta11+Bjk%*%beta22, col="green", lwd=2, lty="dotted")
lines(eta.seq, Bjk%*%beta22, , col="red", lwd=2, lty="dotted")
                                       # should add legend


sp3<-lm(data[,3]~0+Bjk, offset=Bj1%*%beta11+Bjk%*%beta22)
plot(eta.seq, data[,3], xlab = "x", ylab = "y")
beta31<-beta11
beta32<-beta22
beta33<-sp3$coefficients
lines(eta.seq, Bj1%*%beta11+Bjk%*%beta22+Bjk%*%beta33, col="green", lwd=2, lty="dotted")


beta <- c(beta11, beta21, beta22, beta31,beta32, beta33)
#beta<-c(rep(0, 0.5*p*(p+1)*m))

totdraws <- 1000
mu.proc.ar <-array(0,c(n, 0.5*p*(p+1), totdraws))
mu.proc.ar[,1,1] <-Bj1%*%beta11
mu.proc.ar[,2,1] <-Bj1%*%beta21
mu.proc.ar[,3,1] <-Bjk%*%beta22
mu.proc.ar[,4,1] <-Bj1%*%beta31
mu.proc.ar[,5,1] <-Bjk%*%beta32
mu.proc.ar[,6,1] <-Bjk%*%beta33

## beta
# betaj is of size j*m
# beta=c(beta1, ..., betap) is of size m*0.5*(1+p)*p
# betaj=beta[(0.5*j*(j-1)*m+1):(0.5*j*(j+1)*m)]

tau <- as.vector(rep(1,0.5*p*(p+1)))
a.tau <- 1
b.tau <- 0.05

sigma.vec<-c(rep(1,p))
a.sigma <-1
b.sigma <-0.005

totdraws <- 100

source("fcts_plotting.R")
source("fcts_eta_Bar.R")

eta<-matrix(0, ncol=p, nrow=n)
for (j in 1:p) { eta[,j]=seq(0.001,0.99,length.out=n) }

Eta.ar <-array(0,c(n, p, totdraws))
# n rows, p columns, totdraws tables
Eta.ar[,,1] <- eta

sigma.ar<-matrix(0, ncol=p, nrow=totdraws)
beta.ar<-matrix(0, ncol=0.5*p*(p+1)*m, nrow=totdraws)
beta.ar[1,]<-beta

