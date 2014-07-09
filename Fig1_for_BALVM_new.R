f12 <- function( y1, y2,  theta11=0, theta21=-3, theta31=5,
			         	 theta12=0, theta22=5, theta32=3 )
# joint pdf (probability density function) of y1 and y2
# density of bivariate mixture N(theta_1, I)/3 + N(theta_2, I)/3 + N(theta_3, I)/3
# theta_1=c(theta11, theta12);
# theta_2=c(theta21, theta22);
# theta_3=c(theta31, theta32)
# dnorm are used for easier vectorization
{
 f12 <-  dnorm( y1, mean = theta11, sd = 1 )*dnorm( y2, mean = theta12, sd = 1 )/3+
         dnorm( y1, mean =theta21, sd = 1 )*dnorm( y2, mean = theta22, sd = 1 )/3+
         dnorm( y1, mean =theta31, sd = 1 )*dnorm( y2, mean = theta32, sd = 1 )/3
  return( f12 )
}

F1.ex <- function( y1, theta11=0, theta21=-3, theta31=5 )
# y1 marginal cdf (cumulative density function)
{

  f1 <- pnorm( y1, mean = theta11, sd = 1 )/3+
        pnorm( y1, mean = theta21, sd = 1 )/3+
        pnorm( y1, mean = theta31, sd = 1 )/3
  return( f1)
}

f1.ex <- function( y1, theta11=0, theta21=-3, theta31=5)
# marginal pdf of y1 (probability density function)
{
 f1 <- dnorm(y1, mean = theta11, sd = 1)/3+
       dnorm(y1, mean = theta21, sd = 1)/3+
       dnorm(y1, mean = theta31, sd = 1)/3
  return(f1)
}

f2.1 <- function( y2, y1,  theta11=0, theta21=-3, theta31=5,
				          theta12=0, theta22=5, theta32=3)
# conditional pdf of y2|y1
{
 intermed.result <- exp( log( f12( y1, y2,  theta11, theta21, theta31, theta12, theta22, theta32 ) )-
			                   log( f1.ex( y1,  theta11, theta21, theta31 ) ) )
 result <- c( rep( 0, times=length( intermed.result ) ) )
 result[ which( is.finite( intermed.result ) ) ] <- intermed.result[ which( is.finite( intermed.result ) ) ]
 # we can take only finite values as we know that all the infinite values are just purely numerical artefacts
 return( result )
}

F2.1.ex <- function( y2, y1, theta11=0, theta21=-3, theta31=5,
				             theta12=0, theta22=5, theta32=3)
# conditional cdf of y2|y1
# approximate integrattion with quadrature
{
 integrate( f=f2.1, lower=-Inf, upper=y2, y1, theta11, theta21, theta31,
				    theta12, theta22, theta32, subdivisions=10000 ) 
} 

find.y2 <- function( u1, u2,  theta11=0, theta21=-3, theta31=5,
				             theta12=0, theta22=5, theta32=3 )
# conditional quantile function F_{2|1}(U_2|U_1; ...)
# finds y2 given values of marginal cdf F1 in u1 and
# the value of conditional cdf F2.1 in u1 and u2
# using approximate integrattion with quadrature
{
	F1.shift.cdf <- function( y1, u1, theta11, theta21, theta31 ){ F1.ex(y1, theta11, theta21, theta31)-u1 }
	y1.val <- uniroot( f=F1.shift.cdf, interval=c( -9, 9 ), u1=u1, theta11, theta21, theta31 )$root

  F2.1.shift.cdf <- function( y2, y1, u2, theta11=0, theta21=-3, theta31=5,
				                    theta12=0, theta22=5, theta32=3 ){
                                F2.1.ex( y2=y2, y1=y1, theta11, theta21, theta31,
				                        theta12, theta22, theta32 )$value-u2
                        }
	y2 <- uniroot( f=F2.1.shift.cdf, interval=c( -9, 9 ), y1=y1.val, u2=u2,
                theta11, theta21, theta31,
				        theta12, theta22, theta32 )$root
 return(y2)
}


rmixnorm <- function( weights, mus, sample.size=10000 )
# this fct simulate realizations from univariate mixture of normals
# with weights and mus, standard deviations all being one
{
	components <- sample( 1:3, prob=c( weights ), size=sample.size, replace=TRUE )
	sds<-1
#	sds <- c( 1, 1, 1 ) # if want different sd
	samples <- rnorm( n=sample.size, mean=mus[components], sd=sds ) # if want different sd put # sd=sds[components]
	samples
}

F2.1.ex.MC <- function( y2, y1, sample.size=10000,
                        theta11=0, theta21=-3, theta31=5,
				                theta12=0, theta22=5, theta32=3 )
# conditional cdf of y2|y1
# approximate Monte-Carlo integration
{
   w <- c( dnorm( y1, mean=theta11, sd=1 ), dnorm( y1, mean=theta21, sd=1 ), dnorm( y1, mean=theta31, sd=1 ) )
   w <- w/ sum( w )

   y2.samples <- rmixnorm( weights=w, mus=c( theta12, theta22, theta32 ), N=sample.size )

   smpls.smllr.y2 <- 0+ ( y2.samples<= y2 )
   smaller.y2 <- y2.samples[ which( smpls.smllr.y2==1 ) ]
   bigger.y2<-y2.samples[ which( smpls.smllr.y2==0 ) ]
   full.F2.1<-sum( f2.1( y2=y2.samples, y1, theta11, theta21, theta31,
  				    theta12, theta22, theta32 ) )/length( y2.samples )
   if ( length( smaller.y2 )>=length( bigger.y2 ) )
  	F2.1<-sum( f2.1( y2=smaller.y2, y1, theta11, theta21, theta31,
  				    theta12, theta22, theta32) )/
  ( full.F2.1*length( smaller.y2 ) ) else F2.1<-1-sum( f2.1( y2=bigger.y2, y1, theta11, theta21, theta31,
  				    theta12, theta22, theta32 ) )/( full.F2.1* length( bigger.y2 ) )
   if ( F2.1<0 ) F2.1<-0
   F2.1
} 

# F2.1.ex.MC( y2=6, y1=10, sample.size=100000 )
# F2.1.ex( y2=6, y1=10 )$value

find.y2.MC<-function( u1, u2, sample.size=100000, theta11=0, theta21=-3, theta31=5,
				    theta12=0, theta22=5, theta32=3 )
# conditional quantile function F_{2|1}(U_2|U_1; ...)
# so finds y2 given marginal cdf F1, conditional cdf F2.1 and u1, u2
# using approximate Monte-Carlo integration
{
	F1.cdf<-function( y1, u1, theta11, theta21, theta31 ){ F1.ex(y1, theta11, theta21, theta31)-u1 }
	y1.val<-uniroot( f=F1.cdf, interval=c( -9, 9 ), u1=u1, theta11, theta21, theta31 )$root
	# cat( "y1.val=", y1.val, "\n" )
	F2.1.cdf<-function( y2, y1, u2, theta11, theta21, theta31, theta12, theta22, theta32 )
  { F2.1.ex.MC( y2, y1, sample.size=sample.size, theta11, theta21, theta31,
				    theta12, theta22, theta32 )-u2 }
	y2<-uniroot( f=F2.1.cdf, interval=c( -9, 9 ), y1=y1.val, u2=u2, theta11, theta21, theta31,
				    theta12, theta22, theta32 )$root
 return(y2)
}
# to check
 find.y2.MC( u1=0.1, u2=0.9, sample.size=100000, theta11=0, theta21=-3, theta31=5,
				    theta12=0, theta22=5, theta32=3 )
 find.y2( u1=0.1, u2=0.9,  theta11=0, theta21=-3, theta31=5,
				    theta12=0, theta22=5, theta32=3 )


f13<-function( y1, y3, sigma3, theta11=0, theta21=-3, theta31=5,
	theta12=0, theta22=5, theta32=3, c1=0.65, c2=0.5)
# joint pdf of y1 and y3
{
	cy1<-exp( c1*y1-c2 )
	cy1<-cy1/( 1+cy1 )
	theta12<-cy1*theta12; theta22<- cy1*theta22; theta32<-cy1*theta32;
	sd.sqrt<-sqrt( cy1^2+sigma3^2 )

	f13<-dnorm( y1, mean = theta11, sd = 1 )*dnorm( y3, mean = theta12, sd =sd.sqrt )/3+
	dnorm( y1, mean =theta21, sd = 1 )*dnorm( y3, mean = theta22, sd = sd.sqrt )/3+
	dnorm( y1, mean =theta31, sd = 1 )*dnorm( y3, mean = theta32, sd = sd.sqrt )/3
	return( f13 )
}

f23.MC<-function( y2, y3, sigma3, sample.size,  theta11=0, theta21=-3, theta31=5,
	theta12=0, theta22=5, theta32=3, c1=0.65, c2=0.5
 ) 
# joint pdf of y2 and y3
{
	w<-c(1/3, 1/3, 1/3)
	y1.samples<-rmixnorm( weights=w, mus=c( theta12, theta22, theta32 ), N=sample.size )
	cy1<-exp( c1*y1.samples-c2 )
	cy1<-cy1/( 1+cy1 )
	the<-cy1*y2; 

	f23.integrand<-( dnorm( y1.samples, mean = theta11, sd = 1 )*dnorm( y2, mean = theta12, sd =1 )+
		dnorm( y1.samples, mean =theta21, sd = 1 )*dnorm( y2, mean = theta22, sd = 1 )+
		dnorm( y1.samples, mean =theta31, sd = 1 )*dnorm( y3, mean = theta32, sd = 1 ) )*
		dnorm( y3, mean = the, sd = sigma3 )/3
	f23<-sum( f23.integrand )/length( y1.samples )
	return( f23 )
}

# to check
  f23.MC( y2=0, y3=0, sigma3=1, sample.size=1000,  theta11=0, theta21=-3, theta31=5,
	theta12=0, theta22=5, theta32=3, c1=0.65, c2=0.5)

f2.ex<-function( y2, theta12=0, theta22=5, theta32=3 )   
# marginal pdf of y1
{
 f2<-dnorm(y2, mean = theta12, sd = 1)/3+
 dnorm(y2, mean = theta22, sd = 1)/3+
 dnorm(y2, mean = theta32, sd = 1)/3
 return(f2)
}

f3.2<-function( y3, y2, sigma3=1, sample.size, theta11=0, theta21=-3, theta31=5,
	               theta12=0, theta22=5, theta32=3 )
# conditional pdf of y3|y2
{
 intermed.result<-exp( log( f23.MC( y2, y3, sigma3=1, sample.size=sample.size, theta11, theta21, theta31,
	                                   theta12, theta22, theta32 ) )-
                       log( f2.ex( y2, theta12, theta22, theta32 ) ) )
 #result<-c( rep( 0, times=length( intermed.result ) ) )
 #result[ which( is.finite( intermed.result ) ) ]<-intermed.result[ which( is.finite( intermed.result ) ) ]
 # we can take only finite values as we know all the infinite values are just purely numerical artefacts
 result<-intermed.result
 return( result )
}


# to check
f3.2( y3=0, y2=1,sigma3=1, sample.size=1000, theta11=0, theta21=-3, theta31=5,
	theta12=0, theta22=5, theta32=3 )

  #important



F3.2.ex<-function( y3, y2,sigma3=1, theta11=0, theta21=-3, theta31=5,
	theta12=0, theta22=5, theta32=3, sample.size=1000, div=10000 )
# conditional cdf of y3|y2
# approximate integrattion with adaptive quadrature
{
vf3.2 <- Vectorize(f3.2)
 integrate( f=vf3.2, lower=-Inf, upper=y3, y2, sigma3=1, theta11, theta21, theta31,
	                   theta12, theta22, theta32, sample.size=sample.size, subdivisions=div )
}
# to check
F3.2.ex( y3=0, y2=1, sigma3=1, theta11=0, theta21=-3, theta31=5,
	         theta12=0, theta22=5, theta32=3, sample.size=1000, div=10000 )




f3.1<-function( y3, y1,sigma3=1, theta11=0, theta21=-3, theta31=5,
	theta12=0, theta22=5, theta32=3, c1=0.65, c2=0.5, sample.size=1000, div=10000  )
# conditional pdf of y3|y1
{
 intermed.result <- exp( log( f13( y1, y3, sigma3, theta11, theta21, theta31,
	                               theta12, theta22, theta32, c1, c2) )-
                        log( f1.ex( y1, theta11, theta21, theta31 ) ) )
 result<-c( rep( 0, times=length( intermed.result ) ) )
 result[ which( is.finite( intermed.result ) ) ]<-intermed.result[ which( is.finite( intermed.result ) ) ]
 # we can take only finite values as we all the infinite values are just purely numerical artefacts
 return( result )
}

F3.1.ex <- function( y3, y1, sigma3=1, theta11=0, theta21=-3, theta31=5,
	theta12=0, theta22=5, theta32=3, c1=0.65, c2=0.5, div=10000)
# conditional cdf of y3|y1
# approximate integrattion with adaptive quadrature
{
  vf3.1 <- vectorize(f3.1)
 integrate( f=vf3.1, lower=-Inf, upper=y3, y1, sigma3, theta11, theta21, theta31,
	theta12, theta22, theta32, c1, c2, subdivisions=div )
}

F1.cdf <- function(y1, u1, theta11=0, theta21=-3, theta31=5)
  {
      F1.ex(y1, theta11, theta21, theta31)-u1
  }
F3.1.cdf<-function( y3, y1, u3, sigma3=1, theta11=0, theta21=-3, theta31=5,
	theta12=0, theta22=5, theta32=3, c1=0.65, c2=0.5, div=10000  )
  {
    F3.1.ex( y3=y3, y1, sigma3, theta11, theta21, theta31,
             theta12, theta22, theta32, c1, c2, div )$value-u3
  }


find.y3 <- function( u1, u3, sigma3=1, theta11=0, theta21=-3, theta31=5,
	theta12=0, theta22=5, theta32=3, c1=0.65, c2=0.5, div=10000 )
# conditional quantile function F_{3|1}(U_3|U_1; ...)
# finds y3 given marginal cdf F1, conditional cdf F3.1 and u1, u2
# using approximate Monte-Carlo integration
{

	y1.val<-uniroot( f=F1.cdf, interval=c( -9, 9 ), u1=u1, theta11, theta21, theta31)$root
	# cat( "y1.val=", y1.val, "\n" )
	y3 <- uniroot(f=F3.1.cdf, interval=c( -9, 9 ), y1=y1.val, u3=u3,
              sigma3, theta11, theta21, theta31,
              theta12, theta22, theta32, c1, c2, div)$root
 return(y3)
}
#find.y3( 0.1, 0.3 )


# PLOTTING

#postscript("example_density.eps",horizontal=FALSE, paper="special",height=8,width=12)
pdf("example_density.pdf", paper="special",height=8,width=12)


layout(matrix(c(1:6), ncol=3, nrow=2, byrow=TRUE))
par(mar=c(5,5,2,2), mgp=c(2.5, 1, 0))

# contours of density f12 - (a) picture
y1<-seq(from=-6, to=8, by=0.01)
y2<-seq(from=-4, to=8, by=0.01)
z1<-outer(y1,y2, FUN=f12)
contour(y1,y2,z1, method="edge", xlab=expression(y[1]), ylab=expression(y[2]),
	levels=c(0.0000001, 0.00001, 0.0001, 0.001, 0.005,  0.02,  0.05), labcex = 1.5,
	cex.axis=2,cex.lab=2.5, drawlabels=F, main="(a)", cex.main=2.5)

# contours of density f23 - (b) picture

y2<-seq(from=-3, to=7, by=0.1); y3<-seq(from=-4, to=6, by=0.1)
#f23.contours<-matrix(ncol=length( y2 ), nrow=length( y3 ) )

#for (j in 1:length( y2 ) )
#{
#  for (i in 1:length( y3 ) ) f23.contours[i,j]<-f23.MC( y2=y2[i], y3=y3[j], sigma3=1, sample.size=40000 )
#}
contour(y2, y3, f23.contours, method="edge", xlab=expression(y[2]), ylab=expression(y[3]),labcex = 1.5,
	levels=c(0.000001, 0.00001, 0.0001, 0.0003,  0.001, 0.005,  0.02,  0.05), 
	cex.axis=2,cex.lab=2.5, drawlabels=F, main="(b)", cex.main=2.5)

# contours of density f13 - (c) picture

y1<-seq(from=-6, to=8, by=0.1); y3<-seq(from=-4, to=6, by=0.1)
f13.dens<-outer(y1, y3, FUN=f13, sigma3=1 ) 

contour(y1,y3,f13.dens, method="edge", xlab=expression(y[1]), ylab=expression(y[3]),
	levels=c(0.000001, 0.00001, 0.0001, 0.0003,  0.001, 0.005,  0.02,  0.05), 
	labcex = 1.5,
	cex.axis=2,cex.lab=2.5, drawlabels=F, main="(c)", cex.main=2.5)


# constant nu2 (for all possible y2 values) contours as functions of u1 and u2  - (d) picture

u1.vec<-seq(from=0.000001, to=0.999999, by=0.01)
u2.vec<-seq(from=0.000001, to=0.999999, by=0.01)
z12<-matrix(ncol=length(u2.vec), nrow=length(u1.vec))
#computing matrix of nu2 values
for (j in 1:length(u2.vec))
{
  for (i in 1:length(u1.vec)) z12[i,j]<-find.y2(u1=u1.vec[i],  u2=u2.vec[j])
}
contour(u1.vec,u2.vec,z12, method="edge", xlab=expression(u[1]), ylab=expression(u[2]),labcex = 1.5,
	cex.axis=2,cex.lab=2.5, drawlabels=F, main="(d)", cex.main=2.5)

# nu2 function (for one possible y2 values) for u2=0.2 - (e) picture
#u2.1.val<-0.2
#u1.vec<-seq(from=0.0001, to=0.9999, by=0.001)
#nu2.1.vec<-c(); for (i in 1:length(u1.vec)) nu2.1.vec<-c(nu2.1.vec, find.y2(u1=u1.vec[i], u2=u2.1.val))
#
#plot(u1.vec, nu2.1.vec, type="l", xlab=expression(u[1]),
#	ylab=expression(nu[2](u[1], 0.2)), cex.axis=2,cex.lab=2.5, main="(e)", cex.main=2.5)

# constant nu3 (for all possible y1 values) contours as functions of u1 and u3  - (f) picture

u1.vec<-seq(from=0.000001, to=0.999999, by=0.01)
u3.vec<-seq(from=0.000001, to=0.999999, by=0.01)
z13<-matrix(ncol=length(u3.vec), nrow=length(u1.vec))
#computing matrix of nu2 values
for (j in 1:length(u3.vec))
{
  for (i in 1:length(u1.vec)) z13[i,j]<-find.y3(u1=u1.vec[i],  u3=u3.vec[j])
}
contour(u1.vec,u3.vec,z13, method="edge", xlab=expression(u[1]), ylab=expression(u[3]),labcex = 1.5,
	cex.axis=2,cex.lab=2.5, drawlabels=F, main="(e)", cex.main=2.5)
dev.off()
