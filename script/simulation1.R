getwd()

####Exercise 4 Monte Carlo Simulation (1)
###Assessing the properties of a test statistic
library(spdep)

set.seed(1000)  ##生成随机数,使结果一致
x <- rnorm(25)
mean(x)
sd(x)

set.seed(1)
xx <- rnorm(100000)            #rnorm
mean(xx)
sd(xx)

y <- rnorm(1000,4,30)
mean(y)
sd(y)

z <- sort(y)                #sort
z[950]

qz <- qnorm(0.95,4,30)      #qnorm
qz

pz <- pnorm(z[950],4,30)   #pnorm
pz

##4.2.2 Practice

z2 <- rpois(50,2)
summary(z2)
mean(z2)
var(z2)
pdf(file="poisson.pdf")   #生成pdf文件
hist(z2)
dev.off()

###4.3 Basic matrix algebra operations
a<- c(1:5)
a
A <- matrix(1:25,5,5)
A
Z <- c(1:5)
B <- array(Z,dim=c(5,5))
B


###4.4 Creating Spatial weights for a grid layout (cell2nb)

rook4x4 <- cell2nb(4,4)
summary(rook4x4)

rook4x4t <- cell2nb(4,4)
summary(rook4x4t)

##4.4.2 Practice

rook5x10 <- cell2nb(5,10)
summary(rook5x10)

queen5x10 <- cell2nb(5,10,type = "queen")
summary(queen5x10)


###4.5 converting spatial weights to a matrix (nb2mat)

wrks <- nb2mat(rook4x4)
wrks

##4.5.2 Practice: Creating a spatial lag using matrix multiplication
x <- rnorm(50)
w50 <- nb2mat(queen5x10)
w50[1,]

##compare to 2 methods
#method1
wx <- w50%*%x                                 #矩阵运算
wx[1:4]
#method2
wwx <- lag.listw(nb2listw(queen5x10),x)           #nb2listw (邻居矩阵转换成权重矩阵)
wwx[1:4]


###4.6 A simulation experiment
rm(list = ls(all=TRUE))
##4.6.1 computing Moran's I using matrix algebra
set.seed(12345)
y <- rnorm(16)
yz <- y-mean(y)
mean(yz)
#denominator
yy <- crossprod(yz)
yy
#numerator
rook4x4 <- cell2nb(4,4)
wy <- lag.listw(nb2listw(rook4x4),yz)
ywy <- crossprod(yz,wy)
ywy
#computation
mi <- ywy/yy
mi

##4.6.2 the theoretical momoents of Moran's I

rks <- spweights.constants(nb2listw(rook4x4))
rks

eirk4x4 <- -1 / rks$n1
eirk4x4
ei2rk4x4 <- (rks$n^2*rks$S1-rks$n*rks$S2+3*rks$S0^2)/(rks$n1*(rks$n+1)*rks$S0^2)
ei2rk4x4
varirk4x4 <- ei2rk4x4 - eirk4x4^2
varirk4x4
sdirk4x4 <- sqrt(varirk4x4)
sdirk4x4

##4.6.3 Inference
#method1
iz <- (mi-eirk4x4)/sdirk4x4
iz

piz <- pnorm(iz,lower.tail = FALSE)*2
piz
#method2
moran.test(y,nb2listw(rook4x4),randomisation = FALSE,alternative = "two.sided")   #results identical

###4.6.5 the simulation loop /循环

r <- 1000
moran <- vector(mode = "numeric",length = r)
reject <- vector(mode="numeric",length=r)
n <- 16
crit <- 0.05
w4 <- nb2listw(rook4x4)

for (i in 1:r) {
  y <- rnorm(n)
  yz <- y-mean(y)
  yy <- crossprod(yz)
  wy <- lag.listw(w4,yz)
  ywy <- crossprod(yz,wy)
  mi <- ywy/yy
  moran[i] <- mi
  iz <- (mi-eirk4x4)/sdirk4x4
  piz <- pnorm(iz,lower.tail = FALSE)*2
  if (piz<crit){reject[i]<1}
}

#4.6.6 analyzing the results

mean(moran)
sd(moran)
summary(moran)
rejfreq <- sum(reject) / r
rejfreq

varsim <- crit*(1-crit)/r
varsim
sdsim <- sqrt(varsim)
sdsim

simint <- c(crit - 2*sdsim,crit + 2*sdsim)
simint


###4.6.7 Practice:simulating moran's I

sim.moran.rand <- function(rn,cn,rpl,pvalue=0.05)
{
  morans <- vector(mode="numeric",length=rpl)
  reject <- vector(mode="numeric",length=rpl)
  w <- nb2listw(cell2nb(rn,cn))
  n <- rn * cn
  for (i in 1:rpl){
    y <- rnorm(n)
    mi <- moran.test(y,w,alternative="two.sided")
    morans[i] <- mi$statistic
    if (mi$p.value) { reject[i] <- 1}
  }
  rejfreq <- sum(reject) / rpl
  value <- list(morans = morans, rej = rejfreq)
}
resu <- sim.moran.rand(5,5,1000)
resu$rej
pdf(file="moransim.pdf")
plot(resu$morans,main="Moran's I under Null - Randomization",xlab="Reference Distribution",lwd=2,col=2)
dev.off()


