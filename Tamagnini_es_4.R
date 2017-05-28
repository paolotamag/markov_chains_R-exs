#Paolo Tamagnini
#1536242
#paolotamag@gmail.com

#install.packages('plot3D', repos = 'http://cran.r-project.org')
#install.packages('rgl', repos = 'http://cran.r-project.org')

library(plot3D)
library(rgl)

set.seed(123)
Y=c(4,5,4,1,0,4,3,4,0,6,
    + 3,3,4,0,2,6,3,3,5,4,5,3,1,4,4,1,5,5,3,4,2,5,2,2,3,4,2,1,3,2,
    + 1,1,1,1,1,3,0,0,1,0,1,1,0,0,3,1,0,3,2,2,0,1,1,1,0,1,0,1,0,0,
    + 0,2,1,0,0,0,1,1,0,2,2,3,1,1,2,1,1,1,1,2,4,2,0,0,0,1,4,0,0,0,
    + 1,0,0,0,0,0,1,0,0,1,0,0)
n = length(Y)
alfa<-2
beta<-1
a<-2
b<-1
nchain<-11000

post <- function(theta){
  lamb = theta[1]
  ph = theta[2]
  mm = theta[3]
  prod = exp(-mm*lamb)
  prod = exp(-(n-mm)*ph)*prod
  prod = lamb^(sum(Y[1:mm]))*prod
  prod = phi^(sum(Y[(mm+1):n]))*prod
  prod = exp(-beta*lamb)*prod
  prod = lamb^(a-1)*prod
  prod = exp(-beta*ph)*prod
  prod = ph^(a-1)*prod
  return(prod)}

mChainDist<-function(t1, t2, X) {
  n<-length(X)
  #print(n)
  probVector<-rep(0,n-1)
  for (i in 1:n-1) {
    uno<-exp((t2-t1)*i)
    due<-t1^sum(X[1:i])
    tre<-t2^sum(X[(i+1):n]) 
    probVector[i]<-uno*due*tre }
  probVector<-probVector/sum(probVector)
  return(probVector) }

lambdaChain<-rep(0, nchain)
phiChain<-rep(0, nchain)
mChain<-rep(0, nchain)
lambda<-rgamma(1, alfa, beta)
phi<-rgamma(1, a, b)
m<-sample(1:n, 1)


for (s in 1:nchain) {
  lambdaChain[s]<-lambda<-rgamma(1, alfa+sum(Y[1:m]), beta+m)
  phiChain[s]<-phi<-rgamma(1, a+sum(Y[(m+1):n]), b+n-m)
  probVect<-mChainDist(lambda, phi, Y)
  mChain[s]<-m<-sample(1:(n-1), 1, prob=probVect) }

lambdaMean = mean(lambdaChain[1001:nchain])
phiMean = mean(phiChain[1001:nchain])
m = mean(mChain[1001:nchain])
t = table(mChain)/length(mChain)
moda = strtoi(names(t[which(t == max(t))[1]]))
mYearsMode = moda + 1850
mYearsMean = m + 1850
postValue = post(c(lambdaMean,phiMean,moda))

lambdaMean
phiMean
moda
mYearsMode
m
mYearsMean
postValue

lamdaSim = lambdaChain[1001:nchain]

hist(lamdaSim,freq = F,breaks = 75, col = 'orange', main = 'simulation vs full-conditional',xlab='lambda')
xfit<-seq(min(lamdaSim),max(lamdaSim),length=length(lamdaSim)) 
yfit<-dgamma(xfit,alfa+sum(Y[1:moda]),beta+moda) 
lines(xfit, yfit, col="blue", lwd=2)

phiSim = phiChain[1001:nchain]

hist(phiSim,freq = F,breaks = 75, col = 'orange', main = 'simulation vs full-conditional',xlab='phi')
xfit<-seq(min(phiSim),max(phiSim),length=length(phiSim)) 
yfit<-dgamma(xfit,a+sum(Y[(moda+1):n]),b+n-moda) 
lines(xfit, yfit, col="blue", lwd=2)

mSim = mChain[1001:nchain]

tavol = table(mSim)/length(mSim)
tbl <- barplot(tavol,col = 'orange', main = 'simulation vs full-conditional',ylim=c(0,0.30),xlab='m')
listM = c(30,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51)
yfit<-mChainDist(lambdaMean,phiMean,Y)[listM] 
lines(x = tbl, y = yfit, col="blue",lwd=2)
points(x = tbl, y = yfit,lwd=2)


#3d histograms and 2d maps
definition = 50

x <- lamdaSim
y <- phiSim
x_c <- cut(x, definition)
y_c <- cut(y, definition)
z <- table(x_c, y_c)

hist3D(z=z, border="black",xlab='lambda',ylab='phi',zlab='frequency',ticktype = "detailed")

image2D(z=z, border="black",xlab='lambda',ylab='phi')


x <- lamdaSim
y <- mSim
x_c <- cut(x, definition)
y_c <- cut(y, definition)
z <- table(x_c, y_c)

hist3D(z=z,border="black", xlab='lambda',ylab='m',zlab='frequency',ticktype = "detailed")

image2D(z=z, border="black",xlab='lambda',ylab='m')


x <- phiSim
y <- mSim
x_c <- cut(x, definition)
y_c <- cut(y, definition)
z <- table(x_c, y_c)

hist3D(z=z, border="black",xlab='phi',ylab='m',zlab='frequency',ticktype = "detailed")

image2D(z=z, border="black",xlab='phi',ylab='m')


#3d scatterplots

nccz= length(lamdaSim)
loll = nccz - 100
xp <- lamdaSim[loll:nccz]
yp <- phiSim[loll:nccz]
np = length(xp)
Zp<-rep(0, np^2)
Xp<-rep(0, np^2)
Yp<-rep(0, np^2)
zi = 1
cc = 1
for (ip in 1:np) {
  for (jp in 1:np) {
    Zp[zi]=post(c(xp[ip],yp[jp],moda))
    Xp[zi]=xp[ip]
    Yp[zi]=yp[jp] 
    zi = zi + 1}
  cc = cc + 1}
plot3d(x=Xp, y=Yp, z=Zp, type="p", col="red", xlab="lambda", ylab="phi", zlab="post", 
       size=5, lwd=15, box=F)

xp <- lamdaSim[loll:nccz]
yp <- mChain[1001:nchain][loll:nccz]
np = length(xp)
Zp<-rep(0, np^2)
Xp<-rep(0, np^2)
Yp<-rep(0, np^2)
zi = 1
cc = 1
for (ip in 1:np) {
  for (jp in 1:np) {
    Zp[zi]=post(c(xp[ip],phiMean,yp[jp]))
    Xp[zi]=xp[ip]
    Yp[zi]=yp[jp] 
    zi = zi + 1}
  cc = cc + 1}
plot3d(x=Xp, y=Yp, z=Zp, type="p", col="red", xlab="lambda", ylab="m", zlab="post", 
       size=5, lwd=15, box=F)

xp <- phiSim[loll:nccz]
yp <- mChain[1001:nchain][loll:nccz]
np = length(xp)
Zp<-rep(0, np^2)
Xp<-rep(0, np^2)
Yp<-rep(0, np^2)
zi = 1
cc = 1
for (ip in 1:np) {
  for (jp in 1:np) {
    Zp[zi]=post(c(lambdaMean,xp[ip],yp[jp]))
    Xp[zi]=xp[ip]
    Yp[zi]=yp[jp] 
    zi = zi + 1}
  cc = cc + 1}
plot3d(x=Xp, y=Yp, z=Zp, type="p", col="red", xlab="phi", ylab="m", zlab="post", 
       size=5, lwd=15, box=F)

percChange = (phiMean - lambdaMean)/lambdaMean*100
percChange
