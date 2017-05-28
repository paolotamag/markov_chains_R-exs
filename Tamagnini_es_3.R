#Paolo Tamagnini
#1536242
#paolotamag@gmail.com

seme = 123
xzero = 1
#a
set.seed(seme)
mpt<-matrix(c(0,1/2,1/2,5/8,1/8,1/4,2/3,1/3,0),nrow=3,byrow=T)
S=c(1,2,3) 
x0<-xzero #
nsample<-1000
chain<-rep(NA,nsample+1) 
chain[1]<-x0 
for(t in 1:nsample){
  chain[t+1]<-sample(S,size=1,prob=mpt[chain[t],]) }
plot(chain,ylim=c(0,4))

#b
resultB = table(chain)/nsample
print(resultB)


#c
set.seed(seme)
mpt<-matrix(c(0,0.5,0.5,0.625,0.125,0.25,2/3,1/3,0),nrow=3,byrow=T)
S=c(1,2,3) 
x0<-xzero
nsample<-1000
nchains = 500
lastVal<-rep(NA,nchains)
for (k in 1:nchains) {
  chain<-rep(NA,nsample+1) 
  chain[1]<-x0 
  for(t in 1:nsample){
    chain[t+1]<-sample(S,size=1,prob=mpt[chain[t],]) }
  lastVal[k] = chain[nsample+1]
}

plot(lastVal,ylim=c(0,4))
resultC = table(lastVal)/nchains
print(resultC)

#d
M<-matrix(c(-1,5/8,2/3,1/2,-7/8,1/3,1/2,1/4,-1,1,1,1),nrow=4,byrow=T)
c = c(0,0,0,1)
resultPi = qr.solve(M,c)
print(resultPi)

resultPi - resultB
resultPi - resultC

mean(abs(resultPi - resultB))
mean(abs(resultPi - resultC))
