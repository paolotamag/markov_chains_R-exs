#Paolo Tamagnini
#1536242
#paolotamag@gmail.com

#install.packages('msm', repos = 'http://cran.r-project.org')

library(msm)

set.seed(123)
x = c( 1.0,  1.5,  1.5,  1.5, 2.5, 4.0, 5.0, 5.0,  7.0, 8.0,  8.5,  9.0,  9.5, 9.5, 10.0, 12.0, 12.0, 13.0, 13.0, 14.5, 15.5, 15.5, 16.5, 17.0, 22.5, 29.0, 31.5)
y = c(1.80, 1.85, 1.87, 1.77, 2.02, 2.27, 2.15, 2.26, 2.47,2.19, 2.26, 2.40, 2.39, 2.41, 2.50, 2.32, 2.32, 2.43, 2.47, 2.56, 2.65, 2.47, 2.64, 2.56, 2.70, 2.72, 2.57)

n = length(y)

sigmaAlfa2 = 10000
sigmaBeta2 = 10000
a = 0.001
b = 0.001

nchain = 10000

alfaChainOrig<-rep(0, nchain)
betaChainOrig<-rep(0, nchain)
gammaChainOrig<-rep(0, nchain)
tau2ChainOrig<-rep(0, nchain)
Y20ChainOrig<-rep(0, nchain)
Y30ChainOrig<-rep(0, nchain)

alfa<-rnorm(1, mean = 0, sd = sqrt(sigmaAlfa2))

while(!is.finite(alfa) || alfa<1){
  alfa<-rnorm(1, mean = 0, sd = sqrt(sigmaAlfa2)) }

beta<-rnorm(1, mean = 0, sd = sqrt(sigmaBeta2))

while(!is.finite(beta) || beta<1){
  beta<-rnorm(1, mean = 0, sd = sqrt(sigmaBeta2)) }

gamma <- runif(1)

tau2<-1/rgamma(1,0.001, 0.001)

while(!is.finite(tau2)|| tau2<0){
  tau2<-1/rgamma(1,a, b) }

alfaMu <-function(gamma, tau2, beta) {
  num = sigmaAlfa2*sum(beta*gamma^x+y)
  den = (n*sigmaAlfa2+tau2)
  return(num/den) }

betaMu <-function(gamma, tau2, alfa) {
  num = sigmaBeta2*sum((alfa-y)*gamma^x)
  den = (sigmaBeta2*sum(gamma^(2*x))+tau2)
  return(num/den) }

alfaVar <-function(tau2) {
  num = tau2*sigmaAlfa2
  den = sigmaAlfa2*n+tau2
  return(num/den) }

betaVar <-function(tau2,gamma) {
  num = tau2*sigmaBeta2
  den = sigmaBeta2*sum(gamma^(2*x))+tau2
  return(num/den) }



metrWithGibbsUnif <-function(alfa,beta,gamma,tau2) {
  delta=10^(-2)
  h = runif(1,gamma-delta,gamma+delta)
  
  if (h<0 || h>1) {
    numh = 0 }
  
  numh = -sum((y-alfa+beta*h^(x))^2)
  numg = -sum((y-alfa+beta*gamma^(x))^2)
  den = 2*tau2
  
  fnum = exp(numh/den)
  fden = exp(numg/den)
  
  p = min(1,(fnum)/(fden))
  if(!is.finite(p)){
    p = 1  }
  if(runif(1)<p){
    return(h)
  }
  else {
    return(gamma)  }}


for (s in 1:nchain) {
  alfaChainOrig[s]<-alfa<-rtnorm(1, mean = alfaMu(gamma,tau2,beta), sd = sqrt(alfaVar(tau2)),lower = 1, upper = Inf)
  betaChainOrig[s]<-beta<-rtnorm(1, mean = betaMu(gamma,tau2,alfa), sd = sqrt(betaVar(tau2,gamma)),lower = 1, upper = Inf)
  gammaChainOrig[s]<-gamma<-metrWithGibbsUnif(alfa,beta,gamma,tau2)
  tau2ChainOrig[s]<-tau2<-1/rgamma(1,shape = a+n/2, rate=(b + sum((y-alfa+beta*gamma^(x))^2)/2))
  Y20ChainOrig[s]<-rnorm(1,alfa-beta*gamma^20,sqrt(tau2))
  Y30ChainOrig[s]<-rnorm(1,alfa-beta*gamma^30,sqrt(tau2))}

lb = 1001

alfaChain = alfaChainOrig[lb:10000]
betaChain = betaChainOrig[lb:10000]
gammaChain = gammaChainOrig[lb:10000]
tau2Chain = tau2ChainOrig[lb:10000]
Y20Chain = Y20ChainOrig[lb:10000]
Y30Chain = Y30ChainOrig[lb:10000]


alfa = mean(alfaChain)
beta = mean(betaChain)
gamma = mean(gammaChain)
tau2 = mean(tau2Chain)
Y20 = mean(Y20Chain)
Y30 = mean(Y30Chain)

####MLE#########
# alfa = 2.65  #
# beta = 0.96  #
# gamma = 0.87 #
# tau2 = 0.008 #
################

print(alfa)
print(beta)
print(gamma)
print(tau2)
print(Y20)
print(Y30)

lb = 0

alfaChain = alfaChainOrig[lb:10000]
gammaChain = gammaChainOrig[lb:10000]
Y20Chain = Y20ChainOrig[lb:10000]
Y30Chain = Y30ChainOrig[lb:10000]


itVec=seq(1,length(alfaChain))
plot(itVec,alfaChain,type = 'l',main='trace_plot_alpha',xlab = 't',ylab = 'alpha_t')
abline(h=alfa, col='red', lw = 2)


betaChain = betaChainOrig[3:10000]
itVecbeta=seq(1,length(betaChain))

plot(itVecbeta,betaChain,type = 'l',main='trace_plot_beta',xlab = 't',ylab = 'beta_t')
abline(h=beta, col='red', lw = 2)

plot(itVec,gammaChain,type = 'l',main='trace_plot_gamma',xlab = 't',ylab = 'gamma_t')
abline(h=gamma, col='red', lw = 2)


tau2Chain = tau2ChainOrig[4:10000]
itVectau=seq(1,length(tau2Chain))

plot(itVectau,tau2Chain,type = 'l',main='trace_plot_tau2',xlab = 't',ylab = 'tau2_t')
abline(h=tau2, col='red', lw = 2)

plot(itVec,Y20Chain,type = 'l',main='trace_plot_Y20',xlab = 't',ylab = 'Y20_t')
abline(h=Y20, col='red', lw = 2)

plot(itVec,Y30Chain,type = 'l',main='trace_plot_Y30',xlab = 't',ylab = 'Y30_t')
abline(h=Y30, col='red', lw = 2)



lb = 1001

alfaChain = alfaChainOrig[lb:10000]
betaChain = betaChainOrig[lb:10000]
gammaChain = gammaChainOrig[lb:10000]
tau2Chain = tau2ChainOrig[lb:10000]
Y20Chain = Y20ChainOrig[lb:10000]
Y30Chain = Y30ChainOrig[lb:10000]

#fancy:
x = alfaChain
iter = (lb):(length(x)+(lb-1))
runningmeans=cumsum(x)/(1:length(x))
plot(iter,runningmeans,type="l",col='red', main='run_mean_alfa', xlab = 't',ylab = 'I_t',ylim =  c(2.55,2.7))
abline(h=runningmeans[length(runningmeans)], col='blue')
abline(h=alfa, col='green')

x = betaChain
runningmeans=cumsum(x)/(1:length(x))
plot(iter,runningmeans,type="l",col='red', main='run_mean_beta', xlab = 't',ylab = 'I_t', ylim = c(1.02,1.06))
abline(h=runningmeans[length(runningmeans)], col='blue')
abline(h=beta, col='green')

x = gammaChain
runningmeans=cumsum(x)/(1:length(x))
plot(iter,runningmeans,type="l",col='red', main='run_mean_gamma', xlab = 't',ylab = 'I_t', ylim = c(0.8,0.9))
abline(h=runningmeans[length(runningmeans)], col='blue')
abline(h=gamma, col='green')

x = tau2Chain
runningmeans=cumsum(x)/(1:length(x))
plot(iter,runningmeans,type="l",col='red', main='run_mean_tau2', xlab = 't',ylab = 'I_t')
abline(h=runningmeans[length(runningmeans)], col='blue')
abline(h=tau2, col='green')

x = Y20Chain
runningmeans=cumsum(x)/(1:length(x))
plot(iter,runningmeans,type="l",col='red', main='run_mean_Y20', xlab = 't',ylab = 'I_t')
abline(h=runningmeans[length(runningmeans)], col='blue')
abline(h=Y20, col='green')

x = Y30Chain
runningmeans=cumsum(x)/(1:length(x))
plot(iter,runningmeans,type="l",col='red', main='run_mean_Y30', xlab = 't',ylab = 'I_t')
abline(h=runningmeans[length(runningmeans)], col='blue')
abline(h=Y30, col='green')


lb = 1001

alfaChain = alfaChainOrig[lb:10000]
betaChain = betaChainOrig[lb:10000]
gammaChain = gammaChainOrig[lb:10000]
tau2Chain = tau2ChainOrig[lb:10000]
Y20Chain = Y20ChainOrig[lb:10000]
Y30Chain = Y30ChainOrig[lb:10000]

gamma_k <- function(k,vet,Icap) {
  t = length(vet)
  somma = 0
  for (ki in 1:(t-k)) {
    somma = somma + (vet[ki] - Icap)*(vet[ki+k] - Icap) }
  return(somma/(t-k))}


varCap <- function(vet,Icap) {
  somma = 0
  tVar = length(vet)
  print('Summing up all gamma_k, from k = 1 to:')
  print(tVar)
  print('Starting now..')
  print(1)
  for (klol in 1:(tVar-1)) {
    if (klol%%500==0) {
      print(klol) }
    somma = somma + gamma_k(klol,vet,Icap) }
  num = gamma_k(0,vet,Icap) + 2*somma
  return(num/tVar) }

errAlfa = varCap(alfaChain,alfa)

errBeta = varCap(betaChain,beta)

errGamma = varCap(gammaChain,gamma)

errTau2 = varCap(tau2Chain,tau2)

errY20 = varCap(Y20Chain,Y20)

errY30 = varCap(Y30Chain,Y30)

print(alfa)
print(errAlfa)
print(beta)
print(errBeta)
print(gamma)
print(errGamma)
print(tau2)
print(errTau2)
print(Y20)
print(errY20)
print(Y30)
print(errY30)

precY20 = 1 / abs(errY20)
precY30 = 1 / abs(errY30)

print(precY20)
print(precY30)

coeffAlfa = sqrt(abs(errAlfa))/alfa
coeffBeta = sqrt(abs(errBeta))/beta
coeffGamma = sqrt(abs(errGamma))/gamma
coeffTau2 = sqrt(abs(errTau2))/tau2

coeffAlfa
coeffBeta
coeffGamma
coeffTau2

#not required
acf(alfaChain, type = "correlation")
acf(betaChain, type = "correlation")
acf(gammaChain, type = "correlation")
acf(tau2Chain, type = "correlation")


cor_A_B = cor(alfaChain,betaChain)
cor_A_G = cor(alfaChain,gammaChain)
cor_A_T = cor(alfaChain,tau2Chain)
cor_B_G = cor(betaChain,gammaChain)
cor_B_T = cor(betaChain,tau2Chain)
cor_G_T = cor(gammaChain,tau2Chain)

listCor <- list("cor_Alfa_Beta"= abs(cor_A_B), "cor_Alfa_Gamma"= abs(cor_A_G), "cor_Alfa_Tau2"= abs(cor_A_T), 
                "cor_Beta_Gamma"= abs(cor_B_G),"cor_Beta_Tau2"= abs(cor_B_T),"cor_Gamma_Tau2"= abs(cor_G_T))

listCor
print('the best corr has the symbol V')
sapply(listCor, function(X) if(X != max(unlist(listCor))) print('X') else print('V'))

