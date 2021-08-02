setwd("C:/Users/rl02898/Documents/Bayesian Statistics/Final Project")

##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################

data1 <- read.csv("Drug_Overdose_Deaths.csv",header = TRUE,sep = ",")
a <- data1[data1$Year=="2018" & data1$Indicator=="Number of Drug Overdose Deaths",]
b <- as.matrix(a)
c <- b[,-c(2:5,7:12)]
d <- transform(c, Data.Value = as.numeric(Data.Value))
e <- aggregate(Data.Value~State,d,sum)
f <- e[-c(45,53),]           ## Finally done sorting the large data set
Populations <- c(737438,4887871,3013825,7171646,39865590,5695564,
                 3572665,702455,967171,21299325,10519475,1420491,
                 3156145,1754208,12741080,6691878,2911505,4468402,
                 4659978,6902149,6042718,1338404,9995915,5611179,
                 6126452,2986530,1062305,10383620,760077,1929268,
                 1356458,8908520,2095428,3034392,19542209,11689442,
                 3943079,4190713,12807060,1057315,5084127,882235,
                 6770010,29865997,3161105,8517685,626299,7535591,
                 5813568,1805832,577737)
Region <- c(1,5,5,3,2,3,7,6,6,5,5,2,4,1,4,4,4,4,5,7,6,7,4,4,4,
            5,1,5,4,4,7,6,3,2,6,4,3,1,6,7,5,4,5,3,3,6,7,1,4,6,1)
Income <- c(73181,48123,45869,56581,71805,69117,74168,82372,
             62852,52594,56183,77765,58570,52225,62992,54181,
             56422,48375,46145,77385,80776,55277,54909,68388,
             53578,43529,53386,52752,61843,59970,73381,80088,
             46744,58003,64894,54021,50051,60212,59105,63870,
             50570,56521,51340,59206,65358,71535,57513,70979,
             59305,43469,60434)
Health.Care <- c(11064,7281,7408,6452,7549,6804,9859,11944,10254,
                 8076,6587,7299,8200,6927,8262,8300,7651,8004,7815,
                 10559,8602,9531,8055,8871,8107,7646,8221,7264,9851,
                 8412,9589,8859,7214,6714,9778,8712,7627,8044,9258,
                 9551,7311,8933,7372,6998,5982,7556,10190,7913,8702,
                 9462,8320
                 )
Percent.Pop <- (f[,2]/Populations)*100
data2 <- cbind(f,Percent.Pop,Region,Income,Health.Care)
row.names(data2) <- 1:nrow(data2)
colnames(data2) <- c("State","Deaths","Percent.Pop","Region","Income","Health.Care")
data <- data2[,3]


## Let's first look at a frequentist approach to get a feel for the data.
linreg <- lm(data2[,3]~data2[,4]+data2[,5]+data2[,6])
summary(linreg)
linreg1 <- lm(data2[,3]~data2[,4])
summary(linreg1)

pdf(file = "linreg1.pdf",width = 9,height=9)
par(mfrow=c(2,2))
plot(data2[,2],data2[,3])
abline(0.06978,-0.000003112)
plot(data2[,3],data2[,4])
plot(data2[,3],data2[,5])
plot(data2[,3],data2[,6])
dev.off()

posterior.x <- rnorm(1000,(51/0.03*0.058+1/sd(data)*mean(data))/(51/0.03+1/sd(data)),sqrt(1/(51/0.03+1/sd(data))))
posterior.y <- dnorm(posterior.x,(51/0.03*0.058+1/sd(data)*mean(data))/(51/0.03+1/sd(data)),sqrt(1/(51/0.03+1/sd(data))))
pdf(file = "posteriors.pdf",width = 9,height=9)
par(mfrow=c(1,1))
plot(sort(posterior.x),posterior.y[order(posterior.x)],type = "l",col="red")
prior.x <- rnorm(1000,0.058,0.03)
prior.y <- dnorm(prior.x,0.058,0.03)
lines(sort(prior.x),prior.y[order(prior.x)],col="blue")
likelihood.x <- rnorm(1000,mean(data),sd(data))
likelihood.y <- dnorm(likelihood.x,mean(data),sd(data))
lines(sort(likelihood.x),likelihood.y[order(likelihood.x)])
lines(density(data),lty=2)
legend("topright",inset = .03,legend = c("Likelihood","Prior","Posterior","density()"), col = c("black","blue","red","black"),lty = c(1,1,1,2),cex = .8)
abline(v=0.058,col="blue")
abline(v=mean(data))
abline(v=(51/0.03*0.058+1/sd(data)*mean(data))/(51/0.03+1/sd(data)),col="red")
abline(v=0.051,lty=2)
dev.off()



##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################

## Let's see if Monte Carlo Sampling supports our results so far.
## Establish initials
mu <- 0.05
sigma <- 1
delta <- 1
prior.mean <- 0.058
prior.sd <- 0.03
theta.initial <- rnorm(1,mu,1)
N <- 20000
theta.post <- rep(NA, N)
theta.post[1] <- theta.initial

## Metropolis-Hastings Algorithm
for (i in 2:N) {
  theta = rnorm(1,theta.post[i-1],delta)
  num = prod(dnorm(data,theta,sigma))*dnorm(theta,prior.mean,prior.sd)*dnorm(theta.post[i-1],theta,delta)
  denom = prod(dnorm(data,theta.post[i-1],sigma))*dnorm(theta.post[i-1],prior.mean,prior.sd)*dnorm(theta,theta.post[i-1],delta)
  W = num/denom
  alpha = min(1,W)
  U = runif(1,0,1)
  
  if (U <= alpha) {
    theta.post[i] = theta
    print("advance")
  }
  else {
    theta.post[i] <- theta.post[i-1]
    print("maintain")
  }
  print(theta.post[i])
}

## Burn-in
burn <- 7500
post.sample <- theta.post[(burn+1):N]

## Posterior Density vs Algorithm Density
pdf(file = "metropolis.pdf",width = 9,height=9)
plot(sort(posterior.x),posterior.y[order(posterior.x)],type = "l",col="red")
lines(density(post.sample))
dev.off()

## Posterior Summary
mean(post.sample)
sd(post.sample)

## History Plot
pdf(file = "history.pdf",width = 9,height=9)
par(mfrow=c(2,1))
plot(c(1:N), theta.post, type="l")                     # before burn-in
plot(c(1:length(post.sample)), post.sample, type="l")  # after burn-in
dev.off()



##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################

## Now let's construct our posterior via OpenBUGS.
## Create the model
#install.packages("R2OpenBUGS")
library(R2OpenBUGS)
model1 <- function(){
  for (i in 1:51){
    y[i] ~ dnorm(mu,29.33735)
  }
  
  mu ~ dnorm(0.058,33.3333333)
  y.pred ~ dnorm(mu,29.33735)
  prob.y <- step(y.pred-0.1)
}
model.drug1 <- file.path("C:/Users/rl02898/Documents/Bayesian Statistics/Final Project/model.drug1.txt")
write.model(model1, model.drug1)

data.drug1 <- list("y"=data)
inits.drug1 <- function(){list(mu=0.5)}
params.drug1 <- c("mu","y.pred","prob.y")
out.drug1 <- bugs(data=data.drug1, inits=inits.drug1, parameters.to.save=params.drug1,
                    model.file=model.drug1,n.iter=20000, n.burnin=5000, debug=F,n.chains=1)
samples1.ypred <- out.drug1$sims.list$y.pred
samples1.proby <- out.drug1$sims.list$prob.y
samples1.mu <- out.drug1$sims.list$mu
pdf(file = "openbugs.pdf",width = 9,height=9)
par(mfrow=c(1,1))
plot(sort(posterior.x),posterior.y[order(posterior.x)],type = "l",col="red")
lines(density(samples1.mu))   ## The posteriors compare very well, they should.
dev.off()

## 95% Credible Interval for mu
M <- 12500
credible.interval <- function(sample,alpha){
  c(sort(sample)[length(sample)*(alpha/2)],sort(sample)[length(sample)*(1-(alpha/2))])
}

credible.interval(samples1.mu,0.05)


## Highest Posterior Density (HPD Interval)
#install.packages("HDInterval")
library(HDInterval)
hdi(samples1.mu)


## Function for the estimated CDF(a,b)
integrate.sample <- function(sample,a,b){
  if (missing(b)){
    sum(sample < a)/length(sample)
  }
  else {
    if (missing(a)) {
      sum(sample > b)/length(sample)
    }
    else {
      sum(sample < b)/length(sample) - sum(sample < a)/length(sample)
    }
  }
}
# If a is left blank, it will evalute upper tail area.
# If b is left blank, it will evaulate lower tail area.
integrate.sample(samples1.mu,0.00,0.05)

# Returns the approximate HPD interval for any given sample and alpha, for unimodal posterior samples.
interval.HPD <- function(sample,alpha){
  
  pdf <- density(sample,n=1000)
  a <- 0
  q <- (which(pdf$y==max(pdf$y))-5)
  Q <- length(pdf$y)
  
  for (i in q:1){
    for (j in (i+1):Q){
      if (pdf$y[j] <= pdf$y[i]){
        a <- integrate.sample(sample,pdf$x[i],pdf$x[j])
      }
      if (a > (1-alpha)){
        break
      }
      if (pdf$y[j] <= pdf$y[i]){
        break
      }
    }
    if (a > (1-alpha)){
      print("The HDP Interval is:")
      print(pdf$x[i])
      print(pdf$x[j])
      pdf(file = "hdi.pdf",width = 9,height=9)
      plot <- plot(density(sample))
      x <- pdf$x[i:j]
      y <- pdf$y[i:j]
      polygon(x,y,density=50,col='turquoise1')
      z <- c(pdf$x[i],pdf$x[i],pdf$x[j],pdf$x[j])
      w <- c(0,pdf$y[i],pdf$y[j],0)
      polygon(z,w,density=50,col='turquoise1')
      abline(h=pdf$y[i],col="gray34")
      abline(v=pdf$x[i],col="gray34")
      abline(v=pdf$x[j],col="gray34")
      dev.off()
      break
    }
  }
}

interval.HPD(samples1.mu,0.05)
hdi(samples1.mu)                    ## These intervals are very close



##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################

## Now let's run a linear regression on the data.
## Let's start by using a shrinking prior, to find important factors.
## Zellner's g-prior
Y <- data2[,3]
X <- matrix(c(rep(1,length=length(Y)),data2[,4],data2[,5],data2[,6]),ncol = 4)
colnames(X) <- c("Intercept","X1","X2","X3")
taumatrix <- t(X) * X
mean <- c(0,0,0,0)
m <- length(Y)

model2 <- function(){
  for (i in 1:m){
    Y[i] ~ dnorm(mu[i],tau)
    mu[i] <- inprod(beta[],X[i,])
  }
  
  g <- 2 * m
  beta[1:4] ~ dmnorm(mean[],prec[,])
  
  for(i in 1:4){
    for(j in 1:4){
      prec[i,j] <- (1/g)*tau*taumatrix[i,j]
    }
  }
  
  for (i in 1:m) {
    y.pred[i] ~ dnorm(mu[i],tau)
  }
  
  for (i in 1:m) {
    log_cpo[i] <- -0.5*log(tau/6.283)-0.5*pow((Y[i]-mu[i]),2)
  }
  LPML <- sum(log_cpo[1:m])
  
  for (i in 1:m) {
    y[i] ~ dnorm(mu[i],tau)
  }
  
  tau ~ dgamma(0.01,0.01)
}
model.drug2 <- file.path("C:/Users/rl02898/Documents/Bayesian Statistics/Final Project/model.drug2.txt")
write.model(model2, model.drug2)

data.drug2 <- list("Y","X","m","taumatrix","mean")
inits.drug2 <- function(){list(mu=0.5)}
params.drug2 <- c("mu","beta","LPML")
out.drug2 <- bugs(data=data.drug2, inits=inits.drug2, parameters.to.save=params.drug2,
                  model.file=model.drug2,n.iter=20000, n.burnin=5000, debug=F,n.chains=1)
samples2.mu <- out.drug2$sims.list$mu
samples2.beta1 <- out.drug2$sims.list$beta[,1]
samples2.beta2 <- out.drug2$sims.list$beta[,2]
samples2.beta3 <- out.drug2$sims.list$beta[,3]
samples2.beta4 <- out.drug2$sims.list$beta[,4]
samples2.LPML <- out.drug2$sims.list$LPML
mean(samples2.mu)
mean(samples2.beta1)
mean(samples2.beta2)
mean(samples2.beta3)
mean(samples2.beta4)

## How well can this model predict our death rates?
## Let's calculate the residuals.
residual2 <- rep(NA,length(Y))

for (i in 1:length(Y)) {
  y <- mean(samples2.beta1)+mean(samples2.beta2)*data2[i,4]+mean(samples2.beta3)*data2[i,5]+mean(samples2.beta4)*data2[i,6]
  residual2[i] <- Y[i]-y
}
max(abs(residual2))
min(abs(residual2))      ## Second lowest min residual


## Now let's try using a double exponential prior
model3 <- function(){
  for (i in 1:m){
    Y[i] ~ dnorm(mu[i],tau)
    mu[i] <- inprod(beta[],X[i,])
  }
  
  for (i in 1:4) {
    beta[i] ~ ddexp(0,tau)
  }
  
  for (i in 1:m) {
    log_cpo[i] <- -0.5*log(tau/6.283)-0.5*pow((Y[i]-mu[i]),2)
  }
  LPML <- sum(log_cpo[1:m])
  
  tau ~ dnorm(0.03,1)       ## Model was not converging with  tau ~ gamma(0.01,0.01)
}
model.drug3 <- file.path("C:/Users/rl02898/Documents/Bayesian Statistics/Final Project/model.drug3.txt")
write.model(model3, model.drug3)

data.drug3 <- list("Y","X","m")
inits.drug3 <- function(){list(mu=0.05,tau=0.03)}
params.drug3 <- c("mu","beta","LPML")
out.drug3 <- bugs(data=data.drug3, inits=inits.drug3, parameters.to.save=params.drug3,
                  model.file=model.drug3,n.iter=20000, n.burnin=5000, debug=F,n.chains=1)
samples3.mu <- out.drug3$sims.list$mu
samples3.beta1 <- out.drug3$sims.list$beta[,1]
samples3.beta2 <- out.drug3$sims.list$beta[,2]
samples3.beta3 <- out.drug3$sims.list$beta[,3]
samples3.beta4 <- out.drug3$sims.list$beta[,4]
samples3.LPML <- out.drug3$sims.list$LPML
mean(samples3.mu)
mean(samples3.beta1)
mean(samples3.beta2)
mean(samples3.beta3)
mean(samples3.beta4)
residual3 <- rep(NA,length(Y))

for (i in 1:length(Y)) {
  y <- mean(samples3.beta1)+mean(samples3.beta2)*data2[i,4]+mean(samples3.beta3)*data2[i,5]+mean(samples3.beta4)*data2[i,6]
  residual3[i] <- Y[i]-y
}
max(abs(residual3))
min(abs(residual3))


## The double exponential did not seem like a good fit at all. Let's try a flat normal prior to see where it takes us.
model4 <- function(){
  for (i in 1:m){
    Y[i] ~ dnorm(mu[i],tau)
    mu[i] <- inprod(beta[],X[i,])
  }
  
  for (i in 1:4) {
    beta[i] ~ dnorm(0,1e-6)
  }
  
  for (i in 1:m) {
    log_cpo[i] <- -0.5*log(tau/6.283)-0.5*pow((Y[i]-mu[i]),2)
  }
  LPML <- sum(log_cpo[1:m])
  
  tau ~ dgamma(0.01,0.01)
  
  for (i in 1:m) {
    y[i] ~ dnorm(mu[i],tau)
  }
  
}
model.drug4 <- file.path("C:/Users/rl02898/Documents/Bayesian Statistics/Final Project/model.drug4.txt")
write.model(model4, model.drug4)

data.drug4 <- list("Y","X","m")
inits.drug4 <- function(){list(mu=0.05,tau=0.03)}
params.drug4 <- c("mu","beta","LPML")
out.drug4 <- bugs(data=data.drug4, inits=inits.drug4, parameters.to.save=params.drug4,
                  model.file=model.drug4,n.iter=20000, n.burnin=5000, debug=F,n.chains=1)
samples4.mu <- out.drug4$sims.list$mu
samples4.beta1 <- out.drug4$sims.list$beta[,1]
samples4.beta2 <- out.drug4$sims.list$beta[,2]
samples4.beta3 <- out.drug4$sims.list$beta[,3]
samples4.beta4 <- out.drug4$sims.list$beta[,4]
samples4.LPML <- out.drug4$sims.list$LPML
mean(samples4.mu)
mean(samples4.beta1)
mean(samples4.beta2)
mean(samples4.beta3)
mean(samples4.beta4)
residual4 <- rep(NA,length(Y))

for (i in 1:length(Y)) {
  y <- mean(samples4.beta1)+mean(samples4.beta2)*data2[i,4]+mean(samples4.beta3)*data2[i,5]+mean(samples4.beta4)*data2[i,6]
  residual4[i] <- Y[i]-y
}
max(abs(residual4))     ## Lowest max residual
min(abs(residual4))



## Let's try this again by eliminating the last two factors to make this simple linear regression.
## From the frequentist approach, this is what we were predicting to happen.
## We will test whether the region you live in affects the percentage of overdose.
Y <- data2[,3]
X1 <- matrix(c(rep(1,length=length(Y)),data2[,4]),ncol = 2)
colnames(X1) <- c("Intercept","X1")
taumatrix1 <- t(X1) %*% X1
mean1 <- c(0,0)
m <- length(Y)

model5 <- function(){
  for (i in 1:m){
    Y[i] ~ dnorm(mu[i],tau)
    mu[i] <- inprod(beta[],X1[i,])
  }
  
  g <- 2 * m
  beta[1:2] ~ dmnorm(mean1[],prec[,])
  
  for(i in 1:2){
    for(j in 1:2){
      prec[i,j] <- (1/g)*tau*taumatrix1[i,j]
    }
  }
  
  for (i in 1:m) {
    log_cpo[i] <- -0.5*log(tau/6.283)-0.5*pow((Y[i]-mu[i]),2)
  }
  LPML <- sum(log_cpo[1:m])
  
  tau ~ dgamma(0.01,0.01)
  
  for (i in 1:m) {
    y[i] ~ dnorm(mu[i],tau)
  }
}
model.drug5 <- file.path("C:/Users/rl02898/Documents/Bayesian Statistics/Final Project/model.drug5.txt")
write.model(model5, model.drug5)

data.drug5 <- list("Y","X1","m","taumatrix1","mean1")
inits.drug5 <- function(){list(mu=0.5)}
params.drug5 <- c("mu","beta","LPML")
out.drug5 <- bugs(data=data.drug5, inits=inits.drug5, parameters.to.save=params.drug5,
                  model.file=model.drug5,n.iter=20000, n.burnin=5000, debug=F,n.chains=1)
samples5.mu <- out.drug5$sims.list$mu
samples5.beta1 <- out.drug5$sims.list$beta[,1]
samples5.beta2 <- out.drug5$sims.list$beta[,2]
samples5.LPML <- out.drug5$sims.list$LPML
mean(samples5.mu)
mean(samples5.beta1)
mean(samples5.beta2)
pdf(file = "regression.pdf",width = 9,height=9)
plot(data2[,4],data2[,3],type = "p")
abline(mean(samples5.beta1),mean(samples5.beta2))     ## There is clearly an upward trend depending on where you live.
dev.off()
residual5 <- rep(NA,length(Y))

for (i in 1:length(Y)) {
  y <- mean(samples5.beta1)+mean(samples5.beta2)*data2[i,4]
  residual5[i] <- Y[i]-y
}
max(abs(residual5))
min(abs(residual5))


## Now let's try using a double exponential prior
model6 <- function(){
  for (i in 1:m){
    Y[i] ~ dnorm(mu[i],tau)
    mu[i] <- inprod(beta[],X1[i,])
  }
  
  for (i in 1:2) {
    beta[i] ~ ddexp(0,tau)
  }
  
  for (i in 1:m) {
    log_cpo[i] <- -0.5*log(tau/6.283)-0.5*pow((Y[i]-mu[i]),2)
  }
  LPML <- sum(log_cpo[1:m])
  
  tau ~ dnorm(0.03,1)       ## Model was not converging with  tau ~ gamma(0.01,0.01)
}
model.drug6 <- file.path("C:/Users/rl02898/Documents/Bayesian Statistics/Final Project/model.drug6.txt")
write.model(model6, model.drug6)

data.drug6 <- list("Y","X1","m")
inits.drug6 <- function(){list(mu=0.05,tau=0.03)}
params.drug6 <- c("mu","beta","LPML")
out.drug6 <- bugs(data=data.drug6, inits=inits.drug6, parameters.to.save=params.drug6,
                  model.file=model.drug6,n.iter=20000, n.burnin=5000, debug=F,n.chains=1)
samples6.mu <- out.drug6$sims.list$mu
samples6.beta1 <- out.drug6$sims.list$beta[,1]
samples6.beta2 <- out.drug6$sims.list$beta[,2]
samples6.LPML <- out.drug6$sims.list$LPML
mean(samples6.mu)
mean(samples6.beta1)
mean(samples6.beta2)
pdf(file = "regression1.pdf",width = 9,height=9)
plot(data2[,4],data2[,3],type = "p")
abline(mean(samples6.beta1),mean(samples6.beta2))     ## There is clearly an upward trend depending on where you live.
dev.off()
residual6 <- rep(NA,length(Y))

for (i in 1:length(Y)) {
  y <- mean(samples6.beta1)+mean(samples6.beta2)*data2[i,4]
  residual6[i] <- Y[i]-y
}
max(abs(residual6))
min(abs(residual6))


## Let's try a flat normal prior to see where it takes us.
model7 <- function(){
  for (i in 1:m){
    Y[i] ~ dnorm(mu[i],tau)
    mu[i] <- inprod(beta[],X1[i,])
  }
  
  for (i in 1:2) {
    beta[i] ~ dnorm(0,1e-6)
  }
  
  for (i in 1:m) {
    log_cpo[i] <- -0.5*log(tau/6.283)-0.5*pow((Y[i]-mu[i]),2)
  }
  LPML <- sum(log_cpo[1:m])
  
  tau ~ dgamma(0.01,0.01)
  
  for (i in 1:m) {
    y[i] ~ dnorm(mu[i],tau)
  }
  
  for (i in 1:m) {
    y.pred[i] ~ dnorm(mu[i],tau)
  }
  
}
model.drug7 <- file.path("C:/Users/rl02898/Documents/Bayesian Statistics/Final Project/model.drug7.txt")
write.model(model7, model.drug7)

data.drug7 <- list("Y","X1","m")
inits.drug7 <- function(){list(mu=0.05,tau=0.03)}
params.drug7 <- c("mu","beta","LPML","y.pred")
out.drug7 <- bugs(data=data.drug7, inits=inits.drug7, parameters.to.save=params.drug7,
                  model.file=model.drug7,n.iter=20000, n.burnin=5000, debug=F,n.chains=1)
samples7.mu <- out.drug7$sims.list$mu
samples7.beta1 <- out.drug7$sims.list$beta[,1]
samples7.beta2 <- out.drug7$sims.list$beta[,2]
samples7.LPML <- out.drug7$sims.list$LPML
samples7.ypred <- out.drug7$sims.list$y.pred
mean(samples7.mu)
mean(samples7.beta1)
mean(samples7.beta2)
mean(samples7.ypred)
pdf(file = "regression2.pdf",width = 9,height=9)
plot(data2[,4],data2[,3],type = "p")
abline(mean(samples7.beta1),mean(samples7.beta2))     ## There is clearly an upward trend depending on where you live.
dev.off()
residual7 <- rep(NA,length(Y))

for (i in 1:length(Y)) {
  y <- mean(samples7.beta1)+mean(samples7.beta2)*data2[i,4]
  residual7[i] <- Y[i]-y
}
max(abs(residual7))       ## Second lowest max residual
min(abs(residual7))       ## Lowest min residual



##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################

## We just tested six different models:
## Zellner's g-prior with three covariates
## Double exponential with three covariates
## Flat prior with three covariates
## Zellner's g-prior with one factor
## Double exponential with one factor
## Flat prior with one factor
## Let's list the LPML and DIc for each to find a superior model.

## Zellner's g-prior with three factors
## LPML
mean(samples2.LPML)           ## -123.1792
## DIC
out.drug2$DIC                 ## -203.9


## Double exponential with three covariates
## LPML
mean(samples3.LPML)           ## -1258577
## DIC
out.drug3$DIC                 ## 1140


## Flat prior with three covariates
## LPML
mean(samples4.LPML)           ## -121.8963
## DIC
out.drug4$DIC                 ## -202.9


## Zellner's g-prior with one factor
## LPML
mean(samples5.LPML)           ## -121.8831
## DIC
out.drug5$DIC                 ## -204.6


## Double exponential with one factor
## LPML
mean(samples6.LPML)           ## 5.136655
## DIC
out.drug6$DIC                 ## 14.37


## Flat prior with one factor
## LPML
mean(samples7.LPML)           ## -121.5848
## DIC
out.drug7$DIC                 ## -204.6


##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################

## Calculating the LPML for both the g-prior w/ three covariates and the flat prior w/ one covariate
## Do not run!! Takes about 15 minutes
cpo1 = rep(NA, length=length(Y))
cpo2 = rep(NA, length=length(Y))
data.drugLPML1 <- list("Y"=Y[-i],"X"=X[-i,],"m"=(length(Y)-1),"taumatrix"=taumatrix,"mean"=mean)
data.drugLPML2 <- list("Y"=Y[-i],"X1"=X1[-i,],"m"=(length(Y)-1))
params.drugLPML <- c("y")


for (i in 1:m) {
  m1 <- bugs(data=data.drugLPML1, inits=inits.drug2, parameters.to.save=params.drugLPML,
             model.file=model.drug2,n.iter=20000, n.burnin=5000, debug=F,n.chains=1)
  m2 <- bugs(data=data.drugLPML2, inits=inits.drug7, parameters.to.save=params.drugLPML,
             model.file=model.drug7,n.iter=20000, n.burnin=5000, debug=F,n.chains=1)
  cpo1[i] <- mean(m1$sims.list$y)
  cpo2[i] <- mean(m2$sims.list$y)
  print(cpo1[i])
}

LPML_H1 <- sum(log(cpo1))
LPML_H2 <- sum(log(cpo2))
LPML_H1
LPML_H2

## Results
# > LPML_H1
# [1] -142.5956
# > LPML_H2
# [1] -141.9897



##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################

model8 <- function(){
  M ~ dcat(p[1:2])
  
  for (i in 1:m){
    Y[i] ~ dnorm(mu[M,i],tau[M])
    mu[1,i] <- inprod(alpha[],X[i,])
    mu[2,i] <- inprod(beta[],X1[i,])
  }
  
  p[1] <- 0.5
  p[2] <- 0.5
  
  for (i in 1:4) {
    alpha[i] ~ dnorm(0,1e-6)
  }
  
  for (i in 1:2) {
    beta[i] ~ dnorm(0,1e-6)
  }
  
  tau[1] ~ dgamma(0.01,0.01)
  tau[2] ~ dgamma(0.01,0.01)
}
model.sales8 <- file.path("C:/Users/rl02898/Documents/Bayesian Statistics/Final Project/model.drug8.txt")
write.model(model8, model.sales8)

data.sales8 <- list("X","X1","m","Y")
inits.sales8 <- function(){list(tau=0.01,alpha=c(0.1,0.1,0.1,0.1),beta=c(0.1,0.1))}
params.sales8 <- c("mu","M")
out.sales8 <- bugs(data=data.sales8, inits=inits.sales8, parameters.to.save=params.sales8,
                    model.file=model.sales8,n.iter=20000, n.burnin=5000, debug=F,n.chains=1)
samples8.M <- out.sales8$sims.list$M
BF12 = sum(samples8.M==1)/sum(samples8.M==2)
BF12



##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################

## According to these statistics, the flat prior with one factor is a more superior model.
## With this being said, let's find a prediction for each region across the US.
## To calculate this, let's take our superior model paramters and plug in the values.
mean(samples7.beta1)+mean(samples7.beta2)*1
mean(samples7.beta1)+mean(samples7.beta2)*2
mean(samples7.beta1)+mean(samples7.beta2)*3
mean(samples7.beta1)+mean(samples7.beta2)*4
mean(samples7.beta1)+mean(samples7.beta2)*5
mean(samples7.beta1)+mean(samples7.beta2)*6
mean(samples7.beta1)+mean(samples7.beta2)*7





