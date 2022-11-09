rm(list=ls())
library(MBA)
library(fields)
library('coda')
library(leaps)
library(spBayes)
library(geoR)

library(leaps)
library(MBA)
library(fields)


events.data <- read.csv(".\\Datasets\\2017_tamil_nadu_subdistrict.csv", header = TRUE)
test.data <-read.csv(".\\Datasets\\Final_Tamilnadu_2017.csv", header = TRUE)

s.all <- as.matrix(events.data[,3:4])
X <- as.matrix(scale(events.data[,6:13]))
y <- log(events.data[,5])
ns <- nrow(s.all)

X.test <- as.matrix(scale(test.data[,5:12]))
y.test <- log(test.data[,4])
s.test <-as.matrix(test.data[,2:3])


#scatter plot of Y vs. X
pairs(cbind(y,X))

#y <- as.matrix(events.data[,5])
b <- regsubsets(X, y)
rs <- summary(b)
rs$outmat

#Linear models for estiamtion
om1 <- lm(y~X[,6])
om2 <- lm(y~X[,1]+X[,6])
om3 <- lm(y~X[,1]+X[,6]+X[,7])
om4 <- lm(y~X[,1]+X[,2]+X[,6]+X[,7])
om5 <- lm(y~X[,1]+X[,2]+X[,3]+X[,6]+X[,7])
om6 <- lm(y~X[,1]+X[,2]+X[,3]+X[,4]+X[,6]+X[,7])
om7 <- lm(y~X[,1]+X[,2]+X[,3]+X[,4]+X[,6]+X[,7]+X[,8])
om8 <- lm(y~X[,1]+X[,2]+X[,3]+X[,4]+X[,5]+X[,6]+X[,7]+X[,8])


om <- list(om1,om2,om3,om4,om5,om6,om7,om8)
n.y <- length(y)

# number of parameters in models in om
m <- 8
npar <- 6:13
#AIC
AIC<- sapply(1:m, function(x) round(extractAIC(om[[x]],k=2)[2],2))
AIC

#BIC
BIC<- sapply(1:m, function(x) round(extractAIC(om[[x]],k=log(n.y))[2],2))
BIC

#all possible subsets selection with PRESS
myPRESS <- function(x,y,indx){
  m1 <- lm(y~x[,indx])
  press <- sum((m1$residuals/(1-hatvalues(m1)))^2)
  return(press)
}

#PRESS Index
PRESS.indx <- matrix("", nrow = m, ncol = m)
colnames(PRESS.indx) <- c("x1", "x2", "x3", "x4","x5", "x6", "x7","x8")
PRESS <- rep(0,m)
for(i in 1:m)
{
  indx <- combn(m,i)
  n.indx <- ncol(indx)
  tmp <- sapply(1:n.indx, function(k) myPRESS(X, y, indx[,k]))
  PRESS[i] <- round(min(tmp),2)
  PRESS.indx[i, indx[,which.min(tmp)]] <- "*"
}


PRESS.indx
PRESS

#Check if residuals have spatial assosciation
surf <- mba.surf(cbind(s.all,om4$residuals), no.X = 102, no.Y = 102, extend = TRUE)$xyz.est
image.plot(surf, xaxs = "r", yaxs = "r", xlab = "Longitude", ylab = "Latitude")
points(s.all, pch = 20)
contour(surf, add = T)


#empirical variograms
om4.sum <- summary(om4)
max.dist <-   0.5*max(iDist(s.all))
#apply variogram methods
bins <- 100
#Matheron estimator
vario.M <- variog(coords = s.all, data = y, trend = ~X[,1:8],
                  estimator.type = "classical", uvec = (seq(0, max.dist, l = bins )))


#Cressie-Hawkins estimator
vario.CH <- variog(coords = s.all, data = y, trend = ~X[,1:8],
                   estimator.type = "modulus", uvec = (seq(0, max.dist, l = bins )))

plot(vario.M, main = "y")
points(vario.CH$u, vario.CH$v, col = "red")
legend("bottomright", legend=c("Matheron", "Cressie-Hawkins"), pch = c(1,1),
       col = c("black", "red"), cex = 0.6)


#WLS with weights N(h)/gamma^2(h)
fit.vario <- variofit(vario.M, ini.cov.pars = c(om4.sum$sigma^2,
                                                -max.dist/log(0.05)), cov.model = "exponential",
                      minimisation.function = "optim", weights = "cressie", control = list(factr=1e-10, maxit = 500))


fit.vario


#Bayesian Hierarchical model
p <- 5
dis <- iDist(s.all)
max.dist <- 0.5 * max(dis)
min.dist <- min(dis[dis!=0])


#pick the variogram estimate as the starting value.
lower.phi <- -log(0.05)/max.dist
upper.phi <- -log(0.01)/min.dist

tuning <- list("phi"=0.5, "sigma.sq"=0.01, "tau.sq"=0.01)
priors <- list("beta.Norm"=list(rep(0,p), diag(1000000,p)),
               "phi.Unif"=c(lower.phi, upper.phi),
               "sigma.sq.IG"=c(2, fit.vario$cov.pars[1]),  "tau.sq.IG"= c(2,0.23))     #c(2, fit.vario$nugget))
w.phi <- c(0.85, 0.9, 0.95)
phi.init <- lower.phi * w.phi + upper.phi * (1-w.phi)
sigmaSq.init <- fit.vario$cov.pars[1] * c(1, 0.5, 2)
tauSq.init <- fit.vario$nugget * c(1, 2, 0.5)

#amcmc
n.batch <- 800
batch.length <- 20
n.samples <- n.batch*batch.length
n.report <- 20

ptm <- proc.time()
starting <- list("phi"= phi.init[1], "sigma.sq"= sigmaSq.init[1], "tau.sq"=tauSq.init[1])
m.chain1 <- spLM(y ~ X[,c(1,2,6,7)], coords=s.all, starting=starting,
                 tuning=tuning, priors=priors, cov.model="exponential",
                 amcmc=list("n.batch"=n.batch, "batch.length"=batch.length, "accept.rate"=0.43),
                 n.report=n.report)


starting <- list("phi"= phi.init[2], "sigma.sq"= sigmaSq.init[2], "tau.sq"=tauSq.init[2])
m.chain2 <- spLM(y ~ X[,c(1,2,6,7)], coords=s.all, starting=starting,
                 tuning=tuning, priors=priors, cov.model="exponential",
                 amcmc=list("n.batch"=n.batch, "batch.length"=batch.length, "accept.rate"=0.43),
                 n.report=n.report)

starting <- list("phi"= phi.init[3], "sigma.sq"= sigmaSq.init[3], "tau.sq"=tauSq.init[3])
m.chain3 <- spLM(y ~ X[,c(1,2,6,7)], coords=s.all, starting=starting,
                 tuning=tuning, priors=priors, cov.model="exponential",
                 amcmc=list("n.batch"=n.batch, "batch.length"=batch.length, "accept.rate"=0.43),
                 n.report=n.report)

#posterior sample of theta
theta.samps <- mcmc.list(m.chain1$p.theta.samples, m.chain2$p.theta.samples,  m.chain3$p.theta.samples)
plot(theta.samps, density = FALSE, auto.layout = FALSE, ask = par(mfrow=c(1,1)))

#Diagnostics Of MCMC converges:
#"potential scale reduction factor" for each parameter obtained by "gelman.diag"
#approximate convergence is diagnosed when the upper confidence limit is close to 1. 
print(gelman.diag(theta.samps))

gelman.plot(theta.samps, auto.layout = FALSE, ask = par(mfrow=c(1,1)))


#we pick the first 10,000 samples as burn-in.  
burn.in <- 10000
ptm <- proc.time()
m.chain1 <- spRecover(m.chain1, start=burn.in, thin=10, verbose=FALSE)
m.chain2 <- spRecover(m.chain2, start=burn.in, thin=10, verbose=FALSE)
m.chain3 <- spRecover(m.chain3, start=burn.in, thin=10, verbose=FALSE)
proc.time() - ptm

m4.nsp<- bayesLMRef(om4, n.samples)

burn.in <- 10000
m4.nsp.Diag <- spDiag(m4.nsp, start = burn.in, verbose = FALSE)
DIC.nsp <- m4.nsp.Diag$DIC
GP.nsp <- m4.nsp.Diag$GP

m1.sp.Diag <- spDiag(m.chain1, verbose = FALSE)
DIC.sp <- m1.sp.Diag$DIC
GP.sp <- m1.sp.Diag$GP

DIC <- cbind(DIC.nsp, DIC.sp)
colnames(DIC) <- c("non-spatial", "spatial")
GP <- cbind(GP.nsp, GP.sp)
colnames(GP) <- c("non-spatial", "spatial")
print(DIC)

print(GP)

#MSPR
#prediction at s.test with non-spatial model
pred.nsp <- spPredict(m4.nsp, start = burn.in, thin = 10,  pred.covars = cbind(1,X.test[,c(1,2,6,7)]))
#posterior predicted value and prediction interval
pred.nsp.summary <- apply(pred.nsp$p.y.predictive.samples, 1, function(x){quantile(x, prob=c(0.025,0.5,0.975))})
MSPR.nsp <- mean((pred.nsp.summary[2,]-y.test)^2)

#prediction at s.test with spatial model
pred.sp <- spPredict(m.chain1, start = burn.in, thin = 10, pred.coords = s.test, pred.covars = cbind(1,X.test[,c(1,2,6,7)]))

#posterior predicted value and prediction interval
pred.sp.summary<- apply(pred.sp$p.y.predictive.samples, 1, function(x){quantile(x, prob=c(0.025,0.5,0.975))})
MSPR.sp <- mean((pred.sp.summary[2,]-y.test)^2)

print(cbind(MSPR.nsp, MSPR.sp))

#plots for comparing the prediction performance of spatial vs non-spatial models
par(mfrow=c(1,2))
#plots comparing observed vs prediction of non-spatial models
a<-min(c(pred.nsp.summary, pred.sp.summary,y.test))
b<-max(c(pred.nsp.summary, pred.sp.summary,y.test))+8
n.test <- length(y.test)
indx <- order(pred.nsp.summary[2,])
plot(c(1:n.test),y.test[indx], typ="l",ylim=c(a,b),xlab=NA,ylab=NA, main = "non-spatial prediction")
#fitted quantiles 0.025,0.5,0.975
polygon(c(c(1:n.test),rev(c(1:n.test))),c(pred.nsp.summary[1, indx],rev(pred.nsp.summary[3, indx])),col="grey90",border=FALSE)
lines(c(1:n.test),y.test[indx],lty=1,col="black")
lines(c(1:n.test),pred.nsp.summary[1, indx],lty=2,col="grey60")
lines(c(1:n.test),pred.nsp.summary[2, indx],lty=4,col="red")
lines(c(1:n.test),pred.nsp.summary[3, indx],lty=2,col="grey60")
legend(1,b,c("95% credible band","predicted y","obs y"),lty=c(2,4,1),lwd=c(2.5,2.5),col=c("grey60","red","black"));

#plots comparing observed vs prediction of spatial models
indx <- order(pred.sp.summary[2,])
plot(c(1:n.test),y.test[indx], typ="l",ylim=c(a,b),xlab=NA,ylab=NA,
     main = "spatial prediction")
#fitted quantiles 0.025,0.5,0.975
polygon(c(c(1:n.test),rev(c(1:n.test))),c(pred.sp.summary[1, indx],
                                          rev(pred.sp.summary[3, indx])),col="grey90",border=FALSE)
lines(c(1:n.test),y.test[indx],lty=1,col="black")
lines(c(1:n.test),pred.sp.summary[1, indx],lty=2,col="grey60")
lines(c(1:n.test),pred.sp.summary[2, indx],lty=4,col="red")
lines(c(1:n.test),pred.sp.summary[3, indx],lty=2,col="grey60")
legend(1,b,c("95% posterior band","predicted",
             "obs"),lty=c(2,4,1),lwd=c(2.5,2.5),col=c("grey60","red","black"))

#Predict Y at new locations where Xs were observed.
pred.data <- read.csv(".\\Datasets\\Final_Tamilnadu_2018.csv", header = TRUE)
s.pred <- as.matrix(pred.data[,2:3])
X.pred <- as.matrix(scale(pred.data[,5:12]))
y.pred <- log(pred.data[,4])

library(sp)
pred <- spPredict(m.chain1, start = burn.in, thin = 10, pred.coords = s.pred, pred.covars = cbind(1,X.pred[,c(1,2,4,7)]))

pred.summary <- apply(pred$p.y.predictive.samples, 1, function(x){quantile(x, prob=c(0.025,0.5,0.975))})

#District Level Map
ind2 = readRDS(".\\gadm36_IND_2_sp.rds") #Download the subdistrict shape file necessary
tn = (ind2[ind2$NAME_1=="Tamil Nadu",])
spplot(tn,"NAME_1", main = "TN Districts", colorkey =F)

quilt.plot(s.pred[,1],s.pred[,2],exp(pred.summary[2,]), xlab = "Longitude", ylab = "Latitude",main="Posterior Predictive Median")
map(tn ,add=TRUE)

hist(exp(pred.summary[2,]))
#prediction error
pred_table <- as.data.frame(cbind(s.pred, exp(pred.summary[2,])))
MSPR.pred <- mean((exp(pred.summary[2,]-y.pred)^2))
MSPR.pred

rng95 <- pred.summary[3,]-pred.summary[1,]
#District Level Map
ind2 = readRDS(".\\gadm36_IND_2_sp.rds") #Download the subdistrict shape file
tn = (ind2[ind2$NAME_1=="Tamil Nadu",])
spplot(tn,"NAME_1", main = "TN Districts", colorkey =F)

quilt.plot(s.pred[,1],s.pred[,2],rng95, xlab = "Longitude", ylab = "Latitude",main="Posterior Predictive 95% CI range")
map(tn ,add=TRUE)

