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

events.data <- read.csv("./Dataset/2017_India.csv", header = TRUE)

s.all <- as.matrix(events.data[,2:3])
s.districts <- as.matrix(events.data[,1])
y <- log(events.data[,4])
X <- as.matrix(scale(events.data[,5:7]))
ns <- nrow(s.all)

#random splitting the data to training and testing data.(0.7 training, 0.3 testing)
set.seed(1)
s.train.indx <- sample(1:ns, 0.7*ns)
s.train <- as.matrix(s.all[s.train.indx,])
s.train.districts <- as.matrix(s.districts[s.train.indx,])
s.test <- as.matrix(s.all[-s.train.indx,])
s.test.districts <- as.matrix(s.districts[-s.train.indx,])
y.train <- as.matrix(y[s.train.indx])
y.test <- as.matrix(y[-s.train.indx])
X.train <- as.matrix(X[s.train.indx, ])
X.test <- as.matrix(X[-s.train.indx, ])

max.dist <-  0.5 * max(iDist(s.train))

# Fit the training data using REML (get error if number of iteration(maxit) is increased greater than 20: Error in optim(par = c(18.5991373589123, 0), fn = function (pars, fp,:non-finite value supplied by optim )
fit.reml <- likfit(coords = s.train, data = y.train, trend = "cte", ini.cov.pars =
                     c(var(y.train), -max.dist/log(0.05)), cov.model = "exponential",
                   lik.method = "REML", control = list(factr=1e-10, maxit = 500))

print(fit.reml)


# prediction using ordinary krigging for test data using the estimated parameters
y.obj <- as.geodata(cbind(s.train, y.train), coords.col = 1:2, data.col = 3)
pred <- krige.conv(y.obj, coords=s.train,locations=s.test, krige=
                     krige.control(cov.model="matern",cov.pars=fit.reml$cov.pars,
                                   nugget=fit.reml$nugget, kappa = fit.reml$kappa, micro.scale = 0))

plot(y.test, pred$predict)
abline(lsfit(y.test,pred$predict),lty=2,col=1)

MSPR.onlykrigging <- mean((pred$predict-y.test)^2)
print(MSPR.onlykrigging)


#interpolated 2D plot of the prediction standard errors
nsx <- 155
nsy <- 155
#B-spline interpolation
surf <- mba.surf(cbind(s.test,sqrt(pred$krige.var)), no.X = nsx, no.Y = nsy, extend = FALSE)$xyz.est
image.plot(surf, xaxs = "r", yaxs = "r", xlab = "x(ft)", ylab = "y(ft)")
points(s.test, pch = 20)
contour(surf, add = T)

final.preds <- data.frame(District = s.test.districts, Lat =s.test[,1], Long = s.test[,2], Predictions = pred$predict, Truth = y.test[,1])

write.csv(x=final.preds, file="D:/UNL/Second Sem/STAT 832 Intro to Spatial Statistics/CLass Project/SummerProject/KrigPred2018.csv")




pairs(cbind(y,X))


#regression model

om1 <- lm(y.train~X.train[,1]+X.train[,2]+X.train[,3])
#om2 <- lm(y.train~X.train[,1]+X.train[,3])
om2 <- lm(y.train~X.train[,1])
om3 <- lm(y.train~1)


#Check if residuals have spatial assosciation
surf <- mba.surf(cbind(s.train,om2$residuals), no.X = 379, no.Y = 379, extend = TRUE)$xyz.est
image.plot(surf, xaxs = "r", yaxs = "r", xlab = "Long", ylab = "Lat", ylim = rev(range(s.train[,2])))
points(s.train, pch = 20)
contour(surf, add = T)


#empirical variograms
om2.sum <- summary(om2)
max.dist <-   0.5*max(iDist(s.all))
#apply variogram methods
bins <- 100
#Matheron estimator
vario.M <- variog(coords = s.train, data = y.train, trend = ~X.train[,c(1)],
                  estimator.type = "classical", uvec = (seq(0, max.dist, l = bins )))


#Cressie-Hawkins estimator
vario.CH <- variog(coords = s.train, data = y.train, trend = ~X.train[,c(1)],
                   estimator.type = "modulus", uvec = (seq(0, max.dist, l = bins )))

plot(vario.M, main = "y")
points(vario.CH$u, vario.CH$v, col = "red")
legend("bottomright", legend=c("Matheron", "Cressie-Hawkins"), pch = c(1,1),
       col = c("black", "red"), cex = 0.6)


#WLS with weights N(h)/gamma^2(h)
fit.vario <- variofit(vario.CH, ini.cov.pars = c(om2.sum$sigma^2,
                                                 -max.dist/log(0.05)), cov.model = "exponential",
                      minimisation.function = "optim", weights = "cressie", control = list(factr=1e-10, maxit = 500))


fit.vario


#Bayesian Hierarchical model
p <- 2
dis <- iDist(s.train)
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
m.chain1 <- spLM(y.train ~ X.train[,c(1)], coords=s.train, starting=starting,
                 tuning=tuning, priors=priors, cov.model="exponential",
                 amcmc=list("n.batch"=n.batch, "batch.length"=batch.length, "accept.rate"=0.43),
                 n.report=n.report)


starting <- list("phi"= phi.init[2], "sigma.sq"= sigmaSq.init[2], "tau.sq"=tauSq.init[2])
m.chain2 <- spLM(y.train ~ X.train[,c(1)], coords=s.train, starting=starting,
                 tuning=tuning, priors=priors, cov.model="exponential",
                 amcmc=list("n.batch"=n.batch, "batch.length"=batch.length, "accept.rate"=0.43),
                 n.report=n.report)

starting <- list("phi"= phi.init[3], "sigma.sq"= sigmaSq.init[3], "tau.sq"=tauSq.init[3])
m.chain3 <- spLM(y.train ~ X.train[,c(1)], coords=s.train, starting=starting,
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


#we pick the first 12,000 samples as burn-in.  
burn.in <- 12000
ptm <- proc.time()
m.chain1 <- spRecover(m.chain1, start=burn.in, thin=10, verbose=FALSE)
m.chain2 <- spRecover(m.chain2, start=burn.in, thin=10, verbose=FALSE)
m.chain3 <- spRecover(m.chain3, start=burn.in, thin=10, verbose=FALSE)
proc.time() - ptm

m1.nsp<- bayesLMRef(om2, n.samples)

burn.in <- 12000
m1.nsp.Diag <- spDiag(m1.nsp, start = burn.in, verbose = FALSE)
DIC.nsp <- m1.nsp.Diag$DIC
GP.nsp <- m1.nsp.Diag$GP

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
pred.nsp <- spPredict(m1.nsp, start = burn.in, thin = 10,  pred.covars = cbind(1,X.test[,c(3)]))
#posterior predicted value and prediction interval
pred.nsp.summary <- apply(pred.nsp$p.y.predictive.samples, 1, function(x){quantile(x, prob=c(0.025,0.5,0.975))})
MSPR.nsp <- mean((pred.nsp.summary[2,]-y.test)^2)

#prediction at s.test with spatial model
pred.sp <- spPredict(m.chain1, start = burn.in, thin = 10, pred.coords = s.test, pred.covars = cbind(1,X.test[,c(3)]))

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
library(sp)
pred <- spPredict(m.chain1, start = burn.in, thin = 10, pred.coords = s.test, pred.covars = cbind(1,X.test[,c(3)]))

pred.summary <- apply(pred$p.y.predictive.samples, 1, function(x){quantile(x, prob=c(0.025,0.5,0.975))})

hist(exp(pred.summary[2,]))

quilt.plot(s.test[,1],s.test[,2],exp(pred.summary[2,]),main="Predictions")

a <- as.data.frame(cbind(s.test,y.test, exp(pred.summary[2,])))


rng95 <- pred.summary[3,]-pred.summary[1,]
b <- as.data.frame(cbind(s.pred, rng95))
