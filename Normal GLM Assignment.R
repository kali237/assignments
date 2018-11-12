#Normal GLM assignment
library(carData)
library(car)

#setwd(dir="C:\\Users\\kwickens\\Documents\\Regression and Non Parametric\\Normal GLM")
setwd(dir = "C:\\Users\\Kali\\Documents\\7880 Assingments\\Assignment 2 - Normal GLM")

load("Normal_GLM.Rdata")

############### functions ############### 
lik=function(beta1,beta2,y,X1,X2){
  #y consists of n observations
  #X1 is an n-by-p1 matrix, with p1 covariates for each subject
  #X2 is an n-by-p2 matrix, with p2 covariates for each subject
  #beta1 is a vector of length p1 OR a m-by-p1 matrix
  #beta2 is a vector of length p2 OR a m-by-p2 matrix
  if (is.matrix(beta1)){
    beta1=t(beta1)
    beta2=t(beta2)
  }
  L=t(y)%*%X1%*%beta1 - t(y^2)%*%X2%*%beta2/2
  if (is.matrix(beta1)){
    L=L - 0.5*colSums((X1%*%beta1)^2/(X2%*%beta2))+0.5*colSums(log(X2%*%beta2))
    L=t(L)
    colnames(L)="Log-likelihood"
  } else {
    L=L - 0.5*sum((X1%*%beta1)^2/(X2%*%beta2))+0.5*sum(log(X2%*%beta2))
    L=L[1,1]
  }
  return(L-0.5*log(2*pi)*dim(X1)[1]) #returns a vector with m values of the likelihood
}

Dlik=function(beta1,beta2,y,X1,X2){
  #y consists of n observations
  #X1 is an n-by-p1 matrix, with p1 covariates for each subject
  #X2 is an n-by-p2 matrix, with p2 covariates for each subject
  #beta1 is a vector of length p1 OR a m-by-p1 matrix
  #beta2 is a vector of length p2 OR a m-by-p2 matrix
  
  if (is.matrix(beta1)){
    m=dim(beta1)[1]
    L=t(t(rep(1,m)))%*%c(t(y)%*%X1, - t(y^2)%*%X2/2)
    L=L - cbind(t((X1%*%t(beta1))/(X2%*%t(beta2)))%*%X1, -0.5*t((X1%*%t(beta1))^2/(X2%*%t(beta2))^2)%*%X2 -0.5*t(1/(X2%*%t(beta2)))%*%X2)   
  } else {
    L=c(t(y)%*%X1, - t(y^2)%*%X2/2)
    L=L - c(t((X1%*%beta1)/(X2%*%beta2))%*%X1, -0.5*t((X1%*%beta1)^2/(X2%*%beta2)^2)%*%X2 -0.5*t(1/(X2%*%beta2))%*%X2)   
  }
  if (is.matrix(beta1)) {
    colnames(L)=c(colnames(beta1),colnames(beta2))
    rownames(L)=rownames(beta1)} else {
      names(L)=c(names(beta1),names(beta2))
    }
  return(L) #returns m-by-(p1+p2) matrix with the m gradient vectors
}

D2lik=function(beta1,beta2,y,X1,X2){
  #y consists of n observations
  #X1 is an n-by-p1 matrix, with p1 covariates for each subject
  #X2 is an n-by-p2 matrix, with p2 covariates for each subject
  #beta1 is a vector of length p1 (only one parameter!)
  #beta2 is a vector of length p2 (only one parameter!)
  
  p1=length(beta1)
  p2=length(beta2)
  eta1=X1%*%beta1
  eta2=X2%*%beta2
  L=matrix(0,p1+p2,p1+p2)
  for (m in 1:(p1+p2)){
    for(n in 1:(p1+p2)){
      if (m<=p1 & n<=p1){
        L[m,n]=-sum(X1[,m]*X1[,n]/eta2)
      }
      if (m<=p1 & n>p1){
        L[m,n]=sum(eta1*X1[,m]*X2[,n-p1]/eta2^2)
      }
      if (m>p1 & n<=p1){
        L[m,n]=sum(eta1*X2[,m-p1]*X1[,n]/eta2^2)
      }
      if (m>p1 & n>p1){
        L[m,n]=-sum(X2[,m-p1]*X2[,n-p1]*(eta1^2/eta2^3 + 0.5/eta2^2))
      }
    }
  }
  #colnames(L)=c(names(beta1),names(beta2))
  #rownames(L)=c(names(beta1),names(beta2))
  
  return(L) #returns (p1+p2)-by-(p1+p2) second derivative matrix
}

#########################################



#a)

x_2 = x**2

fit = lm(y ~ x + x_2) #yes, very significant
summary(fit)

#just for fun
fit0 = lm(y~x)
summary(fit0)

SSres = sum(fit$res^2)
SSresH0 = sum(fit0$res^2)

f = (SSresH0 - SSres)/SSres * (47/1)
pvalue = 1-pf(f,1,47)

#plot the residuals
plot(fit)
plot(x=x, y = fit$res, ylab = "residuals", xlab = NULL)
plot(x=x_2, y = fit$res, ylab = "residuals", xlab = NULL)
qqplot(fit$residuals)


#b) fit model is a submodel assuming equivariance
#for the \eta_{i1} this is basically scaling the fitted values from the first model and scaling by the variance. Because we are assuming equivariance, we just divide this all by the variance (or the standard error squared)
#for \eta_{i2} we do not want the value of the variance to change based on the data, therefore we want \beta_{21} and \beta_{22} to be 0, and just scale the intercept by the variance

sqrt(SSres/47) #residual standard error, \sigma^hat

sig2_hat = SSres/47 #estimate for the variance

beta1 = fit$coefficients/sig2_hat

beta2 = c(1/sig2_hat,0,0)

X_1 = model.matrix(fit)
X_2 = model.matrix(fit)

lik(beta1,beta2,y,X_1,X_2)
Dlik(beta1,beta2,y,X_1,X_2)
#the first three are 0 because we already found in model a. the last 3 are not, thus we are not at the minimum (MLE), and thus our assumption of equivariance is wrong

#c) stupid optimization shit
#use the beta1 and beta2 and the lik functions

l = lik(beta1,beta2,y,X_1,X_2)
Dl = Dlik(beta1,beta2,y,X_1,X_2)
D2l = D2lik(beta1,beta2,y,X_1,X_2)

D2l_inv = -solve(D2l)

h = D2l_inv%*%Dl

b = c(beta1,beta2)

#stuff from Nick

D2l=function(b)  {D2lik(b[1:3],b[4:6],y,model.matrix(fit),model.matrix(fit))}
Dl=function(b)  {Dlik(b[1:3],b[4:6],y,model.matrix(fit),model.matrix(fit))}
l=function(b)  {lik(b[1:3],b[4:6],y,model.matrix(fit),model.matrix(fit))}
h_m=function(b){-solve(D2l(b))%*%Dl(b)}
opt = optim(b,l,Dl,control= list(fnscale = -1))

summary(opt)
n_beta = opt$par
n_beta1 = n_beta[1:3]
n_beta2 = n_beta[4:6]

lik(n_beta1,n_beta2,y,X_1,X_2) #cool, give -40.25

eta1_1 = X_1%*%n_beta1
eta2_1 = X_1%*%n_beta2

#Newton-Raphson! This now works. It produces a better minimum than 'optim'
h=h_m(b)
del=1
while(del>10^-35){ # 10^-35 is around 100 times the size of single-precision 'machine zero.'
  #In reality this will be satisfied when l(b) and l(b+h) agree to ~16 digits,
  #which will happen much much sooner.
  while(is.nan(l(b+h))||l(b+h)<l(b)){
    h = h/2
  }
  del = l(b+h)-l(b)
  b = b+h
  h = h_m(b)
  b
}
#Load 'b' with the \beta coefficients
#simulate; y = N(\eta1/\eta2, 1/\eta2)

lik(b[1:3],b[4:6],y,X_1,X_2) #gives -38.13

eta1 = model.matrix(fit)%*%b[1:3] #\eta_{1,i} = X_i*\beta_1
eta2 = model.matrix(fit)%*%b[4:6] #\eta_{2,i} = X_i*\beta_2

#plots plots plots!!

#Plot of 'y' values, fitted values from the GLM, and fitted values from the linear model
#Black is the 'y' values, red is the linear model and green is the GLM
plot(y)
points(fit$fitted.values,col="red")
points(eta1/eta2, col="blue")

#Residual plots
#Residuals from the GLM
plot(y - eta1/eta2)
points(sqrt(1/eta2),pch="-", col = "red")
points(-sqrt(1/eta2),pch="-", col = "red")

#Residuals from the linear model
plot(fit$residuals)
points(sqrt(sig2_hat)*replicate(50,1),pch="-", col = "red")
points(sqrt(sig2_hat)*replicate(50,-1),pch="-", col = "red")
#It looks good, since there are lots of points within the +/- 1 stdev bars...
#But actually there are TOO MANY points within the bars! We should only have ~70% between the bars!
#~30% of the points SHOULD be outside the bars, but only 7 of 50 (14%) points are outside the bars.

#d) test the quadratic terms for both \eta_1 and \eta_2. Then test which model is better GLM or standard linear model

eta1_fit = lm(eta1 ~ x + x_2)
summary(eta1_fit)

b_h0 = c(b[1:2],0)

eta1_fit0 = lm(eta1 ~ x)
summary(eta1_fit0)

eta2_fit = lm(eta2 ~ x + x_2)
summary(eta2_fit)


#N(m,s) = m + \sqrt(s)N(0,1)
n  = rnorm(50)
sim = eta1/eta2 + sqrt(1/eta2)*n

