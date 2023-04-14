library(tidyverse)
library(tictoc)


tic()

rm(list=ls()) 		# Clear workspace

set.seed(1) 		# Set seed for random number generator

x <- c(-0.1,0.2,0.7)
B <- c(10,100,1000)

## Question 1

distribution <- c()

# option 1
for (b in B) {
  mean = matrix(0,B,1)
  for (i in 1:b) {
    distribution[i] <- mean(sample(x, 
                            size = length(x), 
                            replace = T))
                }
  emp_cdf <- ecdf(distribution)
  plot(emp_cdf, col=3)
}

# option 2

mean_matrix = matrix(0,B,1)

for (it in 1:B) {
  ## for each sample
  drawn_sample <- sample(x,length(x),replace=TRUE)
  ## write the expectation into the matrix
  mean_matrix[it] = mean(drawn_sample)
}
mean_matrix  %>% as_tibble()  %>% ggplot(aes(V1)) +
  stat_ecdf(geom = "step")




#Question B




library(numDeriv) 	# loads the functions grad and hessian which numerically evaluate the gradient and hessian
library(tictoc)

tic()

rm(list=ls()) 		# Clear workspace

set.seed(1) 		# Set seed for random number generator

n = 300 			# Set the sample size

theta0 = c(1,1) 	# True parameter value

k = length(theta0) 

theta_start = rep(0,k)

# Set working directory
# setwd("...")

# Load the log-likelihood, its derivative, and the hessian
Probit_LL_h <- function (y,x,par) {
  
  Phi = pnorm(x %*% par)
  phi = dnorm(x %*% par)
  
  n = length(y) 
  k = length(par)
  
  # Computing the Hessian
  h = t(x)%*% (matrix(rep((1-y)*(phi/(1-Phi)*(x%*%par)-(phi^2)/(1-Phi)^2),k),nrow=n)*x) + 
    t(x)%*% (matrix(rep((y)*(phi/Phi*(-x%*%par)-(phi^2)/Phi^2),k),nrow=n)*x)
  h = -h
  
  return(h)
}

J_1 <- function(y,x,par) {
  
  Phi = pnorm(x %*% par)
  phi = dnorm(x %*% par)
  
  n = length(y) 
  k = length(par)
  
  g = matrix(rep(y*phi/Phi - (1-y)*phi/(1-Phi),k),nrow=n)*x
  
  f = t(g)%*%g
  
  return(f)
}

J_3 <- function(y,x,par) {
  Phi = pnorm(x %*% par)
  phi = dnorm(x %*% par)
  
  n = length(y) 
  k = length(par)
  
  f = t(x)%*%(matrix(rep(phi/Phi*phi/(1-Phi),k),nrow=n)*x)
  
  return(f)
}

num = 100			# Number of Monte Carlo iterations

theta_hat_vec = matrix(0,num,k)

J_1_inv = matrix(0,k,k)
J_2_inv = matrix(0,k,k)
J_3_inv = matrix(0,k,k)

B = 399

J_inv_boot = matrix(0,k,k)

for (it in 1:num) {
  # Data generating process
  x = cbind(matrix(1,n,1),matrix(rnorm(n,0,1),ncol=1)) 	# regressors
  
  u = rnorm(n,0,1)							 	# error term
  
  y_star = x %*% theta0 + u						# latent "utility"
  
  y = ceiling(y_star/(max(abs(y_star))+0.1))				# observed outcome
  
  dat = data.frame(x,y)
  probit <- glm(y~x-1, data=dat, family = binomial(link = "probit"))
  theta_hat = probit$coefficients
  
  theta_hat_vec[it,1:k] = theta_hat
  
  J_1_inv = J_1_inv + solve(J_1(y,x,theta_hat))
  
  J_2_inv = J_2_inv + solve(Probit_LL_h(y,x,theta_hat))
  
  J_3_inv = J_3_inv + solve(J_3(y,x,theta_hat))
  
  theta_hat_boot = matrix(0,B,k)
  #empty matrice for every column
  
  for (b in 1:B) {
    index_b <- sample(length(y),length(y),replace=TRUE)
    y_b <- y[index_b]
    x_b <- x[index_b,]
    dat_b = data.frame(x_b,y_b)
    probit <- glm(y_b~x_b-1, data=dat_b, family = binomial(link = "probit"))
    theta_hat_boot[b,1:k] = probit$coefficients
  }
  
  J_inv_boot = J_inv_boot + var(theta_hat_boot)
}

# Average of theta_hat over Monte Carlo iterations
colMeans(theta_hat_vec)

# Variance of theta_hat over Monte Carlo iterations
var(theta_hat_vec)

# Average of variance estimate based on J_1
J_1_inv/num

# Average of variance estimate based on J_2
J_2_inv/num

# Average of variance estimate based on J_3
J_3_inv/num

# Average of variance estimate based on J_3
J_inv_boot/num 

toc()






## Question C:

#ML we already defined it before in Probit_LL_h
#bootstrap testing

theta0 = c(1,1) 	# True parameter value
k = length(theta0) 
theta_start = rep(0,k)
num = 1000

# NLS
Probit_NLS <- function(y,x,par) {
  n = length(y)
  Phi = pnorm(x %*% par)
  f = mean((y-Phi)^2)
  return(f)
}

# GMM

Probit_GMM <- function(y,x,par) {
  n = length(y)
  Phi = pnorm(x %*% par)
  m = t((y-Phi)) %*% x
  f = (m %*% t(m))/n^2
  return(f)
}

# Average marginal effect for each estimator
ame_vec_ML = matrix(0,num,1)
ame_vec_NLS = matrix(0,num,1)
ame_vec_GMM = matrix(0,num,1)

for (it in 1:num) {
  x = cbind(matrix(1,n,1),matrix(rnorm((k-1)*n,0,1),ncol=k-1)) 	# regressors
  
  u = (rchisq(n,1)-1)/sqrt(2)		# error term
  
  y_star = x %*% theta0 + u						# latent "utility"
  
  y = ceiling(y_star/(max(abs(y_star))+0.1))				# observed outcome
  
  dat = data.frame(x,y)
  
  # ML
  probit_ML <- glm(y~x-1, data=dat, family = binomial(link = "probit"))
  theta_hat_ML = probit_ML$coefficients
  
  margeff_ML = dnorm(x %*% theta_hat_ML)*theta_hat_ML[2]
  ame_ML = mean(margeff_ML)
  ame_vec_ML[it,1:1] = ame_ML
  
  # NLS
  probit_NLS <- optim(par = theta0, Probit_NLS, y = y, x = x, method = c("BFGS"), control = list(reltol=1e-9), hessian=TRUE)
  theta_hat_NLS = probit_NLS$par
  
  margeff_NLS = dnorm(x %*% theta_hat_NLS)*theta_hat_NLS[2]
  ame_NLS = mean(margeff_NLS)
  ame_vec_NLS[it,1:1] = ame_NLS  
  
  # GMM
  probit_GMM <- optim(par = theta0, Probit_GMM, y = y, x = x, method = c("BFGS"), control = list(reltol=1e-9), hessian=TRUE)
  theta_hat_GMM = probit_GMM$par
  
  margeff_GMM = dnorm(x %*% theta_hat_GMM)*theta_hat_GMM[2]
  ame_GMM = mean(margeff_GMM)
  ame_vec_GMM[it,1:1] = ame_GMM
}

# Average of average marginal effect for each estimator
mean_ame_ML = mean(ame_vec_ML)
mean_ame_NLS = mean(ame_vec_NLS)
mean_ame_GMM = mean(ame_vec_GMM)

#see them
mean_ame_ML
mean_ame_NLS
mean_ame_GMM

#second stage, you can increase large number

large_number = 10000

x = cbind(matrix(1,large_number,1),matrix(rnorm((k-1)*large_number,0,1),ncol=k-1)) 	# regressors
ame_true = dnorm(x %*% theta0) # we omit "* theta_1 [= 1]"

mean_ame_true = mean(ame_true)

m = cbind(mean_ame_ML,mean_ame_NLS,mean_ame_GMM,mean_ame_true)

library(knitr)
library(kableExtra)


rownames(m) = c("Value")
kable(list(m), col.names = c("ML","NLS","GMM","Truth approx."), 
      escape = F, 
      caption = "") %>%
  kable_styling(latex_options = "hold_position")

#to compare it with the well specified model just change u






################### question d #########################?

set.seed(1) 		# Set seed for random number generator

n = 300			# Set the sample size

num = 1000
theta0 = c(1,1) 	# True parameter value

theta_null = c(1,1) 

k = length(theta0) 

theta_start = theta0

# define functions

Probit_LL <- function(y,x,par) {
  
  Phi = pnorm(x %*% par)
  phi = dnorm(x %*% par)
  
  n = length(y) 
  k = length(par)
  
  # Computing the log-likelihood
  f = sum(y*log(Phi)) + sum((1-y)*log(1-Phi))
  f = -f/n
  
  return(f)
}


constrained_Probit_LL <- function(y,x,par_1) {
  
  par_2 = theta_null[2]
  par = c(par_1,par_2)
  Phi = pnorm(x %*% par)
  phi = dnorm(x %*% par)
  
  n = length(y) 
  k = length(par)
  
  # Computing the log-likelihood
  f = sum(y*log(Phi)) + sum((1-y)*log(1-Phi))
  f = -f/n
  
  return(f)
}


Probit_LL_g <- function (y,x,par) {
  
  Phi = pnorm(x %*% par)
  phi = dnorm(x %*% par)
  
  n = length(y) 
  k = length(par)
  
  # Computing the gradient
  g = t(matrix(rep(phi/Phi,k),nrow=n)*x) %*% y - 
    t(matrix(rep(phi/(1-Phi),k),nrow=n)*x) %*% (1-y)
  g = -g
  
  return(g)
}

J_1 <- function(y,x,par) {
  
  Phi = pnorm(x %*% par)
  phi = dnorm(x %*% par)
  
  n = length(y) 
  k = length(par)
  
  g = matrix(rep(y*phi/Phi - (1-y)*phi/(1-Phi),k),nrow=n)*x
  
  f = t(g)%*%g
  
  return(f)
}

# create empty matrices 
J1_vec = matrix(0,num,1)
J2_vec = matrix(0,num,1)
J_2_sum = matrix(rep(0,4),2)

theta_hat_vec = matrix(0,num,k)

Wald_vec = matrix(0,num,1)

reject_Wald = matrix(0,num,1)
reject_Score = matrix(0,num,1)
reject_LR = matrix(0,num,1)

cv = qchisq(.95, df=1) 

for (it in 1:num) {
  # Data generating process
  x = cbind(matrix(1,n,1),matrix(rnorm((k-1)*n,0,1),ncol=k-1)) 	# regressors
  
  u = rnorm(n,0,1)							 	# error term
  
  y_star = x %*% theta0 + u						# latent "utility"
  
  y = ceiling(y_star/(max(abs(y_star))+0.1))				# observed outcome
  
  dat = data.frame(x,y)
  probit <- glm(y~x-1, data=dat, family = binomial(link = "probit"))
  theta_hat = probit$coefficients
  
  theta_hat_vec[it,1:k] = theta_hat
  
  result <- optim(par = theta_hat, Probit_LL, y = y, x = x, method = c("BFGS"), control = list(reltol=1e-9), hessian=TRUE)
  
  
  ###constrained model 
  result_constrained <- optim(par = theta_null[1], constrained_Probit_LL, y = y, x = x, method = c("BFGS"), control = list(reltol=1e-9), hessian=TRUE)
  
  theta_constrained = result_constrained$par 
  theta_transplant = c(theta_constrained,theta_null[2])
  
  #J_2 = (1/n)*(result$hessian)
  #V = solve(J_2)
  #V = solve((1/n)*V)
  J_2 = 1/n*solve(result$hessian)
  
  # Wald 
  temp=theta_hat[2]-theta_null[2]
  Wald = temp^2/J_2[2,2]
  Wald_vec[it]= Wald
  
  
  #temp = theta_hat[2] - theta_null[2]
  
  #Wald = t(temp)%*%V[2,2]%*%temp
  
  if (is.nan(Wald) == 0) {
    if (Wald > cv) {
      reject_Wald[it] = 1}
  }
  
  # Score
  #Score = Probit_LL_g(y,x,theta_null)
  
  #Score = t(Score)[2]%*%solve(J_1(y,x,theta_hat))[2,2]%*%Score[2]
  J1_vec[it] = J_1(y,x,theta_transplant)
  Score = Probit_LL_g(y,x,theta_transplant)
  
  Score = t(Score)%*%solve(J_1(y,x,theta_transplant))%*%Score
  
  if (is.nan(Score) == 0) {
    if (Score > cv) {
      reject_Score[it] = 1}
  }
  
  
  # Likelihood ratio
  LR = 2*(constrained_Probit_LL(y,x,theta_constrained)*n - Probit_LL(y,x,theta_hat)*n)
  
  if (is.nan(LR) == 0) {
    if (LR > cv) {
      reject_LR[it] = 1 }   
  }
}


colMeans(reject_Wald)
colMeans(reject_Score)
colMeans(reject_LR)




##########Question E########################################
#we do the same but changing theta, in this case theta=0.9

set.seed(1) 		# Set seed for random number generator

n = 300			# Set the sample size

num = 1000
theta0 = c(1,1) 	# True parameter value

theta_null = c(0.9,0.9) 

k = length(theta0) 

theta_start = theta0

# define functions

Probit_LL <- function(y,x,par) {
  
  Phi = pnorm(x %*% par)
  phi = dnorm(x %*% par)
  
  n = length(y) 
  k = length(par)
  
  # Computing the log-likelihood
  f = sum(y*log(Phi)) + sum((1-y)*log(1-Phi))
  f = -f/n
  
  return(f)
}


constrained_Probit_LL <- function(y,x,par_1) {
  
  par_2 = theta_null[2]
  par = c(par_1,par_2)
  Phi = pnorm(x %*% par)
  phi = dnorm(x %*% par)
  
  n = length(y) 
  k = length(par)
  
  # Computing the log-likelihood
  f = sum(y*log(Phi)) + sum((1-y)*log(1-Phi))
  f = -f/n
  
  return(f)
}


Probit_LL_g <- function (y,x,par) {
  
  Phi = pnorm(x %*% par)
  phi = dnorm(x %*% par)
  
  n = length(y) 
  k = length(par)
  
  # Computing the gradient
  g = t(matrix(rep(phi/Phi,k),nrow=n)*x) %*% y - 
    t(matrix(rep(phi/(1-Phi),k),nrow=n)*x) %*% (1-y)
  g = -g
  
  return(g)
}

J_1 <- function(y,x,par) {
  
  Phi = pnorm(x %*% par)
  phi = dnorm(x %*% par)
  
  n = length(y) 
  k = length(par)
  
  g = matrix(rep(y*phi/Phi - (1-y)*phi/(1-Phi),k),nrow=n)*x
  
  f = t(g)%*%g
  
  return(f)
}

# create empty matrices 
J1_vec = matrix(0,num,1)
J2_vec = matrix(0,num,1)
J_2_sum = matrix(rep(0,4),2)

theta_hat_vec = matrix(0,num,k)

Wald_vec = matrix(0,num,1)

reject_Wald = matrix(0,num,1)
reject_Score = matrix(0,num,1)
reject_LR = matrix(0,num,1)

cv = qchisq(.95, df=1) 

for (it in 1:num) {
  # Data generating process
  x = cbind(matrix(1,n,1),matrix(rnorm((k-1)*n,0,1),ncol=k-1)) 	# regressors
  
  u = rnorm(n,0,1)							 	# error term
  
  y_star = x %*% theta0 + u						# latent "utility"
  
  y = ceiling(y_star/(max(abs(y_star))+0.1))				# observed outcome
  
  dat = data.frame(x,y)
  probit <- glm(y~x-1, data=dat, family = binomial(link = "probit"))
  theta_hat = probit$coefficients
  
  theta_hat_vec[it,1:k] = theta_hat
  
  result <- optim(par = theta_hat, Probit_LL, y = y, x = x, method = c("BFGS"), control = list(reltol=1e-9), hessian=TRUE)
  
  
  ###constrained model 
  result_constrained <- optim(par = theta_null[1], constrained_Probit_LL, y = y, x = x, method = c("BFGS"), control = list(reltol=1e-9), hessian=TRUE)
  
  theta_constrained = result_constrained$par 
  theta_transplant = c(theta_constrained,theta_null[2])
  
  #J_2 = (1/n)*(result$hessian)
  #V = solve(J_2)
  #V = solve((1/n)*V)
  J_2 = 1/n*solve(result$hessian)
  
  # Wald 
  temp=theta_hat[2]-theta_null[2]
  Wald = temp^2/J_2[2,2]
  Wald_vec[it]= Wald
  
  
  #temp = theta_hat[2] - theta_null[2]
  
  #Wald = t(temp)%*%V[2,2]%*%temp
  
  if (is.nan(Wald) == 0) {
    if (Wald > cv) {
      reject_Wald[it] = 1}
  }
  
  # Score
  #Score = Probit_LL_g(y,x,theta_null)
  
  #Score = t(Score)[2]%*%solve(J_1(y,x,theta_hat))[2,2]%*%Score[2]
  J1_vec[it] = J_1(y,x,theta_transplant)
  Score = Probit_LL_g(y,x,theta_transplant)
  
  Score = t(Score)%*%solve(J_1(y,x,theta_transplant))%*%Score
  
  if (is.nan(Score) == 0) {
    if (Score > cv) {
      reject_Score[it] = 1}
  }
  
  
  # Likelihood ratio
  LR = 2*(constrained_Probit_LL(y,x,theta_constrained)*n - Probit_LL(y,x,theta_hat)*n)
  
  if (is.nan(LR) == 0) {
    if (LR > cv) {
      reject_LR[it] = 1 }   
  }
}


colMeans(reject_Wald)
colMeans(reject_Score)
colMeans(reject_LR)



######################################
###################for theta=1.1



set.seed(1) 		# Set seed for random number generator

n = 300			# Set the sample size

num = 1000
theta0 = c(1,1) 	# True parameter value

theta_null = c(1.1,1.1) 

k = length(theta0) 

theta_start = theta0

# define functions

Probit_LL <- function(y,x,par) {
  
  Phi = pnorm(x %*% par)
  phi = dnorm(x %*% par)
  
  n = length(y) 
  k = length(par)
  
  # Computing the log-likelihood
  f = sum(y*log(Phi)) + sum((1-y)*log(1-Phi))
  f = -f/n
  
  return(f)
}


constrained_Probit_LL <- function(y,x,par_1) {
  
  par_2 = theta_null[2]
  par = c(par_1,par_2)
  Phi = pnorm(x %*% par)
  phi = dnorm(x %*% par)
  
  n = length(y) 
  k = length(par)
  
  # Computing the log-likelihood
  f = sum(y*log(Phi)) + sum((1-y)*log(1-Phi))
  f = -f/n
  
  return(f)
}


Probit_LL_g <- function (y,x,par) {
  
  Phi = pnorm(x %*% par)
  phi = dnorm(x %*% par)
  
  n = length(y) 
  k = length(par)
  
  # Computing the gradient
  g = t(matrix(rep(phi/Phi,k),nrow=n)*x) %*% y - 
    t(matrix(rep(phi/(1-Phi),k),nrow=n)*x) %*% (1-y)
  g = -g
  
  return(g)
}

J_1 <- function(y,x,par) {
  
  Phi = pnorm(x %*% par)
  phi = dnorm(x %*% par)
  
  n = length(y) 
  k = length(par)
  
  g = matrix(rep(y*phi/Phi - (1-y)*phi/(1-Phi),k),nrow=n)*x
  
  f = t(g)%*%g
  
  return(f)
}

# create empty matrices 
J1_vec = matrix(0,num,1)
J2_vec = matrix(0,num,1)
J_2_sum = matrix(rep(0,4),2)

theta_hat_vec = matrix(0,num,k)

Wald_vec = matrix(0,num,1)

reject_Wald = matrix(0,num,1)
reject_Score = matrix(0,num,1)
reject_LR = matrix(0,num,1)

cv = qchisq(.95, df=1) 

for (it in 1:num) {
  # Data generating process
  x = cbind(matrix(1,n,1),matrix(rnorm((k-1)*n,0,1),ncol=k-1)) 	# regressors
  
  u = rnorm(n,0,1)							 	# error term
  
  y_star = x %*% theta0 + u						# latent "utility"
  
  y = ceiling(y_star/(max(abs(y_star))+0.1))				# observed outcome
  
  dat = data.frame(x,y)
  probit <- glm(y~x-1, data=dat, family = binomial(link = "probit"))
  theta_hat = probit$coefficients
  
  theta_hat_vec[it,1:k] = theta_hat
  
  result <- optim(par = theta_hat, Probit_LL, y = y, x = x, method = c("BFGS"), control = list(reltol=1e-9), hessian=TRUE)
  
  
  ###constrained model 
  result_constrained <- optim(par = theta_null[1], constrained_Probit_LL, y = y, x = x, method = c("BFGS"), control = list(reltol=1e-9), hessian=TRUE)
  
  theta_constrained = result_constrained$par 
  theta_transplant = c(theta_constrained,theta_null[2])
  
  #J_2 = (1/n)*(result$hessian)
  #V = solve(J_2)
  #V = solve((1/n)*V)
  J_2 = 1/n*solve(result$hessian)
  
  # Wald 
  temp=theta_hat[2]-theta_null[2]
  Wald = temp^2/J_2[2,2]
  Wald_vec[it]= Wald
  
  
  #temp = theta_hat[2] - theta_null[2]
  
  #Wald = t(temp)%*%V[2,2]%*%temp
  
  if (is.nan(Wald) == 0) {
    if (Wald > cv) {
      reject_Wald[it] = 1}
  }
  
  # Score
  #Score = Probit_LL_g(y,x,theta_null)
  
  #Score = t(Score)[2]%*%solve(J_1(y,x,theta_hat))[2,2]%*%Score[2]
  J1_vec[it] = J_1(y,x,theta_transplant)
  Score = Probit_LL_g(y,x,theta_transplant)
  
  Score = t(Score)%*%solve(J_1(y,x,theta_transplant))%*%Score
  
  if (is.nan(Score) == 0) {
    if (Score > cv) {
      reject_Score[it] = 1}
  }
  
  
  # Likelihood ratio
  LR = 2*(constrained_Probit_LL(y,x,theta_constrained)*n - Probit_LL(y,x,theta_hat)*n)
  
  if (is.nan(LR) == 0) {
    if (LR > cv) {
      reject_LR[it] = 1 }   
  }
}


colMeans(reject_Wald)
colMeans(reject_Score)
colMeans(reject_LR)




### Question G  (F does not require code)
#use only for wald

set.seed(1) 		# Set seed for random number generator

n = 300			# Set the sample size

num = 1000
theta0 = c(1,1) 	# True parameter value

theta_null = c(1,1) 

k = length(theta0) 

theta_start = theta0

##################WALD (The right one)

rm(list=ls()) 		# Clear workspace

set.seed(1) 		# Set seed for random number generator

n = 300 			# Set the sample size

theta0 = c(1,1) # True parameter value

theta_null = c(1,1) 

k = length(theta0) 

theta_start = theta0

# Load the log-likelihood, its derivative, and the hessian
Probit_LL <- function(y,x,par) {
  
  Phi = pnorm(x %*% par)
  phi = dnorm(x %*% par)
  
  n = length(y) 
  k = length(par)
  
  # Computing the log-likelihood
  f = sum(y*log(Phi)) + sum((1-y)*log(1-Phi))
  f = -f
  
  return(f)
}


Probit_GMM <- function(y,x,par) {
  
  n = length(y) 
  k = length(par)
  
  Phi = pnorm(x %*% par)
  
  f = (1/n^2)*(t(t(x)%*%(y-Phi))) %*% (t(x)%*%(y-Phi))
  return(f)
}

Probit_GMM_g <- function(y,x,par) {
  
  n = length(y) 
  k = length(par)
  
  Phi = pnorm(x %*% par)
  phi = dnorm(x %*% par)
  
  
  g = (2/(n^2)) * t(x) %*% (((x%*%t(x)%*%(y-Phi))) * phi)  
  g = -g
  return(g)
}

Probit_LL_g <- function (y,x,par) {
  
  n = length(y) 
  k = length(par)
  Phi = pnorm(x %*% par)
  phi = dnorm(x %*% par)
  
  
  # Computing the gradient
  g =  t(y*phi/Phi-(1-y)*phi/(1-Phi))%*%x
  g = -(1/n)*t(g)
  
  return(g)
}

Probit_LL_h <- function (y,x,par) {
  
  Phi = pnorm(x %*% par)
  phi = dnorm(x %*% par)
  
  n = length(y)
  k = length(par)
  
  # Computing the Hessian
  h = t(x)%*% (matrix(rep((1-y)*(phi/(1-Phi)*(x%*%par)-(phi^2)/(1-Phi)^2),k),nrow=n)*x) +
    t(x)%*% (matrix(rep((y)*(phi/Phi*(-x%*%par)-(phi^2)/Phi^2),k),nrow=n)*x)
  # (1/n) *( sum(y*log(Phi)) + sum((1-y)*log(1-Phi)))
  h = -h
  
  return(h)
}

Probit_Var_GMM <- function(y,x,par) {
  
  n = length(y)
  k = length(par)
  
  Phi = pnorm(x %*% par)
  phi = dnorm(x %*% par)
  
  m = matrix(rep((y-Phi),k),nrow=n)*x
  A = t(x) %*% (matrix(rep(phi,k),nrow=n)*x)
  A_inv = solve(A)
  
  var_hat = A_inv %*% (t(m) %*% m) %*% A_inv
  
  return(var_hat)
}
num = 1000			# Number of Monte Carlo iterations

theta_hat_vec = matrix(0,num,k)

J_1_inv = matrix(0,k,k)
J_2_inv = matrix(0,k,k)
J_3_inv = matrix(0,k,k)

Wald_vec = matrix(0,num,1)

reject_Wald = matrix(0,num,1)
reject_Score = matrix(0,num,1)
reject_LR = matrix(0,num,1)

cv = qchisq(.95, df=k) 

for (it in 1:num) {
  # Data generating process
  x = cbind(matrix(1,n,1),matrix(rnorm((k-1)*n,0,1),ncol=k-1)) 	# regressors
  
  u = rnorm(n,0,1)							 	# error term
  
  y_star = x %*% theta0 + u						# latent "utility"
  
  y = ceiling(y_star/(max(abs(y_star))+0.1))				# observed outcome
  
  result <- optim(par = theta_start, Probit_GMM, y = y, x = x, gr = Probit_GMM_g, method = c("BFGS"), control = list(reltol=1e-9))
  theta_hat = result$par
  
  #dat = data.frame(x,y)
  #probit <- glm(y~x-1, data=dat, family = binomial(link = "probit"))
  #theta_hat = probit$coefficients
  
  theta_hat_vec[it,1:k] = theta_hat
  
  # Wald 
  temp = theta_hat[2] - theta_null[2]
  
  Wald = t(temp)%*%solve(Probit_Var_GMM(y,x,theta_hat))[2,2]%*%temp
  
  Wald_vec[it] = Wald
  
  if (is.nan(Wald) == 0) {
    if (Wald > cv) {
      reject_Wald[it] = 1
    }
  }
  
  
  # LR
  # LR = 2*(Probit_Var_GMM(y,x,theta_null)-Probit_Var_GMM(y,x,theta_hat))
  # 
  # if (is.nan(LR) == 0) {
  #   if (LR > cv) {
  #     reject_LR[it] = 1
  #   }
  # }
  # 
}

hist(Wald_vec,30, freq=FALSE) 
curve(dchisq(x, df = k), from = 0, to = 15, add=TRUE)

mean(reject_Wald)
mean(reject_LR)



######################################à
#probit (wrong wald, use only for qlr)

set.seed(2) 		# Set seed for random number generator

n  = 300		# Set the sample size

num = 1000
theta0 = c(1,1) 	# True parameter value

theta_null = c(1,1) 

k = length(theta0) 

theta_start = theta0

# define functions
Probit_GMM <- function(y,x,par) {
  
  n = length(y) 
  k = length(par)
  
  Phi = pnorm(x %*% par)
  
  f = (1/n^2)*t(t(x)%*%(y - Phi))%*%(t(x)%*%(y - Phi))
  return(f)
}

Probit_GMM_OG <- function(y,x,par) {
  
  n = length(y) 
  k = length(par)
  Phi = pnorm(x %*% par)
  
  f = (1/n)*t(y - Phi)%*%x
  return(t(f))
}



constrained_Probit_GMM <- function(y,x,par_1) {
  par_2 = theta_null[2]
  par = c(par_1,par_2)
  n = length(y) 
  k = length(par)
  
  Phi = pnorm(x %*% par)
  
  f = (1/n^2)*t(t(x)%*%(y - Phi))%*%(t(x)%*%(y - Phi))
  return(f)
}


# Probit_GMM_g <- function(y,x,par) {

#   n = length(y) 
#   k = length(par)

#   Phi = pnorm(x %*% par)
#   phi = dnorm(x %*% par)

#   g = (2/n^2) * t(x) %*% (((x%*%t(x)%*%(y-Phi))) * phi) 
#   g = -g
#   return(g)
# }

Probit_GMM_g <- function(y,x,par) {
  
  n = length(y) 
  k = length(par)
  
  Phi = pnorm(x %*% par)
  phi = dnorm(x %*% par)
  
  f = m = t((y-Phi)) %*% x
  A = t(x) %*% (matrix(rep(phi,k),nrow=n)*x)
  
  f = -2*(A %*% t(m))/n^2
  
  return(f)
}


Probit_Var_GMM <- function(y,x,par) {
  
  n = length(y)
  k = length(par)
  
  Phi = pnorm(x %*% par)
  phi = dnorm(x %*% par)
  
  m = matrix(rep((y-Phi),k),nrow=n)*x
  A = t(x) %*% (matrix(rep(phi,k),nrow=n)*x)
  A_inv = solve(A)
  
  var_hat = A_inv %*% (t(m) %*% m) %*% A_inv
  
  return(var_hat)
}



Probit_GMM_omega <- function(y,x,par) {
  
  n = length(y)
  k = length(par)
  Phi = pnorm(x %*% par)
  phi = dnorm(x %*% par)
  
  a = -(1/n)*sum(phi)
  b = -(1/n)*sum(x[,2]*phi)
  c = -(1/n)*sum(x[,2]^2*phi)
  g = matrix(data=c(c(a,b),c(b,c)),ncol=2)
  G = solve(g)
  
  h = (1/n)*sum((y-Phi)^2)
  i = (1/n)*sum((y-Phi)^2 * x[,2])
  j = (1/n)*sum((y-Phi)^2 * x[,2]^2)
  O = matrix(data=c(c(h,i),c(i,j)), ncol=2)
  return(O)
}


# create empty matrices 
J_2_sum = matrix(rep(0,4),2)

theta_hat_vec = matrix(0,num,k)

Wald_vec = matrix(0,num,1)
LR_vec  = matrix(0,num,1)

reject_Wald = matrix(0,num,1)
reject_LR = matrix(0,num,1)

cv = qchisq(.95, df=1) 


for (it in 1:num) {
  # Data generating process
  x = cbind(matrix(1,n,1),matrix(rnorm((k-1)*n,0,1),ncol=k-1)) 	# regressors
  
  u = rnorm(n,0,1)							 	# error term
  
  y_star = x %*% theta0 + u						# latent "utility"
  
  y = ceiling(y_star/(max(abs(y_star))+0.1))				# observed outcome
  
  dat = data.frame(x,y)
  
  # GMM
  
  result_constrained <- optim(par = theta0[1], constrained_Probit_GMM, y = y, x = x, method = c("BFGS"), control = list(reltol=1e-9), hessian=TRUE)
  theta_hat_GMM_constrained = result_constrained$par
  
  result <- optim(par = theta0, Probit_GMM, y = y, x = x, method = c("BFGS"), control = list(reltol=1e-9), hessian=TRUE)
  theta_hat_GMM = result$par
  
  
  
  # Wald 
  theta_GMM_transplant = c(theta_hat_GMM_constrained,theta_null[2])
  temp = theta_hat_GMM - theta_GMM_transplant
  
  Wald = t(temp) %*% solve(Probit_Var_GMM(y,x,theta_hat_GMM)) %*% temp
  
  
  
  if (is.nan(Wald) == 0) {
    if (Wald > cv) {
      reject_Wald[it] = 1}
  }
  
  
  # likelihood ratio
  
  LR = (t(Probit_GMM_OG(y,x,theta_GMM_transplant))%*% solve(Probit_GMM_omega(y,x,theta_hat_GMM)/n) %*% Probit_GMM_OG(y,x,theta_GMM_transplant)-t(Probit_GMM_OG(y,x,theta_hat_GMM)) %*% solve(Probit_GMM_omega(y,x,theta_hat_GMM)/n)  %*% Probit_GMM_OG(y,x,theta_hat_GMM))
  LR_vec[it] = LR
  
  if (is.nan(LR) == 0) {
    if (LR > cv) {
      reject_LR[it] = 1 }   
  }
}



m = matrix(c(mean(reject_Wald),mean(reject_LR)))
rownames(m) = c("Wald","QLR")
kable(list(m), col.names = c("$\\theta_1 = 1$"), escape = F, 
      caption = "GMM Rejection Probability") %>%
  kable_styling(latex_options = "hold_position")



hist(LR_vec, col="red")
