
library(numDeriv) 	# loads the functions grad and hessian which numerically evaluate the gradient and hessian

rm(list=ls()) 		# Clear workspace

set.seed(1) 		# Set seed for random number generator

n = 300 			# Set the sample size

theta0 = c(1,1) 	# True parameter value

k = length(theta0) 

# Data generating process
x = cbind(matrix(1,n,1),matrix(rnorm((k-1)*n,0,1),ncol=k-1)) 	# regressors

u = rnorm(n,0,1)							 	# error term

y_star = x %*% theta0 + u						# latent "utility"

y = ceiling(y_star/(max(abs(y_star))+0.1))				# observed outcome

# Load the log-likelihood
Probit_LL <- function(y,x,par) {
  
  n = length(y) 
  k = length(par)
  
  Phi = pnorm(x %*% par)
 
  # Computing the log-likelihood
  f =(1/n) *( sum(y*log(Phi)) + sum((1-y)*log(1-Phi)))
  return(-f)
  
}

Probit_LL(y, x, theta0)
# optim without user-specified gradient
result_b <- optim(par = theta0, Probit_LL, y = y, x = x, method = c("BFGS"), control = list(reltol=1e-9), hessian=TRUE)
result_b$par

# Load the derivative of the log-likelihood
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
# Check if the gradient function was correctly programmed by comparing it to a numerical approximation of it
Probit_LL_g(y,x,theta0)
grad(function(u) Probit_LL(y,x,u),theta0)

# optim with user-specified gradient
result_c <- optim(par = theta0, Probit_LL, y = y, x = x, gr = Probit_LL_g, method = c("BFGS"), control = list(reltol=1e-9), hessian=TRUE)
result_c$par

# Load the objective function for NLS
Probit_NLS <- function(y,x,par) {
  
  n = length(y) 
  k = length(par)
  
  Phi = pnorm(x %*% par)
  
  f = (1/n) * (sum((y - Phi)^2))
  
  return(f)
}

# optim without user-specified gradient
result_e <- optim(par = theta0, Probit_NLS, y = y, x = x, method = c("BFGS"), control = list(reltol=1e-30), hessian=TRUE)
result_e$par

# Load the derivative of the NLS objective function
Probit_NLS_g <- function(y,x,par) {
  
  n = length(y) 
  k = length(par)
  
  Phi = pnorm(x %*% par)
  phi = dnorm(x %*% par)
  
  g = (2/n) * t(-phi*(y-Phi)) %*% x
  
  return(g)
}

# Check if the gradient function was correctly programmed by comparing it to a numerical approximation of it
Probit_NLS_g(y,x,theta0)

grad(function(u) Probit_NLS(y,x,u),theta0)

# optim with user-specified gradient
result_f <- optim(par = theta0, Probit_NLS, y = y, x = x, gr = Probit_NLS_g, method = c("BFGS"), control = list(reltol=1e-9), hessian=TRUE)
result_f$par

# Load the objective function for MM
Probit_GMM <- function(y,x,par) {
  
  n = length(y) 
  k = length(par)
  
  Phi = pnorm(x %*% par)
  
  f = (1/n^2)*(t(t(x)%*%(y-Phi))) %*% (t(x)%*%(y-Phi))
  return(f)
}


Probit_GMM(y,x,theta0) 


# optim without user-specified gradient
result_h <- optim(par = theta0, Probit_GMM, y = y, x = x, method = c("BFGS"), control = list(reltol=1e-9), hessian=TRUE)
result_h$par

# Load the derivative of the GMM objective function
Probit_GMM_g <- function(y,x,par) {
  
  n = length(y) 
  k = length(par)
  
  Phi = pnorm(x %*% par)
  phi = dnorm(x %*% par)
  
  
  g = (2/(n^2)) * t(x) %*% (((x%*%t(x)%*%(y-Phi))) * phi)  
  g = -g
  return(g)
}

# Check if the gradient function was correctly programmed by comparing it to a numerical approximation of it
Probit_GMM_g(y,x,theta0)
grad(function(u) Probit_GMM(y,x,u),theta0)

# optim with user-specified gradient
result_i <- optim(par = theta0, Probit_GMM, y = y, x = x, gr = Probit_GMM_g, method = c("BFGS"), control = list(reltol=1e-9), hessian=TRUE)
result_i$par

############################ MONTECARLO PART #######################################

library(tictoc)

tic()

# rm(list=ls()) 		# Clear workspace

set.seed(1) 		# Set seed for random number generator

n = 300 			# Set the sample size

theta0 = c(1,1) 	# True parameter value

k = length(theta0) 

source("Probit_LL.R")
source("Probit_NLS.R")
source("Probit_GMM.R")
# source("Probit_J_1.R"
source("Probit_Sigma_NLS.R")
source("Probit_Var_GMM.R")

Probit_J_1 <- function(y,x,par) {
  n = length(y)
  k = length(par)
  
  score = Probit_LL_g(y,x,par) 
  f = score %*% t(score)
  
  return(f)
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

J_2_ML <- function (y,x,par) {

   # Computing the Hessian
  result_f <- optim(par = par, Probit_LL, y = y, x = x, gr = Probit_LL_g, method = c("BFGS"),
                    control = list(reltol=1e-9), hessian=TRUE)
  hessian <- result_f$hessian

  return(hessian)
}

J_2_ML(y,x,theta0)

num = 1000			# Number of Monte Carlo iterations

theta_hat_ML_vec = matrix(0,num,k)
theta_hat_NLS_vec = matrix(0,num,k)
theta_hat_GMM_vec = matrix(0,num,k)

inside_N_ML = rep(0,num)
inside_N_NLS = rep(0,num)
inside_N_GMM = rep(0,num)

J_1_sum = matrix(rep(0,4),2)
J_2_sum = matrix(rep(0,4),2)
Var_hat_NLS_sum = matrix(rep(0,4),2)
Var_hat_GMM_sum = matrix(rep(0,4),2)

epsilon = 0.1

for (it in 1:num) {
  # Data generating process
  x = cbind(matrix(1,n,1),matrix(rnorm((k-1)*n,0,1),ncol=k-1)) 	# regressors
  
  u = rnorm(n,0,1)							 	# error term
  
  y_star = x %*% theta0 + u						# latent "utility"
  
  y = ceiling(y_star/(max(abs(y_star))+0.1))				# observed outcome
  
  # ML
  result <- optim(par = theta0, Probit_LL, y = y, x = x, method = c("BFGS"), control = list(reltol=1e-9), hessian=TRUE)
  theta_hat_ML = result$par
  
  theta_hat_ML_vec[it,1:k] = theta_hat_ML
  
  if (sqrt(sum((theta_hat_ML - theta0)^2)) < epsilon) {
    inside_N_ML[it] = 1
  }
  
  J_1_sum = J_1_sum + Probit_J_1(y,x,theta_hat_ML) 
  J_2_sum = J_2_sum + J_2_ML(y,x,theta_hat_ML) # it just is the Hessian
  
  # NLS
  result <- optim(par = theta0, Probit_NLS, y = y, x = x, method = c("BFGS"), control = list(reltol=1e-9), hessian=TRUE)
  theta_hat_NLS = result$par
  
  theta_hat_NLS_vec[it,1:k] = theta_hat_NLS
  
  if (sqrt(sum((theta_hat_NLS - theta0)^2)) < epsilon) {
    inside_N_NLS[it] = 1
  }
  
  #	Var_hat_NLS_sum = Var_hat_NLS_sum + ...
  
  # GMM
  result <- optim(par = theta0, Probit_GMM, y = y, x = x, method = c("BFGS"), control = list(reltol=1e-9), hessian=TRUE)
  theta_hat_GMM = result$par
  
  theta_hat_GMM_vec[it,1:k] = theta_hat_GMM
  
  if (sqrt(sum((theta_hat_GMM - theta0)^2)) < epsilon) {
    inside_N_GMM[it] = 1
  }
  
  #	Var_hat_GMM_sum = Var_hat_GMM_sum + ...
}

# Averages - not asked for
colMeans(theta_hat_ML_vec)
colMeans(theta_hat_NLS_vec)
colMeans(theta_hat_GMM_vec)

# Variances
var(theta_hat_ML_vec)
var(theta_hat_NLS_vec)
var(theta_hat_GMM_vec)

# (Estimated) Probability that theta_hat lies inside a neigborhood around theta_0
mean(inside_N_ML)
mean(inside_N_NLS)
mean(inside_N_GMM)














library("ggplot2")
library("latex2exp")
gmm_to_plot <- as.data.frame(theta_hat_GMM_vec[1:num,2]) 
colnames(gmm_to_plot)[1]  <- "theta_2" 
nls_to_plot <- as.data.frame(theta_hat_NLS_vec[1:num,2]) 
colnames(nls_to_plot)[1]  <- "theta_2" 
ML_to_plot <- as.data.frame(theta_hat_ML_vec[1:num,2]) 
colnames(ML_to_plot)[1]  <- "theta_2" 

gmm <- ggplot(gmm_to_plot, aes(x = theta_2)) + # we use the second column as it is the one for
                                                     # theta_1, the coefficient of x
  geom_histogram(fill = "#99d8c9") + 
  stat_function(fun = dnorm, args = list(mean(gmm_to_plot$theta_2), sd = sd(gmm_to_plot$theta_2)), linewidth = 1) +
  geom_vline(xintercept = 1, color = "#2ca25f", size = 1) +
    theme_bw() +
    labs(title = TeX(r"(Density of $\theta_2$ estimated with GMM)"),
      x = TeX(r"($\theta_2$)"),
      y = "Density")
  
nls <- ggplot(nls_to_plot, aes(x = theta_2)) + # we use the second column as it is the one for
  # theta_1, the coefficient of x
  geom_histogram(fill = "#99d8c9") + 
  stat_function(fun = dnorm, args = list(mean(nls_to_plot$theta_2), sd = sd(nls_to_plot$theta_2)), linewidth = 1) +
  geom_vline(xintercept = 1, color = "#2ca25f", size = 1) +
  theme_bw() +
  labs(title = TeX(r"(Density of $\theta_2$ estimated with NLS)"),
       x = TeX(r"($\theta_2$)"),
       y = "Density")

ml <- ggplot(ml_to_plot, aes(x = theta_2)) + # we use the second column as it is the one for
  # theta_1, the coefficient of x
  geom_histogram(fill = "#99d8c9") + 
  stat_function(fun = dnorm, args = list(mean(ml_to_plot$theta_2), sd = sd(ml_to_plot$theta_2)), linewidth = 1) +
  geom_vline(xintercept = 1, color = "#2ca25f", size = 1) +
  theme_bw() +
  labs(title = TeX(r"(Density of $\theta_2$ estimated with ML)"),
       x = TeX(r"($\theta_2$)"),
       y = "Density")

gridExtra::grid.arrange(ml, nls, gmm, nrow = 1)

J_1_sum/num
J_2_sum/num
Var_hat_NLS_sum/num
Var_hat_GMM_sum/num

toc()


##### question n ####
Probit_NLS_sigma <- function(y,x,par) {
  
  n = length(y) 
  k = length(par)
  
  Phi = pnorm(x %*% par)
  phi = dnorm(x %*% par)
  
  sigma= t((2/n) * t(-phi*(y-Phi)) %% x) %% ((2/n) * t(-phi*(y-Phi)) %*% x)
  result_f <- optim(par = theta0, Probit_NLS, y = y, x = x, gr = Probit_NLS_g, method = c("BFGS"), control = list(reltol=1e-9), hessian=TRUE)
  hessian <- result_f$hessian
  
  g= ((hessian)^-1 *sigma *(hessian)^-1)
  
  return(g)
}

Probit_NLS_sigma(y=y, x=x, par=theta0)

############### question r ##############
set.seed(1) 		# Set seed for random number generator

n = 300 			# Set the sample size

theta0 = c(1,1) 	# True parameter value

num = 1000

Asvar <-  matrix(0,num)

for (it in 1:num) {
  # Data generating process
  x = cbind(matrix(1,n,1),matrix(rnorm((k-1)*n,0,1),ncol=k-1)) 	# regressors
  
  u = rnorm(n,0,1)							 	# error term
  
  y_star = x %*% theta0 + u						# latent "utility"
  
  y = ceiling(y_star/(max(abs(y_star))+0.1))		

  result <- optim(par = theta0, Probit_GMM, y = y, x = x, method = c("BFGS"), control = list(reltol=1e-9), hessian=TRUE)
  theta_hat_GMM = result$par
  
Phi = pnorm(x %*% theta_hat_GMM)
phi = dnorm(x %*% theta_hat_GMM)

G <- (1/n)* sum((2/n) * t(x) %*% (((x%*%t(x)%*%(y-Phi))) * phi) )
G_inverse <- G^(-1)
G_trans <- t(G)

Omega <- (1/n)* sum((1/n^2)*(t(t(x)%*%(y-Phi))) %*% (t(x)%*%(y-Phi)) * (1/n^2)*(t(t(x)%*%(y-Phi))) %*% (t(x)%*%(y-Phi)))

Asy <- G_trans * Omega * G

Asvar[it,] = (G_inverse * Omega * G_inverse)/n
}

Asvar <- as.data.frame(Asvar)
options(scipen=0)
mean(Asvar$V1)

# now with theta0
Asvar_theta0 <- as.data.frame(matrix(0,num))

for (it in 1:num) {
  # Data generating process
  x = cbind(matrix(1,n,1),matrix(rnorm((k-1)*n,0,1),ncol=k-1)) 	# regressors
  
  u = rnorm(n,0,1)							 	# error term
  
  y_star = x %*% theta0 + u						# latent "utility"
  
  y = ceiling(y_star/(max(abs(y_star))+0.1))		
 
  Phi = pnorm(x %*% theta0)
  phi = dnorm(x %*% theta0)
  
  G <- (1/n)* sum((2/n) * t(x) %*% (((x%*%t(x)%*%(y-Phi))) * phi) )
  G_inverse <- G^(-1)
  
  Omega <- (1/n)* sum((1/n^2)*(t(t(x)%*%(y-Phi))) %*% (t(x)%*%(y-Phi)) * (1/n^2)*(t(t(x)%*%(y-Phi))) %*% (t(x)%*%(y-Phi)))
  
  Asvar_theta0[it,] = (G_inverse * Omega * G_inverse)/n
}

mean(Asvar_theta0$V1)

