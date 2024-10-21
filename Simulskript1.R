
library(splines)
library(stats4)
library(VGAM)


### Generating random data process.

set.seed(123)
#1) create a dataset with X following a normal but different patterns for the error.
datagen <- function(n, beta, rcdf, meanX, sdX, M=100) {
  d <- length(meanX)
  datasets <- vector("list", M)  # Create a list to store M datasets
  for (m in 1:M) {  # Loop over the number of replicates
    x <- matrix(data = NA, nrow = n, ncol = d)  # Design matrix
    for (j in 1:d)  x[,j] <- rnorm(n,mean=meanX[j],sd=sdX[j])
    eps <- rcdf(n)  # Generate epsilon
    y <- x%*%beta + eps  # Response variable
    datasets[[m]] <- list(Y = y, X = x, B = beta, E = eps)  # Store dataset in the list
  }
  return(datasets)  # Return the list of datasets
}


### Formula for criterion

# Differentiable form
critbis <- function(beta, data, F_eps) {
  n <- length(data$Y)
  residuals <- data$Y - data$X%*%beta  # This gives a vector of length n
  F_eps_residuals <- F_eps(as.vector(residuals))  # Ensure it’s a vector
  p2 <- mean(F_eps_residuals^2)
  p3 <- (1 / (n^2)) * sum(sapply(1:n, function(i) {
    sum(sapply(1:n, function(j) {
      if (i != j) {
        return(max(F_eps_residuals[i], F_eps_residuals[j]))
      } else {
        return(0)
      }
    }))
  }))
  p4 <- (1 / (n^2)) * sum(F_eps_residuals)
  result <- 1/3 + p2 - p3 - p4
  return(result)
}


### Gradient for the criterion.
# Implement the gradient. Returns a vector in R^d and the corresponding L-2 norm.
grad <- function(data, beta, F_eps,f_eps) {
  n <- length(data$Y)
  d <- length(beta)
  residuals <- data$Y - data$X%*%beta
  f_residuals <- f_eps(residuals) #returns a nx1 matrix
  F_residuals <- F_eps(residuals) #returns a nx1 matrix
  term1_weights <- f_residuals * (1 - 2 * F_residuals) # returns a 1xn matrix.
  term1 <- (t(term1_weights)%*%data$X) / n #1xd matrix.
  term2 <- 0
  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j) {
        r_i <- residuals[i]
        r_j <- residuals[j]
        indicator <- as.numeric(r_i > r_j)
        term2 <- term2 + indicator * (f_eps(r_i) * data$X[i,] - f_eps(r_j) * data$X[i,])
      }
    }
  }
  term2 <- term2 / n^2
  result <- term1 + term2
  norm <- sqrt(sum(result^2))
  return(list(G=result,N=norm)) #$G is 1xd matrix.
}


# Gradient descent algorithm.

# Function to build OLS.
ols <- function(data) {
  return(solve(t(data$X)%*%data$X)%*%t(data$X)%*%data$Y)
}

# Gradient with Backtracking Line Search
# Since the initial values are crucial for the performance of the algorithm, take as initial value classic OLS.
# Try to improve the speed 
backtracking<-function(data, cdf, pmf, tol = 1e-3, a = .001, b = .9,max_iteration = 1000){
  k = 1
  beta <- ols(data) #returns a 2x1 matrix.
  betaplus <- beta - a * t(grad(data, beta, cdf,pmf)$G/grad(data, beta, cdf,pmf)$N)
  while (sqrt(sum((beta - betaplus)^2)) > tol && k < max_iteration){
    a <- b * a
    beta <- betaplus
    betaplus <- beta - a * t(grad(data, beta, cdf,pmf)$G/grad(data, beta, cdf,pmf)$N)
    k <- k + 1
  }
  return(list(iteration = k - 1, res = betaplus, compa = cbind(data$B,betaplus,ols(data))))
}


### Simulations for || \beta^*-\beta_0||, || OLS-\beta_0||

# Create a function which returns as a list both || \beta^*-\beta_0|| and || OLS-\beta_0||.
normbeta <- function(data, F_eps, f_eps) {
  return(list(norm1=sqrt(sum((backtracking(data,F_eps, f_eps)$res-data$B)^2)),norm2=sqrt(sum((data$B-ols(data))^2))))
}

# MLE algorithm.
# Define the negative log-likelihood function for Laplace-distributed errors
neg_lllaplace <- function(params, X, Y,b=1) {
  beta <- params # Extract regression coefficients (beta)
  
  # Calculate residuals (errors) from the linear model
  residuals <- Y - X %*% beta
  
  # Negative log-likelihood for Laplace distribution
  n <- length(Y)
  neg_ll <- n * log(2 * b) + sum(abs(residuals) / b)
  
  return(neg_ll)
}

dat1 <- datagen(1500,c(-3,2,0.3,0.8,-6),rlaplace,c(0.4,-4,1,2.3,1.75),c(0.23,3,0.75,1,2))

rescoco <- list()
for (i in 1:100) {
  rescoco[[i]] <- normbeta(dat1[[i]],plaplace,dlaplace)
  print(i)
}

file_path <- "~/polybox - Justine Leclerc (jleclerc@student.ethz.ch)@polybox.ethz.ch (2)/THESIS !/res_simul.xlsx"
wb <- loadWorkbook(file_path)
addWorksheet(wb, "15005laplace")
writeData(wb, sheet = "15005laplace", rescoco)
saveWorkbook(wb, file_path, overwrite = TRUE)

library(optimx)
initial_params <- c(rep(0, ncol(dat1[[1]]$X)))
# Optimize the negative log-likelihood

mle <- lapply(dat1, function(data) optimx(par=initial_params, fn=neg_lllaplace, X=data$X, Y=data$Y, method="BFGS"))

mlebis <- matrix(rep(0,500),nrow=100)
for (i in 1:100){
  mlebis[i,] <-c(mle[[i]]$p1,mle[[i]]$p2,mle[[i]]$p3)
}

normmle <- apply(mlebis,1,function(data) sqrt(sum(data-c(-3,2,0.3,0.8,-6))^2))

file_path <- "~/polybox - Justine Leclerc (jleclerc@student.ethz.ch)@polybox.ethz.ch (2)/THESIS !/res_simul.xlsx"
wb <- loadWorkbook(file_path)
addWorksheet(wb, "15005laplacebis")
writeData(wb, sheet = "15005laplacebis", normmle)
saveWorkbook(wb, file_path, overwrite = TRUE)

dat1 <- datagen(500,c(-3,2,0.3,0.8,-6),rlaplace,c(0.4,-4,1,2.3,1.75),c(0.23,3,0.75,1,2))

rescoco <- list()
for (i in 1:100) {
  rescoco[[i]] <- normbeta(dat1[[i]],plaplace,dlaplace)
  print(i)
}

file_path <- "~/polybox - Justine Leclerc (jleclerc@student.ethz.ch)@polybox.ethz.ch (2)/THESIS !/res_simul.xlsx"
wb <- loadWorkbook(file_path)
addWorksheet(wb, "5005laplace")
writeData(wb, sheet = "5005laplace", rescoco)
saveWorkbook(wb, file_path, overwrite = TRUE)

library(optimx)
initial_params <- c(rep(0, ncol(dat1[[1]]$X)))
# Optimize the negative log-likelihood

mle <- lapply(dat1, function(data) optimx(par=initial_params, fn=neg_lllaplace, X=data$X, Y=data$Y, method="BFGS"))

mlebis <- matrix(rep(0,500),nrow=100)
for (i in 1:100){
  mlebis[i,] <-c(mle[[i]]$p1,mle[[i]]$p2,mle[[i]]$p3)
}

normmle <- apply(mlebis,1,function(data) sqrt(sum(data-c(-3,2,0.3,0.8,-6))^2))

file_path <- "~/polybox - Justine Leclerc (jleclerc@student.ethz.ch)@polybox.ethz.ch (2)/THESIS !/res_simul.xlsx"
wb <- loadWorkbook(file_path)
addWorksheet(wb, "5005laplacebis")
writeData(wb, sheet = "5005laplacebis", normmle)
saveWorkbook(wb, file_path, overwrite = TRUE)

dat1 <- datagen(1500,c(-3,2,0.3,0.8,-6),function(x) rlaplace(x, scale=4),c(0.4,-4,1,2.3,1.75),c(0.23,3,0.75,1,2))

rescoco <- list()
for (i in 1:100) {
  rescoco[[i]] <- normbeta(dat1[[i]],function(x) plaplace(x, scale=4),function(x) dlaplace(x, scale=4))
  print(i)
}

file_path <- "~/polybox - Justine Leclerc (jleclerc@student.ethz.ch)@polybox.ethz.ch (2)/THESIS !/res_simul.xlsx"
wb <- loadWorkbook(file_path)
addWorksheet(wb, "5005laplacesc")
writeData(wb, sheet = "5005laplacesc", rescoco)
saveWorkbook(wb, file_path, overwrite = TRUE)

library(optimx)
initial_params <- c(rep(0, ncol(dat1[[1]]$X)))
# Optimize the negative log-likelihood

mle <- lapply(dat1, function(data) optimx(par=initial_params, fn=neg_lllaplace, X=data$X, Y=data$Y, b=4, method="BFGS"))

mlebis <- matrix(rep(0,500),nrow=100)
for (i in 1:100){
  mlebis[i,] <-c(mle[[i]]$p1,mle[[i]]$p2,mle[[i]]$p3)
}

normmle <- apply(mlebis,1,function(data) sqrt(sum(data-c(-3,2,0.3,0.8,-6))^2))

file_path <- "~/polybox - Justine Leclerc (jleclerc@student.ethz.ch)@polybox.ethz.ch (2)/THESIS !/res_simul.xlsx"
wb <- loadWorkbook(file_path)
addWorksheet(wb, "5005laplacebissc")
writeData(wb, sheet = "5005laplacebissc", normmle)
saveWorkbook(wb, file_path, overwrite = TRUE)