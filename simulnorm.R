library(splines)
library(stats4)
library(VGAM)
library(openxlsx)
library(optimx)


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

# Function to build OLS. => retourne une dx1 matrice.
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
} #compa = matrice dx3.


### Simulations for || \beta^*-\beta_0||, || OLS-\beta_0||

# Create a function which returns as a list both || \beta^*-\beta_0|| and || OLS-\beta_0||.
normbeta <- function(data, F_eps, f_eps) {
  return(list(norm1=sqrt(sum((backtracking(data,F_eps, f_eps)$res-data$B)^2)),norm2=sqrt(sum((data$B-ols(data))^2)),bias1=mean(backtracking(data,F_eps, f_eps)$res)-data$B, bias2=mean(ols(data))-data$B, var1=var(backtracking(data,F_eps, f_eps)$res), var2=var(ols(data))))
}


datares <- function(data, F_eps, f_eps){
  olse <- ols(data)
  crit <- backtracking(data, F_eps, f_eps)$res
  return(cbind.data.frame(olse,crit))
}

trainnorm <- datagen(500, c(0.95,-0.16,1.23,0.37), rnorm, c(3,-1,2,0.5), c(1,1,1,1), M=100)

res5004 <- lapply(trainnorm,function(x) datares(x, pnorm, dnorm))
write.csv(res5004 , "res5004.csv")
system('git add res5004.csv')
system('git commit -m "ajout res"')
system('git push origin main ')

trainnorm <- datagen(200, c(0.95,-0.16,1.23,0.37), rnorm, c(3,-1,2,0.5), c(1,1,1,1), M=100)

res2004 <- lapply(trainnorm,function(x) datares(x, pnorm, dnorm))
write.csv(res2004 , "res2004.csv")
system('git add res2004.csv')
system('git commit -m "ajout res"')
system('git push origin main ')

trainnorm <- datagen(500, c(0.95,-0.16), rnorm, c(3,-1), c(1,1), M=100)

res5002 <- lapply(trainnorm,function(x) datares(x, pnorm, dnorm))
write.csv(res5002 , "res5002.csv")
system('git add res5002.csv')
system('git commit -m "ajout res"')
system('git push origin main ')

# Manip. des résultats.
library(readr)
res2004 <- read_csv("Desktop/Thesis GIT/res2004.csv")
res5002 <- read_csv("Desktop/Thesis GIT/res5002.csv")
res5004 <- read_csv("Desktop/Thesis GIT/res5004.csv")

beta <- c(0.95,-0.16,1.23,0.37)

res2004olse <- as.data.frame(t(as.matrix(res2004 %>% select(starts_with("olse")))))
res2004crit <- as.data.frame(t(as.matrix(res2004 %>% select(starts_with("crit")))))
res5002olse <- as.data.frame(t(as.matrix(res5002 %>% select(starts_with("olse")))))
res5002crit <- as.data.frame(t(as.matrix(res5002 %>% select(starts_with("crit")))))
res5004olse <- as.data.frame(t(as.matrix(res5004 %>% select(starts_with("olse")))))
res5004crit <- as.data.frame(t(as.matrix(res5004 %>% select(starts_with("crit")))))

mse <- function(trueb,data){
  return((1/nrow(data))*sum(apply(data,1,function(x) (x-trueb)^2)))
}

mse(beta,res2004olse)
mse(beta, res2004crit)
