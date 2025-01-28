library(splines)
library(stats4)
library(VGAM)
library(openxlsx)
library(optimx)
library(dplyr)
library(readr)

# Functions for data treatment.
mse <- function(trueb,data){
  return((1/nrow(data))*sum(apply(data,1,function(x) (x-trueb)^2)))
}

dataresbis <- function(data,beta){
  mse <- mse(beta[1:ncol(data)],data)
  var <- apply(data,2,var)
  mean <- apply(data,2,mean)
  return(list(mse=mse, var=var, mean=mean))
}

# Fisher matrix normal.
cltnorm1 <- matrix(nrow=100,ncol=4)
for (i in 1:100) {
  for (j in 1:4) {
  cltnorm1[i,j]<- sqrt(500)*sqrt(xxtnorm1[[i]][j,j])*res2004crit[i,j ]
  }
}

cltnorm2 <- matrix(nrow=100,ncol=4)
for (i in 1:100) {
  for (j in 1:4) {
    cltnorm2[i,j]<- sqrt(500)*sqrt(xxtnorm2[[i]][j,j])*res2004crit[i,j ]
  }
}

cltnorm3 <- matrix(nrow=100,ncol=4)
for (i in 1:100) {
  for (j in 1:2) {
    cltnorm3[i,j]<- sqrt(500)*sqrt(xxtnorm3[[i]][j,j])*res2004crit[i,j ]
  }
}

# Manip. des résultats SIMULNORM.
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



filenorm<-list(res2004olse,
            res2004crit,
            res5004olse,
            res5004crit,
            res5002olse,
            res5002crit
            )

syntnorm<-lapply(filenorm,function(x) dataresbis(x,beta))
msenorm <- as.data.frame(cbind(rep(c("OLSE","beta_hat"),3),as.numeric(c(50,50,125,125,250,250)),as.numeric(unlist(lapply(syntnorm, function(x) x$mse)))))
msenorm$V1 <- as.factor(msenorm$V1)
msenorm$V2 <- as.numeric(msenorm$V2)
msenorm$V3 <- as.numeric(msenorm$V3)
plot(msenorm[,2],msenorm[,3], xlab="Number of observation per dimension",ylab="MSE", col=msenorm[,1],type="b", main="Gaussian errors N(0,1)" )
legend(legend =c("OLSE","beta_hat"), col=mselp[,1],x="topright", pch=1 )


apply(res5002olse,2,var)/apply(res5002crit,2,var)
apply(res5004olse,2,var)/apply(res5004crit,2,var)

plot(res2004crit[,1],res2004olse[,1])
abline(a=0,b=1)

plot(res5004crit[,1],res5004olse[,1])
abline(a=0,b=1)

lm(res2004crit[,1]~res2004olse[,1])
lm(res5004crit[,1]~res5004olse[,1])
plot(lm(res2004crit[,1]~res5004olse[,1]))

hist(res5004crit[,1])
lines(density(res5004crit[,1]))

# Manip. des résultats SIMULLAPLACE.
betalp <- c(0.95,-0.16,1.23,0.37)
xlp <- cbind(c(3,-1,2,0.5), c(1,1,1,1))

reslp2004 <- read_csv("Desktop/Thesis GIT/reslp2004.csv")
reslp5002 <- read_csv("Desktop/Thesis GIT/reslp5002.csv")
reslp5004 <- read_csv("Desktop/Thesis GIT/reslp5004.csv")

reslp2004olse <- as.data.frame(t(as.matrix(reslp2004 %>% select(starts_with("olse")))))
reslp2004crit <- as.data.frame(t(as.matrix(reslp2004 %>% select(starts_with("crit")))))
reslp5002olse <- as.data.frame(t(as.matrix(reslp5002 %>% select(starts_with("olse")))))
reslp5002crit <- as.data.frame(t(as.matrix(reslp5002 %>% select(starts_with("crit")))))
reslp5004olse <- as.data.frame(t(as.matrix(reslp5004 %>% select(starts_with("olse")))))
reslp5004crit <- as.data.frame(t(as.matrix(reslp5004 %>% select(starts_with("crit")))))
reslp5002mle <- as.data.frame(t(as.matrix(reslp5002 %>% select(starts_with("BFGS")))))
reslp5004mle <- as.data.frame(t(as.matrix(reslp5004 %>% select(starts_with("BFGS")))))
reslp2004mle <- as.data.frame(t(as.matrix(reslp2004 %>% select(starts_with("BFGS")))))

filelp<-list(reslp2004olse,
               reslp2004crit,
             reslp2004mle,
               reslp5004olse,
               reslp5004crit, reslp5004mle, reslp5002olse,
             reslp5002crit, reslp5002mle)

syntlp<-lapply(filelp,function(x) dataresbis(x,betalp))
mselp <- as.data.frame(cbind(rep(c("OLSE","beta_hat","MLE"),3),as.numeric(rep(c(50,125,250),each=3)),as.numeric(unlist(lapply(syntlp, function(x) x$mse)))))
mselp$V1 <- as.factor(mselp$V1)
mselp$V2 <- as.numeric(mselp$V2)
mselp$V3 <- as.numeric(mselp$V3)

fisherlp <- solve(diag(nrow=4,ncol = 4) + c(3,-1,2,0.5)%*%t(c(3,-1,2,0.5)))
  
par(mfrow=c(2,2))
hist(reslp2004mle$V2, probability = TRUE)
lines(density(reslp2004mle$V2))
hist(reslp2004olse$V2, probability = TRUE)
lines(density(reslp2004olse$V2))
hist(reslp2004crit$V2, probability = TRUE)
lines(density(reslp2004crit$V2))

plot(density(reslp2004mle$V2), col="red")
lines(density(reslp2004olse$V2), col="blue")
lines(density(reslp2004crit$V2), col="green")

boxplot(cbind(reslp2004mle$V2,reslp2004olse$V2,reslp2004crit$V2))

#lp bis
reslpbis2004 <- read_csv("Desktop/Thesis GIT/reslpbis2004.csv")
reslpbis5002 <- read_csv("Desktop/Thesis GIT/reslpbis5002.csv")
reslpbis5004 <- read_csv("Desktop/Thesis GIT/reslpbis5004.csv")

reslpbis2004olse <- as.data.frame(t(as.matrix(reslpbis2004 %>% select(starts_with("olse")))))
reslpbis2004crit <- as.data.frame(t(as.matrix(reslpbis2004 %>% select(starts_with("crit")))))
reslpbis5002olse <- as.data.frame(t(as.matrix(reslpbis5002 %>% select(starts_with("olse")))))
reslpbis5002crit <- as.data.frame(t(as.matrix(reslpbis5002 %>% select(starts_with("crit")))))
reslpbis5004olse <- as.data.frame(t(as.matrix(reslpbis5004 %>% select(starts_with("olse")))))
reslpbis5004crit <- as.data.frame(t(as.matrix(reslpbis5004 %>% select(starts_with("crit")))))
reslpbis5002mle <- as.data.frame(t(as.matrix(reslpbis5002 %>% select(starts_with("BFGS")))))
reslpbis5004mle <- as.data.frame(t(as.matrix(reslpbis5004 %>% select(starts_with("BFGS")))))
reslpbis2004mle <- as.data.frame(t(as.matrix(reslpbis2004 %>% select(starts_with("BFGS")))))

filelpbis <- list(reslpbis2004olse,
     reslpbis2004crit,
     reslpbis2004mle,
     reslpbis5004olse,
     reslpbis5004crit, reslpbis5004mle, reslpbis5002olse,
     reslpbis5002crit, reslpbis5002mle)

fisherlpbis <- fisherlp*16

syntlpbis<-lapply(filelpbis,function(x) dataresbis(x,betalp))
mselpbis <- as.data.frame(cbind(rep(c("OLSE","beta_hat","MLE"),3),as.numeric(rep(c(50,125,250),each=3)),as.numeric(unlist(lapply(syntlpbis, function(x) x$mse)))))
mselpbis$V1 <- as.factor(mselpbis$V1)
mselpbis$V2 <- as.numeric(mselpbis$V2)
mselpbis$V3 <- as.numeric(mselpbis$V3)

# Calculate global x and y limits
x_range <- range(c(mselp[,2], mselpbis[,2]))
y_range <- range(c(mselp[,3], mselpbis[,3]))

# Plot the first dataset with defined ranges
plot(mselp[,2], mselp[,3], 
     xlab = "Number of observation per dimension", 
     ylab = "MSE", 
     col = mselp[,1], 
     type = "b", 
     main = "Laplace errors L(1) and L(4)",
     xlim = x_range, 
     ylim = y_range, 
     pch = 16) # Solid points for L(1)

# Add the second dataset, keeping colors consistent with mselp[,1]
points(mselpbis[,2], mselpbis[,3], 
       col = mselp[,1], # Use the same color from mselp[,1]
       pch = 17,        # Use a different point style for L(4)
       type = "b")      # Add lines and points for consistency

# Add a legend to differentiate between datasets
legend("topright", 
       legend = c("L(1)", "L(4)","OLSE","beta_hat","MLE"), 
       col = c("black", "black", mselp[,1]), # Use unique colors from mselp[,1]
       pch = c(16, 17,15,15,15),         # Matching point styles
       )                 # Add line styles to legend


par(mfrow=c(2,2))
hist(reslpbis2004mle$V2, probability = TRUE)
lines(density(reslpbis2004mle$V2))
hist(reslpbis2004olse$V2, probability = TRUE)
lines(density(reslpbis2004olse$V2))
hist(reslpbis2004crit$V2, probability = TRUE)
lines(density(reslpbis2004crit$V2))

plot(density(reslpbis2004mle$V2), col="red")
lines(density(reslpbis2004olse$V2), col="blue")
lines(density(reslpbis2004crit$V2), col="green")

boxplot(cbind(reslpbis2004mle$V2,reslpbis2004olse$V2,reslpbis2004crit$V2))

# Manip. des résultats SIMULT.
betat <- c(0.95,-0.16,1.23,0.37)
xt <- cbind(c(3,-1,2,0.5), c(1,1,1,1))

rest2004 <- read_csv("Desktop/Thesis GIT/rest2004.csv")
rest5002 <- read_csv("Desktop/Thesis GIT/rest5002.csv")
rest5004 <- read_csv("Desktop/Thesis GIT/rest5004.csv")

rest2004olse <- as.data.frame(t(as.matrix(rest2004 %>% select(starts_with("olse")))))
rest2004crit <- as.data.frame(t(as.matrix(rest2004 %>% select(starts_with("crit")))))
rest5002olse <- as.data.frame(t(as.matrix(rest5002 %>% select(starts_with("olse")))))
rest5002crit <- as.data.frame(t(as.matrix(rest5002 %>% select(starts_with("crit")))))
rest5004olse <- as.data.frame(t(as.matrix(rest5004 %>% select(starts_with("olse")))))
rest5004crit <- as.data.frame(t(as.matrix(rest5004 %>% select(starts_with("crit")))))
rest5002mle <- as.data.frame(t(as.matrix(rest5002 %>% select(starts_with("BFGS")))))
rest5004mle <- as.data.frame(t(as.matrix(rest5004 %>% select(starts_with("BFGS")))))
rest2004mle <- as.data.frame(t(as.matrix(rest2004 %>% select(starts_with("BFGS")))))

filet<-list(rest2004olse,
             rest2004crit,
            rest2004mle,
             rest5004olse,
             rest5004crit, rest5004mle,rest5002olse,
            rest5002crit,rest5002mle)

syntt<-lapply(filet,function(x) dataresbis(x,betat))
mset <- as.data.frame(cbind(rep(c("OLSE","beta_hat","MLE"),3),as.numeric(rep(c(50,125,250),each=3)),as.numeric(unlist(lapply(syntt, function(x) x$mse)))))
mset$V1 <- as.factor(mset$V1)
mset$V2 <- as.numeric(mset$V2)
mset$V3 <- as.numeric(mset$V3)

plot(mset[,2],mset[,3], xlab="Number of observation per dimension",ylab="MSE", col=mselp[,1],type="p", main="T-distribution errors t(1)")
legend(legend =c("OLSE","beta_hat","MLE"), col=mset[,1],x="topright", pch=1 )

# Calculate global x and y limits
x_range <- range(c(mset[,2], msetbis[,2]))
y_range <- range(c(mset[,3], msetbis[,3]))

# Plot the first dataset with defined ranges
plot(mset[,2], mset[,3], 
     xlab = "Number of observation per dimension", 
     ylab = "MSE", 
     col = mset[,1], 
     type = "b", 
     main = "t-errors df=1 and df=4",
     xlim = x_range, 
     ylim = y_range, 
     pch = 16) # Solid points for L(1)

# Add the second dataset, keeping colors consistent with mset[,1]
points(msetbis[,2], msetbis[,3], 
       col = mset[,1], # Use the same color from mset[,1]
       pch = 17,        # Use a different point style for L(4)
       type = "b")      # Add lines and points for consistency

# Add a legend to differentiate between datasets
legend("topright", 
       legend = c("df=1", "df=4","OLSE","beta_hat","MLE"), 
       col = c("black", "black", mset[,1]), # Use unique colors from mset[,1]
       pch = c(16, 17,15,15,15),         # Matching point styles
)                 # Add line styles to legend



par(mfrow=c(2,2))
hist(rest2004mle$V2, probability = TRUE)
lines(density(rest2004mle$V2))
hist(rest2004olse$V2, probability = TRUE)
lines(density(rest2004olse$V2))
hist(rest2004crit$V2, probability = TRUE)
lines(density(rest2004crit$V2))

plot(density(rest2004mle$V2), col="red")
lines(density(rest2004olse$V2), col="blue")
lines(density(rest2004crit$V2), col="green")

boxplot(cbind(rest2004mle$V2,rest2004olse$V2,rest2004crit$V2))

#t bis
restbis2004 <- read_csv("Desktop/Thesis GIT/restbis2004.csv")
restbis5002 <- read_csv("Desktop/Thesis GIT/restbis5002.csv")
restbis5004 <- read_csv("Desktop/Thesis GIT/restbis5004.csv")

restbis2004olse <- as.data.frame(t(as.matrix(restbis2004 %>% select(starts_with("olse")))))
restbis2004crit <- as.data.frame(t(as.matrix(restbis2004 %>% select(starts_with("crit")))))
restbis5002olse <- as.data.frame(t(as.matrix(restbis5002 %>% select(starts_with("olse")))))
restbis5002crit <- as.data.frame(t(as.matrix(restbis5002 %>% select(starts_with("crit")))))
restbis5004olse <- as.data.frame(t(as.matrix(restbis5004 %>% select(starts_with("olse")))))
restbis5004crit <- as.data.frame(t(as.matrix(restbis5004 %>% select(starts_with("crit")))))
restbis5002mle <- as.data.frame(t(as.matrix(restbis5002 %>% select(starts_with("BFGS")))))
restbis5004mle <- as.data.frame(t(as.matrix(restbis5004 %>% select(starts_with("BFGS")))))
restbis2004mle <- as.data.frame(t(as.matrix(restbis2004 %>% select(starts_with("BFGS")))))

filetbis<-list(restbis2004olse,
            restbis2004crit,
            restbis2004mle,
            restbis5004olse,
            restbis5004crit, restbis5004mle,restbis5002olse,
            restbis5002crit,restbis5002mle)

synttbis<-lapply(filet,function(x) dataresbis(x,betat))
msetbis <- as.data.frame(cbind(rep(c("OLSE","beta_hat","MLE"),3),as.numeric(rep(c(50,125,250),each=3)),as.numeric(unlist(lapply(synttbis, function(x) x$mse)))))
msetbis$V1 <- as.factor(msetbis$V1)
msetbis$V2 <- as.numeric(msetbis$V2)
msetbis$V3 <- as.numeric(msetbis$V3)

plot(msetbis[,2],msetbis[,3], xlab="Number of observation per dimension",ylab="MSE", col=mset[,1],type="p", main="T-distribution errors t(4)")
legend(legend =c("OLSE","beta_hat","MLE"), col=mset[,1],x="topright", pch=1 )

synttbis<-lapply(filetbis,function(x) dataresbis(x,betat))
par(mfrow=c(2,2))


# 1.PLOTS MSE.
par(mfrow=c(1,2))
plot(msenorm[,2],msenorm[,3], xlab="Number of observation per dimension",ylab="MSE", col=msenorm[,1],type="b", main="",pch=16 )
legend(legend =c("OLSE/MLE",expression(hat(beta))), col=mselp[,1],x="topright", pch=16 )
# Calculate global x and y limits
x_range <- range(c(mselp[,2], mselpbis[,2]))
y_range <- range(c(mselp[,3], mselpbis[,3]))

# Plot the first dataset with defined ranges
plot(mselp[,2], mselp[,3], 
     xlab = "Number of observation per dimension", 
     ylab = "MSE", 
     col = mselp[,1], 
     type = "b", 
     main = "",
     xlim = x_range, 
     ylim = y_range, 
     pch = 16) # Solid points for L(1)

# Add the second dataset, keeping colors consistent with mselp[,1]
points(mselpbis[,2], mselpbis[,3], 
       col = mselp[,1], # Use the same color from mselp[,1]
       pch = 17,        # Use a different point style for L(4)
       type = "b")      # Add lines and points for consistency

# Add a legend to differentiate between datasets
legend("topright", 
       legend = c("L(1)", "L(4)","OLSE",expression(hat(beta)),"MLE"), 
       col = c("black", "black", mselp[,1]), # Use unique colors from mselp[,1]
       pch = c(16, 17,15,15,15),         # Matching point styles
)                 # Add line styles to legend

# 2. t-tests ; is the crit that different from OLSE ?
resttestlp <- rbind(t.test(reslp2004olse$V1,reslp2004crit$V1)$conf.int[1:2],t.test(reslp2004olse$V2,reslp2004crit$V2)$conf.int[1:2],
t.test(reslp2004olse$V3,reslp2004crit$V3)$conf.int[1:2],t.test(reslp2004olse$V4,reslp2004crit$V4)$conf.int[1:2],
t.test(reslp5004olse$V1,reslp5004crit$V1)$conf.int[1:2],t.test(reslp5004olse$V2,reslp5004crit$V2)$conf.int[1:2],
t.test(reslp5004olse$V3,reslp2004crit$V3)$conf.int[1:2],t.test(reslp5004olse$V4,reslp2004crit$V4)$conf.int[1:2],
t.test(reslp5002olse$V1,reslp5002crit$V1)$conf.int[1:2],t.test(reslp5002olse$V2,reslp5002crit$V2)$conf.int[1:2])

resttestlpbis <- rbind(t.test(reslpbis2004olse$V1,reslpbis2004crit$V1)$conf.int[1:2],t.test(reslpbis2004olse$V2,reslpbis2004crit$V2)$conf.int[1:2],
t.test(reslpbis2004olse$V3,reslpbis2004crit$V3)$conf.int[1:2],t.test(reslpbis2004olse$V4,reslpbis2004crit$V4)$conf.int[1:2],
t.test(reslpbis5004olse$V1,reslpbis2004crit$V1)$conf.int[1:2],t.test(reslpbis5004olse$V2,reslpbis5004crit$V2)$conf.int[1:2],
t.test(reslpbis5004olse$V3,reslpbis2004crit$V3)$conf.int[1:2],t.test(reslpbis5004olse$V4,reslpbis2004crit$V4)$conf.int[1:2],
t.test(reslpbis5002olse$V1,reslpbis5002crit$V1)$conf.int[1:2],t.test(reslpbis5002olse$V2,reslpbis5002crit$V2)$conf.int[1:2])

resttestnorm <- rbind(t.test(res2004olse$V1,res2004crit$V1)$conf.int[1:2],t.test(res2004olse$V2,res2004crit$V2)$conf.int[1:2],
                   t.test(res2004olse$V3,res2004crit$V3)$conf.int[1:2],t.test(res2004olse$V4,res2004crit$V4)$conf.int[1:2],
                   t.test(res5004olse$V1,res2004crit$V1)$conf.int[1:2],t.test(res5004olse$V2,res5004crit$V2)$conf.int[1:2],
                   t.test(res5004olse$V3,res2004crit$V3)$conf.int[1:2],t.test(res5004olse$V4,res2004crit$V4)$conf.int[1:2],
                   t.test(res5002olse$V1,res5002crit$V1)$conf.int[1:2],t.test(res5002olse$V2,res5002crit$V2)$conf.int[1:2])

resttestlpmle <- rbind(t.test(reslp2004mle$V1,reslp2004crit$V1)$conf.int[1:2],t.test(reslp2004mle$V2,reslp2004crit$V2)$conf.int[1:2],
                      t.test(reslp2004mle$V3,reslp2004crit$V3)$conf.int[1:2],t.test(reslp2004mle$V4,reslp2004crit$V4)$conf.int[1:2],
                      t.test(reslp5004mle$V1,reslp5004crit$V1)$conf.int[1:2],t.test(reslp5004mle$V2,reslp5004crit$V2)$conf.int[1:2],
                      t.test(reslp5004mle$V3,reslp2004crit$V3)$conf.int[1:2],t.test(reslp5004mle$V4,reslp2004crit$V4)$conf.int[1:2],
                      t.test(reslp5002mle$V1,reslp5002crit$V1)$conf.int[1:2],t.test(reslp5002mle$V2,reslp5002crit$V2)$conf.int[1:2])

resttestlpbismle <- rbind(t.test(reslpbis2004mle$V1,reslpbis2004crit$V1)$conf.int[1:2],t.test(reslpbis2004mle$V2,reslpbis2004crit$V2)$conf.int[1:2],
                      t.test(reslpbis2004mle$V3,reslpbis2004crit$V3)$conf.int[1:2],t.test(reslpbis2004mle$V4,reslpbis2004crit$V4)$conf.int[1:2],
                      t.test(reslpbis5004mle$V1,reslpbis2004crit$V1)$conf.int[1:2],t.test(reslpbis5004mle$V2,reslpbis5004crit$V2)$conf.int[1:2],
                      t.test(reslpbis5004mle$V3,reslpbis2004crit$V3)$conf.int[1:2],t.test(reslpbis5004mle$V4,reslpbis2004crit$V4)$conf.int[1:2],
                      t.test(reslpbis5002mle$V1,reslpbis5002crit$V1)$conf.int[1:2],t.test(reslpbis5002mle$V2,reslpbis5002crit$V2)$conf.int[1:2])

sum(unlist(lapply(resttestlp, function(x) sign(x[1])==sign(x[2]))))
sum(unlist(lapply(resttestnorm, function(x) sign(x[1])==sign(x[2]))))
mag <- cbind.data.frame("norm"=as.numeric(unlist(apply(resttestnorm, 1, function(x) abs(x[1]-x[2])))),
                        "lp"=as.numeric(unlist(apply(resttestlp,1, function(x) abs(x[1]-x[2])))),
                      "lpbis"= as.numeric(unlist(apply(resttestlpbis,1, function(x) abs(x[1]-x[2])))),
                      "lpmle"=as.numeric(unlist(apply(resttestlpmle, 1,function(x) abs(x[1]-x[2])))),
                      "lpbismle"=as.numeric(unlist(apply(resttestlpbismle, 1,function(x) abs(x[1]-x[2])))))
# Ca donne pas grand chose.

#Comparaison des densités.

#Densités normales.
par(mfrow=c(1,3))
plot(density(res2004crit$V2),col="blue",main="",ylim=c(0,23))
lines(density(res2004olse$V2),col="green")
legend(x="topright", legend=c(expression(hat(beta)),"OLSE"),lty=1,col=c("blue","green"),cex=0.75)
plot(density(res5004crit$V2),col="blue",main="",ylim=c(0,23))
lines(density(res5004olse$V2),col="green")
legend(x="topright",legend=c(expression(hat(beta)),"OLSE"),lty=1,col=c("blue","green"),cex=0.75)
plot(density(res5002crit$V1),col="blue",main="",ylim=c(0,23))
lines(density(res5002olse$V1),col="green")
legend(x="topright",legend=c(expression(hat(beta)),"OLSE"),lty=1,col=c("blue","green"),cex=0.75)
#Densités Laplace
par(mfrow=c(1,3))
plot(density(reslp2004crit$V2),col="blue",main="",lty=1,xlim=c(-1.5,1.5),ylim=c(0,6))
lines(density(reslp2004olse$V2),col="green",lty=1)
lines(density(reslp2004mle$V2),col="red",lty=1)
lines(density(reslpbis2004crit$V2),col="blue",main="",lty=2)
lines(density(reslpbis2004olse$V2),col="green",lty=2)
lines(density(reslpbis2004mle$V2),col="red",lty=2)
legend(x="topright", legend=c("L(1)","L(4)",expression(hat(beta)),"OLSE","MLE"),lty=c(1,2,1,1,1),col=c("black","blue","green","red"),cex=0.7)

plot(density(reslp5004crit$V2),col="blue",main="",lty=1,xlim=c(-1.5,1.5),ylim=c(0,10))
lines(density(reslp5004olse$V2),col="green",lty=1)
lines(density(reslp5004mle$V2),col="red",lty=1)
lines(density(reslpbis5004crit$V2),col="blue",main="",lty=2)
lines(density(reslpbis5004olse$V2),col="green",lty=2)
lines(density(reslpbis5004mle$V2),col="red",lty=2)
legend(x="topright", legend=c("L(1)","L(4)",expression(hat(beta)),"OLSE","MLE"),lty=c(1,2,1,1,1),col=c("black","blue","green","red"),cex=0.7)

plot(density(reslp5002crit$V2),col="blue",main="",lty=1,xlim=c(-1.5,1.5),ylim=c(0,10))
lines(density(reslp5002olse$V2),col="green",lty=1)
lines(density(reslp5002mle$V2),col="red",lty=1)
lines(density(reslpbis5002crit$V2),col="blue",main="",lty=2)
lines(density(reslpbis5002olse$V2),col="green",lty=2)
lines(density(reslpbis5002mle$V2),col="red",lty=2)
legend(x="topright", legend=c("L(1)","L(4)",expression(hat(beta)),"OLSE","MLE"),lty=c(1,2,1,1,1),col=c("black","blue","green","red"),cex=0.7)

#Boxplots.
par(mfrow=c(3,3))
boxplot(res2004crit$V2,res2004olse$V2,names=c(expression(hat(beta)),"OLSE"),ylim=c(-0.35,max(max(res2004crit$V2),max(res2004olse$V2))))
abline(h=-0.16)
boxplot(res5004crit$V2,res5004olse$V2,names=c(expression(hat(beta)),"OLSE"),ylim=c(-0.35,max(max(res2004crit$V2),max(res2004olse$V2))))
abline(h=-0.16)
boxplot(res5002crit$V2,res5002olse$V2,names=c(expression(hat(beta)),"OLSE"),ylim=c(-0.35,max(max(res2004crit$V2),max(res2004olse$V2))))
abline(h=-0.16)
boxplot(reslp2004crit$V2,reslp2004olse$V2,reslp2004mle$V2,names=c(expression(hat(beta)),"OLSE","MLE"),ylim=c(min(min(reslp2004crit$V2),min(reslp2004olse$V2)),max(max(reslp2004crit$V2),max(reslp2004olse$V2))))
abline(h=-0.16)
boxplot(reslp5004crit$V2,reslp5004olse$V2,reslp5004mle$V2,names=c(expression(hat(beta)),"OLSE","MLE"),ylim=c(min(min(reslp2004crit$V2),min(reslp2004olse$V2)),max(max(reslp2004crit$V2),max(reslp2004olse$V2))))
abline(h=-0.16)
boxplot(reslp5002crit$V2,reslp5002olse$V2,reslp5002mle$V2,names=c(expression(hat(beta)),"OLSE","MLE"),ylim=c(min(min(reslp2004crit$V2),min(reslp2004olse$V2)),max(max(reslp2004crit$V2),max(reslp2004olse$V2))))
abline(h=-0.16)
boxplot(reslpbis2004crit$V2,reslpbis2004olse$V2,reslpbis2004mle$V2,names=c(expression(hat(beta)),"OLSE","MLE"),ylim=c(min(min(reslpbis2004crit$V2),min(reslpbis2004olse$V2)),max(max(reslpbis2004crit$V2),max(reslpbis2004olse$V2))))
abline(h=-0.16)
boxplot(reslpbis5004crit$V2,reslpbis5004olse$V2,reslpbis5004mle$V2,names=c(expression(hat(beta)),"OLSE","MLE"),ylim=c(min(min(reslpbis2004crit$V2),min(reslpbis2004olse$V2)),max(max(reslpbis2004crit$V2),max(reslpbis2004olse$V2))))
abline(h=-0.16)
boxplot(reslpbis5002crit$V2,reslpbis5002olse$V2,reslpbis5002mle$V2,names=c(expression(hat(beta)),"OLSE","MLE"),ylim=c(min(min(reslpbis2004crit$V2),min(reslpbis2004olse$V2)),max(max(reslpbis2004crit$V2),max(reslpbis2004olse$V2))))
abline(h=-0.16)



