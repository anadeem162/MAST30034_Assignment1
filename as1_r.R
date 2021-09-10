TC <- read.csv(file = 'TC.csv', header = FALSE)


TC = as.matrix(TC)


X <- read.csv(file = 'X_df.csv', header = FALSE)


X = as.matrix(X)


### Question 2.3

# number of sources
nsrcs = 6 

# horizontal slice size
x1 = 21

# vertical slice size
x2 = 21

# number of time samples
N = 240



rho_list = c()
avg_MSE_list = c()

# counts number of iterations in following loop up until 21
count = 0

for(rho in seq(0, 1, 0.05)){

  step <- 1/(norm(TC %*% t(TC)) * 1.1)
  thr <- rho*N*step
  Ao <- matrix(0, nsrcs, 1)
  A <- matrix(0, nsrcs, 1)
  Alr <- matrix(0, nsrcs, x1*x2)
  
  for (k in 1:(x1*x2)) {
    A <- Ao+step*(t(TC)%*%(X[,k]-(TC%*%Ao)))
    A <- (1/(1+thr)) * (sign(A)*pmax(replicate(nsrcs, 0), abs(A)-thr))
    
    for (i in 1:10) {
      Ao <- A
      A <- Ao+step * (t(TC)%*%(X[,k]-(TC%*%Ao)))
      A <- (1/(1+thr)) * (sign(A)*pmax(replicate(nsrcs, 0), abs(A)-thr))
    }
    Alr[,k] <- A 

  }

  
  Dlr = X%*%t(Alr)
  
  V = 441
  
  MSE = sum( (X -(Dlr %*% Alr) )^ 2 ) / (N*V)

  
  
  rho_list <- c(rho_list, rho)
  avg_MSE_list <- c(avg_MSE_list, MSE/10)
  
  count = count + 1
  
}


print(rho_list)
print(avg_MSE_list)

png(file="2.3_Avg_MSE_curve.png",
    width=600, height=350)

plot(rho_list, avg_MSE_list,
     main="Average MSE Curve",
     xlab="rho",
     ylab="Average of MSE",
     type = "l",
     lwd=2,
     col="blue")

dev.off()

###install.packages("gplots")
library(gplots)

png(file="2.3_rho_Avg_MSE_Table.png",
    width=600, height=350)

textplot(cbind(rho_list, avg_MSE_list))

dev.off()


##############################################################################################################################


### Question 2.4


rho = 0.25
  
step <- 1/(norm(TC %*% t(TC)) * 1.1)
thr <- rho*N*step
Ao <- matrix(0, nsrcs, 1)
A <- matrix(0, nsrcs, 1)
Alr <- matrix(0, nsrcs, x1*x2)

for (k in 1:(x1*x2)) {
  A <- Ao+step*(t(TC)%*%(X[,k]-(TC%*%Ao)))
  A <- (1/(1+thr)) * (sign(A)*pmax(replicate(nsrcs, 0), abs(A)-thr))
  
  for (i in 1:10) {
    Ao <- A
    A <- Ao+step * (t(TC)%*%(X[,k]-(TC%*%Ao)))
    A <- (1/(1+thr)) * (sign(A)*pmax(replicate(nsrcs, 0), abs(A)-thr))
  }
  Alr[,k] <- abs(A)   # here applied absolute
}

Dlr = X%*%t(Alr)


write.csv(Dlr, file = "Dlr.csv")
write.csv(Alr, file = "Alr.csv")




#########################################################################################

### Question 2.5

Z <- read.csv(file = 'Z.csv', header = FALSE)
#View(Z)

Z = as.matrix(Z)


X <- read.csv(file = 'X_df.csv', header = FALSE)

X = as.matrix(X)

# number of sources
nsrcs = 6 

# horizontal slice size
x1 = 21

# vertical slice size
x2 = 21

# number of time samples
N = 240

V = 441

rho = 0.001

step <- 1/(norm(Z %*% t(Z)) * 1.1)
thr <- rho*N*step
Ao <- matrix(0, nsrcs, 1)
A <- matrix(0, nsrcs, 1)
APCR <- matrix(0, nsrcs, x1*x2)

for (k in 1:(x1*x2)) {
  A <- Ao+step*(t(Z)%*%(X[,k]-(Z%*%Ao)))
  A <- (1/(1+thr)) * (sign(A)*pmax(replicate(nsrcs, 0), abs(A)-thr))
  
  for (i in 1:10) {
    Ao <- A
    A <- Ao+step * (t(Z)%*%(X[,k]-(Z%*%Ao)))
    A <- (1/(1+thr)) * (sign(A)*pmax(replicate(nsrcs, 0), abs(A)-thr))
  }
  APCR[,k] <- abs(A)   # here applied absolute
}

DPCR = X%*%t(APCR)

write.csv(DPCR, file = "DPCR.csv")
write.csv(APCR, file = "APCR.csv")

