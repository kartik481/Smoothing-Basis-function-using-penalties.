################################################################################
############## 3: Smoothing with basis expansions and penalties ################
################################################################################

################################ Contributors ################################## 

##---------------------------  Kartik (s2407270)------------------------------##

#################################### CODE  #####################################

library(MASS)
library(methods)

Basis_mat<- function (x,k,bord){
  
  dk <- diff(range(x))/(k-bord)                ## knot spacing
  knots <- seq(min(x)-dk*bord,by=dk,length=k+bord+1)
  X <- splines::splineDesign(knots,x,ord=bord+1,outer.ok=TRUE)
  return(X)
}


estimate <- function(n,Q,R,D,X,y,k,logsp,ngrid){
## This function estimates the optimal value of lambda by minimizing the gcv 
## value. The function takes argument Q, R (Founded using QR decomposition), X
## (Basis Matrix), y (labels), logsp (the ends of interval to search for lambda)
## ngrid (the of lambda values to try)
  
  
  
  
  S <- t(solve(R)) %*% t(D) %*% D %*% solve(R)
                                             ## Defining the eigen decomposition
                                             ## vector to decompose to speed up  
                                             ## search for optimal lambda
  
  ev <- eigen(S)                        ## getting Eigen values and vector
  A <- diag(ev$values)                  ## Extract the eigen components for A
  U <- ev$vectors                       ## Extracting the eigen vectors for U
  
  vals <- exp(seq(logsp[1],logsp[-1],length.out=ngrid))
                                       ## Creating a grid to search
                                       ## for lambda values to try. Taking 
                                       ## exponential because all values are in
                                       ## log scale to convert back we use exp()
                                       ## exponential function.
  
  
  gcv <- 1e6                           ## initializing the min_gcv value  
                                       ## to very large value in order to 
                                       ## get minimum gcv to specific lambda
  
  lambda <- 0                          ## making global variables in order
  edk <-0                              ## extract it from the loop's scope
                                       ## where edk is effective degrees of 
                                       ## freedom
  
  I <- diag(k)                         ## creating a Identity matrix of 
                                       ## same dimensions as k
  
  for (temp_lambda in vals){           ## iterating vals(i.e. lambda values)
    
    ## Initializing the beta with respective value of lambda stored in vals 
    beta <- solve(R)%*%U%*%solve( I + temp_lambda * A)%*%t(U)%*%t(Q)%*% y
    
    temp_edk <- sum(diag(solve( I + temp_lambda * A)))
                                       ## temporary Effective Degrees of Freedom
                                       ## for that specific lambda
    
    fitted <- X %*% beta               ## fitted values course
    
    sigma_2 <- (t(y-fitted) %*% (y-fitted)) / (n-temp_edk)
                                       ## Calculating the residual variance to 
                                       ## find the temporary gcv in order to 
                                       ## minimize it
    
    temp_gcv <- sigma_2 / (n-temp_edk) ## computing temporary GCV values using  
                                       ## given formula
    
    
    if (temp_gcv < gcv){                     ## comparing the minimum value with
                                             ## previous results
      
      gcv <- temp_gcv                        ## getting the minimum value of the
                                             ## corresponding gcv 
      
                                             ## extracting the parameters from
      lambda <- temp_lambda                  ## the for corresponding gcv i.e.
                                             ## lambda  
      
      edk <- temp_edk                        ## Effective degrees of freedom edk 
      opt_sigma_2 <<-sigma_2                 ## and optimal value of residual
                                             ## variance by minimizing GCV.
    }
  }
  return(c(lambda,gcv,edk,opt_sigma_2))   ## return the optimal parameters 
                                          ## obtained
  
  
}

pspline<- function (x,y,k=20,logsp=c(-5,5),bord=3,pord=2,ngrid=100){
## function arguments are as follows:
## x and y – the vectors of x, y data to smooth. k – the number of basis 
## functions to use logsp – the ends of the interval over which to search for 
## the smoothing parameter (log λ scale). 
## bord – the B-spline order to use: 3 corresponds to cubic.
## pord – the order of difference to use in the penalty. 2 is for the penalty 
## ngrid – the number of smoothing parameter values to try. 
  
  X <- Basis_mat(x,k,bord)                 ## Setting up the basis (X matrix)
  
  n <- length(x)                           ## No. of input values or datapoints
  
  D <- diff(diag(k),differences=pord)      ## D is a k−2 × k matrix of zeroes
                                           ## used to estimate beta
  
  lambda <- 0                              ## initialing as global var.
  
  if (length(logsp)==1){
    lambda <- logsp[1]
  }else{
    QR <- qr(X)                            ## Doing QR decomposition of X vector 
    Q <- qr.Q(QR)                          ## Extarcting the Q matrix
    R <- qr.R(QR)                          ## Extracting the R matrix
    
    optimal <- estimate(n,Q,R,D,X,y,k,logsp,ngrid) 
                                           ## Calling the func estimate to
    lambda <- optimal[1]                   ## to get optimal value of λ 
                                                      
    gcv <<- optimal[2]                     ## Extracting the optimal value of 
    edk <<- optimal[3]                     ## GCV, edk(effective degrees of 
                                           ## freedom)
    sigma_2 <<- optimal[4]                 ## Residual variance                            
    
  }
  beta <- solve(t(X) %*% X + lambda * t(D) %*% D) %*% t(X) %*% y
  fit <- X %*% beta
  return (list(lambda,beta, fit,sigma_2,edk,gcv))
}

data(mcycle)                           ## Getting the mcycle dataset from 
                                       ## mass library

x_train <- mcycle$times                ## Input values that corresponds to time 
                                       ## at an instant

y_train <- mcycle$accel                ## Labels for the train values that 
                                       ## corresponds to the acceleration at 
                                       ## that specific time



attr <- pspline(x_train,y_train,k=20,logsp=c(-5,5),bord=3,pord=2,ngrid=100)

a <-list(x=x_train,y=y_train,k=20,logsp=c(-5,5),bord=3,pord=2,ngrid=100,lambda=attr[[1]]
         ,beta=attr[[2]],fit=attr[[3]],sigma_2=attr[[4]],edk=attr[[5]],gcv=attr[[6]])

class(a) <- "pspline"


print.pspline <- function (m){
  cat("Order",m$bord,"p-spline with order",m$pord,"penalty","\n")
  cat("Effective degrees of freedom:",m$edk,"Coefficients:",length(m$beta),"\n")
  
  r_2 <- 1-(length(m$y)-1)* sigma_2 / sum(t(m$y-mean(m$y)) %*% (m$y-mean(m$y)))
  
  cat("Residual std dev:",sqrt(m$sigma_2),"r-squared:",r_2,"GCV:",m$gcv,"\n")
  
  return(list(gcv,m$edk,r_2))
}

predict.pspline <- function(m,x,se=TRUE){
  Xp <- Basis_mat(x,m$k,m$bord)
  
  D <- diff(diag(m$k),differences=m$pord) 
  
  V <- solve(t(Xp) %*% Xp + m$lambda * t(D) %*%  D) * sigma_2            
  
  std_err <- sqrt(rowSums(Xp * (Xp %*% V)))
  if (se){
    fit <- Xp %*% m$beta
    se <- std_err
    return(list(fit,se)) 
  }else{
    
    predictions <- Xp %*% m$beta 
    return(predictions)
  }
  
}

plot.pspline <- function(m)
{
  layout(matrix(c(1,1,2,3),2,2,byrow=TRUE))   ## splitting the plot window in 
                                              ## areas with custom sizes using 
                                              ## layout function, where 
                                              ## matrix in function specifying 
                                              ## the location of figures.
  
  plot(m$x,m$y,xlab="Data",ylab="Labels",col=7)
  ## Plotting the 1st figure 
  lines(m$x,m$fit,col="red")
  
  se <- predict(m,m$x)[[2]]
  print(se)
  upper<- m$fit + 1.96*se                                # upper 95% conf. band
  lower <- m$fit - 1.96*se                               # lower 95% conf. band
  
  lines(m$x, upper, type="plp", pch="-",col=4)
  lines(m$x, lower, type="plp", pch="-",col=4)
  
  
  residuals <- m$y-m$fit
  plot(m$fit,residuals,xlab="Fitted values",ylab="Residuals",col="red")
  
  qqnorm(residuals, pch = 1, frame= FALSE)
  qqline(residuals, col="steelblue", lwd=2)
  
  
}

imp <-print(a)

pred <-predict(a,x_train)
plot(a)
