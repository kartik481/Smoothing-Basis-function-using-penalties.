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
  
                                             
  
     
  S <- solve(t(R)) %*% t(D) %*% D %*% solve(R)## Defining the eigen decomposition
                                             ## vector to decompose
  ev <- eigen(S)                        ## getting Eigen values and vector
  A <- diag(ev$values)                  ## Extract the components
  U <- ev$vectors                       ## Extracting the vectors
  
  vals <- seq(logsp[1],logsp[-1],length.out=ngrid)## Creating a grid to search
                                                  ## for lambda values to try
 
   gcv <- 1e6                               ## initializing the min_gcv value  
                                             ## to very large value in order to 
                                             ## get minimum
  
  lambda <- 0                                ## making global variable in order
  edk <-0                                    ## extract it from the loop's scope
  
  I <- diag(k)                               ## creating a Indentity matrix of 
                                             ## same dimensions as k
  for (i in c(1:ngrid)){              ## iterating vals(lambda values)
   
    ## Initializing the beta with respective value of lambda stored in vals 
    
    beta <- solve(R) %*% U %*% solve( I + vals[i] * A) %*% t(U) %*% t(Q) %*% y
    
    temp_edk <- sum(diag(solve( I + vals[i] * A)))  ## Effective Degrees of Freedom
    
    fitted <- X %*% beta                         ## fitted values course
    
    sigma_2 <- (t(y-fitted) %*% (y-fitted)) / (n-temp_edk)
                                             ## Calculating the residual variance
    
    temp_gcv <- sigma_2 / (n-temp_edk)       ## computing GCV values using the 
                                             ## given formula
           
                                             
    if (temp_gcv < gcv){                     ## comparing the minimum value with
                                             ## previous results
      
      gcv <- temp_gcv                      ## getting the index of minimum value
                                           ## corresponding to gcv 
                       
                                           ## extracting the parameters from
      lambda <- vals[i]                    ## the for corresponding gcv i.e.
                                           ## lambda and effective degrees of 
                                           ## freedom i.e. edk
      edk <- temp_edk
    }
  }
  return(c(exp(lambda),gcv,edk))     ## return the optimal parameters
                                         
  
}

pspline<- function (x,y,k=20,logsp=c(-5,5),bord=3,pord=2,ngrid=100){
## function with arguments:
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
    
    optimal <- estimate(n,Q,R,D,X,y,k,logsp,ngrid) ## Calling the func estimate to
    lambda <- optimal[1]                           ## to get optimal value of λ &
    gcv <<- optimal[2]
    ed_k <<- optimal[[3]]
    print(lambda)
  }
  beta <- solve(t(X) %*% X + lambda * t(D) %*% D) %*% t(X) %*% y
  mu <- X %*% beta
  sigma_2 <- (t(y-mu) %*% (y-mu))/ (n-edk)
  return (list(lambda,beta, mu,sigma_2,edk,gcv))
}

data(mcycle)                           ## Getting the mcycle dataset from 
                                       ## mass library
x_train <- mcycle$times        ## Input values that corresponds to time 
                                       ## at an instant

y_train <- mcycle$accel         ## Labels for the train values that 
                                       ## corresponds to the acceleration at 
                                       ## that specific time

x_new <- mcycle$times[101:133]
y_new <- mcycle$accel[101:133]

attr <- pspline(x_train,y_train,k=20,logsp=c(-5,5),bord=3,pord=2,ngrid=100)

a <-list(x_train,y_train,k=20,logsp=c(-5,5),bord=3,pord=2,ngrid=100,lambda=attr[[1]]
          ,beta=attr[[2]],mu=attr[[3]],sigma_2=attr[[4]],edk=attr[[5]],gcv=attr[[6]])

class(a) <- "pspline"

estimate_errors<- function(m){
  print(length(m$x))
  e <- rep(1,length.out=length(m$x))
  y_pred <- m$x %*% m$beta + e
  
}

print.pspline <- function (m){
  cat(m$edk)
  cat("Order",m$bord,"p-spline with order",m$pord,"penalty","\n")
  cat("Effective degrees of freedom:",m$edk,"Coefficients:",length(m$beta),"\n")
  cat("Residual std dev:","r-squared:","GCV:",m$gcv,"\n")
  ##return(list(gcv,m$degree,r_2))
}

predict.pspline <- function(m,x,se=TRUE){
  Xp <- Basis_mat(x,m$k,m$bord)
  D <- diff(diag(k),differences=pord) 
  V <- solve(t(Xp) %*% Xp + m$lambda %*% t(D) %%  D) * m$sigma_2            
  std_err <- rowSums(Xp %*% (Xp %*% V))^0.5
  if (se){
    e <- rep(1,length.out=nrow(m$Xp))
    fit <- Xp %*% m$beta + e
    return(list(fit,std_err)) 
  }else{
    e <- rep(1,length.out=nrow(m$Xp))
    predictions <- Xp %*% m$beta + e
    return(predictions)
  }
  
}

plot.pspline <- function(m)
{
  plot(m$x,m$y,xlab="Data",ylab="Labels",type = 'l')
  
}
##estimate(a)
print(a)
#predict(a,x_new)

  