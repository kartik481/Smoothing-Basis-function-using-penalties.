################################################################################
############## 3: Smoothing with basis expansions and penalties ################
################################################################################

################################ Contributor ################################## 
##---------------------------  Kartik (s2407270)------------------------------##

##------------------------------ Problem Description -------------------------##

## In this practical, we have to write the R functions for smoothing x, y data. 
## The model is defined using the equation:- yi =f(xi)+εi, i=1,...,n
## where xi and yi are observed and where as f(xi) is a unknown smooth function
## for which we have to estimate the parameters. To estimate f you approximate 
## it using a basis expansion, f(x) = sum(βjbj(x)), j=1...,k where the bj(x) are 
## where as k is chosen such that it can estimate wide range of function shapes. 
## In this Problem, the B-spline basis functions are evaluated at any x value
## within the range of data. Our main aim of these is to write the linear
## model as y = Xβ + ε where Xij =bj(xi)

## But to avoid over fitting, we are introducing a smoothing penalty. 
## Then the model is penalized by using least square method:   

##         β = argmin |y−Xβ|**2 + λ * β^T *D^T *D *β (T means transpose) 

## where lambda(λ) is the smoothing parameter which we need to estimate
## by using the grid search between the specified intervals by minimizing the 
## gcv(generalized cross validation criterion). But if only one value is given 
## then no searching is done for it. Since the values of lambda are given in
## log scale but can drived taking the exponential of log(lambda). 

## If we find optimal value for lambda then we can easily find our coefficients
## Given that our smooth function estimate can vary from something very wiggly 
## to a simple straight line fit as lambda(λ) increases, it does not make sense 
## to treat its statistical degrees of freedom as k. Instead the effective 
## degrees of freedom, κ = tr{(XT X+ λDT D)−1XT X} is used.
## And from this result we can then caclculate the residual variance given by
##                        sig2(σˆ2) = |y − μ|^2 / (n − κ), 
## By calculating these parameters for our data, we can finally find GCV which  
## is called generalized cross validation criterion. [GCV= σˆ2 / (n − κ)] 
## But computing GCV directly can be very costly and can take O(k^3) operations 
## to find the value of lambda.
## To avoid this, the QR decomposition of our basis is done matrix such
## that X = QR, then after doing that we are decomposing the resultant matrix 
## (U* A* U^T= R^−T* D^T* D * R^−1) by eigen decomposition.
##                          
## By doing this trick we easily find our parameters using only O(k) operation 
## for each new value of lambda:
## The coefficients are then calculated using R^−1 *U *(I + λ*A)^−1* U^T * Q^T*y 
## edk(effective degrees of freedom) =  tr{(I + λA)^−1}
## To find sig2 and gcv the results remains same as in previous method. 

## Finally, the function can be estimated by using the appropriate penalty which
## can be seen in final plots. 
## s.t. estimated values are given by y_est = Xβ + ε



#################################### CODE  #####################################


##library(MASS)
## Used to get the mcycle data from MASS library

Basis_mat<- function (x,k,bord){
## This function is used for seeting up Basis matrix for x vector,
## where the bord is the order of B-spline to use.
  
  dk <- diff(range(x))/(k-bord)        ## knot spacing
  knots <- seq(min(x)-dk*bord,by=dk,length=k+bord+1)
  X <- splines::splineDesign(knots,x,ord=bord+1,outer.ok=TRUE)
                                       ## Spline Design function from R’s built 
                                       ## in splines package to set up matrix
 
   return(X)                           ## Returning the modified matrix
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
  
  ev <- eigen(S)                       ## getting Eigen values and vector
  A <- diag(ev$values)                 ## Extract the eigen components for A
  U <- ev$vectors                      ## Extracting the eigen vectors for U
  
  gcv <- 1e6                           ## initializing the min_gcv value  
                                       ## to very large value in order to 
                                       ## get minimum gcv to specific lambda
  
  lambda <- 0                          ## making global variables in order
  edk <-0                              ## extract it from the loop's scope
                                       ## where edk is effective degrees of 
                                       ## freedom
  
  I <- diag(k)                         ## creating a Identity matrix of 
                                       ## same dimensions as k
  
  if (length(logsp)!=1){  
  ## If the length is not equal to 1 then we have to search for optimal
  ## lambda value in the grid.
    
  vals <- exp(seq(logsp[1],logsp[-1],length.out=ngrid))
                                       ## Creating a grid to search
                                       ## for lambda values to try. Taking 
                                       ## exponential because all values are in
                                       ## log scale to convert back we use exp()
                                       ## exponential function.
    
  for (temp_lambda in vals){           ## iterating vals(i.e. lambda values)
    
    
    beta <- solve(R)%*%U%*%solve( I + temp_lambda * A)%*%t(U)%*%t(Q)%*% y
    ## Initializing the beta with respective value of lambda stored in vals 
    
    temp_edk <- sum(diag(solve( I + temp_lambda * A)))
                                       ## temporary Effective Degrees of Freedom
                                       ## for that specific lambda
    
    fitted <- X %*% beta               ## fitted values course for the specific
                                       ## value in the grid
    
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
      
      sig2 <<- as.numeric(sigma_2)           ## Optimal value of residual
                                             ## variance by minimizing GCV.
    }
  }
  
  V <- (solve(R)%*%U%*%solve( I + lambda * A)%*%t(U) %*% t(solve(R))) * sig2
                                             ## covariance matrix for the 
                                             ## coefficients(β)
  
  return(list(lambda,gcv,edk,sig2,V))        ## return the optimal parameters 
                                             ## obtained
                                       
  }else{
  
    lambda <- exp(logsp)
    ## If only a single value is provided then no searching is done,
    ## we just take the exponential of given logsp value.
      
    beta <- solve(R)%*%U%*%solve( I + lambda * A)%*%t(U)%*%t(Q)%*% y
    ## Calculating the beta with respective value of lambda.
    
    edk <- sum(diag(solve( I + lambda * A))) ## Effective Degrees of Freedom
                                             ## for that specific lambda
                                             
    
    fitted <- X %*% beta                     ## fitted values for the
                                             ## estimated coefficients.
    
    sig2 <- as.numeric((t(y-fitted) %*% (y-fitted)) / (n-edk))
                                             ## Calculating the residual  
                                             ## variance to find the GCV
    
    gcv <- sig2 / (n-edk)                    ## Computing GCV value using  
                                             ## given formula
    
    
    V <- (solve(R)%*%U%*%solve( I + lambda * A)%*%t(U) %*% t(solve(R))) * sig2
                                             ## covariance matrix for the 
                                             ## coefficients(β)
                  
    return(list(lambda,gcv,edk,sig2,V))      ## return the obtained parameters.
  }
  
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
  
  lambda <- 0                              ## initialing lambda as global var.
  
  QR <- qr(X)                              ## Doing QR decomposition of X vector
                                           ## using inbuilt qr() function.
  
  Q <- qr.Q(QR)                            ## Extracting the Q matrix from QR
  
  R <- qr.R(QR)                            ## Extracting the R matrix from QR
  
  if (length(logsp)==1){                   ## if single value (log λ scale) 
                                           ## is given.
    
    optimal <<- estimate(n,Q,R,D,X,y,k,logsp,ngrid)
                                           ## getting the required components
                                           ## in a list named optimal from the
                                           ## estimate function.
    
    lambda <<- optimal[[1]]                ## Getting optimal lambda from list
    
    gcv <<- optimal[[2]]                   ## Getting GCV value
    
    edk <<- optimal[[3]]                   ## Extracting the effective degrees 
                                           ## of freedom
    
    sig2 <<- optimal[[4]]                  ## Extracting the residual variance
    
    V <<- optimal[[5]]                     ## covariance matrix for the 
                                           ## coefficients, β,
    
  }else{
    
    optimal <- estimate(n,Q,R,D,X,y,k,logsp,ngrid) 
                                           ## Calling the function estimate to
    lambda <- optimal[[1]]                 ## to get optimal value of λ 
                                                      
    gcv <<- optimal[[2]]                   ## Extracting the optimal value of 
    edk <<- optimal[[3]]                   ## GCV, edk (effective degrees of 
                                           ## freedom)
    
    sig2 <<- optimal[[4]]                  ## Residual variance 
    
    V <<- optimal[[5]]                     ## Covariance  matrix for coefficient
                                           ## beta
  
  }
  beta <- solve(t(X) %*% X + lambda * t(D) %*% D) %*% t(X) %*% y
                                           ## Calculating the coefficients(beta)
                                           ## for optimal(best) value of lambda
                                          
                                    
  fitted <- X %*% beta                     ## Calculating the course of fitted
                                           ## values
  
  res <- y - fitted                        ## Calculating the residuals values
  
  m <- list(x=x, y=y, k=k, logsp=logsp, bord=bord, pord=pord, ngrid=ngrid,
       lambda=lambda, coef=beta, fitted=fitted, sig2=sig2, edk=edk, gcv=gcv, V=V
       , residuals=res)
                                           ## Creating a list of objects for
                                           ## class pspline 
  
  class(m)<- "pspline"                     ## Created pspline class with the 
                                           ## attributes given in the list.
  
  return(m)                                ## returning the objects of pspline
                                           ## class to caller
}

print.pspline <- function (m){
## This a method function for pspline class. The arguments are the objects of 
## pspline class. 
  
  cat("Order", m$bord, "p-spline with order", m$pord, "penalty","\n")
                                      ## Printing the the spline order and 
                                      ## and penalty
  
  cat("Effective degrees of freedom:",m$edk,"Coefficients:",length(m$coef),"\n")
                                      ## Printing the Effective degrees of 
                                      ## freedom and coefficients
  
  
  r2 <- 1-(length(m$y)-1)* m$sig2 / sum(t(m$y-mean(m$y)) %*% (m$y-mean(m$y)))
                                      ## Calculating the r^2 error using the 
                                      ## given formula
  
  cat("Residual std dev:",sqrt(m$sig2),"r-squared:",r2,"GCV:",m$gcv,"\n")
                                      ## Printing the Residual standard 
                                      ## deviation, r^2 and gcv
  
  invisible(list(m$gcv, m$edk, r2))   ## Silently returning the GCV, edk, r^2
                                      ## using invisible 
                                    
}

predict.pspline <- function(m,x,se=TRUE){
## This function should make predictions from the smooth fit, for new x values 
## within the range of the original data. The arguments are object of pspline
## class(m), new x values(x), se is used for returning specific lists.
  
  Xp <- Basis_mat(x, m$k, m$bord)     ## Setting the Basis matrix using the 
                                      ## original settings.
  
  fit <- Xp %*% m$coef                ## New Fitted values using the 
                                      ## coefficients obtained from the smooth
                                      ## fit
  
  std_err <- rowSums(Xp * (Xp %*% m$V))^0.5
                                      ## Calculating the standard error by 
                                      ## taking the sqrt of covariance matrix.
                                      
  
  if (se){     
    ## If se is TRUE then return the list name fit with corresponding 
    ## se(standard error).
    
    l <- list(fit=fit,se=std_err)     ## creating a named list 
    
    return(l)                         ## returning the list created
    
  }else{
    ## If se is FALSE then returning only new fitted values.

    return(fit)                       ## returning only fitted values.
  }
  
}

plot.pspline <- function(m){
## This function plots the three different types of plots by taking pspline
## class object(m) as a argument. 
  layout(matrix(c(1,1,2,3),2,2,byrow=TRUE))   ## splitting the plot window in 
                                              ## areas with custom sizes using 
                                              ## layout function, where 
                                              ## matrix in function specifying 
                                              ## the location of figures.
  
  plot(m$x, m$y, main=" Smooth function fitted to the Data", 
  xlab="x", ylab="y")          
                                              ## Plotting the 1st figure of 
                                              ## x and y (original data) 
   
  
  
  se <- predict(m,m$x)[[2]]                   ## Getting standard errors from
                                              ## the previous predict.pspline
                                              ## method which are stored at 
                                              ## index 2 of the list.
  
  ul <- m$fitted + 1.96*se                    ## upper 95% credible interval 
   
  ll <- m$fitted - 1.96*se                    ## lower 95% credible interval
  
  lines(m$x,m$fitted,col="red")               ## Plotting estimated smooth
                                              ## function overlaid as line
  
  lines(m$x, ul,type="l", lty=2, col="blue")  ## Plotting the upper 95% upper
                                              ## credible intervals 
                                             
  lines(m$x, ll, type="l", lty=2, col="blue") ## Plotting the upper 95% lowe
                                              ## credible intervals
  
  legend("topleft", legend=c("Fitted smooth function", "95% credible intervals")
  , col=c("red","blue"),lty=1:2 ,cex=0.437, bty="n") 
                                              ## Adding legend to the plot 
                                              ## describing the plot
  
  
  plot(m$fitted, m$residuals, main="Residuals against fitted values", 
  xlab="Fitted values", ylab="Residuals", col="red")
                                              ## plot the model residuals 
                                              ## against fitted values
  
  qqnorm(m$residuals, pch = 1, frame= FALSE)  ## qqplot of the residuals.
  qqline(m$residuals, col="red", lwd=1)       ## Adding line to qqplot
  
  invisible(list(ll,ul,m$x))                  ## returning the list silently
                                              ## using invisible function.
}

data(mcycle)                           ## Getting the mcycle dataset from 
                                       ## mass library

x <- mcycle$times                      ## Input values that corresponds to time 
                                       ## at an instant

y <- mcycle$accel                      ## Labels for the train values that 
                                       ## corresponds to the acceleration at 
                                       ## that specific time



model <- pspline(x, y)
## Storing the results returned from pspline function in a object named model.

print(model)                             ## print is a method of pspline class
                                         ## which is printing details about our
                                         ## fitted model.

set.seed(0)     
## Setting the seed so that data don't randomize every time

x_new <- runif(100,min(x),max(x))        ## New x values within the range of the 
                                         ## original data generated using 
                                         ## uniform distribution

predict(model,x_new)                     ## Making predictions for newly 
                                         ## generated x values using the predict
                                         ## method of pspline

plot(model)                              ## Calling the method plot to plot
                                         ## three required plots
