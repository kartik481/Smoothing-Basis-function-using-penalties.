################################################################################
############## 3: Smoothing with basis expansions and penalties ################
################################################################################

################################ Contributors ################################## 

##---------------------------  Kartik (s2407270)------------------------------##

#################################### CODE  #####################################

library(methods)

Basis_mat<- function (x,k,bord){

  dk <- diff(range(x))/(k-bord)                ## knot spacing
  knots <- seq(min(x)-dk*bord,by=dk,length=k+bord+1)
  X <- splines::splineDesign(knots,x,ord=bord+1,outer.ok=TRUE)
  return(X)
}
#Basis_mat(c(1:3),20,3)


eigen_decom<- function(A){
## functions (A) arguments is the matrix multiplication of given equation and 
## using that we do the eigen decomposition of the resultant vector.
  
  ev <- eigen(A)                        ## getting Eigen values and vector
  A <- diag(ev$values)                  ## Extract the components
  U <- ev$vectors                       ## Extracting the vectors
  
  decom <- cbind(U,A)                   ## Combining vectors in order to return          

  return(decom)                         ## return the vector that contains U and
                                        ## A to find optimal gcv.
}


min_gcv <- function(Q,R,D,X,y,logsp,ngrid){

## This function estimates the optimal value of lambda by minimizing the gcv 
## value. The function takes argument Q, R (Founded using QR decomposition), X
## (Basis Matrix), y (labels), logsp (the ends of interval to search for lambda)
## ngrid (the of lambda values to try)
  
  n <-nrow(X)                                ## length of number of inputs 
                                             ## can be extracted from X vector
                                             ## rows elements as X is N*N vector
  
     
  S <- t(solve(R))%*% t(D) %*% D %*% solve(R)## Defining the eigen decomposition
                                             ## vector to decompose
  decom <- eigen_decom(S)                    ## Function eigen_decom will 
                                             ## decompose the into U, A vectors
  
  U <- decom[1:20,1:20]                      ## Extracting the vector U from 
                                             ## decom matrix
  
  A <- decom[1:20,21:40]                     ## Extracting the vector U from 
                                             ## decom matrix
  
  vals <- seq(logsp[1],logsp[-1],length.out=ngrid)## Creating a grid to search
                                                  ## for lambda values to try
 
   min_gcv <- 1e9                            ## initializing the min_gcv value  
                                             ## to very large value in order to 
                                             ## get minimum
  
  lambda <- 0                                ## making global variable in order
                                             ## extract it from the loop's scope
  
  sigma_2 <- 0                               ## making global variable in order
                                             ## extract it from the loop's scope
  
  I <- diag(nrow(A))                         ## creating a Indentity matrix of 
                                             ## same dimensions as A vector

  for (i in c(1:ngrid)){              ## iterating vals(lambda values)
   
    ## Initializing the beta with respective value of lambda stored in vals 
    beta <- solve(R) %*% U %*% solve( I + vals[i] * A) %*% t(U) %*% t(Q) %*% y
    k <- sum(diag(solve( I + vals[i] * A)))  ## Effective Degrees of Freedom
    mu <- X %*% beta                         ## fitted values course
    sigma_2 <- ((y-mu)**2) / (n-k)           ## Calculating the residual variance
    gcv <- sigma_2 / (n-k)                   ## computing GCV values using the 
                                             ## given formula
   
     optimal<-cbind(gcv,sigma_2)             ## Creating a matrix to get 
                                             ## the minimum value of gcv and 
                                             ## other paramters lambda, sigma**2
    
     if (min(optimal[1]) < min_gcv){         ## comparing the minimum value with
                                             ## previous results
      
      index <- which.min(optimal[1])      ## getting the index of minimum value
                                          ## corresponding to gcv 
      
      min_gcv <- optimal[index,1]         ## extracting the parameters from
      lambda <- vals[i]                   ## the vector corresponding to the
      min_sigma_2=optimal[index,2]        ## index
    }
  
  }
  
  return(c(lambda,sigma_2,min_gcv))      ## return the optimal parameters
  
}

pspline<- function (x,y,k=20,logsp=c(-5,5),bord=3,pord=2,ngrid=100){

  X <- Basis_mat(x,k,bord)
  e <- rep(1,length=length(x))
  D <- diff(diag(k),differences=pord)
  lambda <- 0
  
  if (length(logsp)==1){
    lambda <- logsp[1]
  }else{
    QR <- qr(X)
    Q <- qr.Q(QR)
    R <- qr.R(QR)
    estimate <<- min_gcv(Q,R,D,X,y,logsp,ngrid)
    lambda <- estimate[1]
    sigma_2 <- estimate[2]
    gcv <- estimate[3]
  }

  beta <- solve(t(X) %*% X + lambda * t(D) %*% D)* t(X) *y
  mu <- X %*% beta
  #return (c(beta,mu,sigma_2))
  return(X%*% beta +e )
}
x=c(1:20)
y=sin(x)
f=pspline(x,y)
plot(x,y,type = "l")
plot(x,diag(f),pch=19,col=4,xlab = "Input",ylab = "Predictions",type="l")
