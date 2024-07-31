#install.packages("CVXR")
#install.packages("Rcsdp")
library(CVXR)
library(Rcsdp)

generate_list <- function(r, m, length_out){
  my_list <- numeric(length_out)
  if (r <= length_out) {
    my_list[r] <- 1
  }
  if ((m + r) <= length_out) {
    my_list[m + r] <- -1
  }
  return(my_list)
}
generate_matrix<-function(i,j,N){
  my_matrix <- matrix(0,N,N)
  if(i==j){
    my_matrix[i,j]= 1
  }else{
    my_matrix[i,i]= 1
    my_matrix[j,j]= 1
    my_matrix[i,j]= -1
    my_matrix[j,i]= -1
  }
  return(my_matrix)
}
index_function<-function(n,k,plt=FALSE){
  index_number = rep(0,n)
  for(i in 1:n){
    if(i==1){
      index_sub = sample(2:n,k)
      result = list(cbind(rep(1,k),index_sub))
      I = rep(0,n)
      index_number[index_sub] = I[index_sub]+1
    }else{
      if(index_number[i]<k){
        L = which(index_number<k)
        index_s = L[L>i]
        num = k-index_number[i]
        if(num > length(index_s)){
          index_sub1 = c()
          for(j in 1:length(result)){
            if(length(result[[j]])!= 0){
              if(!(i %in% result[[j]][,2])){
                index_sub1 = append(index_sub1,j)
              }
            }
          }
          index_sub1 = append(index_sub1,index_s)
          index_sub = sample(index_sub1,num)
          result = c(result,list(cbind(rep(i,num),index_sub)))
          index_number[index_sub] = index_number[index_sub] + 1 
        }else{
          index_sub = sample(index_s,num)
          result = c(result,list(cbind(rep(i,num),index_sub)))
          index_number[index_sub] = index_number[index_sub] + 1 
        }
      }else{
        next
      }
    }
  }
  
  Index = matrix(0,nrow=1,ncol=2)
  for(l in 1:length(result)){
    Index = rbind(Index,result[[l]])
  }
  Index =Index[-1,]
  Index = Index[(Index[,1] != Index[,2]) & (Index[,1] < Index[,2]),]
  if(plt==TRUE){
    print(table(c(Index[,1],Index[,2])))
    print(nrow(Index))
  }
  return(Index)
}

b_new_function<-function(d,K,N){
  b_list <- c(1,0,numeric(N))
  for(i in 3:(N+2)){
    b_list[i] = d[i-2]- K[i-2,i-2]
  }
  return(b_list)
}

C_new_function<-function(K,N){
  p = qr(K)$rank
  newC1 <- list(-rep(1,2*N),-matrix(c(0,0,0,0),2,2,byrow=TRUE),-rep(0,p+1))
  return(newC1)
}

A_new_function<-function(d,K,N){
  p = qr(K)$rank
  for(i in 1:(N+2)){
    if(i==1){
      A_sub = list(rep(0,2*N),matrix(c(1,0,0,0),2,2,byrow=TRUE),rep(0,p+1))
      A = list(A_sub)
    }else if(i==2){
      A_sub = list(rep(0,2*N),matrix(c(0,0.5,0.5,0),2,2,byrow=TRUE),-c(1,rep(0,p)))
      A=c(A,list(A_sub))
    }else{
      Sigma = cbind(rep(0,N),2*(svd(K)$u[,1:p] %*% diag(sqrt(svd(K)$d))[1:p,1:p]))
      A_sub = list(-generate_list(r=i, m=N, length_out=2*N),matrix(c(0,0,0,1),2,2,byrow=TRUE),-Sigma[i-2,])
      A=c(A,list(A_sub))
    }
  }
  return(A)
}

gen_kernel_matrix<-function(D,sample,n_size,lambda){
  n= nrow(D)
  if(sample==FALSE){
    m = n*(n-1)/2
    k1=1
    index=c()
    index = matrix(0,nrow=m,ncol=2)
    for(i in 1:(n-1)){ for(j in (i+1):n){index[k1,]=c(i,j) 
    k1=k1+1}}
    b=c()
    for(i in 1:m){
      if(i==1){
        A_sub = list(generate_matrix(index[i,1],index[i,2],n),-generate_list(r=i,m=m,length_out=2*m))
        A=list(A_sub)
        b=append(b,D[index[i,1],index[i,2]])
      }else{
        A_sub = list(generate_matrix(index[i,1],index[i,2],n),-generate_list(r=i,m=m,length_out=2*m))
        A=c(A,list(A_sub))
        b =append(b,D[index[i,1],index[i,2]])
      }
    }
    C = list(-lambda*diag(n),-rep(1,2*m))
    K <- list(type=c("s","l"),size=c(n,2*m))
    result = csdp(C,A,b,K,control = csdp.control(printlevel = 0))
    return(result$X[[1]])
  }else{
    index = index_function(n=n,k=n_size)
    b=c()
    for(i in 1:nrow(index)){
      if(i==1){
        A_sub = list(generate_matrix(index[i,1],index[i,2],n),-generate_list(r=i,m=nrow(index),length_out=2*nrow(index)))
        A=list(A_sub)
        b=append(b,D[index[i,1],index[i,2]])
      }else{
        A_sub = list(generate_matrix(index[i,1],index[i,2],n),-generate_list(r=i,m=nrow(index),length_out=2*nrow(index)))
        A=c(A,list(A_sub))
        b =append(b,D[index[i,1],index[i,2]])
      }
    }
    C = list(-lambda*diag(n),-rep(1,2*nrow(index)))
    K <- list(type=c("s","l"),size=c(n,2*nrow(index)))
    result = csdp(C,A,b,K,control = csdp.control(printlevel = 0))
    return(result$X[[1]])
  }
}
Newkernel_solve <- function(d,Kernel_matrix,n){
  p = qr(Kernel_matrix)$rank
  u <- Variable(2*n)
  Z <- Variable(2,2, PSD= TRUE)
  x <- Variable(p+1)
  Sigma = cbind(rep(0,n),2*(svd(Kernel_matrix)$u[,1:p] %*% diag(sqrt(svd(Kernel_matrix)$d))[1:p,1:p]))
  
  for(i in 1:(n+2)){
    if(i==1){
      constraints <-list(sum(diag(matrix(c(1,0,0,0),2,2,byrow=TRUE) %*% Z))==1)
    }else if(i==2){
      constraints <- c(constraints,
                       list(sum(diag(matrix(c(0,0.5,0.5,0),2,2,byrow=TRUE) %*% Z))-matrix(c(1,rep(0,p)),nrow=1,ncol=p+1,byrow=TRUE) %*% x == 0))
    }else{
      constraints <- c(constraints,list(d[i-2]-Kernel_matrix[i-2,i-2]-sum(diag(matrix(c(0,0,0,1),2,2,byrow=TRUE) %*% Z))+
                                          matrix(Sigma[i-2,],nrow=1,ncol=p+1,byrow=TRUE) %*% x +
                                          matrix(generate_list(r=i-2, m=n, length_out=2*n),nrow=1,ncol=2*n,byrow=TRUE) %*% u == 0))
    }
  }
  constraints <-c(constraints,list(norm2(x[2:(p+1)])<=x[1]))
  constraints <-c(constraints,list(u>=0))
  #constraints <-c(constraints,list(eigen(Z)$values>=0))
  C1 <- matrix(rep(1,2*n),nrow=1,ncol=2*n,byrow=TRUE)
  objective <- Minimize(C1 %*% u)
  problem <- Problem(objective, constraints)
  result  <- psolve(problem)
  c = result$getValue(Z)[2,2]
  b = (svd(Kernel_matrix)$u[,1:p] %*% diag(sqrt(svd(Kernel_matrix)$d))[1:p,1:p]) %*% t(t(result$getValue(x)[-1]))
  K_new = rbind(cbind(Kernel_matrix,b),c(t(b),c))
  re = list(result=result,K_new=K_new)
  return(K_new)
}

Newkernel_solve1<- function(d,K,N){
  p = qr(K)$rank
  newA1 = A_new_function(d=d,K=K,N=N)
  newb1 = b_new_function(d=d,K=K,N=N)
  newC1 = C_new_function(K=K,N=N)
  newK1 <- list(type=c("l","s","l"),size=c(2*N,2,p+1))
  result = csdp(newC1,newA1,newb1,newK1, control = csdp.control(printlevel = 0))
  c<-result$X[[2]][2,2]
  b = (svd(K)$u[,1:p] %*% diag(sqrt(svd(K)$d))[1:p,1:p]) %*% t(t(result$X[[3]][-1]))
  K_new = rbind(cbind(K,b),c(t(b),c))
  return(K_new)
}

#set.seed(1234)

# E.1
#x = rnorm(101,0,1)
#D_all = as.matrix(dist(x))
#D = D_all[1:100,1:100]
#d = D_all[101,1:100]
#n = nrow(D)

## gen 100*100 kernel matrix  D:distance_matrix; sample: yes or no sample; n_size: show times for every i; 
#kernel_matrix = gen_kernel_matrix(D,sample=TRUE,n_size=30,lambda = 1)
#eigen(kernel_matrix)$values

# E.2

## y = sin(3/4pi*x) + x+ e
#sig = 0.01
#x = rnorm(101,0,1)
#y = sin(3*pi*x/4) + x + rnorm(101,0,sig)
#D_all = as.matrix(dist(x))
#D = D_all[1:100,1:100]
#d = D_all[101,1:100]
#n = nrow(D)

## gen kernel matrix
#kernel_matrix = gen_kernel_matrix(D,sample=TRUE,n_size=30,lambda = 1)
#new_kernel_matrix = Newkernel_solve1(d=d, K=kernel_matrix, N=n)
#new_kernel_matrix = Newkernel_solve(d=d,Kernel_matrix=kernel_matrix,n=n)

## gen y_hat by guassian process model
#y_hat = t(as.matrix(new_kernel_matrix[101,1:100],nrow=1,ncol=100)) %*% solve(kernel_matrix + sig^2 * diag(n)) %*% as.matrix(y[1:100],nrow=1,ncol=100)
#error = (y_hat-y[101])^2
#error1 = ((y_hat-y[101])/y[101])^2
#print(c(error, error1))


