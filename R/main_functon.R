library(lme4)
library(numDeriv)
library(optimParallel)
library(nloptr)
library(dplyr)
library(glmmTMB)
library(parallel)

# grad(dummy,x,y=1)
## X is the HiC matrix (only lower triangular)
## P is the corresponding 3D coordinates
## Par is parameters
## Calculate Distance first, then turn back to 3D coordinate

objective_poisson=function(P.V,Par,X,random,method.type)
{
  X=as.matrix(X)
  X[upper.tri(X,diag=TRUE)]=0
  N=nrow(X)
  P=matrix(P.V,ncol=3)
  D.V=as.vector(sqrt(dist(P)))
  D.mat=matrix(1,N,N)
  D.mat[lower.tri(D.mat)]=D.V
  D.mat=t(D.mat)
  D.mat[lower.tri(D.mat)]=D.V
  D.V.nozero=D.V[D.V>0]
  D.V.nozero.min=min(D.V.nozero)
  D.mat[D.mat==0]=D.V.nozero.min/10000
  lambda=exp(Par[1]+random+Par[2]*log(D.mat))
  result=matrix(0,N,N)
  if(method.type=='tRex'|method.type=='tPAM')
  {
  result[lambda<=30]=X[lambda<=30]*log(lambda[lambda<=30])-log(exp(lambda[lambda<=30])-1)
  result[lambda>30]=X[lambda>30]*log(lambda[lambda>30])-lambda[lambda>30]
  result[X==0]=0
  }
  else if(method.type=='PRAM'|method.type=='PAM')
  {result=X*log(lambda)-lambda }
  else {
    stop('method.type should be one of the following: tRex, tPAM, PRAM, PAM')
  }
  output=result[lower.tri(result)]
  # Negative log likelihood
  return(-sum(output))
}

gradient_poisson=function(P.V,Par,X,random,method.type)
{
  N=nrow(X)
  P=matrix(P.V,ncol=3)
  X=as.matrix(X)
  X[upper.tri(X,diag=TRUE)]=0
  x=P[,1]
  y=P[,2]
  z=P[,3]
  D.V=as.vector(sqrt(dist(P)))
  D.mat=matrix(1,N,N)
  D.mat[lower.tri(D.mat)]=D.V
  D.mat=t(D.mat)
  D.mat[lower.tri(D.mat)]=D.V
  D.V.nozero=D.V[D.V>0]
  D.V.nozero.min=min(D.V.nozero)
  D.mat[D.mat==0]=D.V.nozero.min/10000
  lambda=exp(Par[1]+random+Par[2]*log(D.mat))
  coef_m=matrix(0,N,N)
  if(method.type=='tRex'|method.type=='tPAM')
  {
  coef_m[lambda<=30]=(X[lambda<=30]/lambda[lambda<=30]-exp(lambda[lambda<=30])/(exp(lambda[lambda<=30])-1))*Par[2]*lambda[lambda<=30]/(D.mat[lambda<=30]^2)
  coef_m[lambda>30]=(X[lambda>30]/lambda[lambda>30]-1)*Par[2]*lambda[lambda>30]/(D.mat[lambda>30]^2)
  coef_m[X==0]=0
  }
  else if(method.type=='PRAM'|method.type=='PAM')
  {
   coef_m=(X/lambda-1)*Par[2]*lambda/(D.mat^2)  
  }
  else {
    stop('method.type should be one of the following: tRex, tPAM, PRAM, PAM')
    }
  temp=coef_m[lower.tri(coef_m)]
  coef_m=t(coef_m)
  coef_m[lower.tri(coef_m)]=temp
  
  gradient_x=(matrix(rep(x,N),nrow = N,byrow = TRUE)-x)*coef_m
  gradient_y=(matrix(rep(y,N),nrow = N,byrow = TRUE)-y)*coef_m
  gradient_z=(matrix(rep(z,N),nrow = N,byrow = TRUE)-z)*coef_m
  
  gradient_3d=cbind(colSums(gradient_x),colSums(gradient_y),colSums(gradient_z))
  # gradient of the negative log likelihood
  return(-as.vector(gradient_3d))
}  

estimate_parameters=function(P.V,X,method.type)
{
  P=matrix(P.V,ncol=3)
  X=as.matrix(X)
  N=ncol(X)
  X.V=X[lower.tri(X)]
  D.V=log(sqrt(as.vector(dist(P))))
  Wx=c()
  Wy=c()
  for(j in 1:(N-1)) 
  {  
    for(i in (j+1):N)
    {
      Wx=c(Wx,i)
      Wy=c(Wy,j)
    }
  }
  data=data.frame(Y=X.V,D=D.V,Wx=as.numeric(Wx),Wy=as.numeric(Wy))
  if(method.type=='tRex')
  {
  data.R=data[data$Y>0,]
  modelp=glmmTMB(Y~D+(1|Wx)+(1|Wy),data=data.R,family=truncated_poisson(link = "log"))
  models=summary(modelp)
  temp1=ranef(modelp)$cond$Wx
  temp2=ranef(modelp)$cond$Wy
  data1=data.frame(Wx=as.numeric(rownames(temp1)),as.vector(temp1))
  colnames(data1)=c('Wx','random1')
  data1$random1[is.na(data1$random1)]=0
  data2=data.frame(Wy=as.numeric(rownames(temp2)),as.vector(temp2))
  colnames(data2)=c('Wy','random2')
  data2$random2[is.na(data2$random2)]=0
  data=left_join(data,data1,by='Wx')
  data=left_join(data,data2,by='Wy')
  data$random=data$random1+data$random2
  data$random[is.na(data$random)]=0
  random_matrix=matrix(0,N,N)
  random_matrix[lower.tri(random_matrix)]=data$random
  out=list(models$coefficients$cond[,1],random_matrix)
  }
  else if(method.type=='PRAM')
  {
    data.R=data
    modelp=glmmTMB(Y~D+(1|Wx)+(1|Wy),data=data.R,family=poisson(link = "log"))
    models=summary(modelp)
    temp1=ranef(modelp)$cond$Wx
    temp2=ranef(modelp)$cond$Wy
    data1=data.frame(Wx=as.numeric(rownames(temp1)),as.vector(temp1))
    colnames(data1)=c('Wx','random1')
    data1$random1[is.na(data1$random1)]=0
    data2=data.frame(Wy=as.numeric(rownames(temp2)),as.vector(temp2))
    colnames(data2)=c('Wy','random2')
    data2$random2[is.na(data2$random2)]=0
    data=left_join(data,data1,by='Wx')
    data=left_join(data,data2,by='Wy')
    data$random=data$random1+data$random2
    data$random[is.na(data$random)]=0
    random_matrix=matrix(0,N,N)
    random_matrix[lower.tri(random_matrix)]=data$random
    out=list(models$coefficients$cond[,1],random_matrix)
  }
  else if(method.type=='tPAM')
  {
    data.R=data[data$Y>0,]
    modelp=glmmTMB(Y~D,data=data.R,family=truncated_poisson(link = "log"))
    models=summary(modelp)
    random_matrix=matrix(0,N,N)
    out=list(models$coefficients$cond[,1],random_matrix) 
  }
  else if(method.type=='PAM')
  {
    data.R=data
    modelp=glm(Y~D,data=data.R,family = poisson(link = "log"))
    models=summary(modelp)
    random_matrix=matrix(0,N,N)
    out=list(models$coefficients[,1],random_matrix)
    }
  
  else {
    stop('method.type should be one of the following: tRex, tPAM, PRAM, PAM')
  }
  return(out)
}

estimate_3D=function(P.int,Par,X,random,method.type)
{
  P.V=as.vector(P.int)
  out=optim(P.V,fn=objective_poisson,gr=gradient_poisson,Par=Par,X=X,random=random,method.type=method.type,method='BFGS',control= list(trace = 2,reltol=1e-12))
  struct=out$par
  struct3d=matrix(struct,ncol=3)
  return(list(struct3d,out$value))
}  
## estimated 3d coordinate given all the other information
## X is the HiC matrix
## P is the corresponding 3D coordinates
## Ignoring covariate for now
## W is the random effects

trexopt=function(P.int,X,max.iter,tol,method.type)
{
  X=as.matrix(X)
  N=nrow(X)
  P.int=as.vector(P.int)
  # Initial value for P
  neg_loglik=c()
  parameter_result=matrix(NA,max.iter,2)
  parameter_temp=c(1,-3)
  X.V=X[lower.tri(X)]
  random_temp=matrix(rnorm(N*N),N,N)
  for(i in 1:max.iter)
  {
    cat("3D iteration = ", i, "\n")
    struct_temp=estimate_3D(P.int,parameter_temp,X,random_temp,method.type)
    neg_loglik=c(neg_loglik,struct_temp[[2]])
    P.int=as.vector(struct_temp[[1]])
    cat("parameter iteration = ", i, "\n")
    temp_parameters=estimate_parameters(P.int,X,method.type)
    parameter_temp=temp_parameters[[1]]
    random_temp=temp_parameters[[2]]
    if(parameter_temp[2]>0)
    { parameter_temp[2]=runif(1,-4,-0.2) }
    parameter_result[i,]=parameter_temp
    if(i>1)
    {
      if(abs((neg_loglik[i]-neg_loglik[i-1])/neg_loglik[i])<tol)
      break
    }
  }
  return(list(struct=struct_temp[[1]],neg_loglik=neg_loglik,parameter=parameter_temp))
}

### method.type='tRex': truncated poisson, random effect
### method.type='tPAM': truncated Poisson, without random effect
### method.type='PRAM': poisson, with random effect
### method.type='PAM':  Poisson, without random effect


trex_opt=function(P.int=NULL,N.int=3,X,max.iter=50,tol=1e-4,method.type='tRex',parallel=FALSE,cpu=4)
{
X=as.matrix(X)
N=ncol(X)
P.int.list=list()
if(!(method.type %in% c('tRex','tPAM','PAM','PRAM')))
   { stop('method.type should be one of the following: tRex, tPAM, PRAM, PAM') }

if(is.null(P.int))
{
 for(i in 1:N.int)
 {
   P.int.list[[i]]=rnorm(3*N)
 }
}
else
{  
  if(length(P.int)!=N*3) {stop("The dimension of initial 3D structure and Hi-C matrix don't match!")}
  P.int.list[[1]]=P.int
  if(N.int>1)
  {
    for(i in 2:N.int)
    {
      P.int.list[[i]]=rnorm(3*N)
    }
  }
}
if(parallel==TRUE){
result_list=mclapply(P.int.list,trexopt,X,max.iter,tol,method.type,mc.cores = getOption("mc.cores", cpu))
}else{
    result_list=list()
    for(i in 1:length(P.int.list)){
      #cat("ith block",i,"\n")
      P.int.v=P.int.list[[i]]
      result_list[[i]]=trexopt(P.int.v,X,max.iter,tol,method.type)
    }
}  
### Find the optimal initial structure
neg.obj=c()
for(j in 1:length(result_list))
{
  temp=result_list[[j]][[2]]
  temp.m=temp[length(temp)]
  neg.obj=c(neg.obj,temp.m)
}
index=which.min(neg.obj)
return(result_list[[index]])
}


