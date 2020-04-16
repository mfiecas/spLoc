library(MASS)
library(mvnfast)
library(Matrix)

frob=function(X){sum(X^2, na.rm=T)}

mySample <- function(n,x,repl=5000){
  nx=choose(n,x)
  if (nx<repl){repl=nx}
  rowsample <- replicate(repl,sort(sample(n,x,FALSE)),simplify=FALSE)
  while(length(unique(rowsample))<repl){
    rowsample <- unique(rowsample)
    rowsample <- c(rowsample,
                   replicate(repl-length(rowsample),
                             sort(sample(n,x,FALSE)),simplify=FALSE))
  }
  return(do.call(cbind,rowsample))
}

buildGrid=function(n=10, jitter=0.0001){
  coords.x=kronecker(1:n, rep(1,n))/n
  coords.y=kronecker(rep(1,n),1:n)/n
  coords=cbind(coords.x,coords.y)
  coords.jitter=coords+runif(prod(dim(coords)),-jitter,jitter )
  distmat=as.matrix(dist(cbind(coords.x,coords.y)))
  return(list(coords=coords, coords.jitter=coords.jitter, distmat=distmat))
}

buildNNmat=function(nn=1, distMat){
  n=nrow(distMat)
  out=matrix(0, n,n)
  for (i in 1:n){
    ind =which(distMat[,i]<=sort(distMat[,i],partial=nn)[nn])
    out[ind,i]=1/(nn)
  }
  return(out)
}

spLocST.fit=function(y,X,dx.status,time, distMat, n.visits, 
                    n.iter=10,m.iter=20, eps=1e-3, eps2=1e-2, 
                    start=NULL,...){
  if (sum(n.visits)!=nrow(y)) stop("n.visits do not match")
  distMat=distMat^2
  m=ncol(y); p=ncol(X)
  
  iter=0; bool=T
  betamat=NULL; varcompsmat=NULL
  residmat=NULL
  augX=augdx=NULL
  for (j in 1:length(n.visits)){
    for (k in 1:n.visits[j]){
      augX=rbind(augX, X[j,])
      augdx=c(augdx, dx.status[j])
    }
  }
  augX2=cbind(augX,augdx)
  
  beta.old=rep(0, p+1+m)
  varcomps.old=0
  while (bool){
    if (iter==0){
      if (is.null(start)){
        
        beta=as.numeric(updateMeansInit(y,cbind(X,dx.status),n.visits,time,distMat))
        betamat=cbind(betamat,beta)
        resid=matrix(NA,sum(n.visits),m)
        
        for (j in 1:m){
          beta1=beta[1:(p+1)]
          beta2=beta[-(1:(p+1))]
          resid[,j]=y[,j]-augX2%*%beta1-time*beta2[j]
        }
        residmat=cbind(residmat,c(resid))
        residT=t(resid)
        phi.hat.old=phi.hat=5
        psi.hat.old=psi.hat=1
        bool2=T;iter2=0
        while (bool2){
          vcomps=c(covRegC(residT,n.visits, time,distMat, phi.hat, psi.hat));
          phi.hat=optim(c(phi.hat), varcomps=c(psi.hat, vcomps[1],vcomps[2]), variomatSpatial, residT=residT,
                        lower=0.000001, upper=200,distMat=distMat, time=time,nvisits=n.visits,
                        method = "L-BFGS-B")$par[1]
          psi.hat=optim(c(psi.hat), varcomps=c(phi.hat, vcomps[1],vcomps[2]), variomatTemporal, residT=residT,
                        lower=0.000001, upper=1,distMat=distMat, time=time,nvisits=n.visits,
                        method = "L-BFGS-B")$par[1]
          if (iter2==m.iter){bool2=F}
          else if (abs(phi.hat-phi.hat.old)<eps2){bool2=F}
          else if (abs(psi.hat-psi.hat.old)<eps2){bool2=F}
          else{
            iter2=iter2+1
            phi.hat.old=phi.hat
            psi.hat.old=psi.hat
          }
        }
        varcomps=c(phi.hat, psi.hat,vcomps)    
        varcompsmat=cbind(varcompsmat,varcomps)
      }
      else{
        beta.start=fit$betamat
        varcomps.start=fit$varcompsmat
        beta=fit$beta; beta.old=fit$betamat[,ncol(fit$betamat)-1]
        varcomps=fit$varcomps; varcomps.old=fit$varcompsmat[,ncol(fit$varcompsmat)-1]
        phi.hat=varcomps[1]
        psi.hat=varcomps[2]
      }
    }
    else{
      beta=as.numeric(updateMeansST(y,cbind(X,dx.status),n.visits,time,distMat,
                                  phi=varcomps[1],psi=varcomps[2],
                                  sigma_s=varcomps[3], sigma_e=varcomps[4]))
      
      resid=matrix(NA,sum(n.visits),m)
      
      
      for (j in 1:m){ 
        beta1=beta[1:(p+1)]
        beta2=beta[-(1:(p+1))]
        resid[,j]=y[,j]-augX2%*%beta1-time*beta2[j] 
      }
      residmat=cbind(residmat,c(resid))
      
      residT=t(resid)
      bool2=T;iter2=0
      while (bool2){
        vcomps=c(covRegC(residT,n.visits, time,distMat, phi.hat, psi.hat));
        phi.hat=optim(c(phi.hat), varcomps=c(psi.hat, vcomps[1],vcomps[2]), variomatSpatial, residT=residT,
                      lower=0.0001, upper=100,distMat=distMat, time=time,nvisits=n.visits,
                      method = "L-BFGS-B")$par[1]
        psi.hat=optim(c(psi.hat), varcomps=c(phi.hat, vcomps[1],vcomps[2]), variomatTemporal, residT=residT,
                      lower=0.0001, upper=1,distMat=distMat, time=time,nvisits=n.visits,
                      method = "L-BFGS-B")$par[1]
        if (iter2==m.iter){bool2=F}
        else if (abs(phi.hat-phi.hat.old)<eps2){bool2=F}
        else if (abs(psi.hat-psi.hat.old)<eps2){bool2=F}
        else{
          iter2=iter2+1
          phi.hat.old=phi.hat
          psi.hat.old=psi.hat
        }
      }
      varcomps=c(phi.hat,psi.hat,vcomps)
      betamat=cbind(betamat,beta)
      varcompsmat=cbind(varcompsmat, varcomps)
    }
    
    if (iter==n.iter){bool=F}
    else if ( (frob(beta-beta.old)/length(beta))<eps ){bool=F}
    else{ 
      beta.old=beta
      varcomps.old=varcomps
      iter=iter+1
    }
  }
  return(list(beta=as.numeric(beta),varcomps=varcomps,iter=iter,
              resid=resid,betamat=betamat, varcompsmat=varcompsmat,residmat=residmat))
}


spLocST.perm=function(residMat,X,time,dx.status,distMat, n.visits, varcomps,nn.set=NULL, n.perm=1000,alpha=0.05){
  m=nrow(distMat)
  if (is.null(nn.set)){ nn.set=unique(c(2^(c(0:floor(log(m)/log(2)) ) ),m )) }
  nn.set=nn.set[1:5]
  
  n.nn.set=length(nn.set)
  distMat=distMat^2
  
  nnList=matrix(NA, m, n.nn.set*m)
  for (k in 1:n.nn.set){
    if (k==1){
      start=1;end=m
    }else{
      start=(k-1)*m+1;end=k*m
    }
    nnList[,start:end]= (buildNNmat(nn=nn.set[k], distMat))
  }
  
  myPerm=mySample(length(dx.status),sum(dx.status==-1), n.perm)-1
  pval=c(splocst_modelpermC(resid=residMat, distMat=distMat,nnList = nnList,
                           time = time,dxStatus = dx.status,permMat = myPerm,
                           varcomps=varcomps,nvisits = n.visits,alpha=alpha))
  
  return(list(pval=pval[1], thres=pval[2],tstat=matrix(pval[-(1:2)],m),nn.set=nn.set,alpha=alpha))
}

spLocT.fit=function(y,X,dx.status,time, n.visits, n.iter=10,eps=1e-3, distMat=distMat,
                   start=NULL,...){
  m=ncol(y); p=ncol(X)
  
  iter=0; bool=T
  betamat=NULL; varcompsmat=NULL; residmat=NULL
  augX=augdx=NULL
  for (j in 1:length(n.visits)){
    for (k in 1:n.visits[j]){
      augX=rbind(augX, X[j,])
      augdx=c(augdx, dx.status[j])
    }
  }
  augX2=cbind(augX,augdx)
  
  beta.old=rep(0, p+1+m)
  varcomps.old=0
  while (bool){
    if (iter==0){
      beta=as.numeric(updateMeansInit(y,cbind(X,dx.status),n.visits,time,distMat))
      betamat=cbind(betamat,beta)
      resid=matrix(NA,sum(n.visits),m)
      
      q=length(beta)
      
      for (j in 1:m){
        beta1=beta[1:(p+1)]
        beta2=beta[-(1:(p+1))]
        resid[,j]=y[,j]-augX2%*%beta1-time*beta2[j]
      }
      residmat=cbind(residmat, c(t(resid)))
      sigma2=frob(resid)/(prod(dim(resid))-q)
      
      rrt.summ=0
      for (i in 1:length(n.visits)){
        visits=n.visits[i]
        if (i==1){ ind=1:visits}
        else{ ind=(sum(n.visits[1:(i-1)])+1):(sum(n.visits[1:(i-1)])+visits) }
        resid.sub=resid[ind,]
        rrt=tcrossprod(resid.sub);diag(rrt)=0
        rrt.summ=rrt.summ+sum(rrt)/2
      }
      phi=rrt.summ/(sum( n.visits*(n.visits-1)*m/2)-q)/sigma2
      varcomps=c(sigma2, phi)
      varcompsmat=cbind(varcompsmat, varcomps)
    }
    else{
      beta=as.numeric(updateMeansT(y,cbind(X,dx.status),n.visits,time,
                                   sigma_e=varcomps[1], phi=varcomps[2]))
      
      resid=matrix(NA,sum(n.visits),m)
      for (j in 1:m){ 
        beta1=beta[1:(p+1)]
        beta2=beta[-(1:(p+1))]
        resid[,j]=y[,j]-augX2%*%beta1-time*beta2[j] 
      }
      residmat=cbind(residmat, c(t(resid)))
      
      sigma2=frob(resid)/(prod(dim(resid))-q)
      
      rrt.summ=0
      for (i in 1:length(n.visits)){
        visits=n.visits[i]
        if (i==1){ ind=1:visits}
        else{ ind=(sum(n.visits[1:(i-1)])+1):(sum(n.visits[1:(i-1)])+visits) }
        resid.sub=resid[ind,]
        rrt=tcrossprod(resid.sub);diag(rrt)=0
        rrt.summ=rrt.summ+sum(rrt)/2
      }
      phi=rrt.summ/(sum( n.visits*(n.visits-1)*m/2)-q)/sigma2
      
      varcomps=c(sigma2, phi)
      betamat=cbind(betamat,beta)
      varcompsmat=cbind(varcompsmat, varcomps)
    }
    
    if (iter==n.iter){bool=F}
    else if ( (frob(beta-beta.old)/length(beta))<eps ){bool=F}
    else{ 
      beta.old=beta
      varcomps.old=varcomps
      iter=iter+1
    }
  }
  return(list(beta=as.numeric(beta),varcomps=varcomps,iter=iter,
              resid=resid,betamat=betamat, varcompsmat=varcompsmat,residmat=residmat))
}

spLocT.perm=function(residMat,X,time,dx.status,distMat, n.visits, varcomps,nn.set=NULL, n.perm=1000,alpha=0.05/70){
  m=nrow(distMat)
  if (is.null(nn.set)){ nn.set=unique(c(2^(c(0:floor(log(m)/log(2)) ) ),m )) }
  
  n.nn.set=length(nn.set)
  
  nnList=matrix(NA, m, n.nn.set*m)
  for (k in 1:n.nn.set){
    if (k==1){
      start=1;end=m
    }else{
      start=(k-1)*m+1;end=k*m
    }
    nnList[,start:end]= (buildNNmat(nn=nn.set[k], distMat))
  }
  
  myPerm=mySample(length(dx.status),sum(dx.status==-1), n.perm)-1
  
  pval=c(sploct_modelpermC(resid=residMat, distMat=distMat,nnList = nnList,
                          time = time,dxStatus = dx.status,permMat = myPerm,
                          varcomps=varcomps,nvisits = n.visits,alpha=alpha))
  
  return(list(pval=pval[1], thres=pval[2], tstat=matrix(pval[-(1:2)],m), nn.set=nn.set,alpha=alpha))
}

selectRegions=function(permfit,distMat){
  tstat=permfit$tstat
  nn.set=permfit$nn.set
  alpha=permfit$alpha
  thres=permfit$thres
  pval=permfit$pval
  if (pval<alpha){
    m=nrow(distMat)
    out=NULL
    Iset=which(tstat>=thres, arr.ind=T)
    bool=T
    while (bool){
      ind.min=which(tstat[Iset]==max(tstat[Iset]))
      if (length(ind.min)>1){
        ind.min=ind.min[1]
      }
      ind.min.loc=Iset[ind.min,1]; ind.min.nn=nn.set[Iset[ind.min,2]]
      ind.vert=as.numeric(which(distMat[ind.min.loc,]<=sort(distMat[ind.min.loc,])[ind.min.nn]))
      out=c(out, ind.vert)
      
      Iset=Iset[-ind.min,,drop=F]
      include.set=NULL
      if (nrow(Iset)!=0){
        for (j in 1:nrow(Iset)){
          vert=Iset[j,1]; vert.nn=nn.set[Iset[j,2]]
          if (!vert%in%ind.vert){
            sub.neighbors=as.numeric(which(distMat[vert,]<=sort(distMat[vert,])[vert.nn]))
            if (length(setdiff(ind.vert, sub.neighbors))==length(ind.vert)){
              include.set=c(include.set, j)
            }
          }
        }
      }
      
      #update Iset
      if (is.null(include.set)){ bool=F }
      else{    Iset=Iset[include.set,,drop=F]}
    }    
  } else{
    out=NULL
  }
  
  
  return(out)
}
