library(lme4)
library(lmerTest)
library(Matrix)
library(MASS)
library(mutoss)

lmefit=function(ymat, Z, time, dx.status, n.visits){
  p=ncol(ymat)
  n=nrow(ymat)
  n.subj=length(n.visits)
  Subject=rep(paste0("Subj",1:n.subj),n.visits)
  dx.expand=rep(dx.status, n.visits)
  time2=time^2
  residMat=matrix(NA, n, p)
  for (j in 1:p){
    fit=lmer(ymat[,j] ~ Z+time+dx.expand+(time|Subject), 
             control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-5)))
    sigma2=attr(VarCorr(fit),"sc")^2
    residMat[,j]=residuals(fit)/sigma2
  }
  return(residMat)
}


vlme=function(ymat, Z, time, dx.status, n.visits){
  p=ncol(ymat)
  n=nrow(ymat)
  n.subj=length(n.visits)
  Subject=rep(paste0("Subj",1:n.subj),n.visits)
  dx.expand=rep(dx.status, n.visits)
  pvecs=rep(NA, p)
  for (j in 1:p){
    fit=lmer(ymat[,j] ~ Z+time+dx.expand+time*dx.expand+(time|Subject), 
             control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-5)))
    pvecs[j]=anova(fit)[4,6]
  }
  return(pvecs)
}

scoretest.perm=function(residMatLH, residMatRH, 
                        NNmatLH, NNmatRH,
                        time, dx.status, n.visits, n.perm=10000){
  n=nrow(residMatLH)
  p1=ncol(residMatLH)
  p2=ncol(residMatRH)
  ind1=which(residMatLH[1,]!=0)
  ind2=which(residMatRH[1,]!=0)
  if (length(ind1)==0){
    pval=stat=NA
    perms=rep(NA, n.perm)
  } else{
    residMatLH=residMatLH[,ind1,drop=F]
    residMatRH=residMatRH[,ind2,drop=F]
    time.residm1=matrix(NA, length(n.visits),p1)
    time.residm2=matrix(NA, length(n.visits),p2)
    
    for (i in 1:length(n.visits)){
      if (i==1){start=1} else{ start=sum(n.visits[1:(i-1)])+1}
      end=sum(n.visits[1:i])
      time.sub=time[start:end]
      time.residm1[i,]=apply(residMatLH[start:end,]*matrix(rep(time.sub, p1), length(time.sub)),2,sum)
      time.residm2[i,]=apply(residMatRH[start:end,]*matrix(rep(time.sub, p2), length(time.sub)),2,sum)
    }
    
    dx.mat1=matrix(rep(dx.status, p1),length(n.visits))
    dx.mat2=matrix(rep(dx.status, p2),length(n.visits))
    
    score1=as.numeric(crossprod(NNmatLH, apply(dx.mat1*time.residm1,2,sum)))
    score2=as.numeric(crossprod(NNmatRH, apply(dx.mat2*time.residm2,2,sum)))
    
    permMat1=matrix(NA, n.perm, ncol(NNmatLH))
    permMat2=matrix(NA, n.perm, ncol(NNmatRH))
    
    for (perm in 1:n.perm){
      perm.index=sample(length(n.visits))
      permMat1[perm,]=as.numeric(crossprod(NNmatLH, apply(dx.mat1[perm.index,]*time.residm1,2,sum)))
      permMat2[perm,]=as.numeric(crossprod(NNmatRH, apply(dx.mat2[perm.index,]*time.residm2,2,sum)))
    }
    
    for (j in 1:ncol(NNmatLH)){
      permvec=permMat1[,j]
      sdp=sd(permvec)
      permMat1[,j]=((permvec-mean(permvec)) /sdp)
      score1[j]=(score1[j]/sdp)
    }
    
    for (j in 1:ncol(NNmatRH)){
      permvec=permMat2[,j]
      sdp=sd(permvec)
      permMat2[,j]=((permvec-mean(permvec)) /sdp)
      score2[j]=(score2[j]/sdp)
    }
    
    perms=apply( (cbind(permMat1, permMat2))^2,1,max)
    stat=max(c(score1,score2)^2)
    pval=mean(perms>stat )
  }
  
  return(list(scoreLH=score1^2, scoreRH=score2^2,
              permMatLH=permMat1,permMatRH=permMat2,
              stat=stat,perms=perms, pval=pval))
}

ClusterSearch=function(tstat, thres, NNmat){
  sig.ind=which(tstat>thres)
  tstat=tstat[sig.ind]
  NNmat=NNmat[,sig.ind,drop=F]
  
  if (ncol(NNmat)==0){bool=F}else{bool=T}
  
  clust=1
  sig=NULL
  while(bool){
    ind=which(tstat==max(tstat))[1]
    sig.vertices=which(NNmat[,ind]!=0)
    sig=c(sig, sig.vertices)
    out.set=which(apply(NNmat[sig.vertices,,drop=F],2,max)>0)
    
    tstat=tstat[-out.set]
    NNmat=NNmat[,-out.set,drop=F]
    
    print(paste0("Cluster ", clust, " with ", length(sig.vertices), " vertices is selected (total=", length(sig), ")" ))
    clust=clust+1
    
    if (length(tstat)==0){bool=F}
  }
  
  return(sig)
}
