library(Rcpp)
library(RcppArmadillo)

sourceCpp("spLoc_support.cpp")
source("spLoc.R")

set.seed(1131)


##STEP 1: generate simulated data
n.subj=50  #number of subjects
s=100      #number of vertices 
alpha=0.05 #FWER

sigma2=0.3
tau2=0.1     #nugget (iid noise)
phi=1;       #parameter for spatial correlation
psi=0.5      #temporal correlation

gamma=0.1    #effect size
gammamat=matrix(0,10,10) #generate a (10 x 10) grid
gammamat[4:6,4:6]=gamma

gamma.vec=c(gammamat)

grid=buildGrid(n = sqrt(s))
distmat=as.matrix(dist(grid$coords.jitter))
n.visits=sample(2:4, n.subj, replace=T)
n=sum(n.visits)  #total number of visits
dx.status=c(rep(-1, n.subj/2),rep(1,n.subj/2))

X1=runif(n.subj)
X2=rbinom(n.subj, 1,0.5)
X3=X1*X2
X=cbind(1,X1,X2,X3)
y=matrix(NA, n, s)

augX=NULL
augdx=NULL
for (j in 1:length(n.visits)){
  for (k in 1:n.visits[j]){
    augX=rbind(augX, X[j,])
    augdx=c(augdx, dx.status[j])
  }
}

time.vec=NULL
for (i in 1:length(n.visits)){
  visits=n.visits[i]
  time=c(0, (1:(visits-1))*0.3 )
  time.vec=c(time.vec, time)
}

#generate covariates
beta.covariate=c(1,2,2,-1)
beta.dx=rnorm(1)
beta.time=rnorm(s)
beta=c(beta.covariate, beta.dx, beta.time)
Xbeta=augX%*%beta.covariate
Xbetamat=kronecker(Xbeta, t(rep(1,s)))
for (j in 1:s){
  Xbetamat[,j]=Xbetamat[,j]+augdx*beta.dx+time.vec*beta.time[j]
}

resid=matrix(rnorm(n*s, 0, sqrt(tau2)),n)

for (j in 1:length(n.visits)){
  visits=n.visits[j]
  if (j==1) {ind=1:visits}
  else {
    csum=sum(n.visits[1:(j-1)])
    ind= (csum+1):(csum+visits)
  }
  time=time.vec[ind]
  sptemp=t(matrix(mvrnorm(1, rep(0,s*visits),sigma2*kronecker(generateTmat(time,psi), exp(-phi*distmat^2)) ),s) )
  resid[ind,]=resid[ind,]+sptemp
  
  for (v in 1:s){
    resid[ind,v]=resid[ind,v]+time*dx.status[j]*gamma.vec[v]
  }
}

for (j in 1:s){
  y[,j]=Xbetamat[,j]+resid[,j]
}

##Step 2: fitting spLoc-ST and spLoc-T
#fitting spLoc-ST
fit1=spLocST.fit(y,X, dx.status=dx.status,
                time=time.vec,distMat=distmat,n.visits=n.visits, 
                n.iter = 20,m.iter=30,eps = 1e-2,eps2 = 5e-2)

pvecs1=spLocST.perm(residMat=fit1$resid,X = X,dx.status = dx.status, n.visits = n.visits,distMat = distmat,varcomps = fit1$varcomps, time = time.vec,n.perm = 20000, alpha=alpha)
pval1=pvecs1$pval    #p-value for the null
regions1=selectRegions(pvecs1, distmat) #selected signal regions NULL if there's none.

#fitting sploc-T
fit2=spLocT.fit(y=y,X=X,dx.status=dx.status,time=time.vec, n.visits=n.visits, distMat=distmat)
pvecs2=spLocT.perm(residMat=fit2$resid,X=X,time=time.vec,dx.status=dx.status,distMat=distmat, n.visits=n.visits, varcomps=fit2$varcomps,n.perm=20000, alpha=alpha)
pval2=pvecs2$pval    #p-value for the null
regions2=selectRegions(pvecs2, distmat) #selected signal regions. NULL if there's none.

