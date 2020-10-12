source("SpLoc_functions_manuscript.R")

#Specify parameters for LME
tau2=0.5                          #Variance of noise
b.cov=matrix(c(3,0.5,0.5,0.2),2)  #Covariance of random intercept and slope

alpha0=1
alpha1=c(1,-1,0.5)
beta0=0.5
beta1=1

n.subj=50                         #Number of subjects

#Load nearest neighbor information (sparse matrix format in R)
NNmatLH=readRDS("NNmat_LH.rds") #Left hemisphere
NNmatRH=readRDS("NNmat_RH.rds") #Right hemisphere


# Define signals in Design 1
ind.signal.lh=which(NNmatLH[,28000]==1)
ind.signal.rh=which(NNmatRH[,28000]==1)

# Define signals in Design 2
# ind.signal.lh=c(which(NNmatLH[,14000]==1),which(NNmatLH[,15000]==1), which(NNmatLH[,16000]==1))
# ind.signal.rh=c(which(NNmatRH[,14000]==1),which(NNmatRH[,15000]==1), which(NNmatRH[,16000]==1))

# Define signals in Design 3
# ind.signal.lh=c(which(NNmatLH[,6500]==1),which(NNmatLH[,8500]==1), which(NNmatLH[,10500]==1),
#                 which(NNmatLH[,12500]==1),which(NNmatLH[,14500]==1))
# ind.signal.rh=c(which(NNmatRH[,6500]==1),which(NNmatRH[,8500]==1), which(NNmatRH[,10500]==1),
#                 which(NNmatRH[,12500]==1),which(NNmatRH[,14500]==1))

gamma=0.5

###Step 1 : generate simulation data
n.visits=sample(3:4, n.subj, replace = T)
time=NULL 
for (i in 1:n.subj){
  time=c(time, (0:(n.visits[i]-1))/2 ) 
}

X=matrix(rnorm(n.subj*3), n.subj) 
X.expand=X[rep(1:n.subj, n.visits),]
dx.status=c(rep(-1, n.subj/2),rep(1,n.subj/2)) 
dx.expand=rep(dx.status, n.visits)

pred=c(alpha0+X.expand%*%alpha1+dx.expand*beta0+time*beta1)
ymat1=matrix(NA,sum(n.visits),2329) 
ymat2=matrix(NA,sum(n.visits),2332) 

#Generate data for the left and right hemispheres
for (j in 1:2329){ 
  b=mvrnorm(n.subj,c(0,0), b.cov) 
  b.expand=b[rep(1:50, n.visits),]
  epsilon=rnorm(sum(n.visits),0,sqrt(tau2))
  ymat1[,j]=pred+b.expand[,1]+time*b.expand[,2]+epsilon 
}

for (j in 1:2332){
  b=mvrnorm(n.subj,c(0,0), b.cov)
  b.expand=b[rep(1:50, n.visits),]
  epsilon=rnorm(sum(n.visits),0,sqrt(tau2))
  ymat2[,j]=pred+b.expand[,1]+time*b.expand[,2]+epsilon
}

#Add signals
ymat1[,ind.signal.lh]=ymat1[,ind.signal.lh]+matrix(rep(gamma*dx.expand*time, length(ind.signal.lh)), sum(n.visits) )
ymat2[,ind.signal.rh]=ymat2[,ind.signal.rh]+matrix(rep(gamma*dx.expand*time, length(ind.signal.rh)), sum(n.visits) )


###Step 2: Fitting LME and conducing permutation
#Fitting univariate LME model to each vertex
fit1=lmefit(ymat1, X.expand, time, dx.status,n.visits)
fit2=lmefit(ymat2, X.expand, time, dx.status,n.visits)

#Conducting SpLoc
sc=scoretest.perm(fit1,fit2, NNmatLH,NNmatRH, time, dx.status,n.visits,n.perm = 5000)

###Step 3: Select signal clusters
thres=quantile(apply(cbind(sc$permMatLH^2,sc$permMatRH^2),1,max),0.95)  #Thresholds controlling FWER
cl1=ClusterSearch(sc$scoreLH,thres,NNmatLH)  #Detected clusters in the left hemisphere
cl2=ClusterSearch(sc$scoreRH,thres,NNmatRH)  #Detected clusters in the right hemisphere
