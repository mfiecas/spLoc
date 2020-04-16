#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// [[Rcpp::export]]
double frobC(arma::mat X){
  double out=0;
  int p=X.n_rows;
  int q=X.n_cols;
  
  for (int i=0; i<p; ++i){
    for (int j=0; j<q; ++j){
      out=out+pow(X(i,j),2);
    }
  }
  return out;
}

// [[Rcpp::export]]
arma::mat diagC(int a) {
  arma::mat out(a,a); out.fill(0);
  for (int i=0; i<a; ++i){
    out(i,i)=1;
  }
  return out;
}

// [[Rcpp::export]]
double meanprop(arma::vec a, double b){
  double out=0;
  int n=a.size();
  for (int i=0; i<n; ++i){
    if (a(i)>=b){
      out=out+1;
    }
  }
  out=out/(n+1);
  return out;
}

// [[Rcpp::export]]
arma::mat STInvC(arma::mat Us, arma::vec ds, arma::mat T2mat, double sigma2, double tau2){
  arma::mat Ut;
  arma::vec dt;
  arma::mat Vt;
  
  svd_econ(Ut, dt,Vt, T2mat);
  
  arma::mat dst=ds*dt.t();
  dst=1/(dst*sigma2+tau2);
  arma::mat dstmat=diagmat(vectorise(dst));
  
  arma::mat Ust=kron(Ut, Us);
  arma::mat out=Ust*dstmat*Ust.t();
  
  return out;
}

// [[Rcpp::export]]
arma::vec avg_rank(arma::vec x) {
  arma::uvec w = arma::stable_sort_index(x, "descend");
  R_xlen_t sz = x.size();
  arma::vec r(sz);

  for (R_xlen_t n, i = 0; i < sz; i += n) {
    n = 1;
    while (i + n < sz && x[w[i]] == x[w[i + n]]) ++n;
    for (R_xlen_t k = 0; k < n; k++) {
      // r[w[i + k]] = i + (n + 1) / 2.;
      r[w[i + k]] = i + (n ) ;
    }
  }
  return r;
}

// [[Rcpp::export]]
double quantileC(arma::vec Tstatvec, double alpha){
  int n=Tstatvec.size();
  arma::vec avgrank=avg_rank(Tstatvec);
  int thres=floor(n*alpha);
  int index=0;
  for (int i=0; i<n; ++i){
    if (avgrank(i)==thres){
      index=i;
      break;
    }
  }
  double out=Tstatvec(index);
  return out;
}

// [[Rcpp::export]]
arma::mat solveC(arma::mat X){
	arma::mat Xinv=inv_sympd(X);
	return Xinv;
}

// [[Rcpp::export]]
arma::mat solveC2(arma::mat X){
  arma::mat Xinv=X.i();
  return Xinv;
}

// [[Rcpp::export]]
arma::mat solveCstable(arma::mat X){
  int n=X.n_cols;
  arma::mat out(n,n); out.fill(0);
  arma::mat u;
  arma::vec d;
  arma::mat v;
  svd(u, d, v, X);
  
  for (int i=0; i<d.size(); ++i){
    if (d(i)>0){
      out=out+u.col(i)*v.col(i).t()/d(i);
    }
  }
  
  return out;
}

// [[Rcpp::export]]
arma::mat mykron(arma::vec time, int visits){
  arma::mat iden=diagC(visits);
  arma::mat out=kron(time, iden);
  return out;
}

// [[Rcpp::export]]
double computeABA(arma::mat X, arma::vec y){
	double out;
	arma::mat Xy=X*y;
	out=sum(y%Xy);
	return out;
}

// [[Rcpp::export]]
double computeABAstable(arma::vec Bd, arma::mat Bv, arma::vec y){
  double out=0;
  int d=Bv.n_cols;
  arma::vec Bvty=Bv.t()*y;
  for (int i=0; i<d; ++i){
    if (Bd(i)>0){
      out=out+Bd(i)*pow(Bvty(i),2);
    }
  }
  return out;
}

// [[Rcpp::export]]
arma::mat computeABAmat(arma::mat A, arma::mat B){
  arma::mat BA=B*A;
  int p=A.n_cols;
  arma::mat out(p,p);
  for (int i=0; i<p; ++i){
    for (int j=0; j<p;++j){
      out(i,j)=sum(A.col(i)%BA.col(j));
    }
  }
  return out;
}

// [[Rcpp::export]]
arma::mat computeABCmat(arma::mat A, arma::mat B, arma::mat C){
  arma::mat BC=B*C;
  int p=A.n_cols;
  int q=C.n_cols;
  arma::mat out(p,q);
  for (int i=0; i<p; ++i){
    for (int j=0; j<q;++j){
      out(i,j)=sum(A.col(i)%BC.col(j));
    }
  }
  return out;
}


// [[Rcpp::export]]
arma::mat generateTmat(arma::vec time, double psi){
  int visits=time.size();
  arma::mat mat(visits,visits);mat.fill(1);
  arma::mat V=psi*mat+(1-psi)*diagC(visits);
  return V;
}

// [[Rcpp::export]]
arma::mat generateSmat(double phi, arma::mat distmat){
  arma::mat distmat2=exp(-phi*distmat);
  return distmat2;
}

// [[Rcpp::export]]
arma::vec updateMeansST(arma::mat y, arma::mat augX, arma::vec nvisits, arma::vec time,arma::mat distmat,
                      double phi, double psi, double sigma_s,double sigma_e){
  int start;
  int end;
  int p=augX.n_cols;
  int s=y.n_cols;
  int nsubj=nvisits.size();
  arma::mat XtX(p,p); XtX.fill(0);
  arma::mat XtW(p,s); XtW.fill(0);
  arma::mat WtX(s,p); WtX.fill(0);
  arma::mat WtW(s,s); WtW.fill(0);
  arma::vec Xty(p); Xty.fill(0);
  arma::vec Wty(s); Wty.fill(0);
  
  arma::vec time0;
  arma::rowvec subX1;
  arma::mat subX;
  arma::mat subW;
  arma::mat t1mat;
  arma::mat V;
  arma::mat Vinv;
  arma::mat subresid;
  arma::mat subresidvec;
  arma::vec VinvResidvec;
  arma::mat spMat=generateSmat(phi, distmat);
  arma::mat spU;
  arma::vec spD;
  arma::mat spV;
  svd(spU,spD,spV, spMat);
  
  for (int i=0; i<nsubj; ++i){
    int visits=nvisits(i);
    int svisits=s*visits;
    if (i==0){
      start=0;
      end=visits-1;
    } else{
      start=sum(nvisits.subvec(0, i-1));
      end=sum(nvisits.subvec(0,i))-1;
    }
    
    time0=time.subvec(start, end);
    arma::mat mat(svisits,1);mat.fill(1);
    subX1=augX.row(i);
    subX=kron(mat,subX1);
    subW=mykron(time0, s);
    t1mat=generateTmat(time0, psi);
    Vinv=STInvC(spU,spD,t1mat,sigma_s,sigma_e);
    double Vinvsum=accu(Vinv);
    XtX=XtX+Vinvsum*(subX1.t()*subX1); //computeABAmat(subX, Vinv);
    XtW=XtW+computeABCmat(subX, Vinv, subW);
    WtX=XtW.t();
    WtW=WtW+computeABAmat(subW, Vinv);
    
    subresid=y.rows(start,end);
    subresidvec=vectorise(subresid.t());
    VinvResidvec=Vinv*subresidvec;
    Xty=Xty+subX.t()*VinvResidvec;
    Wty=Wty+subW.t()*VinvResidvec;
  }
  
  arma::mat XtXinv=solveC2(XtX);
  arma::mat out22=solveC2(WtW-WtX*XtXinv*XtW);
  arma::mat out11=XtXinv+XtXinv*XtW*out22*WtX*XtXinv;
  arma::mat out12=-(XtXinv*XtW*out22);
  arma::mat out21=-(out22*WtX*XtXinv);
  
  arma::vec beta1=out11*Xty+out12*Wty;
  arma::vec beta2=out21*Xty+out22*Wty;
  arma::mat beta=join_cols(beta1,beta2);
  return beta;
}


// [[Rcpp::export]]
arma::vec updateMeansT(arma::mat y, arma::mat augX, arma::vec nvisits, arma::vec time,
                      double phi,double sigma_e){
  int start;
  int end;
  int p=augX.n_cols;
  int s=y.n_cols;
  int nsubj=nvisits.size();
  arma::mat XtX(p,p); XtX.fill(0);
  arma::mat XtW(p,s); XtW.fill(0);
  arma::mat WtX(s,p); WtX.fill(0);
  arma::mat WtW(s,s); WtW.fill(0);
  arma::vec Xty(p); Xty.fill(0);
  arma::vec Wty(s); Wty.fill(0);
  
  arma::vec time0;
  arma::rowvec subX1;
  arma::mat subX;
  arma::mat subW;
  arma::mat t1mat;
  arma::mat V;
  arma::mat Vinv;
  arma::mat subresid;
  arma::mat subresidvec;
  arma::vec VinvResidvec;

  for (int i=0; i<nsubj; ++i){
    int visits=nvisits(i);
    int svisits=s*visits;
    if (i==0){
      start=0;
      end=visits-1;
    } else{
      start=sum(nvisits.subvec(0, i-1));
      end=sum(nvisits.subvec(0,i))-1;
    }
    
    time0=time.subvec(start, end);
    arma::mat mat1(visits,visits);mat1.fill(1);
    V=generateTmat(time0,phi);

    arma::mat mat(svisits,1);mat.fill(1);
    subX1=augX.row(i);
    subX=kron(mat,subX1);
    subW=mykron(time0, s);
    Vinv=kron(solveC(V),diagC(s));
    double Vinvsum=accu(Vinv);
    XtX=XtX+Vinvsum*(subX1.t()*subX1); 
    XtW=XtW+computeABCmat(subX, Vinv, subW);
    WtX=XtW.t();
    WtW=WtW+computeABAmat(subW, Vinv);
    
    subresid=y.rows(start,end);
    subresidvec=vectorise(subresid.t());
    VinvResidvec=Vinv*subresidvec;
    Xty=Xty+subX.t()*VinvResidvec;
    Wty=Wty+subW.t()*VinvResidvec;
  }
  
  arma::mat XtXinv=solveC2(XtX);
  arma::mat out22=solveC2(WtW-WtX*XtXinv*XtW);
  arma::mat out11=XtXinv+XtXinv*XtW*out22*WtX*XtXinv;
  arma::mat out12=-(XtXinv*XtW*out22);
  arma::mat out21=-(out22*WtX*XtXinv);
  
  arma::vec beta1=out11*Xty+out12*Wty;
  arma::vec beta2=out21*Xty+out22*Wty;
  arma::mat beta=join_cols(beta1,beta2);
  return beta;
}

// [[Rcpp::export]]
arma::vec updateMeansInit(arma::mat y, arma::mat augX, arma::vec nvisits, arma::vec time,arma::mat distmat){
  int start;
  int end;
  int p=augX.n_cols;
  int s=y.n_cols;
  int nsubj=nvisits.size();
  arma::mat subX; subX.fill(0);
  arma::mat XtX(p,p); XtX.fill(0);
  arma::mat XtW(p,s); XtW.fill(0);
  arma::mat WtX(s,p); WtX.fill(0);
  arma::mat WtW(s,s); WtW.fill(0);
  arma::vec Xty(p); Xty.fill(0);
  arma::vec Wty(s); Wty.fill(0);
  
  for (int i=0; i<nsubj; ++i){
    int visits=nvisits(i);
    int svisits=s*visits;
    if (i==0){
      start=0;
      end=visits-1;
    } else{
      start=sum(nvisits.subvec(0, i-1));
      end=sum(nvisits.subvec(0,i))-1;
    }
    
    arma::vec time0=time.subvec(start, end);
    arma::mat mat(svisits,1);mat.fill(1);
    arma::rowvec subX1=augX.row(i);
    arma::mat subX=kron(mat,subX1);
    arma::mat subW=mykron(time0, s);
    
    XtX=XtX+subX.t()*subX;
    XtW=XtW+subX.t()*subW;
    WtX=WtX+subW.t()*subX;
    WtW=WtW+subW.t()*subW;
    
    arma::mat subresid=y.rows(start,end);
    arma::mat subresidvec=vectorise(subresid.t());
    
    Xty=Xty+subX.t()*subresidvec;
    Wty=Wty+subW.t()*subresidvec;
  }
  arma::mat XtXinv=solveC(XtX);
  arma::mat out22=solveC2(WtW-WtX*XtXinv*XtW);
  arma::mat out11=XtXinv+XtXinv*XtW*out22*WtX*XtXinv;
  arma::mat out12=-(XtXinv*XtW*out22);
  arma::mat out21=-(out22*WtX*XtXinv);
  
  arma::vec beta1=out11*Xty+out12*Wty;
  arma::vec beta2=out21*Xty+out22*Wty;
  arma::vec beta=join_cols(beta1,beta2);
  return beta;
}

// [[Rcpp::export]]
arma::vec splocst_modelpermC(arma::mat resid, arma::vec nvisits, arma::vec time,arma::mat distMat,
                            arma::mat nnList,arma::vec varcomps, arma::vec dxStatus, arma::mat permMat, double alpha){
  int nperm=permMat.n_cols;
  int N0=permMat.n_rows;
  double phi=varcomps(0);
  double psi=varcomps(1);
  double sigma_s=varcomps(2);
  double sigma_e=varcomps(3);
  int start;
  int end;
  int s=resid.n_cols;
  int nsubj=nvisits.size();
  int nnm=nnList.n_cols;

  arma::mat ubMat(s, nsubj);
  resid=resid.t();

  arma::mat SpMat=generateSmat(phi, distMat);

  for (int i=0; i<nsubj; ++i){
    int visits=nvisits(i);
    int svisits=s*visits;
    if (i==0){
      start=0; end=visits-1;
    } else{
      start=sum(nvisits.subvec(0, i-1)); end=sum(nvisits.subvec(0,i))-1;
    }

    arma::vec time0=time.subvec(start, end);
    arma::mat mat(svisits,1);mat.fill(1);
    arma::mat subW=mykron(time0, s);

    arma::mat subresidvec=vectorise(resid.cols(start,end));
    arma::mat Vinv=solveCstable(sigma_s*kron(generateTmat(time0,psi), SpMat)+sigma_e*diagC(svisits));

    ubMat.col(i)=subW.t()*Vinv*subresidvec;
  }
  ubMat=ubMat.t();

  arma::vec pvec(nnm); pvec.fill(0);
  arma::vec pPerm(nperm); pPerm.fill(0);
  arma::rowvec Tstat;
  arma::mat out(nperm,s); out.fill(0);
  arma::mat Tstatmat(nperm, nnm); Tstatmat.fill(0);
  arma::rowvec U(s); U.fill(0);
  arma::rowvec Us(s);
  arma::vec permStatus;

  arma::rowvec sumU(s); sumU.fill(0);
  for (int i=0; i<nsubj; ++i){
    U=U+dxStatus(i)*ubMat.row(i);
    sumU=sumU+ubMat.row(i);
  }

  for (int perm=0; perm<nperm; ++perm){
    permStatus=permMat.col(perm);
    Us=sumU;
    for (int i=0; i<N0; ++i){
      Us=Us-2*ubMat.row(permStatus(i));
    }
    out.row(perm)=Us;
  }

  Tstat=U*nnList;
  Tstatmat=out*nnList;

  for (int j=0; j<nnm; ++j){
    double std=stddev(Tstatmat.col(j));
    Tstatmat.col(j)=Tstatmat.col(j)/std;
    Tstat(j)=Tstat(j)/std;
  }

  Tstat=pow(Tstat,2);
  Tstatmat=pow(Tstatmat,2);

  arma::vec Tstatvec(nperm); Tstatvec.fill(0);
  double Tstatmax=Tstat.max();
  for (int perm=0; perm<nperm; ++perm){
    Tstatvec(perm)=Tstatmat.row(perm).max();
  }

  double pval=meanprop(Tstatvec, Tstatmax);
  double tstat=quantileC(Tstatvec, alpha);
  arma::vec out2(nnm+2);
  for (int j=0; j<(nnm+2); ++j){
    if (j==0){ out2(j)=pval; }
    else if (j==1){ out2(j)=tstat;}
    else {
      out2(j)=Tstat(j-2);
    }
  }
  
  return out2;
}

// [[Rcpp::export]]
arma::vec covRegC(arma::mat residT, arma::vec nvisits, arma::vec time,arma::mat distmat,
                        double phi, double psi){
  int start;
  int end;
  int s=residT.n_rows;
  int nsubj=nvisits.size();
  arma::vec out(2);
  arma::mat XtX(2,2); XtX.fill(0);
  arma::vec Xty(2); Xty.fill(0);
  arma::mat spMat=generateSmat(phi, distmat);

  for (int i=0; i<nsubj; ++i){
    int visits=nvisits(i);
    int svisits=s*visits;
    if (i==0){
      start=0;end=visits-1;
    } else{
      start=sum(nvisits.subvec(0, i-1));end=sum(nvisits.subvec(0,i))-1;
    }

    arma::vec time0=time.subvec(start, end);
    arma::vec subresid1=vectorise(residT.cols(start,end));
    arma::vec subresidvec=vectorise(subresid1*subresid1.t());

    arma::mat STmat=kron(generateTmat(time0,psi),spMat);

    arma::vec cormat=vectorise(STmat);
    arma::vec emat=vectorise(diagC(svisits));

    XtX(0,0)=XtX(0,0)+sum(cormat%cormat);
    XtX(0,1)=XtX(0,1)+sum(cormat%emat);
    XtX(1,1)=XtX(1,1)+sum(emat);
    XtX(1,0)=XtX(0,1);

    Xty(0)=Xty(0)+sum(cormat%subresidvec);
    Xty(1)=Xty(1)+sum(emat%subresidvec);
  }

  out=solveC(XtX)*Xty;
  return out;
}

// [[Rcpp::export]]
double variomatSpatial(double phi, arma::vec varcomps,arma::mat residT, arma::mat distMat, arma::vec nvisits, arma::vec time){
  double psi=varcomps(0);
  double sigma2=varcomps(1);
  double tau2=varcomps(2);
  arma::mat Sigma1=generateSmat(phi, distMat);
  arma::mat Sigma;
  arma::vec subresidvec;
  int n=distMat.n_cols;
  int nsubj=nvisits.size();
  int visits;
  int start;
  int end;
  double out=0;
  for (int i=0; i<nsubj; ++i){
    visits=nvisits(i);
    if (i==0){
      start=0;
      end=visits-1;
    } else{
      start=sum(nvisits.subvec(0, i-1));
      end=sum(nvisits.subvec(0,i))-1;
    }
    arma::vec time0=time.subvec(start,end);
    
    subresidvec=vectorise(residT.cols(start,end));
    Sigma=sigma2*kron(generateTmat(time0, psi),Sigma1)+tau2*diagC(n*visits);
    out=out+frobC(subresidvec*subresidvec.t()-Sigma);
  }
  return out;
}

// [[Rcpp::export]]
double variomatTemporal(double psi, arma::vec varcomps,arma::mat residT, arma::mat distMat, arma::vec nvisits, arma::vec time){
  double phi=varcomps(0);
  double sigma2=varcomps(1);
  double tau2=varcomps(2);
  arma::mat Sigma1=generateSmat(phi, distMat);
  arma::mat Sigma;
  arma::vec subresidvec;
  int n=distMat.n_cols;
  int nsubj=nvisits.size();
  int visits;
  int start;
  int end;
  double out=0;
  for (int i=0; i<nsubj; ++i){
    visits=nvisits(i);
    if (i==0){
      start=0;end=visits-1;
    } else{
      start=sum(nvisits.subvec(0, i-1));end=sum(nvisits.subvec(0,i))-1;
    }
    arma::vec time0=time.subvec(start,end);
    subresidvec=vectorise(residT.cols(start,end));
    Sigma=sigma2*kron(generateTmat(time0, psi),Sigma1)+tau2*diagC(n*visits);
    out=out+frobC(subresidvec*subresidvec.t()-Sigma);
  }
  return out;
}

// [[Rcpp::export]]
arma::vec sploct_modelpermC(arma::mat resid, arma::vec nvisits, arma::vec time,arma::mat distMat,
                           arma::mat nnList,arma::vec varcomps, arma::vec dxStatus, arma::mat permMat, double alpha){
  int nperm=permMat.n_cols;
  int N0=permMat.n_rows;
  double phi=varcomps(1);
  int start;
  int end;
  int s=resid.n_cols;
  int nsubj=nvisits.size();
  int nnm=nnList.n_cols;
  
  arma::mat ubMat(s, nsubj);
  resid=resid.t();
  
  for (int i=0; i<nsubj; ++i){
    int visits=nvisits(i);
    int svisits=s*visits;
    if (i==0){
      start=0; end=visits-1;
    } else{
      start=sum(nvisits.subvec(0, i-1)); end=sum(nvisits.subvec(0,i))-1;
    }
    
    arma::vec time0=time.subvec(start, end);
    arma::mat V=generateTmat(time0,phi);
    
    arma::mat mat(svisits,1);mat.fill(1);
    arma::mat subW=mykron(time0, s);
    
    arma::mat subresidvec=vectorise(resid.cols(start,end));
    arma::mat Vinv=kron(solveC(V),diagC(s));
    
    ubMat.col(i)=subW.t()*Vinv*subresidvec;
  }
  ubMat=ubMat.t();
  
  double pval;
  arma::vec pvec(nnm); pvec.fill(0);
  arma::vec pPerm(nperm); pPerm.fill(0);
  arma::rowvec Tstat;
  arma::mat out(nperm,s); out.fill(0);
  arma::mat Tstatmat(nperm, nnm); Tstatmat.fill(0);
  arma::rowvec U(s); U.fill(0);
  arma::rowvec Us(s);
  arma::vec permStatus;
  
  arma::rowvec sumU(s); sumU.fill(0);
  for (int i=0; i<nsubj; ++i){
    U=U+dxStatus(i)*ubMat.row(i);
    sumU=sumU+ubMat.row(i);
  }
  
  for (int perm=0; perm<nperm; ++perm){
    permStatus=permMat.col(perm);
    Us=sumU;
    for (int i=0; i<N0; ++i){
      Us=Us-2*ubMat.row(permStatus(i));
    }
    out.row(perm)=Us;
  }
  
  Tstat=U*nnList;
  Tstatmat=out*nnList;
  
  for (int j=0; j<nnm; ++j){
    double std=stddev(Tstatmat.col(j));
    Tstatmat.col(j)=Tstatmat.col(j)/std;
    Tstat(j)=Tstat(j)/std;
  }
  
  Tstat=pow(Tstat,2);
  Tstatmat=pow(Tstatmat,2);
  
  arma::vec Tstatvec(nperm); Tstatvec.fill(0);
  double Tstatmax=Tstat.max();
  for (int perm=0; perm<nperm; ++perm){
    Tstatvec(perm)=Tstatmat.row(perm).max();
  }
  
  pval=meanprop(Tstatvec, Tstatmax);
  double tstat=quantileC(Tstatvec, alpha);
  arma::vec out2(nnm+2);
  for (int j=0; j<(nnm+2); ++j){
    if (j==0){ out2(j)=pval; }
    else if (j==1){ out2(j)=tstat;}
    else {
      out2(j)=Tstat(j-2);
    }
  }
  return out2;
}

