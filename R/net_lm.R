

###############################
#####  Linear Regression  #####
###############################

LmL0=function(x, y, Omega=NULL, alpha=1.0, lambda=NULL, nlambda=50, rlambda=NULL, wbeta=rep(1,ncol(x)), sgn=rep(1,ncol(x)), nfolds=1, foldid=NULL, ill=TRUE, iL0=TRUE, icutB=TRUE, ncutB=10, ifast=TRUE, isd=FALSE, iysd=FALSE, keep.beta=FALSE, thresh=1e-6, maxit=1e+3) {
  
  N0=nrow(x); p=ncol(x)
  
  sdy=1.0
  if (iysd) {
    sdy=sd(y)
    y=y/sdy
  }
  
  ### Adaptive
  aPen=ifelse(all(wbeta>0), TRUE, FALSE)
  wbeta2=ifelse(wbeta==0,0,1)
  
  ### Lambda path
  if (all(wbeta==0)) {
    lambda=0.0
  }
  
  if (is.null(lambda)) {
    ilambda=1
    if (is.null(rlambda)) {
      rlambda=ifelse(N0>p, 0.0001, 0.01)
    }
    lambda=(rlambda)^(c(0:(nlambda-1))/(nlambda-1))
  } else {
    ilambda=0
    nlambda=length(lambda)
  }
  
  if (is.null(Omega) | p==1) {
    penalty=ifelse(alpha==1, "Lasso", "Enet")
    adaptive=ifelse(any(wbeta!=1), TRUE, FALSE)
  } else {
    penalty=ifelse(alpha==1, "Lasso", "Net")
    adaptive=c(ifelse(any(wbeta!=1), TRUE, FALSE),ifelse(any(sgn!=1), TRUE, FALSE))
    
    Omega=rbind(0,cbind(0,Omega)); sgn1=c(1,sgn) # intercept
    ## Check Omega: positive off-diagonal, zero diagonal
    if (any(diag(Omega)!=0)) {
      diag(Omega)=0
      # cat("Diagonal of Omega was set to all zeros\n")
    }
    if (any(Omega<0)) {
      Omega=abs(Omega)
      # cat("Off-diagonal of Omega was foced to non-negative values\n")
    }
    
    ### Correlation/Adjacency matrix
    if (inherits(Omega, "dgCMatrix")) {
      W=OmegaSC(Omega, sgn1); W$loc=W$loc+1
    } else {
      W=OmegaC(Omega, sgn1); W$loc=W$loc+1
    }
    rm(Omega)
  }
  
  
  #####  Run  #####
  out=switch(penalty,
             "Net"=NetLmC(x, y, alpha, lambda, nlambda, ilambda, wbeta, wbeta2, W$Omega, W$loc, W$nadj, p, N0, thresh, maxit, 1e-5),
             EnetLmC(x, y, alpha, lambda, nlambda, ilambda, wbeta, wbeta2, p, N0, thresh, maxit, 1e-5)
  )
  nlambdai=out$nlambda ## number of lambdas
  if (nlambdai==0)
    return(NULL)
  lambdai=out$lambda[1:nlambdai]
  
  out$Beta[is.na(out$Beta)]=out$BetaSTD[is.na(out$Beta)]
  out$Beta=Matrix(out$Beta[, 1:nlambdai,drop=FALSE]*sdy, sparse=TRUE)
  out$BetaSTD=Matrix(out$BetaSTD[, 1:nlambdai,drop=FALSE]*sdy, sparse=TRUE)
  out$nzero=apply(out$Beta!=0, 2, sum)
  out$flag=out$flag[1:nlambdai]
  out$rsq=out$rsq[1:nlambdai]
  out$a=out$a[1:nlambdai]
  
  
  if (nfolds==1 & is.null(foldid)) {
    
    fit=data.frame(lambda=lambdai, rsq=out$rsq, nzero=out$nzero)
    if (!isd) {
      return(list(a=out$a, Beta=out$Beta, fit=fit, penalty=penalty, adaptive=adaptive, flag=out$flag))
    } else {
      return(list(a=out$a, Beta=out$BetaSTD, fit=fit, penalty=penalty, adaptive=adaptive, flag=out$flag))
    }
    
  } else {
    
    ########################################
    #####  Cross-validation estimates  #####
    
    ###  Split data for cross-validation
    if (is.null(foldid)) {
      foldid=sample(rep(seq(nfolds), length=N0))
    } else {
      nfolds=max(foldid)
    }
    tb=table(foldid)
    N0i=numeric(nfolds); Nf=numeric(nfolds)
    for (i in 1:nfolds) {
      N0i[i]=sum(tb[-i]); Nf[i]=tb[i]
    }
    weighti=as.vector(tapply(rep(1,N0), foldid, sum))
    
    
    outi=list(); cvRSS=matrix(NA, nrow=nfolds, ncol=nlambdai)
    for (i in 1:nfolds) {
      temid=(foldid==i)
      
      outi[[i]]=switch(penalty,
                       "Net"=cvNetLmC(x[!temid, ,drop=FALSE], y[!temid], alpha, lambdai, nlambdai, wbeta, wbeta2, W$Omega, W$loc, W$nadj, N0i[i],  p,  thresh, maxit, x[temid, ,drop=FALSE], y[temid], Nf[i], 1e-5),
                       cvEnetLmC(x[!temid, ,drop=FALSE], y[!temid], alpha, lambdai, nlambdai, wbeta, wbeta2, N0i[i],p, thresh, maxit, x[temid, ,drop=FALSE], y[temid], Nf[i], 1e-5)
      )
      
      outi[[i]]$Beta[is.na(outi[[i]]$Beta)]=outi[[i]]$BetaSTD[is.na(outi[[i]]$Beta)]
      
      
      if (ill) {
        delta2i=outi[[i]]$RSS[1:outi[[i]]$nlambda]/N0i[i]
        cvRSS[i, 1:outi[[i]]$nlambda]=outi[[i]]$RSSp[1:outi[[i]]$nlambda]/delta2i+Nf[i]*log(delta2i) ## for ith fold
      } else {
        cvRSS[i, 1:outi[[i]]$nlambda]=outi[[i]]$RSSp[1:outi[[i]]$nlambda] ## for ith fold
      }
      
    }
    
    cvRSS=cvRSS[, 1:nlambdai,drop=FALSE]
    cvraw=cvRSS/weighti; nfoldi=apply(!is.na(cvraw), 2, sum); #rm(cvRSS) #
    cvm=apply(cvraw, 2, weighted.mean, w=weighti, na.rm=TRUE)
    cvse=sqrt(apply(sweep(cvraw, 2, cvm, "-")^2, 2, weighted.mean, w=weighti, na.rm=TRUE)/(nfoldi-1))
    
    
    indexi=which.min(cvm)
    indexij=which(cvm<=(cvm[indexi]+cvse[indexi]))[1]
    temi=rep("", nlambdai)
    temi[indexi]="*";#temi[indexij]=ifelse(temi[indexij]=="", "*", "***")
    #temCV=data.frame(lambda=lambdai, cvm=cvm, cvse=cvse, nzero=out$nzero, index=temi,stringsAsFactors=FALSE)
    temCV=data.frame(lambda=lambdai, rsq=out$rsq, cvm=cvm, cvse=cvse, nzero=out$nzero, index=temi, stringsAsFactors=FALSE)
    
    
    if (!iL0) {
      if (!keep.beta) {
        if (!isd) {
          return(list(a=out$a[indexi], Beta=out$Beta[, indexi], fit=temCV, lambda.min=lambdai[indexi], penalty=penalty, adaptive=adaptive, flag=out$flag))
        } else {
          return(list(Beta=out$BetaSTD[, indexi], fit=temCV, lambda.min=lambdai[indexi], penalty=penalty, adaptive=adaptive, flag=out$flag))
        }
        
      } else {
        if (!isd) {
          return(list(a=out$a, Beta=out$Beta, fit=temCV, lambda.min=lambdai[indexi], penalty=penalty, adaptive=adaptive, flag=out$flag))
        } else {
          return(list(Beta=out$BetaSTD, fit=temCV, lambda.min=lambdai[indexi], penalty=penalty, adaptive=adaptive, flag=out$flag))
        }
      }
    }
    
    
    
    #####################################
    #####  Cross-validation for L0  #####
    if (icutB) {
      
      ###  cutoff  ###
      il0=ifelse(ifast, indexi, 1)
      cvm=list(); cv.min=rep(NA, nlambdai)
      repeat {
        
        BetaSTD=out$BetaSTD[,il0]
        cut0=sort(unique(c(0,abs(BetaSTD))),decreasing=FALSE); cuti=NULL
        for (i in 1:ncutB) {
          cuti=c(cuti, diff(cut0)/ncutB*i+cut0[-length(cut0)])
        }
        cut0=sort(unique(c(cut0,cuti)))
        
        Betai=matrix(sapply(outi, function(x){x$Beta[, il0,drop=FALSE]}), nrow=p)
        BetaSTDi=matrix(sapply(outi, function(x){x$BetaSTD[, il0,drop=FALSE]}), nrow=p)
        
        cvRSS=matrix(NA, nrow=nfolds, ncol=length(cut0)); i=1
        for (i in 1:nfolds) {
          temid=foldid==i
          Betaj=Betai[, i]; BetaSTDj=BetaSTDi[, i]
          
          outij=cvHardLmC(Betaj, BetaSTDj, cut0, wbeta, x[!temid,,drop=FALSE],y[!temid],N0i[i],p, x[temid,,drop=FALSE],y[temid],Nf[i])
          
          if (ill) {
            delta2i=outij$RSS/N0i[i]
            cvRSS[i,]=outij$RSSp/delta2i+Nf[i]*log(delta2i)
          } else {
            cvRSS[i,]=outij$RSSp
          }
        }
        
        cvraw=cvRSS/weighti; nfoldi=apply(!is.na(cvraw), 2, sum); #rm(cvRSS) #
        cvm[[il0]]=apply(cvraw, 2, weighted.mean, w=weighti, na.rm=TRUE)
        cv.min[il0]=min(cvm[[il0]])
        
        
        il1=which.min(cv.min)
        if (nlambdai==1) break
        if (ifast) {
          if (sum(!is.na(cv.min))>1) {
            if (!is.na(cv.min[pmin(il1+1,nlambdai)]) & !is.na(cv.min[pmax(il1-1,1)])) {
              break
            } else if (is.na(cv.min[pmin(il1+1,nlambdai)])) {
              il0=il1+1
            } else if (is.na(cv.min[pmax(il1-1,1)])) {
              il0=il1-1
            } else {
              il0=which.min(cv.min)
            }
            
          } else if (il0+1<=nlambdai) {
            il0=il0+1
          } else if (il0-1>=1) {
            il0=il0-1
          }
        } else {
          il0=il0+1
          if (il0>nlambdai) break
        }
      }
      
      il0=which.min(cv.min)
      nzero0=which.min(cvm[[il0]])
      
      Beta0=out$Beta[,il0]
      BetaSTD0=out$BetaSTD[,il0]
      
      
      BetaSTD=out$BetaSTD[,il0]
      cut0=sort(unique(c(0,abs(BetaSTD))),decreasing=FALSE); cuti=NULL
      for (i in 1:ncutB) {
        cuti=c(cuti, diff(cut0)/ncutB*i+cut0[-length(cut0)])
      }
      cut0=sort(unique(c(cut0,cuti)))
      
      
      Beta0[abs(BetaSTD0)<=cut0[nzero0] & wbeta>0.0]=0.0
      BetaSTD0[abs(BetaSTD0)<=cut0[nzero0] & wbeta>0.0]=0.0
      
      temCV0=data.frame(lambda=lambdai[il0],cvm=cv.min[il0],nzero=sum(Beta0!=0))
      
      
    } else {
      
      
      ###  number of non-zeros  ###
      il0=ifelse(ifast, indexi, 1)
      cvm=list(); cv.min=rep(NA, nlambdai)
      repeat {
        
        numi=out$nzero[il0]
        Betai=matrix(sapply(outi, function(x){x$Beta[, il0,drop=FALSE]}), nrow=p)
        BetaSTDi=matrix(sapply(outi, function(x){x$BetaSTD[, il0,drop=FALSE]}), nrow=p)
        
        Betao=apply(Betai!=0, 2, sum)
        numi2=pmax(min(max(Betao), numi),1)
        
        cvRSS=matrix(NA, nrow=nfolds, ncol=numi2)
        for (i in 1:nfolds) {
          numj=min(Betao[i], numi); temid=foldid==i
          Betaj=Betai[, i]; BetaSTDj=BetaSTDi[, i]
          
          if (numj==0) {
            outij=cvTrimLmC(c(0.0, 0.0), numj, numi2, c(0, 0), x[!temid,,drop=FALSE],y[!temid],N0i[i],p, x[temid, ,drop=FALSE], y[temid], Nf[i])
            
            if (ill) {
              delta2i=outij$RSS/N0i[i]
              cvRSS[i,]=outij$RSSp/delta2i+Nf[i]*log(delta2i)
            } else {
              cvRSS[i,]=outij$RSSp
            }
            
          } else {
            
            BetaSTDjj=BetaSTDj
            BetaSTDjj[wbeta==0]=max(abs(BetaSTDj))+1
            temo=rank(-abs(BetaSTDjj), ties.method="min")
            
            temo=data.frame(o=temo[which(temo<=numj)], loc=which(temo<=numj))
            temo=temo[order(temo$o), ]
            temo=temo[1:numj,]
            
            outij=cvTrimLmC(Betaj[temo$loc], numj, numi2, temo$loc-1, x[!temid,,drop=FALSE],y[!temid],N0i[i],p, x[temid, ,drop=FALSE], y[temid], Nf[i])
            
            
            if (ill) {
              delta2i=outij$RSS/N0i[i]
              cvRSS[i,]=outij$RSSp/delta2i+Nf[i]*log(delta2i)
            } else {
              cvRSS[i,]=outij$RSSp
            }
          }
        }
        
        cvraw=cvRSS/weighti; nfoldi=apply(!is.na(cvraw), 2, sum); #rm(cvRSS) #
        cvm[[il0]]=apply(cvraw, 2, weighted.mean, w=weighti, na.rm=TRUE)
        temi=cvm[[il0]]
        
        if (aPen) {
          cv.min[il0]=min(temi)
        } else {
          cv.min[il0]=min(temi[sum(wbeta==0):length(temi)])
        }
        
        
        il1=which.min(cv.min)
        if (nlambdai==1) break
        if (ifast) {
          if (sum(!is.na(cv.min))>1) {
            if (!is.na(cv.min[pmin(il1+1,nlambdai)]) & !is.na(cv.min[pmax(il1-1,1)])) {
              break
            } else if (is.na(cv.min[pmin(il1+1,nlambdai)])) {
              il0=il1+1
            } else if (is.na(cv.min[pmax(il1-1,1)])) {
              il0=il1-1
            } else {
              il0=which.min(cv.min)
            }
            
          } else if (il0+1<=nlambdai) {
            il0=il0+1
          } else if (il0-1>=1) {
            il0=il0-1
          }
        } else {
          il0=il0+1
          if (il0>nlambdai) break
        }
      }
      
      il0=which.min(cv.min)
      
      Beta0=out$Beta[,il0]
      BetaSTD0=out$BetaSTD[,il0]
      
      temi=cvm[[il0]]
      if (aPen) {
        cuti=which.min(temi)
      } else {
        cuti=which.min(temi[sum(wbeta==0):length(temi)])+sum(wbeta==0)-1
      }
      
      Beta0j=out$BetaSTD[,il0]
      Beta0j[which(wbeta==0)]=max(abs(Beta0j))+1
      
      if (cuti<length(Beta0j)) {
        Beta0[abs(Beta0j)<=sort(abs(Beta0j),TRUE)[cuti+1]]=0
        BetaSTD0[abs(Beta0j)<=sort(abs(Beta0j),TRUE)[cuti+1]]=0
      }
      
      temCV0=data.frame(lambda=lambdai[il0],cvm=cv.min[il0],nzero=sum(Beta0!=0))
    }
    
    
    if (!keep.beta) {
      
      if (!isd) {
        return(list(a=out$a[indexi], Beta=out$Beta[, indexi], Beta0=Beta0, fit=temCV, fit0=temCV0, lambda.min=lambdai[indexi], lambda.opt=lambdai[il0], penalty=penalty, adaptive=adaptive, flag=out$flag))
      } else {
        return(list(Beta=out$BetaSTD[, indexi], Beta0=BetaSTD0, fit=temCV, fit0=temCV0, lambda.min=lambdai[indexi], lambda.opt=lambdai[il0], penalty=penalty, adaptive=adaptive, flag=out$flag))
      }
      
    } else {
      if (!isd) {
        return(list(a=out$a, Beta=out$Beta, Beta0=Beta0, fit=temCV, fit0=temCV0, lambda.min=lambdai[indexi], lambda.opt=lambdai[il0], penalty=penalty, adaptive=adaptive, flag=out$flag))
      } else {
        return(list(Beta=out$BetaSTD, Beta0=BetaSTD0, fit=temCV, fit0=temCV0, lambda.min=lambdai[indexi], lambda.opt=lambdai[il0], penalty=penalty, adaptive=adaptive, flag=out$flag))
      }
    }
    
  }
}




