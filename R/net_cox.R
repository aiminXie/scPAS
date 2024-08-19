

#######################
#####  Cox model  #####
#######################

CoxL0=function(x, y, Omega=NULL, alpha=1.0, lambda=NULL, nlambda=100, rlambda=NULL, wbeta=rep(1,ncol(x)), sgn=rep(1,ncol(x)), nfolds=1, foldid=NULL, iL0=TRUE, icutB=TRUE, ncutB=10, ifast=TRUE, isd=FALSE, ifastr=TRUE, keep.beta=FALSE, thresh=1e-6, maxit=1e+5) {
  
  N0=nrow(x); p=ncol(x)
  ifastr=as.integer(ifastr)
  
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
  
  
  if (is.null(Omega)) {
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
  prep0=PrepCox(x, y)
  out=switch(penalty,
             "Net"=NetCoxC(prep0$x, prep0$tevent, alpha, lambda, nlambda, ilambda, wbeta, wbeta2, W$Omega, W$loc, W$nadj, prep0$N, prep0$nevent, prep0$nevent1, prep0$loc1, prep0$n, p, N0, thresh, maxit, ifastr),
             EnetCoxC(prep0$x, prep0$tevent, alpha, lambda, nlambda, ilambda, wbeta, wbeta2, prep0$N, prep0$nevent, prep0$nevent1, prep0$loc1, prep0$n, p, N0, thresh, maxit, ifastr)
  )
  nlambdai=out$nlambda
  if (nlambdai==0)
    return(NULL)
  lambdai=out$lambda[1:nlambdai]
  
  out$Beta[is.na(out$Beta)]=out$BetaSTD[is.na(out$Beta)]
  out$Beta=Matrix(out$Beta[, 1:nlambdai,drop=F], sparse=TRUE)
  out$BetaSTD=Matrix(out$BetaSTD[, 1:nlambdai], sparse=TRUE)
  out$nzero=apply(out$Beta!=0, 2, sum)
  out$flag=out$flag[1:nlambdai]
  
  
  if (nfolds==1 & is.null(foldid)) {
    
    fit=data.frame(lambda=lambdai, nzero=out$nzero)
    if (!isd) {
      return(list(Beta=out$Beta, fit=fit, penalty=penalty, adaptive=adaptive, flag=out$flag))
    } else {
      return(list(Beta=out$BetaSTD, fit=fit, penalty=penalty, adaptive=adaptive, flag=out$flag))
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
    tb=table(foldid);N0i=numeric(nfolds)
    for (i in 1:nfolds)
      N0i[i]=sum(tb[-i])
    
    prepk=list()
    for (i in 1:nfolds) {
      temid=which(foldid!=i)
      prepk[[i]]=PrepCox(x[temid, ], y[temid, ])
    }
    weighti=as.vector(tapply(y[, "status"], foldid, sum))
    
    
    outi=list(); cvPL=matrix(NA, nrow=nfolds, ncol=nlambdai)
    for (i in 1:nfolds) {
      outi[[i]]=switch(penalty,
                       "Net"=cvNetCoxC(prepk[[i]]$x, prepk[[i]]$tevent, alpha, lambdai, nlambdai, wbeta, wbeta2, W$Omega, W$loc, W$nadj, prepk[[i]]$N, prepk[[i]]$nevent, prepk[[i]]$nevent1, prepk[[i]]$loc1, prepk[[i]]$n, p, N0i[i], thresh, maxit, 0, prep0$x, prep0$N, prep0$nevent, prep0$nevent1, prep0$loc1, prep0$n),
                       cvEnetCoxC(prepk[[i]]$x, prepk[[i]]$tevent, alpha, lambdai, nlambdai, wbeta, wbeta2, prepk[[i]]$N, prepk[[i]]$nevent, prepk[[i]]$nevent1, prepk[[i]]$loc1, prepk[[i]]$n, p, N0i[i], thresh, maxit, 0, prep0$x, prep0$N, prep0$nevent, prep0$nevent1, prep0$loc1, prep0$n)
      )
      
      outi[[i]]$Beta[is.na(outi[[i]]$Beta)]=outi[[i]]$BetaSTD[is.na(outi[[i]]$Beta)]
      
      cvPL[i, 1:outi[[i]]$nlambda]=outi[[i]]$lf[1:outi[[i]]$nlambda]-outi[[i]]$ll[1:outi[[i]]$nlambda]
    }
    
    temi=apply(abs(cvPL)==Inf,2,sum,na.rm=TRUE)
    if (any(temi>0)) {
      nlambdai=max(min(which(temi>0))-1,1)
    }
    
    cvraw=cvPL/weighti; nfoldi=apply(!is.na(cvraw), 2, sum); #rm(cvPL) #
    cvm=apply(cvraw, 2, weighted.mean, w=weighti, na.rm=TRUE)
    cvse=sqrt(apply(sweep(cvraw, 2, cvm, "-")^2, 2, weighted.mean, w=weighti, na.rm=TRUE)/(nfoldi-1))
    
    cvraw=cvraw[1:nlambdai]
    cvm=cvm[1:nlambdai]
    cvse=cvse[1:nlambdai]
    
    indexi=which.max(cvm)
    indexij=which(cvm>=(cvm[indexi]-cvse[indexi]))[1]
    temi=rep("", nlambdai)
    temi[indexi]="*"# ;temi[indexij]=ifelse(temi[indexij]=="", "*", "***")
    #temCV=data.frame(lambda=lambdai, cvm=cvm, cvse=cvse, nzero=out$nzero, index=temi,stringsAsFactors=FALSE)
    temCV=data.frame(lambda=lambdai[1:nlambdai], cvm=cvm, cvse=cvse, nzero=out$nzero[1:nlambdai], index=temi[1:nlambdai], stringsAsFactors=FALSE)
    temCV$cvm=-temCV$cvm # compatible with LM
    
    if (!iL0) {
      if (!keep.beta) {
        if (!isd) {
          return(list(Beta=out$Beta[, indexi], fit=temCV, lambda.min=lambdai[indexi], penalty=penalty, adaptive=adaptive, flag=out$flag[1:nlambdai]))
        } else {
          return(list(Beta=out$BetaSTD[, indexi], fit=temCV, lambda.min=lambdai[indexi], penalty=penalty, adaptive=adaptive, flag=out$flag[1:nlambdai]))
        }

      } else {
        if (!isd) {
          return(list(Beta=out$Beta[,1:nlambdai], fit=temCV, lambda.min=lambdai[indexi], penalty=penalty, adaptive=adaptive, flag=out$flag[1:nlambdai]))
        } else {
          return(list(Beta=out$BetaSTD[,1:nlambdai], fit=temCV, lambda.min=lambdai[indexi], penalty=penalty, adaptive=adaptive, flag=out$flag[1:nlambdai]))
        }
      }
    }
    
    
    
    #####################################
    #####  Cross-validation for L0  #####
    if (icutB) {
      
      ###  cutoff  ###
      il0=ifelse(ifast, indexi, 1)
      cvm=list(); cv.max=rep(NA, nlambdai)
      repeat {
        
        BetaSTD=out$BetaSTD[,il0]
        cut0=sort(unique(c(0,abs(BetaSTD))),decreasing=FALSE); cuti=NULL
        for (i in 1:ncutB) {
          cuti=c(cuti, diff(cut0)/ncutB*i+cut0[-length(cut0)])
        }
        cut0=sort(unique(c(cut0,cuti)))
        
        Betai=matrix(sapply(outi, function(x){x$Beta[, il0,drop=FALSE]}), nrow=p)
        BetaSTDi=matrix(sapply(outi, function(x){x$BetaSTD[, il0,drop=FALSE]}), nrow=p)
        
        cvPL=matrix(NA, nrow=nfolds, ncol=length(cut0)); i=1
        for (i in 1:nfolds) {
          Betaj=Betai[, i]; BetaSTDj=BetaSTDi[, i]
          
          cvPL[i, ]=cvHardCoxC(Betaj, BetaSTDj, cut0, wbeta, p, prep0$x, prep0$N, prep0$nevent, prep0$nevent1, prep0$loc1, prep0$n, prepk[[i]]$x, prepk[[i]]$N, prepk[[i]]$nevent, prepk[[i]]$nevent1, prepk[[i]]$loc1, prepk[[i]]$n, 0, 1)
        }
       
        
        cvraw=cvPL/weighti;nfoldi=apply(!is.na(cvraw), 2, sum); #rm(cvPL) #
        cvm[[il0]]=apply(cvraw, 2, weighted.mean, w=weighti, na.rm=TRUE)
        temi=cvm[[il0]]
        
        if (aPen) {
          cv.max[il0]=max(temi)
        } else {
          cv.max[il0]=max(temi[sum(wbeta==0):length(temi)])
        }
        
        
        il1=which.max(cv.max)
        if (nlambdai==1) break
        if (ifast) {
          if (sum(!is.na(cv.max))>1) {
            if (!is.na(cv.max[pmin(il1+1,nlambdai)]) & !is.na(cv.max[pmax(il1-1,1)])) {
              break
            } else if (is.na(cv.max[pmin(il1+1,nlambdai)])) {
              il0=il1+1
            } else if (is.na(cv.max[pmax(il1-1,1)])) {
              il0=il1-1
            } else {
              il0=which.max(cv.max)
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
      
      
      il0=which.max(cv.max)
      nzero0=which.max(cvm[[il0]])

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

      temCV0=data.frame(lambda=lambdai[il0],cvm=cv.max[il0],nzero=sum(Beta0!=0))
      temCV0$cvm=-temCV0$cvm # compatible with LM
      
      
    } else {
      
      
      ###  number of non-zeros  ###
      il0=ifelse(ifast, indexi, 1)
      cvm=list(); cv.max=rep(NA, nlambdai)
      repeat {
        
        numi=out$nzero[il0]
        Betai=sapply(outi, function(x){x$Beta[, il0, drop=FALSE]})
        BetaSTDi=sapply(outi, function(x){x$BetaSTD[, il0, drop=FALSE]})
        
        Betao=apply(Betai!=0, 2, sum)
        numi2=pmax(min(max(Betao), numi),1)
        
        cvPL=matrix(NA, nrow=nfolds, ncol=numi2)
        for (i in 1:nfolds) {
          numj=min(Betao[i], numi)
          Betaj=Betai[, i]; BetaSTDj=BetaSTDi[, i]
          
          if (numj==0) {
            
            cvPL[i, ]=cvTrimCoxC(c(0.0, 0.0), numj, numi2, c(0, 0), prep0$x, prep0$N, prep0$nevent, prep0$nevent1, prep0$loc1, prep0$n, prepk[[i]]$x, prepk[[i]]$N, prepk[[i]]$nevent, prepk[[i]]$nevent1, prepk[[i]]$loc1, prepk[[i]]$n, 0, 1)
            
          } else {
            
            BetaSTDjj=BetaSTDj
            BetaSTDjj[wbeta==0]=max(abs(BetaSTDj))+1
            temo=rank(-abs(BetaSTDjj), ties.method="min")
            
            temo=data.frame(o=temo[which(temo<=numj)], loc=which(temo<=numj))
            temo=temo[order(temo$o), ]
            temo=temo[1:numj,]
            
            cvPL[i, ]=cvTrimCoxC(Betaj[temo$loc], numj, numi2, temo$loc-1, prep0$x, prep0$N, prep0$nevent, prep0$nevent1, prep0$loc1, prep0$n, prepk[[i]]$x, prepk[[i]]$N, prepk[[i]]$nevent, prepk[[i]]$nevent1, prepk[[i]]$loc1, prepk[[i]]$n, 0, 1)
          }
          
        }
        
        cvraw=cvPL/weighti;nfoldi=apply(!is.na(cvraw), 2, sum); #rm(cvPL) #
        cvm[[il0]]=apply(cvraw, 2, weighted.mean, w=weighti, na.rm=TRUE)
        temi=cvm[[il0]]
        
        if (aPen) {
          cv.max[il0]=max(temi)
        } else {
          cv.max[il0]=max(temi[sum(wbeta==0):length(temi)])
        }
        
        
        il1=which.max(cv.max)
        if (nlambdai==1) break
        if (ifast) {
          if (sum(!is.na(cv.max))>1) {
            if (!is.na(cv.max[pmin(il1+1,nlambdai)]) & !is.na(cv.max[pmax(il1-1,1)])) {
              break
            } else if (is.na(cv.max[pmin(il1+1,nlambdai)])) {
              il0=il1+1
            } else if (is.na(cv.max[pmax(il1-1,1)])) {
              il0=il1-1
            } else {
              il0=which.max(cv.max)
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
      
      
      il0=which.max(cv.max)
      
      Beta0=out$Beta[,il0]
      BetaSTD0=out$BetaSTD[,il0]
      
      temi=cvm[[il0]]
      if (aPen) {
        cuti=which.max(temi)
      } else {
        cuti=which.max(temi[sum(wbeta==0):length(temi)])+sum(wbeta==0)-1
      }
      
      Beta0j=out$BetaSTD[,il0]
      Beta0j[which(wbeta==0)]=max(abs(Beta0j))+1
      
      if (cuti<length(Beta0j)) {
        Beta0[abs(Beta0j)<=sort(abs(Beta0j),TRUE)[cuti+1]]=0
        BetaSTD0[abs(Beta0j)<=sort(abs(Beta0j),TRUE)[cuti+1]]=0
      }
      
      temCV0=data.frame(lambda=lambdai[il0],cvm=cv.max[il0],nzero=sum(Beta0!=0))
      temCV0$cvm=-temCV0$cvm # compatible with LM
      
    }
    
    
    if (!keep.beta) {
      
      if (!isd) {
        return(list(Beta=out$Beta[, indexi], Beta0=Beta0, fit=temCV, fit0=temCV0, lambda.min=lambdai[indexi], lambda.opt=lambdai[il0], penalty=penalty, adaptive=adaptive, flag=out$flag[1:nlambdai]))
      } else {
        return(list(Beta=out$BetaSTD[, indexi], Beta0=BetaSTD0, fit=temCV, fit0=temCV0, lambda.min=lambdai[indexi], lambda.opt=lambdai[il0], penalty=penalty, adaptive=adaptive, flag=out$flag[1:nlambdai]))
      }
      
    } else {
      
      if (!isd) {
        return(list(Beta=out$Beta[,1:nlambdai], Beta0=Beta0, fit=temCV, fit0=temCV0, lambda.min=lambdai[indexi], lambda.opt=lambdai[il0], penalty=penalty, adaptive=adaptive, flag=out$flag[1:nlambdai]))
      } else {
        return(list(Beta=out$BetaSTD[,1:nlambdai], Beta0=BetaSTD0, fit=temCV, fit0=temCV0, lambda.min=lambdai[indexi], lambda.opt=lambdai[il0], penalty=penalty, adaptive=adaptive, flag=out$flag[1:nlambdai]))
      }
    }
    
    
  }
}





##################################################
#####  Cox: Prepare data for log-likelihood  #####
##################################################

PrepCox=function(x, y){
  
  N0=nrow(x)
  oi=order(y[, "status"], decreasing=TRUE)
  x=x[oi, ];y=y[oi, ]
  oi=order(y[, "time"])
  x=x[oi, ];y=y[oi, ]
  
  ## remove the first censored cases
  i1=which(y[, "status"]==1);mi1=min(i1)-1
  if (mi1!=0) {
    x=x[-c(1:mi1), ];y=y[-c(1:mi1), ]
  }
  ty=y[, "time"];tevent=y[, "status"]
  N=nrow(x);n1=sum(y[, "status"])
  
  dty=duplicated(ty) # ties
  
  ### for calculation of log-likelihood
  if (any(dty)) {
    tevent0=tevent
    tevent0[which(dty)]=0
    
    ievent=cumsum(tevent0);loc1=which(tevent0==1)
    nevent=table(ievent);n=length(unique(ievent))
    nevent1=tapply(tevent==1, ievent, sum)
  } else {
    ievent=cumsum(tevent);loc1=which(tevent==1)
    nevent=table(ievent);n=length(unique(ievent))
    nevent1=rep(1, n)
  }
  
  return(list(x=x, N0=N0, tevent=tevent, N=N, nevent=nevent, nevent1=nevent1, loc1=loc1, n=n))
}




