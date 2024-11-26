
###Bayesian Information Criteria


## aux functions

CovARp <- function(phi,ti) {
  p <- length(phi)
  n <- max(ti)
  if (n==1) Rn <- matrix(1)
  else Rn <- stats::toeplitz(ARMAacf(ar=phi, ma=0, lag.max = n-1))
  rhos <- stats::ARMAacf(ar=phi, ma=0, lag.max = p)[-1]
  Rn <- Rn/(1-sum(rhos*phi))
  return(as.matrix(Rn[ti,ti]))
}


CovDEC <- function(phi1, phi2, ti) {
  ni <- length(ti)
  Rn <- diag(ni)
  if (ni==1) Rn <- matrix(1)
  else {
    for (i in 1:(ni-1)) for (j in (i+1):ni) Rn[i,j] <- phi1^(abs(ti[i]-ti[j])^phi2)
    Rn[lower.tri(Rn)] <- t(Rn)[lower.tri(Rn)]
  }
  return(Rn)
}


#### Slash function

## find pdf of SLASH

dSL<-function(y,mu,sigma,nu){
  f<-function(u) nu*u^(nu-1)*mnormt::dmnorm(y,mu,u^(-1)*sigma)
  resp <- stats::integrate(Vectorize(f),0,1)$value
  return(resp)
}


## find CDF of SLASH

acumSL<-function(y,mu,sigma,nu){
  f<-function(u) mnormt::pmnorm(y,mu,u^(-1)*sigma)*nu*u^(nu-1)
  resp <- stats::integrate(Vectorize(f),0,1)$value
  return(resp)
}

## Conditional CDF of Slash


grafSL<-function(y0,mu0,sigma0,Qc,muc,sigmac,nu){
  f<-function(u) mnormt::dmnorm(y0,mu0,u^(-1)*sigma0)*pmnorm(Qc,muc,u^(-1)*sigmac)*nu*u^(nu-1)
  resp <- stats::integrate(Vectorize(f),0,1)$value
  return(resp)
}



##Likelihood functions for Student-t and Normal models

VerCensLMM <- function(cc, nj, y, x, z, beta1, sigmae, DD, tt, phi=NULL, nu=NULL, distr="Normal", cens.type="left", depstr = "UNC",LI=NULL,LS=NULL){

  GB = mvtnorm::GenzBretz(maxpts = 5e4, abseps = 1e-9, releps = 0)

  m<-length(nj)[1]

  N<-sum(nj)

  p1<-dim(x)[2]

  q1<-dim(z)[2]

  ver<-matrix(0,m,1)



  if(cens.type=="left"){

    LI<-rep(-Inf,length(cc))

    LS<-rep(Inf,length(cc))

    LS[cc==1]<-y[cc==1]

    LI<-as.vector(LI)

    LS<-as.vector(LS)

  }



  if(cens.type=="right"){

    LI<-rep(-Inf,length(cc))

    LI[cc==1]<-y[cc==1]

    LS<-rep(Inf,length(cc))

    LI<-as.vector(LI)

    LS<-as.vector(LS)

  }



  if(cens.type=="both"){

    LI<-LI

    LS<-LS

    LI<-as.vector(LI)

    LS<-as.vector(LS)

  }





  if(distr=="Normal"){



    for (j in 1:m){

      cc1<-cc[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))]

      y1<- y[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))]

      t1<- tt[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))]

      x1<- matrix(x[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),  ],ncol=p1)

      z1<- matrix(z[(sum(nj[1:j-1])+1) : (sum(nj[1:j])) ,  ],ncol=q1)

      LI1<- LI[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))]

      LS1<- LS[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))]

      gammai<- x1%*%beta1

      if (depstr == "AR1" || depstr == "ARp") {
        Gama<- sigmae*CovARp(phi = phi, ti = t1)#MatArpphi(phi,nj[j],)
      } else if (depstr == "UNC") {
        Gama<- sigmae*diag(nj[j])
      } else if (depstr == "DEC" || depstr == "CAR") {
        Gama<- sigmae*CovDEC(phi1 = phi[1], phi2 = phi[2], ti = t1)
      }

      Gama<-(Gama+t(Gama))/2

      Psi<- (Gama+(z1)%*%DD%*%t(z1))

      Psi<- (Psi+t(Psi))/2



      if(sum(cc1)==0){

        ver[j]<- mvtnorm::dmvnorm(c(y1),mean=c(gammai),sigma=Psi)

      }



      if(sum(cc1)>=1){



        if(sum(cc1)==nj[j]){

          ver[j]<- mvtnorm::pmvnorm(LI1, LS1, mean=c(gammai),sigma=Psi,algorithm = GB)

        }



        else {



          muc<-x1[cc1==1,]%*%beta1+Psi[cc1==1,cc1==0]%*%ginv(Psi[cc1==0,cc1==0])%*%(y1[cc1==0]-x1[cc1==0,]%*%beta1)

          Sc<-Psi[cc1==1,cc1==1]-Psi[cc1==1,cc1==0]%*%ginv(Psi[cc1==0,cc1==0])%*%Psi[cc1==0,cc1==1]

          Sc<-(Sc+t(Sc))/2

          So<-(Psi[cc1==0,cc1==0]+t(Psi[cc1==0,cc1==0]))/2

          ver[j]<- mvtnorm::dmvnorm(c(y1[cc1==0]),mean=c(gammai[cc1==0]),sigma=So)*(pmvnorm(LI1[cc1==1],LS1[cc1==1],mean=as.vector(muc),sigma=Sc,algorithm = GB))[1]
        }
      }
    }
  }

  if(distr=="Student"){

    for (j in 1:m){

      cc1<-cc[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))]

      y1<- y[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))]

      t1<- tt[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))]

      x1<- matrix(x[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),  ],ncol=p1)

      z1<- matrix(z[(sum(nj[1:j-1])+1) : (sum(nj[1:j])) ,  ],ncol=q1)

      LI1<- LI[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))]

      LS1<- LS[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))]

      gammai<- x1%*%beta1

      if (depstr == "AR1" || depstr == "ARp") {
        Gama<- sigmae*CovARp(phi = phi, ti = t1)#MatArpphi(phi,nj[j],)
      } else if (depstr == "UNC") {
        Gama<- sigmae*diag(nj[j])
      } else if (depstr == "DEC" || depstr == "CAR") {
        Gama<- sigmae*CovDEC(phi1 = phi[1], phi2 = phi[2], ti = t1)
      }

      Gama<-(Gama+t(Gama))/2

      Psi<- (Gama+(z1)%*%DD%*%t(z1))

      Psi<- (Psi+t(Psi))/2



      if(sum(cc1)==0){

        ver[j]<- tmvtnorm::dtmvt((y1), mean=c(gammai),sigma=Psi,df=nu,lower=rep(-Inf,length(mean)),upper=rep(Inf,length(mean)))

      }



      if(sum(cc1)>=1){



        if(sum(cc1)==nj[j]){

          ver[j]<- tmvtnorm::ptmvt(lowerx=(LI1), upperx=(LS1), mean=c(gammai),  sigma=Psi, df=nu)

        }



        else {



          nu1<-(nu+length(cc1[cc1==0]))

          muc<-gammai[cc1==1]+Psi[cc1==1,cc1==0]%*%ginv(Psi[cc1==0,cc1==0])%*%(y1[cc1==0]-gammai[cc1==0])

          Sc <-Psi[cc1==1,cc1==1]-Psi[cc1==1,cc1==0]%*%ginv(Psi[cc1==0,cc1==0])%*%Psi[cc1==0,cc1==1]

          Sc.1<-(Sc+t(Sc))/2

          Qy1<-t(y1[cc1==0]-gammai[cc1==0])%*%ginv(Psi[cc1==0,cc1==0])%*%(y1[cc1==0]-gammai[cc1==0])

          auxcte<-(nu+as.numeric(Qy1))/(nu+length(cc1[cc1==0]))

          Sc22<-auxcte*Sc.1

          SigmaUi<-Sc22

          SigmaUi.1<-(SigmaUi+t(SigmaUi))/2

          So<-(Psi[cc1==0,cc1==0]+t(Psi[cc1==0,cc1==0]))/2

          ver[j]<- mvtnorm::dmvt((y1[cc1==0]),delta=c(gammai[cc1==0]),sigma=So,df=nu,log=FALSE)*pmvt(lower=LI1[cc1==1], upper=LS1[cc1==1], delta=c(muc), df=nu1, sigma=SigmaUi.1,algorithm = GB)[1]

        }
      }
    }
  }



  if(distr=="Slash") {


    for (j in 1:m ){

      cc1<-cc[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))]

      y1<- y[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))]

      t1<- tt[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))]

      x1<- matrix(x[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),  ],ncol=p1)

      z1<- matrix(z[(sum(nj[1:j-1])+1) : (sum(nj[1:j])) ,  ],ncol=q1)

      LI1<- LI[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))]

      LS1<- LS[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))]


      beta1 <- as.matrix(beta1)

      DD <- as.matrix(DD)

      gammai<- x1%*%beta1

        if (depstr == "AR1" || depstr == "ARp") {
          Gama<- sigmae*CovARp(phi = phi, ti = t1)#MatArpphi(phi,nj[j],)
        } else if (depstr == "UNC") {
          Gama<- sigmae*diag(nj[j])
        } else if (depstr == "DEC" || depstr == "CAR") {
          Gama<- sigmae*CovDEC(phi1 = phi[1], phi2 = phi[2], ti = t1)
        }

        Gama<-(Gama+t(Gama))/2

        Psi<- (Gama+(z1)%*%DD%*%t(z1))

        Psi<- (Psi+t(Psi))/2



        if(sum(cc1)==0){

          ver[j]<-dSL(c(y1),c(gammai),Psi,nu)
        }

        if(sum(cc1)>0){

          if(sum(cc1)==nj[j]){
            muc= x1%*%beta1
            Sc = Psi
            Bi<- diag(sqrt(diag(Sc)),sum(cc1),sum(cc1))
            Ri<- solve(Bi)%*%Sc%*%solve(Bi)
            ai<- solve(Bi)%*%(y1-muc)
            ver[j]<-acumSL(c(ai),rep(0,sum(cc1)),Ri,nu)
          }

          else {
            muc<-x1[cc1==1,]%*%beta1+Psi[cc1==1,cc1==0]%*%solve(Psi[cc1==0,cc1==0])%*%(y1[cc1==0]-x1[cc1==0,]%*%beta1)
            Sc<- Psi[cc1==1,cc1==1]-Psi[cc1==1,cc1==0]%*%solve(Psi[cc1==0,cc1==0])%*%Psi[cc1==0,cc1==1]

            ver[j]<-grafSL(c(y1[cc1==0]),c(gammai[cc1==0]),Psi[cc1==0,cc1==0],y1[cc1==1],c(muc),Sc,nu)


          }
        }

      }
    }



  return(ver)

}

# distr = Normal, Student, Slash

##Model comparison criterion


criteria<-function(cc, nj, y, x, z, tt, espac=5, stanobj, distr="Normal", depstr = "UNC",
                    cens.type="left",LI=NULL, LS=NULL){

  if(distr == "Normal"){

    if (depstr == "AR1") {
      cadeia = list(BetaF =rstan::extract(stanobj, "beta")[[1]],
                  sigmaeF = rstan::extract(stanobj, "sigma2")[[1]],
                  DD = rstan::extract(stanobj, "D1"),
                  phiF = rstan::extract(stanobj, "phi")[[1]])
    } else if (depstr == "ARp") {
      cadeia = list(BetaF =rstan::extract(stanobj, "beta")[[1]],
                  sigmaeF = rstan::extract(stanobj, "sigma2")[[1]],
                  DD = rstan::extract(stanobj, "D1"),
                  phiF = rstan::extract(stanobj, "phi")[[1]])
    } else if (depstr == "UNC") {
      cadeia = list(BetaF =rstan::extract(stanobj, "beta")[[1]],
                  sigmaeF = rstan::extract(stanobj, "sigma2")[[1]],
                  DD = rstan::extract(stanobj, "sigmab2")[[1]])
    } else if (depstr == "DEC") {
      cadeia = list(BetaF =rstan::extract(stanobj, "beta")[[1]],
                  sigmaeF = rstan::extract(stanobj, "sigma2")[[1]],
                  DD = rstan::extract(stanobj, "D1"),
                  phiF = cbind(rstan::extract(stanobj, "phi1")[[1]],
                               rstan::extract(stanobj, "phi2")[[1]]))
    } else if (depstr == "CAR") {
      cadeia = list(BetaF =rstan::extract(stanobj, "beta")[[1]],
                    sigmaeF = rstan::extract(stanobj, "sigma2")[[1]],
                    DD = rstan::extract(stanobj, "D1"),
                    phiF = cbind(rstan::extract(stanobj, "phi1")[[1]],
                                 rep(1,length(rstan::extract(stanobj, "phi1")[[1]]))))
    }
  }

  if (distr == "Student" || distr == "Slash"){
    if (depstr == "AR1") {
      cadeia = list(BetaF =rstan::extract(stanobj, "beta")[[1]],
                    sigmaeF = rstan::extract(stanobj, "sigma2")[[1]],
                    DD = rstan::extract(stanobj, "D1"),
                    phiF = rstan::extract(stanobj, "phi")[[1]],
                    nuF = rstan::extract(stanobj, "nu")[[1]])
    } else if (depstr == "ARp") {
      cadeia = list(BetaF =rstan::extract(stanobj, "beta")[[1]],
                    sigmaeF = rstan::extract(stanobj, "sigma2")[[1]],
                    DD = rstan::extract(stanobj, "D1"),
                    phiF = rstan::extract(stanobj, "phi")[[1]],
                    nuF = rstan::extract(stanobj, "nu")[[1]])
    } else if (depstr == "UNC") {
      cadeia = list(BetaF =rstan::extract(stanobj, "beta")[[1]],
                    sigmaeF = rstan::extract(stanobj, "sigma2")[[1]],
                    DD = rstan::extract(stanobj, "sigmab2")[[1]],
                    nuF = rstan::extract(stanobj, "nu")[[1]])
    } else if (depstr == "DEC") {
      cadeia = list(BetaF =rstan::extract(stanobj, "beta")[[1]],
                    sigmaeF = rstan::extract(stanobj, "sigma2")[[1]],
                    DD = rstan::extract(stanobj, "D1"),
                    phiF = cbind(rstan::extract(stanobj, "phi1")[[1]],
                                 rstan::extract(stanobj, "phi2")[[1]]),
                    nuF = rstan::extract(stanobj, "nu")[[1]])
    } else if (depstr == "CAR") {
      cadeia = list(BetaF =rstan::extract(stanobj, "beta")[[1]],
                    sigmaeF = rstan::extract(stanobj, "sigma2")[[1]],
                    DD = rstan::extract(stanobj, "D1"),
                    phiF = cbind(rstan::extract(stanobj, "phi1")[[1]],
                                 rep(1,length(rstan::extract(stanobj, "phi1")[[1]]))),
                    nuF = rstan::extract(stanobj, "nu")[[1]])
    }
  }

  m<-length(nj)

  if (missing(tt)) tt <- unlist(lapply(nj,seq_len))



  if(cens.type=="left"){

    LI<-rep(-Inf,length(cc))

    LS<-rep(Inf,length(cc))

    LS[cc==1]<-y[cc==1]

    LI<-as.vector(LI)

    LS<-as.vector(LS)

  }



  if(cens.type=="right"){

    LI<-rep(-Inf,length(cc))

    LI[cc==1]<-y[cc==1]

    LS<-rep(Inf,length(cc))

    LI<-as.vector(LI)

    LS<-as.vector(LS)

  }



  if(cens.type=="both"){

    LI<-LI

    LS<-LS

    LI<-as.vector(LI)

    LS<-as.vector(LS)

  }


  if (depstr == "UNC") {
    ntot <-length(cadeia$DD)

    DDs <- Reduce("+", cadeia$DD)/ntot

    DDs <- diag(DDs,dim(z)[2],dim(z)[2])
  }

  else {

  ntot<- dim(cadeia$DD$D1)[1]

  DDs <- matrix(0, nrow = dim(cadeia$DD$D1)[2], ncol = dim(cadeia$DD$D1)[3])

    for(i in 1:ntot){

    DDs <- DDs + cadeia$DD$D1[i,,]

    }

  DDs <- DDs/ntot
  DDs <- (DDs+t(DDs))/2

  }



    if (depstr=="AR1") {
      phis <- mean(cadeia$phiF)
    }  else if (depstr == "UNC") {
      phis = NULL
    } else if (depstr == "DEC"||depstr == "ARp" || depstr == "CAR") {
      phis <- apply(cadeia$phiF,2,mean)
    }

  if (distr == "Normal"){

    ver<-VerCensLMM(cc = cc, nj = nj, y = y, x = x, z = z, beta1 = apply(cadeia$BetaF,2,mean), sigmae = mean(cadeia$sigmaeF),
                    DD = DDs , tt = tt, phi = phis,nu =NULL,distr = distr, cens.type = cens.type ,depstr = depstr, LI = LI, LS = LS)

    iter<-floor(ntot/espac)

    Loglikaux<-matrix(0,m,iter)

    CPOaux<-matrix(0, m,iter)


    for(k in 1:iter){

      i<-espac*k

      betasF<-cadeia$BetaF[i,]

      sigmasF<-cadeia$sigmaeF[i]

      if (depstr == "UNC"){
        DDF<- cadeia$DD[i]

        DDF <- diag(DDF,dim(z)[2],dim(z)[2])
      }

      else {
        DDF<-cadeia$DD$D1[i,,]
      }
      DDF<-(DDF+t(DDF))/2

      if (depstr == "AR1") {
        phisF<-cadeia$phiF[i]
      } else if (depstr == "ARp"){
        phisF<-cadeia$phiF[i,]
      } else if (depstr == "UNC") {
        phisF<-NULL
      } else if (depstr == "DEC" || depstr == "CAR") {
        phisF<-cadeia$phiF[i,]
      }

      Loglikaux[,k]<- VerCensLMM(cc = cc, nj = nj,y =  y, x = x, z = z, beta1 = betasF,sigmae = sigmasF, DD = DDF, tt = tt,
                                 phi = phisF, nu = NULL, distr = distr,cens.type =  cens.type,depstr = depstr, LI = LI, LS = LS)

      CPOaux[,k]<-1/Loglikaux[,k]

    }

    Np<-length(c(betasF,sigmasF, DDF[upper.tri(DDF, diag = T)], phisF))

  }


  if (distr == "Student" || distr == "Slash"){



    ver <- VerCensLMM(cc = cc, nj = nj, y = y, x = x, z = z,
                         beta1 = apply(cadeia$BetaF,2,mean), sigmae =  mean(cadeia$sigmaeF),
                         DD = DDs , tt = tt, phi = phis, nu = ceiling(mean(cadeia$nuF)),
                         distr = distr, cens.type = cens.type, LI = LI, LS = LS)

    iter<-floor(ntot/espac)

    Loglikaux<-matrix(0,m,iter)

    CPOaux<-matrix(0, m,iter)


    for(k in 1:iter){

      i<-espac*k

      betasF<-cadeia$BetaF[i,]

      sigmasF<-cadeia$sigmaeF[i]

      if (depstr == "UNC"){
        DDF<- cadeia$DD[i]

        DDF <- diag(DDF,dim(z)[2],dim(z)[2])
      }

      else {
        DDF<-cadeia$DD$D1[i,,]
      }

      DDF<-(DDF+t(DDF))/2

      if (depstr == "AR1") {
        phisF <- cadeia$phiF[i]
      } else if (depstr == "ARp"){
        phisF <- cadeia$phiF[i,]
      } else if (depstr == "UNC") {
        phisF <- NULL
      } else if (depstr == "DEC" || depstr == "CAR") {
        phisF <- cadeia$phiF[i,]
      }

      nu1<-cadeia$nuF[i]

      Loglikaux[,k]<- VerCensLMM(cc = cc, nj = nj,y =  y, x = x, z = z, beta1 = betasF,sigmae = sigmasF, DD = DDF, tt = tt,
                                 phi = phisF, nu = ceiling(nu1), distr = distr,cens.type =  cens.type,depstr = depstr, LI = LI, LS = LS)

      CPOaux[,k]<-1/Loglikaux[,k]

    }

    Np<-length(c(betasF,sigmasF, DDF[upper.tri(DDF, diag = T)], phisF, nu1))

  }


  CPO<- sum(log(1/(apply(CPOaux,1,mean))))

  DIC<- 2*mean(apply(-2*log(Loglikaux),2,sum))-sum(2*log(ver))

  EAIC<- mean(apply(-2*log(Loglikaux),2,sum))+2*Np

  EBIC<- mean(apply(-2*log(Loglikaux),2,sum))+log(sum(nj))*Np

  IL<- -log(1/(apply(CPOaux,1,mean)))+ apply(log(Loglikaux),1,mean)

  return(list("LPML" = CPO, "DIC" = DIC, "EAIC" = EAIC, "EBIC" = EBIC, "KL" = IL))

}



