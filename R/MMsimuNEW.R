

is.wholenumber <- function(x, tol1 = .Machine$double.eps^0.5)  abs(x - round(x)) < tol1


matrix.sqrt <- function(A){
  if (length(A)==1) return(sqrt(A))
  else{
    sva <- svd(A)
    if (min(sva$d)>=0) {
      Asqrt <- sva$u%*%diag(sqrt(sva$d))%*%t(sva$v) # svd e decomposi??o espectral
      if (all(abs(Asqrt%*%Asqrt-A)<1e-4)) return(Asqrt)
      else stop("Matrix square root is not defined/not real")
    }
    else stop("Matrix square root is not defined/not real")
  }
}


# Transformation function: phi to pi
tphitopi <- function(phit) {
  p <- length(phit)
  Phi <- matrix(0,ncol=p,nrow=p)
  Phi[p,] <- phit
  if (p>1) {
    for (k in p:2) {
      for (i in 1:(k-1)) {
        Phi[k-1,i] <- (Phi[k,i] + Phi[k,k]*Phi[k,k-i])/(1-Phi[k,k]^2)
      }
    }
    return(diag(Phi))
  }
  else return(phit)
}

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

# Corr matrix of comp symmetry (CS)
CovCS <- function(phi, n) {
  if (n==1) Rn <- matrix(1)
  else Rn <- toeplitz(c(1,rep(phi,n-1)))
  return(Rn)
}



errorVar<- function(times,object=NULL,sigma2=NULL,depStruct=NULL,phi=NULL) {
  if((!is.null(object))&&(!inherits(object,c("SMSN","SMN")))) stop("object must inherit from class SMSN or SMN")
  if (is.null(object)&&is.null(depStruct)) stop("object or depStruct must be provided")
  if (is.null(object)&&is.null(sigma2)) stop("object or sigma2 must be provided")
  if (is.null(depStruct)) depStruct<-object$depStruct
  if (depStruct=="CI") depStruct = "UNC"
  if (depStruct!="UNC" && is.null(object)&is.null(phi)) stop("object or phi must be provided")
  if (!(depStruct %in% c("UNC","ARp","CS","DEC","CAR1"))) stop("accepted depStruct: UNC, ARp, CS, DEC or CAR1")
  if (is.null(sigma2)) sigma2<-object$estimates$sigma2
  if (is.null(phi)&&depStruct!="UNC") phi<-object$estimates$phi
  if (depStruct=="ARp" && (any(!is.wholenumber(times))|any(times<=0))) stop("times must contain positive integer numbers when using ARp dependency")
  if (depStruct=="ARp" && any(tphitopi(phi)< -1|tphitopi(phi)>1)) stop("AR(p) non stationary, choose other phi")
  #
  if (depStruct=="UNC") var.out<- sigma2*diag(length(times))
  if (depStruct=="ARp") var.out<- sigma2*CovARp(phi,times)
  if (depStruct=="CS") var.out<- sigma2*CovCS(phi,length(times))
  if (depStruct=="DEC") var.out<- sigma2*CovDEC(phi[1],phi[2],times)
  if (depStruct=="CAR1") var.out<- sigma2*CovDEC(phi,1,times)
  var.out
}

rsmsn.lmm <- function(time1,x1,z1,sigma2,D1,beta,lambda,depStruct="UNC",phi=NULL,distr="sn",nu=NULL) {
  if (length(D1)==1 && !is.matrix(D1)) D1=as.matrix(D1)
  q1 = nrow(D1)
  p = length(beta)
  if (ncol(as.matrix(x1))!=p) stop("incompatible dimension of x1/beta")
  if (ncol(as.matrix(z1))!=q1) stop ("incompatible dimension of z1/D1")
  if (length(lambda)!=q1) stop ("incompatible dimension of lambda/D1")
  if (!is.matrix(D1)) stop("D must be a matrix")
  if ((ncol(D1)!=q1)|(nrow(D1)!=q1)) stop ("wrong dimension of D")
  if (length(sigma2)!=1) stop ("wrong dimension of sigma2")
  if (sigma2<=0) stop("sigma2 must be positive")
  Sig <- errorVar(time1,depStruct = depStruct,sigma2=sigma2,phi=phi)
  #
  if (distr=="ssl") distr<-"ss"
  if (!(distr %in% c("sn","st","ss","scn"))) stop("Invalid distribution")
  if (distr=="sn") {ui=1; c.=-sqrt(2/pi)}
  if (distr=="st") {ui=rgamma(1,nu/2,nu/2); c.=-sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)}
  if (distr=="ss") {ui=rbeta(1,nu,1); c.=-sqrt(2/pi)*nu/(nu-.5)}
  if (distr=="scn") {ui=ifelse(runif(1)<nu[1],nu[2],1);
  c.=-sqrt(2/pi)*(1+nu[1]*(nu[2]^(-.5)-1))}
  #if (all(lambda==0)) c.=0
  delta = lambda/as.numeric(sqrt(1+t(lambda)%*%(lambda)))
  Delta = matrix.sqrt(D1)%*%delta
  Gammab = D1 - Delta%*%t(Delta)
  Xi = matrix(x1,ncol=p)
  Zi = matrix(z1,ncol=q1)
  Beta = matrix(beta,ncol=1)
  ti = c.+abs(rnorm(1,0,ui^-.5))
  bi = t(rmvnorm(1,Delta*ti,sigma=ui^(-1)*Gammab))
  Yi = t(rmvnorm(1,Xi%*%Beta+Zi%*%bi,sigma=ui^(-1)*Sig))
  if (all(Xi[,1]==1)) Xi = Xi[,-1]
  if (all(Zi[,1]==1)) Zi = Zi[,-1]
  return(data.frame(time=time1,y=Yi,x=Xi,z=Zi))
}

MMsimu.mod=function(m,x,z,tt,nj,beta,sigmae,D,phi,struc,typeModel,percCensu,nivel.Censu,cens.type,nu){

  y <- rep(0,length(tt))

  subject_log <- 0

  for(i in 1:m) {

    if(typeModel == "Normal"){

      data_i <-rsmsn.lmm(tt[(subject_log + 1) : (subject_log + nj[i])],
                       x[(subject_log + 1) : (subject_log + nj[i]),],
                       z[(subject_log + 1) : (subject_log + nj[i]),],
                       sigmae, D, beta, rep(0,dim(D)[1]), depStruct = struc,
                       phi = phi, distr = "sn", nu = NULL)

    }

    if(typeModel == "Student") {

      data_i <-rsmsn.lmm(tt[(subject_log + 1) : (subject_log + nj[i])],
                         x[(subject_log + 1) : (subject_log + nj[i]),],
                         z[(subject_log + 1) : (subject_log + nj[i]),],
                         sigmae, D, beta, rep(0,dim(D)[1]), depStruct = struc,
                         phi = phi, distr = "st", nu = nu)
    }

    if(typeModel == "CN") {
      data_i <-rsmsn.lmm(tt[(subject_log + 1) : (subject_log + nj[i])],
                         x[(subject_log + 1) : (subject_log + nj[i]),],
                         z[(subject_log + 1) : (subject_log + nj[i]),],
                         sigmae, D, beta, rep(0,dim(D)[1]), depStruct = struc,
                         phi = phi, distr = "scn", nu = nu)
    }

    if(typeModel == "Slash") {
      data_i <-rsmsn.lmm(tt[(subject_log + 1) : (subject_log + nj[i])],
                         x[(subject_log + 1) : (subject_log + nj[i]),],
                         z[(subject_log + 1) : (subject_log + nj[i]),],
                         sigmae, D, beta, rep(0,dim(D)[1]), depStruct = struc,
                         phi = phi, distr = "ss", nu = nu)
    }


    y[(subject_log + 1):(subject_log + nj[i])] <- data_i[,2]

    subject_log <- subject_log + nj[i]
  }

  yy=y
  y_cc=y
  cc=rep(0,length(y_cc))

  if(!is.null(percCensu))
  {
    if(percCensu!=0)
    {
      if(cens.type=="left") {
        aa=sort(y, decreasing = FALSE)
        bb=aa[1:(percCensu*sum(nj))]
        cutof<-bb[percCensu*sum(nj)]
        cc=matrix(1,sum(nj),1)*(y< cutof)
        y[cc==1]=cutof
        y_cc=y
      }

      if(cens.type=="right") {
        aa=sort(y, decreasing = TRUE)
        bb=aa[1:(percCensu*sum(nj))]
        cutof<-bb[percCensu*sum(nj)]
        cc=matrix(1,sum(nj),1)*(y> cutof)
        y[cc==1]=cutof
        y_cc=y
      }

      if(cens.type=="interval") {
        aa=sort(y, decreasing = FALSE)
        bbi=aa[1:(percCensu*sum(nj)*0.5)]
        aa=sort(y, decreasing = TRUE)
        bbs=aa[1:(percCensu*sum(nj)*0.5)]
        cutofi<-bbi[percCensu*sum(nj)*0.5]
        cutofs<-bbs[percCensu*sum(nj)*0.5]
        cci=matrix(1,sum(nj),1)*(y< cutofi)
        y[cci==1]=cutofi
        ccs=matrix(1,sum(nj),1)*(y>cutofs)
        y[ccs==1]=cutofs
        y_cc=y
        cc=cci+ccs
      }
    }
  }

  if(!is.null(nivel.Censu))
  {
    if(length(nivel.Censu)==1)
    {

      if(cens.type=="left") {
        cutof<-nivel.Censu
        cc=matrix(1,sum(nj),1)*(y< cutof)
        y[cc==1]=cutof
        y_cc=y
      }

      if(cens.type=="right") {
        cutof<-nivel.Censu
        cc=matrix(1,sum(nj),1)*(y> cutof)
        y[cc==1]=cutof
        y_cc=y
      }

    }

    if(length(nivel.Censu)>1)
    {
      if(cens.type=="interval") {
        cutofi<-nivel.Censu[1]
        cutofs<-nivel.Censu[2]
        cci=matrix(1,sum(nj),1)*(y< cutofi)
        y[cci==1]=cutofi
        ccs=matrix(1,sum(nj),1)*(y>cutofs)
        y[ccs==1]=cutofs
        y_cc=y
        cc=cci+ccs
      }
    }
  }

  # return(list(cc=cc, y_cc=y_cc, y_origin = yy))
  return(list(cc=cc, y_cc=y_cc))
}




