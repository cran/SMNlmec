




#' @title Bayesian Censored Mixed-Effects Models with Damped Exponential Correlation Structures for Scale Mixture of Normal distributions error
#' @import rstan
#' @import StanHeaders
#' @import MASS
#' @import tmvtnorm
#' @import mnormt
#' @import stats
#' @importFrom mvtnorm pmvnorm dmvnorm rmvnorm dmvt
#' @importFrom LaplacesDemon dmvn is.positive.definite as.inverse logdet
#' @importFrom TruncatedNormal pmvt mvNcdf mvNqmc
#' @importFrom numDeriv jacobian
#' @import methods
#' @description This function fits left, right censored mixed-effects linear model, with scale mixture of normal distribution errors, using the Stan. It returns estimates, standard errors and LPML, AIC, BIC and DIC.
#' @param ID Vector \code{N x 1} of the ID of the data set, specifying the ID for each measurement.
#' @param x_set Design matrix of the fixed effects of order \code{N x p}.
#' @param z_set Design matrix of the random effects of order \code{N x d}.
#' @param tt Vector \code{N x 1} with the time the measurements were made, where \code{N} is the total number of measurements for all individuals. Default it's considered regular times.
#' @param y_complete Vector \code{N x 1} of the complete responses.
#' @param censor_vector Vector \code{N x 1} of the indicator vector of censored responses.
#' @param dist Distribution of the random effects and random error. Available options are \code{Normal}, \code{Student} and \code{Slash}.
#' @param struc Structure of the correlation structure. Available options are \code{UNC}, \code{DEC}, \code{CAR}.
#' @param direction Direction of censoring type. Available options are \code{left} and \code{right}.
#' @param thin_num A positive integer specifying the period for saving samples. The default is 5. See more details in rstan::stan().
#' @param chains_num A positive integer specifying the number of chains generating by rstan::stan(). The default is 3.
#' @param iter_num A positive integer specifying the number of iterations for each chain (including warmup). The default is 5000.
#' @param burn_percen A percentage of the warm-up iterations in each chain the Stan. The default is 0.1.
#' @param seed_set A random seed. The default is NULL.
#' @param adapt_delta_set A parameter to control the sampler's behavior. The default is 0.8. See rstan::stan() for more details.
#' @return Return a S4 class SMNlmecfit object. Using function \code{SMNlmec.summary()} to obtain the estimation of parameters and model selection criteria. The SMNlmecfit include:
#' \item{stan_object}{A stanfit object from rstan::stan().}
#' \item{model_criteria}{A list includes LPML, DIC, EAIC, EBIC, K-L divergence.}
#' \item{dist_set}{The setting of distribution of the stan model.}
#' \item{struc_set}{The setting of correlation structure of the stan model.}
#' @references Kelin Zhong, Fernanda L. Schumacher, Luis M. Castro and Victor H. Lachos. Bayesian analysis of censored linear mixed-effects models for heavy-tailed  irregularly observed repeated measures. Statistics in Medicine, 2025. doi:10.1002/sim.10295
#' @examples
#' \donttest{
#' require(rstan)
#' require(StanHeaders)
#' require(MASS)
#' require(tmvtnorm)
#' require(mvtnorm)
#' require(mnormt)
#'
#' data("UTIdata_sub")
#' data1 <- UTIdata_sub
#' y1 <- c(log10(data1$RNA))
#' cc <- (data1$RNAcens==1)+0
#' y_com<-as.numeric(y1)
#' rho_com<-as.numeric(cc)
#' x <- cbind(
#'  (data1$Fup==0)+0,
#'  (data1$Fup==1)+0,
#'  (data1$Fup==3)+0,
#'  (data1$Fup==6)+0,
#'  (data1$Fup==9)+0,
#'  (data1$Fup==12)+0,
#'  (data1$Fup==18)+0,
#'  (data1$Fup==24)+0
#'  )
#' z <- matrix(rep(1, length(y1)), ncol=1)
#'
#' UTI_T_DEC <- SMNlmec.est(ID = data1$Patid, x_set = x, z_set = z,
#'                          tt = data1$Fup, y_complete = y_com,
#'                          censor_vector = rho_com, dist = "Student",
#'                          struc = "DEC", direction = "left",
#'                          thin_num = 1, chains_num = 1, iter_num = 3000,
#'                          burn_percen = 0.1, seed_set = 9955,
#'                          adapt_delta_set = 0.8)
#'
#' SMNlmec.summary(UTI_T_DEC)
#' }
#'
#' @export


SMNlmec.est <- function(ID, x_set, z_set, tt, y_complete,
                        censor_vector, dist = "Normal",
                        struc = "UNC", direction = "left",
                        thin_num = 1, chains_num = 1, iter_num = 3000,
                        burn_percen = 0.1, seed_set = NULL,
                        adapt_delta_set = 0.8) {

  y_com <- y_complete
  rho_com <- censor_vector
  x <- x_set
  z <- z_set

  subjects <- unique(ID)
  cluster <- c(match(ID,subjects))
  m <- length(subjects)
  N_com <-length(cluster)
  cens_N <- sum(rho_com)
  obs_N <- N_com - cens_N

  l_set <- dim(x)[2]
  q1_set <- dim(z)[2]
  ycen <- y_com[rho_com == 1]
  yobs <- y_com[rho_com == 0]

  nj<- rep(0,m)
  for(j in 1:m) {
    nj[j] <- sum(cluster == j)
  }

  ind_set<-numeric()
  log_nj<-0

  for(i in 1:length(nj)) {
    for(j in 1:nj[i]){
      ind_set[log_nj+j] <- i
    }
    log_nj <- log_nj + nj[i]
  }

  cens_nj <- numeric()
  log_nj <- 0
  cens_count <- 0

  for(i in 1:length(nj)) {
    for(j in 1:nj[i]){
      if(rho_com[log_nj+j] == 1) {
        cens_count <- cens_count+1
      }
    }
    log_nj <- log_nj + nj[i]
    cens_nj[i] <- cens_count
    cens_count <- 0
  }


  if(dist == "Normal") {
    if(struc == "UNC")   {
      if(direction == "left") {

        stan_file <- system.file("stan", "lmm-UNC-N-censored-left.stan", package = "SMNlmec")

        stan_obj <- rstan::stan(file= stan_file,
                         data = list(N_complete = N_com, N_cens=cens_N, N_obs= obs_N,
                                     n = m, l = l_set,yobs=yobs,ycen=ycen,x=x,
                                     ind=ind_set,rho=rho_com),
                         thin = thin_num, chains = chains_num,
                         iter = iter_num, warmup = iter_num*burn_percen,
                         seed = seed_set, control = list(adapt_delta=adapt_delta_set))

        SMNlmec_criteria <- criteria(cc = rho_com, nj = nj, y = y_com, x = x, z = z,
                              tt = tt, espac=5,
                              stanobj = stan_obj, distr="Normal",
                              depstr = "UNC", cens.type="left",LI=NULL, LS=NULL)

      }

      if(direction == "right") {

        stan_file <- system.file("stan", "lmm-UNC-N-censored-right.stan", package = "SMNlmec")

        stan_obj <- rstan::stan(file= stan_file,
                         data = list(N_complete = N_com, N_cens=cens_N, N_obs= obs_N,
                                     n = m, l = l_set,yobs=yobs,ycen=ycen,x=x,
                                     ind=ind_set,rho=rho_com),
                         thin = thin_num, chains = chains_num,
                         iter = iter_num, warmup = iter_num*burn_percen,
                         seed = seed_set, control = list(adapt_delta=adapt_delta_set))

        SMNlmec_criteria <- criteria(cc = rho_com, nj = nj, y = y_com, x = x, z = z,
                                     tt = tt, espac=5,
                                     stanobj = stan_obj, distr="Normal",
                                     depstr = "UNC", cens.type="right",LI=NULL, LS=NULL)

      }
    }

    if(struc == "DEC")   {
      if(direction == "left") {

        stan_file <- system.file("stan", "lmm-DEC-N-censored-left.stan", package = "SMNlmec")

        stan_obj <- rstan::stan(file= stan_file,
                              data = list(N_complete = N_com, N_cens = cens_N, ycen=ycen,
                                          n=m, l=l_set, q1 = q1_set,
                                          x=x, z=z, timevar = tt,
                                          y_complete = y_com, rho = rho_com,
                                          njvec =nj, cens_nj = cens_nj,
                                          ind = ind_set),
                              thin = thin_num, chains = chains_num,
                              iter = iter_num, warmup = iter_num*burn_percen,
                              seed = seed_set, control = list(adapt_delta=adapt_delta_set))


        SMNlmec_criteria <- criteria(cc = rho_com, nj = nj, y = y_com, x = x, z = z,
                              tt = tt, espac=5,
                              stanobj = stan_obj, distr="Normal",
                              depstr = "DEC", cens.type="left",LI=NULL, LS=NULL)
      }

      if(direction == "right") {

        stan_file <- system.file("stan", "lmm-DEC-N-censored-right.stan", package = "SMNlmec")

        stan_obj <- rstan::stan(file= stan_file,
                         data = list(N_complete = N_com, N_cens = cens_N, ycen=ycen,
                                     n=m, l=l_set, q1 = q1_set,
                                     x=x, z=z, timevar = tt,
                                     y_complete = y_com, rho = rho_com,
                                     njvec =nj, cens_nj = cens_nj,
                                     ind = ind_set),
                         thin = thin_num, chains = chains_num,
                         iter = iter_num, warmup = iter_num*burn_percen,
                         seed = seed_set, control = list(adapt_delta=adapt_delta_set))


        SMNlmec_criteria <- criteria(cc = rho_com, nj = nj, y = y_com, x = x, z = z,
                                     tt = tt, espac=5,
                                     stanobj = stan_obj, distr="Normal",
                                     depstr = "DEC", cens.type="right",LI=NULL, LS=NULL)
      }
    }

    if(struc == "CAR")   {
      if(direction == "left") {

        stan_file <- system.file("stan", "lmm-CAR-N-censored-left.stan", package = "SMNlmec")

        stan_obj <- rstan::stan(file= stan_file,
                         data = list(N_complete = N_com, N_cens = cens_N, ycen=ycen,
                                     n=m, l=l_set, q1 = q1_set,
                                     x=x, z=z, timevar = tt,
                                     y_complete = y_com, rho = rho_com,
                                     njvec =nj, cens_nj = cens_nj,
                                     ind = ind_set),
                         thin = thin_num, chains = chains_num,
                         iter = iter_num, warmup = iter_num*burn_percen,
                         seed = seed_set, control = list(adapt_delta=adapt_delta_set))


        SMNlmec_criteria <- criteria(cc = rho_com, nj = nj, y = y_com, x = x, z = z,
                                     tt = tt, espac=5,
                                     stanobj = stan_obj, distr="Normal",
                                     depstr = "CAR", cens.type="left",LI=NULL, LS=NULL)
      }

      if(direction == "right") {

        stan_file <- system.file("stan", "lmm-CAR-N-censored-right.stan", package = "SMNlmec")

        stan_obj <- rstan::stan(file= stan_file,
                         data = list(N_complete = N_com, N_cens = cens_N, ycen=ycen,
                                     n=m, l=l_set, q1 = q1_set,
                                     x=x, z=z, timevar = tt,
                                     y_complete = y_com, rho = rho_com,
                                     njvec =nj, cens_nj = cens_nj,
                                     ind = ind_set),
                         thin = thin_num, chains = chains_num,
                         iter = iter_num, warmup = iter_num*burn_percen,
                         seed = seed_set, control = list(adapt_delta=adapt_delta_set))


        SMNlmec_criteria <- criteria(cc = rho_com, nj = nj, y = y_com, x = x, z = z,
                                     tt = tt, espac=5,
                                     stanobj = stan_obj, distr="Normal",
                                     depstr = "CAR", cens.type="right",LI=NULL, LS=NULL)

      }
    }

  }

  if(dist == "Student") {
    if(struc == "UNC")   {
      if(direction == "left") {

        stan_file <- system.file("stan", "lmm-UNC-t-censored-left.stan", package = "SMNlmec")

        stan_obj <- rstan::stan(file= stan_file,
                         data = list(N_complete = N_com, N_cens=cens_N, N_obs= obs_N,
                                     n = m, l = l_set,yobs=yobs,ycen=ycen,x=x,
                                     ind=ind_set,rho=rho_com),
                         thin = thin_num, chains = chains_num,
                         iter = iter_num, warmup = iter_num*burn_percen,
                         seed = seed_set, control = list(adapt_delta=adapt_delta_set))

        SMNlmec_criteria <- criteria(cc = rho_com, nj = nj, y = y_com, x = x, z = z,
                              tt = tt, espac=5,
                              stanobj = stan_obj, distr="Student",
                              depstr = "UNC", cens.type="left",LI=NULL, LS=NULL)
      }

      if(direction == "right") {

        stan_file <- system.file("stan", "lmm-UNC-t-censored-right.stan", package = "SMNlmec")

        stan_obj <- rstan::stan(file= stan_file,
                         data = list(N_complete = N_com, N_cens=cens_N, N_obs= obs_N,
                                     n = m, l = l_set,yobs=yobs,ycen=ycen,x=x,
                                     ind=ind_set,rho=rho_com),
                         thin = thin_num, chains = chains_num,
                         iter = iter_num, warmup = iter_num*burn_percen,
                         seed = seed_set, control = list(adapt_delta=adapt_delta_set))

        SMNlmec_criteria <- criteria(cc = rho_com, nj = nj, y = y_com, x = x, z = z,
                                     tt = tt, espac=5,
                                     stanobj = stan_obj, distr="Student",
                                     depstr = "UNC", cens.type="right",LI=NULL, LS=NULL)
      }
    }

    if(struc == "DEC")   {
      if(direction == "left") {

        stan_file <- system.file("stan", "lmm-DEC-t-censored-left.stan", package = "SMNlmec")

        stan_obj <- rstan::stan(file= stan_file,
                              data = list(N_complete = N_com, N_cens = cens_N, ycen=ycen,
                                          n=m, l=l_set, q1 = q1_set,
                                          x=x, z=z, timevar = tt,
                                          y_complete = y_com, rho = rho_com,
                                          njvec =nj, cens_nj = cens_nj,
                                          ind = ind_set),
                              thin = thin_num, chains = chains_num,
                              iter = iter_num, warmup = iter_num*burn_percen,
                              seed = seed_set, control = list(adapt_delta=adapt_delta_set))

        SMNlmec_criteria <- criteria(cc = rho_com, nj = nj, y = y_com, x = x, z = z,
                              tt = tt, espac=5,
                              stanobj = stan_obj, distr="Student",
                              depstr = "DEC", cens.type="left",LI=NULL, LS=NULL)
      }

      if(direction == "right") {

        stan_file <- system.file("stan", "lmm-DEC-t-censored-right.stan", package = "SMNlmec")

        stan_obj <- rstan::stan(file= stan_file,
                         data = list(N_complete = N_com, N_cens = cens_N, ycen=ycen,
                                     n=m, l=l_set, q1 = q1_set,
                                     x=x, z=z, timevar = tt,
                                     y_complete = y_com, rho = rho_com,
                                     njvec =nj, cens_nj = cens_nj,
                                     ind = ind_set),
                         thin = thin_num, chains = chains_num,
                         iter = iter_num, warmup = iter_num*burn_percen,
                         seed = seed_set, control = list(adapt_delta=adapt_delta_set))

        SMNlmec_criteria <- criteria(cc = rho_com, nj = nj, y = y_com, x = x, z = z,
                                     tt = tt, espac=5,
                                     stanobj = stan_obj, distr="Student",
                                     depstr = "DEC", cens.type="right",LI=NULL, LS=NULL)
      }
    }

    if(struc == "CAR")   {
      if(direction == "left") {

        stan_file <- system.file("stan", "lmm-CAR-t-censored-left.stan", package = "SMNlmec")

        stan_obj <- rstan::stan(file= stan_file,
                              data = list(N_complete = N_com, N_cens = cens_N, ycen=ycen,
                                          n=m, l=l_set, q1 = q1_set,
                                          x=x, z=z, timevar = tt,
                                          y_complete = y_com, rho = rho_com,
                                          njvec =nj, cens_nj = cens_nj,
                                          ind = ind_set),
                              thin = thin_num, chains = chains_num,
                              iter = iter_num, warmup = iter_num*burn_percen,
                              seed = seed_set, control = list(adapt_delta=adapt_delta_set))

        SMNlmec_criteria <- criteria(cc = rho_com, nj = nj, y = y_com, x = x, z = z,
                              tt = tt, espac=5,
                              stanobj = stan_obj, distr="Student",
                              depstr = "CAR", cens.type="left",LI=NULL, LS=NULL)
      }

      if(direction == "right") {

        stan_file <- system.file("stan", "lmm-CAR-t-censored-right.stan", package = "SMNlmec")

        stan_obj <- rstan::stan(file= stan_file,
                         data = list(N_complete = N_com, N_cens = cens_N, ycen=ycen,
                                     n=m, l=l_set, q1 = q1_set,
                                     x=x, z=z, timevar = tt,
                                     y_complete = y_com, rho = rho_com,
                                     njvec =nj, cens_nj = cens_nj,
                                     ind = ind_set),
                         thin = thin_num, chains = chains_num,
                         iter = iter_num, warmup = iter_num*burn_percen,
                         seed = seed_set, control = list(adapt_delta=adapt_delta_set))

        SMNlmec_criteria <- criteria(cc = rho_com, nj = nj, y = y_com, x = x, z = z,
                                     tt = tt, espac=5,
                                     stanobj = stan_obj, distr="Student",
                                     depstr = "CAR", cens.type="right",LI=NULL, LS=NULL)
      }
    }

  }

  if(dist == "Slash") {
    if(struc == "UNC")   {
      if(direction == "left") {

        stan_file <- system.file("stan", "lmm-UNC-SL-censored-left.stan", package = "SMNlmec")

        stan_obj <- rstan::stan(file= stan_file,
                               data = list(N_complete = N_com, N_cens=cens_N, N_obs= obs_N,
                                           n = m, l = l_set,yobs=yobs,ycen=ycen,x=x,
                                           ind=ind_set,rho=rho_com),
                               thin = thin_num, chains = chains_num,
                               iter = iter_num, warmup = iter_num*burn_percen,
                               seed = seed_set, control = list(adapt_delta=adapt_delta_set))

        SMNlmec_criteria <- criteria(cc = rho_com, nj = nj, y = y_com, x = x, z = z,
                               tt = tt, espac=5,
                               stanobj = stan_obj, distr="Slash",
                               depstr = "UNC", cens.type="left",LI=NULL, LS=NULL)

      }

      if(direction == "right") {

        stan_file <- system.file("stan", "lmm-UNC-SL-censored-right.stan", package = "SMNlmec")

        stan_obj <- rstan::stan(file= stan_file,
                         data = list(N_complete = N_com, N_cens=cens_N, N_obs= obs_N,
                                     n = m, l = l_set,yobs=yobs,ycen=ycen,x=x,
                                     ind=ind_set,rho=rho_com),
                         thin = thin_num, chains = chains_num,
                         iter = iter_num, warmup = iter_num*burn_percen,
                         seed = seed_set, control = list(adapt_delta=adapt_delta_set))

        SMNlmec_criteria <- criteria(cc = rho_com, nj = nj, y = y_com, x = x, z = z,
                                     tt = tt, espac=5,
                                     stanobj = stan_obj, distr="Slash",
                                     depstr = "UNC", cens.type="right",LI=NULL, LS=NULL)

      }
    }

    if(struc == "DEC")   {
      if(direction == "left") {

        stan_file <- system.file("stan", "lmm-DEC-SL-censored-left.stan", package = "SMNlmec")

        stan_obj <- rstan::stan(file= stan_file,
                               data = list(N_complete = N_com, N_cens = cens_N, ycen=ycen,
                                           n=m, l=l_set, q1 = q1_set,
                                           x=x, z=z, timevar = tt,
                                           y_complete = y_com, rho = rho_com,
                                           njvec =nj, cens_nj = cens_nj,
                                           ind = ind_set),
                               thin = thin_num, chains = chains_num,
                               iter = iter_num, warmup = iter_num*burn_percen,
                               seed = seed_set, control = list(adapt_delta=adapt_delta_set))

        SMNlmec_criteria <- criteria(cc = rho_com, nj = nj, y = y_com, x = x, z = z,
                               tt = tt, espac=5,
                               stanobj = stan_obj, distr="Slash",
                               depstr = "DEC", cens.type="left",LI=NULL, LS=NULL)
      }

      if(direction == "right") {

        stan_file <- system.file("stan", "lmm-DEC-SL-censored-right.stan", package = "SMNlmec")

        stan_obj <- rstan::stan(file= stan_file,
                         data = list(N_complete = N_com, N_cens = cens_N, ycen=ycen,
                                     n=m, l=l_set, q1 = q1_set,
                                     x=x, z=z, timevar = tt,
                                     y_complete = y_com, rho = rho_com,
                                     njvec =nj, cens_nj = cens_nj,
                                     ind = ind_set),
                         thin = thin_num, chains = chains_num,
                         iter = iter_num, warmup = iter_num*burn_percen,
                         seed = seed_set, control = list(adapt_delta=adapt_delta_set))

        SMNlmec_criteria <- criteria(cc = rho_com, nj = nj, y = y_com, x = x, z = z,
                                     tt = tt, espac=5,
                                     stanobj = stan_obj, distr="Slash",
                                     depstr = "DEC", cens.type="right",LI=NULL, LS=NULL)
      }
    }

    if(struc == "CAR")   {
      if(direction == "left") {

        stan_file <- system.file("stan", "lmm-CAR-SL-censored-left.stan", package = "SMNlmec")

        stan_obj <- rstan::stan(file= stan_file,
                               data = list(N_complete = N_com, N_cens = cens_N, ycen=ycen,
                                           n=m, l=l_set, q1 = q1_set,
                                           x=x, z=z, timevar = tt,
                                           y_complete = y_com, rho = rho_com,
                                           njvec =nj, cens_nj = cens_nj,
                                           ind = ind_set),
                               thin = thin_num, chains = chains_num,
                               iter = iter_num, warmup = iter_num*burn_percen,
                               seed = seed_set, control = list(adapt_delta=adapt_delta_set))

        SMNlmec_criteria <- criteria(cc = rho_com, nj = nj, y = y_com, x = x, z = z,
                               tt = tt, espac=5,
                               stanobj = stan_obj, distr="Slash",
                               depstr = "CAR", cens.type="left",LI=NULL, LS=NULL)
      }

      if(direction == "right") {

        stan_file <- system.file("stan", "lmm-CAR-SL-censored-right.stan", package = "SMNlmec")

        stan_obj <- rstan::stan(file= stan_file,
                         data = list(N_complete = N_com, N_cens = cens_N, ycen=ycen,
                                     n=m, l=l_set, q1 = q1_set,
                                     x=x, z=z, timevar = tt,
                                     y_complete = y_com, rho = rho_com,
                                     njvec =nj, cens_nj = cens_nj,
                                     ind = ind_set),
                         thin = thin_num, chains = chains_num,
                         iter = iter_num, warmup = iter_num*burn_percen,
                         seed = seed_set, control = list(adapt_delta=adapt_delta_set))

        SMNlmec_criteria <- criteria(cc = rho_com, nj = nj, y = y_com, x = x, z = z,
                                     tt = tt, espac=5,
                                     stanobj = stan_obj, distr="Slash",
                                     depstr = "CAR", cens.type="right",LI=NULL, LS=NULL)
      }
    }
  }

  SMNlmec.est_object <- SMNlmecfit.creator(stan_obj, SMNlmec_criteria, dist, struc)

  return(SMNlmec.est_object)

}






