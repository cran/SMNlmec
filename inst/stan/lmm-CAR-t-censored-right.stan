functions{
	matrix vector_to_symmat(vector x, int n) {
	  matrix[n, n] m;
	  int k;
	  k = 1;
	  for (j in 1:n) {
		for (i in 1:j) {
		  m[i, j] = x[k];
		  if (i != j) {
			m[j, i] = m[i, j];
		  }
		  k = k + 1;
		}
	  }
	  return m;
	}

	matrix CAR_matrix(real phi1_set, int nj, vector tt) {
     matrix[nj, nj] mat;
     real lag;

     for(i in 1:nj){
       for(j in 1:nj){
         lag = fabs(tt[i]-tt[j]);
         mat[i,j] = pow(phi1_set, pow(lag,1));
       }
     }


     return mat;
   }

     vector nj_cens_mu(vector v_nj,vector c_nj,int nj_obs,int nj_cens, int nj_n, vector Mu_set, matrix Sigma_set){
    vector[nj_cens] mu_hat;

    matrix[nj_obs,nj_obs] var22;
    matrix[nj_cens,nj_obs] var12;
    vector[nj_obs] y_nj_obs;
    vector[nj_obs] mu_obs;
    vector[nj_cens] mu_cens;
    int cens_ind[nj_cens];
    int obs_ind[nj_obs];
    int l;
    int m;

    l = 1;
    while(l <= nj_cens){
      for(j in 1:nj_n){
        if(c_nj[j] == 1) {
          cens_ind[l] = j;
          l = l + 1;
        }
      }
    }

    m = 1;
    while(m <= nj_obs){
      for(j in 1:nj_n){
        if(c_nj[j] == 0) {
          obs_ind[m] = j;
          m = m + 1;
        }
      }
    }

    for(i in 1:nj_obs){
      y_nj_obs[i] = v_nj[obs_ind[i]];
      mu_obs[i] = Mu_set[obs_ind[i]];
    }

    for(i in 1:nj_cens){
      mu_cens[i] = Mu_set[cens_ind[i]];
    }

    for(i in 1:nj_obs){
      for(j in 1:nj_cens){
        var12[j,i] = Sigma_set[cens_ind[j],obs_ind[i]];
      }
    }

    for(i in 1:nj_obs){
      for(j in 1:nj_obs){
        var22[i,j] = Sigma_set[obs_ind[i],obs_ind[j]];
      }
    }

    mu_hat = mu_cens + var12*inverse(var22)*(y_nj_obs - mu_obs);

    return mu_hat;
  }

  matrix nj_cens_sigma(vector c_nj, int nj_obs,int nj_cens, int nj_n, matrix Sigma_set){
    matrix[nj_cens,nj_cens] sigma_hat;
    matrix[nj_obs,nj_obs] var22;
    matrix[nj_cens,nj_cens] var11;
    matrix[nj_cens,nj_obs] var12;
    matrix[nj_obs,nj_cens] var21;
    int cens_ind[nj_cens];
    int obs_ind[nj_obs];
    int l;
    int m;

    l = 1;
    while(l <= nj_cens){
      for(j in 1:nj_n){
        if(c_nj[j] == 1 ) {
          cens_ind[l] = j;
          l = l + 1;
        }
      }
    }

    m = 1;
    while(m <= nj_obs){
      for(j in 1:nj_n){
        if(c_nj[j] == 0 ) {
          obs_ind[m] = j;
          m = m + 1;
        }
      }
    }



    for(i in 1:nj_cens){
      for(j in 1:nj_cens){
        var11[i,j] = Sigma_set[cens_ind[i],cens_ind[j]];
      }
    }

    for(i in 1:nj_obs){
      for(j in 1:nj_obs){
        var22[i,j] = Sigma_set[obs_ind[i],obs_ind[j]];
      }
    }



    for(i in 1:nj_obs){
      for(j in 1:nj_cens){
        var12[j,i] = Sigma_set[cens_ind[j],obs_ind[i]];
        var21[i,j] = Sigma_set[obs_ind[i],cens_ind[j]];
      }
    }


    sigma_hat = var11 - var12*inverse(var22)*var21;

    return sigma_hat;
  }

  vector Mu_subset(int njloc, int nj_n, vector Mu_set){
    vector[nj_n] Mu_temp;

      for(k in 1:nj_n){
        Mu_temp[k] = Mu_set[njloc + k];
      }

      return Mu_temp;
  }

  matrix Sigma_subset(int njloc, int nj_n, matrix Sigma_set){
      matrix[nj_n,nj_n] Sigma_temp;

       for(k in 1:nj_n){
          for(l in 1:nj_n){
            Sigma_temp[k,l] = Sigma_set[njloc + k, njloc + l];
          }
      }
      return Sigma_temp;
    }


  vector Mu_yhat(vector Mu_c, int nj, int njloc){
    vector[nj] Mu_temp_hat;

    Mu_temp_hat = Mu_c[njloc + 1: njloc + nj];

    return Mu_temp_hat;
  }
}

data {
     int<lower=1> N_complete;
     int<lower=0> N_cens;

     int<lower=1> n; // total number of subjects
     int<lower=1> l; // column of x
     int<lower=1> q1; // column of z

     matrix[N_complete,l] x;
     matrix[N_complete,q1] z;
     vector[N_complete] timevar;
     vector[N_cens] ycen;

     vector[N_complete] y_complete;
     vector[N_complete] rho; // censored indicator
     int njvec[n]; // number of observation in each subject
     int cens_nj[n]; // censored number in each subject

     int<lower=1,upper=n> ind[N_complete];
}

transformed data{
  int N_obs;

  N_obs = N_complete - N_cens;
}

parameters {
      vector[l] beta;

		  cholesky_factor_corr[q1] Lcorr;// cholesky factor (L_u matrix for D1R)
		  vector<lower=0>[q1] ddsqrt;

      real<lower=0> sigmae;
		  real<lower=0.0001, upper=.9999> phi1;

      matrix[q1, n] etavec;
      vector<lower=0>[n] uvec;
		  real<lower=2> nu;
      vector<lower = ycen>[N_cens] y_cens;

}


transformed parameters {
  matrix[q1, n] bvec;
  cov_matrix[q1] D1; // VCV matrix
  real<lower=0> sigma2;

  sigma2 = sigmae * sigmae;

  D1 = quad_form_diag(multiply_lower_tri_self_transpose(Lcorr), ddsqrt); // quad_form_diag: diag_matrix(ddsqrt) * d1R * diag_matrix(ddsqrt)
  {
    matrix[q1,q1] dL = diag_pre_multiply(ddsqrt, Lcorr);
	for (j in 1:n){
	  bvec[,j] = dL*inv_sqrt(uvec[j]) * etavec[,j];
    }
  }
}


model {

  beta ~ normal(0,100);
  sigmae ~ student_t(4,0,5);

  ddsqrt ~ student_t(4,0,5);
  Lcorr ~ lkj_corr_cholesky(2.0);
  to_vector(etavec) ~ std_normal();

  phi1 ~ beta(1,1);

  nu ~ student_t(4,0,5);//cauchy(0, 2.5)

  uvec ~ gamma(nu/2,nu/2);

  {
	vector[N_complete] yhat;

	int njloc = 0;
  int cens_node = 0;
  int count_nj = 0;

  int obs_temp;
  real sigma_sqrt;

	for (i in 1:N_complete){
	  yhat[i] = x[i]*beta+row(z,i)*col(bvec,ind[i]);
	}


	for (j in 1:n){
		{
		  matrix[njvec[j],njvec[j]] Sigma_Ri = CAR_matrix(phi1,njvec[j],timevar[(njloc+1):(njloc + njvec[j])]);
		  matrix[njvec[j],njvec[j]] Sigma = sigma2 *  inv(uvec[j]) * Sigma_Ri;
		  vector[njvec[j]] Mu = yhat[njloc + 1:njloc + njvec[j]];

      if(cens_nj[j] == 0){
        njloc = njloc + njvec[j];
      }

      else if(cens_nj[j] == 1){

          obs_temp = njvec[j] - cens_nj[j];
          sigma_sqrt = sqrt(nj_cens_sigma(rho[(njloc + 1):(njloc + njvec[j])],obs_temp,cens_nj[j],njvec[j], Sigma)[1,1]);
          y_cens[cens_node + 1] ~ normal(nj_cens_mu(y_complete[(njloc + 1):(njloc + njvec[j])],rho[(njloc + 1):(njloc + njvec[j])],obs_temp,cens_nj[j],njvec[j],Mu,Sigma)[1],sigma_sqrt);
          cens_node = cens_node + cens_nj[j];
          njloc = njloc + njvec[j];
        }

      else if(cens_nj[j] == njvec[j]){

       y_cens[(cens_node + 1):(cens_node + njvec[j])] ~ multi_normal(Mu,Sigma);
       cens_node = cens_node + cens_nj[j];
       njloc = njloc + njvec[j];
      }

      else {

        obs_temp = njvec[j] - cens_nj[j];
        y_cens[(cens_node + 1):(cens_node + cens_nj[j])] ~ multi_normal(nj_cens_mu(y_complete[(njloc + 1):(njloc + njvec[j])],rho[(njloc + 1):(njloc + njvec[j])],obs_temp,cens_nj[j],njvec[j],Mu,Sigma),nj_cens_sigma(rho[(njloc + 1):(njloc + njvec[j])],obs_temp,cens_nj[j],njvec[j], Sigma));
        cens_node = cens_node + cens_nj[j];
        njloc = njloc + njvec[j];
        }

      y_complete[(count_nj+1):(count_nj+njvec[j])] ~ multi_normal(yhat[(count_nj+1):(count_nj+njvec[j])], Sigma);
		  count_nj += njvec[j];
      }
	  }
	}
}
