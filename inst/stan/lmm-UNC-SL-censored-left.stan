data {
  int<lower=0> N_complete;   // total observations
  int<lower=0> N_cens;   // total censuras
  int<lower=0> N_obs;   // total obs
  int<lower=0> n;   // total subjects
  int<lower=0> l;    // number covariates
  vector[N_obs] yobs;      // y-obs
  vector[N_cens] ycen;      // y-cens
  vector[N_complete] rho;     //indicador de censuras
  int<lower=1,upper=n> ind[N_complete];
  matrix[N_complete,l] x;  //covariates
}

parameters {
  vector[l]  beta;
  vector[n]  b;
  real<lower=0> sigma;
  real<lower=0> sigmab;
  vector<upper = ycen>[N_cens] y_cens;
  real<lower=1> nu;      //slash
  vector<lower = 0, upper = 1>[n] w; // slash and Student
 }

transformed parameters{
    real<lower=0> sigma2;
    real<lower=0> sigmab2;
    sigma2 = sigma^2;
    sigmab2 = sigmab^2;
}

model {
   
    for (r in 1:l){
  beta[r] ~ normal(0, 100); //regression parameters
     }
  
   sigma ~ student_t(4,0,5);
   sigmab ~ student_t(4,0,5);

   nu ~  student_t(4,0,5); 
 
    {
     int cens_node = 0;
   int obs_node = 0;
    vector[N_complete] mu;
    vector[N_complete] media;  
   vector[N_obs] mediay;
   vector[N_cens] mediac;
    vector[N_obs] sigmao;
    vector[N_cens] sigmac;
    mu = x*beta;
 		for(i in 1 : n) {
    w[i] ~ beta(nu,1);   //Slash
    b[i] ~ normal(0, inv_sqrt(w[i])*sigmab);
    	}
     
  
for(i in 1 : N_complete){
 media[i] = mu[i]+b[ind[i]];
  if(rho[i] == 0){
  obs_node = obs_node+1;
  mediay[obs_node] = media[i];
  sigmao[obs_node] = sigma*inv_sqrt(w[ind[i]]);
  yobs[obs_node]  ~ normal(mediay[obs_node], sigmao[obs_node]);
  }
  else{
  cens_node = cens_node+1;
  mediac[cens_node] = media[i];
  sigmac[cens_node] = sigma*inv_sqrt(w[ind[i]]);
  y_cens[cens_node] ~ normal(mediac[cens_node], sigmac[cens_node]);
  }
 }
      }
}
