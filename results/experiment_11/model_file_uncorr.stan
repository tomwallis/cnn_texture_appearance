// This Stan code was generated with the R package 'brms'. 
// We recommend generating the data with the 'make_standata' function. 
functions { 
} 
data { 
  int<lower=1> N;  // total number of observations 
  int Y[N];  // response variable 
  int<lower=1> K;  // number of population-level effects 
  matrix[N, K] X;  // centered population-level design matrix 
  vector[K] X_means;  // column means of X before centering 
  // data for group-specific effects of ID 1 
  int<lower=1> J_1[N]; 
  int<lower=1> N_1; 
  int<lower=1> M_1; 
  vector[N] Z_1_1;  
  vector[N] Z_1_2;  
  vector[N] Z_1_3;  
  vector[N] Z_1_4;  
  vector[N] Z_1_5;  
  vector[N] Z_1_6;  
  vector[N] Z_1_7;  
  vector[N] Z_1_8;  
  vector[N] Z_1_9;  
  vector[N] Z_1_10;  
  vector[N] Z_1_11;  
  vector[N] Z_1_12;  
  // data for group-specific effects of ID 2 
  int<lower=1> J_2[N]; 
  int<lower=1> N_2; 
  int<lower=1> M_2; 
  vector[N] Z_2_1;  
  vector[N] Z_2_2;  
  vector[N] Z_2_3;  
  vector[N] Z_2_4;  
  vector[N] Z_2_5;  
  vector[N] Z_2_6;  
  vector[N] Z_2_7;  
  vector[N] Z_2_8;  
  vector[N] Z_2_9;  
  vector[N] Z_2_10;  
  vector[N] Z_2_11;  
  vector[N] Z_2_12;  
  int prior_only;  // should the likelihood be ignored? 
} 
transformed data { 
} 
parameters { 
  vector<lower=-10,upper=10>[K] b;  // population-level effects 
  real temp_Intercept;  // temporary Intercept 
  vector<lower=0>[M_1] sd_1;  // group-specific standard deviations 
  vector[N_1] z_1[M_1];  // unscaled group-specific effects 
  vector<lower=0>[M_2] sd_2;  // group-specific standard deviations 
  vector[N_2] z_2[M_2];  // unscaled group-specific effects 
  // parameters to store prior samples 
  real<lower=-10,upper=10> prior_b; 
  real<lower=0> prior_sd_1; 
  real<lower=0> prior_sd_2; 
} 
transformed parameters { 
  // group-specific effects 
  vector[N_1] r_1_1; 
  vector[N_1] r_1_2; 
  vector[N_1] r_1_3; 
  vector[N_1] r_1_4; 
  vector[N_1] r_1_5; 
  vector[N_1] r_1_6; 
  vector[N_1] r_1_7; 
  vector[N_1] r_1_8; 
  vector[N_1] r_1_9; 
  vector[N_1] r_1_10; 
  vector[N_1] r_1_11; 
  vector[N_1] r_1_12; 
  // group-specific effects 
  vector[N_2] r_2_1; 
  vector[N_2] r_2_2; 
  vector[N_2] r_2_3; 
  vector[N_2] r_2_4; 
  vector[N_2] r_2_5; 
  vector[N_2] r_2_6; 
  vector[N_2] r_2_7; 
  vector[N_2] r_2_8; 
  vector[N_2] r_2_9; 
  vector[N_2] r_2_10; 
  vector[N_2] r_2_11; 
  vector[N_2] r_2_12; 
  r_1_1 = sd_1[1] * (z_1[1]); 
  r_1_2 = sd_1[2] * (z_1[2]); 
  r_1_3 = sd_1[3] * (z_1[3]); 
  r_1_4 = sd_1[4] * (z_1[4]); 
  r_1_5 = sd_1[5] * (z_1[5]); 
  r_1_6 = sd_1[6] * (z_1[6]); 
  r_1_7 = sd_1[7] * (z_1[7]); 
  r_1_8 = sd_1[8] * (z_1[8]); 
  r_1_9 = sd_1[9] * (z_1[9]); 
  r_1_10 = sd_1[10] * (z_1[10]); 
  r_1_11 = sd_1[11] * (z_1[11]); 
  r_1_12 = sd_1[12] * (z_1[12]); 
  r_2_1 = sd_2[1] * (z_2[1]); 
  r_2_2 = sd_2[2] * (z_2[2]); 
  r_2_3 = sd_2[3] * (z_2[3]); 
  r_2_4 = sd_2[4] * (z_2[4]); 
  r_2_5 = sd_2[5] * (z_2[5]); 
  r_2_6 = sd_2[6] * (z_2[6]); 
  r_2_7 = sd_2[7] * (z_2[7]); 
  r_2_8 = sd_2[8] * (z_2[8]); 
  r_2_9 = sd_2[9] * (z_2[9]); 
  r_2_10 = sd_2[10] * (z_2[10]); 
  r_2_11 = sd_2[11] * (z_2[11]); 
  r_2_12 = sd_2[12] * (z_2[12]); 
} 
model { 
  vector[N] eta; 
  eta = X * b + temp_Intercept; 
  for (n in 1:N) { 
    eta[n] = eta[n] + r_1_1[J_1[n]] * Z_1_1[n] + r_1_2[J_1[n]] * Z_1_2[n] + r_1_3[J_1[n]] * Z_1_3[n] + r_1_4[J_1[n]] * Z_1_4[n] + r_1_5[J_1[n]] * Z_1_5[n] + r_1_6[J_1[n]] * Z_1_6[n] + r_1_7[J_1[n]] * Z_1_7[n] + r_1_8[J_1[n]] * Z_1_8[n] + r_1_9[J_1[n]] * Z_1_9[n] + r_1_10[J_1[n]] * Z_1_10[n] + r_1_11[J_1[n]] * Z_1_11[n] + r_1_12[J_1[n]] * Z_1_12[n] + r_2_1[J_2[n]] * Z_2_1[n] + r_2_2[J_2[n]] * Z_2_2[n] + r_2_3[J_2[n]] * Z_2_3[n] + r_2_4[J_2[n]] * Z_2_4[n] + r_2_5[J_2[n]] * Z_2_5[n] + r_2_6[J_2[n]] * Z_2_6[n] + r_2_7[J_2[n]] * Z_2_7[n] + r_2_8[J_2[n]] * Z_2_8[n] + r_2_9[J_2[n]] * Z_2_9[n] + r_2_10[J_2[n]] * Z_2_10[n] + r_2_11[J_2[n]] * Z_2_11[n] + r_2_12[J_2[n]] * Z_2_12[n]; 
  } 
  // prior specifications 
  b ~ student_t(3,0,2.5); 
  sd_1 ~ student_t(3,0,2.5); 
  z_1[1] ~ normal(0, 1); 
  z_1[2] ~ normal(0, 1); 
  z_1[3] ~ normal(0, 1); 
  z_1[4] ~ normal(0, 1); 
  z_1[5] ~ normal(0, 1); 
  z_1[6] ~ normal(0, 1); 
  z_1[7] ~ normal(0, 1); 
  z_1[8] ~ normal(0, 1); 
  z_1[9] ~ normal(0, 1); 
  z_1[10] ~ normal(0, 1); 
  z_1[11] ~ normal(0, 1); 
  z_1[12] ~ normal(0, 1); 
  sd_2 ~ student_t(3,0,2.5); 
  z_2[1] ~ normal(0, 1); 
  z_2[2] ~ normal(0, 1); 
  z_2[3] ~ normal(0, 1); 
  z_2[4] ~ normal(0, 1); 
  z_2[5] ~ normal(0, 1); 
  z_2[6] ~ normal(0, 1); 
  z_2[7] ~ normal(0, 1); 
  z_2[8] ~ normal(0, 1); 
  z_2[9] ~ normal(0, 1); 
  z_2[10] ~ normal(0, 1); 
  z_2[11] ~ normal(0, 1); 
  z_2[12] ~ normal(0, 1); 
  // likelihood contribution 
  if (!prior_only) { 
    Y ~ bernoulli_logit(eta); 
  } 
  // additionally draw samples from priors 
  prior_b ~ student_t(3,0,2.5); 
  prior_sd_1 ~ student_t(3,0,2.5); 
  prior_sd_2 ~ student_t(3,0,2.5); 
} 
generated quantities { 
  real b_Intercept;  // population-level intercept 
  b_Intercept = temp_Intercept - dot_product(X_means, b); 
} 