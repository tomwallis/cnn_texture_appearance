// generated with brms 1.4.0
functions { 
} 
data { 
  int<lower=1> N;  // total number of observations 
  int Y[N];  // response variable 
  int<lower=1> K;  // number of population-level effects 
  matrix[N, K] X;  // population-level design matrix 
  // data for group-level effects of ID 1 
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
  int<lower=1> NC_1; 
  // data for group-level effects of ID 2 
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
  int<lower=1> NC_2; 
  int prior_only;  // should the likelihood be ignored? 
} 
transformed data { 
  int Kc; 
  matrix[N, K - 1] Xc;  // centered version of X 
  vector[K - 1] means_X;  // column means of X before centering 
  Kc = K - 1;  // the intercept is removed from the design matrix 
  for (i in 2:K) { 
    means_X[i - 1] = mean(X[, i]); 
    Xc[, i - 1] = X[, i] - means_X[i - 1]; 
  } 
} 
parameters { 
  vector<lower=-5,upper=5>[Kc] b;  // population-level effects 
  real temp_Intercept;  // temporary intercept 
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations 
  matrix[M_1, N_1] z_1;  // unscaled group-level effects 
  // cholesky factor of correlation matrix 
  cholesky_factor_corr[M_1] L_1; 
  vector<lower=0>[M_2] sd_2;  // group-level standard deviations 
  matrix[M_2, N_2] z_2;  // unscaled group-level effects 
  // cholesky factor of correlation matrix 
  cholesky_factor_corr[M_2] L_2; 
  // parameters to store prior samples
  real<lower=-5,upper=5> prior_b;
  real<lower=0> prior_sd_1;
  real<lower=0> prior_sd_2;
} 
transformed parameters { 
  // group-level effects 
  matrix[N_1, M_1] r_1; 
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
  // group-level effects 
  matrix[N_2, M_2] r_2; 
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
  r_1 = (diag_pre_multiply(sd_1, L_1) * z_1)'; 
  r_1_1 = r_1[, 1];  
  r_1_2 = r_1[, 2];  
  r_1_3 = r_1[, 3];  
  r_1_4 = r_1[, 4];  
  r_1_5 = r_1[, 5];  
  r_1_6 = r_1[, 6];  
  r_1_7 = r_1[, 7];  
  r_1_8 = r_1[, 8];  
  r_1_9 = r_1[, 9];  
  r_1_10 = r_1[, 10];  
  r_1_11 = r_1[, 11];  
  r_1_12 = r_1[, 12];  
  r_2 = (diag_pre_multiply(sd_2, L_2) * z_2)'; 
  r_2_1 = r_2[, 1];  
  r_2_2 = r_2[, 2];  
  r_2_3 = r_2[, 3];  
  r_2_4 = r_2[, 4];  
  r_2_5 = r_2[, 5];  
  r_2_6 = r_2[, 6];  
  r_2_7 = r_2[, 7];  
  r_2_8 = r_2[, 8];  
  r_2_9 = r_2[, 9];  
  r_2_10 = r_2[, 10];  
  r_2_11 = r_2[, 11];  
  r_2_12 = r_2[, 12];  
} 
model { 
  vector[N] eta; 
  eta = Xc * b + temp_Intercept; 
  for (n in 1:N) { 
    eta[n] = eta[n] + (r_1_1[J_1[n]]) * Z_1_1[n] + (r_1_2[J_1[n]]) * Z_1_2[n] + (r_1_3[J_1[n]]) * Z_1_3[n] + (r_1_4[J_1[n]]) * Z_1_4[n] + (r_1_5[J_1[n]]) * Z_1_5[n] + (r_1_6[J_1[n]]) * Z_1_6[n] + (r_1_7[J_1[n]]) * Z_1_7[n] + (r_1_8[J_1[n]]) * Z_1_8[n] + (r_1_9[J_1[n]]) * Z_1_9[n] + (r_1_10[J_1[n]]) * Z_1_10[n] + (r_1_11[J_1[n]]) * Z_1_11[n] + (r_1_12[J_1[n]]) * Z_1_12[n] + (r_2_1[J_2[n]]) * Z_2_1[n] + (r_2_2[J_2[n]]) * Z_2_2[n] + (r_2_3[J_2[n]]) * Z_2_3[n] + (r_2_4[J_2[n]]) * Z_2_4[n] + (r_2_5[J_2[n]]) * Z_2_5[n] + (r_2_6[J_2[n]]) * Z_2_6[n] + (r_2_7[J_2[n]]) * Z_2_7[n] + (r_2_8[J_2[n]]) * Z_2_8[n] + (r_2_9[J_2[n]]) * Z_2_9[n] + (r_2_10[J_2[n]]) * Z_2_10[n] + (r_2_11[J_2[n]]) * Z_2_11[n] + (r_2_12[J_2[n]]) * Z_2_12[n]; 
  } 
  // prior specifications 
  b ~ normal(0, 2.5); 
  sd_1 ~ cauchy(0, 2.5); 
  L_1 ~ lkj_corr_cholesky(1); 
  to_vector(z_1) ~ normal(0, 1); 
  sd_2 ~ cauchy(0, 2.5); 
  L_2 ~ lkj_corr_cholesky(1); 
  to_vector(z_2) ~ normal(0, 1); 
  // likelihood contribution 
  if (!prior_only) { 
    Y ~ bernoulli_logit(eta); 
  } 
  // additionally draw samples from priors
  prior_b ~ normal(0,2.5);
  prior_sd_1 ~ cauchy(0,2.5);
  prior_sd_2 ~ cauchy(0,2.5);
} 
generated quantities { 
  real b_Intercept;  // population-level intercept 
  corr_matrix[M_1] Cor_1; 
  vector<lower=-1,upper=1>[NC_1] cor_1; 
  corr_matrix[M_2] Cor_2; 
  vector<lower=-1,upper=1>[NC_2] cor_2; 
  real prior_cor_1; 
  real prior_cor_2; 
  b_Intercept = temp_Intercept - dot_product(means_X, b); 
  // take only relevant parts of correlation matrix 
  Cor_1 = multiply_lower_tri_self_transpose(L_1); 
  cor_1[1] = Cor_1[1,2]; 
  cor_1[2] = Cor_1[1,3]; 
  cor_1[3] = Cor_1[2,3]; 
  cor_1[4] = Cor_1[1,4]; 
  cor_1[5] = Cor_1[2,4]; 
  cor_1[6] = Cor_1[3,4]; 
  cor_1[7] = Cor_1[1,5]; 
  cor_1[8] = Cor_1[2,5]; 
  cor_1[9] = Cor_1[3,5]; 
  cor_1[10] = Cor_1[4,5]; 
  cor_1[11] = Cor_1[1,6]; 
  cor_1[12] = Cor_1[2,6]; 
  cor_1[13] = Cor_1[3,6]; 
  cor_1[14] = Cor_1[4,6]; 
  cor_1[15] = Cor_1[5,6]; 
  cor_1[16] = Cor_1[1,7]; 
  cor_1[17] = Cor_1[2,7]; 
  cor_1[18] = Cor_1[3,7]; 
  cor_1[19] = Cor_1[4,7]; 
  cor_1[20] = Cor_1[5,7]; 
  cor_1[21] = Cor_1[6,7]; 
  cor_1[22] = Cor_1[1,8]; 
  cor_1[23] = Cor_1[2,8]; 
  cor_1[24] = Cor_1[3,8]; 
  cor_1[25] = Cor_1[4,8]; 
  cor_1[26] = Cor_1[5,8]; 
  cor_1[27] = Cor_1[6,8]; 
  cor_1[28] = Cor_1[7,8]; 
  cor_1[29] = Cor_1[1,9]; 
  cor_1[30] = Cor_1[2,9]; 
  cor_1[31] = Cor_1[3,9]; 
  cor_1[32] = Cor_1[4,9]; 
  cor_1[33] = Cor_1[5,9]; 
  cor_1[34] = Cor_1[6,9]; 
  cor_1[35] = Cor_1[7,9]; 
  cor_1[36] = Cor_1[8,9]; 
  cor_1[37] = Cor_1[1,10]; 
  cor_1[38] = Cor_1[2,10]; 
  cor_1[39] = Cor_1[3,10]; 
  cor_1[40] = Cor_1[4,10]; 
  cor_1[41] = Cor_1[5,10]; 
  cor_1[42] = Cor_1[6,10]; 
  cor_1[43] = Cor_1[7,10]; 
  cor_1[44] = Cor_1[8,10]; 
  cor_1[45] = Cor_1[9,10]; 
  cor_1[46] = Cor_1[1,11]; 
  cor_1[47] = Cor_1[2,11]; 
  cor_1[48] = Cor_1[3,11]; 
  cor_1[49] = Cor_1[4,11]; 
  cor_1[50] = Cor_1[5,11]; 
  cor_1[51] = Cor_1[6,11]; 
  cor_1[52] = Cor_1[7,11]; 
  cor_1[53] = Cor_1[8,11]; 
  cor_1[54] = Cor_1[9,11]; 
  cor_1[55] = Cor_1[10,11]; 
  cor_1[56] = Cor_1[1,12]; 
  cor_1[57] = Cor_1[2,12]; 
  cor_1[58] = Cor_1[3,12]; 
  cor_1[59] = Cor_1[4,12]; 
  cor_1[60] = Cor_1[5,12]; 
  cor_1[61] = Cor_1[6,12]; 
  cor_1[62] = Cor_1[7,12]; 
  cor_1[63] = Cor_1[8,12]; 
  cor_1[64] = Cor_1[9,12]; 
  cor_1[65] = Cor_1[10,12]; 
  cor_1[66] = Cor_1[11,12]; 
  // take only relevant parts of correlation matrix 
  Cor_2 = multiply_lower_tri_self_transpose(L_2); 
  cor_2[1] = Cor_2[1,2]; 
  cor_2[2] = Cor_2[1,3]; 
  cor_2[3] = Cor_2[2,3]; 
  cor_2[4] = Cor_2[1,4]; 
  cor_2[5] = Cor_2[2,4]; 
  cor_2[6] = Cor_2[3,4]; 
  cor_2[7] = Cor_2[1,5]; 
  cor_2[8] = Cor_2[2,5]; 
  cor_2[9] = Cor_2[3,5]; 
  cor_2[10] = Cor_2[4,5]; 
  cor_2[11] = Cor_2[1,6]; 
  cor_2[12] = Cor_2[2,6]; 
  cor_2[13] = Cor_2[3,6]; 
  cor_2[14] = Cor_2[4,6]; 
  cor_2[15] = Cor_2[5,6]; 
  cor_2[16] = Cor_2[1,7]; 
  cor_2[17] = Cor_2[2,7]; 
  cor_2[18] = Cor_2[3,7]; 
  cor_2[19] = Cor_2[4,7]; 
  cor_2[20] = Cor_2[5,7]; 
  cor_2[21] = Cor_2[6,7]; 
  cor_2[22] = Cor_2[1,8]; 
  cor_2[23] = Cor_2[2,8]; 
  cor_2[24] = Cor_2[3,8]; 
  cor_2[25] = Cor_2[4,8]; 
  cor_2[26] = Cor_2[5,8]; 
  cor_2[27] = Cor_2[6,8]; 
  cor_2[28] = Cor_2[7,8]; 
  cor_2[29] = Cor_2[1,9]; 
  cor_2[30] = Cor_2[2,9]; 
  cor_2[31] = Cor_2[3,9]; 
  cor_2[32] = Cor_2[4,9]; 
  cor_2[33] = Cor_2[5,9]; 
  cor_2[34] = Cor_2[6,9]; 
  cor_2[35] = Cor_2[7,9]; 
  cor_2[36] = Cor_2[8,9]; 
  cor_2[37] = Cor_2[1,10]; 
  cor_2[38] = Cor_2[2,10]; 
  cor_2[39] = Cor_2[3,10]; 
  cor_2[40] = Cor_2[4,10]; 
  cor_2[41] = Cor_2[5,10]; 
  cor_2[42] = Cor_2[6,10]; 
  cor_2[43] = Cor_2[7,10]; 
  cor_2[44] = Cor_2[8,10]; 
  cor_2[45] = Cor_2[9,10]; 
  cor_2[46] = Cor_2[1,11]; 
  cor_2[47] = Cor_2[2,11]; 
  cor_2[48] = Cor_2[3,11]; 
  cor_2[49] = Cor_2[4,11]; 
  cor_2[50] = Cor_2[5,11]; 
  cor_2[51] = Cor_2[6,11]; 
  cor_2[52] = Cor_2[7,11]; 
  cor_2[53] = Cor_2[8,11]; 
  cor_2[54] = Cor_2[9,11]; 
  cor_2[55] = Cor_2[10,11]; 
  cor_2[56] = Cor_2[1,12]; 
  cor_2[57] = Cor_2[2,12]; 
  cor_2[58] = Cor_2[3,12]; 
  cor_2[59] = Cor_2[4,12]; 
  cor_2[60] = Cor_2[5,12]; 
  cor_2[61] = Cor_2[6,12]; 
  cor_2[62] = Cor_2[7,12]; 
  cor_2[63] = Cor_2[8,12]; 
  cor_2[64] = Cor_2[9,12]; 
  cor_2[65] = Cor_2[10,12]; 
  cor_2[66] = Cor_2[11,12]; 
  // additionally draw samples from priors 
  prior_cor_1 = lkj_corr_rng(2,1)[1, 2]; 
  prior_cor_2 = lkj_corr_rng(2,1)[1, 2]; 
} 