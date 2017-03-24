// This Stan code was generated with the R package 'brms'. 
// We recommend generating the data with the 'make_standata' function. 
functions { 
} 
data { 
  int<lower=1> N;  // total number of observations 
  int Y[N];  // response variable 
  int<lower=1> K_eta;  // number of population-level effects 
  matrix[N, K_eta] X_eta;  // population-level design matrix 
  // data for group-specific effects of ID 1 
  int<lower=1> J_1[N]; 
  int<lower=1> N_1; 
  int<lower=1> M_1; 
  vector[N] Z_1_eta_1; 
  vector[N] Z_1_eta_2; 
  vector[N] Z_1_eta_3; 
  vector[N] Z_1_eta_4; 
  vector[N] Z_1_eta_5; 
  vector[N] Z_1_eta_6; 
  vector[N] Z_1_eta_7; 
  vector[N] Z_1_eta_8; 
  vector[N] Z_1_eta_9; 
  vector[N] Z_1_eta_10; 
  vector[N] Z_1_eta_11; 
  vector[N] Z_1_eta_12; 
  int<lower=1> NC_1; 
  // data for group-specific effects of ID 2 
  int<lower=1> J_2[N]; 
  int<lower=1> N_2; 
  int<lower=1> M_2; 
  vector[N] Z_2_eta_1; 
  vector[N] Z_2_eta_2; 
  vector[N] Z_2_eta_3; 
  vector[N] Z_2_eta_4; 
  vector[N] Z_2_eta_5; 
  vector[N] Z_2_eta_6; 
  vector[N] Z_2_eta_7; 
  vector[N] Z_2_eta_8; 
  vector[N] Z_2_eta_9; 
  vector[N] Z_2_eta_10; 
  vector[N] Z_2_eta_11; 
  vector[N] Z_2_eta_12; 
  int<lower=1> NC_2; 
  int prior_only;  // should the likelihood be ignored? 
} 
transformed data { 
} 
parameters { 
  vector<lower=-5,upper=5>[K_eta] b_eta;  // population-level effects 
  vector<lower=0>[M_1] sd_1;  // group-specific standard deviations 
  matrix[M_1, N_1] z_1;  // unscaled group-specific effects 
  // cholesky factor of correlation matrix 
  cholesky_factor_corr[M_1] L_1; 
  vector<lower=0>[M_2] sd_2;  // group-specific standard deviations 
  matrix[M_2, N_2] z_2;  // unscaled group-specific effects 
  // cholesky factor of correlation matrix 
  cholesky_factor_corr[M_2] L_2; 
  // parameters to store prior samples 
  real<lower=-5,upper=5> prior_b_eta; 
  real<lower=0> prior_sd_1; 
  real<lower=0> prior_sd_2; 
} 
transformed parameters { 
  // group-specific effects 
  matrix[N_1, M_1] r_1; 
  vector[N_1] r_1_eta_1; 
  vector[N_1] r_1_eta_2; 
  vector[N_1] r_1_eta_3; 
  vector[N_1] r_1_eta_4; 
  vector[N_1] r_1_eta_5; 
  vector[N_1] r_1_eta_6; 
  vector[N_1] r_1_eta_7; 
  vector[N_1] r_1_eta_8; 
  vector[N_1] r_1_eta_9; 
  vector[N_1] r_1_eta_10; 
  vector[N_1] r_1_eta_11; 
  vector[N_1] r_1_eta_12; 
  // group-specific effects 
  matrix[N_2, M_2] r_2; 
  vector[N_2] r_2_eta_1; 
  vector[N_2] r_2_eta_2; 
  vector[N_2] r_2_eta_3; 
  vector[N_2] r_2_eta_4; 
  vector[N_2] r_2_eta_5; 
  vector[N_2] r_2_eta_6; 
  vector[N_2] r_2_eta_7; 
  vector[N_2] r_2_eta_8; 
  vector[N_2] r_2_eta_9; 
  vector[N_2] r_2_eta_10; 
  vector[N_2] r_2_eta_11; 
  vector[N_2] r_2_eta_12; 
  r_1 = (diag_pre_multiply(sd_1, L_1) * z_1)'; 
  r_1_eta_1 = r_1[, 1];  
  r_1_eta_2 = r_1[, 2];  
  r_1_eta_3 = r_1[, 3];  
  r_1_eta_4 = r_1[, 4];  
  r_1_eta_5 = r_1[, 5];  
  r_1_eta_6 = r_1[, 6];  
  r_1_eta_7 = r_1[, 7];  
  r_1_eta_8 = r_1[, 8];  
  r_1_eta_9 = r_1[, 9];  
  r_1_eta_10 = r_1[, 10];  
  r_1_eta_11 = r_1[, 11];  
  r_1_eta_12 = r_1[, 12];  
  r_2 = (diag_pre_multiply(sd_2, L_2) * z_2)'; 
  r_2_eta_1 = r_2[, 1];  
  r_2_eta_2 = r_2[, 2];  
  r_2_eta_3 = r_2[, 3];  
  r_2_eta_4 = r_2[, 4];  
  r_2_eta_5 = r_2[, 5];  
  r_2_eta_6 = r_2[, 6];  
  r_2_eta_7 = r_2[, 7];  
  r_2_eta_8 = r_2[, 8];  
  r_2_eta_9 = r_2[, 9];  
  r_2_eta_10 = r_2[, 10];  
  r_2_eta_11 = r_2[, 11];  
  r_2_eta_12 = r_2[, 12];  
} 
model { 
  vector[N] eta_eta; 
  vector[N] eta; 
  eta_eta = X_eta * b_eta; 
  for (n in 1:N) { 
    eta_eta[n] = eta_eta[n] + r_1_eta_1[J_1[n]] * Z_1_eta_1[n] + r_1_eta_2[J_1[n]] * Z_1_eta_2[n] + r_1_eta_3[J_1[n]] * Z_1_eta_3[n] + r_1_eta_4[J_1[n]] * Z_1_eta_4[n] + r_1_eta_5[J_1[n]] * Z_1_eta_5[n] + r_1_eta_6[J_1[n]] * Z_1_eta_6[n] + r_1_eta_7[J_1[n]] * Z_1_eta_7[n] + r_1_eta_8[J_1[n]] * Z_1_eta_8[n] + r_1_eta_9[J_1[n]] * Z_1_eta_9[n] + r_1_eta_10[J_1[n]] * Z_1_eta_10[n] + r_1_eta_11[J_1[n]] * Z_1_eta_11[n] + r_1_eta_12[J_1[n]] * Z_1_eta_12[n] + r_2_eta_1[J_2[n]] * Z_2_eta_1[n] + r_2_eta_2[J_2[n]] * Z_2_eta_2[n] + r_2_eta_3[J_2[n]] * Z_2_eta_3[n] + r_2_eta_4[J_2[n]] * Z_2_eta_4[n] + r_2_eta_5[J_2[n]] * Z_2_eta_5[n] + r_2_eta_6[J_2[n]] * Z_2_eta_6[n] + r_2_eta_7[J_2[n]] * Z_2_eta_7[n] + r_2_eta_8[J_2[n]] * Z_2_eta_8[n] + r_2_eta_9[J_2[n]] * Z_2_eta_9[n] + r_2_eta_10[J_2[n]] * Z_2_eta_10[n] + r_2_eta_11[J_2[n]] * Z_2_eta_11[n] + r_2_eta_12[J_2[n]] * Z_2_eta_12[n]; 
    // compute non-linear predictor 
    eta[n] = 0.333333333 + 0.666666667 * inv_logit(eta_eta[n]); 
  } 
  // prior specifications 
  b_eta ~ normal(0, 2.5); 
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
  prior_b_eta ~ normal(0,2.5); 
  prior_sd_1 ~ cauchy(0,2.5); 
  prior_sd_2 ~ cauchy(0,2.5); 
} 
generated quantities { 
  corr_matrix[M_1] Cor_1; 
  vector<lower=-1,upper=1>[NC_1] cor_1; 
  corr_matrix[M_2] Cor_2; 
  vector<lower=-1,upper=1>[NC_2] cor_2; 
  real prior_cor_1; 
  real prior_cor_2; 
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
  prior_cor_1 = lkj_corr_rng(2,1)[1,2]; 
  prior_cor_2 = lkj_corr_rng(2,1)[1,2]; 
} 