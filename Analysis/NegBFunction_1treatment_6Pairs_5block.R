## fitting Neg Bin model without generation time, without the error term for lambda
## FEB 2021
## JINLIN

Write_Stan_NegB_1layer_6pairs_5blocks <- function(PreparedData1,Name,
                                     AlphaLowerBound = TRUE,
                                     AlphaSD= 1){
  
  # Function to Fit a Stan, using raw count data and a Negative-Binomial error distrubtion
  
  #############
  ### Core Start
  #########
  
  Base<- paste0('data {
  int<lower=0> N;
  int y[N];// Observed Count of emergences
}
parameters {
vector<lower=0>[5] mur;   // 5 growth terms for 5 species. Fixed Lower bound
vector',ifelse(AlphaLowerBound,'<lower=0>',''),'[17] A;  //17 A terms (as only some of the inter-pairs combination are measured). Fixed Lower bound
real<lower=0> phi;  // divergence
vector<lower=0>[3] r1;  // 3 growth terms of each of the 3 blocks for species 1
vector<lower=0>[5] r2;  // 5 growth terms of each of the 5 blocks for species 2 (PAL)
vector<lower=0>[5] r3;  // 5 growth terms of each of the 5 blocks for species 3 (PAN)
vector<lower=0>[3] r4;  // 3 growth terms of each of the 3 blocks for species 4
vector<lower=0>[3] r5;  // 3 growth terms of each of the 3 blocks for species 5
real<lower=0> sigmar;  // variance among blocks for r, assuming the same sigma for all species
}
model {
mur[1:5] ~ normal(20, 10);  // as from the previous PAN-PAL pair, their growth rates are around 20.
r1[1:3] ~ normal(mur[1], sigmar);
r2[1:5] ~ normal(mur[2], sigmar);
r3[1:5] ~ normal(mur[3], sigmar);
r4[1:3] ~ normal(mur[4], sigmar);
r5[1:3] ~ normal(mur[5], sigmar);
sigmar ~ gamma(2,0.1);  // recommended by the rstan to avoid prior problem, but I do not know how it helps
A[1:17] ~ normal(0, ',AlphaSD,');
phi ~ cauchy(0,10);  // prior of phi
') 
  
  GQ<-'generated quantities {
  vector[N] log_lik; 
  vector[N] y_sim;\n'          # Start of generated quantitities section
  
  
  
  d=PreparedData1
  n=nrow(d)
  
  r_sel <- d$r.ind
  blockN <- d$blockN
  
  R <-  data.frame(N0_i = d$FocalSp_Den,
                   N0_j = d$OtherSp_Den,
                   Aii_2 = d$Aii.ind,  # index of Aii
                   Aij_2 = d$Aij.ind)  # index of Aij
  
  
  for(i in 1:n){
    GrowthCore = paste0( 'r', r_sel[i], '[',blockN[i],']/(1+',R$N0_i[i],'*A[',R$Aii_2[i],']')  # NB needs brackets of denominator closing!
    if(R$Aij_2[i] >0){
      GrowthCore <- paste0(GrowthCore, '+',R$N0_j[i],'*A[',R$Aij_2[i],  ']')    # Add on other species  
    }
    Base <- paste0(Base, 'y[',i,'] ~  neg_binomial_2(',R$N0_i[i],'*(',GrowthCore,') )',', phi);\n')
    GQ <- paste0(GQ, 'y_sim[',i,'] = neg_binomial_2_rng( (',R$N0_i[i],'*( ', GrowthCore, ')  )','), phi);\n',
                 'log_lik[',i,'] = neg_binomial_2_lpmf(y[',i,']|',R$N0_i[i],'*(', GrowthCore,') )',', phi);\n')
  }
  
  Both <- paste0(Base,'\n} \n', GQ, '\n}')
  writeLines(Both, paste0('StanModels/BuiltModel_NegB_',Name,'.stan'))
  return(paste0('Function Finished. Output: <BuiltModel_NegB_',Name,'.stan>'))
}