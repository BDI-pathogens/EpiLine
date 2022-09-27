functions {
   vector gamma_lpdf_v( vector x, vector alpha, vector beta )
  {
    return( alpha .* log( beta ) - lgamma( alpha ) + ( alpha-1) .* log(x ) - beta .* x );
  }
  
  int ceil_integer(real x, int min_val, int max_val)
  {
    int val = min_val;
    if( x > max_val )
      reject( "x must be less than max_val" );

    while( val < x )
      val = val + 1;

    return val;
  }
}

data {
  int<lower=1> t_rep;           // total days of reporting
  int<lower=1> t_symptom_pre;   // max time of reporting after symptoms
  int reported[t_rep];          // number of report cases by day
  int n_ll;
  int ll_report[n_ll];
  int ll_symptoms[n_ll];
  int<lower=0> ll_N[n_ll];
  real prior_alpha_min;
  real prior_beta_min;
  real prior_alpha_max;
  real prior_beta_max;
  real prior_alpha_gp_sd_max;
  real prior_beta_gp_sd_max;
  real prior_log_symptoms0_min;
  real prior_log_symptoms0_max;
  real prior_r_0_min;
  real prior_r_0_max;
  real prior_r_gp_sd_max;
  real prior_phi_od_max;
  int<lower=1,upper=t_rep + t_symptom_pre > hyper_gp_period_r;
  int<lower=1,upper=t_rep + t_symptom_pre > hyper_gp_period_dist;
}

transformed data {
  int t_max            = t_rep + t_symptom_pre;
  int t_rep_symptoms   = t_symptom_pre + 1;
  real t_rep_symptoms_a[t_rep_symptoms];
  vector[t_rep] t_rep_ones_v;
  matrix[t_rep_symptoms,t_max] line_list;
  int t_r_max = ceil_integer( ( 1.0 * t_max ) / hyper_gp_period_r, 1, t_max );
  int t_t_r_map[ t_max ];
  real prior_r_gp_sd_max_adj;
  int t_dist_max = ceil_integer( ( 1.0 * t_max ) / hyper_gp_period_dist, 1, t_max );
  int t_t_dist_map[ t_max ];
  real prior_alpha_gp_sd_max_adj;
  real prior_beta_gp_sd_max_adj;
  int sdxx;
  int ddxx;
  int tau;
  int<lower=1,upper=t_rep> rdx_conv_min_idx[t_max];
  int<lower=1,upper=t_rep> rdx_conv_max_idx[t_max];
  int<lower=1,upper=t_rep_symptoms> ddx_conv_min_idx[t_max];
  int<lower=1,upper=t_rep_symptoms> ddx_conv_max_idx[t_max];

  // NOTE ON INDEXING
  // reported cases are indexed by:                   rdx = 1..t_rep
  // symptoms are indexed by:                         sdx = 1..t_max
  // the symptom and report time are the same when:   sdx = rdx + t_symptom_pre
  // symptom-report time indexed by:                  ddx = 1..t_rep_sypmtoms
  // actual symptom-report time is:                   tau = ddx - 1 

  // useful vectors for convolution
  for( ddx in 1:t_rep_symptoms ) {
    t_rep_symptoms_a[ ddx ] = ddx;
  }
  t_rep_ones_v = rep_vector( 1, t_rep );

  // convert line-list to relative time from symptom time
  line_list = rep_matrix( 0, t_rep_symptoms, t_max );
  for( idx in 1:n_ll ) {
    tau  = ll_report[ idx ] - ll_symptoms[ idx];
    ddxx = tau + 1;
    ddxx = max( min( ddxx, t_rep_symptoms ), 1 );
    sdxx = max( min( ll_symptoms[idx] + t_symptom_pre, t_max ), 1 );
    line_list[ ddxx, sdxx ] = ll_N[ idx ];
  }

  // gp at a lower frequency, adjust prior and create time conversion maps
  prior_r_gp_sd_max_adj      = prior_r_gp_sd_max / sqrt( hyper_gp_period_r );
  prior_alpha_gp_sd_max_adj  = prior_alpha_gp_sd_max / sqrt( hyper_gp_period_dist );
  prior_beta_gp_sd_max_adj  = prior_beta_gp_sd_max / sqrt( hyper_gp_period_dist );
  for( sdx in 1:t_max ) {
    t_t_r_map[ sdx ]    = ceil_integer( ( 1.0 * sdx ) / hyper_gp_period_r, 1, t_max );
    t_t_dist_map[ sdx ] = ceil_integer( ( 1.0 * sdx ) / hyper_gp_period_dist, 1, t_max );
  }
  
  // When calculating the expected number of reports on each day, we take the
  // convolution of the report-symptom time distribution and the number of 
  // symptomatic. Given the the report-symptom time distribution changes on 
  // each symptomatic date, it is easiest for each day we calculate the number 
  // of symptomatic to distribute them between the report groups.
  for( sdx in 1:t_max ){
    rdx_conv_min_idx[sdx] = max( 1, sdx - t_rep_symptoms + 1 );
    rdx_conv_max_idx[sdx] = min( sdx, t_rep );
    ddx_conv_min_idx[sdx] = max( 1, t_rep_symptoms + 1 - sdx );
    ddx_conv_max_idx[sdx] = min( t_rep_symptoms,t_rep_symptoms + t_rep - sdx );
  }
}

parameters {
  real<lower=prior_log_symptoms0_min,upper=prior_log_symptoms0_max> log_symptoms0;
  real<lower=prior_r_0_min,upper=prior_r_0_max> r_0;
  real<lower=0,upper=prior_r_gp_sd_max_adj> r_gp_sd;
  real<lower=-3*prior_r_gp_sd_max_adj,upper=3*prior_r_gp_sd_max_adj> r_gp[t_r_max];
  real<lower=-3*prior_alpha_gp_sd_max_adj,upper=3*prior_alpha_gp_sd_max_adj> alpha_gp[t_dist_max];
  real<lower=-3*prior_beta_gp_sd_max_adj,upper=3*prior_beta_gp_sd_max_adj> beta_gp[t_dist_max];
  real<lower=0,upper=prior_alpha_gp_sd_max_adj> alpha_gp_sd;
  real<lower=0,upper=prior_beta_gp_sd_max_adj> beta_gp_sd;
  real<lower=0,upper=prior_phi_od_max> phi_od;
  real<lower=prior_alpha_min,upper=prior_alpha_max> alpha0;
  real<lower=prior_beta_min,upper=prior_beta_max> beta0;
}

transformed parameters {
  vector[t_max] symptoms;
  vector[t_max] r;
  real<lower=prior_alpha_min,upper=prior_alpha_max> alpha[t_max];
  real<lower=prior_beta_min,upper=prior_beta_max> beta[t_max];
  vector[t_rep] reported_intensity;
  matrix[t_rep_symptoms,t_max] rep_symptoms_lpdf;
  real dist_cdf[t_rep_symptoms];

  // the number of symptoms cases follows a log normal process
  for( idx in 1:t_max )
    r[idx] = r_gp[ t_t_r_map[ idx ]];
  r[1]  = r_0;
  r     = cumulative_sum(r);
  symptoms = exp( log_symptoms0 + cumulative_sum( r ) );

  // the distribution parameters follow a normal process
  alpha[1]  = alpha0;
  beta[1]  = beta0;
  for( sdx in 2:t_max ) {
    alpha[sdx] = fmax( alpha[sdx-1] + alpha_gp[ t_t_dist_map[ sdx ]], prior_alpha_min );
    beta[sdx]  = fmax( beta[sdx-1] + beta_gp[ t_t_dist_map[ sdx ]], prior_beta_min );
  }

  reported_intensity = rep_vector( 0, t_rep );
  for( sdx in 1:t_max ){
    // calculate the report-symptom distribution for each symptom time
    dist_cdf[1] = gamma_cdf( 1, alpha[sdx], beta[sdx] );
    rep_symptoms_lpdf[ 1, sdx ] = dist_cdf[1];
    for( idx in 2:t_rep_symptoms ) {
       dist_cdf[idx] = gamma_cdf( idx, alpha[sdx], beta[sdx] );
       rep_symptoms_lpdf[ idx, sdx ] = dist_cdf[idx] - dist_cdf[idx-1];
    }

    // take the convoluation with the symptoms to get expected report
    reported_intensity[ rdx_conv_min_idx[sdx]:rdx_conv_max_idx[sdx]] += symptoms[ sdx ] * 
      rep_symptoms_lpdf[ ddx_conv_min_idx[sdx]:ddx_conv_max_idx[sdx], sdx ];
      
    // when fitting the line list data need to take in to account many observations are truncated
    rep_symptoms_lpdf[ , sdx ] = log( rep_symptoms_lpdf[ , sdx ] ) - log( sum( rep_symptoms_lpdf[ ddx_conv_min_idx[sdx]:ddx_conv_max_idx[sdx], sdx ] ) );
  }
}

model {
  // the number of cases with symptoms follows a log normal process
  r_gp      ~ normal( 0, r_gp_sd );
  alpha_gp  ~ normal( 0, alpha_gp_sd );
  beta_gp   ~ normal( 0, beta_gp_sd );

  // the expected number of
  reported ~ neg_binomial_2( reported_intensity, t_rep_ones_v / phi_od );

  // probability of line list
  target += sum( rep_symptoms_lpdf .* line_list );
}

