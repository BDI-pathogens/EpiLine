functions {
  vector johnson_su_lpd( vector x, vector gamma, vector delta, vector lambda, vector xi, int t_rep_symptoms )
  {
    vector[t_rep_symptoms] z1 = ( x - xi ) ./ lambda;
    vector[t_rep_symptoms] z2 = gamma + delta .* asinh( z1 );
    return log( delta ) - log( lambda ) - 0.5 * ( 1.837877 + log( 1 + z1 .* z1 ) + z2 .* z2 );
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
  int<lower=1> t_max;                                      // total days to calculate new symptoms
  int<lower=1,upper=t_max> t_rep_min;                      // first day of reporting
  int<lower=t_rep_min,upper=t_max> t_rep_max;              // final day of reporting
  int<lower=1,upper=t_rep_min-1> t_rep_symptoms_max;       // max time of reporting after symptoms
  int<lower=t_rep_max-t_max+1,upper=0> t_rep_symptoms_min; // min time of reporting before symptoms
  int reported[t_rep_max-t_rep_min+1];                     // number of report cases by day
  int n_ll;
  int ll_report[n_ll];
  int ll_symptoms[n_ll];
  int<lower=0> ll_N[n_ll];
  real prior_gamma_min;
  real prior_delta_min;
  real prior_xi_min;
  real prior_lambda_min;
  real prior_gamma_max;
  real prior_delta_max;
  real prior_xi_max;
  real prior_lambda_max;
  real prior_xi_gp_sd_max;
  real prior_lambda_gp_sd_max;
  real prior_gamma_gp_sd_max;
  real prior_delta_gp_sd_max;
  real prior_log_symptoms0_min;
  real prior_log_symptoms0_max;
  real prior_r_0_min;
  real prior_r_0_max;
  real prior_r_gp_sd_max;
  real prior_phi_od_max;
  int<lower=1,upper=t_max> hyper_gp_period_r;
  int<lower=1,upper=t_max> hyper_gp_period_dist;
}

transformed data {
  int t_rep            = t_rep_max - t_rep_min + 1;
  int t_rep_symptoms   = t_rep_symptoms_max - t_rep_symptoms_min + 1;
  int t_rep_offset_min = t_rep_min - t_rep_symptoms_max;
  int t_rep_offset_max = t_rep_min - t_rep_symptoms_min;
  vector[t_rep_symptoms] t_rep_symptoms_v;
  vector[t_rep] t_rep_ones_v;
  matrix[t_rep_symptoms,t_rep] line_list;
  int rep_symptoms_offset;
  int t_r_max = ceil_integer( ( 1.0 * t_max ) / hyper_gp_period_r, 1, t_max );
  int t_t_r_map[ t_max ];
  real prior_r_gp_sd_max_adj;
  int t_dist_max = ceil_integer( ( 1.0 * t_rep ) / hyper_gp_period_dist, 1, t_rep );
  int t_t_dist_map[ t_rep ];
  real prior_xi_gp_sd_max_adj;
  real prior_lambda_gp_sd_max_adj;
  real prior_gamma_gp_sd_max_adj;
  real prior_delta_gp_sd_max_adj;

  // useful vectors for convolution
  t_rep_symptoms_v = t_rep_symptoms_max - cumulative_sum( rep_vector( 1, t_rep_symptoms ) ) + 1;
  t_rep_ones_v     = rep_vector( 1, t_rep );

  // convert line-list to relative time from report time
  line_list = rep_matrix( 0, t_rep_symptoms, t_rep );
  for( idx in 1:n_ll ) {
    rep_symptoms_offset = t_rep_symptoms_max + 1 - ( ll_report[ idx ] - ll_symptoms[ idx ] );
    rep_symptoms_offset = max( min( rep_symptoms_offset, t_rep_symptoms ), 1 );
    line_list[ rep_symptoms_offset, ll_report[ idx ] - t_rep_min + 1 ] = ll_N[ idx ];
  }

  // gp at a lower frequency, adjust prior and create time conversion maps
  prior_r_gp_sd_max_adj      = prior_r_gp_sd_max / sqrt( hyper_gp_period_r );
  prior_xi_gp_sd_max_adj     = prior_xi_gp_sd_max / sqrt( hyper_gp_period_dist );
  prior_lambda_gp_sd_max_adj = prior_lambda_gp_sd_max / sqrt( hyper_gp_period_dist );
  prior_gamma_gp_sd_max_adj  = prior_gamma_gp_sd_max / sqrt( hyper_gp_period_dist );
  prior_delta_gp_sd_max_adj  = prior_delta_gp_sd_max / sqrt( hyper_gp_period_dist );
  for( t in 1:t_max ) {
    t_t_r_map[ t ]    = ceil_integer( ( 1.0 * t ) / hyper_gp_period_r, 1, t_max );
  }
  for( t in 1:t_rep ) {
    t_t_dist_map[ t ] = ceil_integer( ( 1.0 * t ) / hyper_gp_period_dist, 1, t_max );
  }
}

parameters {
  real<lower=prior_log_symptoms0_min,upper=prior_log_symptoms0_max> log_symptoms0;
  real<lower=prior_r_0_min,upper=prior_r_0_max> r_0;
  real<lower=0,upper=prior_r_gp_sd_max_adj> r_gp_sd;
  real<lower=-3*prior_r_gp_sd_max_adj,upper=3*prior_r_gp_sd_max_adj> r_gp[t_r_max];
  real<lower=-3*prior_xi_gp_sd_max_adj,upper=3*prior_xi_gp_sd_max_adj> xi_gp[t_dist_max];
  real<lower=-3*prior_lambda_gp_sd_max_adj,upper=3*prior_lambda_gp_sd_max_adj> lambda_gp[t_dist_max];
  real<lower=-3*prior_gamma_gp_sd_max_adj,upper=3*prior_gamma_gp_sd_max_adj> gamma_gp[t_dist_max];
  real<lower=-3*prior_delta_gp_sd_max_adj,upper=3*prior_delta_gp_sd_max_adj> delta_gp[t_dist_max];
  real<lower=0,upper=prior_xi_gp_sd_max_adj> xi_gp_sd;
  real<lower=0,upper=prior_lambda_gp_sd_max_adj> lambda_gp_sd;
  real<lower=0,upper=prior_gamma_gp_sd_max_adj> gamma_gp_sd;
  real<lower=0,upper=prior_delta_gp_sd_max_adj> delta_gp_sd;
  real<lower=0,upper=prior_phi_od_max> phi_od;
  real<lower=prior_gamma_min,upper=prior_gamma_max> gamma0;
  real<lower=prior_delta_min,upper=prior_delta_max> delta0;
  real<lower=prior_lambda_min,upper=prior_lambda_max> lambda0;
  real<lower=prior_xi_min,upper=prior_xi_max> xi0;
}

transformed parameters {
  vector[t_max] symptoms;
  vector[t_max] r;
  real<lower=prior_xi_min,upper=prior_xi_max> xi[t_rep];
  real<lower=prior_lambda_min,upper=prior_lambda_max> lambda[t_rep];
  real<lower=prior_gamma_min,upper=prior_gamma_max> gamma[t_rep];
  real<lower=prior_delta_min,upper=prior_delta_max> delta[t_rep];
  vector[t_rep] reported_intensity;
  matrix[t_rep_symptoms,t_rep] rep_symptoms_lpdf;

  // the number of symptoms cases follows a log normal process
  for( idx in 1:t_max )
    r[idx] = r_gp[ t_t_r_map[ idx ]];
  r[1]  = r_0;
  r     = cumulative_sum(r);
  symptoms = exp( log_symptoms0 + cumulative_sum( r ) );

  // the distribution parameters follow a normal process
  xi[1]     = xi0;
  lambda[1] = lambda0;
  gamma[1]  = gamma0;
  delta[1]  = delta0;
  for( idx in 2:t_rep ) {
    xi[idx]     = xi[idx-1] + xi_gp[ t_t_dist_map[ idx ]];
    lambda[idx] = lambda[idx-1] + lambda_gp[ t_t_dist_map[ idx ]];
    gamma[idx]  = gamma[idx-1] + gamma_gp[ t_t_dist_map[ idx ]];
    delta[idx]  = delta[idx-1] + delta_gp[ t_t_dist_map[ idx ]];
  }

  // take the convoluation with the daily symptoms-report distribution
  for( t in 1:t_rep )
    rep_symptoms_lpdf[ , t ] = johnson_su_lpd( t_rep_symptoms_v, rep_vector( gamma[t], t_rep_symptoms ), rep_vector( delta[t], t_rep_symptoms ),
        rep_vector( lambda[t], t_rep_symptoms ),rep_vector( xi[t], t_rep_symptoms ), t_rep_symptoms );

 for( t in 1:t_rep )
    reported_intensity[ t ] = dot_product( symptoms[ (t + t_rep_offset_min ):(t + t_rep_offset_max )], exp( rep_symptoms_lpdf[ , t ]));
}

model {
  // the number of cases with symptoms follows a log normal process
  r_gp      ~ normal( 0, r_gp_sd );
  xi_gp     ~ normal( 0, xi_gp_sd );
  lambda_gp ~ normal( 0, lambda_gp_sd );
  gamma_gp  ~ normal( 0, gamma_gp_sd );
  delta_gp  ~ normal( 0, delta_gp_sd );

  // the expected number of
  reported ~ neg_binomial_2( reported_intensity, t_rep_ones_v / phi_od );

  // probability of line list
  target += sum( rep_symptoms_lpdf .* line_list );
}

