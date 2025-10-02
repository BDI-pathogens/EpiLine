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
  int<lower=1> t_rep;           // total days of reporting
  int<lower=1> t_symptom_pre;   // max time of reporting after symptoms
  int<lower=1> t_symptom_post;  // max time of reporting before symptoms
  array[t_rep] int reported;          // number of report cases by day
  int n_ll;
  array[n_ll] int ll_report;
  array[n_ll] int ll_symptoms;
  array[n_ll] int<lower=0> ll_N;
  int<lower=0,upper=t_rep+t_symptom_pre+t_symptom_post> t_static_dist; // freeze delay distribution after this time (-1=no freeze)
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
  int<lower=1,upper=t_rep + t_symptom_pre + t_symptom_post> hyper_gp_period_r;
  int<lower=1,upper=t_rep + t_symptom_pre + t_symptom_post> hyper_gp_period_dist;
}

transformed data {
  int t_max            = t_rep + t_symptom_pre + t_symptom_post;
  int t_rep_symptoms   = t_symptom_pre + t_symptom_post + 1;
  vector[t_rep_symptoms] t_rep_symptoms_v;
  vector[t_rep] t_rep_ones_v;
  matrix[t_rep_symptoms,t_max] line_list;
  int t_r_max = ceil_integer( ( 1.0 * t_max ) / hyper_gp_period_r, 1, t_max );
  array[ t_max ] int t_t_r_map;
  real prior_r_gp_sd_max_adj;
  int t_dist_max = ceil_integer( ( 1.0 * t_max ) / hyper_gp_period_dist, 1, t_max );
  int t_static_dist_adj = t_static_dist;
  array[ t_max ] int t_t_dist_map;
  real prior_xi_gp_sd_max_adj;
  real prior_lambda_gp_sd_max_adj;
  real prior_gamma_gp_sd_max_adj;
  real prior_delta_gp_sd_max_adj;
  int sdxx;
  int ddxx;
  int tau;
  array[t_max] int<lower=1,upper=t_rep> rdx_conv_min_idx;
  array[t_max] int<lower=1,upper=t_rep> rdx_conv_max_idx;
  array[t_max] int<lower=1,upper=t_rep_symptoms> ddx_conv_min_idx;
  array[t_max] int<lower=1,upper=t_rep_symptoms> ddx_conv_max_idx;

  // NOTE ON INDEXING
  // reported cases are indexed by:                   rdx = 1..t_rep
  // symptoms are indexed by:                         sdx = 1..t_max
  // the symptom and report time are the same when:   sdx = rdx + t_symptom_pre
  // symptom-report time indexed by:                  ddx = 1..t_rep_sypmtoms
  // actual symptom-report time is:                   tau = ddx - 1 - t_symptom_post

  // useful vectors for convolution
  for( ddx in 1:t_rep_symptoms )
    t_rep_symptoms_v[ ddx ] = -t_symptom_post + ddx - 1;
  t_rep_ones_v = rep_vector( 1, t_rep );

  // convert line-list to relative time from symptom time
  line_list = rep_matrix( 0, t_rep_symptoms, t_max );
  for( idx in 1:n_ll ) {
    tau  = ll_report[ idx ] - ll_symptoms[ idx];
    ddxx = tau + t_symptom_post + 1;
    ddxx = max( min( ddxx, t_rep_symptoms ), 1 );
    sdxx = max( min( ll_symptoms[idx] + t_symptom_pre, t_max ), 1 );
    line_list[ ddxx, sdxx ] = ll_N[ idx ];
  }

  // gp at a lower frequency, adjust prior and create time conversion maps
  prior_r_gp_sd_max_adj      = prior_r_gp_sd_max / sqrt( hyper_gp_period_r );
  prior_xi_gp_sd_max_adj     = prior_xi_gp_sd_max / sqrt( hyper_gp_period_dist );
  prior_lambda_gp_sd_max_adj = prior_lambda_gp_sd_max / sqrt( hyper_gp_period_dist );
  prior_gamma_gp_sd_max_adj  = prior_gamma_gp_sd_max / sqrt( hyper_gp_period_dist );
  prior_delta_gp_sd_max_adj  = prior_delta_gp_sd_max / sqrt( hyper_gp_period_dist );
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
  
  // if have a static period in the distributions
  if( t_static_dist > 0 )
  {
    if( t_static_dist == 1 ) {
      t_dist_max = 0;
    } else {
      t_dist_max = ceil_integer( ( 1.0 * t_static_dist ) / hyper_gp_period_dist, 1, t_static_dist);
    }
  } else {
    t_static_dist_adj = t_max;
  }
}

parameters {
  real<lower=prior_log_symptoms0_min,upper=prior_log_symptoms0_max> log_symptoms0;
  real<lower=prior_r_0_min,upper=prior_r_0_max> r_0;
  real<lower=0,upper=prior_r_gp_sd_max_adj> r_gp_sd;
  array[t_r_max] real<lower=-3*prior_r_gp_sd_max_adj,upper=3*prior_r_gp_sd_max_adj> r_gp;
  array[t_dist_max] real<lower=-3*prior_xi_gp_sd_max_adj,upper=3*prior_xi_gp_sd_max_adj> xi_gp;
  array[t_dist_max] real<lower=-3*prior_lambda_gp_sd_max_adj,upper=3*prior_lambda_gp_sd_max_adj> lambda_gp;
  array[t_dist_max] real<lower=-3*prior_gamma_gp_sd_max_adj,upper=3*prior_gamma_gp_sd_max_adj> gamma_gp;
  array[t_dist_max] real<lower=-3*prior_delta_gp_sd_max_adj,upper=3*prior_delta_gp_sd_max_adj> delta_gp;
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
  array[t_max] real<lower=prior_xi_min,upper=prior_xi_max> xi;
  array[t_max] real<lower=prior_lambda_min,upper=prior_lambda_max> lambda;
  array[t_max] real<lower=prior_gamma_min,upper=prior_gamma_max> gamma;
  array[t_max] real<lower=prior_delta_min,upper=prior_delta_max> delta;
  vector[t_rep] reported_intensity;
  matrix[t_rep_symptoms,t_max] rep_symptoms_lpdf;

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
  for( sdx in 2:t_static_dist_adj ) {
    xi[sdx]     = xi[sdx-1] + xi_gp[ t_t_dist_map[ sdx ]];
    lambda[sdx] = lambda[sdx-1] + lambda_gp[ t_t_dist_map[ sdx ]];
    gamma[sdx]  = gamma[sdx-1] + gamma_gp[ t_t_dist_map[ sdx ]];
    delta[sdx]  = delta[sdx-1] + delta_gp[ t_t_dist_map[ sdx ]];
  }
  for( sdx in (t_static_dist_adj+1):t_max ){
    xi[sdx]     = xi[sdx-1];
    lambda[sdx] = lambda[sdx-1];
    gamma[sdx]  = gamma[sdx-1];
    delta[sdx]  = delta[sdx-1];
  }

  reported_intensity = rep_vector( 0, t_rep );
  for( sdx in 1:t_max ){
    // calculate the report-symptom distribution for each symptom time
    rep_symptoms_lpdf[ , sdx ] = johnson_su_lpd( t_rep_symptoms_v, rep_vector( gamma[sdx], t_rep_symptoms ), rep_vector( delta[sdx], t_rep_symptoms ),
        rep_vector( lambda[sdx], t_rep_symptoms ),rep_vector( xi[sdx], t_rep_symptoms ), t_rep_symptoms );
  
    // take the convoluation with the symptoms to get expected report
    reported_intensity[ rdx_conv_min_idx[sdx]:rdx_conv_max_idx[sdx]] += symptoms[ sdx ] * 
      exp( rep_symptoms_lpdf[ ddx_conv_min_idx[sdx]:ddx_conv_max_idx[sdx], sdx ]);
      
    // when fitting the line list data need to take in to account many observations are truncated
    rep_symptoms_lpdf[ , sdx ] -= log( sum( exp( rep_symptoms_lpdf[ ddx_conv_min_idx[sdx]:ddx_conv_max_idx[sdx], sdx ]) ) );
  }
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

