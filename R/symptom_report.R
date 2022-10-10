.eps = 1e-6

##################################################################/
#  Name: symptom_report.fit.class 
###################################################################/
symptom_report.fit.class <- R6Class( 
  "symptom_report.fit.class",
  private = list(
    .stan_raw = NULL,
    .stan_extract = NULL, 
    .stan_params = NULL,
    .stan_data = NULL,
    .fitted_data = NULL,
    .report_date = NULL,
    
    .staticReturn = function( val, name )
    {
      if( !is.null( val ) )
        stop( sprintf( "cannot set %s", name ) )
      
      privateName <- sprintf( ".%s", name )
      return( private[[ privateName ]])
    }
  ),
  active  = list(
    stan_raw     = function( val = NULL ) private$.staticReturn( val, "stan_raw" ),
    stan_extract = function( val = NULL ) private$.staticReturn( val, "stan_extract" ),
    stan_params  = function( val = NULL ) private$.staticReturn( val, "stan_params" ),
    stan_data    = function( val = NULL ) private$.staticReturn( val, "stan_data" ),
    fitted_data  = function( val = NULL ) private$.staticReturn( val, "fitted_data" ),
    report_date  = function( val = NULL ) private$.staticReturn( val, "report_date" )
  ),
  public  = list(
    ##################################################################/
    #  Name: initialize
    ###################################################################/
    initialize = function( stan_raw, stan_data, fitted_data, stan_params, report_date ) {
      private$.stan_raw     <- stan_raw
      private$.stan_extract <- extract( stan_raw )
      private$.stan_data    <- stan_data
      private$.fitted_data  <- fitted_data
      private$.stan_params  <- stan_params
      private$.report_date  <- report_date
    },
    ##################################################################/
    #  Name: plot.symptoms
    ###################################################################/
    plot.symptoms = function( 
      show = TRUE,
      simulation = NULL
    ) 
    { 
      t_rep       <- self$stan_data$t_rep    
      time_offset <- self$stan_data$t_symptom_pre
      t_max       <- t_rep + time_offset + self$stan_data$t_symptom_post
      t_rep_min   <- self$stan_data$t_symptom_pre + 1
      t_rep_max   <- t_rep + time_offset
      extract     <- self$stan_extract
      report_date <- self$report_date
      reported    <- self$fitted_data$reported

      dates   <- (1:t_max) - time_offset
      t_rep_0 <- 0
      if( !is.null( report_date ) ) {
        dates   <- report_date + dates 
        t_rep   <- report_date + t_rep 
        t_rep_0 <- report_date + t_rep_0
      } 
      
      p1 <- plot_ly(
        x = dates,
        y = colQuantiles( extract$symptoms, probs = 0.025 ),
        type = "scatter",
        mode = "lines",
        line = list( color = rgb(0,0,0.5), width= 0 ),
        showlegend = FALSE
      ) %>%
        add_trace( y = colQuantiles(extract$symptoms, probs = 0.975), fill = "tonexty", fillcolor = "rgba(0,0,0.5,0.3)", showlegend = TRUE, name = "CI 5%-95%") %>%
        add_trace( y = colMedians(extract$symptoms ), line = list( width = 5 ), name = "posterior", showlegend = TRUE ) %>%
        add_trace( y = c( rep(NA,t_rep_min-1), reported,rep(NA,t_max-t_rep_max) ), mode = "markers", showlegend = TRUE, name = "reported data" ) %>%
        layout( 
          xaxis  = list( title = list( text = "date" ) ),
          yaxis  = list( title = list( text = "number of people" ) ),
          legend = list( x = 0.05 ),
          shapes = list(
            list( x0=t_rep_0, x1 = t_rep_0, y0 = 0, y1 = 1, type = "line", yref = "paper", marker = list(dash = "dot", width = 1)),
            list( x0=t_rep, x1 = t_rep, y0 = 0, y1 = 1, type = "line", yref = "paper", marker = list(dash = "dot", width = 1))))
      
      if( !is.null( simulation ) ) {
        if( is.R6( simulation ) & class( simulation )[1] == "symptom_report.simulation.class" ) {
          if( length( simulation$symptom ) != length( dates) )
            stop( "simulation time period not consistent with that of fit" )
          p1 <- p1 %>% add_trace( y = simulation$symptom, mode = "markers", name = "symptoms (simulation)", showlegend = TRUE)
          
        } else
          stop( "simulation must be of type symptom_report.simulation.class" )
      }
      
      if( show ) show( p1 )
      return( p1 )
    },
    ##################################################################/
    #  Name: plot.r
    ###################################################################/
    plot.r = function( show = TRUE, simulation = NULL ) 
    {
      t_rep       <- self$stan_data$t_rep    
      time_offset <- self$stan_data$t_symptom_pre
      t_max       <- t_rep + time_offset + self$stan_data$t_symptom_post
      t_rep_min   <- self$stan_data$t_symptom_pre + 1
      t_rep_max   <- t_rep + time_offset
      extract     <- self$stan_extract
      report_date <- self$report_date
      reported    <- self$fitted_data$reported
        
      dates   <- (1:t_max) - time_offset
      t_rep_0 <- 0
      if( !is.null( report_date ) ) {
        dates   <- report_date + dates 
        t_rep   <- report_date + t_rep 
        t_rep_0 <- report_date + t_rep_0
      } 
     
      p1 <- plot_ly(
        x = dates,
        y = colQuantiles( extract$r, probs = 0.05 ),
        type = "scatter",
        mode = "lines",
        line = list( color = rgb(0,0,0.5), width= 0 ),
        showlegend = FALSE
      ) %>%
        add_trace( y = colQuantiles( extract$r, probs = 0.95), fill = "tonexty", fillcolor = "rgba(0,0,0.5,0.3)", showlegend = TRUE, name = "CI 5%-95%") %>%
        add_trace( y = colMedians(extract$r ), line = list( width = 5 ), name = "posterior", showlegend = TRUE ) %>%
        layout( 
          legend = list( x = 0.05 ),
          shapes = list(
            list( x0=t_rep_0, x1 = t_rep_0, y0 = 0, y1 = 1, type = "line", yref = "paper", line = list(dash = "dot", width = 1)),
            list( x0=t_rep, x1 = t_rep, y0 = 0, y1 = 1, type = "line", yref = "paper", line = list(dash = "dot", width = 1))),
          xaxis  = list( title = list( text = "date" ) ),
          yaxis  = list( title = list( text= "r(t)" ) )
        )
      
      if( !is.null( simulation ) ) {
        if( is.R6( simulation ) & class( simulation )[1] == "symptom_report.simulation.class" ) {
          if( length( simulation$r ) != length( dates) )
            stop( "simulation time period not consistent with that of fit" )
          p1 <- p1 %>% add_trace( y = simulation$r, mode = "markers", name = "r (simulation)", showlegend = TRUE)
        } else
          stop( "simulation must be of type symptom_report.simulation.class" )
      }
      if( show ) show( p1 )
      return( p1 )
    },
    ##################################################################/
    #  Name: plot.symptom_report.dist
    ###################################################################/
    plot.symptom_report.dist = function( 
      show = TRUE, 
      simulation = NULL,
      date1 = NULL,
      date2 = NULL
    ) 
    {
      t_rep       <- self$stan_data$t_rep    
      time_offset <- self$stan_data$t_symptom_pre
      t_max       <- t_rep + time_offset + self$stan_data$t_symptom_post
      t_rep_min   <- self$stan_data$t_symptom_pre + 1
      t_rep_max   <- t_rep + time_offset
      extract     <- self$stan_extract
      report_date <- self$report_date
      
      if( !is.null( date1 ) ) 
        t_rep_min <- t_rep_min + as.numeric( date1 - report_date)
      if( !is.null( date2 ) ) 
        t_rep_max <- t_rep_min + as.numeric( date2 - date1)
      
      dist_days   <- (- self$stan_data$t_symptom_post):(self$stan_data$t_symptom_pre)
      dist_from_t <- function( t_dist ) {
        dist <- data.table( 
          xi     = extract$xi[,t_dist], 
          lambda = extract$lambda[,t_dist], 
          gamma  = extract$gamma[,t_dist], 
          delta  = extract$delta[,t_dist],
          dummy  = 1 )
        
        dist <- data.table( x = dist_days, dummy = 1)[ dist, on= "dummy", allow.cartesian = TRUE]
        dist[ , pdf := .djsu( x, xi, lambda, gamma, delta )]
        dist <- dist[ , .(pdf025 = quantile(pdf, probs = 0.025), pdf50 = median(pdf), pdf975=quantile(pdf, probs = 0.975)), by = "x"]
        return( list( dist = dist, date =  report_date + t_dist - time_offset - 1) )  
      }
      dist1 = dist_from_t( t_rep_min )
      dist2 = dist_from_t( t_rep_max )
      
      p1 = plot_ly(
        dist1$dist,
        x = ~x,
        y = ~pdf025,
        type = "scatter",
        mode = "lines",
        line = list( color = rgb(0,0,0.5), width= 0 ),
        name = "CI 5%-95%",
        showlegend = FALSE
      ) %>%
        add_trace( y = ~pdf975, fill = "tonexty", fillcolor = "rgba(0,0,0.5,0.3)") %>%
        add_trace( y = ~pdf50, line = list( width = 5 ), name = dist1$date, showlegend = TRUE ) %>%
        add_trace( data = dist2$dist, line = list( color = rgb(0,0.5,0), width= 0 )) %>%
        add_trace( data = dist2$dist, y = ~pdf975, fill = "tonexty", fillcolor = "rgba(0,0.5,0,0.3)", name = "CI 5%-95%") %>%
        add_trace( data = dist2$dist, y = ~pdf50, line = list( color = rgb(0,0.5,0),  width = 5 ), name = dist2$date, showlegend = TRUE ) %>%
      layout(
        legend = list( x = 0.8),
        xaxis  = list( title = list( text = "symptom-report date interval")),
        yaxis  = list( title = list( text = "posterior density"))
      )
      
      if( !is.null( simulation ) ) {
        if( is.R6( simulation ) & class( simulation )[1] == "symptom_report.simulation.class" ) {
          if( t_rep != simulation$t_rep )
            stop( "simulation time period not consistent with that of fit" )
          if( ( self$stan_data$t_symptom_pre != simulation$t_symptom_pre ) | ( self$stan_data$t_symptom_post != simulation$t_symptom_post ) )
            stop( "simulation t_sympton_pre and _post not consistent with that used in fit")
          xi     <- simulation$dist_xi
          lambda <- simulation$dist_lambda
          gamma  <- simulation$dist_gamma
          delta  <- simulation$dist_delta
          sim1   <- .djsu( dist_days, xi[t_rep_min], lambda[t_rep_min], gamma[t_rep_min], delta[t_rep_min] )
          sim2   <- .djsu( dist_days, xi[t_rep_max], lambda[t_rep_max], gamma[t_rep_max], delta[t_rep_max] )
          p1 <- p1 %>% add_trace( y = sim1, mode = "markers", marker = list( color = rgb(0,0,0.5) ), name = sprintf( "%s (sim)", dist1$date ), showlegend = TRUE )
          p1 <- p1 %>% add_trace( y = sim2, mode = "markers", marker = list( color = rgb(0,0.5,0) ), name = sprintf( "%s (sim)", dist2$date ), showlegend = TRUE )
        } else
          stop( "simulation must be of type symptom_report.simulation.class" )
      }

      if( show ) show( p1 )
      return( p1 )
    },
    ##################################################################/
    #  Name: plot.symptom_report.quantiles
    ###################################################################/
    plot.symptom_report.quantiles= function( 
      show = TRUE, 
      simulation = NULL,
      quantiles =  c( 0.05,0.25, 0.5, 0.75, 0.95)
    ) {
      if( length( quantiles ) > 7 ) 
        stop( "maximum allowable number of quantiles is 7")
      
      t_rep       <- self$stan_data$t_rep    
      time_offset <- self$stan_data$t_symptom_pre
      t_max       <- t_rep + time_offset + self$stan_data$t_symptom_post
      t_rep_min   <- self$stan_data$t_symptom_pre + 1
      t_rep_max   <- t_rep + time_offset
      extract     <- self$stan_extract
      report_date <- self$report_date
      n_dates     <- t_rep_max - t_rep_min + 1
      dates       <- report_date + 0:(n_dates - 1)
      
      qs_05 <- matrix( nrow = n_dates, ncol = length( quantiles) )
      qs_50 <- matrix( nrow = n_dates, ncol = length( quantiles) )
      qs_95 <- matrix( nrow = n_dates, ncol = length( quantiles) )
      for( tdx in 1:n_dates ) {
        t      <- t_rep_min + tdx - 1
        xi     <- extract$xi[,t]
        lambda <- extract$lambda[,t] 
        gamma  <- extract$gamma[,t]
        delta  <- extract$delta[,t]
        for( qdx in 1:length( quantiles ) ) {
          qs <- .qjsu( quantiles[ qdx], xi, lambda, gamma, delta ) 
          qs_05[ tdx, qdx ] <- quantile( qs, probs = 0.05 )
          qs_50[ tdx, qdx ] <- quantile( qs, probs = 0.5 )
          qs_95[ tdx, qdx ] <- quantile( qs, probs = 0.95 )
        }
      }
      
      colour_template = c(
        "rgba(0,0,0.5,%.3f)",
        "rgba(0,0.5,0,%.3f)",
        "rgba(0.5,0,0,%.3f)",
        "rgba(0,0.5,0.5,%.3f)",
        "rgba(0.5,0,0.5,%.3f)",
        "rgba(0,0.5,0,%.3f)",
        "rgba(0.5,0.5,0.5,%.3f)"
      )
      
      p1 <- plot_ly( type = "scatter", mode = "lines", showlegend = FALSE, line = list( width = 0 ) )
      for( qdx in length( quantiles ):1) {
        p1 <- p1 %>% add_trace( x = dates, y = qs_05[ , qdx ] )
        p1 <- p1 %>% add_trace( x = dates, y = qs_95[ , qdx ], fill = "tonexty", 
                                fillcolor = sprintf( colour_template[ qdx], 0.3) )
        p1 <- p1 %>% add_trace( x = dates, y = qs_50[ , qdx ], name = sprintf( "%s percentile", quantiles[qdx] * 100 ),
                                line = list( width = 5, color = sprintf( colour_template[ qdx], 1 ) ), showlegend = TRUE )
      }
      p1 <- p1 %>% layout(
        legend = list( x = 0.75),
        xaxis  = list( title = list( text = "symptom date")),
        yaxis  = list( title = list( text = "symptom-report interval"))
      )

      if( !is.null( simulation ) ) {
        if( is.R6( simulation ) & class( simulation )[1] == "symptom_report.simulation.class" ) {
          if( t_rep != simulation$t_rep )
            stop( "simulation time period not consistent with that of fit" )
          if( ( self$stan_data$t_symptom_pre != simulation$t_symptom_pre ) | ( self$stan_data$t_symptom_post != simulation$t_symptom_post ) )
            stop( "simulation t_sympton_pre and _post not consistent with that used in fit")
          
          for( tdx in 1:n_dates ) {
            t      <- t_rep_min + tdx - 1
            xi     <- simulation$dist_xi[t]
            lambda <- simulation$dist_lambda[t]
            gamma  <- simulation$dist_gamma[t]
            delta  <- simulation$dist_delta[t]
            for( qdx in 1:length( quantiles ) ) {
              qs_50[ tdx, qdx] <- .qjsu( quantiles[ qdx], xi, lambda, gamma, delta ) 
            }
          }
          for( qdx in length( quantiles ):1) {
            p1 <- p1 %>% add_trace( x = dates, y = qs_50[ , qdx ], mode = "markers",
                                    marker = list( color = sprintf( colour_template[ qdx], 1 ) ))
          } 
          
        } else
          stop( "simulation must be of type symptom_report.simulation.class" )
      }
         
      if( show ) show( p1 )
        return( p1 )
    }
  )
)

##################################################################/
#  Name: symptom_report.simulation.class 
###################################################################/
symptom_report.simulation.class <- R6Class( 
  "symptom_report.simulation.class",
  private = list(
    .reported = NULL,
    .r        = NULL,
    .symptom  = NULL,
    .linelist = NULL,
    .dist_xi     = NULL,
    .dist_lambda = NULL,
    .dist_gamma  = NULL,
    .dist_delta  = NULL,
    .t_rep          = NULL,
    .t_symptom_pre  = NULL,
    .t_symptom_post = NULL,
   
    .staticReturn = function( val, name )
    {
      if( !is.null( val ) )
        stop( sprintf( "cannot set %s", name ) )
      
      privateName <- sprintf( ".%s", name )
      return( private[[ privateName ]])
    }
  ),
  
  active  = list(
    reported = function( val = NULL ) private$.staticReturn( val, "reported" ),
    r        = function( val = NULL ) private$.staticReturn( val, "r" ),
    symptom  = function( val = NULL ) private$.staticReturn( val, "symptom" ),
    linelist = function( val = NULL ) private$.staticReturn( val, "linelist" ),
    dist_xi     = function( val = NULL ) private$.staticReturn( val, "dist_xi" ),
    dist_lambda = function( val = NULL ) private$.staticReturn( val, "dist_lambda" ),
    dist_gamma  = function( val = NULL ) private$.staticReturn( val, "dist_gamma" ),
    dist_delta  = function( val = NULL ) private$.staticReturn( val, "dist_delta" ),
    t_rep          = function( val = NULL ) private$.staticReturn( val, "t_rep" ),
    t_symptom_pre  = function( val = NULL ) private$.staticReturn( val, "t_symptom_pre" ),
    t_symptom_post = function( val = NULL ) private$.staticReturn( val, "t_symptom_post" )
  ),
  public  = list(
    ##################################################################/
    #  Name: initialize
    ###################################################################/
    initialize = function( 
      reported,
      r,
      symptom,
      linelist,
      dist_xi,
      dist_lambda,
      dist_gamma,
      dist_delta,
      t_rep,
      t_symptom_pre,
      t_symptom_post
    ) 
    {
      private$.reported <- reported 
      private$.r        <- r  
      private$.symptom  <- symptom 
      private$.linelist <- linelist
      private$.dist_xi     <- dist_xi
      private$.dist_lambda <- dist_lambda
      private$.dist_gamma  <- dist_gamma
      private$.dist_delta  <- dist_delta
      private$.t_rep          <- t_rep
      private$.t_symptom_pre  <- t_symptom_pre
      private$.t_symptom_post <- t_symptom_post
    }
  )
)

##################################################################/
#  Name: symptom_report.simulator
#
#  Description: Runs a simulation of the symptom_report model and
#  returns the line list of cases by day
# 
#  Arguments:
#  t_max          - total time for which data is reported
#  t_symptom_pre  - maximum time before report of onset of symptoms
#  t_symptom_post - maximum time after report of onset of symptoms
#  symptom_0      - initial number of symptomatic individuals
#  r              - daily r(t) at start of simulation (single value or vector of length t_rep, or t_rep+t_symptom_pre+t_symptom_post)
#  dist_xi        - symptom-report distribution xi parameter (Johnson SU)
#  dist_lambda    - symptom-report distribution lambda parameter (Johnson SU)
#  dist_gamma     - symptom-report distribution gamma parameter (Johnson SU)
#  dist_delta     - symptom-report distribution delta parameter (Johnson SU)
#  report_var_over_mean - the variance over mean for negative binomial distribution of daily cases given expected cases
###################################################################/
symptom_report.simulator <- function(
  t_rep          = 50,
  t_symptom_pre  = 30,
  t_symptom_post = 5,
  symptom_0      = 2,
  r              = 0.05,
  dist_xi        = 2,
  dist_lambda    = 2.5,
  dist_gamma     = -2,
  dist_delta     = 1.5,
  report_var_over_mean = 2,
  linelist_coverage = 1
) 
{
  # run the simulation of symptomatic cases for an extended period of time
  # so we can estimate the number of reported cases in the required range
  t_rep_min <- t_symptom_pre
  t_rep_max <- t_rep - t_symptom_post
  t_max     <- t_rep + t_symptom_pre + t_symptom_post 
  if( t_rep < 1 )
    stop( "t_rep must be positive")
  if( t_symptom_pre < 10 )
    stop( "t_symptom_pre must be at least 10")
  if( t_symptom_post < 1 )
    stop( "t_symptom_post must be at least 1")
  
  # check r is acceptable
  if( length( r ) == 1 ) {
    r <- rep( r, t_max )
  } else if( length( r ) == t_rep ) {
    r <- c( rep( r[1], t_symptom_pre ), r, rep( r[1], t_symptom_pre[ t_rep ] ) )
  } else if( length( r ) != t_max )
    stop( "r must be of length 1 or t_rep or t_max")

  # daily symptomic cases is given by 
  symptom <- round( symptom_0 * exp( cumsum( r ) ) )
  
  # check the length of the distribution paramters
  if( length( dist_xi ) == 1 )     dist_xi     <- rep( dist_xi, t_max )
  if( length( dist_lambda ) == 1 ) dist_lambda <- rep( dist_lambda, t_max )
  if( length( dist_gamma ) == 1 )  dist_gamma  <- rep( dist_gamma, t_max )
  if( length( dist_delta ) == 1 )  dist_delta  <- rep( dist_delta, t_max )
  if( length( dist_xi ) != t_max | length( dist_lambda ) != t_max | length( dist_gamma ) != t_max | length( dist_delta ) != t_max  )
    stop( "dist_xxxxx paramters must be of length 1 or t_max" )
  
  # get the expected number reported each day
  report <- rep(0, t_rep)
  for( sdx in 1:t_max ) {
    for( tau in (-t_symptom_post):t_symptom_pre ) {
      rdx = sdx + tau - t_symptom_pre;
      if( rdx < 1 | rdx > t_rep )
        next;
      report[ rdx ] <- report[ rdx ] + symptom[sdx] * .djsu( tau, dist_xi[sdx], dist_lambda[sdx], dist_gamma[sdx], dist_delta[sdx])
    } 
  }
  
  # the number of actual report cases is assumed to be a negative binomial based on this
  if( report_var_over_mean == 1 ) {
    report <- rpois( t_rep, report );
  } else if( report_var_over_mean > 1 ) {
    report <- rnbinom( t_rep, size = report / ( report_var_over_mean - 1 ), mu = report )
  } else 
    stop( "variance over mean must be greater or eqaual to one")
 
  # generate a global linelist
  ll_symptom  <- rep( 1:t_max, symptom )
  ll_gamma    <- rep( dist_gamma, symptom )
  ll_delta    <- rep( dist_delta, symptom)
  ll_xi       <- rep( dist_xi, symptom )
  ll_lambda   <- rep( dist_lambda, symptom )
  ll_report   <- ll_symptom + round( .rjsu( length( ll_symptom), ll_xi, ll_lambda, ll_gamma, ll_delta ) )
  linelist    <- data.table( report = ll_report, symptom = ll_symptom )
  
  linelist[ , report  := report - t_symptom_pre ]
  linelist[ , symptom := symptom - t_symptom_pre ]
  linelist = linelist[ report > 0 & report <= t_rep ]
  
  if( linelist_coverage <= 0 | linelist_coverage > 1 )
    stop( "linelist_coverage miust be (0,1)")
  linelist = linelist[ which( rbinom( linelist[ ,.N ], 1, linelist_coverage ) == 1 ) ]
  
  sim = symptom_report.simulation.class$new( 
    report,
    r,
    symptom,
    linelist,
    dist_xi,
    dist_lambda,
    dist_gamma,
    dist_delta,
    t_rep, 
    t_symptom_pre,
    t_symptom_post
  )
  return( sim ) 
}


##################################################################/
#  Name: .fit_data_from_csv 
#
#  Description:  read the fit data from csv
##################################################################/
.fit_data_from_csv <- function( file_reported, file_linelist )
{
  if( !file.exists( file_reported ) )
    stop( "file_reported does note exist")
  dt   <- fread( file_reported )
  cols <- names( dt )
  if( length( cols ) != 2 | cols[1] != "date" | cols[2] != "reported" )
    stop( "file_reported must contain 2 columns with header date,reported")
  dt[ , date := as.Date( date )]
  dt <- dt[ order( date ) ]
  dt[ , reported := as.integer( reported ) ]
  if( dt[ , max( is.na( date )) ] | dt[ , max( is.null( date )) ] )
    stop( "a date is misformed in the file_reported, should be of form 2022-05-01")
  if( dt[ , max( is.na( reported )) ] | dt[ , max( is.null( reported )) ] | dt[ , min( reported ) ] < 0 )
    stop( "a reported values is misformed in file_reported, should be positive integers")
  # fill in any missing dates
  report_date <- dt[ 1, date ]
  max_date    <- tail( dt[ , date ], 1 )
  dt_clean    <- data.table( date = report_date + ( 0:as.integer(max_date-report_date)) )
  dt_clean    <- dt[ dt_clean, on = "date" ]   
  dt_clean[ , reported := ifelse( is.na( reported ), 0 , reported ) ]
  reported <- dt_clean[ , reported]
  
  if( !file.exists( file_linelist ) )
    stop( "file_linelist does note exist")
  dt   <- fread( file_linelist )
  cols <- names( dt )
  if( length( cols ) != 2 | cols[1] != "report" | cols[2] != "symptom" )
    stop( "file_linelist must contain 2 columns with header report,symptom")
  dt[ , report := as.Date( report )]
  dt[ , symptom := as.Date( symptom )]
  if( dt[ , max( is.na( report)) ] | dt[ , max( is.null( symptom )) ] )
    stop( "a date is misformed in the file_linelist, should be of form 2022-05-01 (remove any entries where symptom date is unknown)")
  
  # remove points outside of the reporting range
  outside_points <- dt[ report < report_date | report > max_date, .N]
  if( outside_points > 0 ) {
    warning( sprintf( "%s report dates in the linelist are outside the reported total data, removing", outside_points ) )
    dt <- dt[ report >= report_date & report <= max_date ]
  }
  
  linelist_symptom <- as.integer( dt[ , symptom ] - report_date ) + 1 
  linelist_report  <- as.integer( dt[ , report ] - report_date ) + 1  
  
  return( list( report_date = report_date, reported = reported, linelist_symptom = linelist_symptom, linelist_report = linelist_report ))
}

##################################################################/
#  Name: symptom_report.fit
#
#  Description: Estimates the parameters in the model from the daily
#  total number of reported cases and the line list of symptom-report pairs
# 
#  Arguments:
#  reported         - vector of number of reported cases by day
#  linelist_symptom - vector of symptom dates (1 per person)
#  linelist_report  - vector of symptom dates (1 per person)
#  t_symptom_pre    - maximum time before report of onset of symptoms
#  t_symptom_post   - maximum time after report of onset of symptoms
###################################################################/
symptom_report.fit <- function(
  reported         = NULL,
  linelist_symptom = NULL,
  linelist_report  = NULL,
  t_symptom_pre  = 30,
  t_symptom_post = 5,
  mcmc_n_samples = 1e2,
  mcmc_n_chains  = 1,
  t_static_dist  = NA,
  prior_gamma_min  = -20,
  prior_delta_min  = -20,
  prior_xi_min     = -20,
  prior_lambda_min = 0.1,
  prior_gamma_max  = 20,
  prior_delta_max  = 20,
  prior_xi_max     = 20,
  prior_lambda_max = 20,
  prior_xi_gp_sd_max      = 0.6,
  prior_lambda_gp_sd_max  = 0.4,
  prior_gamma_gp_sd_max   = 0.4,
  prior_delta_gp_sd_max   = 0.2,
  prior_r_0_min           = -0.1,
  prior_r_0_max           = 0.2,
  prior_r_gp_sd_max       = 0.055,
  prior_phi_od_max        = 1,  
  prior_log_symptoms0_min = log( 0.5 ),
  prior_log_symptoms0_max = log( 50 ),
  hyper_gp_period_r    = 2,
  hyper_gp_period_dist = 2,
  report_date = NULL,
  file_reported = NULL,
  file_linelist = NULL
)
{
    # if data is stored in csv then use it
    if( !is.null( file_reported ) ) 
    {
      data <- .fit_data_from_csv( file_reported, file_linelist )
      report_date      <- data$report_date
      reported         <- data$reported 
      linelist_report  <- data$linelist_report
      linelist_symptom <- data$linelist_symptom
    }
    
    if( is.null( reported ) | is.null( linelist_symptom) | is.null( linelist_report) )
      stop( "must specify reported, linelist_symptom and linelist_report OR file_reported and file_linelist")
    
    # calculate the times require for the simultion
    t_rep <- length( reported )
    t_max <- t_rep + t_symptom_post + t_symptom_pre
  
    # if all dynamic the set static time to 0 else add extra period (i.e. t_static_dist is relative to initial symptoms)
    t_static_dist = ifelse( is.na( t_static_dist ), 0, t_static_dist )  
    
    # calculate totals and symptom-report pairs
    if( length( linelist_symptom ) != length( linelist_report ) )
      stop( "linelist of symptom and report dates must be equal" )
    linelist <- data.table( report = linelist_report, symptom = linelist_symptom )
 
    # enforce the max/min times
    linelist[ , symptom := pmin( pmax( symptom, report - t_symptom_pre ), report + t_symptom_post ) ]
    linelist <- linelist[ , .N, by = c("report", "symptom" ) ][ order(report,symptom) ]
    
    # prepare the data for Stan
    data <- list(
      t_rep          = t_rep,
      t_symptom_pre  = t_symptom_pre,
      t_symptom_post = t_symptom_post,
      reported    = reported,
      n_ll        = linelist[, .N ],
      ll_report   = linelist[ , report ],
      ll_symptoms = linelist[ , symptom ],
      ll_N        = linelist[ , N], 
      t_static_dist = t_static_dist,
      prior_gamma_min  = prior_gamma_min ,
      prior_delta_min  = prior_delta_min,
      prior_xi_min     = prior_xi_min,
      prior_lambda_min = prior_lambda_min,
      prior_gamma_max  = prior_gamma_max,
      prior_delta_max  = prior_delta_max,
      prior_xi_max     = prior_xi_max,
      prior_lambda_max = prior_lambda_max,
      prior_xi_gp_sd_max      = prior_xi_gp_sd_max,
      prior_lambda_gp_sd_max  = prior_lambda_gp_sd_max,
      prior_gamma_gp_sd_max   = prior_gamma_gp_sd_max,
      prior_delta_gp_sd_max   = prior_delta_gp_sd_max,
      prior_log_symptoms0_min = prior_log_symptoms0_min,
      prior_log_symptoms0_max = prior_log_symptoms0_max,
      prior_r_0_min        = prior_r_0_min,
      prior_r_0_max        = prior_r_0_max,
      prior_r_gp_sd_max    = prior_r_gp_sd_max,
      prior_phi_od_max     = prior_phi_od_max,
      hyper_gp_period_r    = hyper_gp_period_r,
      hyper_gp_period_dist = hyper_gp_period_dist
    )
    
    # initialise some of the data in the chain  
    # 1. moment match the Johnson SU parameters for the entire line list
    ll_mean <- mean( linelist[ , rep( report - symptom, N)])
    ll_var  <- var( linelist[ , rep( report - symptom, N)])
    ll_skew <- skewness( linelist[ , rep( report - symptom, N)])
    
    if( prior_delta_min != prior_delta_max ) {
      delta <- max( min( runif( 1,1.5,2), prior_delta_max ), prior_delta_min )
    } else {
      delta <- prior_delta_min  
      data$prior_delta_min <- delta - .eps
      data$prior_delta_max <- delta + .eps
    }

    if( prior_gamma_min != prior_gamma_max ) {
      t_skew <- function( gamma ) return( .jsu.skewness( 1, 1, gamma, delta ) - ll_skew )
      if( t_skew( prior_gamma_min ) * t_skew( prior_gamma_max ) > 0 ) {
        if( ll_skew > 0 )
          gamma <- prior_gamma_min
        else
          gamma <- prior_gamma_max 
      } else
        gamma  <- uniroot( t_skew, lower = prior_gamma_min, upper = prior_gamma_max )$root
      gamma  <- pmin( pmax( gamma, prior_gamma_min + .eps ), prior_gamma_max - .eps )
    } else {
      gamma <- prior_gamma_min  
      data$prior_gamma_min <- gamma - .eps
      data$prior_gamma_max <- gamma + .eps
    }
 
    if( prior_lambda_min != prior_lambda_max ) {
      t_var  <- function( lambda ) return( .jsu.var( 1, lambda, gamma, delta ) - ll_var )
      lambda <- uniroot( t_var, lower = prior_lambda_min, upper = prior_lambda_max)$root
      lambda <- pmin( pmax( lambda, prior_lambda_min + .eps ), prior_lambda_max - .eps )
    } else {
      lambda <- prior_lambda_min  
      data$prior_lambda_min <- lambda - .eps
      data$prior_lambda_max <- lambda + .eps
    }
    
    if( prior_xi_min != prior_xi_max ) {
      t_mean <- function( xi )return( .jsu.mean( xi, lambda, gamma, delta ) - ll_mean )
      xi     <- uniroot( t_mean, lower = prior_xi_min, upper = prior_xi_max)$root
      xi     <- pmin( pmax( xi, prior_xi_min + .eps ), prior_xi_max - .eps )
    } else {
      xi <- prior_xi_min  
      data$prior_xi_min <- xi - .eps
      data$prior_xi_max <- xi + .eps
    }
  
    # 2.initialize r(0)=0 and for the GPs to be constant
    if( t_static_dist == 1 ) {
      n_dist_gp <- 0
    } else {
      t_dist_max <- ifelse( t_static_dist, t_static_dist, t_max )
      n_dist_gp  <- ceiling( t_dist_max / data$hyper_gp_period_dist )
    }
       
    init_func <- function(x) {
      return( list(
        r_gp = rep( 0, ceiling( t_max / data$hyper_gp_period_r ) ),
        r0 = 0.02,
        gamma0 = gamma,
        delta0 = delta,
        lambda0 = lambda,
        xi0 = xi,
        xi_gp     = rep( 0, n_dist_gp ),
        lambda_gp = rep( 0, n_dist_gp ),
        gamma_gp  = rep( 0, n_dist_gp ),
        delta_gp  = rep( 0, n_dist_gp )
      ) ) } 
    
    # get Stan model and sample
    model <- model.symptom_report.stan()
    init  <- lapply( 1:mcmc_n_chains, init_func )
    raw   <-sampling(
      model,
      chains =mcmc_n_chains,
      cores = min( mcmc_n_chains, 3 ),
      iter = mcmc_n_samples,
      data = data,
      init = init,
      pars = c( "symptoms", "r", "delta", "gamma", "lambda", "xi", "r_gp_sd")
    )
    
    stan_params <- list( n_chains = mcmc_n_chains, samples = mcmc_n_samples )
    fitted_data <- list( reported = reported, linelist_report = linelist_report, linelist_symptom = linelist_symptom )
    fit <- symptom_report.fit.class$new( raw, data, fitted_data, stan_params, report_date ) 
    return( fit )
}