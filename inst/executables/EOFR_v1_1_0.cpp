#include <TMB.hpp>
#include <Eigen/Eigenvalues>

// Function for detecting NAs
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

// Input:  n_g, n_f, n_t, logical_flag, Ags_ij, Ags_x, Mat_sc
template<class Type>                                                                                        //
array<Type> project_knots( int n_g, int n_c, int n_f, array<Type> Mat_scf, matrix<int> A_ij, vector<Type> A_x ){
  array<Type> Mat_gcf(n_g, n_c, n_f);
  Mat_gcf.setZero();
  for( int c=0; c<n_c; c++ ){
  for( int Arow=0; Arow<A_ij.rows(); Arow++ ){
  for( int f=0; f<n_f; f++ ){
    int g = A_ij(Arow,0);
    int s = A_ij(Arow,1);
    Mat_gcf(g,c,f) += A_x(Arow) * Mat_scf(s,c,f);
  }}}
  return Mat_gcf;
}

// Space time
template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace R_inla;
  using namespace Eigen;
  using namespace density;

  // Settings
  DATA_SCALAR(Cross_correlation);
  DATA_INTEGER(Constrain_orthogonality);

  // Dimensions
  DATA_INTEGER(n_i);         // Number of physical observations (stacked across all categories and years)
  DATA_INTEGER(n_j);         // Number of biological observations (stacked across species)
  DATA_INTEGER(n_l);         // Number of calibration fields
  DATA_INTEGER(n_s);         // Number of "strata" (i.e., vectices in SPDE mesh)
  DATA_INTEGER(n_g);         // Number of extrapolation-grid cells
  DATA_INTEGER(n_t);         // Number of time-indices
  DATA_INTEGER(n_c);         // Number of categories for physical response
  DATA_INTEGER(n_p);         // Number of biological response variables
  DATA_INTEGER(n_f);         // Number of spatial factors

  // Physical data
  DATA_VECTOR(B_i);       	// Response for each physical observation
  DATA_IVECTOR(c_i);         // Category for each observation
  DATA_IVECTOR(t_i);          // Time-indices (year, season, etc.) for each observation
  DATA_IVECTOR(l_i);          // Calibration index for each observation (NA is no calibration)

  // Biological data
  DATA_VECTOR(Y_j);          // Response for each biological observation
  DATA_MATRIX(X_jk);         // Predictors for each biological observation
  DATA_IVECTOR(p_j);         // Category for each biological observation
  DATA_IVECTOR(t_j);          // Time-index for each biological observation

  // Aniso objects
  DATA_STRUCT(spde_aniso,spde_aniso_t);

  // Projection matrices from knots s to data i or extrapolation-grid cells x
  DATA_IMATRIX( Ais_ij );
  DATA_VECTOR( Ais_x );
  DATA_IMATRIX( Ags_ij );
  DATA_VECTOR( Ags_x );

  // Shared parameters
  PARAMETER_MATRIX(lambda_tf);

  // Physical parameters
  PARAMETER_VECTOR(ln_H_input); // Anisotropy parameters
  PARAMETER(logkappa);             // log-spatial decorrelation rate
  PARAMETER_MATRIX(alpha_ct);       // Year effect
  PARAMETER_ARRAY(epsiloninput_scf);   // Annual variation
  PARAMETER_VECTOR(ln_sigma_c);     // residual variance for physical variable
  PARAMETER_VECTOR(delta_l);        // Mean intercalibration for each calibration-level l
  PARAMETER_VECTOR(ln_sigma_l);      // SD for spatial variation in intercalibration for each calibration-level l
  PARAMETER_MATRIX(deltainput_sl);   // Zero-centered spatial intercalibration for each calibration-level l

  // Biological parameters
  PARAMETER_VECTOR(beta0_p);  // Intercept
  PARAMETER_VECTOR(beta_k);  // slope
  PARAMETER_VECTOR(gamma_p);  // Slope, response to EOF index
  PARAMETER_VECTOR(ln_sigma_p); // residual variance for biological variable

  ////////////////////////
  // Preparatory bookkeeping
  ////////////////////////

  // Indices -- i=Observation; t=Year; c=Category; p=Dynamic-covariate
  int i,j,k,l,t,c,s,f,f2,f3;
  
  // Objective function
  vector<Type> jnll_comp(4);
  // Slot 0 -- Physical data
  // Slot 1 -- Spatial process
  // Slot 2 -- Biological data
  // Slot 3 -- Intercalibration process
  jnll_comp.setZero();
  Type jnll = 0;                

  // Derived parameters
  Type Range_raw;
  Range_raw = sqrt(8) / exp( logkappa );   // Range = approx. distance @ 10% correlation
  vector<Type> sigma_c( ln_sigma_c.size() );
  sigma_c = exp( ln_sigma_c );
  vector<Type> sigma_p( ln_sigma_p.size() );
  sigma_p = exp( ln_sigma_p );
  vector<Type> sigma_l( ln_sigma_l.size() );
  sigma_l = exp( ln_sigma_l );

  // Anisotropy elements
  matrix<Type> H(2,2);
  H(0,0) = exp(ln_H_input(0));
  H(1,0) = ln_H_input(1);
  H(0,1) = ln_H_input(1);
  H(1,1) = (1+ln_H_input(1)*ln_H_input(1)) / exp(ln_H_input(0));

  // Constrain orthogonality
  matrix<Type> L_tf( n_t, n_f );
  L_tf = lambda_tf;
  if( (Constrain_orthogonality==true) & (n_f>=2) ){
    vector<Type> tmpsum( n_f-1 );
    matrix<Type> tmpmat( n_f-1, n_f-1 );
    vector<Type> tmpvec( n_f-1 );
    // Loop through columns
    for( f=1; f<n_f; f++ ){
      tmpsum.setZero();
      for( t=f; t<n_t; t++ ){
      for( f2=0; f2<f; f2++ ){
        tmpsum(f2) += L_tf(t,f2) * L_tf(t,f);
      }}
      tmpmat.setIdentity();
      for( f2=0; f2<f; f2++ ){
      for( f3=0; f3<f; f3++ ){
        tmpmat(f2,f3) = L_tf(f2,f3);
      }}
      tmpvec = -1 * tmpsum.matrix().transpose() * tmpmat.inverse();
      for( f2=0; f2<f; f2++ ){
        L_tf(f2,f) = tmpvec(f2);
      }
    }
  }

  // Biological covariates
  vector<Type> eta_j( n_j );
  eta_j.setZero();
  for( j=0; j<n_j; j++ ){
  for( k=0; k<X_jk.cols(); k++ ){
    eta_j(j) += X_jk(j,k) * beta_k(k);
  }}

  ////////////////////////
  // Random effects
  //  Spatial and spatio-temporal variation
  // (Should be first to optimize speedup involving normalization of GMRFs, but must be after interactions)
  ////////////////////////

  // Random field probability
  Eigen::SparseMatrix<Type> Q( n_s, n_s );
  Q = Q_spde(spde_aniso, exp(logkappa), H);
  GMRF_t<Type> gmrf_Q;
  gmrf_Q = GMRF( Q );
  Type logtau = log( 1 / (exp(logkappa) * sqrt(4*M_PI)) );

  // Likelihood of epsiloninput
  for(f=0; f<n_f; f++){
  for(c=0; c<n_c; c++){
    jnll_comp(1) += gmrf_Q( epsiloninput_scf.col(f).col(c) );
  }}

  // Likelihood of deltainput_sl
  for(l=0; l<n_l; l++){
    jnll_comp(3) += gmrf_Q( deltainput_sl.col(l) );
  }

  // Transformation of epsiloninput_scf
  array<Type> epsilon_sct(n_s, n_c, n_t);
  epsilon_sct.setZero();
  for(s=0; s<n_s; s++){
  for(t=0; t<n_t; t++){
  for(f=0; f<n_f; f++){
  for(c=0; c<n_c; c++){
    // PDF for first year of autoregression
    epsilon_sct(s,c,t) += epsiloninput_scf(s,c,f)/exp(logtau) * L_tf(t,f);
  }}}}

  // Projection for epsilon_sct to data
  vector<Type> epsilon_i(n_i);
  epsilon_i.setZero();
  for( int Arow=0; Arow<Ais_ij.rows(); Arow++ ){
    i = Ais_ij(Arow,0);
    s = Ais_ij(Arow,1);
    epsilon_i(i) += Ais_x(Arow) * epsilon_sct(s,c_i(i),t_i(i));
  }

  // Transformation for deltainput_sl
  matrix<Type> delta_sl( deltainput_sl.rows(), deltainput_sl.cols() );
  delta_sl.setZero();
  for(s=0; s<n_s; s++){
  for(l=0; l<n_l; l++){
    delta_sl(s,l) += deltainput_sl(s,l)/exp(logtau) * sigma_l(l);
  }}

  // Projection for delta_sl to data
  vector<Type> delta_i(n_i);
  delta_i.setZero();
  for( int Arow=0; Arow<Ais_ij.rows(); Arow++ ){
    i = Ais_ij(Arow,0);
    s = Ais_ij(Arow,1);
    if( (l_i(i)>=0) & (l_i(i)<n_l) ){
      delta_i(i) += delta_l(l_i(i)) + ( Ais_x(Arow) * delta_sl(s,l_i(i)) );
    }
  }
  REPORT( delta_sl );
  REPORT( delta_i );
  REPORT( sigma_l );

  ////////////////////////
  // Likelihood for physical data
  ////////////////////////

  // Derived quantities
  vector<Type> Bhat_i(n_i);
  vector<Type> ln_prob_i(n_i);
  ln_prob_i.setZero();

  // Likelihood contribution from observations
  for(i=0; i<n_i; i++){
    // Linear predictors
    Bhat_i(i) = alpha_ct(c_i(i),t_i(i)) + epsilon_i(i) + delta_i(i);
    if( !isNA(B_i(i)) ){
      // Likelihood for delta-models with continuous positive support
      ln_prob_i(i) = dnorm(B_i(i), Bhat_i(i), sigma_c(c_i(i)), true);
    }
  }

  ////////////////////////
  // Likelihood for biological data
  ////////////////////////

  // Derived quantities
  vector<Type> Yhat_j(n_j);
  vector<Type> ln_prob_j(n_j);
  ln_prob_j.setZero();

  // Likelihood contribution from observations
  for(j=0; j<n_j; j++){
    // Linear predictors
    Yhat_j(j) = beta0_p(p_j(j)) + eta_j(j) + gamma_p(p_j(j)) * L_tf(t_j(j),0) * Cross_correlation;
    if( !isNA(Y_j(j)) ){
      // Likelihood for delta-models with continuous positive support
      ln_prob_j(j) = dnorm(Y_j(j), Yhat_j(j), sigma_p(p_j(j)), true);
    }
  }


  /////////////////////////
  //  Joint likelihood
  /////////////////////////

  jnll_comp(0) = -1 * (ln_prob_i).sum();
  jnll_comp(2) = -1 * (ln_prob_j).sum();
  jnll = jnll_comp.sum();

  ////////////////////////
  // Diagnostic outputs
  ////////////////////////

  REPORT( Q );
  REPORT( Bhat_i );
  REPORT( ln_prob_i );
  REPORT( sigma_c );
  REPORT( Yhat_j );
  REPORT( ln_prob_j );
  REPORT( sigma_p );
  REPORT( epsilon_sct );
  REPORT( H );
  REPORT( Range_raw );
  REPORT( jnll_comp );
  REPORT( jnll );
  REPORT( L_tf );
  REPORT( lambda_tf );
  REPORT( eta_j );

  ADREPORT( Yhat_j );

  // Calculate value of vactors at extrapolation-grid cells (e.g., for use when visualizing estimated or rotated factor estimates)
  array<Type> epsiloninput_gcf( n_g, n_c, n_f );
  epsiloninput_gcf = project_knots( n_g, n_c, n_f, epsiloninput_scf, Ags_ij, Ags_x );
  REPORT( epsiloninput_gcf );

  // Projection from epsilon to observations
  array<Type> epsilon_gct(n_g, n_c, n_t);
  epsilon_gct = project_knots( n_g, n_c, n_t, epsilon_sct, Ags_ij, Ags_x );
  REPORT( epsilon_gct );

  return jnll;
}
