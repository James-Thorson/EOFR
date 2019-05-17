#' Calculate parameter inputs for TMB
#'
#' \code{Param_Fn} generates the \code{parameters} input for \code{TMB::MakeADFun}
#'
#' @param DataList list outputted from \code{EOFR::make_data}
#' @inheritParams make_data


#' @export
make_parameters <-
function( Version, DataList, Rank="Expanded" ){

  # Local function to make a random array
  rarray = function( dim, mean=0, sd=0.01 ) array( rnorm(prod(dim),mean=mean,sd=sd), dim=dim)
  seq_pos <- function( length.out, from=1 ){
    seq(from=from, to=length.out, length.out=max(length.out,0))
  }

  #######################
  # Make Parameters for each version
  #######################

  if(Version%in%c("EOFR_v1_0_0")){
    Return = list("lambda_tf"=rarray(dim=c(DataList$n_t,DataList$n_f),sd=0.2), "ln_H_input"=c(0,0), "logkappa"=log(0.9),
      "alpha_ct"=rarray(dim=c(DataList$n_c,DataList$n_t)), "epsiloninput_scf"=rarray(dim=c(DataList$n_s,DataList$n_c,DataList$n_f)),
      "ln_sigma_c"=rep(1,DataList$n_c), "beta0_p"=rep(0,DataList$n_p), "beta_k"=rep(0,ncol(DataList$X_jk)), "gamma_p"=rep(0,DataList$n_p),
      "ln_sigma_p"=rep(0,DataList$n_p) )
  }

  #######################
  # Fill in values that are shared across versions
  #######################

  # Restrictions on loadings
  if( tolower(Rank) == "expanded" ){
    if( DataList$Cross_correlation==TRUE ){
      Fix = cbind( rep(1,DataList$n_t) %o% rep(FALSE,DataList$n_p), upper.tri(Return[["lambda_tf"]][,-seq_pos(DataList$n_p),drop=FALSE]) )
      Return[["lambda_tf"]][ ifelse(Fix==1,TRUE,FALSE) ] = 0
    }else{
      Return[["lambda_tf"]][upper.tri(Return[["lambda_tf"]])] = 0
    }
  }
  if( tolower(Rank) == "reduced" ){
    Return[["lambda_tf"]][upper.tri(Return[["lambda_tf"]])] = 0
  }
  if( tolower(Rank) == "full" ){
    # Nothing needed
  }

  # Error messages
  if( any(sapply(Return, FUN=function(num){any(is.na(num))})) ) stop("Some parameter is NA")

  # Return tagged list
  return( Return )
}
