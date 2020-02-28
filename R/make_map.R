#' @export
make_map <-
function( DataList, TmbParams, Rank="Expanded", Aniso=TRUE, intercept_structure="none",
  sigmac_structure="separate" ){

  # Local functions
  fix_value <- function( fixvalTF ){
    vec = rep(0,length(fixvalTF))
    if(sum(fixvalTF)>0) vec[which(fixvalTF==1)] = NA
    if(sum(!fixvalTF)>0) vec[which(!is.na(vec))] = 1:sum(!is.na(vec))
    vec = factor( vec ) 
    return( vec )
  }
  seq_pos <- function( length.out, from=1 ){
    seq(from=from, to=length.out, length.out=max(length.out,0))
  }

  # Create tagged-list in TMB format for fixing parameters
  Map = list()

  # Restrictions on loadings
  Map[["lambda_tf"]] = matrix( 1:(DataList$n_t*DataList$n_f), ncol=DataList$n_f )
  if( tolower(Rank) == "expanded" ){
    if( DataList$Cross_correlation==TRUE ){
      Fix = cbind( rep(1,DataList$n_t) %o% rep(FALSE,DataList$n_p), upper.tri(Map[["lambda_tf"]][,-seq_pos(DataList$n_p),drop=FALSE]) )
      Map[["lambda_tf"]][ ifelse(Fix==1,TRUE,FALSE) ] = NA
    }else{
      Map[["lambda_tf"]][upper.tri(Map[["lambda_tf"]])] = NA
    }
  }
  if( tolower(Rank) == "reduced" ){
    Map[["lambda_tf"]][upper.tri(Map[["lambda_tf"]])] = NA
  }
  if( tolower(Rank) == "full" ){
    # Nothing needed
  }
  Map[["lambda_tf"]] = as.factor(Map[["lambda_tf"]])

  # Restrictions on cross-correlation
  if( DataList$Cross_correlation==FALSE ){
    Map[["gamma_p"]] = as.factor(rep(NA,DataList$n_p))
  }

  # Restriction on intercepts
  if( tolower(intercept_structure) == "category" ){
    Map[["alpha_ct"]] = as.factor( 1:DataList$n_c %o% rep(1,DataList$n_t) )
  }
  if( tolower(intercept_structure) == "none" ){
    Map[["alpha_ct"]] = as.factor( array(NA, dim=c(DataList$n_c,DataList$n_t)) )
  }
  if( tolower(intercept_structure) == "separate" ){
    # Nothing needed to do
  }

  # Restrictions on ln_sigma_c
  if( tolower(sigmac_structure) == "identical" ){
    Map[["ln_sigma_c"]] = as.factor( rep(1,DataList$n_c) )
  }

  # Anisotropy
  if( Aniso==FALSE ){
    Map[["ln_H_input"]] = as.factor(rep(NA,length(TmbParams$ln_H_input)))
  }

  # Return
  return(Map)
}

