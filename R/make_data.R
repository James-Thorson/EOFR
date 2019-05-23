
#' Build data input for VAST model
#'
#' \code{make_data} builds a tagged list of data inputs used by TMB for running the model
#'
#' @param b_i Sampled biomass for each observation i
#' @param a_i Sampled area for each observation i
#' @param c_i Category (e.g., species, length-bin) for each observation i (running from 0 to number of categories minus 1)
#' @param t_i vector specifying time for each observation, e.g., year
#' @param l_i vector specifying which calibration factor is applied to a given observation  (running from 0 to number of categories minus 1); use NA if not calibrated.
#' @param Version a version number (see example for current default).
#' @param spatial_list tagged list of locatoinal information from , i.e., from \code{FishStatsUtils::make_spatial_info}

#' @return Tagged list containing inputs to function \code{VAST::Build_TMB_Fn()}

#' @export
make_data <-
function( Version, B_i, Y_j, c_i, t_i, n_f, t_j, p_j, spatial_list, l_i=c_i, Cross_correlation=TRUE,
  Constrain_orthogonality=FALSE, CheckForErrors=TRUE, X_jk=matrix(0,nrow=length(Y_j),ncol=0) ){

  # Rescale tprime_iz to start at 0
  tmin = min( c(t_i,t_j), na.rm=TRUE)
  tmax = max( c(t_i,t_j), na.rm=TRUE)
  tprime_i = t_i - tmin
  tprime_j = t_j - tmin

  # Determine dimensions
  n_i = length(B_i)
  n_j = length(Y_j)
  n_l = ifelse( all(is.na(l_i)), 0, max(l_i,na.rm=TRUE) + 1 )
  n_t = tmax - tmin + 1
  n_c = max(c_i, na.rm=TRUE) + 1
  n_p = max(p_j, na.rm=TRUE) + 1
  n_g = spatial_list$n_g

  ###################
  # Check for bad data entry
  ###################

  if( CheckForErrors==TRUE ){
    if( n_c!=length(unique(na.omit(as.vector(c_i)))) ) stop("n_c doesn't equal the number of levels in c_i")
  }

  # Check for wrong dimensions
  if( CheckForErrors==TRUE ){
    if( any(c(length(B_i),length(c_i),length(tprime_i))!=n_i) ) stop("b_i, c_i, or tprime_i doesn't have length n_i")
  }

  ###################
  # Check for incompatibilities amongst versions
  ###################

  ###################
  # Check for incompatible settings
  ###################


  ###################
  # switch defaults if necessary
  ###################

  ###################
  # Output tagged list
  ###################

  # CMP_xmax should be >100 and CMP_breakpoint should be 1 for Tweedie model
  Return = NULL
  if(Version%in%c("EOFR_v1_0_0")){
    Return = list( "Cross_correlation"=Cross_correlation, "Constrain_orthogonality"=Constrain_orthogonality, "n_i"=n_i, "n_j"=n_j, "n_s"=spatial_list$MeshList$anisotropic_spde$n.spde, "n_g"=n_g, "n_t"=n_t, "n_c"=n_c, "n_p"=n_p, "n_f"=n_f, "B_i"=B_i, "c_i"=c_i, "t_i"=tprime_i, "Y_j"=Y_j, "X_jk"=X_jk, "p_j"=p_j, "t_j"=tprime_j, "spde_aniso"=list(), "Ais_ij"=cbind(spatial_list$A_is@i,spatial_list$A_is@j), "Ais_x"=spatial_list$A_is@x, "Ags_ij"=cbind(spatial_list$A_gs@i,spatial_list$A_gs@j), "Ags_x"=spatial_list$A_gs@x )
  }
  if(Version%in%c("EOFR_v1_1_0")){
    Return = list( "Cross_correlation"=Cross_correlation, "Constrain_orthogonality"=Constrain_orthogonality, "n_i"=n_i, "n_j"=n_j, "n_l"=n_l, "n_s"=spatial_list$MeshList$anisotropic_spde$n.spde, "n_g"=n_g, "n_t"=n_t, "n_c"=n_c, "n_p"=n_p, "n_f"=n_f, "B_i"=B_i, "c_i"=c_i, "t_i"=tprime_i, "l_i"=l_i, "Y_j"=Y_j, "X_jk"=X_jk, "p_j"=p_j, "t_j"=tprime_j, "spde_aniso"=list(), "Ais_ij"=cbind(spatial_list$A_is@i,spatial_list$A_is@j), "Ais_x"=spatial_list$A_is@x, "Ags_ij"=cbind(spatial_list$A_gs@i,spatial_list$A_gs@j), "Ags_x"=spatial_list$A_gs@x )
  }
  if( is.null(Return) ) stop("`Version` provided does not match the list of possible values")
  if( "spde_aniso" %in% names(Return) ) Return[['spde_aniso']] = list("n_s"=spatial_list$MeshList$anisotropic_spde$n.spde, "n_tri"=nrow(spatial_list$MeshList$anisotropic_mesh$graph$tv), "Tri_Area"=spatial_list$MeshList$Tri_Area, "E0"=spatial_list$MeshList$E0, "E1"=spatial_list$MeshList$E1, "E2"=spatial_list$MeshList$E2, "TV"=spatial_list$MeshList$TV-1, "G0"=spatial_list$MeshList$anisotropic_spde$param.inla$M0, "G0_inv"=INLA::inla.as.dgTMatrix(solve(spatial_list$MeshList$anisotropic_spde$param.inla$M0)) )

  # Return
  return( Return )
}
