# EOF_regression
Package to conduct EOF regression

### Example script

The EOFR package can be run e.g., using the following script.  While long, the majority of the code is either settings, compiling data, or plotting, and each block of code is labeled:

```R
devtools::install_github("james-thorson/FishStatsUtils", ref="2.7.0" )
devtools::install_github("james-thorson/EOFR", ref="development" )

# Set working directory with write access on machine
setwd("C:/Users/James.Thorson/Desktop/Work files/AFSC/2020-08 -- Frederic Maps query about EOFR")

library(EOFR)

###########################
# Get settings
###########################

# Global settings
n_x = 50   # Specify number of knots for predictive process
n_f = 3
Species = "cod"

# Directory
Date = Sys.Date()
  RunFile = DateFile = paste0(getwd(),'/',Date,'/')
  dir.create(RunFile,recursive=TRUE)

Version = get_latest_version( package="EOFR" )
fine_scale = TRUE
Region = "Other"
strata.limits = data.frame('STRATA'="All_areas")
Aniso = FALSE

########################
# Load data
########################

data( EOFR_example )

# Unpack contents
Y_j = EOFR_example$stock_recruit_data[,'ln_recruits_per_spawning_biomass']
X_jk = matrix( EOFR_example$stock_recruit_data[,'spawning_biomass'], ncol=1)
t_j = EOFR_example$stock_recruit_data[,'year']
Data_Geostat = EOFR_example$physical_data

##########################
# Run model
##########################

# Generate grid
Extrapolation_List = make_extrapolation_info( Region=Region, strata.limits=strata.limits,
  observations_LL=cbind("Lon"=Data_Geostat[,"longitude"],"Lat"=Data_Geostat[,"latitude"]), input_grid=Grid, grid_dim_km=c(15,15),
  maximum_distance_from_sample=40 )

# Make spatial info
Spatial_List = make_spatial_info( grid_size_km=1000, n_x=n_x, Method="Mesh", Lon=Data_Geostat[,'longitude'],
  Lat=Data_Geostat[,'latitude'], Extrapolation_List=Extrapolation_List, DirPath=DateFile, Save_Results=TRUE,
  fine_scale=TRUE )

# Build data
TmbData = make_data("Version"=Version, "n_f"=n_f, "B_i"=Data_Geostat[,'response'], "Y_j"=Y_j,  "X_jk"=X_jk,
  "c_i"=as.numeric(Data_Geostat[,'variable'])-1, "p_j"=rep(0,length(Y_j)),
  "t_j"=match(t_j, sort(unique(c(Data_Geostat[,'year'],t_j)))),
  "t_i"=match(Data_Geostat[,'year'], sort(unique(c(Data_Geostat[,'year'],t_j)))),
  "l_i"=as.numeric(Data_Geostat[,'measurement'])-1, "spatial_list"=Spatial_List )

# Build object
TmbList = make_model("TmbData"=TmbData, "RunDir"=DateFile, "Version"=Version, "Aniso"=Aniso, "spatial_list"=Spatial_List )
Obj = TmbList[["Obj"]]

# Optimize
Opt = TMBhelper::fit_tmb( obj=Obj, getsd=TRUE, newtonsteps=1, savedir=RunFile,
  control=list(eval.max=10000,iter.max=10000,trace=1) )
# H = optimHess( par=Opt$par, fn=Obj$fn, gr=Obj$gr )

# Summarize
Report = Obj$report()
ParHat = Obj$env$parList(Opt$par)

####################
# Plots
####################

# Get region-specific settings for plots
MapDetails_List = make_map_info( "Region"=Region, "spatial_list"=Spatial_List, "Extrapolation_List"=Extrapolation_List )
MapDetails_List$Cex = 2
MapDetails_List$MappingDetails = list( "world2Hires", NULL )

# Anisotropy
Report$Range_raw1 = Report$Range_raw2 = Report$Range_raw
plot_anisotropy(FileName=paste0(RunFile,"Aniso.png"), Report=Report, TmbData=list(Options_vec=c(Aniso=Aniso),Options=NA) )

# Indices
Rot = rotate_factors(Cov_jj=NULL, L_pj=Report$L_tf, Psi=aperm(Report$epsiloninput_gcf,c(1,3,2)), RotationMethod="PCA", testcutoff=1e-4)
for(zI in 1:2){
  ThorsonUtilities::save_fig( file=paste0(RunFile,"Factor_indices-",c("Original","Rotated")[zI]), width=4, height=4 )
    if(zI==1) Mat_tf = Report$L_tf
    if(zI==2) Mat_tf = Rot$L_pj_rot
    par( mar=c(2,2,0,0), mgp=c(1.5,0.5,0), tck=-0.02 )
    matplot( y=Mat_tf, x=sort(unique(Data_Geostat[,'year'])), type="l", col=c("black","blue","red"), lty="solid", lwd=2 )
    legend("left", legend=paste0("Factor ",1:n_f), fill=c("black","blue","red"), bty="n" )
  dev.off()
}

# Plot recruitment correlation
ThorsonUtilities::save_fig( file=paste0(RunFile,"Recruitment_correlation"), width=4, height=8 )
  par( mfrow=c(2,1), mar=c(3,3,0,0), mgp=c(1.5,0.5,0), tck=-0.02 )
  #plot( x=ParHat$L_tf[TmbData$t_j+1,1], y=Report$Yhat_j, ylab="", xlab="" )
  Y = cbind(TmbData$Y_j,Report$Yhat_j)
  matplot( x=t_j, y=Y, type="l", lwd=2, col=c("black","red"), ylab="log( Recruits / Spawners )", xlab="" )
  legend( "topleft", bty="n", legend=paste0("Cor: ",formatC(cor(Y[,1],Y[,2]),format="f",digits=2)) )
  plot( x=Y[,2], y=Y[,1], ylab="Observed recruits", xlab="Predicted recruits", xlim=range(Y), ylim=range(Y) )
  abline(a=0, b=1, lwd=2, lty="dotted")
dev.off()

# maps
for( cI in 1:TmbData$n_c ){
  for(zI in 1:2){
    ThorsonUtilities::save_fig( file=paste0(RunFile,"Maps-",c("Original","Rotated")[zI],"-",cI), width=8, height=8 )
      if(zI==1) Mat_sf = Report$epsiloninput_gcf[,cI,]
      if(zI==2) Mat_sf = Rot$Psi_rot[,,cI]
      par( mfrow=c(2,2), mar=c(2,2,0,0), mgp=c(1.5,0.5,0), tck=-0.02 )
      for( fI in 1:TmbData$n_f ){
        plot_variable( map_list=MapDetails_List, xlab="", ylab="", Y_gt=Mat_sf[,fI,drop=FALSE], Format="", add=TRUE )
        if(fI==1) mtext( side=3, text="Spatial map", line=0.5)
        axis(4)
        if(fI==3) axis(1)
      }
      mtext( side=1, text="Longitude", line=2)
      mtext( side=4, text="Latitude", line=1.5, outer=TRUE)
    dev.off()
  }
}
```
