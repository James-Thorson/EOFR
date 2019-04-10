# EOF_regression
Package to conduct EOF regression

### Example script

The EOFR package can be run e.g., using the following script.  While long, the majority of the code is either settings, compiling data, or plotting, and each block of code is labeled:

```R

devtools::install_github("james-thorson/FishStatsUtils", ref="development" )
devtools::install_github("james-thorson/EOFR", ref="development" )

setwd("D:/UW Hideaway (SyncBackFree)/Collaborations/2019 -- EOF regression")

library(TMB)               # Can instead load library(TMBdebug)
library(EOFR)

###########################
# Get settings
###########################

# Global settings
n_x = 50   # Specify number of stations (a.k.a. "knots")
n_f = 3
#Cross_correlation = FALSE
Rank_expanded = TRUE
Species = "cod"
use_REML = TRUE
intercept_structure = "category"
Cross_correlation = TRUE

# Directory
Date = Sys.Date()
  RunFile = DateFile = paste0(getwd(),'/',Date,'/')
  dir.create(RunFile,recursive=TRUE)

Version = get_latest_version( package="EOFR" )
fine_scale = TRUE
Region = "Other"
strata.limits = data.frame('STRATA'="All_areas")
Aniso = FALSE
Constrain_orthogonality = FALSE

########################
# Load data
########################

# Load stock-recruit data
SpeciesCode = 'PCODEBS2016'
R_tc = read.csv( paste0("Data/SARA download/R_tc.csv") )
S_tc = read.csv( paste0("Data/SARA download/S_tc.csv") )
Y_j = log( R_tc[,SpeciesCode] / S_tc[,SpeciesCode] )
X_jk = matrix( S_tc[,SpeciesCode], ncol=1 ) / 1e6
t_j = R_tc[,'X']

# Load EBS bottom temperature data
BottomTemp = NULL
for( t in 1982:2016 ){
  Temp = read.table( file=paste0("Data/EBST/temp",t,".csv") )
  BottomTemp = rbind( BottomTemp, cbind("Year"=t,Temp))
}
Data_Geostat = data.frame( "spp"=1, "Year"=BottomTemp[,"Year"], "Catch_KG"=BottomTemp[,"V3"], "AreaSwept_km2"=1, "Vessel"=0, "Lat"=BottomTemp[,"V2"], "Lon"=BottomTemp[,"V1"] )
Data_Geostat[,'Lon'] = Data_Geostat[,'Lon'] - 360

# Load surface temperature
DF = read.csv( paste0("Data/SST/SST_DF.csv") )
Data2 = data.frame( "spp"=2, "Year"=DF[,"Year"], "Catch_KG"=DF[,"SST"], "AreaSwept_km2"=1, "Vessel"=0, "Lat"=DF[,"Latitude"], "Lon"=DF[,"Longitude"] )
Data2[,'Lon'] = Data2[,'Lon'] - 360
Data2[,'spp'] = 2

# Reduce down surface data to same footprint as EBS data
NN = RANN::nn2( query=Data2[,c("Lat","Lon")], data=Data_Geostat[,c("Lat","Lon")], k=1 )
Data2 = Data2[ which(NN$nn.dists <= 1), ]
Data2 = Data2[ which(Data2[,'Year'] %in% unique(Data_Geostat[,'Year'])), ]

# Combine
Data_Geostat = rbind( Data_Geostat, Data2 )

##########################
# Run model
##########################

# Generate grid
Extrapolation_List = make_extrapolation_info( Region=Region, strata.limits=strata.limits, observations_LL=Data_Geostat[,c("Lon","Lat")], input_grid=Grid, grid_dim_km=c(15,15), maximum_distance_from_sample=40 )

# Make spatial info
Spatial_List = make_spatial_info( grid_size_km=1000, n_x=n_x, Method="Mesh", Lon=Data_Geostat[,'Lon'], Lat=Data_Geostat[,'Lat'], Extrapolation_List=Extrapolation_List, DirPath=DateFile, Save_Results=TRUE, fine_scale=TRUE )

# Build data
TmbData = make_data("Version"=Version, "n_f"=n_f, "B_i"=Data_Geostat[,'Catch_KG'], "Y_j"=Y_j,  "X_jk"=X_jk,
  "c_i"=as.numeric(Data_Geostat[,'spp'])-1, "p_j"=rep(0,length(Y_j)), "t_i"=as.numeric(Data_Geostat[,'Year']), "t_j"=t_j,
  "spatial_list"=Spatial_List, "Cross_correlation"=Cross_correlation, "Constrain_orthogonality"=Constrain_orthogonality )

# Build object
  # dyn.load( paste0(DateFile,"/",TMB::dynlib(Version)) )
TmbList = make_model("TmbData"=TmbData, "RunDir"=DateFile, "Version"=Version, "TmbDir"=TmbDir, "use_REML"=use_REML,
  "Aniso"=Aniso, "spatial_list"=Spatial_List, "Rank_expanded"=Rank_expanded, "intercept_structure"=intercept_structure )
Obj = TmbList[["Obj"]]

# Optimize
Opt = TMBhelper::Optimize( obj=Obj, getsd=TRUE, newtonsteps=1, savedir=RunFile,
  control=list(eval.max=10000,iter.max=10000,trace=1) )

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
Rot = Rotate_Fn(Cov_jj=NULL, L_pj=Report$L_tf, Psi=aperm(Report$epsiloninput_gcf,c(1,3,2)), RotationMethod="PCA", testcutoff=1e-4)
for(zI in 1:2){
  ThorsonUtilities::save_fig( file=paste0(RunFile,"Factor_indices-",c("Original","Rotated")[zI]), width=4, height=4 )
    if(zI==1) Mat_tf = Report$L_tf
    if(zI==2) Mat_tf = Rot$L_pj_rot
    par( mar=c(2,2,0,0), mgp=c(1.5,0.5,0), tck=-0.02 )
    matplot( y=Mat_tf, x=sort(unique(Data_Geostat[,'Year'])), type="l", col=c("black","blue","red"), lty="solid", lwd=2 )
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
Cex = c(0.6,0.85,1.05)[3]
Col = colorRampPalette(colors=c("darkblue","blue","lightblue","lightgreen","yellow","orange","red"))(1000)
for( cI in 1:TmbData$n_c ){
  for(zI in 1:2){
    ThorsonUtilities::save_fig( file=paste0(RunFile,"Maps-",c("Original","Rotated")[zI],"-",cI), width=8, height=8 )
      if(zI==1) Mat_sf = Report$epsiloninput_gcf[,cI,]
      if(zI==2) Mat_sf = Rot$Psi_rot[,,cI]
      par( mfrow=c(2,2), mar=c(2,2,0,0), mgp=c(1.5,0.5,0), tck=-0.02 )
      for( fI in 1:TmbData$n_f ){
        PlotMap_Fn( MappingDetails=MapDetails_List[["MappingDetails"]], xlab="", ylab="", Mat=Mat_sf[,fI,drop=FALSE], Format="", add=TRUE, PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], FileName=paste0(RunFile,"Factor_maps--","Omega2"), Year_Set="", Rotate=MapDetails_List[["Rotate"]], zone=MapDetails_List[["Zone"]], pch=15, Cex=Cex, mfrow=c(1,1), Legend=MapDetails_List[["Legend"]], plot_legend_fig=TRUE, land_color="grey", Col=Col)
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
