#' Create polygons from point locations that describe polygon edges
#'
#' Useful where faults are stored as points in a dataframe (e.g. UCERF 3 fault geometries).
#' See LineMaker for line alternative.
#' @param rownum row number for fault of interest
#' @param cols columns in dataframe to be used in creating polygons
#' @param Geom fault geometry stored in a dataframe
#' @return returns SpatialPolygonsDataFrame object
#' @keywords polygon
#' @export
#' @examples
#' To run with UCERF3 geometry, applying to all faults and skipping the first 8 columns:
#'
#' NFaults <- length(UCERFFaultGeom$Name)
#' #x <- seq(1, NFaults,1)
#' ### Apply faultgeometry function to all faults, make list of fault polygons
#' FaultList <- lapply(7:NFaults, PolygonMaker, Geom=UCERFFaultGeom, cols=9:82)
#' ### Make list into SpatialPolygons (sp faff)
#' FaultPolys <- SpatialPolygons(lapply(FaultList, function(x){x@polygons[[1]]}))
#'

PolygonMaker <- function(rownum, cols, Geom){
  ### Polygons are stored in a row with alternating lat/lon columns and NAs when no further points are necessary
  ### First column is fault name, keep this for polygon ID
  PolyName <- as.character(Geom$Name[rownum])
  ### as.numeric is probably not necessary here, because if it is you might have factors and you should be worried about if this will mess everything up...
  PolyCoords <- as.numeric(Geom[rownum, cols])
  ### Drop NAs (these are just here so num of columns is consistent)
  PolyCoordsOnly <- PolyCoords[is.na(PolyCoords) == FALSE]
  ### Extract alternating lats and lons
  PolyLats <- PolyCoordsOnly[seq(1, length(PolyCoordsOnly), 2)]
  PolyLons <- PolyCoordsOnly[seq(2, length(PolyCoordsOnly), 2)]
  ### Convert lats/lons into Polygon object
  ### NB: matrix needs c - if you don't explicitly hand it a list it will only use one dataset and just ignore the other
  Poly <- Polygon(matrix(c(PolyLons, PolyLats), nrow=length(PolyLats)))
  ### Convert Polygon object to Polygons object, attach name of Polygon/fault
  P2 <- Polygons(list(Poly), ID=PolyName)
  ### Attach CRS to polygon
  P3 <- SpatialPolygons((list(P2)),  proj4string = CRS("+proj=longlat +datum=WGS84"))
  ###  - Should get this into projected CR so we are working in km/m rather than degrees - can use sp_transform but you'll need to figure out suitable projection
  return(P3)
}


#' Create line objects from point locations that describe points of the line
#'
#' Useful where faults are stored as points in a dataframe (e.g. UCERF 3 fault geometries).
#' See PolygonMaker for polygon alternative.
#' @param rownum row number for fault of interest
#' @param Geom fault geometry stored in a dataframe
#' @param cols columns in dataframe to be used in creating polygons
#' @return returns SpatialLinesDataFrame object
#' @keywords line
#' @export
#' @examples
#' To run with UCERF3 geometry, applying to all faults and skipping the first 7 rows:
#' NFaults <- length(UCERFFaultGeom$Name)
#' FaultLineList <- lapply(1:NFaults, LineMaker, Geom=UCERFFaultGeom, cols=9:82)
#' FaultLines <- SpatialLines(lapply(FaultLineList, function(x){x@lines[[1]]}))

LineMaker <- function(rownum, Geom, cols){
  ### Polygons are stored in a row with alternating lat/lon columns and NAs when no further points are necessary
  ### First column is fault name, keep this for polygon ID
  FaultName <- as.character(Geom$Name[rownum])
  ### as.numeric is probably not necessary here, because if it is you might have factors and you should be worried about if this will mess everything up...
  FaultCoords <- as.numeric(Geom[rownum,cols])
  ### Drop NAs (these are just here so num of columns is consistent)
  FaultCoordsOnly <- FaultCoords[is.na(FaultCoords) == FALSE]
  ### Extract alternating lats and lons
  PolyLats <- FaultCoordsOnly[seq(1, length(FaultCoordsOnly), 2)]
  PolyLons <- FaultCoordsOnly[seq(2, length(FaultCoordsOnly), 2)]
  ### Convert lats/lons into lines
  ### NB: matrix needs c - if you don't explicitly hand it a list it will only use one dataset and just ignore the other
  Line <- Line(matrix(c(PolyLons, PolyLats), nrow=length(PolyLats)))
  ### Convert Polygon object to Polygons object, attach name of Polygon/fault
  L2 <- Lines(list(Line), ID=FaultName)
  ### Attach CRS to polygon
  L3 <- SpatialLines((list(L2)),  proj4string = CRS("+proj=longlat +datum=WGS84"))
  ###  - Should get this into projected CR so we are working in km/m rather than degrees - can use sp_transform but you'll need to figure out suitable projection
  return(L3)
}


#' Calculate on and off fault events and create a distance from fault map, which can be saved to specified directory
#'
#' This function identifies points that occur on and off-fault, and returns a distance from fault map.
#' @param Points spatial points dataframe of events/covering map area, should be in EPSG CRS.
#' @param FaultMap SpatialPolygonsDataFrame with fault polygons to be used. Should contain a column called `Fault` that names individual fault segments (numbers will also work - this is used to identify which fault an event is on)
#' @param Res Resolution of output in m. Default is 2500m.
#' @param saveDir Optional save directory for output, if supplied raster will be saved as rds
#' @param EPSG specify EPSG code, defaults to California 3310, see https://spatialreference.org/ref/epsg/ for other areas
#' @param xlim limits for x-axis of fault-distance raster
#' @param ylim limits for the y-axis of fault-distance raster
#' @return returns raster map of fault distances in metres
#' @keywords faultdistance
#' @export
#' @examples
#' Fd <- FaultDistMap(UCERFCat, FTDF2)
#'

FaultDistMap <- function(Points, FaultMap, Res=2500, saveDir = 0, EPSG=3310, xlim, ylim){
  PolyTest2 <- over(Points, FaultMap)
  Points$Fault <- PolyTest2$Name
  OnFault <- subset(Points, is.na(EQs$Fault) == FALSE)
  OffFault <- subset(Points, is.na(EQs$Fault) == TRUE)
  
  FLDF <- st_as_sf(FaultMap)
  FLDF <- st_transform(FLDF, EPSG)
  OF <- st_as_sf(OffFault)
  OF <- st_transform(OF, EPSG)
  
  dd2 <- gDistance(as(FLDF, "Spatial"), as(OF, "Spatial"), byid = TRUE)
  min_dist <- apply(dd2, 1, min)
  
  grdF3= expand.grid(x=seq(from=xlim[1], to= xlim[2], by=0.1), y=seq(from=ylim[1], to=ylim[2], by=0.1))
  coordinates(grdF3) = c("x", "y")
  proj4string(grdF3) = CRS(proj4string(FaultMap))
  # grid expand 'by' arguments don't seem to affect raster, so set resolution in raster layer.
  r2 = raster(grdF3, resolution= 0.05)
  #rm(grdF3)
  # ## Can potentially increase this resolution a bit but it takes a while
  r2 <- projectRaster(r2, crs=crs(FLDF), res=Res)
  dr2 <- gDistance(as(FLDF, "Spatial"), as(r2, "SpatialPoints"), byid = TRUE)
  r2[] <- apply(dr2, 1, min)
  #rm(dr2)
  #ggplot()+ gg(r2)+scale_fill_gradientn(colors=RColorBrewer::brewer.pal(9, "YlOrRd"))  +  ggtitle("Distance to fault (m) - resolution 2.5km")+ labs(x="", y="") + theme_classic() + coord_equal()
  # + gg(as(OF, "Spatial"))
  #r2km <- r2
  r2km <- setValues(r2, r2[]/1000)
  proj_rast <- projectRaster(r2km, crs=CRS(proj4string(FaultMap)))
  FaultDist <- as(proj_rast, 'SpatialPixelsDataFrame')
  if (saveDir != 0){
    saveRDS(FaultDist, file=paste0(saveDir,"FaultDist2pt5kmProj.rds"))
  }
  return(FaultDist)
}


#' Apply uniform buffer to fault polygons
#'
#' Simple function to apply uniform buffer to fault polygons using st_buffer, clean any intersections with clgeo_Clean, and crop large datasets to specifiec boundary.
#' @param Geom Fault geometry as SpatialPolygonsDataFrame
#' @param width width of buffer to be applied in metres
#' @param bdy boundary for cropping. Must be supplied as SpatialPolygon. Default is not to crop.
#' @param EPSG EPSG code for area of application. Defaults to California (3310), see https://spatialreference.org/ref/epsg/ for other areas
#' @return returns SpatialPolygonsDataFrame with specified buffer
#' @keywords buffer
#' @export
#' @examples
#' bdy <- SpatialPolygons(list(Polygons(list(Polygon(boundary[[1]]$loc)), ID="P1")), proj4string=CRS(proj4string(FTDF2)))
#' UnifBuff <- UniformBuffer(FTDF2, 5000, bdy=bdy)


UniformBuffer<- function(Geom, width, bdy = 0, EPSG = 3310){
  #FPr <- st_transform(Geom, EPSG)
  FBuff <- st_buffer(Geom, width)
  length(which(st_is_valid(FBuff) == FALSE))
  #FB <-  st_transform(FBuff, "+proj=longlat +datum=WGS84")
  FB <- as(Geom, "Spatial")
  
  #rm(QFBuff, QFB, QFPr, QFaults)
  ## have to clean up polygons to avoid one weird self-intersection which doesn't show up in sf
  FBc <- clgeo_Clean(FB)
  ### If you've skipped the cleaning step, this may crash Rstudio if you have any invalid geometry. Don't say I didn't warn you!
  ### But actually 'crop' just seems to throw an error, whereas gIntersection just killed everything.
  ### Can also use boundary polygon if you have one, but I'm not sure where this would be necessary?
  if (bdy != 0){
    FBc <- crop(FBc, bdy)
  }
  return(FBc)
}

#' Dip-dependent fault polygon buffer function
#'
#' Apply dip-dependent buffers to fault polygons. Will apply dip in metres, so check CRS!
#' Uses UCERF3 buffers: fault buffer scales with fault dip - 0km at 50 dip to 12km at 90 dip
#' @param PolygonsData a dataset of polygons stored as a spatial polygons object.
#' @param Dip a vector of dips. Must be the same length as PolygonsData
#' @return returns new polygons, buffered according to dip and stored as a SpatialPolygonsDataFrame
#' @keywords buffers
#' @export
#' @examples
#' FaultPolyBuffer <- BufferFunction(FL_T, FL_T$Dip)
#' Use st_transform to transform to EPSG (3310 for state of California, see https://spatialreference.org/ref/epsg/ and vice-versa

DipBufferFunction <- function(PolygonsData, Dip){
  DipAng <- vector(mode="numeric", length=nrow(PolygonsData))
  BufferSize <- vector(mode="numeric", length=nrow(PolygonsData))
  # Units of dist are m
  #BufferSize <- DipAng*116.7 - IF you start buffering at 0 (degrees and km) and extend to 15km at 90 degrees
  ## Here start the buffer if dip > 50 and buffer to 12km at dip=90
  ## But something really dodgy happens if buffersize=0 so set to 0.1 (that's like 10cm, should not make any difference to data!)
  for (i in 1:nrow(PolygonsData)){
    DipAng[i] <- Dip[i]
    if (DipAng[i] > 50){
      BufferSize[i] <- (DipAng[i]-50)*0.3
    }
    else {BufferSize[i] = 0.1}
  }
  #print(BufferSize)
  BufferPoly <- st_buffer(PolygonsData, dist=BufferSize*1000)
  return(BufferPoly)
}


#' Make a SpatialPolygons object from a dataframe 
#'
#' This function converts from a dataframe of x and y coordinates to a SpatialPolygons object
#' This is basically a nested series of conversions from a dataframe to a Polygon to a Polygons list to a SpatialPolygons object.
#' @param coords a dataframe of x and y coordinates describing a spatial polygon
#' @param crs desired crs for the object
#' @return returns SpatialPolygons object
#' @keywords SpatialPolygons
#' @export
#' @examples
#' RELM_poly <- read.delim("RELMTestingPolygon.txt", sep="", header = FALSE)
#' RELM <- spatial_poly_from_df(RELM_poly, CRS(SRS_string = wkt))

spatial_poly_from_df <- function(coords, crs){
  SpatialPolygons((list( Polygons(list(Polygon(coords)), ID="A1"))),  proj4string = crs)
  }


#' Project inlabru model output to uniform grid for CSEP testing 
#'
#' This function takes in an inlabru fit object and model description and parameters required to generate a pyCSEP gridded-forecast
#' @param lgcp_fit an inlabru fit object
#' @param lgcp_model a description of an inlabru model
#' @param b_poly_latlon a boundary polygon in lat-lon coordinates for pyCSEP testing
#' @param dh grid spacing for pyCSEP forecast (in lat/lon)
#' @param mag_min minimum magnitude for the forecast
#' @param mag_max maximum magnitude in forecast 
#' @param b_est estimated b-value to be used in forecast
#' @param mesh inlabru mesh object used for the fitted model
#' @return returns a data frame consistent with those required by pyCSEP
#' @export
csep_grid_wrapper <- function(lgcp_fit, lgcp_model, b_poly_latlon, dh, mag_min, mag_max, b_est, mesh){
  ## Set up magnitude breaks
  mag.break <- seq(mag_min, mag_max, by=0.1)
  ## Find bbox of boundary polygon
  b <- as.data.frame(bbox(b_poly_latlon))
  ### Make uniform grid based on limits + 0.5*spacing to get midpoints expected in RELM polygon
  mid_grd <- expand.grid(x=seq(from=(b$min[1]-(0.5*dh)), to= (b$max[1]+(0.5*dh)), by=dh), y=seq(from=(b$min[2]-(0.5*dh)), to=(b$max[2] +  (0.5*dh)), by=dh))
  
  # Make spatial object
  coordinates(mid_grd) = ~x+y
  proj4string(mid_grd) <- CRS(proj4string(b_poly_latlon))
  # Only keep midpoints in the polygon
  pts_grd_latlon <- crop(mid_grd, b_poly_latlon)
  pts_grd_km <- spTransform(pts_grd_latlon, crs_Cal_km)
  
  ## Predict at mesh locations
  Pr_Num <- predict(lgcp_fit, mesh, lgcp_model)
  ## Project to  grid
  proj <- INLA::inla.mesh.project(mesh, pts_grd_km)
  ## Get values at grid midpoints
  N <- as.vector(proj$A %*% Pr_Num$mean)
  
  ## a value from N scaled by area of grid cells
  a2 <- log(N*(1/(dh^2)))/log(10)
  ## number of bins we want
  n_mb <- length(mag.break)
  
  ## Distribute events/bin over different magnitude bins according to GR
  bmt <- vector()
  for(i in 1:length(pts_grd_km)){
    Ns <-  a2[i] - b_est*(mag.break - 4.95)
    M_D <- c(abs(diff(10^Ns)), 10^Ns[n_mb])
    bmt <- c(bmt, M_D)
  }
  
  ## lower latitude values for grid cells
  lats_1 <- rep((as.numeric(pts_grd_latlon$y) - 0.5*dh), each=n_mb)
  ## Upper values
  lats_2 <- lats_1 + dh
  
  ## Same for lons
  lons_1 <- rep((as.numeric(pts_grd_latlon$x) - 0.5*dh), each=n_mb)
  lons_2 <- lons_1 + dh
  
  ## Set up depth bins (not currently used)
  D1 <- rep(0, length(lats_2))
  D2 <- rep(30, length(lats_2))
  
  ## Set up magnitude bins (upper bounds)
  M2 <- mag.break[2:n_mb]
  M2[n_mb] <- 10
  
  mags_1 <- rep(mag.break, length(pts_grd_latlon$x))
  mags_2 <- rep(M2, length(pts_grd_latlon$x))
  
  ## RELM flag, completely unused but expected in input.
  flag <- rep(1, length(lats_1))
  
  ### Make dataframe of results
  DFT_bm <- as.data.frame(cbind(lons_1, lons_2, lats_1, lats_2, D1, D2, mags_1, mags_2, bmt, flag))
  
  # Sort by lon (probably no longer necessary!)
  grd_forecast <- DFT_bm %>% arrange(lons_1)
  return(grd_forecast)
}

#' Generate catalogue-type forecast from inlabru fitted intensity posterior sample
#'
#' This function takes in posterior intensity sample from the inlabru `generate' function (or a mean/mode/percentile if you feel so inclined)
#' This sample is then used to create an intensity field from which we sample points with a rejection sampler
#' @param loglambda log intensity values from a generate call or a fitted model object
#' @param bdy a boundary polygon over which to generate samples
#' @param mesh inlabru mesh object used for the fitted model
#' @param crs coordinate reference system to be used - must be consistent with mesh, bdy and loglambda
#' @param n_events number of events to sample for each catalogue. Currently assumes this is a Poisson rate and randomly samples exact number
#' @param b_val a b-value to be used in assigning magnitudes
#' @param m_min minimum magnitude in forecast model, used in magnitude creation by TapGRM
#' @return returns a data frame of a synthetic catalogue
#' @export
point_sampler <- function(loglambda, bdy, mesh,  crs, num_events, b_val, m_min){
  ## Number of events for single catalogue from a poisson distribution with lambda = num_events
  cat_n <- rpois(1, num_events)
  ## Find maximum log intensity, this will be used in the rejection sampler
  loglambda_max <- max(loglambda)
  
  ## Set up a spatialpoints dataframe for our results
  samp.points <- SpatialPoints(data.frame(x = 0, y = 0))
  samp.points$mags <- 0
  samp.points <- samp.points[-1,]
  proj4string(samp.points) <- crs
  
  ## Number of points to try sampling
  num <- 0
  n1 <- 300000
  n2 <- 5000000
  ## To sample the correct number of points, keep going until the num >= cat_n
  while (num < cat_n){
    pts <- spsample(bdy, n1, "random")
    #pts <- spTransform(points, crs)
    proj <- INLA::inla.mesh.project(mesh, pts)
    ## calculate ratio of lambda to total lambda
    lambda_ratio <- exp(as.vector(proj$A %*% loglambda) - loglambda_max)
    ## Keep points where the ratio is greater than random variable (0-1)
    keep <- proj$ok & (runif(n1) <= lambda_ratio)
    kept <- pts[keep]
    
    ## If we didn't get any events, run again with more sampled points
    ## Unfortunately this does sometimes happen, rejection samplers are inefficient!
    while (length(kept) == 0){
      
      pts <- spsample(bdy, n2, "random")
      
      proj <- INLA::inla.mesh.project(mesh, pts)
      lambda_ratio <- exp(as.vector(proj$A %*% loglambda) - loglambda_max)
      keep <- proj$ok & (runif(n2) <= lambda_ratio)
      kept <- pts[keep]
    }
    kept$mags <- rep(0, length(kept))
    
    ## Add any new points to existing samp.points
    samp.points <- rbind(samp.points, kept)
    num <- length(samp.points)
    
  }
  
  ## Keep exactly cat_n points, choose these randomly from all of the points we've kept so far
  kp <- sample(seq(1, length(samp.points), by=1), cat_n, replace=FALSE)
  samp.points <- samp.points[kp,]
  ## Get magnitudes for this catalogue from tapered GR distribution
  ## Assumes corner magnitude = 8, should update to set this in input params
  samp.points$mags <- TapGRM(cat_n, b_val, 8, m_min)
  
  return(samp.points)
}


#' Sample magnitudes from a tapered Gutenberg-Richter distribution
#'
#' This function takes a number of required magnitudes, a b-value and corner- and minimum-magnitudes and returns sampled magnitudes consistent with a tapered-GR distribution 
#' Uses method of Vere-Jones et al 2001 (https://doi.org/10.1046/j.1365-246x.2001.01348.x)
#' @param n number of magnitudes to sample
#' @param b b-value for magnitude sampling
#' @param corner_mag a corner magnitude value
#' @param m_min minimum magnitude for magnitude creation
#' @return returns n magnitudes from a TGR distribution
#' @export
TapGRM <- function(n, b, corner_mag, m_min){
  ## Convert Mc to moment and b to beta
  Mt <- 10**((3/2)*m_min + 9.1)
  Mc <- 10**((3/2)*corner_mag + 9.1)
  beta <- (2/3)*b
  
  R1 <- runif(n, 0, 1)
  R2 <- runif(n, 0, 1)
  
  M1 <- Mt*R1**(-1/beta)
  M2 <- Mt - Mc*log(R2)
  
  # Pick minimum value of two sampled options
  mom <- pmin(M1, M2)
  
  ## convert moments back to magnitudes
  ms <- (2/3)*(log10(mom) - 9.1)
  
  return(ms)
}
