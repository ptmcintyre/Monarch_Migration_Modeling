get_model_extent_from_data <- function(element_occurrences, hucs, add_hucs = NULL, remove_hucs = NULL){
  element_occurrences_buffered <- st_buffer(element_occurrences, dist = 5000) %>% 
    st_transform(crs = st_crs(hucs))
  eo_huc_intersection <- st_intersects(hucs, element_occurrences_buffered, sparse = TRUE)
  model_extent <- hucs[which((purrr::map(eo_huc_intersection, length) %>% unlist()) > 0), ]
  if (!is.null(hucs_to_remove)){
    model_extent <- model_extent %>% 
      dplyr::filter(!(HUC10 %in% remove_hucs))
  }
  if (!is.null(hucs_to_add)){
    model_extent <- model_extent %>% 
      rbind(hucs %>% dplyr::filter(HUC10 %in% add_hucs))
  }
  return(model_extent)
}

refine_model_extent_from_MIRTL <- function(cutecode = cutecode, current_extent = species_hucs, hucs){
  
  hucs_to_add <- NULL
  hucs_to_remove <- NULL
  
  polygons_to_add <- googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/17Sy7ujWzAQWA9XRcA1xg8n4BPvg50oIp-z9dE_po8RU/edit#gid=0") %>% 
    dplyr::filter(species == cutecode, map == "extent"| map == "species", shape_type == "drawn Polygon", suggested_action == "add to extent")
  
  if (nrow(polygons_to_add) > 0){
    hucs_to_add <- purrr::map(1:nrow(polygons_to_add), function(poly){
      pol <- st_polygon(
        list(
          cbind(
            polygons_to_add[poly, ]$longitude %>% unlist() %>% strsplit(";") %>% unlist() %>% as.numeric(), 
            polygons_to_add[poly, ]$latitude %>% unlist() %>% strsplit(";") %>% unlist() %>% as.numeric()
          )
        )
      )
      pol <- st_sfc(pol, crs=4326)
      pol <- pol %>% st_transform(crs = st_crs(hucs))
      pol_huc_intersection <- st_intersects(hucs, pol, sparse = TRUE)
      new_hucs <- hucs[which((purrr::map(pol_huc_intersection, length) %>% unlist()) > 0), ]
      new_hucs
    }) %>% 
      bind_rows() %>% 
      distinct(.keep_all = TRUE)
  }
  
  polygons_to_remove <- read_sheet("https://docs.google.com/spreadsheets/d/17Sy7ujWzAQWA9XRcA1xg8n4BPvg50oIp-z9dE_po8RU/edit#gid=0") %>% 
    dplyr::filter(species == cutecode, map == "extent", shape_type == "drawn Polygon", suggested_action == "remove from extent")
  
  if (nrow(polygons_to_remove) > 0){
    hucs_to_remove <- purrr::map(1:nrow(polygons_to_remove), function(poly){
      pol <- st_polygon(
        list(
          cbind(
            polygons_to_remove[poly, ]$longitude %>% unlist() %>% strsplit(";") %>% unlist() %>% as.numeric(), 
            polygons_to_remove[poly, ]$latitude %>% unlist() %>% strsplit(";") %>% unlist() %>% as.numeric()
          )
        )
      )
      pol <- st_sfc(pol, crs=4326)
      pol <- pol %>% st_transform(crs = st_crs(hucs))
      pol_huc_intersection <- st_intersects(hucs, pol, sparse = TRUE)
      new_hucs <- hucs[which((purrr::map(pol_huc_intersection, length) %>% unlist()) > 0), ]
      new_hucs
    }) %>% 
      bind_rows() %>% 
      distinct(.keep_all = TRUE)
  }
  
  final_hucs <- c(setdiff(current_extent$HUC10, hucs_to_remove$HUC10), hucs_to_add$HUC10) %>% unique()
  
  final_extent <- hucs %>% dplyr::filter(HUC10 %in% final_hucs)
  
  return(final_extent)
}

map_model_extent <- function(){
  tmap_mode("plot")
  
  states <- USAboundaries::us_states() %>% 
    dplyr::filter(!(state_name %in% c("Hawaii", "Alaska", "Puerto Rico"))) %>% 
    st_set_crs(4326)
  
  #' #### Panel A: Model predictions with inset
  set_defaults(map_service = "esri", map_type = "world_imagery")
  #' ##### Model predictions
  model_basemap <- basemap_raster(rbind(species_hucs, huc10_layer[huc10_layer$HUC10 %in% hucs_to_remove, ]) %>% st_bbox() %>% st_as_sfc())
  
  model_extent_map = 
    tm_shape(model_basemap) + 
    tm_rgb(alpha = 0.9) +
    tm_layout(frame = FALSE, legend.show = TRUE, legend.text.color = "black", legend.title.color = "black", bg.color = "transparent", legend.text.size = 0.6, legend.title.size = 0.6, legend.outside = TRUE, legend.outside.position = "bottom") +
    tm_shape(shp = species_hucs) + tm_polygons(col = grey(0.5), alpha = 1, zindex = 2000, title.col = "Model extent of current model version")
  
  if (length(hucs_to_add) > 0){
    model_extent_map <- model_extent_map + tm_shape(shp = huc10_layer[huc10_layer$HUC10 %in% hucs_to_add, ]) + tm_polygons(col = "deepskyblue", alpha = 1, zindex = 4000, title.col = "Areas added to previous model extent")
  }
  if (length(hucs_to_remove) > 0){
    model_extent_map <- model_extent_map +   tm_shape(shp = huc10_layer[huc10_layer$HUC10 %in% hucs_to_remove, ]) + tm_polygons(col = "tomato", alpha = 1, zindex = 3000, title.col = "Areas removed from previous model extent")
  }
  model_extent_map <- model_extent_map +
    tm_add_legend("fill", labels = c("Areas kept from previous model extent", "Areas added to previous model extent", "Areas removed from previous model extent"), col = c(grey(0.5), "deepskyblue", "tomato"))
  return(model_extent_map)  
}

group_observations_to_polygons <- function(presence_observations = presence_observations){
  # use Accuracy column to buffer
  ptd_pol <- st_buffer(presence_observations, dist = presence_observations$accuracy)
  # merge overlapping polys
  ptd_pol_m <- st_union(ptd_pol, by_feature = FALSE)
  # convert to single-part, add attribute table
  ptd_pol_m = st_cast(ptd_pol_m, "POLYGON")
  ptd_df <- data.frame("pol_id" = c(1:length(ptd_pol_m)))
  ptd_pol_m <- st_set_geometry(ptd_df, ptd_pol_m)
  
  # get back date (using the most recent obs date for each polygon)
  ptp_j <- st_join(ptd_pol_m, ptd_pol[,c("obsdate","species_cd")], join = st_intersects, left = TRUE)
  ptp_a <- aggregate(ptp_j, by = list(ptp_j$pol_id), FUN = function(x){ max(x)})
  geomName <- attr(ptp_a, "sf_column")
  ptp <- ptp_a[,c("pol_id","obsdate","species_cd",geomName)]
  ptp <- cbind(ptp, RA = "very high")
  
  # Buffer by those units (m) to get group based on separation distance ("EO"), merge them
  # 500 meter default (means 1km sep distance)
  pt_buff1km <- st_buffer(ptp, 500)
  pt_grps <- st_union(pt_buff1km, by_feature = FALSE)
  
  # convert to single-part, add attribute table
  pt_grps = st_cast(pt_grps, "POLYGON")
  ptdf <- data.frame("group_id" = c(1:length(pt_grps)))
  pt_grps <- st_set_geometry(ptdf, pt_grps)
  
  # attribute points with group ID
  ptd_grps <- suppressWarnings(st_intersection(ptp, pt_grps))
  
  return(ptd_grps)
}

get_points_from_presence_polygons <- function(presPolys = presence_polygons, presObs = NULL) {
  
  #check for proper column names. If no error from next code block, then good to go
  #presPolys$RA <- presPolys$SFRACalc
  desiredCols <- c("UID", "GROUP_ID", "SPECIES_CD", "RA", "OBSDATE")
  
  presPolys$GROUP_ID <- as.numeric(factor(presPolys$EO_SEQ_UID, levels = unique(presPolys$EO_SEQ_UID)))
  presPolys$SPECIES_CD <- paste0(substr(presence_polygons$GNAME, 1, 1), substr(strsplit(presence_polygons$GNAME[1], " ")[[1]][2], 1, 4), "_", presence_polygons$EGT_ID)
  presPolys$UID <- presPolys$EO_SEQ_UID
  presPolys <- presPolys %>% 
    dplyr::mutate(
      UID = EO_SEQ_UID,
      GROUP_ID = seq(EO_SEQ_UID),
      SPECIES_CD = taxon_code,
      RA = ifelse(!is.na(EST_REP_AC), EST_REP_AC, "Medium"),
      OBSDATE = LASTOBS_D
    )
  
  # if ("FALSE" %in% c(desiredCols %in% names(presPolys))) {
  #   stop(paste0(
  #     "Column(s) are missing or incorrectly named: ",
  #     paste(desiredCols[!desiredCols %in% shpColNms], collapse = ", ")
  #   ))
  # } else {
  #   print("Required columns are present")
  # }
  # 
  # if (any(is.na(presPolys[, c("UID", "GROUP_ID", "SPECIES_CD", "RA")]))) {
  #   stop(
  #     "The columns 'UID','GROUP_ID','SPECIES_CD', and 'RA' (SFRACalc) cannot have NA values."
  #   )
  # }
  
  #pare down columns
  presPolys <- presPolys[, desiredCols]
  
  ## Merge with presence_observation_polygons, if provided
  if (!is.null(presObs)){
    names(presObs)[!names(presObs) == "geometry"] <- toupper(names(presObs)[!names(presObs) == "geometry"])
    names(presObs)[1] <- "UID"
    presObs <- presObs[, desiredCols]
    presObs$SPECIES_CD <- presPolys$SPECIES_CD[1]
    presObs$GROUP_ID <- presObs$GROUP_ID + max(presPolys$GROUP_ID)
    presObs$OBSDATE <- paste(purrr::map(strsplit(presObs$OBSDATE, "/"), 3), purrr::map(strsplit(presObs$OBSDATE, "/"), 1), purrr::map(strsplit(presObs$OBSDATE, "/"), 2), sep = "-") %>% as.Date() %>% as.character()
    
    presPolys <- rbind(presPolys, presObs)
  }
  
  # set date/year column to [nearest] year, rounding when day is given
  presPolys$OBSDATE <- as.character(presPolys$OBSDATE)
  presPolys$date <- NA
  for (d in 1:length(presPolys$OBSDATE)) {
    dt <- NA
    do <- presPolys$OBSDATE[d]
    if (grepl("^[0-9]{4}.{0,3}$", do)) {
      # year only formats
      dt <- as.numeric(substring(do, 1, 4))
    } else {
      if (grepl("^[0-9]{4}[-|/][0-9]{1,2}[-|/][0-9]{1,2}", do)) {
        # ymd formats
        try(dt <- as.Date(do), silent = TRUE)
      } else if (grepl("^[0-9]{1,2}[-|/][0-9]{1,2}[-|/][0-9]{4}", do)) {
        # mdy formats
        try(dt <- as.Date(do, format = "%m/%d/%Y"), silent = TRUE)
      }
      # if still no match, or if failed
      if (is.na(dt)) {
        if (grepl("[0-9]{4}", do)) {
          # use first 4-digit sequence as year
          dt <- regmatches(do, regexpr("[0-9]{4}", do))
          if (as.integer(dt) < 1900 |
              as.integer(dt) > format(Sys.Date(), "%Y")) {
            # years before 1900 and after current year get discarded
            dt <- Sys.Date()
          } else {
            dt <- as.Date(paste0(dt, "-01-01"))
          }
        }
        # put additional date formats here
      }
      # give up and assign current date
      if (is.na(dt)) {
        dt <- Sys.Date()
      }
      dt <-
        round(as.numeric(format(dt, "%Y")) + (as.numeric(format(dt, "%j")) / 365.25))
    }
    presPolys$date[d] <- dt
  }
  desiredCols <- c(desiredCols, "date")
  
  # What's the name of the geometry column? Change it to geometry if needed
  # presPolys <- backu
  # backu <- presPolys
  
  if (!attr(presPolys, "sf_column") == "geometry") {
    oldName <- attr(presPolys, "sf_column")
    #add new col
    names(presPolys)[names(presPolys) == oldName] <- "geometry"
    attr(presPolys, "sf_column") = "geometry"
  }
  
  # explode multi-part polys ----
  # suppressWarnings(shp_expl <- st_cast(presPolys, "POLYGON"))
  
  shp_expl <- presPolys
  
  #add some columns (explode id and area)
  shp_expl <- cbind(shp_expl,
                    EXPL_ID = 1:nrow(shp_expl),
                    AREAM2 = st_area(shp_expl))
  names(shp_expl$geometry) <- NULL # temp fix until sf patches error
  
  # just in case convert to lower
  names(shp_expl) <- tolower(names(shp_expl))
  
  #calculate Number of points for each poly, stick into new field
  shp_expl$PolySampNum <-
    round(400 * ((2 / (1 + exp(
      -(as.numeric(shp_expl$aream2) / 900 + 1) * 0.004
    ))) - 1))
  #make a new field for the design, providing a stratum name
  shp_expl <-
    cbind(shp_expl, "stratum" = paste("poly_", shp_expl$expl_id, sep = ""))
  
  # sample must be equal or larger than the RA sample size in the random forest model
  shp_expl$ra <- factor(tolower(as.character(shp_expl$ra)))
  
  # QC step: are any records attributed with values other than these?
  raLevels <- c("very high", "high", "medium", "low", "very low")
  if ("FALSE" %in% c(shp_expl$ra %in% raLevels)) {
    stop("at least one record is not attributed with RA appropriately")
  } else {
    print("RA levels attributed correctly")
  }
  
  #EObyRA <- unique(shp_expl[,c("expl_id", "group_id","ra")])
  shp_expl$minSamps[shp_expl$ra == "very high"] <- 5
  shp_expl$minSamps[shp_expl$ra == "high"] <- 4
  shp_expl$minSamps[shp_expl$ra == "medium"] <- 3
  shp_expl$minSamps[shp_expl$ra == "low"] <- 2
  shp_expl$minSamps[shp_expl$ra == "very low"] <- 1
  
  shp_expl$finalSampNum <-
    ifelse(
      shp_expl$PolySampNum < shp_expl$minSamps,
      shp_expl$minSamps,
      shp_expl$PolySampNum
    )
  
  ranPts <- st_sample(shp_expl, size = shp_expl$finalSampNum * 2)
  ranPts.sf <- st_sf(ranPts)
  names(ranPts.sf) <- "geometry"
  st_geometry(ranPts.sf) <- "geometry"
  
  ranPts.joined <- st_join(ranPts.sf, shp_expl)
  
  # check for polys that didn't get any points
  polysWithNoPoints <-
    shp_expl[!shp_expl$expl_id %in% ranPts.joined$expl_id, ]
  if (nrow(polysWithNoPoints) > 0) {
    stop("One or more polygons didn't get any points placed in them.")
  }
  
  #  remove extras using straight table work
  
  # this randomly assigns digits to each point by group (stratum) then next row only takes
  # members in group that are less than target number of points
  rndid <-
    with(ranPts.joined, ave(
      expl_id,
      stratum,
      FUN = function(x) {
        sample.int(length(x))
      }
    ))
  ranPts.joined2 <-
    ranPts.joined[rndid <= ranPts.joined$finalSampNum, ]
  
  ranPts.joined <- ranPts.joined2
  rm(rndid, ranPts.joined2)
  
  #check for cases where sample smaller than requested
  # how many points actually generated?
  ptCount <- table(ranPts.joined$expl_id)
  targCount <- shp_expl[order(shp_expl$expl_id), "finalSampNum"]
  overUnderSampled <- ptCount - targCount$finalSampNum
  
  #positive vals are oversamples, negative vals are undersamples
  print(table(overUnderSampled))
  # If you get large negative values then there are many undersamples and
  # exploration might be worthwhile
  
  names(ranPts.joined) <- tolower(names(ranPts.joined))
  
  colsToKeep <- c("stratum", tolower(desiredCols))
  ranPts.joined <- ranPts.joined[, colsToKeep]
  
  return(ranPts.joined)
  
}

generate_aquatic_occurrence_inputs <- function(mod_ext = model_extent,
                                               presence_polys = presence_polygons,
                                               presence_observation_polys = presence_observation_polygons,
                                               max_stream_level = 4
){
  
  ## Load and process NHDFlowline network data
  ### Load all stream reaches within the model extent's bounding box for efficiency
  model_reaches <- st_read("S://aquaticEnvVars_NHD/US_flowline_network.gdb", layer = "NHDFlowline_Network") %>% 
    st_zm() %>% 
    st_crop(y = st_bbox(mod_ext))
  ### Keep only stream reaches that overlap the model extent
  model_reaches <- model_reaches[st_intersects(model_extent, model_reaches %>% dplyr::select(COMID), sparse = TRUE) %>% unlist(), ]
  
  ## Identify presence reaches
  ### From EO/SF data
  presence_reaches_IDs <- st_intersects(presence_polys, model_reaches %>% dplyr::select(COMID), sparse = TRUE) %>% unlist()
  ### From presence observations
  presence_observations_reaches_IDs <- st_intersects(presence_observation_polys %>% st_buffer(dist = 10), model_reaches %>% dplyr::select(COMID), sparse = TRUE) %>% unlist()
  ### Combine all presence reaches
  presence_reaches_IDs <- c(presence_reaches_IDs, presence_observations_reaches_IDs)
  ## Isolate presence reaches
  presence_reaches <- model_reaches[presence_reaches_IDs, ] %>%
    dplyr::distinct(COMID, .keep_all = TRUE)
  
  ## Group presence reaches to identify validation groups
  ### Buffer by EO separation distance (500m equals the default separation distance of 1km)
  presence_reaches_buffered <- st_buffer(presence_reaches, dist = 500)
  ### Merge overlapping polygons
  presence_reaches_buffered_m <- st_union(presence_reaches_buffered, by_feature = FALSE)
  ### Convert to single-part polygons
  presence_reaches_buffered_m <- st_cast(presence_reaches_buffered_m, "POLYGON")
  ### Add attribute table
  presence_reaches_buffered_m_df <- data.frame("pol_id" = c(1:length(presence_reaches_buffered_m)))
  presence_reaches_buffered_m <- st_set_geometry(presence_reaches_buffered_m_df, presence_reaches_buffered_m)
  ### Aggregate presence reaches that intersect with each other 
  presence_reaches_aggregated <- aggregate(presence_reaches_buffered_m, by = list(presence_reaches_buffered_m_df$pol_id), FUN = function(x){ max(x)})
  ### Add group_id for the identified stream reach group to presence_reaches object
  presence_reaches$group_id <- st_intersects(presence_reaches, presence_reaches_aggregated, sparse = TRUE) %>% unlist()
  
  ## Filter by maximum stream level
  presence_reaches <- presence_reaches %>% 
    dplyr::filter(StreamLeve <= max_stream_level)
  
  ## Get background reaches
  ### Identify reaches adjacent to presence reaches
  presence_adjacent_reaches <- st_intersects(presence_reaches %>% st_buffer(dist = 500), model_reaches %>% dplyr::select(COMID), sparse = TRUE) %>% unlist()
  ### Set background reaches
  background_reaches <- model_reaches[-presence_adjacent_reaches, ] 
  ### Sample subset of background reaches if too many
  if (nrow(background_reaches) > (nrow(presence_reaches) * 10)){ #PJM changed to 12
    background_reaches <- background_reaches %>% 
      dplyr::sample_n(nrow(presence_reaches) * 10) %>% #PJM changed to 12
      dplyr::distinct(COMID, .keep_all = TRUE) 
  }
  
  out <- list(presence_reaches = presence_reaches,
              background_reaches = background_reaches,
              model_reaches = model_reaches
  )
  
  return(out)
}

plot_model_inputs <- function(presence_polys = presence_polygons, 
                              presence_pts = presence_points,
                              background_pts = background_points,
                              mod_extent = model_extent,
                              presence_observation_polys = NULL,
                              proportion_background_points = 0.1
){
  
  tmap_mode("view")
  
  m <- tm_basemap(leaflet::providers$Esri.WorldTopoMap, alpha = 0.7) +
    tm_shape(mod_extent) + tm_polygons(col = grey(0.5), alpha = 0.4, border.col = "black", lwd = 2) +
    tm_shape(sample_n(background_pts, nrow(background_pts) * proportion_background_points)) + tm_dots(col = "tomato") +
    tm_shape(presence_polys) + tm_borders(col = "blue") +
    tm_shape(presence_pts) + tm_dots(col = "blue") +
    tm_layout()
  
  if (!is.null(presence_observation_polys)) m <- m + tm_shape(presence_observation_polys) + tm_borders(col = "blue")
  
  # print(m)
  
  return(m)
  
}

run_models <- function(presence_pts = presence_points,
                       background_pts = background_points,
                       species_prdictrs = env_predictor_set,
                       algorithms = c("rf", "xgb")
){
  
  # algorithms <- match.arg(algorithms)
  
  # Clean up input data
  presence_pts <- cbind(presence_pts, pres=1)
  #presence_pts$group<-as.character(presence_pts$group) #pjm added
  background_pts$stratum <- "pseu-a"
  background_pts <- cbind(background_pts, group_id="pseu-a", pres=0, ra="high", species_cd="background")
  background_pts$group_id <- factor(background_pts$group_id)
  
  # Bind presence and background points
  colList <- c("species_cd", "group_id", "pres", "stratum", "ra", species_prdictrs)
  df.full <- rbind(presence_pts %>% dplyr::select(colList), background_pts %>% dplyr::select(colList)) %>% 
    st_set_geometry(value = NULL) %>% 
    as.data.frame() %>% 
    dplyr::mutate(
      stratum = factor(stratum),
      group_id = factor(group_id),
      pres = factor(pres),
      ra = factor(ra),
      species_cd = factor(species_cd)
    ) %>% 
    dplyr::filter(complete.cases(.))
  
  #df.full$ra<-"high" #PJM hack to see if this fixes issue with sampsize
  
  # Set up model validation routine
  ## Get number of polygons
  numPys <-  nrow(table(presence_pts$stratum))
  ## Get number of EOs
  numEOs <- nrow(table(presence_pts$group_id))
  ## Set maximum number of EOs to jacknife over
  eoCountLimit <- 20
  ## Make sample size groupings
  group <- vector("list")
  group$colNm <- ifelse(numEOs < eoCountLimit,"stratum","group_id")
  group$JackknType <- ifelse(numEOs < eoCountLimit,"polygon","spatial grouping")
  if(numEOs < eoCountLimit) { #pjm changed from < to >
    group$vals <- unique(presence_pts$stratum)
  } else {
    group$vals <- unique(presence_pts$group_id)
  }
  ## Set sample size by EO accuracy
  EObyRA <- unique(df.full[,c(group$colNm,"ra")])
  EObyRA$sampSize[EObyRA$ra == "very high"] <- 5 #PJM changed from 5
  EObyRA$sampSize[EObyRA$ra == "high"] <- 4 #PJM changed from 4
  EObyRA$sampSize[EObyRA$ra == "medium"] <- 3  #PJM changed from 3
  EObyRA$sampSize[EObyRA$ra == "low"] <- 2
  EObyRA$sampSize[EObyRA$ra == "very low"] <- 1
  ## Handle cases where more than one RA is assigned per EO. Handle it here by taking max value
  EObySS <- aggregate(EObyRA$sampSize, by=list(EObyRA[,group$colNm]), max)
  # Set the background pts to the sum of the EO samples
  names(EObySS) <- c(group$colNm,"sampSize")
  EObySS$sampSize[EObySS[group$colNm] == "pseu-a"] <- sum(EObySS[!EObySS[group$colNm] == "pseu-a", "sampSize"])
  # Set sample sizes
  sampSizeVec <- EObySS$sampSize
  names(sampSizeVec) <- as.character(EObySS[,group$colNm])
  rm(EObySS, EObyRA)
  # Reset sample sizes to number of points, when it is smaller than desired sample size
  # This is only relevant when complete.cases may have removed some points from an already-small set of points
  totPts <- table(df.full[,group$colNm])
  totPts[names(totPts)%in%names(sampSizeVec)]
  totPts<-totPts[totPts>0] #PJM added to solve wierd index problem with extra groups
  for (i in names(sampSizeVec)) if (sampSizeVec[i] > totPts[i]) sampSizeVec[i] <- totPts[i]
  rm(totPts)
  # if colList gets modified, 
  # also modify the locations for the independent and dependent variables, here
  depVarCol <- 3
  indVarCols <- c(6:length(colList))
  
  # Run models
  
  ## Prepare output
  model_output <- vector("list", length = length(algorithms))
  names(model_output) <- algorithms
  
  ## Random Forest
  if ("rf" %in% algorithms){
    
    algo <- "rf"
    rf_model_version <- paste(presence_pts$species_cd[1], algo, gsub("-", "", Sys.Date() %>% as.character()), sep = "_")
    rf_model_creation_date <- Sys.Date() %>% as.character()
    
    # tune mtry ----
    #subset df to ten times pres for getting mtry  ## to increase speed
    rowCounts <- table(df.full$pres)
    if(rowCounts["0"] > (10 * rowCounts["1"])){
      tunePres <- df.full[df.full$pres == 1, ]
      tuneAbs <- df.full[df.full$pres == 0, ]
      tuneAbs <- tuneAbs[sample(nrow(tuneAbs), nrow(tunePres)* 10),] ##PJM changed to 12
      df.tune <- rbind(tunePres, tuneAbs)
      rm(tuneAbs, tunePres)
    } else {
      df.tune <- df.full
    }
    # run through mtry twice
    # very small numbers of polys can make a perfect prediction, especially
    # when subsetted like above. If error (prediction error == 0) results from 
    # subsetting, catch error and try tuning with full set
    mtry <- tryCatch(
      {
        x <- tuneRF(df.tune[,indVarCols],
                    y=df.tune[,depVarCol],
                    ntreeTry = 300, stepFactor = 2, mtryStart = 6,
                    strata = df.full[,group$colNm], sampsize = sampSizeVec, replace = TRUE)
        newTry <- x[x[,2] == min(x[,2]),1]
        y <- tuneRF(df.tune[,indVarCols],
                    y=df.tune[,depVarCol],
                    ntreeTry = 300, stepFactor = 1.5, mtryStart = max(newTry),
                    strata = df.full[,group$colNm], sampsize = sampSizeVec, replace = TRUE)
        
        mtry <- max(y[y[,2] == min(y[,2]),1])
        rm(x,y, df.tune, newTry)
        mtry
      }, 
      error=function(cond) {
        message("Can't tune with subset, using full set.")
        df.tune <- df.full
        x <- tuneRF(df.tune[,indVarCols],
                    y=df.tune[,depVarCol],
                    ntreeTry = 300, stepFactor = 2, mtryStart = 6,
                    strata = df.full[,group$colNm], sampsize = sampSizeVec, replace = TRUE)
        newTry <- x[x[,2] == min(x[,2]),1]
        y <- tuneRF(df.tune[,indVarCols],
                    y=df.tune[,depVarCol],
                    ntreeTry = 300, stepFactor = 1.5, mtryStart = max(newTry),
                    strata = df.full[,group$colNm], sampsize = sampSizeVec, replace = TRUE)
        
        mtry <- max(y[y[,2] == min(y[,2]),1])
        rm(x,y, df.tune, newTry)
        mtry
      }
    )
    rm(rowCounts)
    
    ###
    # Remove the least important env vars ----
    ##
    
    ntrees <- 1000
    numCores <- 16 #PJM changed from 5       
    
    # do all randomForest calls in parallel, starting with this one
    cl <- makeCluster(numCores)   
    registerDoParallel(cl)
    
    # pass libPath to workers
    x <- clusterCall(cl, function(x) .libPaths(x), .libPaths())
    
    treeSubs <- ceiling(ntrees/numCores)
    
    rf.find.envars <- foreach(ntree = rep(treeSubs,numCores), .combine = randomForest::combine, 
                              .packages = 'randomForest', .multicombine = TRUE) %dopar% {
                                randomForest(df.full[,indVarCols],
                                             y=df.full[,depVarCol],
                                             importance=TRUE,
                                             ntree=ntree,
                                             mtry=mtry,
                                             strata = df.full[,group$colNm],
                                             sampsize = sampSizeVec, replace = TRUE,
                                             norm.votes = TRUE)
                              }
    
    
    impvals <- importance(rf.find.envars, type = 1)
    OriginalNumberOfEnvars <- length(impvals)
    
    # # first remove the bottom of the correlated vars
    # corrdEVs <- corrdEVs[tolower(corrdEVs$gridName) %in% row.names(impvals),]
    # if(nrow(corrdEVs) > 0 ){
    #   for(grp in unique(corrdEVs$correlatedVarGroupings)){
    #     vars <- tolower(corrdEVs[corrdEVs$correlatedVarGroupings == grp,"gridName"])
    #     imp.sub <- impvals[rownames(impvals) %in% vars,, drop = FALSE]
    #     suppressWarnings(varsToDrop <- imp.sub[!imp.sub == max(imp.sub),, drop = FALSE])
    #     impvals <- impvals[!rownames(impvals) %in% rownames(varsToDrop),,drop = FALSE]
    #   }
    #   rm(vars, imp.sub, varsToDrop)
    # }
    
    # set the percentile, here choosing above 35% percentile
    envarPctile <- 0.50
    y <- quantile(impvals, probs = envarPctile)
    impEnvVars <- impvals[impvals > y,]
    subsetNumberofEnvars <- length(impEnvVars)
    rm(y)
    # which columns are these, then flip the non-envars to TRUE
    impEnvVarCols <- names(df.full) %in% names(impEnvVars)
    impEnvVarCols[1:5] <- TRUE
    # subset!
    rf.df.full <- df.full[,impEnvVarCols]
    # reset the indvarcols object
    indVarCols <- c(6:length(names(rf.df.full)))
    
    rm(impvals, impEnvVars, impEnvVarCols)
    
    ##
    # code above is for removing least important env vars
    ##
    
    # prep for validation loop ----
    #now that entire set is cleaned up, split back out to use any of the three DFs below
    df.in2 <- subset(rf.df.full,pres == "1")
    df.abs2 <- subset(rf.df.full, pres == "0")
    df.in2$stratum <- factor(df.in2$stratum)
    df.abs2$stratum <- factor(df.abs2$stratum)
    df.in2$group_id <- factor(df.in2$group_id)
    df.abs2$group_id <- factor(df.abs2$group_id)
    df.in2$pres <- factor(df.in2$pres)
    df.abs2$pres <- factor(df.abs2$pres)
    
    #reset the row names, needed for random subsetting method of df.abs2, below
    row.names(df.in2) <- 1:nrow(df.in2)
    row.names(df.abs2) <- 1:nrow(df.abs2)
    
    #reduce the number of trees if group$vals has more than 30 entries #commented out to be parallel with aquatic
    #this is for validation
    # if(length(group$vals) > 30) {
    # 	ntrees <- 500
    # } else {
    # 	ntrees <- 750
    # }
    ntrees <- 1000
    
    # reduce the number of validation loops if more than 50. 50 is plenty!
    # randomly draw to get the validation set.
    if(length(group$vals) > 50) {
      group$vals <- sample(group$vals, size = 50)
      group$vals <- factor(group$vals)
    } 
    
    ##initialize the Results vectors for output from the jackknife runs
    trRes <- vector("list",length(group$vals))
    names(trRes) <- group$vals[]
    evSet <- vector("list",length(group$vals))
    names(evSet) <- group$vals[]	   
    evRes <- vector("list",length(group$vals))
    names(evRes) <- group$vals[]
    t.f <- vector("list",length(group$vals))
    names(t.f) <- group$vals[]
    t.ctoff <- vector("list",length(group$vals))
    names(t.ctoff) <- group$vals[]
    v.rocr.rocplot <- vector("list",length(group$vals))
    names(v.rocr.rocplot) <- group$vals[]
    v.rocr.auc <- vector("list",length(group$vals))
    names(v.rocr.auc) <- group$vals[]
    v.y <- vector("list",length(group$vals))
    names(v.y) <- group$vals[]
    v.kappa <- vector("list",length(group$vals))
    names(v.kappa) <- group$vals[]
    v.tss <- vector("list",length(group$vals))
    names(v.tss) <- group$vals[]
    v.OvAc <- vector("list",length(group$vals))
    names(v.OvAc) <- group$vals[]
    t.importance <- vector("list",length(group$vals))
    names(t.importance) <- group$vals[]
    t.rocr.pred <- vector("list",length(group$vals))
    names(t.rocr.pred) <- group$vals[]
    v.rocr.pred <- vector("list",length(group$vals))
    names(v.rocr.pred) <- group$vals[]
    
    #######
    ## This is the validation loop. ----
    ## it creates a model for all-but-one group (EO, polygon, or group),
    ## tests if it can predict that group left out,
    ## then moves on to another group, cycling though all groups
    ## Validation stats in tabular form are the final product.
    #######
    
    # calculate the number of trees to send to each core
    treeSubs <- ceiling(ntrees/numCores)
    
    if(length(group$vals)>1){
      for(i in 1:length(group$vals)){
        # Create an object that stores the select command, to be used by subset.
        trSelStr <- parse(text=paste(group$colNm[1]," != '", group$vals[[i]],"'",sep=""))
        evSelStr <- parse(text=paste(group$colNm[1]," == '", group$vals[[i]],"'",sep=""))
        # apply the subset. do.call is needed so selStr can be evaluated correctly
        trSet <- do.call("subset",list(df.in2, trSelStr))
        evSet[[i]] <- do.call("subset",list(df.in2, evSelStr))
        # use sample to grab a random subset from the background points
        BGsampSz <- nrow(evSet[[i]])
        #### TEMPORARY DURING TESTING ####
        if(BGsampSz > nrow(df.abs2)) BGsampSz <- nrow(df.abs2)/2
        evSetBG <- df.abs2[sample(nrow(df.abs2), BGsampSz , replace = FALSE, prob = NULL),]
        # get the other portion for the training set
        TrBGsamps <- attr(evSetBG, "row.names") #get row.names as integers
        trSetBG <-  df.abs2[-TrBGsamps,]  #get everything that isn't in TrBGsamps
        # join em, clean up
        trSet <- rbind(trSet, trSetBG)
        trSet[,group$colNm] <- factor(trSet[,group$colNm])
        evSet[[i]] <- rbind(evSet[[i]], evSetBG)
        
        ssVec <- sampSizeVec[!names(sampSizeVec) == group$vals[[i]]]
        # re-calc pseudo-absence samples to match input training samples if you can
        tot_pseua <- ifelse(sum(ssVec[!names(ssVec) %in% "pseu-a"]) > nrow(trSetBG),
                            nrow(trSetBG), 
                            sum(ssVec[!names(ssVec) %in% "pseu-a"]))
        ssVec["pseu-a"] <- tot_pseua
        
        rm(trSelStr, evSelStr, trSetBG, evSetBG, TrBGsamps, BGsampSz )
        
        trRes[[i]] <- foreach(ntree = rep(treeSubs,numCores), .combine = randomForest::combine, 
                              .packages = 'randomForest', .multicombine = TRUE) %dopar%
          randomForest(trSet[,indVarCols],y=trSet[,depVarCol],
                       importance=TRUE,mtry=mtry,ntree = ntree,
                       strata = trSet[,group$colNm], sampsize = ssVec, replace = TRUE
          )
        
        # run a randomForest predict on the validation data
        evRes[[i]] <- predict(trRes[[i]], evSet[[i]], type="prob")
        # use ROCR to structure the data. Get pres col of evRes (= named "1")
        v.rocr.pred[[i]] <- prediction(evRes[[i]][,"1"],evSet[[i]]$pres)
        # extract the auc for metadata reporting
        v.rocr.auc[[i]] <- performance(v.rocr.pred[[i]], "auc")@y.values[[1]]
        cat("finished run", i, "of", length(group$vals), "\n")
      }
      
      # restructure validation predictions so ROCR will average the figure
      v.rocr.pred.restruct <- v.rocr.pred[[1]]
      #send in the rest
      for(i in 2:length(v.rocr.pred)){
        v.rocr.pred.restruct@predictions[[i]] <- v.rocr.pred[[i]]@predictions[[1]]
        v.rocr.pred.restruct@labels[[i]] <- v.rocr.pred[[i]]@labels[[1]]
        v.rocr.pred.restruct@cutoffs[[i]] <- v.rocr.pred[[i]]@cutoffs[[1]]
        v.rocr.pred.restruct@fp[[i]] <- v.rocr.pred[[i]]@fp[[1]]
        v.rocr.pred.restruct@tp[[i]] <- v.rocr.pred[[i]]@tp[[1]]
        v.rocr.pred.restruct@tn[[i]] <- v.rocr.pred[[i]]@tn[[1]]
        v.rocr.pred.restruct@fn[[i]] <- v.rocr.pred[[i]]@fn[[1]]
        v.rocr.pred.restruct@n.pos[[i]] <- v.rocr.pred[[i]]@n.pos[[1]]
        v.rocr.pred.restruct@n.neg[[i]] <- v.rocr.pred[[i]]@n.neg[[1]]
        v.rocr.pred.restruct@n.pos.pred[[i]] <- v.rocr.pred[[i]]@n.pos.pred[[1]]
        v.rocr.pred.restruct@n.neg.pred[[i]] <- v.rocr.pred[[i]]@n.neg.pred[[1]]
      }
      
      # run a ROC performance with ROCR
      v.rocr.rocplot.restruct <- performance(v.rocr.pred.restruct, "tpr","fpr")
      # send it to perf for the averaging lines that follow
      rf.perf <- v.rocr.rocplot.restruct
      rm(v.rocr.rocplot.restruct)
      ## for infinite cutoff, assign maximal finite cutoff + mean difference
      ## between adjacent cutoff pairs  (this code is from ROCR)
      if (length(rf.perf@alpha.values)!=0) rf.perf@alpha.values <-
        lapply(rf.perf@alpha.values,
               function(x) { isfin <- is.finite(x);
               x[is.infinite(x)] <-
                 (max(x[isfin]) +
                    mean(abs(x[isfin][-1] -
                               x[isfin][-length(x[isfin])])));
               x[is.nan(x)] <- 0.001; #added by tgh to handle vectors length 2
               x})
      
      for (i in 1:length(rf.perf@x.values)) {
        ind.bool <- (is.finite(rf.perf@x.values[[i]]) & is.finite(rf.perf@y.values[[i]]))
        if (length(rf.perf@alpha.values) > 0)
          rf.perf@alpha.values[[i]] <- rf.perf@alpha.values[[i]][ind.bool]
        rf.perf@x.values[[i]] <- rf.perf@x.values[[i]][ind.bool]
        rf.perf@y.values[[i]] <- rf.perf@y.values[[i]][ind.bool]
      }
      rf.perf.sampled <- rf.perf
      
      # create a list of cutoffs to interpolate off of
      alpha.values <- rev(seq(min(unlist(rf.perf@alpha.values)),
                              max(unlist(rf.perf@alpha.values)),
                              length=max(sapply(rf.perf@alpha.values, length))))
      # interpolate by cutoff, values for y and x
      for (i in 1:length(rf.perf.sampled@y.values)) {
        rf.perf.sampled@x.values[[i]] <-
          approxfun(rf.perf@alpha.values[[i]],rf.perf@x.values[[i]],
                    rule=2, ties=mean)(alpha.values)
        rf.perf.sampled@y.values[[i]] <-
          approxfun(rf.perf@alpha.values[[i]], rf.perf@y.values[[i]],
                    rule=2, ties=mean)(alpha.values)
      }
      
      ## compute average curve
      rf.perf.avg <- rf.perf.sampled
      rf.perf.avg@x.values <- list(rowMeans( data.frame( rf.perf.avg@x.values)))
      rf.perf.avg@y.values <- list(rowMeans( data.frame( rf.perf.avg@y.values)))
      rf.perf.avg@alpha.values <- list( alpha.values )
      
      for(i in 1:length(group$vals)){
        ### get threshold
        # get MTP: minimum training presence (minimum votes recieved [probability]
        # for any training point)
        allVotesPrespts <- trRes[[i]]$votes[,"1"][trRes[[i]]$y == 1]
        MTP <- min(allVotesPrespts/numCores)
        # calculations fail if MTP = 0 so if it does, fall back to maxSSS
        if(MTP == 0) {
          # max sensitivity plus specificity (maxSSS per Liu et al 2016)
          # create the prediction object for ROCR. Get pres col from y, prediction from votes (=named "1")
          pred <- prediction(trRes[[i]]$votes[,"1"],trRes[[i]]$y)
          sens <- performance(pred,"sens")
          spec <- performance(pred,"spec")
          sss <- data.frame(cutSens = unlist(sens@x.values),sens = unlist(sens@y.values),
                            cutSpec = unlist(spec@x.values), spec = unlist(spec@y.values))
          sss$sss <- with(sss, sens + spec)
          maxSSS <- sss[which.max(sss$sss),"cutSens"]/numCores
          cutval.rf <- c(1-maxSSS, maxSSS)
          names(cutval.rf) <- c("0","1")
        } else {
          cutval.rf <- c(1-MTP, MTP)
          names(cutval.rf) <- c("0","1")
        }
        
        #apply the cutoff to the validation data
        v.rf.pred.cut <- predict(trRes[[i]], evSet[[i]],type="response", cutoff=cutval.rf)
        #make the confusion matrix
        v.y[[i]] <- table(observed = evSet[[i]][,"pres"],
                          predicted = v.rf.pred.cut)
        #add estimated accuracy measures
        v.y[[i]] <- cbind(v.y[[i]],
                          "accuracy" = c(v.y[[i]][1,1]/sum(v.y[[i]][1,]), v.y[[i]][2,2]/sum(v.y[[i]][2,])))
        #add row, col names
        rownames(v.y[[i]])[rownames(v.y[[i]]) == "0"] <- "background/abs"
        rownames(v.y[[i]])[rownames(v.y[[i]]) == "1"] <- "known pres"
        colnames(v.y[[i]])[colnames(v.y[[i]]) == "0"] <- "pred. abs"
        colnames(v.y[[i]])[colnames(v.y[[i]]) == "1"] <- "pred. pres"
        print(v.y[[i]])
        #Generate kappa statistics for the confusion matrices
        v.kappa[[i]] <- Kappa(v.y[[i]][1:2,1:2])
        #True Skill Statistic
        v.tss[[i]] <- v.y[[i]][2,3] + v.y[[i]][1,3] - 1
        #Overall Accuracy
        v.OvAc[[i]] <- (v.y[[i]][[1,1]]+v.y[[i]][[2,2]])/sum(v.y[[i]][,1:2])
        ### importance measures ###
        #count the number of variables
        n.var <- nrow(trRes[[i]]$importance)
        #get the importance measures (don't get GINI coeff - see Strobl et al. 2006)
        imp <- importance(trRes[[i]], class = NULL, scale = TRUE, type = NULL)
        imp <- imp[,"MeanDecreaseAccuracy"]
        #get number of variables used in each forest
        used <- varUsed(trRes[[i]])
        names(used) <- names(imp)
        t.importance[[i]] <- data.frame("meanDecreaseAcc" = imp,
                                        "timesUsed" = used )
      } #close loop
      
      #housecleaning
      rm(trSet, evSet)
      
      #average relevant validation/summary stats
      # Kappa - wieghted, then unweighted
      K.w <- unlist(v.kappa, recursive=TRUE)[grep("Weighted.value",
                                                  names(unlist(v.kappa, recursive=TRUE)))]
      Kappa.w.summ <- data.frame("mean"=mean(K.w), "sd"=sd(K.w),"sem"= sd(K.w)/sqrt(length(K.w)))
      K.unw <- unlist(v.kappa, recursive=TRUE)[grep("Unweighted.value",
                                                    names(unlist(v.kappa, recursive=TRUE)))]
      Kappa.unw.summ <- data.frame("mean"=mean(K.unw), "sd"=sd(K.unw),"sem"= sd(K.unw)/sqrt(length(K.unw)))
      #AUC - area under the curve
      auc <- unlist(v.rocr.auc)
      auc.summ <- data.frame("mean"=mean(auc), "sd"=sd(auc),"sem"= sd(auc)/sqrt(length(auc)))
      #TSS - True skill statistic
      tss <- unlist(v.tss) 
      tss.summ <- data.frame("mean"=mean(tss), "sd"=sd(tss),"sem"= sd(tss)/sqrt(length(tss)))
      #Overall Accuracy
      OvAc <- unlist(v.OvAc)
      OvAc.summ <- data.frame("mean"=mean(OvAc), "sd"=sd(OvAc),"sem"= sd(OvAc)/sqrt(length(OvAc)))
      #Specificity and Sensitivity
      v.y.flat <- abind(v.y,along=1)  #collapsed confusion matrices
      v.y.flat.sp <- v.y.flat[rownames(v.y.flat)=="background/abs",]
      v.y.flat.sp <- as.data.frame(v.y.flat.sp, row.names = 1:length(v.y.flat.sp[,1]))
      specif <- v.y.flat.sp[,"pred. abs"]/(v.y.flat.sp[,"pred. abs"] + v.y.flat.sp[,"pred. pres"])   #specificity
      specif.summ <- data.frame("mean"=mean(specif), "sd"=sd(specif),"sem"= sd(specif)/sqrt(length(specif)))
      v.y.flat.sn <- v.y.flat[rownames(v.y.flat)=="known pres",]
      v.y.flat.sn <- as.data.frame(v.y.flat.sn, row.names = 1:length(v.y.flat.sn[,1]))
      sensit <- v.y.flat.sn[,"pred. pres"]/(v.y.flat.sn[,"pred. pres"] + v.y.flat.sn[,"pred. abs"])    #sensitivity
      sensit.summ <- data.frame("mean"=mean(sensit), "sd"=sd(sensit),"sem"= sd(sensit)/sqrt(length(sensit)))
      
      rf_model_summary_table <- data.frame(metric=c("Weighted Kappa", "Unweighted Kappa", "AUC",
                                                    "TSS", "Overall Accuracy", "Specificity",
                                                    "Sensitivity"),
                                           metric_mn=c(Kappa.w.summ$mean, Kappa.unw.summ$mean,auc.summ$mean,
                                                       tss.summ$mean, OvAc.summ$mean, specif.summ$mean,
                                                       sensit.summ$mean),
                                           metric_sd=c(Kappa.w.summ$sd, Kappa.unw.summ$sd,auc.summ$sd,
                                                       tss.summ$sd, OvAc.summ$sd, specif.summ$sd,
                                                       sensit.summ$sd),
                                           metric_sem=c(Kappa.w.summ$sem, Kappa.unw.summ$sem,auc.summ$sem,
                                                        tss.summ$sem, OvAc.summ$sem, specif.summ$sem,
                                                        sensit.summ$sem))
      rf_model_summary_table
    } else {
      cat("Only one polygon, can't do validation", "\n")
      cutval <- NA
    }
    
    # increase the number of trees for the full model
    ntrees <- 2000
    treeSubs <- ceiling(ntrees/numCores)
    
    ####
    #   run the full model ----
    ####
    cat("... creating full model \n")
    
    rf.full <- foreach(ntree = rep(treeSubs,numCores), .combine = randomForest::combine, 
                       .packages = 'randomForest', .multicombine = TRUE) %dopar% {
                         randomForest(rf.df.full[,indVarCols],
                                      y=rf.df.full[,depVarCol],
                                      importance=TRUE,
                                      ntree=ntree,
                                      mtry=mtry,
                                      strata = rf.df.full[,group$colNm],
                                      sampsize = sampSizeVec, replace = TRUE,
                                      norm.votes = TRUE)
                       }
    
    # # write out input data
    # # connect to DB ..
    # db <- dbConnect(SQLite(),dbname=nm_db_file)
    # 
    # # write model input data to database before any other changes made
    rf_model_inputs_table <- data.frame(model_run_name = rf_model_version,
                                        algorithm = algo,
                                        datetime = as.character(Sys.time()),
                                        feat_count = length(unique(df.in2$stratum)),
                                        feat_grp_count = length(unique(df.in2$group_id)),
                                        jckn_grp_column = group$colNm,
                                        jckn_grp_type = group$JackknType,
                                        mn_grp_subsamp = mean(sampSizeVec[!names(sampSizeVec) == "pseu-a"]),
                                        min_grp_subsamp = min(sampSizeVec[!names(sampSizeVec) == "pseu-a"]),
                                        max_grp_subsamp = max(sampSizeVec[!names(sampSizeVec) == "pseu-a"]),
                                        tot_obs_subsamp = sum(sampSizeVec[!names(sampSizeVec) == "pseu-a"]),
                                        tot_bkgd_subsamp = sum(sampSizeVec[names(sampSizeVec) == "pseu-a"]),
                                        obs_count = nrow(df.in2),
                                        bkgd_count = nrow(df.abs2)
    )
    # dbExecute(db, paste0("DELETE FROM tblModelInputs where table_code = '", baseName, 
    #                      "' and algorithm = '", algo, "';")) # remove any previously prepped dataset entry
    # dbWriteTable(db, "tblModelInputs", tblModelInputs, append = TRUE)
    # 
    # # write validation data
    rf_model_summary_table <- cbind("model_run_name" = rep(rf_model_version, nrow(rf_model_summary_table)),
                                    "model_creation_date" = rf_model_creation_date,
                                    "algorithm" = rep(algo, nrow(rf_model_summary_table)),
                                    rf_model_summary_table)
    # 
    # dbExecute(db, paste0("DELETE FROM tblModelResultsValidationStats where model_run_name = '", model_run_name, 
    #                      "' and algorithm = '", algo, "';")) # remove any previously prepped dataset entry
    # dbWriteTable(db, "tblModelResultsValidationStats", summ.table, append = TRUE)
    # 
    # 
    # dbDisconnect(db)
    # rm(db)
    
    ####
    # Importance measures ----
    ####
    #get the importance measures (don't get GINI coeff - see Strobl et al. 2006)
    f.imp <- importance(rf.full, class = NULL, scale = TRUE, type = NULL)
    f.imp <- f.imp[,"MeanDecreaseAccuracy"]
    
    # db <- dbConnect(SQLite(),dbname=nm_db_file)  
    # # get importance data, set up a data frame
    EnvVars <- data.frame(gridName = names(f.imp), impVal = f.imp, fullName="", stringsAsFactors = FALSE)
    # #set the query for the following lookup, note it builds many queries, equal to the number of vars
    # SQLquery <- paste("SELECT gridName, fullName FROM lkpEnvVars WHERE gridName COLLATE NOCASE in ('", paste(EnvVars$gridName,sep=", "),
    # 					"'); ", sep="")
    # #cycle through all select statements, put the results in the df
    # for(i in 1:length(EnvVars$gridName)) {
    #   try(EnvVars$fullName[i] <- as.character(dbGetQuery(db, statement = SQLquery[i])[,2]))
    # }
    
    ###
    # partial plot data ----
    ###
    #get the order for the importance charts
    ord <- order(EnvVars$impVal, decreasing = TRUE)[1:length(indVarCols)]
    if(length(ord) > 15){
      pPlotListLen <- 15
    } else {
      pPlotListLen <- length(ord)
    }
    
    cat("... calculating partial plots \n")
    
    ### subsample, grouped by pres/abs, to speed up partial plots
    ppPres <- rf.df.full[rf.df.full$pres == 1, ]
    ppAbs <- rf.df.full[rf.df.full$pres == 0, ]
    ppPresSamp <- min(c(nrow(ppPres), 6000)) # take all pres samples, or 6000, whichever is less
    ppPresSamp <- sample(1:nrow(ppPres), size = round(ppPresSamp), replace = FALSE)
    ppPresSamp <- ppPres[ppPresSamp,]
    ppAbsSamp <- min(c(nrow(ppAbs), 6000)) # take all abs samples, or 6000, whichever is less
    ppAbsSamp <- sample(1:nrow(ppAbs), size = round(ppAbsSamp), replace = FALSE)
    ppAbsSamp <- ppAbs[ppAbsSamp,]
    
    ppPreddata <- rbind(ppPresSamp, ppAbsSamp)
    
    # run partial plots in parallel
    # curvars = names(f.imp[ord])[1:pPlotListLen]
    curvars = names(f.imp)[1:pPlotListLen]
    rf.pPlots <- foreach(i = iter(curvars), .packages = 'randomForest') %dopar% {
      do.call("partialPlot", list(x = rf.full, pred.data = ppPreddata[,indVarCols],
                                  x.var = i,
                                  which.class = 1,
                                  plot = FALSE))
    }
    
    #fill in names
    names(rf.pPlots) <- c(1:pPlotListLen)
    for(i in 1:length(rf.pPlots)){
      rf.pPlots[[i]]$gridName <- curvars[[i]]
      rf.pPlots[[i]]$fname <- EnvVars$fullName[ord[i]]
    }
    rm(ppPres, ppAbs, ppPresSamp, ppAbsSamp, ppPreddata)
    
    # clear out registered clusters as caret is unhappy with them ...
    unregister <- function() {
      env <- foreach:::.foreachGlobals
      rm(list=ls(name=env), pos=env)
    }
    unregister()
    
    stopCluster(cl)
    rm(cl)
    
    # # tblModelResultsVarsUsed
    rf_var_importance <- data.frame(model_run_name = rf_model_version,
                                    algorithm = algo,
                                    inFinalModel = 0,
                                    EnvVars[c("gridName","impVal")]
    )
    rf_var_importance$inFinalModel[!is.na(rf_var_importance$impVal)] <- 1
    # dbWriteTable(db, "tblModelResultsVarsUsed", varImpDB, append = TRUE)
    # dbDisconnect(db)
    # 
    # closeAllConnections()
    # 
    # message(paste0("Saved rdata file: '", model_run_name , "'."))
    
    ### this isn't needed but keeping for a bit in one script in case we decide otherwise
    #serious cleanup
    # if(isNamespaceLoaded("doParallel")) unloadNamespace("doParallel")
    # if(isNamespaceLoaded("foreach")) unloadNamespace("foreach")
    # if(isNamespaceLoaded("iterators")) unloadNamespace("iterators")
    # if(isNamespaceLoaded("parallel")) unloadNamespace("parallel")
    # if(isNamespaceLoaded("RSQLite")) unloadNamespace("RSQLite")
    # if(isNamespaceLoaded("ROCR")) unloadNamespace("ROCR")
    # if(isNamespaceLoaded("vcd")) unloadNamespace("vcd")
    # if(isNamespaceLoaded("abind")) unloadNamespace("abind")
    # if(isNamespaceLoaded("foreign")) unloadNamespace("foreign")
    # if(isNamespaceLoaded("randomForest")) unloadNamespace("randomForest")
    
    # Package output
    model_output[["rf"]] <- list(
      df = rf.df.full,
      full_model = rf.full,
      model_inputs_table = rf_model_inputs_table,
      model_summary_table = rf_model_summary_table,
      var_importance = rf_var_importance,
      partial_plots = rf.pPlots
    )
    
  }
  
  if ("xgb" %in% algorithms){
    
    algo <- "xgb"
    xgb_model_version <- paste(presence_pts$species_cd[1], algo, gsub("-", "", Sys.Date() %>% as.character()), sep = "_")
    xgb_model_creation_date <- Sys.Date() %>% as.character()
    
    # File: 3_createModel.R
    # Purpose: to create the random forest model. This includes:
    # - validate using leave-one-out jackknifing
    # - create initial model to remove poorest performing env vars
    # - create a final model using all presence points, stratify by EO using RA
    # - build partial plots of top performing env vars for metadata output
    
    # having parallel package loaded messes up caret
    #detach("package:doParallel", unload=TRUE)
    
    
    ###
    # Remove the least important env vars ----
    ##
    
    
    df.full.xgb <- xgb.DMatrix(as.matrix(df.full[,indVarCols]), 
                               label=as.integer(as.character(df.full$pres)))
    
    param <- list(
      booster = "gbtree",
      max_depth = 5, 
      eta = 0.3, 
      verbose = 0, 
      nthread = 4,
      objective = "binary:logistic", 
      eval_metric = "auc")
    
    xgb.find.envars <- xgb.train(params=param, data = df.full.xgb,
                                 nrounds = 20)
    
    xgb.impvals <- xgb.importance(model=xgb.find.envars)
    
    #xgb already removed a bunch of env vars?
    OriginalNumberOfEnvars <- nrow(xgb.impvals)
    
    # first remove the bottom of the correlated vars
    # corrdEVs <- corrdEVs[tolower(corrdEVs$gridName) %in% xgb.impvals$Feature,]
    # if(nrow(corrdEVs) > 0 ){
    #   for(grp in unique(corrdEVs$correlatedVarGroupings)){
    #     vars <- tolower(corrdEVs[corrdEVs$correlatedVarGroupings == grp,"gridName"])
    #     imp.sub <- xgb.impvals[xgb.impvals$Feature %in% vars,, drop = FALSE]
    #     suppressWarnings(varsToDrop <- imp.sub[!imp.sub$Gain == max(imp.sub$Gain),, drop = FALSE])
    #     xgb.impvals <- xgb.impvals[!xgb.impvals$Feature %in% varsToDrop$Feature,,drop = FALSE]
    #   }
    #   rm(vars, imp.sub, varsToDrop)
    # }
    
    # set the percentile, here choosing above 25% percentile
    envarPctile <- 0.25
    y <- quantile(xgb.impvals$Gain, probs = envarPctile)
    impEnvVars <- xgb.impvals[xgb.impvals$Gain > y,]
    subsetNumberofEnvars <- nrow(impEnvVars)
    rm(y)
    # which columns are these, then flip the non-envars to TRUE
    impEnvVarCols <- names(df.full) %in% impEnvVars$Feature
    impEnvVarCols[1:5] <- TRUE
    # subset!
    xgb.df.full <- df.full[,impEnvVarCols]
    # reset the indvarcols object
    indVarCols <- c(6:length(names(xgb.df.full)))
    
    rm(xgb.impvals, impEnvVars, impEnvVarCols)
    
    
    
    # outPth <- file.path(loc_model, ElementNames$Code,"outputs","ensemble")
    # dir.create(outPth, showWarnings = FALSE)
    
    if(length(group$vals) > 50){
      kf <- 50
    } else {
      kf <- length(group$vals)
    }
    
    # randomly assign group ids to background points so background pts
    # can be drawn into KFolds
    
    # vector to split up represents row indices
    rowInd <- seq(1, nrow(xgb.df.full[xgb.df.full$pres == "0",]),by=1)
    # first, randomly assign based on number of pres records in each group
    presCounts <- table(xgb.df.full[,group$colNm])
    presCounts <- presCounts[!names(presCounts) == "pseu-a"]
    # start by assigning a base number to avoid huge skews
    drawAmt <- ifelse(presCounts>5,5,presCounts)
    # now get weights for distributing the rest
    numToDistribute <- length(rowInd) - sum(drawAmt)
    weights <- presCounts/sum(presCounts)
    drawAmt <- drawAmt + round(numToDistribute*weights)
    # rounding might keep this from coming out even, so fix if it didn't
    if(sum(drawAmt) > length(rowInd)){
      reduceAmt <- sum(drawAmt) - length(rowInd)
      # take one away from each item till it balances
      drawAmt[1:reduceAmt] <- drawAmt[1:reduceAmt] - 1
    }
    if(sum(drawAmt) < length(rowInd)){
      increaseAmt <- length(rowInd) - sum(drawAmt)
      if(increaseAmt < length(drawAmt)){
        # add one to each item as needed
        drawAmt[1:increaseAmt] <- drawAmt[1:increaseAmt] + 1
      } else {
        # find the biggest set and add
        drawAmt[drawAmt == max(drawAmt)] <- max(drawAmt) + increaseAmt
      }
    }
    # now randomly place absence rows into groups for k-folding
    sampVec <- integer(0)
    for(i in 1:length(drawAmt)){
      x <- sample(rowInd, drawAmt[i])
      names(x) <- rep(names(drawAmt[i]),drawAmt[i])
      rowInd <- rowInd[!rowInd %in% x]
      sampVec <- c(sampVec, x)
    }
    rm(rowInd, presCounts, drawAmt, numToDistribute, weights, i, x)
    
    #sampVec2 <- unlist(y, use.names = TRUE)
    #sampVec <- c(sampVec, sampVec2)
    sampVec <- sampVec[order(sampVec)]
    # this defines group assignment for each df.full row (assumed pres=1 is first)
    fullSampVec <- c(as.character(xgb.df.full[xgb.df.full$pres == "1",group$colNm]), names(sampVec))
    xgb.df.full$stratum <- fullSampVec
    
    # subsample down to a reasonable number 
    
    # subsample function. x is data frame. ssvec is sampSizeVec (not including pseu-a)
    # this will get the df down to the size indicated by sampSizeVec, 
    # and possibly could be used in a recipe inside the caret model (as a downsample)
    # but I can't figure that out yet. 
    # use it here to downsample before modeling
    # subSampByGp <- function(x, ssvec) {
    #   df.sub = x[FALSE,]
    #   for(i in 1:length(ssvec)){
    #     cls <- names(ssvec[i])
    #     toDrawFrom <- x[x$stratum == cls, ]
    #     draws <- sample(1:nrow(toDrawFrom), ssvec[[i]])
    #     drawnSamps <- toDrawFrom[draws,]
    #     df.sub <- rbind(df.sub, drawnSamps)
    #   }
    #   df.sub
    # }
    
    subSampByGp <- function(x, ssvec, gpColName) {
      df.sub = x[FALSE,]
      for(i in 1:length(ssvec)){
        cls <- names(ssvec[i])
        if(cls == "pseu-a"){
          toDrawFrom <- x[x[,"pres"] == 0, ]
        } else {
          toDrawFrom <- x[x[,gpColName] == cls & x[,"pres"] == 1, ]
        }
        if(ssvec[[i]] > nrow(toDrawFrom)) { #lame
          ssvec[[i]] <- nrow(toDrawFrom)
        }
        draws <- sample(1:nrow(toDrawFrom), ssvec[[i]], replace = FALSE)  
        drawnSamps <- toDrawFrom[draws,]
        df.sub <- rbind(df.sub, drawnSamps)
      }
      df.sub
    }
    
    #df.full.s <- subSampByGp(df.full, sampSizeVec[-grep("pseu-a", names(sampSizeVec))], group$colNm)
    xgb.df.full.s <- subSampByGp(xgb.df.full, sampSizeVec, group$colNm)
    ## TODO sample background down to a reasonable number ##### research it?
    #df.full.s <- rbind(df.full.s, df.full[df.full$pres == 0,])
    
    # define the folds
    folds <- groupKFold(xgb.df.full.s$stratum, k = kf) 
    
    #library(pROC)
    # xgbfitControl <- trainControl(
    #   method = "LGOCV",  #leave group out cross validation
    #   index = folds,
    #   number = 1,
    #   summaryFunction = twoClassSummary,
    #   classProbs = TRUE
    #   )
    
    
    
    xgbfitControl <- trainControl(
      method = "repeatedcv",
      index = folds,
      number = 1,
      summaryFunction = twoClassSummary,
      classProbs = TRUE,
      savePredictions = TRUE
    )
    
    # caret seems to need this
    xgb.df.full.s$pres <- as.character(xgb.df.full.s$pres)
    xgb.df.full.s[xgb.df.full.s$pres == "0","pres"] <- "abs"
    xgb.df.full.s[xgb.df.full.s$pres == "1","pres"] <- "pres"
    xgb.df.full.s$pres <- as.factor(xgb.df.full.s$pres)
    
    # xgbGrid <-  expand.grid(
    #   nrounds = c(5, 10, 20),
    #   max_depth = c(2, 5, 9), 
    #   eta = 1,
    #   gamma = 0,
    #   colsample_bytree = 1,
    #   min_child_weight = 1,
    #   subsample = 1
    # )
    #modRecipe <- recipe(pres ~ ., data = df.full.s
    #)
    
    # run validation with caret
    xgbFit1 <- train(pres ~ ., data = xgb.df.full.s[,c(depVarCol,indVarCols)],
                     method = "xgbTree",
                     trControl = xgbfitControl,
                     metric = "ROC")
    
    # flip *back* for raw xgb
    xgb.df.full.s$pres <- as.character(xgb.df.full.s$pres)
    xgb.df.full.s[xgb.df.full.s$pres == "abs","pres"] <- "0"
    xgb.df.full.s[xgb.df.full.s$pres == "pres","pres"] <- "1"
    xgb.df.full.s$pres <- as.factor(xgb.df.full.s$pres)
    
    df.full.s.xgb <- xgb.DMatrix(as.matrix(xgb.df.full[,indVarCols]), 
                                 label=as.integer(as.character(xgb.df.full$pres)))
    
    ####
    #   run the full model ----
    ####
    
    xgb.full <- xgb.train(params=xgbFit1$finalModel$params, 
                          data = df.full.s.xgb,
                          nrounds = 100)
    
    
    # # write out input data
    # # connect to DB ..
    # db <- dbConnect(SQLite(),dbname=nm_db_file)
    # 
    # # write model input data to database before any other changes made
    xgb_model_inputs_table <- data.frame(model_run_name = xgb_model_version,
                                         algorithm = algo,
                                         datetime = as.character(Sys.time()),
                                         feat_count = length(unique(df.in2$stratum)),
                                         feat_grp_count = length(unique(df.in2$group_id)),
                                         jckn_grp_column = group$colNm,
                                         jckn_grp_type = group$JackknType,
                                         mn_grp_subsamp = mean(sampSizeVec[!names(sampSizeVec) == "pseu-a"]),
                                         min_grp_subsamp = min(sampSizeVec[!names(sampSizeVec) == "pseu-a"]),
                                         max_grp_subsamp = max(sampSizeVec[!names(sampSizeVec) == "pseu-a"]),
                                         tot_obs_subsamp = sum(sampSizeVec[!names(sampSizeVec) == "pseu-a"]),
                                         tot_bkgd_subsamp = nrow(xgb.df.full.s[xgb.df.full.s$pres == "0",]),
                                         obs_count = nrow(df.in2),
                                         bkgd_count = nrow(df.abs2)
    )
    # dbExecute(db, paste0("DELETE FROM tblModelInputs where table_code = '", baseName, 
    #                      "' and algorithm = '", algo, "';")) # remove any previously prepped dataset entry
    # dbWriteTable(db, "tblModelInputs", tblModelInputs, append = TRUE)
    # # write validation data
    # 
    xgbTuneOutput <- xgbFit1$results[as.numeric(rownames(xgbFit1$finalModel$tuneValue)),]
    # 
    # 
    xgb_model_summary_table <- as.data.frame(cbind("model_run_name" = rep(xgb_model_version, 4),
                                                   "model_creation_date" = xgb_model_creation_date,
                                                   "algorithm" = rep(algo, 4),
                                                   "metric" = c("AUC","Sensitivity","Specificity","TSS"),
                                                   "metric_mn" = c(xgbTuneOutput$ROC,
                                                                   xgbTuneOutput$Sens,
                                                                   xgbTuneOutput$Spec,
                                                                   xgbTuneOutput$Sens + xgbTuneOutput$Spec - 1),
                                                   "metric_sd" = c(xgbTuneOutput$ROCSD,
                                                                   xgbTuneOutput$SensSD,
                                                                   xgbTuneOutput$SpecSD,
                                                                   NA)
    ))
    # 
    # dbExecute(db, paste0("DELETE FROM tblModelResultsValidationStats where model_run_name = '", model_run_name, 
    #                      "' and algorithm = '", algo, "';")) # remove any previously prepped dataset entry
    # dbWriteTable(db, "tblModelResultsValidationStats", xgb.summ.table, append = TRUE)
    # 
    # dbDisconnect(db)
    # rm(db)
    
    ####
    # Importance measures ----
    ####
    xgb.impvals <- xgb.importance(model=xgb.full)
    g.imp <- xgb.impvals[,c("Feature","Gain")]
    
    #QC: are all features in impvals?
    #if not, add them as zero (assume model dropped them)
    if(FALSE %in% (xgb.full$feature_names %in% g.imp$Feature)){
      toAdd <- xgb.full$feature_names[!xgb.full$feature_names %in% g.imp$Feature]
      g.imp <- rbind(g.imp, data.frame("Feature" = toAdd, "Gain" = rep(0, length(toAdd))))
    }
    
    # db <- dbConnect(SQLite(),dbname=nm_db_file)
    # # get importance data, set up a data frame
    xgb.EnvVars <- data.frame(gridName = g.imp$Feature, impVal = g.imp$Gain, fullName="", stringsAsFactors = FALSE)
    # #set the query for the following lookup, note it builds many queries, equal to the number of vars
    # SQLquery <- paste("SELECT gridName, fullName FROM lkpEnvVars WHERE gridName COLLATE NOCASE in ('", paste(xgb.EnvVars$gridName,sep=", "),
    # 					"'); ", sep="")
    #cycle through all select statements, put the results in the df
    # for(i in 1:length(xgb.EnvVars$gridName)) {
    #   try(xgb.EnvVars$fullName[i] <- as.character(dbGetQuery(db, statement = SQLquery[i])[,2]))
    # }
    ##clean up
    #dbDisconnect(db)
    
    # ###
    # # partial plot data ----
    # ###
    xgb.df.full.s$pres <- as.integer(as.character(xgb.df.full.s$pres))
    
    #get the order for the importance charts
    ord <- order(xgb.EnvVars$impVal, decreasing = TRUE)[1:length(indVarCols)]
    if(length(ord) > 9){
      pPlotListLen <- 9
    } else {
      pPlotListLen <- length(ord)
    }
    
    xgb.pPlots <- xgb.plot.shap(data = as.matrix(xgb.df.full.s[,indVarCols]), 
                                model = xgb.full, 
                                features = as.character(xgb.EnvVars$gridName[1:pPlotListLen]),
                                target_class = 1,
                                plot = FALSE)
    
    xgb.pPlots$fullNames <- xgb.EnvVars$fullName
    
    # write model metadata to db
    # tblModelResultsVarsUsed
    xgb_var_importance <- data.frame(model_run_name = xgb_model_version,
                                     algorithm = algo,
                                     inFinalModel = 0,
                                     xgb.EnvVars[c("gridName","impVal")])
    xgb_var_importance$inFinalModel[!is.na(xgb_var_importance$impVal)] <- 1
    # dbWriteTable(db, "tblModelResultsVarsUsed", varImpDB, append = TRUE)
    # dbDisconnect(db)
    # 
    # closeAllConnections()
    
    # Package output
    model_output[["xgb"]] = list(
      df = xgb.df.full,
      full_model = xgb.full,
      model_outputs_table = xgb_model_inputs_table,
      model_summary_table = xgb_model_summary_table,
      var_importance = xgb_var_importance,
      partial_plots = xgb.pPlots
    )
    
  }
  
  return(model_output)
}

project_models <- function(model_outputs = species_model_objects, 
                           pred_metadata = env_predictor_metadata
){
  
  algorithms <- names(model_outputs)
  model_projection_output <- vector("list", length = length(algorithms))
  names(model_projection_output) <- algorithms
  env_predictors <- NULL
  env_predictor_names <- NULL
  
  if ("rf" %in% algorithms){
    
    if (!is.null(env_predictors)){
      env_predictor_names <- setdiff(names(model_outputs$rf$full_model$forest$xlevels), env_predictor_names)
    } else {
      env_predictor_names <- names(model_outputs$rf$full_model$forest$xlevels)
    }
    
    env_predictors <- c(env_predictors, env_predictor_files[grep(paste0(env_predictor_names, collapse = "|"), env_predictor_files)] %>% 
                          furrr::future_map(function(r){
                            env_raster <- raster::raster(r) %>% 
                              raster::crop(y = as_Spatial(model_extent)) #%>%  #pjm turned off- already masked for full
                              #raster::mask(mask = as_Spatial(model_extent))
                            env_raster
                          }))
    
    rf_env_stack <- env_predictors %>% raster::stack()
    
    rf_raster_map <- raster::predict(object=rf_env_stack, model=model_outputs$rf$full_model, type = "prob", index=2)
    rf_predictions <- data.frame(presabs = as.numeric(as.character(model_outputs$rf$full_model$y)), prob = (model_outputs$rf$full_model$votes/ncol(model_outputs$rf$full_model$forest$ndbigtree))[, 2])
    
    model_projection_output[["rf"]] <- list(
      predictions_map_continuous = rf_raster_map,
      predictions_df = rf_predictions
    )
    
  }
  
  if ("xgb" %in% algorithms){
    
    if (!is.null(env_predictors)){
      env_predictor_names <- setdiff(species_model_objects$xgb$full_model$feature_names, env_predictor_names)
    } else {
      env_predictor_names <- species_model_objects$xgb$full_model$feature_names
    }
    
    # the 'handle' might go awry
    xgb.full <- xgb.Booster.complete(species_model_objects$xgb$full_model)
    
    # deal with how xgb wants input data
    predfun <- function(m, d){
      colus <- length(d)
      y <- matrix(unlist(d), ncol = colus, byrow = FALSE)
      names(y) <- names(d)
      predict(m,y)
    }
    
    env_predictors <- c(env_predictors, env_predictor_files[grep(paste0(env_predictor_names, collapse = "|"), env_predictor_files)] %>% 
                          furrr::future_map(function(r){
                            env_raster <- raster::raster(r) %>% 
                              raster::crop(y = as_Spatial(model_extent)) %>% 
                              raster::mask(mask = as_Spatial(model_extent))
                            env_raster
                          })) 
    
    xgb_env_stack <- env_predictors[unlist(purrr::map(env_predictors, names)) %in% species_model_objects$xgb$full_model$feature_names] %>% 
      raster::stack()
    
    xgb_raster_map <- raster::predict(object = xgb_env_stack, model = xgb.full, fun=predfun, index = 2)
    
    #needs to be re-created ... contains an external pointer 
    xgb.df.full.s <- species_model_objects$xgb$df
    df.full.s.xgb <- xgb.DMatrix(as.matrix(xgb.df.full.s[,6:ncol(xgb.df.full.s)]), 
                                 label=as.integer(as.character(xgb.df.full.s$pres)))
    xgb_predictions <- predict(species_model_objects$xgb$full_model, df.full.s.xgb)
    xgb_predictions <- as.data.frame(cbind(xgb.df.full.s[,c("pres","group_id","stratum")], "prob" = xgb_predictions))
    
    model_projection_output[["xgb"]] <- list(
      predictions_map_continuous = xgb_raster_map,
      predictions_df = xgb_predictions
    )
  }
  
  return(model_projection_output)
  
}

project_models_aquatic <- function(model_outputs = species_model_objects, 
                                   pred_metadata = env_predictor_metadata
){
  
  algorithms <- names(model_outputs)
  model_projection_output <- vector("list", length = length(algorithms))
  names(model_projection_output) <- algorithms
  
  if ("rf" %in% algorithms){
    #' # Project species habitat models across model extent
    #' ## Attribute background reaches
    #' ### Write query
    species_model_prediction_query <- paste0("SELECT * FROM ", "background_reaches","_att WHERE COMID IN ('", paste(model_reaches$COMID, collapse = "','"),"')") 
    #' ### Run query
    species_model_prediction_data <- dbGetQuery(aquatic_predictors_db, species_model_prediction_query) %>% dplyr::filter(complete.cases(.))
    #' ## Generate prediction
    species_model_predictions <- predict(object = species_model_objects$rf$full_model, newdata = species_model_prediction_data, type = "prob")
    #' ## Bind prediction data
    species_model_predictions <- cbind(species_model_prediction_data, prob = species_model_predictions[, which(colnames(species_model_predictions) == "1")])
    species_model_predictions <- model_reaches %>% dplyr::select(COMID) %>% left_join(species_model_predictions %>% dplyr::select(COMID, prob), by = "COMID")
    #' ## Generate predictions data.frame
    predictions_df <- data.frame(presabs = as.numeric(as.character(species_model_objects$rf$full_model$y)), prob = (species_model_objects$rf$full_model$votes/ncol(species_model_objects$rf$full_model$forest$ndbigtree))[, 2])
    
    model_projection_output[["rf"]] <- list(
      predictions_map_continuous = species_model_predictions,
      predictions_df = predictions_df
    )
  }
  
  return(model_projection_output)
  
}

threshold_predictions <- function(species_predictions = species_model_predictions,
                                  domain = c("terrestrial", "aquatic"),
                                  sensitivity_level = 1,
                                  specificity_level = 1,
                                  set_thresholds = NULL){
  
  domain <- match.arg(domain)
  algorithms <- names(species_predictions)
  
  categorical_predictions <- purrr::map(algorithms, function(algo){
    
    predictions_df <- species_predictions[[algo]]$predictions_df
    roc_curve <- pROC::roc(predictions_df$presabs, predictions_df$prob, quiet = T)
    roc_curve$thresholds <- roc_curve$thresholds[-c(1, length(roc_curve$thresholds))]
    a <- pROC::smooth(roc_curve, method='binormal', n=1024)
    
    if (is.null(set_thresholds)){
      thresholds <- c(roc_curve$thresholds[which(roc_curve$sensitivities >= sensitivity_level) %>% tail(1)],
                      rev(roc_curve$thresholds)[findInterval((1-a$specificities)[which.max(a$specificities+a$sensitivities-1)], rev(1-roc_curve$specificities))],
                      roc_curve$thresholds[which(roc_curve$specificities >= specificity_level) %>% head(1)]
      ) %>% 
        sort()
    } else {
      thresholds <- set_thresholds
    }
    
    thresholds <- thresholds %>% set_names(c("low", "medium", "high")) #(c("low", "medium", "high")) PJM changed, dropped high
    
    thresholds.xy <- c(low.x = 1 - roc_curve$specificities[roc_curve$thresholds == thresholds["low"]],
                       low.y = roc_curve$sensitivities[roc_curve$thresholds == thresholds["low"]],
                       medium.x = 1 - roc_curve$specificities[roc_curve$thresholds == thresholds["medium"]],
                       medium.y = roc_curve$sensitivities[roc_curve$thresholds == thresholds["medium"]],
                       high.x = 1 - roc_curve$specificities[roc_curve$thresholds == thresholds["high"]],
                       high.y = roc_curve$sensitivities[roc_curve$thresholds == thresholds["high"]]
    )
    
    if (domain == "terrestrial"){
      predictions_map_categorical <- species_predictions[[algo]]$predictions_map_continuous
      predictions_map_categorical[] <- ifelse(predictions_map_categorical[] >= thresholds["high"], 1,
                                              ifelse(predictions_map_categorical[] >= thresholds["medium"] & predictions_map_categorical[] < thresholds["high"], 2,
                                                     ifelse(
                                                       predictions_map_categorical[] >= thresholds["low"] & predictions_map_categorical[] < thresholds["medium"], 3, 4
                                                     )
                                              )
      ) %>% as.factor() 
    }
    
    if (domain == "aquatic"){
      predictions_map_categorical <- species_predictions[[algo]]$predictions_map_continuous
      predictions_map_categorical <- predictions_map_categorical %>% 
        dplyr::mutate(cat = ifelse(prob >= thresholds["high"], 1,
                                   ifelse(prob >= thresholds["medium"] & prob < thresholds["high"], 2,
                                          ifelse(
                                            prob >= thresholds["low"] & prob < thresholds["medium"], 3, 4
                                          )
                                   )
        ) %>% as.factor()
        ) %>% 
        dplyr::select(-prob)
      
    }
    
    out <- list(
      roc_curve = roc_curve,
      thresholds = thresholds,
      thresholds.xy = thresholds.xy,
      predictions_map_categorical = predictions_map_categorical
    )
    
    out
  }
  ) %>% set_names(algorithms)
  
  return(categorical_predictions) 
  
}

rasterize_aquatic_predictions <- function(species_predictions = species_model_predictions, 
                                          aggregation_coefficient = 2){
  
  algorithms <- names(species_predictions)
  
  model_extent_raster <- raster::raster("S://terrest_envV_2/clim_bio1_AnnMnTmp_x100_int_30m.tif") %>% 
    raster::crop(y = as_Spatial(model_extent)) %>% 
    raster::mask(mask = as_Spatial(model_extent))
  
  model_extent_raster_agg <- raster::aggregate(model_extent_raster, aggregation_coefficient)
  
  raster_predictions <- purrr::map(algorithms, function(algo){
    
    algo_predictions <- species_predictions[[algo]]$predictions_map_continuous
    
    species_raster <- raster::rasterize(algo_predictions %>% dplyr::select(prob) %>% dplyr::mutate(prob = as.numeric(as.character(prob))), model_extent_raster_agg, field = "prob", fun=mean, na.rm = TRUE)
    
    species_raster
    
  }) %>% set_names(algorithms)
  
  return(raster_predictions) 
  
}

validation_map <- function(predictions_to_map = species_model_predictions_categorical$rf$predictions_map_categorical,
                           presence_polys = presence_polygons, 
                           presence_pts = presence_points,
                           presence_observation_polys = NULL){
  
  tmap_mode("view")
  
  if ("RasterLayer" %in% class(predictions_to_map)){
    predictions_to_map <- predictions_to_map %>% 
      raster::mask(mask = predictions_to_map, maskvalue = 4) %>% 
      raster::trim()
    
    model_palette <- c("#F40000", "#EB8E04", "#EFEC02")
    names(model_palette) <- c("High", "Medium", "Low")
    
    m <- tm_basemap(leaflet::providers$Esri.WorldTopoMap, alpha = 0.7) +
      tm_shape(presence_polys) + tm_borders(col = "blue") +
      tm_shape(presence_pts) + tm_dots(col = "blue") +
      tm_shape(shp = predictions_to_map) + tm_raster(palette = model_palette, alpha = 0.5, style = "fixed", breaks = c(1, 2, 3, 4), labels = c("High", "Medium", "Low"), title = "Habitat probability", zindex = 2000) +
      tm_layout()
    
    if (!is.null(presence_observation_polys)) m <- m + tm_shape(presence_observation_polys) + tm_borders(col = "blue")
  }
  
  if ("sf" %in% class(predictions_to_map)){
    
    predictions_to_map <- predictions_to_map %>% 
      dplyr::filter(cat != 4) %>% 
      dplyr::mutate(cat = as.numeric(as.character(cat)))
    
    model_palette <- c("#F40000", "#EB8E04", "#EFEC02")
    names(model_palette) <- c("High", "Medium", "Low")
    
    m <- tm_basemap(leaflet::providers$Esri.WorldTopoMap, alpha = 0.7) +
      tm_shape(presence_polys) + tm_borders(col = "blue") +
      tm_shape(presence_pts) + tm_dots(col = "blue") +
      tm_shape(shp = predictions_to_map) + 
      tm_lines(col = "cat", palette = model_palette, alpha = 1, style = "fixed", breaks = c(1, 2, 3, 4), labels = c("High", "Medium", "Low"), title.col = "Habitat probability", zindex = 2000) +
      tm_layout()
    
    if (!is.null(presence_observation_polys)) m <- m + tm_shape(presence_observation_polys) + tm_borders(col = "blue")
  }
  
  
  return(m)
  
}

map_predictions <- function(predictions_to_map = species_model_predictions_categorical$rf$predictions_map_categorical){
  
  tmap_mode("plot")
  
  states <- USAboundaries::us_states() %>% 
    dplyr::filter(!(state_name %in% c("Hawaii", "Alaska", "Puerto Rico"))) %>% 
    st_set_crs(4326)
  
  if ("RasterLayer" %in% class(predictions_to_map)){
    ## Get predictions extent
    predictions_to_map <- predictions_to_map %>% 
      raster::mask(mask = predictions_to_map, maskvalue = 4) %>% 
      raster::trim()
    model_palette <- c("#F40000", "#EB8E04", "#EFEC02")
    names(model_palette) <- c("High", "Medium", "Low")
    model_bbox <- st_bbox(predictions_to_map %>% raster::extend(c(300, 300)))
    model_bbox <- model_bbox %>% st_as_sfc()
    model_bbox_WGS84 <- st_bbox(predictions_to_map %>% raster::extend(c(300, 300)) %>% projectRaster(crs = "+init=EPSG:4326"))
    model_bbox_WGS84 <- model_bbox_WGS84 %>% st_as_sfc()
    
    #' #### Panel A: Model predictions with inset
    set_defaults(map_service = "esri", map_type = "world_imagery")
    #' ##### Model predictions
    model_basemap <- basemap_raster(model_bbox)
    
    model_predictions_map = 
      tm_shape(model_basemap) + 
      tm_rgb(alpha = 0.9) +
      tm_shape(states[st_intersects(states, model_bbox_WGS84, sparse = FALSE),]) +
      tm_borders(col = "white", lwd = 1, zindex = 500, alpha = 0.5) +
      tm_layout(frame = FALSE, legend.show = TRUE, legend.text.color = "white", legend.title.color = "white", bg.color = "transparent", legend.text.size = 1, legend.title.size = 1) +
      tm_shape(shp = predictions_to_map) + 
      tm_raster(palette = model_palette, alpha = 1, style = "fixed", breaks = c(1, 2, 3, 4), labels = c("High", "Medium", "Low"), title = "Habitat probability", zindex = 2000)
  }
  
  if ("sf" %in% class(predictions_to_map)){
    ## Get predictions extent
    predictions_to_map <- predictions_to_map %>% 
      dplyr::filter(cat != 4) %>% 
      dplyr::mutate(cat = as.numeric(as.character(cat)))
    model_palette <- c("#F40000", "#EB8E04", "#EFEC02")
    names(model_palette) <- c("High", "Medium", "Low")
    model_bbox <- st_bbox(predictions_to_map)
    model_bbox <- model_bbox %>% st_as_sfc()
    model_bbox_WGS84 <- st_bbox(predictions_to_map %>% st_transform(crs = "+init=EPSG:4326")) %>% st_as_sfc()
    
    #' #### Panel A: Model predictions with inset
    set_defaults(map_service = "esri", map_type = "world_imagery")
    #' ##### Model predictions
    model_basemap <- basemap_raster(model_bbox)
    
    model_predictions_map = 
      tm_shape(model_basemap) + 
      tm_rgb(alpha = 0.9) +
      tm_shape(states[st_intersects(states, model_bbox_WGS84, sparse = FALSE),]) +
      tm_borders(col = "white", lwd = 1, zindex = 500, alpha = 0.5) +
      tm_layout(frame = FALSE, legend.show = TRUE, legend.text.color = "white", legend.title.color = "white", bg.color = "transparent", legend.text.size = 0.7, legend.title.size = 1) +
      tm_shape(shp = predictions_to_map) + 
      tm_lines(col = "cat", palette = model_palette, alpha = 1, style = "fixed", breaks = c(1, 2, 3, 4), labels = c("High", "Medium", "Low"), title.col = "Habitat probability", zindex = 2000)
  }
  
  model_predictions_grob <- tmap_grob(model_predictions_map)
  
  #' ##### Inset map
  inset_map = 
    tm_shape(shp = states[st_intersects(states, model_bbox_WGS84, sparse = FALSE),]) + tm_polygons(border.col = "black", col = "white", lwd = 3) +
    tm_shape(st_centroid(states[st_intersects(states, model_bbox_WGS84, sparse = FALSE),])) + tm_text(text = "state_abbr") +
    tm_shape(shp = model_bbox_WGS84) + tm_polygons(border.col = "red", col = "white", alpha = 0, lwd = 5) +
    tm_layout(frame = FALSE, outer.margins = c(2, 0, 0, .1)) 
  inset_map_grob <- tmap_grob(inset_map)
  
  habitat_map <- ggarrange(inset_map_grob, ggplot() + theme_void(), model_predictions_grob, ncol = 3, nrow = 1, widths = c(.2, .11, .69)) +
    theme(plot.margin = unit(c(20, 20, 20, 20), "points"))
  
  return(habitat_map)
  
}

map_model_refinement <- function(algo = "rf"){
  
  m2 <- map_predictions(predictions_to_map = species_model_predictions_categorical[[algo]]$predictions_map_categorical) 
  
  m2 <- m2 +
    annotate("text", x = 0.2, y = 1, label = paste("Refined model", species_model_objects[[algo]]$model_inputs_table$model_run_name[1]))
  
  if (domain == "terrestrial"){
    m1 <- map_predictions(predictions_to_map = paste0(rootdir, "/output/", taxon_code, "/", parent_model, "/", parent_model, "_categorical.tif") %>% raster())
    
  }
  
  if (domain == "aquatic"){
    m1 <- map_predictions(predictions_to_map = paste0(rootdir, "/output/", taxon_code, "/", parent_model, "/", parent_model, "_categorical.shp") %>% st_read())
    
  }
  
  m1 <- m1 +
    annotate("text", x = 0.2, y = 1, label = paste("Preliminary model", parent_model))
  
  
  return(ggarrange(m1, ggplot() + theme_void(), m2, ncol = 3, widths = c(.475, .05, .475)))
  
}

calculate_overall_confidence_score <- function(
  validation_score = d_validation_confidence_id, 
  review_score = d_review_confidence_id,
  methods_score = d_input_methods_confidence_id
){
  
  overall_confidence_score <- NULL
  if (validation_score == "high" & methods_score == "high" & review_score == "high") overall_confidence_score <- "high"
  if (validation_score == "high" & methods_score == "high" & review_score != "high") overall_confidence_score <- "medium"
  if (validation_score == "high" & methods_score == "medium" & review_score != "low") overall_confidence_score <- "medium"
  if (validation_score == "high" & methods_score == "medium" & review_score == "low") overall_confidence_score <- "low"
  if (validation_score == "high" & methods_score == "low" & review_score == "high") overall_confidence_score <- "medium"
  if (validation_score == "high" & methods_score == "low" & review_score != "high") overall_confidence_score <- "low"
  if (validation_score == "medium" & methods_score == "high" & review_score != "low") overall_confidence_score <- "medium"
  if (validation_score == "medium" & methods_score == "high" & review_score == "low") overall_confidence_score <- "low"
  if (validation_score == "medium" & methods_score == "medium" & review_score != "low") overall_confidence_score <- "medium"
  if (validation_score == "medium" & methods_score == "medium" & review_score == "low") overall_confidence_score <- "low"
  if (validation_score == "medium" & methods_score == "low") overall_confidence_score <- "low"
  if (validation_score == "low" & methods_score == "high" & review_score == "high") overall_confidence_score <- "medium"
  if (validation_score == "low" & methods_score == "medium") overall_confidence_score <- "low"
  if (validation_score == "low" & methods_score == "low") overall_confidence_score <- "low"
  
  return(overall_confidence_score)
  
}

generate_outputs <- function(domain = c("terrestrial", "aquatic"), pdf_only = FALSE){
  
  domain <- match.arg(domain)
  
  #' ## Run through all model runs
  purrr::map(1:length(species_model_objects), function(mod){
    
    #' ## Identify algorithm
    algorithm = names(species_model_objects)[mod]
    
    #' ## Generate map model predictions
    species_habitat_map <- map_predictions(predictions_to_map = species_model_predictions_categorical[[algorithm]]$predictions_map_categorical)
    
    #' ## Generate variable importance table
    variable_importance_table <- species_model_objects[[algorithm]]$full_model$importance[, "MeanDecreaseAccuracy"] %>% stack() %>% set_names("importance", "variable") %>% 
      dplyr::filter(importance >= quantile(importance, .5)) %>% 
      dplyr::left_join(env_predictor_metadata %>% dplyr::select(varName, fullName, description), by = c("variable" = "varName"))
    variable_importance_table <- variable_importance_table %>% 
      dplyr::mutate(fullName = factor(fullName, levels = variable_importance_table %>% dplyr::arrange(importance) %>% dplyr::pull(fullName)))
    
    #' ## Generate Review Ratings plot, if there are reviews for this model
    if (nrow(overall_feedback_table) > 0){
      ratings <- overall_feedback_table %>% 
        dplyr::pull(Rating) 
      p2 <- ggplot() + 
        geom_histogram(data = ratings %>% as.data.frame(), 
                       aes(x = .), 
                       bins = 5, 
                       binwidth = .5,
                       fill = "orange"
        ) + 
        xlim(1, 5) + 
        coord_flip() +
        theme(axis.text = element_text(size = 12), 
              axis.title.x = element_text(size = 12),
              axis.title.y = element_text(size = 12),
              panel.background = element_blank(),
              axis.ticks = element_blank()
        ) +
        xlab("Star Rating") +
        ylab("Number of Reviews")
      p1 <- ggplot() +
        theme_void() +
        geom_text(aes(3, 4.3, label = "Median Rating"), size = 6) +
        geom_text(aes(3, 3.5, label = paste0(ceiling(median(ratings)), "/5")), size = 12) +
        geom_text(aes(3, 2.9, label = paste0("(", length(ratings), " reviews)")), size = 6) +
        ylim(1, 5) +
        xlim(1, 5)
      reviews_p <- ggarrange(p1, p2, widths = c(0.3, 0.7))
    }
    
    #' ## Generate ROC thresholds plot
    roc_p <- data.frame(Sensitivity = species_model_predictions_categorical[[algorithm]]$roc_curve$sensitivities, Specificity = 1 - species_model_predictions_categorical[[algorithm]]$roc_curve$specificities) %>% 
      ggplot(aes(x = Specificity, y = Sensitivity)) +
      geom_line(size = 1) +
      theme_bw() +
      theme(panel.grid = element_blank()) +
      geom_point(data = data.frame(x = species_model_predictions_categorical[[algorithm]]$thresholds.xy["low.x"], y = species_model_predictions_categorical[[algorithm]]$thresholds.xy["low.y"]), aes(x = x, y = y), size = 3, color = "#EFEC02") +
      geom_point(data = data.frame(x = species_model_predictions_categorical[[algorithm]]$thresholds.xy["medium.x"], y = species_model_predictions_categorical[[algorithm]]$thresholds.xy["medium.y"]), aes(x = x, y = y), size = 3, color = "#EB8E04") +
      geom_point(data = data.frame(x = species_model_predictions_categorical[[algorithm]]$thresholds.xy["high.x"], y = species_model_predictions_categorical[[algorithm]]$thresholds.xy["high.y"]), aes(x = x, y = y), size = 3, color = "#F40000") +
      geom_vline(xintercept = species_model_predictions_categorical[[algorithm]]$thresholds.xy["low.x"], linetype = "dashed", color = "#EFEC02") +
      geom_vline(xintercept = species_model_predictions_categorical[[algorithm]]$thresholds.xy["medium.x"], linetype = "dashed", color = "#EB8E04") +
      geom_vline(xintercept = species_model_predictions_categorical[[algorithm]]$thresholds.xy["high.x"], linetype = "dashed", color = "#F40000") +
      geom_abline(intercept = 0, slope = 1, color = grey(0.6)) +
      xlab("False Positive Rate") +
      ylab("True Positive Rate")
    
    #' ## Specify metadata objects 
    common_name = species_info$primaryCommonName
    scientific_name = species_info$scientificName
    grank = species_info$roundedGRank
    element_unique_id = species_info$uniqueId
    model_version = species_model_objects[[algorithm]]$model_summary_table$model_run_name[1]
    algorithm = names(species_model_objects)[1]
    algorithm_name = "Random Forest"
    model_creation_date = paste0(substr(model_version, nchar(model_version)-7, nchar(model_version)-4), "-", substr(model_version, nchar(model_version)-3, nchar(model_version)-2), "-", substr(model_version, nchar(model_version)-1, nchar(model_version)))
    model_iteration = 1
    d_input_methods_confidence_id = c("low", "medium", "high")[2]
    d_validation_confidence_id = ifelse(species_model_objects[[algorithm]]$model_summary_table %>% dplyr::filter(metric == "AUC") %>% dplyr::pull(metric_mn) < 0.8, "low",
                                        ifelse(species_model_objects[[algorithm]]$model_summary_table %>% dplyr::filter(metric == "AUC") %>% dplyr::pull(metric_mn) >= 0.9, "high", "medium")
    )
    d_review_confidence_id = ifelse(nrow(overall_feedback_table) == 0, "low",
                                    ifelse(ceiling(median(ratings)) >= 4, "high", 
                                           ifelse(ceiling(median(ratings)) >= 3 & ceiling(median(ratings)) < 4, "medium")
                                    )
    )
    #d_review_confidence_id ="low" #PJM added
    presence_data_quality = "Acceptable"
    absence_background_data_quality = "Acceptable"
    predictor_relevance_quality = "Acceptable"
    spatio_temporal_data_alignment_quality = "Acceptable"
    algorithm_choice_quality = "Acceptable"
    sensitivity_analysis_quality = "Acceptable"
    d_overall_confidence_id = calculate_overall_confidence_score()
    contact = "gio@natureserve.org"
    model_thresholds = round(species_model_predictions_categorical[[algorithm]]$thresholds, 5)
    variable_importance = species_model_objects[[algorithm]]$var_importance
    partial_plots = species_model_objects[[algorithm]]$partial_plots
    model_inputs_table = species_model_objects[[algorithm]]$model_inputs_table
    model_summary_table = species_model_objects[[algorithm]]$model_summary_table
    model_thresholds.xy = species_model_predictions_categorical[[algorithm]]$thresholds.xy
    period = "current"
    organization = "NatureServe"
    
    #' # Save Output
    #' ## Create output directory
    dir.create(paste0(rootdir, "output/", taxon_code))
    dir.create(paste0(rootdir, "output/", taxon_code, "/", model_version))
    #' ## Knit metadata
    #' ### Set working directory
    setwd(paste0(rootdir, "metadata"))
    knit2pdf(input = "SHM_metadata.rnw", output = paste0(model_version, ".tex"))   
    if (file.exists(paste0(model_version, ".pdf"))) file.remove(paste(model_version, c(".log", ".aux", ".out", ".tex"), sep=""))
    if (file.exists(paste0(model_version, ".pdf"))) file.copy(from = paste0(model_version, ".pdf"), to = paste0(rootdir, "output/", taxon_code, "/", model_version, "/", model_version, ".pdf"))  
    setwd(rootdir)
    
    if (isFALSE(pdf_only)){
      #' ## Save maps
      if (domain == "terrestrial"){
        #' ### Continuous raster
        writeRaster(species_model_predictions[[algorithm]]$predictions_map_continuous, filename = paste0(rootdir, "/", "output/", taxon_code, "/", model_version, "/", model_version, "_continuous.tif"), format="GTiff", overwrite = TRUE)
        #' ### Categorical raster
        r_thres <- species_model_predictions_categorical[[algorithm]]$predictions_map_categorical
        #' #### Edit attribute table
        names(levels(r_thres)[[1]])[2] <- c("hab_value")
        levels(r_thres)[[1]]$hab_class <- c("high probability", "medium probability", "low probability", "non-habitat")
        writeRaster(r_thres, filename = paste0(rootdir, "/", "output/", taxon_code, "/", model_version, "/", model_version, "_categorical.tif"), format="GTiff", datatype='INT4S', overwrite = TRUE)
        write.dbf(levels(r_thres)[[1]], file = paste0(rootdir, "/", "output/", taxon_code, "/", model_version, "/", model_version, "_categorical.tif.vat.dbf"))
      } 
      
      if (domain == "aquatic"){
        #' ### Continuous shapefile
        st_write(species_model_predictions[[algorithm]]$predictions_map_continuous, dsn = paste0(rootdir, "/", "output/", taxon_code, "/", model_version, "/", model_version, "_continuous.shp"))
        #' ### Continuous raster
        #writeRaster(species_model_predictions_rasters[[algorithm]], filename = paste0(rootdir, "/", "output/", taxon_code, "/", model_version, "/", model_version, "_continuous.tif"), format="GTiff", overwrite = TRUE)
        #' ### Categorical shapefile
        st_write(species_model_predictions_categorical[[algorithm]]$predictions_map_categorical, dsn = paste0(rootdir, "/", "output/", taxon_code, "/", model_version, "/", model_version, "_categorical.shp"))
      }
      
      #' ### Package folder for MRT upload
      #' #### Create output folder
      dir.create(paste0(rootdir, "output/", taxon_code, "/", model_version, "/forMRT"))
      dir.create(paste0(rootdir, "output/", taxon_code, "/", model_version, "/forMRT/rasterToUpload"))
      #' #### Include metadata pdf
      file.copy(from = paste0(rootdir, "output/", taxon_code, "/", model_version, "/", model_version, ".pdf"), to = paste0(rootdir, "output/", taxon_code, "/", model_version, "/forMRT/", model_version, ".pdf"))  
      #' #### Include coninuous predictions raster
      file.copy(from = paste0(rootdir, "output/", taxon_code, "/", model_version, "/", model_version, "_continuous.tif"), to = paste0(rootdir, "output/", taxon_code, "/", model_version, "/forMRT/rasterToUpload/", model_version, "_continuous.tif"))
      #' #### Create additional files for upload to MRT
      write_csv(species_hucs %>% dplyr::select(HUC10) %>% st_set_geometry(value = NULL) %>% as.data.frame() %>% dplyr::transmute(HUC_10 = HUC10), paste0(rootdir, "output/", taxon_code, "/", model_version, "/forMRT/modeling-extent.csv"))
      
      #' #### Create rjson file for upload to MRT
      model_json <- list(modelVersion = model_version,
                         iteration = as.numeric(model_iteration),
                         iterationNote = "null",
                         thresholdLow = as.numeric(model_thresholds["low"]),
                         thresholdMid = as.numeric(model_thresholds["medium"]),
                         thresholdHigh = as.numeric(model_thresholds["high"]),
                         elementGlobalId = species_info$elementGlobalId,
                         cuteCode = model_version,
                         scientificName = scientific_name,
                         commonName = common_name,
                         taxonomicGroup = species_info$speciesGlobal$informalTaxonomy$level3,
                         algorithm = algorithm,
                         locationUseClass = ""
      )
      RJSONIO::toJSON(model_json, pretty = TRUE) %>% write(paste0(rootdir, "output/", taxon_code, "/", model_version, "/forMRT/modelRunInfo.json"))
      
      #' Save species R script
      file.copy(from = paste0(rootdir, "SHM_workflow.R"), to = paste0(rootdir, "output/", taxon_code, "/", model_version, "/", model_version, "-SHM_workflow.R"))
      #' Save .RData file
      save.image(paste0(rootdir, "output/", taxon_code, "/", model_version, "/", model_version, ".RData"))
    }
  })
  
}
