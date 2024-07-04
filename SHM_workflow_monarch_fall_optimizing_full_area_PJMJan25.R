#' ---
#' title: NatureServe Species Habitat Modeling Workflow
#' ---

#' # Set up model run (Update fields as necessary) ##########
#' ## Specify root directory
#'
scriptName<-"SHM_workflow_monarch_fall_optimizing_full_area_PJMJan25.R"

rootdir <- paste(getwd(),"/",sep="") # Update root directory
#' ## Create useful objects
species_name <- "Danaus plexipus" # Update with scientific name of focal taxon
element_unique_id <- "ELEMENT_GLOBAL.2.118751" # Update with element global ID of focal taxon
#' ## Get species info from NatureServe Explorer API
species_info <- if(!is.na(element_unique_id)){
  natserv::ns_id(uid = element_unique_id)
} else {
  natserv::ns_search_spp(text_adv = list(searchToken = species_name, matchAgainst = "scientificName", operator = "equals"))$results
}
#' ## Specify taxon code
taxon_code <- paste0(substr(species_info$scientificName, 1, 1), substr(strsplit(species_info$scientificName[1], " ")[[1]][2], 1, 4), "_", species_info$elementGlobalId)
#' ## Select domain for the model: terrestrial or aquatic
domain <- "terrestrial" # select "terrestrial" if modeling terrestrial species or "aquatic" if modeling aquatic species
BLD_type <- "SF" # select "EO" if building model from element occurrences or "SF" if building model from source features
# Identify parent model if this is a model iteration. 
#Change to "" if this is a first-generation mode
parent_model <- ""    #"Nlewi_100341_rf_20220329"

#############################################################

#' # Load model reviews for parent model, if any ############
#' ## Read in MRT tables from AGOL into R
#' ### Check ArcGIS Pro product is installed 
library(arcgisbinding)
arcgisbinding::arc.check_product()
#' ### Check that the connection is successful
arcgisbinding::arc.check_portal()

## Read in MRT tables
### Read in OverallFeedbackRaster table

overall_feedback_table <- arc.open("https://services.arcgis.com/EVsTT4nNRCwmHNyb/arcgis/rest/services/OverallFeedbackRaster/FeatureServer/0") %>% 
  arc.select() %>% 
  dplyr::filter(ModelVersion == parent_model) 
### Read in OverallFeedbackRaster table
detailed_feedback_table <- arc.open("https://services.arcgis.com/EVsTT4nNRCwmHNyb/arcgis/rest/services/DetailedFeedbackRaster/FeatureServer/0") %>% 
  arc.select() %>% 
  dplyr::filter(ModelVersion == parent_model)

#############################################################

#' # Load and create data inputs  ###########################
#' ## Load nationwide HUC10 data layer
huc10_layer <- st_read("S:/other_spatial_layers/_data/other_spatial/feature/HUC10.shp")
#' ## Load EO data
presence_polygons <- st_zm(st_read(paste0(rootdir, "data/monarchs/monarch_template_8_10.shp"), quiet = TRUE)) %>% 
  dplyr::filter(GNAME == species_name) %>% 
  dplyr::filter(EOSF == BLD_type) %>% 
  st_transform(st_crs(huc10_layer))

#' ## Load observation data, if any
external_data<-"yes" # "yes" if there is addtional data
if (external_data=="yes"){
  presence_observations <- read_csv(paste0(rootdir, "data/monarchs/MergedMonarchs_202212_PJM_select_ras.csv"))
    presence_observations<-presence_observations%>%
    #sample_n(1000) %>% #pjm added for testing/bias
    dplyr::filter(coastal!=1)%>%
    #dplyr::filter(fall_coast!=1)%>% # excludes records within 20km of wintering sites
    dplyr::filter(Month>7 & Month<11)%>%
    #dplyr::filter(NV_AZ_win!=1)%>%
    #dplyr::filter(excl_east!=1)%>%
    #dplyr::filter(spr_outl!=1)%>%
    dplyr::filter(Year>2009)%>%
    dplyr::filter(CNDDB_win!=1)%>%   
    dplyr::filter(wRockies==1)%>%  
    #dplyr::filter(!is.na(RASTERVALU_ext))%>% 
    #dplyr::filter(orig_Individuals==1| is.na(orig_Individuals))%>%
    dplyr::filter(!Life_Stage%in%c("Larva","Egg" , "Pupa", "Larva/Egg" , "Pupa/Larva",   "Pupa/Larva/Egg","Egg/larvae" ,"<Null>", "Unknown"))%>%
    dplyr::filter(!is.na(Life_Stage))%>%
    sample_frac(.7,replace = FALSE, weight =bias )%>%
        st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
    st_transform(st_crs(huc10_layer))
  presence_observation_polygons <- group_observations_to_polygons(presence_observations = presence_observations)
}else if (external_data=="no"){
  presence_observations<- NULL  
  presence_observations_polygons<- NULL  
}
#' ## Generate model extent from data
#' ### Specify if given hucs need to be added
#if there is MIRTL input 

#' hucs_to_add <- detailed_feedback_table %>% dplyr::filter(StatusType == 1) %>% dplyr::pull(HUCID)
#' #' ### Specify if given hucs need to be removed
hucs_to_remove <- detailed_feedback_table %>% dplyr::filter(StatusType == 2) %>% dplyr::pull(HUCID)
#' 
sp_hucs<-read_csv(paste0(rootdir, "data/monarchs/monarchs_wRockies_huc10_spr_west.csv"))
#sp_hucs<-read_csv(paste0(rootdir, "data/monarchs/monarchs_wRockies_huc10_noCI.csv"))
hucs_to_add<-sp_hucs$HUC10
#' ### Identify hucs for model extent from data plus any additions/removals
species_hucs <- get_model_extent_from_data(element_occurrences = presence_observation_polygons, #PJM changed from presence_polygons, 
           hucs = huc10_layer,
           add_hucs = hucs_to_add,
           remove_hucs = hucs_to_remove
)


###### If the model inputs been reviewed, following code will include that
# You would need access to applied data science team googledrive or change the googlesheets link in the functions
cutecode<- paste(strsplit(species_name," ")[[1]][1],strsplit(species_name," ")[[1]][2],sep='_')
mirtl_reviewed<- "no"
if (mirtl_reviewed =="yes"){
  species_hucs <-refine_model_extent_from_MIRTL(cutecode =cutecode, current_extent = species_hucs, hucs=huc10_layer)
}else{
  species_hucs<-species_hucs
}
#' ### Map model extent
model_extent_map <- map_model_extent()
model_extent_map
#' ## Specify potential environmental predictors
#' ### Load relevant environmental predictor information
#' #### Terrestrial species workflow
if (domain == "terrestrial"){
  # Specify path
  environmental_predictors_path <-   "P:/PJM/terrestrialEnvVars_330_crop/"  #"S:/terrestrialEnvVars_330/" #"S:/terrest_envV_2/"
  # Get file names
  env_predictor_files <- list.files(environmental_predictors_path, pattern = ".tif$", recursive = TRUE, full.names = TRUE)
  # Load metadata for environmental predictor set
  env_predictor_metadata <- read_csv(paste0(rootdir, "data/environmental-predictor-library-metadata-", domain, ".csv")) %>% 
    dplyr::mutate(varName = gsub(".tif", "", fileName))
}
#' #### Aquatic species workflow
if (domain == "aquatic"){
  # Specify path
  environmental_predictors_path <- "S://aquaticEnvVars_NHD/background.sqlite"
  # Load metadata for environmental predictor set
  env_predictor_metadata <- read_csv(paste0(rootdir, "data/environmental-predictor-library-metadata-", domain, ".csv")) %>% 
    dplyr::mutate(varName = gridName)
  # Identify environmental predictors
  env_predictor_set <- env_predictor_metadata %>% dplyr::pull(varName)
}
#' ## Clean up
#' ### Remove nationwide HUC10 data layer to save space
rm(huc10_layer)
#############################################################

#' ##### NOTE: This is where Model Inputs Review Tool will get spun up

#' # Process model inputs ###################################
#' ## Dissolve species hucs to get model extent
model_extent <- st_union(species_hucs)
#' ## Process species occurrence data
#' ## Terrestrial species workflow
if (domain == "terrestrial"){
  #' ### Sample presence points from presence polygons
  #' #### Sample presence points from presence polygons
  presence_points <- get_points_from_presence_polygons(presPolys = presence_polygons, presObs = presence_observation_polygons)
  #' ### Sample background points
  #' #### Generate random points from study area
  background_points <- st_sample(model_extent, size = 100000) ### Set number of points to generate PJM changed to 250000
  background_points <- st_sf(fid = 1:length(background_points), geometry = background_points)
  background_points <- st_join(background_points, species_hucs, join = st_intersects)[c("fid", "HUC10")]
  names(background_points) <- c("fid", "huc10", "geometry")
  #' ### Update background points
  #' #### Remove background points close to presence points
  #' ##### Buffer presence points
  presence_points_buffered <- st_buffer(presence_points, dist = 3000) #pjm- charnged from 100 to make more sparse
  #' ##### Identify overlapping background points
  coincident_points <- unlist(st_contains(presence_points_buffered, background_points, sparse = TRUE))
  #' ##### Update background_points object
  background_points <- background_points[setdiff(background_points$fid, coincident_points), ] 
  #' #### Subset background points to appropriate size
  #' ##### Use the greater of 10 x number of presence points
  background_points <- background_points %>% sample_n(nrow(presence_points) * 8) #pjm changed from 10
  print(paste0("number of background points: ", nrow(background_points)))
} else if (domain == "aquatic"){ #' ## Aquatic species workflow
  
  aquatic_occurrence_inputs <- generate_aquatic_occurrence_inputs()
  presence_reaches <- aquatic_occurrence_inputs$presence_reaches
  background_reaches <- aquatic_occurrence_inputs$background_reaches
  model_reaches <- aquatic_occurrence_inputs$model_reaches
}
#' ## Process environmental predictor data
#' ## Attribute presence and background points
#' ## Terrestrial species workflow
if (domain == "terrestrial"){
  #' ### Tailor environmental predictors for this species
  #' #### By exclusion
  
  #exclude_vars = ""
  exclude_vars = c("z3impsur1", "z3impsur10", "z3impsur1c", "impsur1", "impsur10", "impsur100",   "RivED330", 
  "FwPond330", "ForWet330", "HerED330", "MilkMax330")
  #, "FWEmer330", "z3upslwr10",  "z3upslpwr1",  "z3upslop10",  "z3upslope1",  "z3upslpc10",  "z3upslpcl1",
                     #"z3upslfl10",  "z3upslpfl1")
  
  
  # e.g. exclude_vars = c("nlcdopn1", "nlcdopn10", "nlcdopn100", "impsur1", "impsur10", "impsur100",
  #  "ntm_1_01", "ntm_1_02", "ntm_1_06", "ntm_1_08", "ntm_1_09", "ntm_2_01",
  #  "ntm_2_02", "ntm_2_05", "ntm_2_06", "ntm_3_01", "ntm_3_03", "ntm_3_09",
  #  "ntm_3_12", "ntm_4_01", "ntm_4_02", "ntm_4_03", "ntm_4_05", "ntm_4_06",
  #  "ntm_5_01", "ntm_6_01", "ntm_6_02", "ntm_6_03", "ntm_6_04", "nlcdshb1",
  #  "nlcdshb10", "nlcdshb100", "distocean"
  # )
  if (exclude_vars[1]!= ""){   #pjm added [1] since was getting indexing error
    env_predictor_files <- env_predictor_files[!grepl(paste0(env_predictor_metadata %>% dplyr::filter(gridName %in% exclude_vars) %>% dplyr::pull(fileName), collapse = "|"), env_predictor_files)]
  }
  #' #### By inclusion
  include_vars = ""
  # include_vars =c("z3clicwd", "z3bulkdens",  "RivED330",  "OthED330",  "HerED330",  "FwPond330",  "ForWet330",
  # "FWEmer330",  "LakeED330",  "NectCt330")
  # 
  
  # include_vars = c("z3clay", "NectCt330", "FWEmer330", "z3poresize","z3soil_ph", "z3clicwd",
  #  "z3cligdd",  "z3clibio10","z3clibio11",  "z3clibio12",
  #  "z3clibio16",  "z3clibio17",  "z3clibio1",  "z3distocea", "z3diststrp","z3flowacc",
  #  "z3geocarb", "z3geoncarb", "z3geoalk",  "z3geosilic", "z3geoextru",   "z3geocollu",   "z3geoclay",
  #  "z3geoloam",   "z3geoctill",   "z3geocoutw",   "z3geofsed",   "z3geopeat",   "z3geosdune",
  #  "z3geoloess",   "z3geosalsd",   "z3geofcsed",   "z3geoccsed",   "z3canopy1",   "z3canopy10",
  #  "z3canopy1c",   "z3nlcdcf1c",   "z3nlcdcf10",   "z3nlcdcfr1",   "z3nlcddf1c",   "z3nlcddf10",
  #  "z3nlcddfr1",   "z3nlcdopn1",   "z3nlcdop10",   "z3nlcdop1c",   "z3nlcdshb1",   "z3nlcdsh10",
  #  "z3nlcdsh1c",   "z3nlcdwat1",   "z3nlcdwa10",   "z3nlcdwa1c",   "z3nlcdww1c",   "z3nlcdww10",
  #  "z3nlcdwwt1",   "z3dnwiffw",   "z3dnwifemw",   "z3dnwisemw",   "z3tp001x1k",   "z3tp010x1k",
  #  "z3tp100x1k",   "z3bulkdens")
  # 
  if (include_vars[1]!= ""){    #pjm added [1] since was getting indexing error
    env_predictor_files <- env_predictor_files[grepl(paste0(env_predictor_metadata %>% dplyr::filter(gridName %in% include_vars) %>% dplyr::pull(fileName), collapse = "|"), env_predictor_files)]
  }
  #' ### Extract environmental values
  #' #### identify environmental predictor files
  environmental_values <- env_predictor_files %>% 
    furrr::future_map(function(r){
      env_raster <- stars::read_stars(r) 
      names(env_raster) <- gsub(".tif", "", names(env_raster))
      presence_values <- stars::st_extract(env_raster, presence_points) %>% 
        sf::st_set_geometry(value = NULL)
      background_values <- stars::st_extract(env_raster, background_points) %>% 
        sf::st_set_geometry(value = NULL)
      list(presence_values = presence_values, background_values = background_values)
    }) 
  #' ### Run environmental predictor intercorrelation routine
  environmental_values_vif <- purrr::map(environmental_values, 2) %>% bind_cols() %>% usdm::vifcor(.9)
  #' ### Isolate predictors for the focal species
  env_predictor_set <- environmental_values_vif@results$Variables
  #' ### Attribute presence and background points
  presence_points <- cbind(presence_points, purrr::map(environmental_values, 1) %>% bind_cols() %>% dplyr::select(env_predictor_set))
  background_points <- cbind(background_points, purrr::map(environmental_values, 2) %>% bind_cols() %>% dplyr::select(env_predictor_set))
} else if (domain == "aquatic"){   #' ## Aquatic species workflow 
  #' ## Attribute presence reaches
  aquatic_predictors_db <- dbConnect(SQLite(), dbname = "S://aquaticEnvVars_NHD/background.sqlite")
  #' ### Write query
  presence_reaches_query <- paste0("SELECT * FROM ", "background_reaches", "_att WHERE COMID IN ('", paste(presence_reaches$COMID, collapse = "','"),"')") 
  #' ### Get query
  presence_reaches_env <- dbGetQuery(aquatic_predictors_db, presence_reaches_query)
  #' ### Merge data
  presence_reaches <- presence_reaches %>%
    dplyr::select(COMID, GNIS_ID, GNIS_NAME, LENGTHKM, REACHCODE, FLOWDIR, WBAREACOMI, FTYPE, FCODE, group_id) %>% 
    dplyr::left_join(presence_reaches_env, by = "COMID") %>% 
    dplyr::mutate(VBL_RL_R = as.numeric(VBL_RL_R),
  VBA_RWA_R = as.numeric(VBA_RWA_R),
  exp = as.numeric(exp),
  JulAug_tempC = as.numeric(JulAug_tempC),
  Flow_cfs = as.numeric(Flow_cfs),
  mean_diam = as.numeric(mean_diam),
  DCI_up = as.numeric(DCI_up),
  DOR = as.numeric(DOR)
    )
  #' ## Attribute background reaches
  #' ### Write query
  background_reaches_query <- paste0("SELECT * FROM ", "background_reaches","_att WHERE COMID IN ('", paste(background_reaches$COMID, collapse = "','"),"')") 
  #' ### Get query
  background_reaches_env <- dbGetQuery(aquatic_predictors_db, background_reaches_query)
  #' ### Merge data
  background_reaches <- background_reaches %>%
    dplyr::select(COMID, GNIS_ID, GNIS_NAME, LENGTHKM, REACHCODE, FLOWDIR, WBAREACOMI, FTYPE, FCODE) %>% 
    dplyr::left_join(background_reaches_env, by = "COMID") %>% 
    dplyr::mutate(VBL_RL_R = as.numeric(VBL_RL_R),
  VBA_RWA_R = as.numeric(VBA_RWA_R),
  exp = as.numeric(exp),
  JulAug_tempC = as.numeric(JulAug_tempC),
  Flow_cfs = as.numeric(Flow_cfs),
  mean_diam = as.numeric(mean_diam),
  DCI_up = as.numeric(DCI_up),
  DOR = as.numeric(DOR)
    )
  #' ### Update presence_reaches object
  presence_reaches <- presence_reaches %>% 
    dplyr::mutate(stratum = as.numeric(factor(COMID, levels = unique(COMID))),
  uid = group_id,
  species_cd = paste0(substr(presence_polygons$GNAME[1], 1, 1), substr(strsplit(presence_polygons$GNAME[1], " ")[[1]][2], 1, 4), "_", presence_polygons$EGT_ID[1]),
  ra = "medium",
  obsdate = substr(Sys.Date(), 1, 4),
  date = obsdate
    )
  #' ### Tailor environmental predictors for this species
  #' #### By exclusion
  #' ##### Remove variables with variation = 0 and a few Pct variables which do not work
  exclude_vars = c(names(which(colSums(presence_reaches[env_predictor_set] %>% st_set_geometry(value = NULL), na.rm = TRUE) == 0)), "RdDensCatRp100", "PctImp2011Cat", "PctImp2011Ws", "PctImp2011CatRp100")
  env_predictor_set <- setdiff(env_predictor_set, exclude_vars)
  #' #### By inclusion
  include_vars = ""
  if (include_vars != "") env_predictor_set <- intersect(env_predictor_set, include_vars)
  
  #' ## Run environmental predictor intercorrelation routine
  #' ### Isolate environmental predictor variable names
  environmental_values_vif <- rbind(presence_reaches %>% dplyr::select(env_predictor_set) %>% st_set_geometry(value = NULL), background_reaches %>% dplyr::select(env_predictor_set) %>% st_set_geometry(value = NULL)) %>% usdm::vifcor(.9)
  #' ### Isolate predictors for the focal species
  env_predictor_set <- environmental_values_vif@results$Variables
  #' ### Attribute presence and background points
  presence_points <- presence_reaches %>% dplyr::select(-environmental_values_vif@excluded)
  background_points <- background_reaches %>% dplyr::select(-environmental_values_vif@excluded)
}
#############################################################

#' # Build species habitat models ###########################
species_model_objects <- run_models(presence_pts = presence_points, 
    background_pts = background_points, 
    species_prdictrs = env_predictor_set, 
    algorithms = "rf"
)
#############################################################

#' # Map model predictions  #################################
#' ## Project species habitat models across model extent
if (domain == "terrestrial"){
  species_model_predictions <- project_models()
} else if (domain == "aquatic"){
  species_model_predictions <- project_models_aquatic()
}
#' ## Threshold model predictions
species_model_predictions_categorical <- threshold_predictions(species_predictions = species_model_predictions, 
               domain = domain,
               sensitivity_level = .95, # set desired proportion of observed presences correctly predicted by low threshold 
               specificity_level = .95 # set desired proportion of observed (pseudo)absences correctly predicted by high threshold
)
#' ## Rasterize aquatic model predictions
if (domain == "aquatic") species_model_predictions_rasters <- rasterize_aquatic_predictions(predictions = species_model_predictions$rf$predictions_map_continuous, aggregation_coefficient = 1)
#' ## Visual check to validate model predictions before generating outputs
validation_map()
#############################################################


#' # Save output
library(ggpubr)
d_input_methods_confidence_id = c("low", "medium", "high")[2]
algorithm="rf"
d_validation_confidence_id = ifelse(species_model_objects[[algorithm]]$model_summary_table %>% dplyr::filter(metric == "AUC") %>% dplyr::pull(metric_mn) < 0.8, "low",
    ifelse(species_model_objects[[algorithm]]$model_summary_table %>% dplyr::filter(metric == "AUC") %>% dplyr::pull(metric_mn) >= 0.9, "high", "medium")
)
d_review_confidence_id = ifelse(nrow(overall_feedback_table) == 0, "low",
                ifelse(ceiling(median(ratings)) >= 4, "high", 
       ifelse(ceiling(median(ratings)) >= 3 & ceiling(median(ratings)) < 4, "medium")
                )
)





#' # Save output
generate_outputs(domain = domain)

#file.create("C:/_automation/ShutMeDown.txt")