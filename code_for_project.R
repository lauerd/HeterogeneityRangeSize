# Daniel Lauer
# 1 April 2020 - 8 February 2023

# Paper title: Habitat and not topographic heterogeneity constrains the range sizes of African mammals

##########################################################################################################################
##########################################################################################################################

# PART I: INTEGRATION OF DATA AND DATA CLEANING

# Set the working directory to the folder containing the data for this project. Also set a random seed to be used
# throughout the script, such that the effects of randomness are removed from all analyses:

setwd('./Data')
Random_seed <- 756; set.seed(Random_seed)

# Generate and read in a list of mammal species that occur only within Africa (hereon called Africa-only species):

  # Create a polygon of Africa's borders, taken from polygons of the world's countries provided by Natural Earth.
  # Process those borders by discerning whether or not they should include Madagascar, and by creating a re-projected
  # version of them for the "Africa_only" function below:

  library(raster); library(sp) # For general spatial analyses, to be used throughout.
  Africa <- shapefile('Shapefile_Africa/Africa.shp') # Shapefile of whole planet.
  Africa <- Africa[Africa@data$CONTINENT == 'Africa', ] # Subset to just Africa.
  plot(Africa, col = 'cornsilk1') # Plot to ensure it looks like Africa.
  
  No_madagascar <- TRUE # Indicate whether or not Madagascar should be included in Africa's borders.
  if (No_madagascar) { # If Madagascar is not to be included...
    Africa <- Africa[Africa@data$SOVEREIGNT != 'Madagascar', ] # Remove Madagascar from "Africa".
    plot(Africa, col = 'cornsilk1') } # Plot to ensure Madagascar was removed.
  
  Africa_samprast <- raster('Shapefile_Africa/Africa_SampleRaster.tif') # Sample raster of species for use below.
  Africa_diffproj <- spTransform(Africa, CRSobj = Africa_samprast@crs) # Re-project "Africa" to CRS of sample raster.
  plot(Africa_diffproj, col = 'cornsilk1') # Plot to ensure re-projection of "Africa" was successful.
  # Africa_diffproj <- rasterize(Africa_diffproj, Africa_samprast, getCover = TRUE) # Convert...
    # ..."Africa_diffproj" to a raster with the same dimensions as "Africa_samprast", and have each of the raster's cells
    # hold the % with which it overlaps with the original "Africa_diffproj" polygon. This will allow partially-
    # overlapping raster cells of species' ranges to be included within Africa's borders when assessing species in
    # "Africa_only" below, such that finer range detail can be attained. Save the result as a file and read it in next.
  # save(Africa_diffproj, file = 'Africa_diffproj.RData') # Save to circumvent inconsequential garbage-collection errors.
  load('Africa_diffproj.RData') # Load it/read it in.
  Africa_diffproj[Africa_diffproj == 0] <- NA # Convert all cells with 0% overlap to NA for the "mask" function below.
  plot(Africa_diffproj) # Plot to ensure rasterization was successful.
  
  # Define a function that takes a mammal species' raster range map as input, and outputs a TRUE/FALSE outcome dictating
  # if the species occurs only within Africa or not:
  
  Africa_only <- function(Range_map) { # For a given species' range map...
    
    # Calculate the sum of the map's pixels, then constrain the map to the borders of Africa, and perform the sum
    # calculation again. If the second divided by the first calculation is >= 0.9 (90%), such that >= 90% of the
    # associated species' range occurs within Africa, then the species can be said to likely occur only within Africa, 
    # as the 90% threshold accounts for map uncertainty. If <90%, the species likely occurs outside of Africa:
    
    Tmp <- cellStats(Range_map, sum) # First sum calculation.
    Range_map <- mask(Range_map, Africa_diffproj) # Constrain the range map to the borders of Africa.
    Tmp2 <- cellStats(Range_map, sum) # Second sum calculation.
    if (Tmp == 0 | Tmp2 / Tmp < 0.9) { return(FALSE) } else { return(TRUE) }} # Does species map only occur in Africa?
  
  # Apply "Africa_only" to every species in the "Current" folder (a folder downloaded from Phylacine's Dryad repository
  # and removed thereafter to save space). Save the resultant species list as a ".txt" file, denoting in the name if
  # Madagascar was included or not. Read the list in (this code block was run once and commented out, with its outputs
  # saved, to save time and space):
  
  # setwd('./Current') # Set the working directory to the appropriate folder.
  # Africa_only_log <- c() # Initiate a TRUE/FALSE vector indicating if a species is in Africa only or not.
  # for (Species in 1:length(list.files())) { # For each species in the folder...
  #   Tmp1 <- raster(list.files()[Species]); Africa_only_log[Species] <- Africa_only(Tmp1) # Apply "Africa_only".
  #   print(paste('Species', Species, 'complete')) }; rm(Species, Tmp1) # Print a record of this loop's progress.
  # setwd('../') # Restore the working directory to its previous location.
  # write.table(gsub('\\.tif', '', list.files('./Current')[Africa_only_log]), file = ifelse(No_madagascar,
  #   'AfricaMammalsNames_ByMe_NoMadagascar.txt', 'AfricaMammalsNames_ByMe_V2.txt'), row.names = FALSE,
  #   col.names = FALSE) # Save list.
  Mammals_names <- as.character(unlist(unname(read.table(list.files()[grepl(ifelse(No_madagascar, 'NoMadagascar',
    'Names_ByMe_V2'), list.files())])))) # Read list in.

# Read in and process a dataset from Phylacine to produce binomial (genus/species) synonyms for all Africa-only species:
    
Mammals_synonyms <- read.csv(list.files()[grepl('Synonyms', list.files())]) # Read in the synonymy dataset.
Mammals_synonyms <- Mammals_synonyms[match(Mammals_names, Mammals_synonyms$Binomial.1.2), ] # Subset to Africa-only.
Mammals_synonyms$Binomial1 <- Mammals_synonyms$Binomial.1.2 # Change name of binomial column.
Counter <- 2 # Initiate a counter that will be used to make additional binomial columns in "Mammals_synonyms".
for (NewCol in c('.1.1', '.1.0$', 'Elton', 'IUCN')) { # For each genus/species pair in "Mammals_synonyms"...
  Mammals_synonyms[paste('Binomial', Counter, sep = '')] <- apply(Mammals_synonyms[, grep(NewCol, colnames(
    Mammals_synonyms))], 1, paste, collapse = '_') # Add a new column representing the binomial of the pair.
  Counter <- Counter + 1 }; rm(NewCol, Counter) # Update the counter.
Mammals_synonyms <- Mammals_synonyms[, grep('Binomial', colnames(Mammals_synonyms))[-1]] # Subset.  
  
# Generate and read in a dataset containing each Africa-only species' select bioclimatic conditions from WorldClim 
# (this code block was run once and commented out, with its outputs saved, to save time and space):

  # # Read in the previous version of this dataset made in previous versions of this script file:
  # 
  # Mammals_bioclim_orig <- read.csv(list.files()[grep('Worldclim_V1', list.files())])
  # 
  # # Upload the libraries used for reading in and analyzing spatial files:
  # 
  # library(sf) # For using the "read_sf" and "st_as_sf" (if necessary) functions.
  # library(exactextractr) # For using the "exact_extract" function below.
  # 
  # # Using each raster file containing a select bioclimatic variable, determine the realizations of that variable for
  # # each Africa-only species from each's polygon range map (using their IUCN names, since the polygons are from IUCN):
  # 
  # BioClim_res <- 50 # Resolution, in kilometers, at which to collect bioclimatic data.
  # BioClim_varnums <- c(1, 4, 12, 15) # Variable numbers, with "c(1, 12)" for means and "c(4, 15)" for seasonalities.
  # BioClim_varnames <- paste(c('MeanTemp', 'TempSeas', 'MeanPrecip', 'PrecipSeas'), '_', as.character(BioClim_res),
  #   'Km', sep = '') # All variable names, with their resolutions included in their names.
  # names(BioClim_varnames) <- Bioclim_varnums # Associate all variable numbers with their corresponding names.
  # BioClim_vars <- sapply(BioClim_varnums, function(Num) { raster(list.files()[grep(paste('bio_', Num, '.tif',
  #   sep = ''), list.files())]) }) # Read in the variables of interest, based on the numbers in "BioClim_varnums".
  # 
  # if (!(BioClim_res %in% c(1, 5))) { # If the climate data resolution is not set to be either 1 or 5 kilometers...
  #   BioClim_vars <- lapply(BioClim_vars, function(Ras) { # For each raster pertaining to a bioclimatic variable...
  #     Ras <- aggregate(Ras, BioClim_res / 5) }) } # Aggregate the 5-km raster to the resolution of interest.
  # 
  # Mammals_bioclim <- data.frame(matrix(ncol = length(BioClim_varnums) + 1, nrow = length(Mammals_names), dimnames =
  #   list(NULL, c('Binomial', BioClim_varnames[as.character(BioClim_varnums)])))) # Create a dataframe to hold values.
  # Mammals_bioclim$Binomial <- Mammals_synonyms$Binomial5 # Update dataframe's first column to be species' IUCN names.
  # 
  # for (Var in 1:length(BioClim_vars)) { # For each bioclimatic variable of interest...
  # 
  #   # Reference the raster file for the bioclimatic variable (its CRS is appropriate for downstream use already):
  # 
  #   Tmp <- BioClim_vars[[Var]]
  # 
  #   # Reference the column name of "Mammals_bioclim" that refers to the bioclimatic variable:
  # 
  #   Col <- unname(BioClim_varnames[as.character(BioClim_varnums[Var])])
  # 
  #   for (Species in 1:nrow(Mammals_bioclim)) { # For each mammal species...
  # 
  #     # If the species' name is in "Mammals_bioclim_orig", and if "Col" is a column name in that dataset, such that
  #     # the species' data for the variable of interest has already been collected...
  # 
  #     if (Mammals_bioclim$Binomial[Species] %in% Mammals_bioclim_orig$Binomial & # Condition 1.
  #       Col %in% colnames(Mammals_bioclim_orig)) { # Condition 2.
  #       Mammals_bioclim[Species, Col] <- Mammals_bioclim_orig[match(Mammals_bioclim$Binomial[Species],
  #         Mammals_bioclim_orig$Binomial), Col] } # Record the data that has previously been collected.
  # 
  #     # Or else, if the species' name has an associated shapefile/polygon in "MammalShapefiles" (a folder downloaded
  #     # from IUCN and removed thereafter to save space)...
  # 
  #     else if (Mammals_bioclim$Binomial[Species] %in% unique(gsub('\\.(.*)', '', list.files('./MammalShapefiles')))) {
  # 
  #       # Read in the shapefile ("sf") of the species and ensure that its CRS matches that of the raster file:
  # 
  #       Tmp2 <- read_sf(paste('MammalShapefiles/', Mammals_bioclim$Binomial[Species], '.shp', sep = ''))
  #       Tmp2 <- st_transform(Tmp2, crs = Tmp@crs)
  # 
  #       # Extract the mean raster file values across the cells of the raster that overlap with the species' shapefile.
  #       # Record those values in "Mammals_bioclim":
  # 
  #       Mammals_bioclim[Species, Col] <- exact_extract(Tmp, Tmp2, 'mean') }
  # 
  #     else { next } # If the species does not have a shapefile, then skip it.
  # 
  #     # Print a record of this loop's progress:
  # 
  #     print(paste('BioClim variable', as.character(Var), 'species', as.character(Species), 'complete', sep = ' ')) }}
  # 
  # rm(Var, Species, Tmp, Col, Tmp2)
  # 
  # # Save "Mammals_bioclim" as a CSV file. If the file already exists, due to the addition of previous data to it, add
  # # the data in "Mammals_bioclim" to that CSV. Then upload and process the updated CSV file for future use:
  # 
  # if ('AfricaMammalsBioclim_FromWorldclim_V2.csv' %in% list.files()) { # If the file already exists...
  #   Mammals_bioclim_prev <- read.csv(list.files()[grepl('Worldclim_V2', list.files())]) # Read it in.
  #   Mammals_bioclim <- cbind(Mammals_bioclim_prev, Mammals_bioclim[, -1]) } # Add new data to it.
  # write.csv(Mammals_bioclim, file = 'AfricaMammalsBioclim_FromWorldclim_V2.csv', row.names = F) # Save updated file.
  Mammals_bioclim <- read.csv(list.files()[grepl('Worldclim_V2', list.files())]) # Read in the updated file.

# Similarly, generate and read in a dataset containing other relevant variables of each Africa-only species from various
# data sources (this code block was run once per var and commented out, with its outputs saved, to save time and space):
  
  # # Read in the previous version of this dataset made in previous versions of this script file:
  # 
  # Mammals_other_orig <- read.csv(list.files()[grep('FromMisc_V1', list.files())])
  # 
  # # Upload the libraries used for reading in and analyzing spatial files:
  # 
  # library(sf) # For using the "read_sf" and "st_as_sf" (if necessary) functions.
  # library(exactextractr) # For using the "exact_extract" function below.
  # library(rgeos) # For using the "gArea" and "gIntersection" functions below.
  # library(geosphere) # For using the "centroid" function below.
  # 
  # # Using the raster/polygon file of the variable of interest, determine the realizations of that variable for
  # # each Africa-only species from each's polygon range map (using their IUCN names, since the polygons are from IUCN):
  # 
  #   # Uncomment if addressing the geographic range of each species (source = IUCN):
  # 
  #   # Other_var_name <- 'GeogRangeIUCN' # Variable name.
  #   # rm(Other_var) # This does not involve an "Other_var".
  # 
  #   # Uncomment if addressing habitat heterogeneity, represented by the mean coefficient of variation of EVI across
  #   # each species' range (source = EarthEnv):
  # 
  #   # Other_var_name <- 'HetHab_VarEVI_1Km' # Variable name at 1-km resolution.
  #   # Other_var <- raster('cv_01_05_1km_uint16.tif') # Raster for the variable at 1-km resolution, OR...
  #   Other_var_name <- 'HetHab_VarEVI_5Km_Updated' # Variable name at 5-km resolution.
  #   Other_var <- aggregate(raster('cv_01_05_1km_uint16.tif'), 5) # Raster for the variable at 5-km resolution, OR...
  #   # Other_var_name <- 'HetHab_VarEVI_50Km' # Variable name at 50-km resolution.
  #   # Other_var <- aggregate(raster('cv_01_05_25km_uint16.tif'), 2) # Raster for the variable at 50-km resolution.
  # 
  #   # Uncomment if addressing topographic heterogeneity, represented by the mean ruggedness of the terrain (e.g. the
  #   # mean degree of variation in terrain slope) across each species' range (source = EarthEnv):
  # 
  #   # Other_var_name <- 'HetTopo_Rugged_1Km' # Variable name at 1-km resolution.
  #   # Other_var <- raster('tri_1KMmn_GMTEDmd.tif') # Raster for the variable at 1-km resolution, OR...
  #   # Other_var_name <- 'HetTopo_Rugged_5Km' # Variable name at 5-km resolution.
  #   # Other_var <- raster('tri_5KMmn_GMTEDmd.tif') # Raster for the variable at 5-km resolution, OR...
  #   # Other_var_name <- 'HetTopo_Rugged_50Km' # Variable name at 50-km resolution.
  #   # Other_var <- raster('tri_50KMmn_GMTEDmd.tif') # Raster for the variable at 50-km resolution.
  # 
  #   # Uncomment if addressing the mean human population density across each species' range (source = CIESIN):
  # 
  #   # Other_var_name <- 'HumPopDen_1Km' # Variable name at 1-km resolution.
  #   # Other_var <- raster('gpw_v4_population_density_rev11_2020_30_sec.tif') # Raster for the variable, OR...
  #   # Other_var_name <- 'HumPopDen_5Km' # Variable name at 5-km resolution.
  #   # Other_var <- raster('gpw_v4_population_density_rev11_2020_2pt5_min.tif') # Raster for the variable, OR...
  #   # Other_var_name <- 'HumPopDen_50Km' # Variable name at 50-km resolution.
  #   # Other_var <- aggregate(raster('gpw_v4_population_density_rev11_2020_2pt5_min.tif'), 10) # Raster for the var.
  # 
  #   # Uncomment if addressing the mid-latitude of the centroid of each species' range (source = IUCN):
  # 
  #   # Other_var_name <- 'MidLat' # Variable name.
  #   # rm(Other_var) # This does not involve an "Other_var".
  # 
  #   # Uncomment if addressing the percentage of each species' range that overlaps with PAs (source = WDPA):
  # 
  #   # Other_var_name <- 'PctPAOverlap' # Variable name.
  #   # Other_var <- list() # To hold the three separate shapefiles that the WDPA uses to store all of its PA info.
  #   # for (File in 1:3) { # For each of the three shapefiles provided by the WDPA...
  #   #   Other_var[[File]] <- shapefile(paste('./WDPA_WDOECM_wdpa_shp', as.character(File-1),
  #   #     '/WDPA_WDOECM_wdpa_shp-polygons.shp', sep = '')) }; rm(File) # Read it in.
  #   # Other_var <- lapply(Other_var, function(Item) { # For each of the three files...
  #   #   Item <- Item[Item@data$IUCN_CAT %in% c('Ia', 'Ib', 'II', 'III', 'IV'), ] # Subset to strictest PA categories.
  #   #   Item <- crop(Item, Africa) }) # Crop to just the extent of Africa.
  #   # Other_var <- bind(Other_var[[1]], Other_var[[2]], Other_var[[3]]) # Combine all files into a single layer.
  # 
  #   # Uncomment if addressing the primary biome/ecoregion in which each species occurs (source = Olson et al., 2001):
  # 
  #   # Other_var_name <- 'PrimaryBiome' # Variable name.
  #   # Other_var <- shapefile('./official/wwf_terr_ecos.shp') # Shapefile for the variable.
  # 
  # if (exists('Other_var')) { # If "Other_var" is an existent variable...
  #   if (grepl('Polygon', class(Other_var))) { # If the variable is represented by a polygon...
  #     Other_var <- crop(Other_var, Africa) }} # Clip the polygon to just the extent of Africa, if not done already.
  # 
  # Mammals_other <- data.frame(matrix(ncol = 2, nrow = length(Mammals_names), dimnames = list(NULL, c(
  #   'Binomial', Other_var_name)))) # Create a dataframe to hold values of interest.
  # Mammals_other$Binomial <- Mammals_synonyms$Binomial5 # Update dataframe's first column to be species' IUCN names.
  # 
  # for (Species in 1:nrow(Mammals_other)) { # For each mammal species...
  # 
  #   # If the species' name is in "Mammals_other_orig", and if "Other_var_name" is a column name in that dataset,
  #   # such that the species' data for the variable of interest has already been collected...
  # 
  #   if (Mammals_other$Binomial[Species] %in% Mammals_other_orig$Binomial & # Condition 1.
  #     Other_var_name %in% colnames(Mammals_other_orig)) { # Condition 2.
  #     Mammals_other[Species, Other_var_name] <- Mammals_other_orig[match(Mammals_other$Binomial[Species],
  #       Mammals_other_orig$Binomial), Other_var_name] } # Record the data that has previously been collected.
  # 
  #   # Or else, if the species' name has an associated shapefile/polygon in "MammalShapefiles" (a folder downloaded
  #   # from IUCN and removed thereafter to save space)...
  # 
  #   else if (Mammals_other$Binomial[Species] %in% unique(gsub('\\.(.*)', '', list.files('./MammalShapefiles')))) {
  # 
  #     # Read in the shapefile ("sf") of the species and ensure that its CRS matches that of the raster/polygon
  #     # file representing the variable of interest:
  # 
  #     Tmp <- read_sf(paste('MammalShapefiles/', Mammals_other$Binomial[Species], '.shp', sep = ''))
  #     Tmp <- st_transform(Tmp, crs = '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')
  # 
  #     # Extract information addressing the variable of interest within the species' range:
  # 
  #     if (exists('Other_var')) { # If "Other_var" is an existent variable...
  #       if (grepl('Raster', class(Other_var))) { # If the variable is represented by a raster...
  #         Mammals_other[Species, Other_var_name] <- exact_extract(Other_var, Tmp, 'mean') }} # Extract
  #           # the mean raster file values across the cells of the raster that overlap with the species' shapefile.
  # 
  #     if (Other_var_name == 'GeogRangeIUCN') { # If the variable is "GeogRangeIUCN"...
  #       Mammals_other[Species, Other_var_name] <- area(Tmp) / 1000000 } # Divide because in lat/long projection.
  # 
  #     if (Other_var_name == 'MidLat') { # If the variable is "MidLat"...
  #       Mammals_other[Species, Other_var_name] <- centroid(Tmp)[2] }
  # 
  #     if (Other_var_name == 'PctPAOverlap') { # If the variable is "PctPAOverlap"...
  #       tryCatch({ # Calculate percent PA overlap, accounting for a situation in which there is no overlap...
  #       if (is.null(raster::intersect(Other_var, Tmp))) { Mammals_other[Species, Other_var_name] <- 0 } else {
  #         Mammals_other[Species, Other_var_name] <- sum(area(raster::intersect(Other_var, Tmp))) / area(Tmp) }},
  #         error = function(E) { # For Loxodonta africana, do this instead...
  #           Mammals_other[Species, Other_var_name] <- gArea(gIntersection(Tmp, Other_var)) / gArea(Tmp) }) }
  # 
  #     if (Other_var_name == 'PrimaryBiome') { # If the variable is "PrimaryBiome"...
  #       Tmp2 <- unique(Other_var@data$BIOME)[unique(Other_var@data$BIOME) <= 14] # Unique biome IDs.
  #       tryCatch({ # Determine which biome ID has the most overlap with the species' range...
  #       Mammals_other[Species, Other_var_name] <- Tmp2[which.max(sapply(Tmp2, function(Biome) {
  #         if (is.null(raster::intersect(Other_var[Other_var@data$BIOME == Biome, ], Tmp))) { return(0) } else {
  #           return(sum(area(raster::intersect(Other_var[Other_var@data$BIOME == Biome, ], Tmp)))) }}))] },
  #           error = function(E) { # For Loxodonta africana, do this instead...
  #             Mammals_other[Species, Other_var_name] <- Tmp2[which.max(sapply(Tmp2, function(Biome) {
  #               if (is.null(gIntersection(Other_var[Other_var@data$BIOME == Biome, ], Tmp))) { return(0) } else {
  #                 return(gArea(gIntersection(Other_var[Other_var@data$BIOME == Biome, ], Tmp))) }}))] }) }}
  # 
  #   else { next } # If the species does not have a shapefile, then skip it.
  # 
  #   # Print a record of this loop's progress:
  # 
  #   print(paste('Species', as.character(Species), 'complete', sep = ' ')) }; rm(Species, Tmp, Tmp2)
  # 
  # # Save "Mammals_other" as a CSV file. If the file already exists, due to the addition of previous data to it, add
  # # the data in "Mammals_other" to that CSV. Then upload the updated CSV file for future use:
  # 
  # if ('AfricaMammalsOther_FromMisc_V2.csv' %in% list.files()) { # If the file already exists...
  #   Mammals_other_prev <- read.csv(list.files()[grepl('FromMisc_V2', list.files())]) # Read it in.
  #   Mammals_other <- cbind(Mammals_other_prev, Mammals_other[, -1]) # Add new data to it.
  #   colnames(Mammals_other)[ncol(Mammals_other)] <- Other_var_name } # Update new data's column name.
  # write.csv(Mammals_other, file = 'AfricaMammalsOther_FromMisc_V2.csv', row.names = FALSE) # Save updated file.
  Mammals_other <- read.csv(list.files()[grepl('FromMisc_V2', list.files())]) # Read in updated file.

# Generate maps of select features represented in "Mammals_other", as they appear across all of Africa (this code block
# is commented out to save the space of holding map files in memory - download files and produce maps as needed):
  
library(DescTools) # For various functions below, including the "BoxCox" function.

# # Map <- shapefile('./official/wwf_terr_ecos.shp'); Map_var <- 'Biomes' # If plotting the biomes/ecoregions in Africa.
# Map <- raster('cv_01_05_1km_uint16.tif'); Map_var <- 'Hab_Xkm' # If plotting habitat heterogeneity at X-km res.
# # Map <- raster('tri_1KMmn_GMTEDmd.tif'); Map_var <- 'Topo_Xkm' # If plotting topographic heterogeneity at X-km res.
# 
# if (grepl('Raster', class(Map))) { # If the map is a raster...
# 
#   # Process the map and determine its settings for plotting it:
# 
#   Map <- crop(Map, Africa) # If the CRS of the map is the same as that of "Africa", crop it to "Africa"'s extent.
#   Map <- mask(Map, Africa) # Clip the cropped map such that it encompasses only the shape of "Africa".
#   # Map <- raster::aggregate(Map, 50) # Potentially aggregate the map to a coarser resolution.
# 
#   Hab <- grepl('Hab', Map_var) # Indicate if the map refers to habitat heterogeneity or not.
#   Map <- Map * ifelse(Hab, 0.0001, 1) # For habitat heterogeneity, multiply by 0.0001 based on EarthEnv metadata.
# 
#   Transform <- '(ln)' # OR "(Box-Cox)", depending upon which transformation of the map's values will be performed.
#   if (Transform == '(Box-Cox)') { values(Map)[!is.na(values(Map))] <- BoxCox(values(Map)[!is.na(values(Map))],
#     lambda = BoxCoxLambda(values(Map)[!is.na(values(Map))], method = 'guerrero', lower = -5, upper = 5)) }
#   if (Transform == '(ln)') { Map <- log(Map); values(Map)[is.infinite(values(Map))] <- NA } # Remove infinite values.
# 
#   Col_no <- 255 # Choose the number of colors that will be used in plotting the map.
#   Cols <- colorRampPalette(c('peachpuff', ifelse(Hab, 'cadetblue4', 'brown4')))(Col_no) # Choose color scheme.
#   Breaks <- unname(quantile(values(Map), seq(0, 1, length.out = Col_no + 1), na.rm = TRUE)) # Breaks for map colors.
#   for (Val in which(duplicated(Breaks))) { Breaks[Val] <- Breaks[Val-1] + 0.000001 } # Prevent duplicate breaks.
# 
#   # Plot the map, using the settings above:
# 
#   par(mar = c(1.6,1.3,1,1)) # Adjust plot margins.
#   plot(Map, axes = TRUE, box = TRUE, cex.axis = 3.5, mgp = c(3,2,0), # Plot map. Maybe remove axes and box.
#     col = Cols, breaks = Breaks, # Set plot colors, and scale them based on distribution of "Map" values.
#     # legend.args = list(text = ifelse(Hab, paste('Coef Var EVI', Transform), paste('TRI', Transform)),
#     #   side = 4, line = 5, font = 2, cex = 2), # Format legend title.
#     axis.args = list(cex.axis = 3.5, at = c(minValue(Map), maxValue(Map)), labels=c('Few\nHabitats','Many\nHabitats')))
#   }; rm(Hab, Transform, Col_no, Breaks, Val) # "c('Few\nHabitats','Many\nHabitats')" or "c('Smooth', 'Rugged')".
# 
# pdf('./Map_Biomes.pdf') # Initiate a PDF file to hold the map made below, which is faster than R rendering the map.
# if (grepl('Polygon', class(Map))) { # If the map is a polygon...
#   Map <- crop(Map, Africa) # If the CRS of the map is the same as that of "Africa", crop it to "Africa"'s extent.
#   if (Map_var == 'Biomes') { # If a map addressing the biomes in Africa is being made...
#     Map@data$BIOME <- ifelse(Map@data$BIOME %in% c(1:6, 12, 14), 'Forest', ifelse(Map@data$BIOME %in% 7:11,
#       'Savanna', 'Desert')) # Convert the map's "BIOME" column into fewer discrete categories.
#     Map@data$COLORS <- ifelse(Map@data$BIOME == 'Forest', 'darkolivegreen4', ifelse(Map@data$BIOME == 'Savanna',
#       'orange', 'peachpuff')) # Make colors for the map based upon the categories in the "BIOME" column.
#     plot(Map, col = Map@data$COLORS, border = Map@data$COLORS) # Plot the map, using the specified colors.
#     legend('topright', legend = unique(Map@data$BIOME), fill = unique(Map@data$COLORS), bty = 'n') }}
# dev.off() # Hand plotting control back to R.

# Analyze habitat vs. topographic heterogeneity across Africa to assess if the two are correlated (this code block is
# commented out to save the space of holding map files in memory - download files and perform analyses as needed):  

# Map_hab <- raster('cv_01_05_1km_uint16.tif') # Habitat heterogeneity map at 1-km resolution.
# Map_topo <- raster('tri_1KMmn_GMTEDmd.tif') # Topographic heterogeneity map at 1-km resolution.
# 
# for (Het_map in c('Map_hab', 'Map_topo')) { # For each map...
#   Map <- eval(parse(text = Het_map)) # Create a variable, "Map", the is equal to the current map.
#   Map <- crop(Map, Africa) # Crop the map to the extent of "Africa".
#   Map <- mask(Map, Africa) # Mask the cropped map to the shape of "Africa".
#   Map <- Map * ifelse(grepl('hab', Het_map), 0.0001, 1) # Multiply the "hab" map by 0.0001 based on EarthEnv metadata.
#   Map <- log(Map); values(Map)[is.infinite(values(Map))] <- NA # Log-transform map values and remove resultant Inf.
#   # Map <- raster::aggregate(Map, 50) # Potentially aggregate map to a coarser resolution.
#   assign(Het_map, Map) }; rm(Het_map, Map) # Update the map to be equal to "Map", now that "Map" is processed.
# 
# Map_valindices <- which(!is.na(values(Map_hab)) & !is.na(values(Map_topo))) # Map pixel indices with paired values.
# 
# library(spatialRF) # For assessing spatial autocorrelation and using the "moran" function.
# 
# Map_valindices_sampler <- function(Num) { # Function to sample spatially-unautocorrelated pixels from "Map_valindices".
#   Samp <- sample(Map_valindices, Num) # Randomly sample "Num" raster pixel indices from "Map_valindices".
#   Coordinates <- xyFromCell(Map_hab, cell = Samp) # Get the lon/lat coordinates of the raster pixels in the sample.
#   Distance_matrix <- pointDistance(Coordinates, lonlat = TRUE) # Produce matrix of distances between the sample pixels.
#   Distance_matrix[upper.tri(Distance_matrix)] <- t(Distance_matrix)[upper.tri(Distance_matrix)] # Make matrix full.
#   Moran_hab <- moran(x = Map_hab[Samp], distance.matrix = Distance_matrix, verbose = FALSE) # Moran's I for hab.
#   Moran_topo <- moran(x = Map_topo[Samp], distance.matrix = Distance_matrix, verbose = FALSE) # Moran's I for topo.
#   if (abs(Moran_hab$test$moran.i) > 0.1 | abs(Moran_topo$test$moran.i) > 0.1) { # If sample has autocorrelation...
#     print('Sample exhibits spatial autocorrelation. Trying again with a new sample') # Provide a message.
#     Map_valindices_sampler(Num) } else { # Recursively try again with a new sample. Otherwise...
#     return(Samp) }} # Return the sample with minimal spatial autocorrelation.
# 
# Map_valindices_samp <- Map_valindices_sampler(1000) # Apply the above function to obtain a proper sample of pixels.
# 
# library(ggplot2) # For making plot visualizations.
# 
# Map_df <- data.frame(Hab = values(Map_hab)[Map_valindices_samp], Topo = values(Map_topo)[Map_valindices_samp],
#   Lat = abs(xyFromCell(Map_hab, Map_valindices_samp)[,2])) # Data frame with variables for the plots below.
# 
# ggplot(Map_df, aes(Hab, Topo)) + # Select data for baseline plot of habitat vs. topographic heterogeneity values.
#   geom_point(size = 3) + geom_smooth(method = 'lm', col = 'deepskyblue3', lwd = 3) + # Scatterplot w/ regression line.
#   labs(x = 'Habitat Heterogeneity', y = 'Topographic Heterogeneity') + # Set axis labels.
#   coord_cartesian(xlim = c(-6, -1), ylim = c(-1, 5)) + # Set axis range limits.
#   scale_x_continuous(labels = function(Lab) { sprintf('%.1f', Lab) }) + # Set X axis tick labels to have one decimal.
#   scale_y_continuous(labels = function(Lab) { sprintf('%.1f', Lab) }) + # Set Y axis tick labels to have one decimal.
#   theme_classic() + # Make the plot have a classic-looking theme.
#   theme(axis.text = element_text(size = 40), axis.text.x = element_text(hjust = 0.8), axis.title = element_text(
#     size = 40), axis.title.x = element_text(margin = margin(t = 20)), axis.title.y = element_text(margin = margin(
#     r = 25))) # Format axes.
# 
# Map_df$Resids <- lm(Map_df$Topo ~ Map_df$Hab)$residuals # Add column of residuals of plot above to "Map_df".
# 
# ggplot(Map_df, aes(Lat, Resids)) + # Select data for plot of residuals vs. abs latitude to assess latitude's affect.
#   geom_point(size = 3) + geom_smooth(method = 'loess', col = 'deepskyblue3', lwd = 3) + # Scatterplot with LOESS curve.
#   geom_hline(yintercept = 0, col = 'darkgoldenrod3', lty = 1, size = 2) + # Add horizontal line at Y = 0.
#   labs(x = 'Latitude (abs Degrees)', y = 'Residuals') + # Set axis labels.
#   coord_cartesian(xlim = c(0, 38), ylim = c(-2, 3.5)) + # Set axis range limits.
#   scale_x_continuous(labels = function(Lab) { sprintf('%.1f', Lab) }) + # Set X axis tick labels to have one decimal.
#   scale_y_continuous(labels = function(Lab) { sprintf('%.1f', Lab) }) + # Set Y axis tick labels to have one decimal.
#   theme_classic() + # Make the plot have a classic-looking theme.
#   theme(axis.text = element_text(size = 40), axis.title = element_text(size = 40), axis.title.x = element_text(
#     margin = margin(t = 20))) # Format axes.
# 
# HeteroCor <- cor(Map_df$Hab, Map_df$Topo, method = 'spearman') # Hab vs. topo correlation.
# ResidsCor <- cor(Map_df$Lat, Map_df$Resids, method = 'spearman') # Resids vs. lat correlation.
# 
# HeteroCor_iter <- c() # Vector to hold multiple random iterations of "HeteroCor".
# for (Iter in 1:100) { # Repeat the measurement of "HeteroCor" many times to generate a distribution of it...
#   Samp <- Map_valindices_sampler(1000) # Obtain a proper sample of raster pixels, then run the calculation on them...
#   HeteroCor_iter[Iter] <- cor(values(Map_hab)[Samp], values(Map_topo)[Samp], method = 'spearman') }; rm(Iter, Samp)
# 
# par(mar = c(7.75,7.75,1,1)) # Adjust plot margin settings.
# hist(HeteroCor_iter, main = '', xlab = '', ylab = '', cex.axis = 2.75, mgp = c(3,1.5,0), col = 'deepskyblue3',
#   breaks = 20) # Historgram of the values in "Heterocor_iter". 
# mtext('Correlation Value', side = 1, line = 5, font = 1, cex = 4) # X-axis label.
# mtext('Frequency', side = 2, line = 5, font = 1, cex = 4) # Y-axis label.

# Determine the resolution at which the below analyses will be done, and based off of that, select which columns from
# "Mammals_bioclim" and "Mammals_other" should be kept. Also make final edits to both dataframes:

Analyses_res <- 1 # Resolution at which the analyses to follow will be done (1, 5, or 50 kilometers).

for (Dataset in c('Mammals_bioclim', 'Mammals_other')) { # For each dataset...
  Dat <- eval(parse(text = Dataset)) # Create a variable, "Dat", the is equal to the dataset.
  Dat <- Dat[, sort(c(setdiff(1:ncol(Dat), grep('Km', colnames(Dat))), grep(paste(as.character(Analyses_res), 'Km',
    sep = ''), colnames(Dat))))] # Keep columns that either don't have a 'Km' resolution or that have the desired one.
  colnames(Dat)[grep('Km', colnames(Dat))] <- gsub('_[^_]+$', '', colnames(Dat)[grep('Km', colnames(Dat))]) # Nix 'Km'.
  assign(Dataset, Dat) }; rm(Dataset, Dat) # Update the dataset to be equal to "Dat", now that "Dat" is processed.

Mammals_other$PrimaryBiome <- as.factor(Mammals_other$PrimaryBiome) # Change to factor.
Mammals_other$HetHab_VarEVI <- Mammals_other$HetHab_VarEVI * 0.0001 # Make this change based on EarthEnv metadata.

for (Name in c('Hab', 'Topo')) { # For each of these phrases, which appear in column names...
  colnames(Mammals_other)[grep(Name, colnames(Mammals_other))] <- gsub('_.*', '', colnames(Mammals_other)[grep(Name,
    colnames(Mammals_other))]) }; rm(Name) # Remove the part of the respective column name that refers to a metric.

# Read in datasets with features of Africa-only species from Phylacine and the literature:

Mammals_paleo <- read.table(list.files()[grepl('Paleo', list.files())], header = TRUE, sep = '\t')
Mammals_ranges <- read.csv(list.files()[grepl('Ranges_FromPhylacine', list.files())])
Mammals_speeds <- read.csv(list.files()[grepl('Locomotion', list.files())])
Mammals_traits <- read.csv(list.files()[grepl('Traits_FromPhylacine', list.files())])

# Update the binomial columns of select datasets to include underscores and to have the proper column name:

colnames(Mammals_speeds)[1] <- 'Binomial'; Mammals_speeds$Binomial <- gsub(' ', '_', Mammals_speeds$Binomial)
colnames(Mammals_paleo)[colnames(Mammals_paleo) == 'Species'] <- 'Binomial'

# Update the names in the various datasets to be consistent with those in Phylacine, using the synonyms provided by
# Phylacine to do so:

Mammals_bioclim$Binomial <- as.character(Mammals_bioclim$Binomial) # Convert to character.
Mammals_other$Binomial <- as.character(Mammals_other$Binomial) # Convert to character.
Mammals_paleo$Binomial <- as.character(Mammals_paleo$Binomial) # Convert to character.
Mammals_synonyms$Binomial1 <- as.character(Mammals_synonyms$Binomial1) # Convert to character.

for (Dataset in ls()[grepl('bioclim$|other$|paleo|speeds', ls())]) { # For each dataset...
  Dat <- eval(parse(text = Dataset)) # Create a variable, "Dat", the is equal to the dataset.
  Dat[, grep('Binomial', colnames(Dat), ignore.case = TRUE)] <- unname( # Modify the "Binomial" column of the dataset.
    sapply(Dat[, grep('Binomial', colnames(Dat), ignore.case = TRUE)], function(Name) { # For each name in the col...
      Tmp <- which(Mammals_synonyms == Name, arr.ind = TRUE) # Which row of "Mammals_synonyms" has the name.
      if (nrow(Tmp) > 0) { # If the name in the column is found in "Mammals_synonyms"...
        return(Mammals_synonyms$Binomial1[Tmp[1,1]]) } else { # Make the name = that in the first col of "synonyms".
          return(Name) }})) # Or else, if the name is not found in "synonyms", keep the original name.
  assign(Dataset, Dat) }; rm(Dataset, Dat) # Update the dataset to be equal to "Dat" with the updated names.

# Update "Mammals_speeds" by subsetting it down to just Africa-only species and finding the max of their speed values:

Mammals_speeds <- Mammals_speeds[Mammals_speeds$Taxon == 'mammal', ]
Mammals_speeds <- aggregate(Mammals_speeds$Speed, by = list(Mammals_speeds$Binomial), FUN = max)
colnames(Mammals_speeds) <- c('Binomial', 'MaxSpeed')
Mammals_speeds_orig <- Mammals_speeds # Keep an original copy of this data before subsetting for later.
Mammals_speeds <- Mammals_speeds[Mammals_speeds$Binomial %in% Mammals_names, ]

# Assign dietary categories to all Africa-only species, based on their diet percentage data in Phylacine:

Mammals_traits$TrophicLevel <- ifelse(Mammals_traits$Diet.Plant == 100, 'Herbivore', ifelse(
  Mammals_traits$Diet.Plant == 0, 'Carnivore', 'Omnivore'))

# Reduce the Phylacine dataframes in size such that they only contain Africa-only species:

Reduced_dataframes <- lapply(list(Mammals_traits, Mammals_ranges), function(Dat) { 
  Dat <- Dat[which(Dat[, grep('Binomial', colnames(Dat), ignore.case = TRUE)] %in% Mammals_names), ] })
Mammals_traits <- Reduced_dataframes[[1]]; Mammals_ranges <- Reduced_dataframes[[2]]; rm(Reduced_dataframes)

# Create a dataframe that contains features of interest for Africa-only species, compiled from different dataframes:

Mammals_final <- cbind(Mammals_traits[, grep('Binomial|Order|Family|Genus|Species', colnames(Mammals_traits))],
  Mammals_traits$Mass.g, Mammals_traits$TrophicLevel, Mammals_bioclim[, -grep('Binomial', colnames(Mammals_bioclim))],
  Mammals_other[, -grep('Binomial|MidLat', colnames(Mammals_other))], Mammals_ranges$Number.Cells.Current.Range)
colnames(Mammals_final) <- c('Binomial', 'Order', 'Family', 'Genus', 'Species', 'BodyMass', 'TrophicLevel',
  colnames(Mammals_bioclim)[-grep('Binomial', colnames(Mammals_bioclim))], colnames(Mammals_other)[-grep(
  'Binomial|MidLat', colnames(Mammals_other))], 'GeogRangePhy') # Give the selected columns clean names.

# Add the data of interest from "Mammals_paleo" into "Mammals_final":

Mammals_final$Binomial <- as.character(Mammals_final$Binomial) # Convert to character.
Mammals_final$SuitChange <- Mammals_paleo$DeltaSuit[match(Mammals_final$Binomial, Mammals_paleo$Binomial)]

# Perform unit conversions on specific columns in "Mammals_final":

Mammals_final$BodyMass <- Mammals_final$BodyMass / 1000 # Convert mass from g to kg.
Mammals_final$GeogRangePhy <- Mammals_final$GeogRangePhy * 96.5 * 96.5 # Convert range sizes into km^2.

# Determine (visually and statistically) if "GeogRangePhy" and "GeogRangeIUCN" in "Mammals_final" are correlated:

Par_orig <- par() # Keep a record of initial plot settings.
par(mar = c(7,7,3,3)) # Adjust plot margins.

plot(Mammals_final$GeogRangeIUCN, Mammals_final$GeogRangePhy, xlab = '', ylab = '', cex.axis = 2, cex = 2, pch = 19)
mtext('Range from IUCN (Sq Km)', side = 1, line = 4, font = 2, cex = 2) # Format the X axis label.
mtext('Range from Phylacine (Sq Km)', side = 2, line = 4, font = 2, cex = 2) # Format the Y axis label.
abline(lm(GeogRangePhy ~ GeogRangeIUCN, data = Mammals_final), col = 'red', lwd = 5) # Add regression line.

GeogRangeCor <- cor(Mammals_final$GeogRangeIUCN, Mammals_final$GeogRangePhy, use = 'complete.obs')

# Remove one of the two "GeogRange" columns, and name the remaining one "GeogRange":

Mammals_final$GeogRange <- Mammals_final$GeogRangeIUCN # Choose which column will be "GeogRange" column.
Mammals_final$GeogRange[which(is.na(Mammals_final$GeogRange))] <- Mammals_final$GeogRangePhy[which(is.na(
  Mammals_final$GeogRange))] # Replace missing values of range with values from other column.
Mammals_final <- Mammals_final[, -grep('IUCN|Phy', colnames(Mammals_final))] # Remove old columns.

# Determine locomotor characteristics of all Africa-only mammals remaining in "Mammals_final" in various ways from
# "Mammals_speeds" and from the literature:

  # Foot position --> categorize each mammal species as digitigrade, plantigrade, unguligrade, or winged:

  library(tibble) # For using the "add_column" function below.
  
  if (!('Locom_FootPos' %in% colnames(Mammals_final))) { # If the "Locom_FootPos" column does not yet exist...
    Mammals_final <- add_column(Mammals_final, Locom_FootPos = NA, .before = 'MeanTemp') } # Initiate and place it.
  Mammals_final$Locom_FootPos[Mammals_final$Family %in% c('Anomaluridae', 'Canidae', 'Dipodidae', 'Erinaceidae',
    'Felidae', 'Hippopotamidae', 'Herpestidae', 'Hyaenidae', 'Leporidae', 'Macroscelididae', 'Nesomyidae',
    'Orycteropodidae', 'Procaviidae', 'Rhinocerotidae', 'Thryonomyidae')] <- 'D' # Digitigrade taxa.
  Mammals_final$Locom_FootPos[Mammals_final$Family %in% c('Bathyergidae', 'Cercopithecidae', 'Chrysochloridae',
    'Cricetidae', 'Ctenodactylidae', 'Elephantidae', 'Galagidae', 'Gliridae', 'Hominidae', 'Hystricidae', 'Leporidae',
    'Lorisidae', 'Manidae', 'Muridae', 'Mustelidae', 'Nandiniidae', 'Pedetidae', 'Petromuridae', 'Procaviidae',
    'Sciuridae', 'Soricidae', 'Spalacidae', 'Tenrecidae', 'Viverridae')] <- 'P' # Plantigrade taxa.
  Mammals_final$Locom_FootPos[Mammals_final$Family %in% c('Bovidae', 'Equidae', 'Giraffidae', 'Suidae',
    'Tragulidae')] <- 'U' # Unguligrade taxa.
  Mammals_final$Locom_FootPos[Mammals_final$Order == 'Chiroptera'] <- 'W' # Winged taxa.
  
  Mammals_final$Locom_FootPos[Mammals_final$Species %in% c('barbarus', 'civetta', 'erythropus')] <- 'D' # Modifications.
  Mammals_final$Locom_FootPos <- as.factor(Mammals_final$Locom_FootPos) # Change to factor.
  
  # Travel speed --> calculate the max speed of each species in kilometers per hour. Do this by first deriving speed
  # values based upon equations, and then replacing those values with known speeds from the literature where possible.
  # Some speed values are random numbers, so perform a sensitivity analysis to ensure that different iterations of
  # random speed values don't greatly impact speed's correlations with species' body masses and range sizes:
  
  if (!('Locom_Speed' %in% colnames(Mammals_final))) { # If the "Locom_Speed" column does not yet exist...
    Mammals_final <- add_column(Mammals_final, Locom_Speed = NA, .after = 'Locom_FootPos') } # Initiate and place it.
  
  Tmp <- which(Mammals_final$Order != 'Chiroptera') # Row indices of all non-volant Africa-only mammals.
  Tmp2 <- setdiff(1:nrow(Mammals_final), Tmp) # Row indices of all volant Africa-only mammals.
  Mammals_final$Locom_Speed[Tmp] <- 25.5 * Mammals_final$BodyMass[Tmp]^0.26 * (1 - exp(-22 *
    Mammals_final$BodyMass[Tmp]^-0.6)) # Equation for deriving speeds in km/hr for non-volant species.
  Mammals_final$Locom_Speed[Tmp2] <- 142.8 * Mammals_final$BodyMass[Tmp2]^0.24 * (1 - exp(-2.4 *
    Mammals_final$BodyMass[Tmp2]^-0.72)) # Equation for deriving speeds in km/hr for volant species.
  
  Locom_rangecor <- c() # Vector to hold correlations of locomotor speed with range size for different random seeds.
  Locom_masscor <- c() # Vector to hold correlations of locomotor speed with body mass for different random seeds.
  Counter <- 1 # Counter to log values in "Locom_rangecor" and "Locom_masscor" in the loop below.
  for (Seed in c(round(runif(20, 1, 2000)), Random_seed)) { # For each of a number of seeds, including "Random_seed"...

    set.seed(Seed) # Set the seed.
    for (G in unique(Mammals_final$Genus)) { # For each unique mammalian genus...
      if (G %in% gsub('_.*', '', Mammals_speeds_orig$Binomial)) { # If the genus is in "Mammals_speeds_orig", set the
        # speeds of its species to random numbers in between the genus' upper and lower speed bounds in that dataset.
        Mammals_final$Locom_Speed[Mammals_final$Genus == G] <- runif(sum(Mammals_final$Genus == G),
          min(Mammals_speeds_orig$MaxSpeed[gsub('_.*', '', Mammals_speeds_orig$Binomial) == G]),
          max(Mammals_speeds_orig$MaxSpeed[gsub('_.*', '', Mammals_speeds_orig$Binomial) == G])) }}; rm(G)
    
    Locom_speed_est <- list( # Initiate list of taxa paired with their bounds for estimated speeds, derived from lit.
      list('Muridae', c(6, 14)), # Miscellaneous rats and mice before specific taxa are addressed next.
      list(c('Acomys', 'Dendromus', 'Graphiurus', 'Hybomys', 'Hylomyscus', 'Lemniscomys', 'Lophuromys', 'Mastomys',
        'Mus', 'Petromyscus', 'Praomys', 'Saccostomus', 'Steatomys'), c(11, 14)), # Mice.
      list(c('Aethomys', 'Arvicanthis', 'Beamys', 'Cricetomys', 'Dasymys', 'Desmomys', 'Grammomys', 'Malacomys',
        'Otomys', 'Pedetes', 'Pelomys', 'Stenocephalemys', 'Thallomys'), c(6, 10)), # Rats and springhares.
      list('Bathyergidae', c(2, 5)), # Mole rats.
      list(c('Cephalophus', 'Neotragus', 'Oryx', 'Philantomba', 'Raphicerus', 'Redunca'), c(60, 70)), # Misc. bovids.
      list(c('Cercopithecidae', 'Hominidae'), c(35, 40)), list('Papio', c(45, 45)), # Monkeys.
      list('Chrysochloridae', c(0.5, 1.2)), # Moles.
      list(c('Emballonuridae', 'Miniopteridae'), c(40, 50)), list(c('Hipposideridae', 'Rhinolophidae'), c(30, 40)),
        list('Nycteridae', c(20, 30)), list('Molossidae', c(61, 96.5)), list('Pteropodidae', c(45, 56)),
        list('Tadarida', c(96.5, 160)), list('Vespertilionidae', c(55, 64)), # Bats.
      list('Galagidae', c(14, 16)), # Bush babies.
      list('Genetta', c(31, 33)), # Genets.
      list(c('Gerbilliscus', 'Gerbillurus', 'Gerbillus', 'Taterillus'), c(8.8, 12.8)), # Gerbils.
      list(c('Hylochoerus', 'Potamochoerus'), c(50, 55)), # Wild pigs.
      list('Kobus', c(49, 51)), # Waterbuck.
      list('Herpestidae', c(31, 33)), # Mongooses.
      list('Macroscelididae', c(19.4, 23.6)), # Elephant shrews.
      list('Manidae', c(5, 5)), # Pangolins.
      list('Mustelidae', c(24, 26)), list(c('Aonyx', 'Hydrictis'), c(9, 10)), # Weasels/badgers/otters/ferrets/etc.
      list('Nanger', c(80, 83)), list('Tragelaphus', c(69, 70)), # Antelopes.
      list(c('Parahyaena', 'Proteles'), c(50, 50)), # Hyenas.
      list('Procaviidae', c(27, 29)), # Hyraxes.
      list('Sciuridae', c(16, 21)), list(c('Anomalurus', 'Idiurus', 'Zenkerella'), c(16, 21)), # Squirrels.
      list('Soricidae', c(0.5, 0.75)), list(c('Congosorex', 'Surdisorex', 'Sylvisorex'), c(1.75, 2))) # Shrews.
    for (G in 1:length(Locom_speed_est)) { # For each taxon in "Locom_speed_est", set the speeds of its species to
      # random numbers in between the taxon's upper and lower speed bounds, bounds derived from the literature.
      Idx <- Mammals_final$Order %in% Locom_speed_est[[G]][[1]] | Mammals_final$Family %in% Locom_speed_est[[G]][[1]] |
        Mammals_final$Genus %in% Locom_speed_est[[G]][[1]] # Logical row indices of the taxon in "Mammals_final".
      Mammals_final$Locom_Speed[Idx] <- runif(sum(Idx), Locom_speed_est[[G]][[2]][1], Locom_speed_est[[G]][[2]][2]) }
    
    Locom_speed_known <- list(list('gnou', 90), list('oryx', 70), list('patas', 55), list('quagga', 70),
      list('taurinus', 80), list('thomsonii', 81)) # Species names paired with known speeds from the literature.
    for (S in 1:length(Locom_speed_known)) { # For each species, set its speed to its known one from the literature.
      Mammals_final$Locom_Speed[Mammals_final$Species == Locom_speed_known[[S]][[1]]] <- Locom_speed_known[[S]][[2]] }
    Mammals_final$Locom_Speed[Mammals_final$Genus == 'Nanger' & Mammals_final$Species == 'granti'] <- 81
    
    Mammals_final$Locom_Speed[match(Mammals_speeds$Binomial, paste(Mammals_final$Genus, Mammals_final$Species,
      sep = '_'))] <- Mammals_speeds$MaxSpeed # Set the known speeds of species recorded in "Mammals_speeds".
    
    Locom_rangecor[Counter] <- cor(log(Mammals_final$Locom_Speed), log(Mammals_final$GeogRange), use = 'complete.obs')
    Locom_masscor[Counter] <- cor(log(Mammals_final$Locom_Speed), log(Mammals_final$BodyMass), use = 'complete.obs')
    Counter <- Counter + 1 } # Update the counter.
  rm(Tmp, Tmp2, Counter, Seed, G, Idx, S)
  set.seed(Random_seed) # Recalibrate seed to "Random_seed".

# Before moving forward, now is the time to split "Mammals_final" into training and testing subsets. The reason is that
# certain procedures below (e.g., imputation, scaling) use all data in a given dataframe, and thus creating the subsets
# later would lead to data leakage between them. Therefore, split "Mammals_final" into the subsets now. Follow the
# protocol of Roberts et al., 2017, to make the split while accounting for the phylogenetic non-independence of species:

  # Read in necessary libraries, for here and/or in general:  
    
  library(caret) # For general machine learning processes moving forward, in case such models are built.
  library(phytools); library(ape) # For generally conducting phylogenetic analyses.
  library(cluster) # For using the "pam" function below.

  # Sub-sample 50 random phylogenetic trees out of the 1,000 provided by Phylacine, and prune each down to Africa-only
  # species (this code block was run once and commented out, with its outputs saved, to save time and space). Select one:

  # Mammals_phylo <- read.nexus('Complete_phylogeny.nex') # Read in all 1,000 trees.
  # Mammals_phylo <- Mammals_phylo[sample(1:length(Mammals_phylo), 50, replace = FALSE)] # Sub-sample.
  # Mammals_phylo <- lapply(Mammals_phylo, keep.tip, tip = Mammals_final$Binomial) # Prune down to Africa-only species.
  # class(Mammals_phylo) <- 'multiPhylo' # Change from list to multiphylo object.
  # writeNexus(Mammals_phylo, 'AfricaMammalsPhylo_FromPhylacine_V2.nex') # Save sub-sampled trees.
  Mammals_phylo <- read.nexus(list.files()[grep('Phylo_FromPhylacine_V2', list.files())]) # Read in.
  Phylosig_treenum <- 13 # Tree number used here for the data split, as well as for the GLS models built later.

  # Using the selected tree, calculate the phylogenetic distances in branch lengths between species. Use the resulting
  # distance matrix as input to a function that splits all species into two clusters based on their distances, with
  # species that are closer together (more evolutionarily related) in the same cluster:
  
  Phylosig_distmat <- cophenetic.phylo(Mammals_phylo[[Phylosig_treenum]]) # Distance matrix.
  ML_datasplit <- pam(Phylosig_distmat, k = 2, cluster.only = TRUE) # Splitting/clustering function.
  ML_trainset <- Mammals_final[ML_datasplit == 1, ]; ML_testset <- Mammals_final[ML_datasplit != 1, ] # Apply split.
  
# Perform a random forest imputation on "Mammals_final", "ML_trainset", and "ML_testset" to fill in their missing values:

  # Import the relevant library for random forest imputation:

  library(missForest)

  # Create a function that performs the imputation per dataframe provided:

  Imputer <- function(Df) { # "Df" = the provided dataframe.
    
    # Create a subset of "Df" that contains the columns that will be used for imputing. "GeogRange" is excluded, because
    # that is the outcome variable of interest, and any variables with few unique values are also excluded, as they
    # provide little model value. Of the taxonomic variables, only "Order" is kept, as the others have too many factor
    # levels for "missForest" to handle:
    
    Impute_subset <- Df[,setdiff(colnames(Df),c('Binomial', 'Family', 'Genus', 'Species', 'Locom_FootPos', 'GeogRange'))]

    # Drop unused levels from all categorical variables in "Impute_subset":
  
    for (Col in colnames(Impute_subset)[sapply(Impute_subset, is.factor)]) { # For each categorical variable...
      Impute_subset[, Col] <- droplevels(Impute_subset[, Col]) }; rm(Col)
    
    # Perform the imputation on "Impute_subset", using the default hyperparameters of the "missForest" function:
    
    Impute_subsetfull <- missForest(Impute_subset, maxiter = 3, variablewise = TRUE)
    
    # Integrate the imputed columns from "Impute_subsetfull" into "Df", keeping a desired column order. Then return "Df":
    
    Df_full <- cbind(Impute_subsetfull$ximp$Order, Df[, -c(which(colnames(Df) %in% colnames(Impute_subsetfull$ximp)),
      grep('GeogRange|Locom_FootPos', colnames(Df)))], Impute_subsetfull$ximp[, -which(colnames(
      Impute_subsetfull$ximp) == 'Order')], Df$GeogRange)
    if (!('Locom_FootPos' %in% colnames(Df_full))) { Df_full <- add_column(Df_full, Locom_FootPos = Df$Locom_FootPos,
      .before = 'Locom_Speed') } # Add back in the "Locom_FootPos" column.
    return(Df_full) } # Return the fully imputed and processed dataframe.
  
  # Apply the "Imputer" function to each dataframe of interest:
  
  Mammals_finalfull <- Imputer(Mammals_final); ML_trainset <- Imputer(ML_trainset); ML_testset <- Imputer(ML_testset)
  
# Perform final modification and cleaning steps on "Mammals_finalfull", "ML_trainset", and "ML_testset" before saving 
# "Mammals_finalfull" as a CSV file and moving forward with data analysis:

Cleaner <- function(Df) { # Create a function that performs these cleaning steps per dataframe.
  Df <- Df[, c(2,1,3:13,16,14:15,17:ncol(Df))] # Change order of dataframe columns as needed.
  colnames(Df)[grep('Order|GeogRange', colnames(Df))] <- c('Order', 'GeogRange') # Cleaner column names.
  Df$PctPAOverlap[Df$PctPAOverlap > 1] <- 1 # Two values of this column are > 1. Fix that, as this is a % column.
  Df <- Df[order(Df$Order, Df$Family, Df$Genus, Df$Species), ] # Order rows based upon taxonomic names of species.
  return(Df) } # Return the cleaned dataframe.
  
Mammals_finalfull <- Cleaner(Mammals_finalfull); ML_trainset <- Cleaner(ML_trainset); ML_testset <- Cleaner(ML_testset)
  
# Save "Mammals_finalfull" as a CSV file:

# write.csv(Mammals_finalfull, 'AfricaMammalsAttributesAndRanges_ByMe_V2.csv', row.names = FALSE)

##########################################################################################################################
##########################################################################################################################

# PART II: DATA MANIPULATION AND EXPLORATION

# Before performing any data exploration on "Mammals_finalfull", remove and/or change certain of its rows and columns.
# Do the same for "ML_trainset" and "ML_testset", such that all dataframes remain consistent:

Rowcol_changer <- function(Df) { # Create a function that performs row and column manipulations per dataframe.

  # Simplify the categories in "PrimaryBiome" down to three, easy-to-comprehend categories:
  
  Df$PrimaryBiome <- as.factor(ifelse(Df$PrimaryBiome %in% c(1:6, 12, 14), 'Forest', ifelse(Df$PrimaryBiome %in% 7:11,
    'Savanna', 'Desert')))
  
  # Simplify the categories in "Locom_FootPos" to refer to species that are vs. are not plantigrade, as the other foot
  # positions exhibit similar species range sizes and thus can be combined:
  
  Df$Locom_FootPos <- factor(ifelse(Df$Locom_FootPos == 'P', 'P', 'NP'))
  
  return(Df) } # Return the processed dataframe.

Mammals_finalfull <- Rowcol_changer(Mammals_finalfull)
ML_trainset <- Rowcol_changer(ML_trainset); ML_testset <- Rowcol_changer(ML_testset)

# Visualize the distributions of all numeric/continuous variables in "Mammals_finalfull":

Mammals_finalfull_num <- Mammals_finalfull[, sapply(Mammals_finalfull, is.numeric)] # Dataframe of numeric variables.
Mammals_finalfull_nump <- Mammals_finalfull_num[, -which(colnames(Mammals_finalfull_num) == 'GeogRange')] # Predictors.

par(mfrow = c(4,3), mar = c(6,6,2,1)) # Set the plot panel to have X rows and Y columns, and set plot margins.
for (Var in colnames(Mammals_finalfull_num)) { # For each variable in "Mammals_finalfull_num"...
  
  # Option 1 - make a histogram of the variable:
  
  # hist(Mammals_finalfull_num[, Var], xlab = '', ylab = '', main = Var, cex.main = 1.5, cex.axis = 1.5, col = 'cyan3')
  
  # Option 2 - make a density plot of the variable, with a polygon included to fill in the plot with color:
  
  # plot(density(Mammals_finalfull_num[, Var]), xlab = '', ylab = '', main = Var, cex.main = 1.5, cex.axis = 1.5)
  # polygon(density(Mammals_finalfull_num[, Var]), col = 'cadetblue3')
  
  # Option 3 - make a Q-Q plot of the variable, which shows the correlation between quantiles generated from the
  # distribution of the variable and those generated from a given normal distribution:
  
  qqnorm(Mammals_finalfull_num[, Var], xlab = '', ylab = '', main = gsub('.*_', '', Var), cex.main = 2, cex.axis = 1.5, 
    cex = 1.5, pch = 19)
  qqline(Mammals_finalfull_num[, Var], col = 'cyan3', lwd = 3)
  
  mtext('Theoretical', side = 1, line = 3, font = 1, cex = 1.5) # X-axis label (OR 'Value').
  mtext('Actual', side = 2, line = 3, font = 1, cex = 1.5) }; rm(Var) # Y-axis label (OR 'Frequency').

# Statistically determine the numeric/continuous variables in "Mammals_finalfull" that have non-normal distributions: 

Pval_method <- 'fdr' # Method used to correct/adjust p-values, to avoid type I errors in multiple comparisons.
Alpha <- 0.05 # Significance level used for various analyses below.

suppressWarnings({ # Remove mention of inconsequential warnings here that are by-products of any K-S test.
  Non_normal <- list() # List to hold K-S test p-values for each numeric variable.
  for (Col in 1:ncol(Mammals_finalfull_num)) { # For each variable, get mean p-value for ensemble of random K-S tests...
    Non_normal[[Col]] <- mean(replicate(100, { ks.test(rnorm(nrow(Mammals_finalfull_num), mean = mean(
      Mammals_finalfull_num[, Col]), sd = sd(Mammals_finalfull_num[, Col])), Mammals_finalfull_num[, Col])$p.value }))
    } }); rm(Col)

Non_normal <- p.adjust(unlist(Non_normal), method = Pval_method) # Correct the p-values because of multiple comparisons.
Non_normal_pvals <- Non_normal # Store p-values before modifying "Non_normal" further.
Non_normal <- colnames(Mammals_finalfull_num)[which(Non_normal < Alpha)] # Which variables are non-normal.
Non_normal_pred <- Non_normal[-which(Non_normal == 'GeogRange')] # Which predictors are non-normal.

# Transform those non-normal numeric/continuous variables in "Mammals_finalfull", "ML_trainset", and "ML_testset" to make
# them normally distributed using a Box-Cox transformation, estimating its "lambda" parameter automatically: 

BoxCox_transformer <- function(Df, Cols) { # Create a function that performs this transformation per dataframe.
  Df[, Cols] <- sapply(Df[, Cols], function(Col) { # For each numeric/continuous column in the dataframe...
    BoxCox(Col, lambda = BoxCoxLambda(Col, method = 'guerrero', lower = -5, upper = 5)) }) # Perform the transform.
  return(Df) } # Return the transformed dataframe.

Mammals_finalfull_num <- BoxCox_transformer(Mammals_finalfull_num, Non_normal) # Apply function.
Mammals_finalfull_nump[, Non_normal_pred] <- Mammals_finalfull_num[, Non_normal_pred] # Integrate into "nump".
Mammals_finalfull[, Non_normal] <- Mammals_finalfull_num[, Non_normal] # Integrate into "Mammals_finalfull".

ML_trainset <- BoxCox_transformer(ML_trainset, Non_normal); ML_testset <- BoxCox_transformer(ML_testset, Non_normal)

# Z-score scale all numeric/continuous predictors for each of "Mammals_finalfull", "ML_trainset", and "ML_testset":

Mammals_finalfull_nump <- scale(Mammals_finalfull_nump)
Mammals_finalfull_num[, colnames(Mammals_finalfull_nump)] <- Mammals_finalfull_nump # Integrate into "num".
Mammals_finalfull[, colnames(Mammals_finalfull_nump)] <- Mammals_finalfull_nump # Integrate into full dataset.

ML_trainset[, colnames(Mammals_finalfull_nump)] <- scale(ML_trainset[, colnames(Mammals_finalfull_nump)])
ML_testset[, colnames(Mammals_finalfull_nump)] <- scale(ML_testset[, colnames(Mammals_finalfull_nump)])

# Make plots that explore links between and distributions of all numeric/continuous features in "Mammals_finalfull":

library(psych) # For using the "pairs.panels" function below.
par(mar = c(7,7,2,2)) # Adjust plot margins.

Mammals_finalfull_pairspanels <- Mammals_finalfull[, c(colnames(Mammals_finalfull_nump), 'GeogRange')] # Features.
colnames(Mammals_finalfull_pairspanels) <- gsub('.*_', '', colnames(Mammals_finalfull_pairspanels)) # Modify colnames.

try(silent = TRUE, # Suppress an inconsequential error message about the size of the plot window.
pairs.panels( # Initiate the function.
  Mammals_finalfull_pairspanels, # Select numeric columns with "GeogRange" last, as was set up above.
  method = 'spearman', # Select correlation coefficient as the visual method.
  density = TRUE, ellipses = FALSE, # Set which statistics will be displayed.
  hist.col = 'steelblue3', cex.cor = 5, cex.labels = 1.1, cex.axis = 1.5, pch = 19, gap = 0.4)) # Format plot.

# For each of any categorical variables of interest found in "Mammals_finalfull", make three kinds of plots: 1) a boxplot
# of the geographic range sizes of species per category, 2) a scatterplot showing how the relationship between range size
# vs. a continuous variable of interest may or may not change across categories, and 3) an interaction plot showing how
# the relationship between range size vs. the categorical variable may change depending upon another categorical variable.
# Conduct statistical tests to determine if the slopes in the scatterplots of item 2 are significantly different from 
# each other. The purpose of items 2 and 3 is to identify candidate interaction terms for the models below:

library(ggplot2) # For making plot visualizations.

Overall_catvars <- sapply(Mammals_finalfull[, (which(colnames(Mammals_finalfull) == 'Species') + 1):ncol(
  Mammals_finalfull)], is.factor) # All categorical variables indices.
Overall_catvars <- names(Overall_catvars)[Overall_catvars] # All categorical variables names.

Trait_plot_colors <- c('tan3', 'deepskyblue3', 'plum4') # Colors to use in the plots below.
Trait_plot_1 <- list(); Trait_plot_2 <- list(); Trait_plot_3 <- list(); Trait_stats <- list() # Lists to hold outputs.
Counter <- 1 # Counter to store items in lists through the loop.

for (Char in Overall_catvars) { # For each categorical variable...
  
  Trait_plot_1[[Counter]] <- ggplot(Mammals_finalfull, aes_string(x = Char, y = 'GeogRange')) + # Select the data.
    geom_boxplot(outlier.size = 3, lwd = 2, fill = 'deepskyblue3') + # Set to be a boxplot with specific settings.
    labs(x = gsub('.*_', '', Char), y = 'Range Size (Box-Cox Sq Km)') + # Set axis labels.
    theme_classic() + # Make the plot have a classic-looking theme.
    theme(axis.text = element_text(size = 30), axis.title = element_text(size = 35), axis.title.x = element_text(
      margin = ggplot2::margin(t = 20)), axis.title.y = element_text(margin = ggplot2::margin(r = 20))) # Format axes.
  
  Var_interest_cont <- 'TempSeas' # Or another continuous variable of interest to analyze below.
  
  Trait_plot_2[[Counter]] <- ggplot(Mammals_finalfull, aes_string(Var_interest_cont, 'GeogRange', col = Char)) + # Data.
    geom_point(size = 3) + geom_smooth(method = 'lm', size = 3) + # Set to be a scatterplot with a regression line.
    labs(x = paste(Var_interest_cont, '(Box-Cox Units)'), y = 'Range Size (Box-Cox Sq Km)', col = gsub(
      '(?!^)(?=[[:upper:]])', '\n', gsub('.*_', '', Char), perl = TRUE)) + # Set axis and legend labels, using regex.
    scale_color_manual(labels = as.character(sort(unique(Mammals_finalfull[, Char]))), values = Trait_plot_colors[
      1:length(unique(Mammals_finalfull[, Char]))]) + # Set the scatter plot color scheme.
    theme_classic() + # Make the plot have a classic-looking theme.
    theme(axis.text = element_text(size = 40), axis.title = element_text(size = 40), axis.title.x = element_text(
      margin = ggplot2::margin(t = 20)), axis.title.y = element_text(margin = ggplot2::margin(r = 20)), legend.text = 
      element_text(size = 40), legend.title = element_text(size = 40, face = 'bold')) # Format axes and legend.
  
  Trait_stats[[Counter]] <- summary(aov(as.formula(paste('GeogRange ~', Var_interest_cont, '*', Char)),
    data = Mammals_finalfull)) # Run statistical test for determining if slopes in "Trait_plot_2" differ significantly.
  Trait_stats[[Counter]] <- Trait_stats[[Counter]][[1]]$`Pr(>F)`[3] # Extract p-value from statistical test.
  
  Var_interest_cat <- 'PrimaryBiome' # Or another categorical variable of interest to analyze below.
  
  Trait_plot_3[[Counter]] <- ggplot(Mammals_finalfull, aes_string(Var_interest_cat, 'GeogRange', col = Char)) + # Data.
    geom_smooth(aes_string(group = Char), method = 'lm', size = 3) + # Set to be a line plot.
    labs(x = Var_interest_cat, y = 'Range Size (Box-Cox Sq Km)', col = gsub('(?!^)(?=[[:upper:]])', '\n', gsub(
      '.*_', '', Char), perl = TRUE)) + # Set axis and legend labels, using regex.
    scale_color_manual(labels = as.character(sort(unique(Mammals_finalfull[, Char]))), values = Trait_plot_colors[
      1:length(unique(Mammals_finalfull[, Char]))]) + # Set the scatter plot color scheme.
    theme_classic() + # Make the plot have a classic-looking theme.
    theme(axis.text = element_text(size = 40), axis.title = element_text(size = 40), axis.title.x = element_text(
      margin = ggplot2::margin(t = 20)), axis.title.y = element_text(margin = ggplot2::margin(r = 20)), legend.text = 
      element_text(size = 40), legend.title = element_text(size = 40, face = 'bold')) # Format axes and legend.
  
  Counter <- Counter + 1 }; rm(Char, Counter, Var_interest_cont, Var_interest_cat) # Update the plot counter.

Trait_stats <- p.adjust(Trait_stats, method = Pval_method) # Adjust the p-values because of multiple comparisons.
# Trait_plot_1; Trait_plot_2; Trait_plot_3 # Print all plots on separate panels.

# Plot boxplots of all numeric/continuous variables in "Mammals_finalfull" to check for outliers and to explore:

par(mar = c(9,5,2.1,2.1), mfrow = c(1,1)) # Adjust plot margins.

Mammals_forboxplot <- as.list(Mammals_finalfull[, colnames(Mammals_finalfull_num)]) # List of variables.
Mammals_forboxplot[[which(colnames(Mammals_finalfull_num) == 'GeogRange')]] <- scale(
  Mammals_forboxplot[[which(colnames(Mammals_finalfull_num) == 'GeogRange')]]) # Scale range variable.

boxplot(Mammals_forboxplot, xaxt = 'n', yaxt = 'n', pch = 19, col = 'darkolivegreen3', lwd = 3) # Make the plot.
axis(1, at = 1:(ncol(Mammals_finalfull_num)), labels = FALSE) # Format the X axis ticks.
text(1:(ncol(Mammals_finalfull_num)), par('usr')[3] - 0.6, labels = colnames(Mammals_finalfull_num), 
  pos = 2, offset = 0, xpd = TRUE, srt = 30, cex = 1.5) # Format the X axis tick labels.
axis(2, cex.axis = 2, font = 1) # Format the Y axis ticks.
mtext('Continuous Feature', side = 1, line = 8, cex = 2, font = 2) # Format the X axis label.
mtext('Transformed and Scaled Value', side = 2, line = 3, cex = 2, font = 2) # Format the Y axis label.

# Check for multicollinearity between select predictors in "Mammals_finalfull". To do so, create a diagnostic model
# formula that does not include the variables in "Mammals_finalfull" that are taxonomic or that will not be modeled
# together. Exclude different combinations of those latter variables to see VIF results for each on its own:

Diagnostic_formula <- paste('GeogRange ~', paste(colnames(Mammals_finalfull)[-unname(which(colnames(Mammals_finalfull)
  %in% c('Binomial', 'Order', 'Family', 'Genus', 'Species', 'GeogRange')))], collapse = ' + '))
library(car) # For using the "vif" function below.
VIF <- vif(lm(Diagnostic_formula, data = Mammals_finalfull)) # Calculate variance inflation factor per predictor.

# Assign categories to the predictors, based upon their reference to taxonomic, biological, or environmental factors:

Predictors <- colnames(Mammals_finalfull)
Predictor_catIDs <- c('Taxonomic', 'BioMass', 'BioDiet', 'BioLocom', 'EnviroOther', 'EnviroHet', 'Response')
Predictor_catNames <- c('Binomial|Order|Family|Genus|Species', 'Body', 'TrophicLevel', 'Locom',
  'Temp|Precip|HumPopDen|PA|MidLat|Biome|Suit', 'Het', 'GeogRange')

for (Idx in 1:length(Predictor_catNames)) { # For each category...
  names(Predictors)[grep(Predictor_catNames[Idx], Predictors)] <- Predictor_catIDs[Idx] }; rm(Idx) # Assign predictors.

# Check for a phylogenetic signal in each numeric predictor in "Mammals_finalfull" by calculating each's Pagel's lambda:
  
  # For each predictor, calculate Pagel's lambda from three random trees in "Mammals_phylo", compiled previously (this 
  # code block was run once and commented out, with its outputs saved, to save time and space):
  
  # Phylosig_pagel <- c() # Initiate a vector that will hold Pagel's lambda values per predictor.
  # for (Pred in 1:ncol(Mammals_finalfull_nump)) { # For each predictor...
  #   Tmp <- Mammals_finalfull[, colnames(Mammals_finalfull_nump)[Pred]] # Variable representing the predictor's values.
  #   names(Tmp) <- Mammals_finalfull$Binomial # Names of "Tmp" as scientific names.
  #   Phylosig_pagel[Pred] <- mean(sapply(sample(1:length(Mammals_phylo), 3, replace = FALSE), function(Tree) {
  #     phylosig(Mammals_phylo[[Tree]], Tmp, method = 'lambda')$lambda }), na.rm = TRUE) # Calculate lambda from trees.
  #   print(paste('Predictor', Pred, 'complete')) }; rm(Pred, Tmp) # Print a record of this loop's progress.
  # write.table(Phylosig_pagel, 'AfricaMammalsPagel_ByMe_V2.txt', row.names = FALSE, col.names = FALSE) # Save values.
  Phylosig_pagel <- read.table(list.files()[grep('Pagel_ByMe_V2', list.files())]) # Read in Pagel's lambda values.
  rownames(Phylosig_pagel) <- colnames(Mammals_finalfull_nump) # Assign values to appropriate predictors.
  
  # Make a bar plot of the calculated Pagel's lambda values per predictor:
  
  Phylosig_pagel$Category <- names(Predictors)[match(rownames(Phylosig_pagel), Predictors)] # Predictor categories.
  Phylosig_pagel$Color <- ifelse(grepl('Bio', Phylosig_pagel$Category), 'orange3', 'cadetblue3') # Colors for bars.
  
  rownames(Phylosig_pagel) <- gsub('.*_', '', rownames(Phylosig_pagel)) # Clean up the names of the predictors.
  Phylosig_pagel <- Phylosig_pagel[rev(order(Phylosig_pagel[,1])), ] # Order the rows by Pagel's lambda value.
  
  par(mar = c(11,8,2.1,2.1)) # Adjust plot margins.
  Pagel_bar <- barplot(Phylosig_pagel[,1], xaxt = 'n', pch = 19, col = Phylosig_pagel$Color, cex.axis = 3) # Make plot.
  text(Pagel_bar, 0, labels = rownames(Phylosig_pagel), pos = 2, offset = 0, xpd = TRUE, srt = 90, cex = 2) # X ticks.
  mtext("Pagel's Lambda", side = 2, line = 5, cex = 3) # Format the Y axis label.

##########################################################################################################################
##########################################################################################################################

# PART III: MACHINE LEARNING AND MODELING

# Upload the library or libraries that will be used for machine learning and model building:

library(nlme) # For the "gls" (phylogenetic generalized least squares) model.

# Make necessary data type conversions in "Mammals_finalfull", "ML_trainset", and "ML_testset" before proceeding:

Converter <- function(Df) { # Create a function that performs data type conversions per dataframe.
  Df <- as.data.frame(unclass(Df)) # Convert all character columns to factor.
  return(Df) } # Return the processed dataframe.
  
Mammals_finalfull <- Converter(Mammals_finalfull)
ML_trainset <- Converter(ML_trainset); ML_testset <- Converter(ML_testset)

# Using an AIC framework, determine the formula that will be used in the models below:

  # Create an initial control formula with select terms in "Mammals_finalfull", excluding the terms of interest:

  Predictors_interest <- Predictors[names(Predictors) %in% c('EnviroHet')] # Terms of interest.

  ML_formula_control <- as.formula(paste('GeogRange ~ ', # Response term.
    paste(setdiff(Predictors[!(names(Predictors) %in% c('Taxonomic', 'Response'))], Predictors_interest),
    collapse = ' + '))) # Feature selection - exclude terms of interest.
  
  # Run models for all possible combinations of formula terms in "ML_formula_control". For each model, output its
  # results, including an AIC value. Whichever model has the lowest AIC value fits the data the best:
  
  library(MuMIn) # For using the "dredge" function below.

  ML_formula_options <- dredge(lm(ML_formula_control, data = Mammals_finalfull, na.action = na.fail), rank = 'AICc')
  
  # Create a vector of all possible formula terms that could be selected by the final formula below:
  
  ML_formula_termoptions <- colnames(ML_formula_options)[2:(which(colnames(ML_formula_options) == 'df') - 1)]
  
  # Choose the final formula by choosing the terms that were used together in the model with the lowest AIC, or in the
  # models with an AIC within a magnitude of X of the lowest. Add in interaction terms to the formula as appropriate. 
  # The formula will also be used as a baseline to add in and analyze terms of interest:
  
  ML_formula_terms <- ML_formula_termoptions[sapply(ML_formula_options[ML_formula_options$delta < 2,
    ML_formula_termoptions], function(Col) { !all(is.na(Col)) })]
  
  ML_formula_terms_interaction <- c() # Vector to hold potential interaction terms.
  for (Term in c('BodyMass:Locom_FootPos', 'Locom_FootPos:MeanTemp', 'OccursInPA:PrecipSeas', 'PrecipSeas:TempSeas', 
    'PrecipSeas:PrimaryBiome', 'PrimaryBiome:TempSeas', 'PrimaryBiome:SuitChange')) { # For each interaction term...
    if (strsplit(Term, ':')[[1]][1] %in% ML_formula_terms & strsplit(Term, ':')[[1]][2] %in% ML_formula_terms) {
      ML_formula_terms_interaction <- c(ML_formula_terms_interaction, Term) }} # Include the term if its components exist.
    
  ML_formula <- paste('GeogRange ~', paste(c(ML_formula_terms, ML_formula_terms_interaction), collapse = ' + '))
  ML_formula_orig <- ML_formula # Keep a record of the initial formula, before changes to it are possibly made.
  
  # Based on the performance of the various models below, perform feature selection by adding or removing terms from
  # "ML_formula", such that an updated, optimized version of "ML_formula" is used for all models:
  
  # ML_formula <- gsub('TrophicLevel \\+ ', '', ML_formula)
  
# Run models on all data in "Mammals_finalfull" using "ML_formula". Run a control model containing no terms of interest,
# then iteratively add one term of interest at a time into "ML_formula" and run the resulting models. Obtain and 
# visualize outputs from each model, and analyze how they compare to each other:

  # Prepare technical elements needed for model building:
  
  library(FSA) # For using the "lrt" function below.
  AIC_performance <- c(); AIC_vals <- c(); AIC_pvals <- list() # Vectors and list to hold metrics outputted by the models.
  Counter <- 1 # Counter for the "Char" loop below.
  
  rownames(Mammals_finalfull) <- Mammals_finalfull$Binomial # Needed for GLS.
  rownames(ML_trainset) <- ML_trainset$Binomial # Perform same manipulation here.
  rownames(ML_testset) <- ML_testset$Binomial # Perform same manipulation here.

  Phylosig_pagel_value <- median(Phylosig_pagel$V1) # Value used in GLS models. Try different values in ML to optimize.

  # Build and analyze the models:
  
  for (Char in unname(c('Control', Predictors_interest))) { # For the control or model containing a term of interest...
    
    # Build and evaluate the models that include an added term of interest:
    
    if (Char != 'Control') { # If the control model is not being addressed...
      
      Model_notcontrol <- gls( # Run a generalized least squares model that accounts for phylogeny.
        as.formula(paste(ML_formula, '+', Char, # Add the term of interest to "ML_formula".
          ifelse(grepl('NA', Char), paste('+', Char, ':MeanPrecip'), ifelse(grepl('NA', Char), paste('+', Char,
          ':MeanTemp'), '')))), # Possibly include the term in an interaction (or if not, put "NA" here).
        data = Mammals_finalfull, # Choose the dataset to be used for the model.
        correlation = corPagel(value = Phylosig_pagel_value, # Account for phylogeny, setting Pagel's lambda and...
          phy = keep.tip(Mammals_phylo[[Phylosig_treenum]], rownames(Mammals_finalfull)), # ...using "Mammals_phylo".
          fixed = FALSE), # Leave the "corPagel" "value" above unfixed to be optimized, or fixed to save time.
        method = 'ML') # Employ maximum likelihood optimization.
      
      AIC_pvals[[Counter]] <- lrt(Model_notcontrol, com = Model_control) # Statistically compare the model to the control.
      AIC_pvals[[Counter]] <- AIC_pvals[[Counter]][length(AIC_pvals[[Counter]])] } # Pull out p-value from comparison.
    
    # Build the control model that does not include an added term of interest:
    
    else { # If the control model is being addressed...
      
      Model_control <- gls(as.formula(ML_formula), data = Mammals_finalfull, correlation = corPagel(value =
        Phylosig_pagel_value, phy = keep.tip(Mammals_phylo[[Phylosig_treenum]], rownames(Mammals_finalfull)),
        fixed = FALSE), method = 'ML') # Perform the same GLS model as above, but do not add any term of interest.
      
      AIC_pvals[[Counter]] <- 1 } # A p-value is not relevant to the control, so fill this first list value with a "1".
    
    # Perform further evaluations on the built model:
    
    if (exists('Model_notcontrol')) { Model <- Model_notcontrol } else { Model <- Model_control } # General model var.
    AIC_performance[Counter] <- cor(fitted(Model), Mammals_finalfull$GeogRange)^2 # R^2 of model predictions.
    AIC_vals[Counter] <- AICc(Model) # AIC value of the model.
    
    # Plot the coefficients of the model's terms:
    
    par(mar = c(20, 12, 0, 0)) # Adjust plot margins.
    Tmp <- Model$coefficients # Temporary variable for the coefficients.
    Tmp <- Tmp[!(grepl('Intercept', names(Tmp)))] # Remove intercept coefficient.
    Tmp <- Tmp[rev(order(abs(Tmp)))] # Order coefficients from largest to smallest in magnitude.
    
    Tmp2 <- summary(Model)$tTable # Table with coefficient p-values.
    Tmp2 <- cbind(Tmp2, pvalue_adj = p.adjust(Tmp2[, 'p-value'], method = Pval_method)) # Add adjusted p-values to "Tmp2".
    Tmp2 <- Tmp2[match(names(Tmp), rownames(Tmp2)), ] # Match "Tmp2" rows to the coefficients in "Tmp".
    print(c(Char, Tmp2[grepl(Char, rownames(Tmp2)), c('Value', 'p-value')])) # Print out coef and p-value of interest.
    Tmp2_orig <- which(Tmp2[, 'p-value'] < Alpha) # Indices of coefficients with significant p-values before adjustment.
    Tmp2_adj <- which(Tmp2[, 'pvalue_adj'] < Alpha) # Indices with significant p-values after adjustment.

    Tmp3 <- intervals(Model)$coef # Temporary variable for the coefficients' 95% confidence intervals.
    Tmp3 <- Tmp3[match(names(Tmp), rownames(Tmp3)), ] # Match confidence intervals to the coefficients in "Tmp".
    
    Tmp4 <- quantile(c(0, max(abs(Tmp), na.rm = TRUE)), 0.1) # Position of first stars in barplots (see below).
    Tmp5 <- quantile(c(0, max(abs(Tmp), na.rm = TRUE)), 0.3) # Position of second stars (see below).
    
    Bar <- barplot(Tmp, ylim = c(-max(abs(Tmp), na.rm = TRUE), max(abs(Tmp), na.rm = TRUE) + 1), xaxt = 'n',
      cex.axis = 5, col = ifelse(grepl(paste(Predictors_interest, collapse = '|'), names(Tmp)), 'cornflowerblue',
      'peachpuff3')) # Make the plot, color-coding the bars by predictors of interest vs. not.
    text(Bar, 0, labels = gsub('Locom_', '', names(Tmp)), pos = 2, offset = 0, xpd = TRUE, srt = 90, cex = 2) # Ticks.
    text(Bar[Tmp2_orig], Tmp4, '*', cex = 5) # Put star by bars with significant coefficient p-values before adjustment.
    text(Bar[Tmp2_adj], Tmp5, '*', cex = 5) # Put another star by bars with significant p-values after adjustment.
    # arrows(Bar, Tmp3[, 'lower'], Bar, Tmp3[, 'upper'], lwd = 2, angle = 90, code = 3, length = 0.05) # Conf intervals.
    # mtext('Attribute', side = 1, line = 10, cex = 5, font = 2) # Format X axis label.
    mtext('Coefficient', side = 2, line = 8, cex = 5) # Format Y axis label.
    
    Counter <- Counter + 1 } # Update counter.
  
  rm(Char, Counter, Model_notcontrol, Model_control, Model, Tmp, Tmp2, Tmp2_orig, Tmp2_adj, Tmp3, Tmp4)
  
  # Prepare the AIC values and associated p-values for visualization:
  
  AIC_pvals <- unlist(AIC_pvals) # Unlist the p-values representing comparisons of non-control to control models.
  AIC_pvals_orig <- AIC_pvals # Save before adjusting below.
  AIC_pvals <- p.adjust(AIC_pvals, method = Pval_method) # Adjust the p-values because of multiple comparisons.
  names(AIC_vals) <- unname(c('Control', gsub('.*_', '', Predictors_interest))) # Match terms to their AICs.
  names(AIC_pvals) <- unname(c('Control', gsub('.*_', '', Predictors_interest))) # Match terms to their AIC p-values.
  AIC_vals <- sort(AIC_vals, decreasing = TRUE) # Sort AIC values from highest to lowest.
  AIC_pvals <- AIC_pvals[match(names(AIC_vals), names(AIC_pvals))] # Sort AIC p-values in the same order.
  AIC_vals <- AIC_vals - AIC_vals[names(AIC_vals) == 'Control'] # Convert AICs to absolute diff from control.
  
  # Visualize the AIC values:
  
  par(mar = c(5,12,1,1), mfrow = c(1,1)) # Adjust plot margins.
  AIC_bar <- barplot(AIC_vals, ylim = c(min(AIC_vals) - 10, max(AIC_vals) + 10), xaxt = 'n', cex.axis = 5,
    col = ifelse(grepl('Control', names(AIC_vals)), 'coral3', ifelse(grepl(paste(gsub('.*_', '', Predictors_interest),
      collapse = '|'), names(AIC_vals)), 'cornflowerblue', 'peachpuff3'))) # Color AIC value bars by group.
  AIC_bar_names <- names(AIC_vals) # Designate the names associated with each bar.
  AIC_bar_names[grep('Hab', AIC_bar_names)] <- 'Habitat'; AIC_bar_names[grep('Topo', AIC_bar_names)] <- 'Topographic'
  text(AIC_bar, 0, labels = AIC_bar_names, pos = 2, offset = 0, xpd = TRUE, srt = 90, cex = 5) # Format X ticks.
  text(AIC_bar[names(AIC_vals) != 'Control' & AIC_pvals < Alpha], 5, '*', cex = 10) # Put stars by bars with sig p-vals.
  # mtext('Attribute', side = 1, line = 3, cex = 5, font = 2) # Format X axis label.
  mtext('Delta AICc Value', side = 2, line = 8, cex = 4.5, font = 2) # Format Y axis label ("Value" or "Percent").

# Run models as above, but this time using a train/test machine learning framework. Train the models on "ML_trainset"
# and test them on "ML_testset". Do so for a control model with no terms of interest included, as well as for additional
# models that each have one term included:

ML_trained <- list() # List to hold trained models.
ML_tested <- list() # List to hold predictions made by the trained models.
ML_evaluated <- list() # List to hold RMSE, R^2, and MAE values stemming from evaluating models' predictions.
ML_abline <- list() # List to hold linear regression lines comparing predicted vs. actual response values.

ML_terms <- c('Control', unname(Predictors_interest)) # Terms defining each model iteration below.
# ML_terms <- 'Control' # (Uncomment this if just running a control model iteration).

for (Term in ML_terms) { # For each term referring to a modeling iteration...
  
  # Update "ML_formula" to include "Term" both on its own and potentially as an interaction with another term (or if no
  # interaction, then "NA"), unless "Term" is referring to the "Control" model, in which case leave "ML_formula" as-is:
  
  if (Term == 'Control') { ML_formula_tmp <- as.formula(ML_formula) } else { ML_formula_tmp <- as.formula(paste(
    ML_formula, ' + ', Term, ifelse(grepl('NA', Term), paste('+', Term, ':MeanPrecip'), ifelse(grepl('NA', Term),
    paste('+', Term, ':MeanTemp'), '')), sep = '')) }

  # Train a model using "ML_formula_tmp" on "ML_trainset", similar to above with all data in "Mammals_finalfull":
  
  ML_trained[[Term]] <- gls(ML_formula_tmp, # Formula to use.
    data = ML_trainset, # Choose the training dataset.
    correlation = corPagel(value = Phylosig_pagel_value, phy = keep.tip(Mammals_phylo[[Phylosig_treenum]],
      rownames(ML_trainset)), fixed = FALSE), # Account for phylogeny using "Mammals_phylo" (try different trees).
    method = 'ML') # Employ maximum likelihood optimization.

  # Make predictions on "ML_testset" with the trained model. Possibly post-process the predictions if they exhibit bias:
  
  ML_tested[[Term]] <- predict(ML_trained[[Term]], ML_testset) * 0.8

  # Evaluate the model by calculating the model's RMSE, R^2, and MAE, as well as by creating a regression plot with
  # predicted model values versus actual values:
  
  Col <- 'GeogRange' # Which column to use for evaluation.
  ML_evaluated[[Term]] <- postResample(pred = ML_tested[[Term]], obs = ML_testset[, Col]) # Perform the evaluation.
  
  par(mar = c(9.5, 9.5, 2.5, 1)) # Adjust plot margins.
  plot(ML_tested[[Term]], ML_testset[, Col], xlab = '', ylab = '', xlim = range(c(ML_testset[, Col], 
    ML_tested[[Term]]), na.rm = TRUE), ylim = range(c(ML_testset[, Col], ML_tested[[Term]]), na.rm = TRUE), 
    pch = 19, cex.axis = 2.75, cex = 3, mgp = c(3,2,0)) # Make the basic plot, equalizing the ranges of the two axes.
  mtext('Predicted (Box-Cox Sq Km)', side = 1, line = 6, cex = 3.5) # Format the X axis label.
  mtext('Actual (Box-Cox Sq Km)', side = 2, line = 6, cex = 3.5) # Format the Y axis label.
  
  lines(x = c(-100,100), y = c(-100,100), col = 'blue', lwd = 10, lty = 'dashed') # Add a diagonal line to the plot.
  
  ML_abline[[Term]] <- lm(ML_testset[, Col] ~ ML_tested[[Term]]) # Perform linear regression on the above plot.
  abline(ML_abline[[Term]], col = 'red', lwd = 10) # Plot the resultant regression line on the plot.
  
  mtext(paste('R-Squared =', as.character(round(summary(ML_abline[[Term]])$adj.r.squared, 2))), side = 1, line = -28,
    cex = 2.75, font = 1, at = range(c(ML_testset[, Col], ML_tested[[Term]]))[1], adj = 0) # R^2 value in plot.
  mtext(paste('Y = ', as.character(round(ML_abline[[Term]]$coefficients[2], 2)), 'X + ', as.character(round(
    ML_abline[[Term]]$coefficients[1], 2)), sep = ''), side = 1, line = -25, cex = 2.75, font = 1, at = range(c(
    ML_testset[, Col], ML_tested[[Term]]))[1], adj = 0) }; rm(Term, Col) # Equation in plot.

# As was done with AIC above, conduct statistical tests to determine if the performance of the models in the object
# "ML_trained" differ from its first indexed model, which represents the control model:

ML_modelcomp <- sapply(2:length(ML_trained), function(Model) { # For each model that is not the control model...
  lrt(ML_trained[[Model]], com = ML_trained[['Control']]) }) # Compare it to the control.
ML_modelcomp <- as.data.frame(ML_modelcomp); ML_modelcomp <- ML_modelcomp[nrow(ML_modelcomp), ] # Pull out only p-vals.
colnames(ML_modelcomp) <- Predictors_interest # Give names to p-values to associate them with their models.

# Make a plot to compare select values in "ML_evaluated":

ML_evaluated_rsq <- sapply(ML_evaluated, '[[', 'Rsquared') # Pull out the metric of interest from "ML_evaluated".
ML_evaluated_rsq <- rev(sort(ML_evaluated_rsq)) # Sort the metric values from highest to lowest.

Floorer <- function(Val, Digits) { round(Val - 5 * 10^(-Digits - 1), Digits) } # Function to round a number down.
Ceilinger <- function(Val, Digits) { round(Val + 5 * 10^(-Digits - 1), Digits) } # Function to round a number up.

par(mar = c(4, 9.5, 2.5, 1)) # Adjust plot margins ("c(20, 13, 2, 2)" if x-axis text is vertical).
names(ML_evaluated_rsq)[grep('Hab', names(ML_evaluated_rsq))] <- 'Habitat' # Better name for plot bar.
names(ML_evaluated_rsq)[grep('Topo', names(ML_evaluated_rsq))] <- 'Topographic' # Better name for plot bar.
Bar <- barplot(ML_evaluated_rsq, log = 'y', ylim = c(Floorer(min(ML_evaluated_rsq, na.rm = T), 2), Ceilinger(max(
  ML_evaluated_rsq, na.rm = T), 2)), xaxt = 'n', cex.axis = 3.5, col = ifelse(grepl('Control', names(ML_evaluated_rsq)),
  'coral3', 'cornflowerblue')) # Make bar plot, and color it as above.
# text(Bar, min(ML_evaluated_rsq) * 0.99, labels = names(ML_evaluated_rsq), pos = 2, offset = 0, xpd = TRUE,
#   srt = 90, cex = 2) # Format X tick labels (if vertical).
text(Bar, Floorer(min(ML_evaluated_rsq, na.rm = T), 2) - 0.0025, labels = names(ML_evaluated_rsq), xpd = T, cex = 3.5,
  font = 1) # Format X tick labels (if horizontal).
mtext('R-Squared', side = 2, line = 6, cex = 4.5) # Y label.

# Perform similar modeling analyses as above, except delineated by different taxonomic groups of species for comparison:

  # Delineate the species in "Mammals_finalfull" into groups, based on their taxonomic affiliations:

  Mammals_groups <- # Initiate a variable to hold groups. Each newly indented line below refers to a group.
    ifelse(Mammals_finalfull$Order %in% c('Afrosoricida', 'Hyracoidea', 'Macroscelidea', 'Proboscidea',
      'Tubulidentata'), 'Afrotheres', # This line and the one above it are for one group.
    ifelse(Mammals_finalfull$Order == 'Chiroptera', 'Bats',
    ifelse(Mammals_finalfull$Order %in% c('Carnivora', 'Pholidota'), 'Ferae',
    ifelse(Mammals_finalfull$Order %in% c('Rodentia', 'Lagomorpha'), 'Glires',
    ifelse(Mammals_finalfull$Order == 'Eulipotyphla', 'Hedgehogs and Shrews',
    ifelse(Mammals_finalfull$Order == 'Primates', 'Primates',
    ifelse(Mammals_finalfull$Order %in% c('Cetartiodactyla', 'Perissodactyla'), 'Ungulates',
    'Other'))))))) # All other species are miscellaneous.
  
  # Create a formula to use for models below that do not include problematic terms, i.e., terms that are uniform across
  # all species in a given group of mammals:
  
  ML_formula_remove <- 'FootPos' # Phrase(s) for the terms to remove.
  ML_formula_terms_interaction_groups <- ML_formula_terms_interaction[!(grepl(ML_formula_remove,
    ML_formula_terms_interaction))] # Remove relevant interaction terms.
  ML_formula_terms_groups <- ML_formula_terms[!(grepl(ML_formula_remove, ML_formula_terms))] # Remove relevant lone terms.
  ML_formula_groups <- paste('GeogRange ~', paste(c(ML_formula_terms_groups, ML_formula_terms_interaction_groups),
    collapse = ' + ')) # Create a new model formula after appropriate terms removed.
  
  # Create a function that builds a "gls" model for each group in "Mammals_groups", and returns its key outputs:
  
  ML_outputs_groups_names <- c('R-squared', 'Coef_val', 'Coef_stderr', 'Coef_tval', 'Coef_pval', 'Coef_pval_adj',
    'Coef_lowint', 'Coef_upint') # Names of outputs to be collected below.
  
  Modeler <- function(Group) { # Initiate the function with the input being the name of the group of interest.
    Model_outputs <- list(); Counter <- 1 # List to hold key model outputs and counter for the loop below.
    for (Term in ML_terms[!(grepl('Control', ML_terms))]) { # For each term referring to a modeling iteration...
      
      # Build and optimize the model for the group and term of interest:
      
      ML_formula_tmp <- as.formula(paste(ML_formula_groups, ' + ', Term, sep = '')) # Set the model iteration's formula.
      if (Group == 'All') { Group_data <- Mammals_finalfull } else { # Select all species OR...
        Group_data <- Mammals_finalfull[Mammals_groups == Group, ] } # ...select only the spp in the group of interest.
      if (grepl('Shrews', Group)) { Group_data <- Group_data[-which(Group_data$PrimaryBiome == 'Desert'), ] }
        # For the "Eulipotyphla" group, just one species is in a desert. This singleton causes model error - remove it.
      if (grepl('Bats', Group)) { Group_data <- Group_data[-which(Group_data$TrophicLevel == 'Omnivore'), ] }
        # For the "Bats" group, just one species is an omnivore. This singleton causes model error - remove it.
      Model <- gls(ML_formula_tmp, data = Group_data, # Build the model with only the species in the group.
        correlation = corPagel(value = Phylosig_pagel_value, phy = keep.tip(Mammals_phylo[[Phylosig_treenum]],
          rownames(Group_data)), fixed = FALSE), # Below, set up model optimization parameters to avoid coding error.
        control = glsControl(maxIter = 1000, msMaxIter = 1000, msTol = 1e-3, tolerance = 1e-3, opt = ifelse(grepl(
          'Ungulates|Bats', Group) & Analyses_res != 50, 'optim', 'nlminb'))) # "Optim" helps certain models converge.
      
      # Return key outputs from the built model:
      
      R_squared <- cor(Group_data$GeogRange, Model$fitted)^2 # Obtain the r-squared of the model.
      Model_table <- summary(Model)$tTable # Obtain the model's summary table, from which further outputs will be taken.
      Model_table <- cbind(Model_table, p.adjust(Model_table[, 'p-value'], method = Pval_method)) # Add corrected p-vals.
      Coef_stats <- Model_table[rownames(Model_table) == Term, ] # Coefficient information for the term "Term".
      Coef_ints <- intervals(Model)$coef; Coef_ints <- unname(Coef_ints[rownames(Coef_ints) == Term, c('lower', 'upper')])
      Model_outputs[[Counter]] <- c(R_squared, Coef_stats, Coef_ints) # Compile the key model outputs obtained above.
      names(Model_outputs[[Counter]]) <- ML_outputs_groups_names # Output names.
      Counter <- Counter + 1 }; rm(Term, Counter) # Update the counter.
    
    return(Model_outputs) } # Return the key model outputs.
  
  # Apply the above function to build models for mammalian groups and obtain model outputs (this code block was run once
  # and commented out, with its outputs saved, to save time and space):
  
  # ML_outputs_groups <- lapply(c('All', unique(Mammals_groups)), Modeler) # Apply function.
  # names(ML_outputs_groups) <- c('All', unique(Mammals_groups)) # Give models associated names.
  # save(ML_outputs_groups, file = paste('ML_outputs_groups_', as.character(Analyses_res), '-Km', '.RData', sep = ''))
  load(paste('ML_outputs_groups_', as.character(Analyses_res), '-Km', '.RData', sep = '')) # Load the outputs.
  
  # Make barplots of the model outputs in "ML_outputs_groups":
  
  par(mfrow = c(2,4), mar = c(1,3.5,3,1)) # Adjust plot layout settings.
  lapply(sort(names(ML_outputs_groups[!(grepl('All|Other', names(ML_outputs_groups)))])), function(Group) {
    Coefs <- sapply(ML_outputs_groups[[Group]], '[[', 'Coef_val') # Get coefficients.
    Lower <- sapply(ML_outputs_groups[[Group]], '[[', 'Coef_lowint') # Get coefficient lower confidence limit.
    Upper <- sapply(ML_outputs_groups[[Group]], '[[', 'Coef_upint') # Get coefficient upper confidence limit.
    Sig <- !(Lower < 0 & Upper > 0) # Which coefficients differ significantly from 0?
    Sig_adj <- sapply(ML_outputs_groups[[Group]], '[[', 'Coef_pval_adj') <= 0.05 # Differ significantly after correction?
    Bar <- barplot(Coefs, ylim = c(-5.5, 2), cex.axis = 3.5, xaxt = 'n', col = c('deepskyblue3', 'orange3')) # Bar plot.
    abline(h = 0, lwd = 2); box() # Add horizontal line at y-axis value of 0, as well as a box.
    arrows(Bar, Lower, Bar, Upper, lwd = 4, angle = 90, code = 3, length = 0.05) # Add error bars.
    # text(2, -2.6, labels = paste0('p = ', round(sapply(ML_outputs_groups[[Group]], '[[', 'Coef_pval_adj')[1], 3)),
    #   cex = 3, font = 2) # Put statistical significance label.
    if (TRUE %in% Sig) { text(Bar[Sig], Coefs[Sig], adj = 0, labels = '*', cex = 5) } # Put stars by bars != 0.
    if (TRUE %in% Sig_adj) { text(Bar[Sig_adj], Coefs[Sig_adj], adj = 0, labels = '**', cex = 5) } }) # Stars by adj sig.
