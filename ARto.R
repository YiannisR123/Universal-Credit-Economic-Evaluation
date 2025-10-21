library(cmdstanr)
library(readxl)
library(haven)
library(dplyr )
'%!in%' = function(x,y)!('%in%'(x,y))
data = list()
data_h = list()
wavesletters = letters[1:12]
setwd("C:/Users/Yiannis/Desktop/UKDA-6931-spss/spss/spss28/ukhls")
for(i in seq_along(wavesletters)){
  data[[i]] = read_sav(paste0(wavesletters[i],"_indresp_protect.sav"))
  data[[i]] = data[[i]][order(data[[i]]$pidp),]
  data_h[[i]] = read_sav(paste0(wavesletters[i],"_hhresp_protect.sav"))
  message(i)
}
setwd("C:/Users/Yiannis/Desktop/UKDA-6670-spss/spss/spss28/ukhls")
data_lsoa = list()
for(i in seq_along(wavesletters)){
  data_lsoa[[i]] = read_sav(paste0(wavesletters[i],"_lsoa01_protect.sav"))
  message(i)
}

library(stringr)
library(spdep)
data.spatial.dir = setwd("C:/Users/Yiannis/Desktop/UKDA-6670-spss/spss/spss28/ukhls")

# spatial polygons at lsoa and lad levels
poly.lsoa <- sf::st_read(dsn = data.spatial.dir, layer = 'LSOA_2011_EW_BGC_V3')
A <- sf::st_read(dsn = data.spatial.dir, layer = 'Wards_December_2017_FCB_in_Great_Britain')
B <- sf::st_read(dsn = data.spatial.dir, layer = 'Wards_December_2016_FCB_in_Great_Britain')
B = B[,-c(3,4)]
E <- sf::st_read(dsn = data.spatial.dir, layer = 'Wards_December_2017_FEB_in_Great_Britain')
colnames(A) = colnames(poly.lsoa)
colnames(B) = colnames(poly.lsoa)
colnames(E) = colnames(poly.lsoa)

# england only
poly.lsoa = rbind(poly.lsoa,A,B,E)
poly.lsoa.england <- poly.lsoa
poly.lsoa.england <- poly.lsoa.england[!duplicated(poly.lsoa.england$LSOA11CD),] 

# make a-mats
lsoa.mat <- spdep::poly2nb(as(poly.lsoa.england, 'Spatial'))
lsoa.mat <- spdep::nb2mat(lsoa.mat, zero.policy = TRUE)
colnames(lsoa.mat) <- rownames(lsoa.mat) <- paste0('lsoa_', 1:dim(lsoa.mat)[1])

# info dataframe to keep track 
lsoa.names <- data.frame(LSOA11CD = poly.lsoa.england$LSOA11CD,
                         LSOA11NM = poly.lsoa.england$LSOA11NM,
                         Internal = rownames(lsoa.mat)) %>%
  mutate(space_id = 1:n())
gc()
#Individuals id to keep across years, filter out retired, disabled, and chronically sick individuals
{
  IdcFilterOut = c()
  IdcFilterOut <- unlist(lapply(seq_along(wavesletters), function(i) {
    data[[i]]$pidp[
      data[[i]][[paste0(wavesletters[i], "_age_dv")]] > 64 |
        data[[i]][[paste0(wavesletters[i], "_age_dv")]] < 16 |
        data[[i]][[paste0(wavesletters[i], "_jbstat")]] %in% c(4, 8) |
        data[[i]][[paste0(wavesletters[i], "_intdaty_dv")]] < 0
    ]
  }))
  IdcToKeep = c()
  IdcToKeep <- unlist(lapply(seq_along(wavesletters), function(i) {
    data[[i]]$pidp
  }))
  IdcFilterOut = IdcFilterOut %>% unique() %>% sort()
  IdcToKeep =  IdcToKeep %>% unique() %>% sort()
  IdcToKeepUnq = c()
  IdcToKeepUnq = IdcToKeep[c(IdcToKeep %!in% IdcFilterOut)]
}

#Extract variables, from indresp and hhresp, link them based on individuals across years
{
  SFS_mat = matrix(NA, nrow = length(IdcToKeepUnq), ncol = length(wavesletters))
  # Use outer loop and Map for efficiency
  for (j in seq_along(wavesletters)) {
    # Extract pidp and finnnow data for the current wave
    pidp_current = data[[j]]$pidp
    finnnow_current = data[[j]][[paste0(wavesletters[j], "_finnow")]]
    m_current = data[[j]][[paste0(wavesletters[j], "_intdaty_dv")]]-2008
    matched_indices = match(pidp_current,IdcToKeepUnq)
    matched_indices_id = matched_indices[!is.na(matched_indices)]
    finnnow_current = finnnow_current[!is.na(matched_indices)]
    m_current = m_current[!is.na(matched_indices)]
    
    matched_indices_id = matched_indices_id[!( m_current>=13)]
    finnnow_current = finnnow_current[!(m_current>=13)]
    m_current = m_current[!( m_current>=13)]
    
    # Populate the matrix column for the current wave
    SFS_mat[matched_indices_id, m_current] = finnnow_current
    message(paste("Wave", j, "completed"))
  }
  SFS_mat[SFS_mat < 0] = NA
  
  # Identify rows where all non-NA values are in the unwanted set
  rows_to_keep <- apply(is.na(SFS_mat),1,mean) != 1
  # Filter the matrix and the vector
  SFS_mat_filtered <- SFS_mat[rows_to_keep, ]
  IdcToKeepUnq <- IdcToKeepUnq[rows_to_keep]
  # Initialize the list of matrices
  vars <- c("SEX", "AGE", "ETH", "MT", "HEQ", "IS", "EA", "GR",
            "UC","GHQ","LS","PSU","STRATA")
  mat_list <- setNames(
    lapply(vars, function(x) matrix(NA, nrow = length(IdcToKeepUnq), ncol = length(wavesletters))),
    vars
  )
  
  # Define a mapping for variable suffixes
  var_suffixes <- c(
    SEX = "_sex", AGE = "_age_dv", ETH = "_ethn_dv",
    MT = "_mastat_dv", HEQ = "_qfhigh_dv",
    IS = "_sclfsat2", EA = "_jbstat", GR = "_gor_dv",
    UC = "_benbase4",GHQ = "_scghq2_dv",
    LS ="_sclfsato",
    PSU = "_psu", STRATA = "_strata"
  )
  
  # Loop over waves
  for (j in seq_along(wavesletters)) {
    
    
    # Populate matrices dynamically
    for (var in vars) {
      if(var != "UC"){
        wave <- wavesletters[j]
        pidp_current = data[[j]]$pidp
        finnnow_current = data[[j]][[paste0(wave, var_suffixes[var])]]
        m_current = data[[j]][[paste0(wavesletters[j], "_intdaty_dv")]]-2008
        matched_indices = match(pidp_current,IdcToKeepUnq)
        matched_indices_id = matched_indices[!is.na(matched_indices)]
        finnnow_current = finnnow_current[!is.na(matched_indices)]
        m_current = m_current[!is.na(matched_indices)]
        
        matched_indices_id = matched_indices_id[!(m_current>=13)]
        finnnow_current = finnnow_current[!( m_current>=13)]
        m_current = m_current[!(m_current>=13)]
        
        
        mat_list[[var]][matched_indices_id, m_current] <- finnnow_current
      }else if(var == "UC" & j >= 6){
        wave <- wavesletters[j]
        pidp_current = data[[j]]$pidp
        finnnow_current = data[[j]][[paste0(wave, var_suffixes[var])]]
        m_current = data[[j]][[paste0(wavesletters[j], "_intdaty_dv")]]-2008
        matched_indices = match(pidp_current,IdcToKeepUnq)
        matched_indices_id = matched_indices[!is.na(matched_indices)]
        finnnow_current = finnnow_current[!is.na(matched_indices)]
        m_current = m_current[!is.na(matched_indices)]
        
        matched_indices_id = matched_indices_id[!(m_current>=13)]
        finnnow_current = finnnow_current[!(m_current>=13)]
        m_current = m_current[!(m_current>=13)]
        
        
        mat_list[[var]][matched_indices_id, m_current] <- finnnow_current
      }
    }
    message("Wave ", j, " completed")
  }
  # Unpack variables from mat_list
  list2env(mat_list, envir = .GlobalEnv)
  
  matrix_list <- lapply(c("LSOA"), function(x) {
    matrix(NA, nrow = length(IdcToKeepUnq), ncol = length(wavesletters))
  })
  names(matrix_list) <- c( "LSOA")
  
  # Loop through each wave
  for (j in seq_along(wavesletters)) {
    wave_letter <- wavesletters[j]
    
    # Extract necessary columns for this wave
    pidp_wave <- data[[j]]$pidp
    hidp_wave <- data[[j]][[paste0(wave_letter, "_hidp")]]
    
    hidp_h_wave <- data_lsoa[[j]][[paste0(wave_letter, "_hidp")]]
    lsoa_h_wave <- data_lsoa[[j]][[paste0(wave_letter, "_lsoa01")]]
    
    # Match and filter indices
    match_indices <- match(IdcToKeepUnq,pidp_wave)
    hid = hidp_wave[match_indices[!is.na(match_indices)]]
    
    match_indices1 <- match(pidp_wave,IdcToKeepUnq)
    match_indices1 =  match_indices1[!is.na(match_indices1)]
    #valid_matches <- !is.na(match_indices)
    #hid <- ifelse(valid_matches, hidp_wave[match_indices[valid_matches]], NA)
    # Filter household data by matching IDs
    #match_indices_h <- match(hidp_h_wave,hid[valid_matches],)
    match_indices_h <- match(hid,hidp_h_wave)
    
    #matrix_list$LSOA[valid_matches, j] <- lsoa_h_wave[match_indices_h]
    matrix_list$LSOA[match_indices1, j] <- lsoa_h_wave[match_indices_h]
    message("Wave ", j, " completed")
  }
  # Extract matrices from the list
  list2env(matrix_list, envir = .GlobalEnv)
}    

#rest of benefits
{
  # Initialize the list of matrices
  vars <- c("CH",  "ESA")
  mat_list <- setNames(
    lapply(vars, function(x) matrix(NA, nrow = length(IdcToKeepUnq), ncol = length(wavesletters))),
    vars
  )
  
  # Define a mapping for variable suffixes
  var_suffixes <- c(CH = "_benctc",  ESA = "_nfh33")
  
  # Loop over waves
  for (j in seq_along(wavesletters)) {
    
    
    # Populate matrices dynamically
    for (var in vars) {
      if(var != "ESA"){
        wave <- wavesletters[j]
        pidp_current = data[[j]]$pidp
        finnnow_current = data[[j]][[paste0(wave, var_suffixes[var])]]
        m_current = data[[j]][[paste0(wavesletters[j], "_intdaty_dv")]]-2008
        matched_indices = match(pidp_current,IdcToKeepUnq)
        matched_indices_id = matched_indices[!is.na(matched_indices)]
        finnnow_current = finnnow_current[!is.na(matched_indices)]
        m_current = m_current[!is.na(matched_indices)]
        
        matched_indices_id = matched_indices_id[!(m_current<0 | m_current>=13)]
        finnnow_current = finnnow_current[!(m_current<0 | m_current>=13)]
        m_current = m_current[!(m_current<0 | m_current>=13)]
        mat_list[[var]][matched_indices_id, m_current] <- finnnow_current
      }else if(var == "ESA" & j > 1){
        wave <- wavesletters[j]
        pidp_current = data[[j]]$pidp
        finnnow_current = data[[j]][[paste0(wave, var_suffixes[var])]]
        m_current = data[[j]][[paste0(wavesletters[j], "_intdaty_dv")]]-2008
        matched_indices = match(pidp_current,IdcToKeepUnq)
        matched_indices_id = matched_indices[!is.na(matched_indices)]
        finnnow_current = finnnow_current[!is.na(matched_indices)]
        m_current = m_current[!is.na(matched_indices)]
        
        matched_indices_id = matched_indices_id[!(m_current<0 | m_current>=13)]
        finnnow_current = finnnow_current[!(m_current<0 | m_current>=13)]
        m_current = m_current[!(m_current<0 | m_current>=13)]
        mat_list[[var]][matched_indices_id, m_current] <- finnnow_current        
      }
    }
    message("Wave ", j, " completed")
  }
  # Unpack variables from mat_list
  list2env(mat_list, envir = .GlobalEnv)
  
  CH[CH %in% c(-9,-8,-7,-2,-1)] = NA#CH[CH %in% c(-9,-8,-2,-1)] = NA
  ESA[ESA %in% c(-10,-9, -8, -7, -2, -1)] = NA#ESA[ESA %in% c(-9, -8, -2, -1)] = NA
  CH[CH == 2] = 0 # does not receive
  ESA[ESA == 2] = 0#does not receive
  
  # Initialize the list of matrices
  vars <- c("W")
  mat_list <- setNames(
    lapply(vars, function(x) matrix(NA, nrow = length(IdcToKeepUnq), ncol = length(wavesletters))),
    vars
  )
  
  # Define a mapping for variable suffixes
  for(j in seq_along(wavesletters)){
    if(j <= 5){
      var_suffixes <- c(W = "_btype6")
      
      
      # Populate matrices dynamically
      for (var in vars) {
        wave <- wavesletters[j]
        pidp_current = data[[j]]$pidp
        finnnow_current = data[[j]][[paste0(wave, var_suffixes[var])]]
        m_current = data[[j]][[paste0(wavesletters[j], "_intdaty_dv")]]-2008
        matched_indices = match(pidp_current,IdcToKeepUnq)
        matched_indices_id = matched_indices[!is.na(matched_indices)]
        finnnow_current = finnnow_current[!is.na(matched_indices)]
        m_current = m_current[!is.na(matched_indices)]
        
        matched_indices_id = matched_indices_id[!(m_current<0 | m_current>=13)]
        finnnow_current = finnnow_current[!(m_current<0 | m_current>=13)]
        m_current = m_current[!(m_current<0 | m_current>=13)]
        mat_list[[var]][matched_indices_id, m_current] <- finnnow_current
      }
    }else{
      var_suffixes <- c(W = "_nfh20")
      
      # Populate matrices dynamically
      for (var in vars) {
        wave <- wavesletters[j]
        pidp_current = data[[j]]$pidp
        finnnow_current = data[[j]][[paste0(wave, var_suffixes[var])]]
        m_current = data[[j]][[paste0(wavesletters[j], "_intdaty_dv")]]-2008
        matched_indices = match(pidp_current,IdcToKeepUnq)
        matched_indices_id = matched_indices[!is.na(matched_indices)]
        finnnow_current = finnnow_current[!is.na(matched_indices)]
        m_current = m_current[!is.na(matched_indices)]
        
        matched_indices_id = matched_indices_id[!(m_current<0 | m_current>=13)]
        finnnow_current = finnnow_current[!(m_current<0 | m_current>=13)]
        m_current = m_current[!(m_current<0 | m_current>=13)]
        mat_list[[var]][matched_indices_id, m_current] <- finnnow_current
      }  
    }
    message("Wave ", j, " completed")
  }
  list2env(mat_list, envir = .GlobalEnv)
  
  W[W %in% c(-10, -9, -8, -7, -2, -1)] = NA
  W[W == 2] = 0 
  
  
  
  
  # Initialize the list of matrices
  vars <- c("HB")
  mat_list <- setNames(
    lapply(vars, function(x) matrix(NA, nrow = length(IdcToKeepUnq), ncol = length(wavesletters))),
    vars
  )
  
  # Define a mapping for variable suffixes
  for(j in seq_along(wavesletters)){
    if(j <= 5){
      var_suffixes <- c(HB = "_benhou1")
      
      # Populate matrices dynamically
      for (var in vars) {
        wave <- wavesletters[j]
        pidp_current = data[[j]]$pidp
        finnnow_current = data[[j]][[paste0(wave, var_suffixes[var])]]
        m_current = data[[j]][[paste0(wavesletters[j], "_intdaty_dv")]]-2008
        matched_indices = match(pidp_current,IdcToKeepUnq)
        matched_indices_id = matched_indices[!is.na(matched_indices)]
        finnnow_current = finnnow_current[!is.na(matched_indices)]
        m_current = m_current[!is.na(matched_indices)]
        
        matched_indices_id = matched_indices_id[!(m_current<0 | m_current>=13)]
        finnnow_current = finnnow_current[!(m_current<0 | m_current>=13)]
        m_current = m_current[!(m_current<0 | m_current>=13)]
        mat_list[[var]][matched_indices_id, m_current] <- finnnow_current
      }
    }else{
      var_suffixes <- c(HB = "_nfh22")
      
      # Populate matrices dynamically
      for (var in vars) {
        wave <- wavesletters[j]
        pidp_current = data[[j]]$pidp
        finnnow_current = data[[j]][[paste0(wave, var_suffixes[var])]]
        m_current = data[[j]][[paste0(wavesletters[j], "_intdaty_dv")]]-2008
        matched_indices = match(pidp_current,IdcToKeepUnq)
        matched_indices_id = matched_indices[!is.na(matched_indices)]
        finnnow_current = finnnow_current[!is.na(matched_indices)]
        m_current = m_current[!is.na(matched_indices)]
        
        matched_indices_id = matched_indices_id[!(m_current<0 | m_current>=13)]
        finnnow_current = finnnow_current[!(m_current<0 | m_current>=13)]
        m_current = m_current[!(m_current<0 | m_current>=13)]
        mat_list[[var]][matched_indices_id, m_current] <- finnnow_current
      }  
    }
    message("Wave ", j, " completed")
  }
  list2env(mat_list, envir = .GlobalEnv)
  
  HB[HB %in% c(-10, -9, -8, -7, -2, -1)] = NA
  HB[HB == 2] = 0 
  
  
  # Initialize the list of matrices
  vars <- c("JSA")
  mat_list <- setNames(
    lapply(vars, function(x) matrix(NA, nrow = length(IdcToKeepUnq), ncol = length(wavesletters))),
    vars
  )
  
  # Define a mapping for variable suffixes
  for(j in seq_along(wavesletters)){
    if(j <= 5){
      var_suffixes <- c(JSA = "_benunemp1")
      
      # Populate matrices dynamically
      for (var in vars) {
        wave <- wavesletters[j]
        pidp_current = data[[j]]$pidp
        finnnow_current = data[[j]][[paste0(wave, var_suffixes[var])]]
        m_current = data[[j]][[paste0(wavesletters[j], "_intdaty_dv")]]-2008
        matched_indices = match(pidp_current,IdcToKeepUnq)
        matched_indices_id = matched_indices[!is.na(matched_indices)]
        finnnow_current = finnnow_current[!is.na(matched_indices)]
        m_current = m_current[!is.na(matched_indices)]
        
        matched_indices_id = matched_indices_id[!(m_current<0 | m_current>=13)]
        finnnow_current = finnnow_current[!(m_current<0 | m_current>=13)]
        m_current = m_current[!(m_current<0 | m_current>=13)]
        mat_list[[var]][matched_indices_id, m_current] <- finnnow_current
      }
    }else{
      var_suffixes <- c(JSA = "_nfh16")
      
      
      # Populate matrices dynamically
      for (var in vars) {
        wave <- wavesletters[j]
        pidp_current = data[[j]]$pidp
        finnnow_current = data[[j]][[paste0(wave, var_suffixes[var])]]
        m_current = data[[j]][[paste0(wavesletters[j], "_intdaty_dv")]]-2008
        matched_indices = match(pidp_current,IdcToKeepUnq)
        matched_indices_id = matched_indices[!is.na(matched_indices)]
        finnnow_current = finnnow_current[!is.na(matched_indices)]
        m_current = m_current[!is.na(matched_indices)]
        
        matched_indices_id = matched_indices_id[!(m_current<0 | m_current>=13)]
        finnnow_current = finnnow_current[!(m_current<0 | m_current>=13)]
        m_current = m_current[!(m_current<0 | m_current>=13)]
        mat_list[[var]][matched_indices_id, m_current] <- finnnow_current
      }  
    }
    message("Wave ", j, " completed")
  }
  list2env(mat_list, envir = .GlobalEnv)
  
  JSA[JSA %in% c(-10, -9, -8, -7, -2, -1)] = NA
  JSA[JSA == 2] = 0 
  
  
  vars <- c("ISS")
  mat_list <- setNames(
    lapply(vars, function(x) matrix(NA, nrow = length(IdcToKeepUnq), ncol = length(wavesletters))),
    vars
  )
  
  # Define a mapping for variable suffixes
  for(j in seq_along(wavesletters)){
    if(j <= 5){
      var_suffixes <- c(ISS = "_btype2")
      
      
      # Populate matrices dynamically
      for (var in vars) {
        wave <- wavesletters[j]
        pidp_current = data[[j]]$pidp
        finnnow_current = data[[j]][[paste0(wave, var_suffixes[var])]]
        m_current = data[[j]][[paste0(wavesletters[j], "_intdaty_dv")]]-2008
        matched_indices = match(pidp_current,IdcToKeepUnq)
        matched_indices_id = matched_indices[!is.na(matched_indices)]
        finnnow_current = finnnow_current[!is.na(matched_indices)]
        m_current = m_current[!is.na(matched_indices)]
        
        matched_indices_id = matched_indices_id[!(m_current<0 | m_current>=13)]
        finnnow_current = finnnow_current[!(m_current<0 | m_current>=13)]
        m_current = m_current[!(m_current<0 | m_current>=13)]
        mat_list[[var]][matched_indices_id, m_current] <- finnnow_current
      }
    }else{
      var_suffixes <- c(ISS = "_nfh15")
      
      
      # Populate matrices dynamically
      for (var in vars) {
        wave <- wavesletters[j]
        pidp_current = data[[j]]$pidp
        finnnow_current = data[[j]][[paste0(wave, var_suffixes[var])]]
        m_current = data[[j]][[paste0(wavesletters[j], "_intdaty_dv")]]-2008
        matched_indices = match(pidp_current,IdcToKeepUnq)
        matched_indices_id = matched_indices[!is.na(matched_indices)]
        finnnow_current = finnnow_current[!is.na(matched_indices)]
        m_current = m_current[!is.na(matched_indices)]
        
        matched_indices_id = matched_indices_id[!(m_current<0 | m_current>=13)]
        finnnow_current = finnnow_current[!(m_current<0 | m_current>=13)]
        m_current = m_current[!(m_current<0 | m_current>=13)]
        mat_list[[var]][matched_indices_id, m_current] <- finnnow_current
      }  
    }
    message("Wave ", j, " completed")
  }
  list2env(mat_list, envir = .GlobalEnv)
  
  ISS[ISS %in% c(-10, -9, -8, -7, -2, -1)] = NA
  ISS[ISS == 2] = 0 
  
  calculate_Ben <- function(...) {
    matrices <- list(...)
    
    # Create an array from the matrices (3D array: 61283 x 10 x 6)
    mat_array <- array(unlist(matrices), dim = c(nrow(matrices[[1]]), ncol(matrices[[1]]), length(matrices)))
    
    # Apply sum logic to each cell across all matrices
    Ben <- apply(mat_array, c(1, 2), function(x) {
      if (all(is.na(x))) {
        return(NA)
      } else {
        # Sum ignoring NAs
        sum_x <- sum(x, na.rm = TRUE)
        return(ifelse(sum_x > 0, 1, 0))
      }
    })
    
    return(Ben)
  }
  
  Benefits_legwel = calculate_Ben(ISS,JSA,HB,W,CH,ESA)
  
  Ben_hist = matrix(NA,dim(UC)[1],dim(UC)[2])
  Ben_hist = Benefits_legwel
  Ben_hist[UC == 1] = 2
}
#replace with NA, missing value or any negative value corresponding to no observation
{
  SEX[SEX %in% c(-9, -1, -2)] = NA
  HEQ[HEQ %in% c(-9, -8)] = NA
  HEQ[HEQ == 96] = 17
  EA[EA %in% c(-9, -8, -2, -1)] = NA
  EA[EA == 97] = 14
  IS[IS %in% c(-10, -9, -8, -7, -2, -1)] = NA
  GR[GR == -9] = NA
  GR[GR == 12] = NA # remove north ireland because in the UCclaims they do not exist
  ETH[ETH == -9] = NA
  ETH[ETH == 97] = 18
  MT[MT %in% c(-9,-8,-2,-1)] = NA
  GHQ[GHQ<0]= NA
  LS[LS<0] = NA
  
  id_1 = LS %in% c(1,2,3) #dissat
  id_2 = LS %in% c(4,5,6,7)
  LS[id_1] = 1
  LS[id_2] = 0
  
  SFS_mat_filtered[which(SFS_mat_filtered <= 3)] = 0
  SFS_mat_filtered[which(SFS_mat_filtered>3)] = 1
  
  GHQ[GHQ <= 3] = 0
  GHQ[GHQ>3] = 1
  
  GHQ_SFS = matrix(NA,dim(GHQ)[1],dim(GHQ)[2])
  GHQ_SFS[GHQ == 0 & SFS_mat_filtered == 0 & LS == 0] = 1 #AR
  GHQ_SFS[GHQ == 1 ] = 2 #MD
  GHQ_SFS[GHQ == 0 & SFS_mat_filtered == 1 & LS == 0] = 3 #FD
  GHQ_SFS[GHQ == 0 & SFS_mat_filtered == 0 & LS == 1] = 4 #LDS
  GHQ_SFS[GHQ == 0 & SFS_mat_filtered == 1 & LS == 1] = 5 #FD & LDS
  
  YY = matrix(NA,dim(GHQ)[1],dim(GHQ)[2])
  for(i in 1:dim(GHQ)[1]){
    for(tt in 2:dim(GHQ)[2]){
      if(GHQ_SFS[i,tt] == 1  & is.na(GHQ_SFS[i,tt])!=TRUE){
        YY[i,tt] = 1
      }else if(GHQ_SFS[i,tt] == 2  & is.na(GHQ_SFS[i,tt])!=TRUE){
        YY[i,tt] = 2
      }else if(GHQ_SFS[i,tt] == 3  & is.na(GHQ_SFS[i,tt])!=TRUE){
        YY[i,tt] = 3
      }else if(GHQ_SFS[i,tt] == 4  & is.na(GHQ_SFS[i,tt])!=TRUE){
        YY[i,tt] = 4
      }else if(GHQ_SFS[i,tt] == 5  & is.na(GHQ_SFS[i,tt])!=TRUE){
        YY[i,tt] = 5
      }
    }
  }
  # 
  # G_F = matrix(NA,dim(GHQ)[1],dim(GHQ)[2])
  # G_F[GHQ == 0 & SFS_mat_filtered == 0 & LS == 0] = 1
  # G_F[!(GHQ == 0 & SFS_mat_filtered == 0 & LS == 0)] = 0
}
gc()
#remove any row which has only NA
{
  matrices_to_filter <- list(
    #G_F = G_F,
    SEX = SEX,
    AGE = AGE,
    HEQ = HEQ,
    #EA = EA,
    IS = IS,
    ETH = ETH,
    MT = MT,
    GR = GR,
    PSU = PSU,
    STRATA = STRATA,
    LSOA = LSOA,
    Ben_hist = Ben_hist,
    GHQ_SFS = GHQ_SFS#,
    #YY = YY
  )
  vars_to_check <- names(matrices_to_filter)
  # Iteratively filter rows where all values in the specified variables are NA
  for (var in vars_to_check) {
    rows_with_only_NA <- rowSums(!is.na(matrices_to_filter[[var]])) == 0
    matrices_to_filter <- lapply(matrices_to_filter, function(mat) mat[!rows_with_only_NA, ])
    IdcToKeepUnq <- IdcToKeepUnq[!rows_with_only_NA]
  }
  # Unpack filtered matrices back to their original variables
  list2env(matrices_to_filter, envir = .GlobalEnv)
}

#group the levels of variables
{
  grp_1 = c(11,16,17) #Below GCSE and othe
  grp_2 = c(7,8,9,10,12,13,14,15) #GCSE, A-level or equivalent
  grp_3 = c(1,2,3,4,5,6) #Degree or higher
  idc_grp_1 = which(HEQ %in% grp_1, arr.ind =  TRUE)
  idc_grp_2 = which(HEQ %in% grp_2, arr.ind =  TRUE)
  idc_grp_3 = which(HEQ %in% grp_3, arr.ind =  TRUE)
  HEQ[idc_grp_1] = 1
  HEQ[idc_grp_2] = 2
  HEQ[idc_grp_3] = 3
  
  # grp_1 = c(1,2,4:14) #Employed and other
  # grp_2 = c(3) #Unemployed
  grp_1 = c(1,2,5) #Employed and other
  grp_2 = c(3,6,7) #Unemployed
  grp_3 = c(8:14)
  idc_grp_1 = which(EA %in% grp_1, arr.ind =  TRUE)
  idc_grp_2 = which(EA %in% grp_2, arr.ind =  TRUE)
  idc_grp_3 = which(EA %in% grp_3, arr.ind =  TRUE)
  EA[idc_grp_1] = 1
  EA[idc_grp_2] = 2
  EA[idc_grp_3] = NA
  
  grp_1 = c(1,2,4)#c(1,2,4,5,6,7,8) #White
  grp_2 = c(5,6,7,8)#c(9,10,11,12,13) #Mixed
  grp_3 = c(9,10,11,12,13)#c(14,15,16) #ASIAN
  grp_4 = c(14,15,16)#c(17) #Black
  grp_5 = c(17,18)#c(18) #Other
  
  idc_grp_1 = which(ETH %in% grp_1, arr.ind =  TRUE)
  idc_grp_2 = which(ETH %in% grp_2, arr.ind =  TRUE)
  idc_grp_3 = which(ETH %in% grp_3, arr.ind =  TRUE)
  idc_grp_4 = which(ETH %in% grp_4, arr.ind =  TRUE)
  idc_grp_5 = which(ETH %in% grp_5, arr.ind =  TRUE)
  
  ETH[idc_grp_1] = 1
  ETH[idc_grp_2] = 2
  ETH[idc_grp_3] = 3
  ETH[idc_grp_4] = 4
  ETH[idc_grp_5] = 5
  
  grp_1 = c(1,4,5,6,7,8,9)#c(1,2)#Single
  grp_2 = c(2,3,10)#c(3,4,5)#Relationship
  idc_grp_1 = which(MT %in% grp_1, arr.ind =  TRUE)
  idc_grp_2 = which(MT %in% grp_2, arr.ind =  TRUE)
  MT[idc_grp_1] = 1
  MT[idc_grp_2] = 2
  
  grp_1 = c(1,2,3)#dissatisfied
  #grp_2 = c(4)#
  grp_2 = c(4,5,6,7)#not dissatisfied, can be satisfied or not dissatisfied or satisfied
  
  idc_grp_1 = which(IS %in% grp_1, arr.ind =  TRUE)
  idc_grp_2 = which(IS %in% grp_2, arr.ind =  TRUE)
  #idc_grp_3 = which(IS %in% grp_3, arr.ind =  TRUE)
  
  IS[idc_grp_1] = 1
  IS[idc_grp_2] = 2
  #IS[idc_grp_3] = 3
  
  
  
  idc_grp_1 = which(AGE >= 16 & AGE < 25) 
  idc_grp_2 = which(AGE >= 25 & AGE < 35) 
  idc_grp_3 = which(AGE >= 35 & AGE < 45) 
  idc_grp_4 = which(AGE >= 45 & AGE < 55) 
  idc_grp_5 = which(AGE >= 55 & AGE < 65) 
  AGE[idc_grp_1] = 1
  AGE[idc_grp_2] = 2
  AGE[idc_grp_3] = 3
  AGE[idc_grp_4] = 4
  AGE[idc_grp_5] = 5
  
  
  
}

gc()
LSOA_id = LSOA
for(i in seq_along(lsoa.names$LSOA11CD)){
  LSOA_id[which(LSOA == lsoa.names$LSOA11CD[i])] = lsoa.names$space_id[i]
}
LSOA_id_2 = LSOA_id
for(i in 1:12){
  LSOA_id_2[grep("^S",LSOA_id[,i]),i] = NA
  LSOA_id_2[grep("^E",LSOA_id[,i]),i] = NA
  LSOA_id_2[grep("^W",LSOA_id[,i]),i] = NA
  LSOA_id_2[grep("^9",LSOA_id[,i]),i] = NA
}
df = LSOA_id_2
df_numeric <- as.data.frame(apply(df, 2, as.numeric))
df_numeric[df_numeric<0] = NA
df_numeric = as.matrix(df_numeric)

matrices_to_filter <- list(
  #G_F = G_F,
  SEX = SEX,
  AGE = AGE,
  HEQ = HEQ,
  #EA = EA,
  IS = IS,
  ETH = ETH,
  MT = MT,
  GR = GR,
  PSU = PSU,
  STRATA = STRATA,
  df_numeric = df_numeric,
  Ben_hist = Ben_hist,
  GHQ_SFS = GHQ_SFS#,
  #YY = YY 
)
idc_to_rem = which(apply(df_numeric,1,sum,na.rm=TRUE) != 0)
rows_with_only_NA <- idc_to_rem 
matrices_to_filter <- lapply(matrices_to_filter, function(mat) mat[rows_with_only_NA, ])
IdcToKeepUnq <- IdcToKeepUnq[rows_with_only_NA]
# Unpack filtered matrices back to their original variables
list2env(matrices_to_filter, envir = .GlobalEnv)
RR = lsoa.mat[sort(unique(as.vector(df_numeric))),sort(unique(as.vector(df_numeric)))]
RR[RR!=0] = 1
df_unq = sort(unique(as.vector(df_numeric)))
df_n = 1:length(df_unq)
for(i in seq_along(df_n)){
  df_numeric[which(df_numeric == df_unq[i])] = i
}

new_matrix <- matrix(NA, nrow = nrow(GHQ_SFS), ncol = ncol(GHQ_SFS))
for (i in 1:nrow(GHQ_SFS)) {
  true_positions <- which(!is.na(GHQ_SFS[i,]))
  new_matrix[i, 1:length(true_positions)] <- true_positions
}

# Recompute n_T for the current state of SFS_mat_filtered
n_T <- apply(GHQ_SFS, 1, function(row) sum(!is.na(row)))

# Function to generate dummy matrices for a given matrix
generate_dummies <- function(matrix_data) {
  # Get the unique values in the matrix, excluding NA
  unique_values <- unique(matrix_data[!is.na(matrix_data)])
  
  # Initialize a list to store the dummy matrices
  dummy_matrices <- list()
  
  # Loop through each unique value to create a dummy matrix
  for (value in unique_values) {
    # Create a dummy matrix: 1 if value matches, 0 otherwise
    dummy_matrix <- ifelse(matrix_data == value, 1, 0)
    
    # Store the dummy matrix in the list with a label
    dummy_matrices[[paste0("dummy_", value)]] <- dummy_matrix
  }
  
  return(dummy_matrices)
}

Ben_hist_unem = matrix(NA,dim(GR)[1],dim(GR)[2])
Ben_hist_unem[Ben_hist == 0 ] = 1
Ben_hist_unem[!(Ben_hist == 0) ] = 0
# Apply the function to each matrix (SEX, HEQ, EA, IS, ETH, MT, SBF)
SEX_dummies <- generate_dummies(SEX)
HEQ_dummies <- generate_dummies(HEQ)
IS_dummies <- generate_dummies(IS)
ETH_dummies <- generate_dummies(ETH)
MT_dummies <- generate_dummies(MT)
EA_dummies <- generate_dummies(EA)
AGE_dummies <- generate_dummies(AGE)
GR_dummies <- generate_dummies(GR)
GHQ_SFS_dummies <- generate_dummies(GHQ_SFS)


M = dim(GR)[1]  

treat = unique(which(Ben_hist == 2, arr.ind = TRUE)[,1])
contr_2 = which(apply(Ben_hist, 1, sum,na.rm=TRUE) == 0)
contr = (1:M)[-c(treat,contr_2)]


Intervention = matrix(NA,dim(GR)[1],dim(GR)[2])
for(r in 1:M){
  for(tt in 1:12){
    if(Ben_hist[r,tt] == 2 & !is.na(Ben_hist[r,tt]) ){
      Intervention[r,tt:12] = 1
    }
    if(Ben_hist[r,tt] != 2 & !is.na(Ben_hist[r,tt])){
      Intervention[r,tt] = 0
    }
  }
}
Intervention_dummies <- generate_dummies(Intervention)

Exposure = matrix(NA,dim(GR)[1],dim(GR)[2])
Exposure[treat,1:12] = 1
Exposure[contr,1:12] = 2 # no ben
Exposure[contr_2,1:12] = 3

Exposure_dummies <- generate_dummies(Exposure)
Exposure_dummies_t <- Exposure_dummies
for(i in 1:M){
  for(t in 1:12){
    Exposure_dummies_t$dummy_1[i,t] = Exposure_dummies_t$dummy_1[i,t]*t
    Exposure_dummies_t$dummy_2[i,t] = Exposure_dummies_t$dummy_2[i,t]*t
    Exposure_dummies_t$dummy_3[i,t] = Exposure_dummies_t$dummy_3[i,t]*t
  }
}
T_mat_plus = matrix(0,M,12)
for(i in 1:M){
  u = 1
  for(t in 1:12){
    if(Intervention[i,t] == 1 & !is.na(Intervention[i,t])){
      T_mat_plus[i,t] = u
      u = u + 1
    }
  }
}
# Exposure_dummies_t_plus <- Exposure_dummies
# for(i in 1:M){
#   for(t in 1:12){
#     Exposure_dummies_t_plus$dummy_1[i,t] = Exposure_dummies_t_plus$dummy_1[i,t]*T_mat_plus[i,t]
#     Exposure_dummies_t_plus$dummy_2[i,t] = Exposure_dummies_t_plus$dummy_2[i,t]*T_mat_plus[i,t]
#     Exposure_dummies_t_plus$dummy_3[i,t] = Exposure_dummies_t_plus$dummy_3[i,t]*T_mat_plus[i,t]
#   }
# }

Intervention_dummies_t_plus <- Intervention_dummies
for(i in 1:M){
  for(t in 1:12){
    Intervention_dummies_t_plus$dummy_1[i,t] = Intervention_dummies_t_plus$dummy_1[i,t]*T_mat_plus[i,t]
    Intervention_dummies_t_plus$dummy_2[i,t] = Intervention_dummies_t_plus$dummy_2[i,t]*T_mat_plus[i,t]
  }
}


gc()
X_matrices <- array(NA, dim = c(M, length(wavesletters), 33))  # M individuals, 21 predictors, TT time points

# Assign dummy variables for each explanatory variable to appropriate positions in the array
X_matrices[,,1] <- SEX_dummies$dummy_2[1:M, ]  # X12 = SEX_dummy_2
X_matrices[,,2] <- HEQ_dummies$dummy_2[1:M, ]  # X32 = HEQ_dummy_2
X_matrices[,,3] <- HEQ_dummies$dummy_3[1:M, ]  # X33 = HEQ_dummy_3
#X_matrices[,,4] <- EA_dummies$dummy_2[1:M, ]   # X42 = EA_dummy_2
X_matrices[,,4] <- IS_dummies$dummy_2[1:M, ]   # X52 = IS_dummy_2
#X_matrices[,,5] <- IS_dummies$dummy_3[1:M, ]  # X53 = IS_dummy_3
X_matrices[,,5] <- ETH_dummies$dummy_2[1:M, ] # X62 = ETH_dummy_2
X_matrices[,,6] <- ETH_dummies$dummy_3[1:M, ] # X63 = ETH_dummy_3
X_matrices[,,7] <- ETH_dummies$dummy_4[1:M, ] # X64 = ETH_dummy_4
X_matrices[,,8] <- ETH_dummies$dummy_5[1:M, ] # X65 = ETH_dummy_5
X_matrices[,,9] <- MT_dummies$dummy_2[1:M, ]  # X72 = MT_dummy_2
X_matrices[,,10] <- AGE_dummies$dummy_2
X_matrices[,,11] <- AGE_dummies$dummy_3
X_matrices[,,12] <- AGE_dummies$dummy_4
X_matrices[,,13] <- AGE_dummies$dummy_5
X_matrices[,,14] <- Intervention_dummies$dummy_1
X_matrices[,,15] <- Exposure_dummies$dummy_1
X_matrices[,,16] <- Exposure_dummies$dummy_2
X_matrices[,,17] <- Exposure_dummies_t$dummy_1
X_matrices[,,18] <- Exposure_dummies_t$dummy_2
# X_matrices[,,21] <- Exposure_dummies_t_plus$dummy_1
# X_matrices[,,22] <- Exposure_dummies_t_plus$dummy_2
# X_matrices[,,23] <- Exp_Int_dummies$dummy_1
# X_matrices[,,24] <- Exp_Int_dummies$dummy_2
# X_matrices[,,25] <- Exp_Int_dummies$dummy_4
X_matrices[,,19] <- Intervention_dummies_t_plus$dummy_1
X_matrices[,,20] <- GR_dummies$dummy_2
X_matrices[,,21] <- GR_dummies$dummy_3
X_matrices[,,22] <- GR_dummies$dummy_4
X_matrices[,,23] <- GR_dummies$dummy_5
X_matrices[,,24] <- GR_dummies$dummy_6
X_matrices[,,25] <- GR_dummies$dummy_7
X_matrices[,,26] <- GR_dummies$dummy_8
X_matrices[,,27] <- GR_dummies$dummy_9
X_matrices[,,28] <- GR_dummies$dummy_10
X_matrices[,,29] <- GR_dummies$dummy_11
X_matrices[,,30] <- GHQ_SFS_dummies$dummy_2
X_matrices[,,31] <- GHQ_SFS_dummies$dummy_3
X_matrices[,,32] <- GHQ_SFS_dummies$dummy_4
X_matrices[,,33] <- GHQ_SFS_dummies$dummy_5


T_ind = matrix(NA,M*12,2)
u = 1
for(i in 1:M){
  for(t in 2:12){
    if(t != 1){
      if(mean(is.na(X_matrices[i,t,]),na.rm=TRUE) == 0 & 
         !is.na(GHQ_SFS[i,t]) & !is.na(GHQ_SFS[i,t-1]) & !is.na(df_numeric[i,t]) & 
         !is.na(PSU[i,t]) &
         !is.na(STRATA[i,t]) ){
        T_ind[u,] = c(i,t)
        u = u + 1
      }
    }else{
      if(mean(is.na(X_matrices[i,t,]),na.rm=TRUE) == 0 & 
         !is.na(GHQ_SFS[i,t]) & !is.na(df_numeric[i,t]) & !is.na(PSU[i,t]) &
         !is.na(STRATA[i,t])){
        T_ind[u,] = c(i,t)
        u = u + 1
      }
    }
  }
}
T_ind = T_ind[which(apply(T_ind,1,sum,na.rm=TRUE)!=0),]
Y = c()
u = 1
for(i in 1:dim(T_ind)[1]){
  Y[u] = GHQ_SFS[T_ind[i,1],T_ind[i,2]]
  u = u + 1
}
X = matrix(NA,dim(X_matrices)[3],length(Y))
u = 1
for(i in 1:dim(T_ind)[1]){
  X[,u] = X_matrices[T_ind[i,1],T_ind[i,2],]
  u = u + 1
}
r = c()
u = 1
for(i in 1:dim(T_ind)[1]){
  r[u] = df_numeric[T_ind[i,1],T_ind[i,2]]
  u = u + 1
}
ps = c()
u = 1
for(i in 1:dim(T_ind)[1]){
  ps[u] = PSU[T_ind[i,1],T_ind[i,2]]
  u = u + 1
}
ps_n = 1:length(unique(ps))
ps_unique = sort(unique(ps))
u = 1
ps_trans = rep(NA,length(ps))
for(i in seq_along(ps_n)){
  ps_trans[which(ps == ps_unique[u])] = i
  u = u + 1
}
PSU_trans = matrix(NA,M,12)
u = 1
for(i in 1:dim(T_ind)[1]){
  PSU_trans[T_ind[i,1],T_ind[i,2]] = ps_trans[u]
  u = u + 1
}
strt = c()
u = 1
for(i in 1:dim(T_ind)[1]){
  strt[u] = STRATA[T_ind[i,1],T_ind[i,2]]
  u = u + 1
}
str_n = 1:length(unique(strt))
str_unique = sort(unique(strt))
u = 1
str_trans = rep(NA,length(strt))
for(i in seq_along(str_n)){
  str_trans[which(strt == str_unique[u])] = i
  u = u + 1
}
STR_trans = matrix(NA,M,12)
u = 1
for(i in 1:dim(T_ind)[1]){
  STR_trans[T_ind[i,1],T_ind[i,2]] = str_trans[u]
  u = u + 1
}
individuals = c()
u = 1
for(i in 1:dim(T_ind)[1]){
  individuals[u] = T_ind[i,1]
  u = u + 1
  
}
T_mat = c()
u = 1
for(i in 1:dim(T_ind)[1]){
  T_mat[u] = T_ind[i,2]
  u = u + 1
}
T_mat_plus_vec = c()
u = 1
for(i in 1:dim(T_ind)[1]){
  T_mat_plus_vec[u] = T_mat_plus[T_ind[i,1],T_ind[i,2]]
  u = u + 1
}

# RR = matrix(c(0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1,
#               1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1,
#               1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0,
#               0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0,
#               0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0,
#               0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0,
#               0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0,
#               0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0,
#               0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0,
#               0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0,
#               1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0), 11, 11)

N_nodes = dim(RR)[1]
N_edges = sum(RR[upper.tri(RR)])
AAAA = as.vector(apply(RR,1,sum))
node_1 = c()
node_2 = c()   
for(i in 1:N_nodes){
  if(AAAA[i] != 0){
    if(length(which(RR[i,(i+1):N_nodes] == 1 ))!=1){
      node_1 = c(node_1,rep(i,length(which(RR[i,(i+1):N_nodes] == 1 ))))
      node_2 = c(node_2,((i+1):N_nodes)[which(RR[i,(i+1):N_nodes] == 1 )])
    }else{
      node_1 = c(node_1,i)
      node_2 = c(node_2,((i+1):N_nodes)[which(RR[i,(i+1):N_nodes] == 1 )])
    }
  }
}

GHQ_SFS[is.na(GHQ_SFS)] = 1
X_matrices[is.na(X_matrices)] = 9999
df_numeric[is.na(df_numeric)] = 0
PSU_trans[is.na(PSU_trans)] = 9999
STR_trans[is.na(STR_trans)] = 9999
data_stan = list(
  K = 5,
  NT = dim(T_ind)[1],
  T_ind = T_ind,
  M = M,
  TT = 12,
  Y = GHQ_SFS,
  n_C = 29,
  X_matrices = X_matrices[,,1:29],
  n_C_t = 4,
  X_matrices_t = X_matrices[,,30:33],
  R = max(unique(r),na.rm=TRUE)+2,
  region = df_numeric,
  PSu = PSU_trans,
  STr = STR_trans,
  Years_plus = T_mat_plus,
  PSU_V = length(unique(ps_trans)),
  STRATA_V = length(unique(str_trans)),
  N_nodes = N_nodes,
  N_edges = N_edges,
  node_1 = node_1,
  node_2 = node_2
)
file <- file.path(cmdstan_path(), "examples", "MarkovModelSFMB",
                  "ResultsUCvsLWvsNo","PosteriorResults", "Results_unif_prior",  "multFDLDS.stan")
mod <- cmdstan_model(file)
gc()
fit_F_ITS = mod$variational(data = data_stan,   save_latent_dynamics = TRUE, draws = 1000)
fit_F_ITS$save_object(file = "outcomes_fit.rds")


b0 = fit_F_ITS$draws()[,8:12]
a0 = fit_F_ITS$draws()[,3:7]
a_coefficients = array(NA,c(1000,29,5))
for(i in 1:29){
  for(k in 1:(5)){
    a_coefficients[,i,k] = fit_F_ITS$draws()[,which(colnames(fit_F_ITS$draws())==paste0("a_coefficients[",i,",",k,"]"))]
    
    
  }
}
a_coefficients_t = array(NA,c(1000,4,5))
for(i in 1:4){
  for(k in 1:(5)){
    a_coefficients_t[,i,k] = fit_F_ITS$draws()[,which(colnames(fit_F_ITS$draws())==paste0("a_coefficients_t[",i,",",k,"]"))]


  }
}
d = array(NA,c(1000,12,5))
for(i in 1:12){
  for(k in 1:(5)){
    d[,i,k] = fit_F_ITS$draws()[,which(colnames(fit_F_ITS$draws())==paste0("d[",i,",",k,"]"))]
    
    
  }
}
g = array(NA,c(1000,16239,5))
for(i in 1:16239){
  for(k in 1:(5)){
    g[,i,k] = fit_F_ITS$draws()[,which(colnames(fit_F_ITS$draws())==paste0("g[",i,",",k,"]"))]
  }
}
e = array(NA,c(1000,M,5))
for(i in 1:M){
  for(k in 1:(5)){
    e[,i,k] = fit_F_ITS$draws()[,which(colnames(fit_F_ITS$draws())==paste0("e[",i,",",k,"]"))]
  }
}
psi = array(NA,c(1000,4683,5))
for(i in 1:4683){
  for(k in 1:(5)){
    psi[,i,k] = fit_F_ITS$draws()[,which(colnames(fit_F_ITS$draws())==paste0("psi[",i,",",k,"]"))]
  }
}
k = array(NA,c(1000,1751,5))
for(i in 1:1751){
  for(j in 1:(5)){
    k[,i,j] = fit_F_ITS$draws()[,which(colnames(fit_F_ITS$draws())==paste0("k[",i,",",j,"]"))]
  }
}

#all FtoF
{
  UC_all = array(NA,c(1000,12,5))
  for(tt in 2:12){
    idc_1 = treat
    idc_1 = idc_1[which(GHQ_SFS[idc_1,tt-1] == 1 & idc_1 %in% T_ind[T_ind[,2]==tt,1])] 
    exp_softmax = array(NA,c(1000,length(idc_1),5))
    for(j in 1:(5)){
      if(length(idc_1)!=1){
      linear_pred = a_coefficients[,,j] %*% t(X_matrices[idc_1,tt,1:29]) + 
        a_coefficients_t[,,j] %*% t(X_matrices[idc_1,tt-1,30:33]) + 
        d[,tt,j] + g[,df_numeric[idc_1,tt],j] + 
        e[,idc_1,j] + as.vector(b0[,j])*tt  +
        as.vector(a0[,j]) + psi[,PSU_trans[idc_1,tt],j] + 
        k[,STR_trans[idc_1,tt],j]
      exp_softmax[,,j] = exp(linear_pred)
      }else{
        linear_pred = a_coefficients[,,j] %*% (X_matrices[idc_1,tt,1:29]) + 
          a_coefficients_t[,,j] %*% (X_matrices[idc_1,tt-1,30:33]) + 
          as.vector(d[,tt,j]) + g[,df_numeric[idc_1,tt],j] + 
          e[,idc_1,j] + as.vector(b0[,j])*tt  +
          as.vector(a0[,j]) + psi[,PSU_trans[idc_1,tt],j] + 
          k[,STR_trans[idc_1,tt],j]
        exp_softmax[,,j] = exp(linear_pred)
      }
    }
    est_prob = array(NA,c(1000,length(idc_1),5))
    for(i in 1:length(idc_1)){
      for(j in 1:5){
        est_prob[,i,j] =  exp_softmax[,i,j]/(exp_softmax[,i,1] + exp_softmax[,i,2] +
                                               exp_softmax[,i,3] + exp_softmax[,i,4] +
                                               exp_softmax[,i,5])
      }
    }
    UC_all[,tt,] = apply(est_prob,c(1,3),mean)
    message(tt)
  }
  saveRDS(UC_all, file = "AR_UC_all_years.RDS")
  
 
  LW_all = array(NA,c(1000,12,5))
  for(tt in 2:12){
    idc_1 = contr
    idc_1 = idc_1[which(GHQ_SFS[idc_1,tt-1] == 1 & idc_1 %in% T_ind[T_ind[,2]==tt,1])] 
    exp_softmax = array(NA,c(1000,length(idc_1),5))
    for(j in 1:(5)){
      if(length(idc_1)!=1){
      linear_pred = a_coefficients[,,j] %*% t(X_matrices[idc_1,tt,1:29]) +
        a_coefficients_t[,,j] %*% t(X_matrices[idc_1,tt-1,30:33]) + 
        as.vector(d[,tt,j]) + g[,df_numeric[idc_1,tt],j] + 
        e[,idc_1,j] + as.vector(b0[,j])*tt  +
        as.vector(a0[,j]) +
        psi[,PSU_trans[idc_1,tt],j] + 
        k[,STR_trans[idc_1,tt],j]
      exp_softmax[,,j] = exp(linear_pred)
      }else{
        linear_pred = a_coefficients[,,j] %*% (X_matrices[idc_1,tt,1:29]) +
          a_coefficients_t[,,j] %*% (X_matrices[idc_1,tt-1,30:33]) + 
          as.vector(d[,tt,j]) + g[,df_numeric[idc_1,tt],j] + 
          e[,idc_1,j] + as.vector(b0[,j])*tt  +
          as.vector(a0[,j]) +
          psi[,PSU_trans[idc_1,tt],j] + 
          k[,STR_trans[idc_1,tt],j]
        exp_softmax[,,j] = exp(linear_pred)
      }
    }
    est_prob = array(NA,c(1000,length(idc_1),5))
    for(i in 1:length(idc_1)){
      for(j in 1:5){
        est_prob[,i,j] =  exp_softmax[,i,j]/(exp_softmax[,i,1] + exp_softmax[,i,2] +
                                               exp_softmax[,i,3] + exp_softmax[,i,4] +
                                               exp_softmax[,i,5])
      }
    }
    LW_all[,tt,] = apply(est_prob,c(1,3),mean)
    message(tt)
  }
  
  saveRDS(LW_all, file = "AR_LW_all_years.RDS")
  

  No_all = array(NA,c(1000,12,5))
  for(tt in 2:12){
    idc_1 = contr_2
    idc_1 = idc_1[which(GHQ_SFS[idc_1,tt-1] == 1 & idc_1 %in% T_ind[T_ind[,2]==tt,1])]
    exp_softmax = array(NA,c(1000,length(idc_1),5))
    for(j in 1:(5)){
      if(length(idc_1)!=1){
      linear_pred = a_coefficients[,,j] %*% t(X_matrices[idc_1,tt,1:29]) + 
        a_coefficients_t[,,j] %*% t(X_matrices[idc_1,tt-1,30:33]) + 
        as.vector(d[,tt,j]) + g[,df_numeric[idc_1,tt],j] + 
        e[,idc_1,j] + as.vector(b0[,j])*tt  +
        as.vector(a0[,j]) +psi[,PSU_trans[idc_1,tt],j] + 
        k[,STR_trans[idc_1,tt],j]
      exp_softmax[,,j] = exp(linear_pred)
      }else{
        linear_pred = a_coefficients[,,j] %*% (X_matrices[idc_1,tt,1:29]) + 
          a_coefficients_t[,,j] %*% (X_matrices[idc_1,tt-1,30:33]) + 
          as.vector(d[,tt,j]) + g[,df_numeric[idc_1,tt],j] + 
          e[,idc_1,j] + as.vector(b0[,j])*tt  +
          as.vector(a0[,j]) +psi[,PSU_trans[idc_1,tt],j] + 
          k[,STR_trans[idc_1,tt],j]
        exp_softmax[,,j] = exp(linear_pred)
      }
    }
    est_prob = array(NA,c(1000,length(idc_1),5))
    for(i in 1:length(idc_1)){
      for(j in 1:5){
        est_prob[,i,j] =  exp_softmax[,i,j]/(exp_softmax[,i,1] + exp_softmax[,i,2] +
                                               exp_softmax[,i,3] + exp_softmax[,i,4] +
                                               exp_softmax[,i,5])
      }
    }
    No_all[,tt,] = apply(est_prob,c(1,3),mean)
    message(tt)
  }
  saveRDS(No_all, file = "AR_No_all_years.RDS")
  
  
}

plot(apply(UC_all[,,5],2,mean),ylim=c(0,1))
points(apply(LW_all[,,5],2,mean),col='red')
points(apply(No_all[,,5],2,mean),col='green')

#gender FtoF
{
  AA = array(NA,c(1000,12,5))
  for(tt in 2:12){
    idc_1 = treat
    idc_1 = idc_1[which(apply(SEX[idc_1,],1,mean,na.rm=TRUE) == 1 & GHQ_SFS[idc_1,tt-1] == 1 & idc_1 %in% T_ind[T_ind[,2]==tt,1])] 
    exp_softmax = array(NA,c(1000,length(idc_1),5))
    for(j in 1:(5)){
      if(length(idc_1)!=1){
        linear_pred = a_coefficients[,,j] %*% t(X_matrices[idc_1,tt,1:29]) + 
          a_coefficients_t[,,j] %*% t(X_matrices[idc_1,tt-1,30:33]) + 
          as.vector(d[,tt,j]) + g[,df_numeric[idc_1,tt],j] + 
          e[,idc_1,j] + as.vector(b0[,j])*tt  +
          as.vector(a0[,j]) +psi[,PSU_trans[idc_1,tt],j] + 
          k[,STR_trans[idc_1,tt],j]
        exp_softmax[,,j] = exp(linear_pred)
      }else{
        linear_pred = a_coefficients[,,j] %*% (X_matrices[idc_1,tt,1:29]) + 
          a_coefficients_t[,,j] %*% (X_matrices[idc_1,tt-1,30:33]) + 
          as.vector(d[,tt,j]) + g[,df_numeric[idc_1,tt],j] + 
          e[,idc_1,j] + as.vector(b0[,j])*tt  +
          as.vector(a0[,j]) +psi[,PSU_trans[idc_1,tt],j] + 
          k[,STR_trans[idc_1,tt],j]
        exp_softmax[,,j] = exp(linear_pred)
      }
    }
    est_prob = array(NA,c(1000,length(idc_1),5))
    for(i in 1:length(idc_1)){
      for(j in 1:5){
        est_prob[,i,j] =  exp_softmax[,i,j]/(exp_softmax[,i,1] + exp_softmax[,i,2] +
                                               exp_softmax[,i,3] + exp_softmax[,i,4] +
                                               exp_softmax[,i,5])
      }
    }
    AA[,tt,] = apply(est_prob,c(1,3),mean)
    message(tt)
  }
  
  AA1 = array(NA,c(1000,12,5))
  for(tt in 2:12){
    idc_1 = treat
    idc_1 = idc_1[which(apply(SEX[idc_1,],1,mean,na.rm=TRUE) == 2 &GHQ_SFS[idc_1,tt-1] == 1& idc_1 %in% T_ind[T_ind[,2]==tt,1])] 
    exp_softmax = array(NA,c(1000,length(idc_1),5))

    for(j in 1:(5)){
     if(length(idc_1)!=1){
        linear_pred = a_coefficients[,,j] %*% t(X_matrices[idc_1,tt,1:29]) + 
          a_coefficients_t[,,j] %*% t(X_matrices[idc_1,tt-1,30:33]) + 
          as.vector(d[,tt,j]) + g[,df_numeric[idc_1,tt],j] + 
          e[,idc_1,j] + as.vector(b0[,j])*tt  +
          as.vector(a0[,j]) +psi[,PSU_trans[idc_1,tt],j] + 
          k[,STR_trans[idc_1,tt],j]
        exp_softmax[,,j] = exp(linear_pred)
     }else{
        linear_pred = a_coefficients[,,j] %*% (X_matrices[idc_1,tt,1:29]) + 
          a_coefficients_t[,,j] %*% (X_matrices[idc_1,tt-1,30:33]) + 
          as.vector(d[,tt,j]) + g[,df_numeric[idc_1,tt],j] + 
          e[,idc_1,j] + as.vector(b0[,j])*tt  +
          as.vector(a0[,j]) +psi[,PSU_trans[idc_1,tt],j] + 
          k[,STR_trans[idc_1,tt],j]
        exp_softmax[,,j] = exp(linear_pred)
     }
    }
    est_prob = array(NA,c(1000,length(idc_1),5))
    for(i in 1:length(idc_1)){
      for(j in 1:5){
        est_prob[,i,j] =  exp_softmax[,i,j]/(exp_softmax[,i,1] + exp_softmax[,i,2] +
                                               exp_softmax[,i,3] + exp_softmax[,i,4] +
                                               exp_softmax[,i,5])
      }
    }
    AA1[,tt,] = apply(est_prob,c(1,3),mean)
    message(tt)
  }
  
  saveRDS(AA, file = "AR_UC_all_years_M.RDS")
  saveRDS(AA1, file = "AR_UC_all_years_F.RDS")
  
  BB = array(NA,c(1000,12,5))
  for(tt in 2:12){
    idc_1 = contr
    idc_1 = idc_1[which(apply(SEX[idc_1,],1,mean,na.rm=TRUE) == 1 & GHQ_SFS[idc_1,tt-1] == 1 & idc_1 %in% T_ind[T_ind[,2]==tt,1])] 
    exp_softmax = array(NA,c(1000,length(idc_1),5))
    for(j in 1:(5)){
      if(length(idc_1)!=1){
        linear_pred = a_coefficients[,,j] %*% t(X_matrices[idc_1,tt,1:29]) + 
          a_coefficients_t[,,j] %*% t(X_matrices[idc_1,tt-1,30:33]) + 
          as.vector(d[,tt,j]) + g[,df_numeric[idc_1,tt],j] + 
          e[,idc_1,j] + as.vector(b0[,j])*tt  +
          as.vector(a0[,j]) +psi[,PSU_trans[idc_1,tt],j] + 
          k[,STR_trans[idc_1,tt],j]
        exp_softmax[,,j] = exp(linear_pred)
      }else{
        linear_pred = a_coefficients[,,j] %*% (X_matrices[idc_1,tt,1:29]) + 
          a_coefficients_t[,,j] %*% (X_matrices[idc_1,tt-1,30:33]) + 
          as.vector(d[,tt,j]) + g[,df_numeric[idc_1,tt],j] + 
          e[,idc_1,j] + as.vector(b0[,j])*tt  +
          as.vector(a0[,j]) +psi[,PSU_trans[idc_1,tt],j] + 
          k[,STR_trans[idc_1,tt],j]
        exp_softmax[,,j] = exp(linear_pred)
      }
    }
    est_prob = array(NA,c(1000,length(idc_1),5))
    for(i in 1:length(idc_1)){
      for(j in 1:5){
        est_prob[,i,j] =  exp_softmax[,i,j]/(exp_softmax[,i,1] + exp_softmax[,i,2] +
                                               exp_softmax[,i,3] + exp_softmax[,i,4] +
                                               exp_softmax[,i,5])
      }
    }
    BB[,tt,] = apply(est_prob,c(1,3),mean)
    message(tt)
  }
  
  BB1 = array(NA,c(1000,12,5))
  for(tt in 2:12){
    idc_1 = contr
    idc_1 = idc_1[which(apply(SEX[idc_1,],1,mean,na.rm=TRUE) == 2 & GHQ_SFS[idc_1,tt-1] == 1 & idc_1 %in% T_ind[T_ind[,2]==tt,1])] 
    exp_softmax = array(NA,c(1000,length(idc_1),5))
    
    for(j in 1:(5)){
      if(length(idc_1)!=1){
        linear_pred = a_coefficients[,,j] %*% t(X_matrices[idc_1,tt,1:29]) + 
          a_coefficients_t[,,j] %*% t(X_matrices[idc_1,tt-1,30:33]) + 
          as.vector(d[,tt,j]) + g[,df_numeric[idc_1,tt],j] + 
          e[,idc_1,j] + as.vector(b0[,j])*tt  +
          as.vector(a0[,j]) +psi[,PSU_trans[idc_1,tt],j] + 
          k[,STR_trans[idc_1,tt],j]
        exp_softmax[,,j] = exp(linear_pred)
      }else{
        linear_pred = a_coefficients[,,j] %*% (X_matrices[idc_1,tt,1:29]) + 
          a_coefficients_t[,,j] %*% (X_matrices[idc_1,tt-1,30:33]) + 
          as.vector(d[,tt,j]) + g[,df_numeric[idc_1,tt],j] + 
          e[,idc_1,j] + as.vector(b0[,j])*tt  +
          as.vector(a0[,j]) +psi[,PSU_trans[idc_1,tt],j] + 
          k[,STR_trans[idc_1,tt],j]
        exp_softmax[,,j] = exp(linear_pred)
      }
    }
    est_prob = array(NA,c(1000,length(idc_1),5))
    for(i in 1:length(idc_1)){
      for(j in 1:5){
        est_prob[,i,j] =  exp_softmax[,i,j]/(exp_softmax[,i,1] + exp_softmax[,i,2] +
                                               exp_softmax[,i,3] + exp_softmax[,i,4] +
                                               exp_softmax[,i,5])
      }
    }
    BB1[,tt,] = apply(est_prob,c(1,3),mean)
    message(tt)
  }
  saveRDS(BB, file = "AR_LW_all_years_M.RDS")
  saveRDS(BB1, file = "AR_LW_all_years_F.RDS")
  
  CC = array(NA,c(1000,12,5))
  for(tt in 2:12){
    idc_1 = contr_2
    idc_1 = idc_1[which(apply(SEX[idc_1,],1,mean,na.rm=TRUE) == 1 & GHQ_SFS[idc_1,tt-1] == 1 & idc_1 %in% T_ind[T_ind[,2]==tt,1])] 
    exp_softmax = array(NA,c(1000,length(idc_1),5))
    
    for(j in 1:(5)){
      if(length(idc_1)!=1){
        linear_pred = a_coefficients[,,j] %*% t(X_matrices[idc_1,tt,1:29]) + 
          a_coefficients_t[,,j] %*% t(X_matrices[idc_1,tt-1,30:33]) + 
          as.vector(d[,tt,j]) + g[,df_numeric[idc_1,tt],j] + 
          e[,idc_1,j] + as.vector(b0[,j])*tt  +
          as.vector(a0[,j]) +psi[,PSU_trans[idc_1,tt],j] + 
          k[,STR_trans[idc_1,tt],j]
        exp_softmax[,,j] = exp(linear_pred)
      }else{
        linear_pred = a_coefficients[,,j] %*% (X_matrices[idc_1,tt,1:29]) + 
          a_coefficients_t[,,j] %*% (X_matrices[idc_1,tt-1,30:33]) + 
          as.vector(d[,tt,j]) + g[,df_numeric[idc_1,tt],j] + 
          e[,idc_1,j] + as.vector(b0[,j])*tt  +
          as.vector(a0[,j]) +psi[,PSU_trans[idc_1,tt],j] + 
          k[,STR_trans[idc_1,tt],j]
        exp_softmax[,,j] = exp(linear_pred)
      }
    }
    est_prob = array(NA,c(1000,length(idc_1),5))
    for(i in 1:length(idc_1)){
      for(j in 1:5){
        est_prob[,i,j] =  exp_softmax[,i,j]/(exp_softmax[,i,1] + exp_softmax[,i,2] +
                                               exp_softmax[,i,3] + exp_softmax[,i,4] +
                                               exp_softmax[,i,5])
      }
    }
    CC[,tt,] = apply(est_prob,c(1,3),mean)
    message(tt)
  }
  
  CC1 = array(NA,c(1000,12,5))
  for(tt in 2:12){
    idc_1 = contr_2
    idc_1 = idc_1[which(apply(SEX[idc_1,],1,mean,na.rm=TRUE) == 2 &GHQ_SFS[idc_1,tt-1] == 1 & idc_1 %in% T_ind[T_ind[,2]==tt,1])] 
    exp_softmax = array(NA,c(1000,length(idc_1),5))
    
    for(j in 1:(5)){
      if(length(idc_1)!=1){
        linear_pred = a_coefficients[,,j] %*% t(X_matrices[idc_1,tt,1:29]) + 
          a_coefficients_t[,,j] %*% t(X_matrices[idc_1,tt-1,30:33]) + 
          as.vector(d[,tt,j]) + g[,df_numeric[idc_1,tt],j] + 
          e[,idc_1,j] + as.vector(b0[,j])*tt  +
          as.vector(a0[,j]) +psi[,PSU_trans[idc_1,tt],j] + 
          k[,STR_trans[idc_1,tt],j]
        exp_softmax[,,j] = exp(linear_pred)
      }else{
        linear_pred = a_coefficients[,,j] %*% (X_matrices[idc_1,tt,1:29]) + 
          a_coefficients_t[,,j] %*% (X_matrices[idc_1,tt-1,30:33]) + 
          as.vector(d[,tt,j]) + g[,df_numeric[idc_1,tt],j] + 
          e[,idc_1,j] + as.vector(b0[,j])*tt  +
          as.vector(a0[,j]) +psi[,PSU_trans[idc_1,tt],j] + 
          k[,STR_trans[idc_1,tt],j]
        exp_softmax[,,j] = exp(linear_pred)
      }
    }
    est_prob = array(NA,c(1000,length(idc_1),5))
    for(i in 1:length(idc_1)){
      for(j in 1:5){
        est_prob[,i,j] =  exp_softmax[,i,j]/(exp_softmax[,i,1] + exp_softmax[,i,2] +
                                               exp_softmax[,i,3] + exp_softmax[,i,4] +
                                               exp_softmax[,i,5])
      }
    }
    CC1[,tt,] = apply(est_prob,c(1,3),mean)
    message(tt)
  }
  saveRDS(CC, file = "AR_No_all_years_M.RDS")
  saveRDS(CC1, file = "AR_No_all_years_F.RDS")  
}

#ethn FtoF
{
  
  
  DD = array(NA,c(1000,12,5))
  for(tt in 2:12){
    idc_1 = treat
    idc_1 = idc_1[which(apply(ETH[idc_1,],1,mean,na.rm=TRUE) == 1 & GHQ_SFS[idc_1,tt-1] == 1& idc_1 %in% T_ind[T_ind[,2]==tt,1])] 
    exp_softmax = array(NA,c(1000,length(idc_1),5))
    
    for(j in 1:(5)){
     if(length(idc_1)!=1){
        linear_pred = a_coefficients[,,j] %*% t(X_matrices[idc_1,tt,1:29]) + 
          a_coefficients_t[,,j] %*% t(X_matrices[idc_1,tt-1,30:33]) + 
          as.vector(d[,tt,j]) + g[,df_numeric[idc_1,tt],j] + 
          e[,idc_1,j] + as.vector(b0[,j])*tt  +
          as.vector(a0[,j]) +psi[,PSU_trans[idc_1,tt],j] + 
          k[,STR_trans[idc_1,tt],j]
        exp_softmax[,,j] = exp(linear_pred)
     }else{
        linear_pred = a_coefficients[,,j] %*% (X_matrices[idc_1,tt,1:29]) + 
         a_coefficients_t[,,j] %*% (X_matrices[idc_1,tt-1,30:33]) + 
          as.vector(d[,tt,j]) + g[,df_numeric[idc_1,tt],j] + 
          e[,idc_1,j] + as.vector(b0[,j])*tt  +
          as.vector(a0[,j]) +psi[,PSU_trans[idc_1,tt],j] + 
          k[,STR_trans[idc_1,tt],j]
        exp_softmax[,,j] = exp(linear_pred)
     }
    }
    est_prob = array(NA,c(1000,length(idc_1),5))
    for(i in 1:length(idc_1)){
      for(j in 1:5){
        est_prob[,i,j] =  exp_softmax[,i,j]/(exp_softmax[,i,1] + exp_softmax[,i,2] +
                                               exp_softmax[,i,3] + exp_softmax[,i,4] +
                                               exp_softmax[,i,5])
      }
    }
    DD[,tt,] = apply(est_prob,c(1,3),mean)
    message(tt)
  }
  
  DD1 = array(NA,c(1000,12,5))
  for(tt in 2:12){
    idc_1 = treat
    idc_1 = idc_1[which((apply(ETH[idc_1,],1,mean,na.rm=TRUE) == 3 |
                           apply(ETH[idc_1,],1,mean,na.rm=TRUE) == 4 ) & GHQ_SFS[idc_1,tt-1] == 1& idc_1 %in% T_ind[T_ind[,2]==tt,1])] 
    exp_softmax = array(NA,c(1000,length(idc_1),5))
    
    for(j in 1:(5)){
      if(length(idc_1)!=1){
        linear_pred = a_coefficients[,,j] %*% t(X_matrices[idc_1,tt,1:29]) + 
          a_coefficients_t[,,j] %*% t(X_matrices[idc_1,tt-1,30:33]) + 
          as.vector(d[,tt,j]) + g[,df_numeric[idc_1,tt],j] + 
          e[,idc_1,j] + as.vector(b0[,j])*tt  +
          as.vector(a0[,j]) +psi[,PSU_trans[idc_1,tt],j] + 
          k[,STR_trans[idc_1,tt],j]
        exp_softmax[,,j] = exp(linear_pred)
      }else{
        linear_pred = a_coefficients[,,j] %*% (X_matrices[idc_1,tt,1:29]) + 
          a_coefficients_t[,,j] %*% (X_matrices[idc_1,tt-1,30:33]) + 
          as.vector(d[,tt,j]) + g[,df_numeric[idc_1,tt],j] + 
          e[,idc_1,j] + as.vector(b0[,j])*tt  +
          as.vector(a0[,j]) +psi[,PSU_trans[idc_1,tt],j] + 
          k[,STR_trans[idc_1,tt],j]
        exp_softmax[,,j] = exp(linear_pred)
      }
    }
    est_prob = array(NA,c(1000,length(idc_1),5))
    for(i in 1:length(idc_1)){
      for(j in 1:5){
        est_prob[,i,j] =  exp_softmax[,i,j]/(exp_softmax[,i,1] + exp_softmax[,i,2] +
                                               exp_softmax[,i,3] + exp_softmax[,i,4] +
                                               exp_softmax[,i,5])
      }
    }
    DD1[,tt,] = apply(est_prob,c(1,3),mean)
    message(tt)
  }
  
  saveRDS(DD, file = "AR_UC_all_years_White.RDS")
  saveRDS(DD1, file = "AR_UC_all_years_AsBl.RDS")
  
  FF = array(NA,c(1000,12,5))
  for(tt in 2:12){
    idc_1 = contr
    idc_1 = idc_1[which(apply(ETH[idc_1,],1,mean,na.rm=TRUE) == 1 & GHQ_SFS[idc_1,tt-1] == 1 & idc_1 %in% T_ind[T_ind[,2]==tt,1])] 
    exp_softmax = array(NA,c(1000,length(idc_1),5))
    
    for(j in 1:(5)){
     if(length(idc_1)!=1){
        linear_pred = a_coefficients[,,j] %*% t(X_matrices[idc_1,tt,1:29]) + 
          a_coefficients_t[,,j] %*% t(X_matrices[idc_1,tt-1,30:33]) + 
          as.vector(d[,tt,j]) + g[,df_numeric[idc_1,tt],j] + 
          e[,idc_1,j] + as.vector(b0[,j])*tt  +
          as.vector(a0[,j]) +psi[,PSU_trans[idc_1,tt],j] + 
          k[,STR_trans[idc_1,tt],j]
        exp_softmax[,,j] = exp(linear_pred)
     }else{
        linear_pred = a_coefficients[,,j] %*% (X_matrices[idc_1,tt,1:29]) + 
          a_coefficients_t[,,j] %*% (X_matrices[idc_1,tt-1,30:33]) + 
          as.vector(d[,tt,j]) + g[,df_numeric[idc_1,tt],j] + 
          e[,idc_1,j] + as.vector(b0[,j])*tt  +
          as.vector(a0[,j]) +psi[,PSU_trans[idc_1,tt],j] + 
          k[,STR_trans[idc_1,tt],j]
        exp_softmax[,,j] = exp(linear_pred)
     }
    }
    est_prob = array(NA,c(1000,length(idc_1),5))
    for(i in 1:length(idc_1)){
      for(j in 1:5){
        est_prob[,i,j] =  exp_softmax[,i,j]/(exp_softmax[,i,1] + exp_softmax[,i,2] +
                                               exp_softmax[,i,3] + exp_softmax[,i,4] +
                                               exp_softmax[,i,5])
      }
    }
    FF[,tt,] = apply(est_prob,c(1,3),mean)
    message(tt)
  }
  
 
  FF1 = array(NA,c(1000,12,5))
  for(tt in 2:12){
    idc_1 = contr
    idc_1 = idc_1[which((apply(ETH[idc_1,],1,mean,na.rm=TRUE) == 3 |
                           apply(ETH[idc_1,],1,mean,na.rm=TRUE) == 4 ) & GHQ_SFS[idc_1,tt-1] == 1 & idc_1 %in% T_ind[T_ind[,2]==tt,1])] 
    exp_softmax = array(NA,c(1000,length(idc_1),5))
    
    for(j in 1:(5)){
      if(length(idc_1)!=1){
        linear_pred = a_coefficients[,,j] %*% t(X_matrices[idc_1,tt,1:29]) + 
          a_coefficients_t[,,j] %*% t(X_matrices[idc_1,tt-1,30:33]) + 
          as.vector(d[,tt,j]) + g[,df_numeric[idc_1,tt],j] + 
          e[,idc_1,j] + as.vector(b0[,j])*tt  +
          as.vector(a0[,j]) +psi[,PSU_trans[idc_1,tt],j] + 
          k[,STR_trans[idc_1,tt],j]
        exp_softmax[,,j] = exp(linear_pred)
      }else{
        linear_pred = a_coefficients[,,j] %*% (X_matrices[idc_1,tt,1:29]) + 
          a_coefficients_t[,,j] %*% (X_matrices[idc_1,tt-1,30:33]) + 
          as.vector(d[,tt,j]) + g[,df_numeric[idc_1,tt],j] + 
          e[,idc_1,j] + as.vector(b0[,j])*tt  +
          as.vector(a0[,j]) +psi[,PSU_trans[idc_1,tt],j] + 
          k[,STR_trans[idc_1,tt],j]
        exp_softmax[,,j] = exp(linear_pred)
      }
    }
    est_prob = array(NA,c(1000,length(idc_1),5))
    for(i in 1:length(idc_1)){
      for(j in 1:5){
        est_prob[,i,j] =  exp_softmax[,i,j]/(exp_softmax[,i,1] + exp_softmax[,i,2] +
                                               exp_softmax[,i,3] + exp_softmax[,i,4] +
                                               exp_softmax[,i,5])
      }
    }
    FF1[,tt,] = apply(est_prob,c(1,3),mean)
    message(tt)
  }
  
  saveRDS(FF, file = "AR_LW_all_years_White.RDS")
  saveRDS(FF1, file = "AR_LW_all_years_AsBl.RDS")
  
  GG = array(NA,c(1000,12,5))
  for(tt in 2:12){
    idc_1 = contr_2
    idc_1 = idc_1[which(apply(ETH[idc_1,],1,mean,na.rm=TRUE) == 1 & GHQ_SFS[idc_1,tt-1] == 1 & idc_1 %in% T_ind[T_ind[,2]==tt,1])] 
    exp_softmax = array(NA,c(1000,length(idc_1),5))
    
    for(j in 1:(5)){
     if(length(idc_1)!=1){
        linear_pred = a_coefficients[,,j] %*% t(X_matrices[idc_1,tt,1:29]) + 
         a_coefficients_t[,,j] %*% t(X_matrices[idc_1,tt-1,30:33]) + 
          as.vector(d[,tt,j]) + g[,df_numeric[idc_1,tt],j] + 
          e[,idc_1,j] + as.vector(b0[,j])*tt  +
          as.vector(a0[,j]) +psi[,PSU_trans[idc_1,tt],j] + 
          k[,STR_trans[idc_1,tt],j]
        exp_softmax[,,j] = exp(linear_pred)
     }else{
        linear_pred = a_coefficients[,,j] %*% (X_matrices[idc_1,tt,1:29]) + 
          a_coefficients_t[,,j] %*% (X_matrices[idc_1,tt-1,30:33]) + 
          as.vector(d[,tt,j]) + g[,df_numeric[idc_1,tt],j] + 
          e[,idc_1,j] + as.vector(b0[,j])*tt  +
          as.vector(a0[,j]) +psi[,PSU_trans[idc_1,tt],j] + 
          k[,STR_trans[idc_1,tt],j]
        exp_softmax[,,j] = exp(linear_pred)
     }
    }
    est_prob = array(NA,c(1000,length(idc_1),5))
    for(i in 1:length(idc_1)){
      for(j in 1:5){
        est_prob[,i,j] =  exp_softmax[,i,j]/(exp_softmax[,i,1] + exp_softmax[,i,2] +
                                               exp_softmax[,i,3] + exp_softmax[,i,4] +
                                               exp_softmax[,i,5])
      }
    }
    GG[,tt,] = apply(est_prob,c(1,3),mean)
    message(tt)
  }
  
  
  
  GG1 = array(NA,c(1000,12,5))
  for(tt in 2:12){
    idc_1 = contr_2
    idc_1 = idc_1[which((apply(ETH[idc_1,],1,mean,na.rm=TRUE) == 3 |
                           apply(ETH[idc_1,],1,mean,na.rm=TRUE) == 4 ) & GHQ_SFS[idc_1,tt-1] == 1 & idc_1 %in% T_ind[T_ind[,2]==tt,1])] 
    exp_softmax = array(NA,c(1000,length(idc_1),5))
    
    for(j in 1:(5)){
    if(length(idc_1)!=1){
        linear_pred = a_coefficients[,,j] %*% t(X_matrices[idc_1,tt,1:29]) + 
          a_coefficients_t[,,j] %*% t(X_matrices[idc_1,tt-1,30:33]) + 
          as.vector(d[,tt,j]) + g[,df_numeric[idc_1,tt],j] + 
          e[,idc_1,j] + as.vector(b0[,j])*tt  +
          as.vector(a0[,j]) +psi[,PSU_trans[idc_1,tt],j] + 
          k[,STR_trans[idc_1,tt],j]
        exp_softmax[,,j] = exp(linear_pred)
    }else{
        linear_pred = a_coefficients[,,j] %*% (X_matrices[idc_1,tt,1:29]) + 
          a_coefficients_t[,,j] %*% (X_matrices[idc_1,tt-1,30:33]) + 
          as.vector(d[,tt,j]) + g[,df_numeric[idc_1,tt],j] + 
          e[,idc_1,j] + as.vector(b0[,j])*tt  +
          as.vector(a0[,j]) +psi[,PSU_trans[idc_1,tt],j] + 
          k[,STR_trans[idc_1,tt],j]
        exp_softmax[,,j] = exp(linear_pred)
    }
    }
    est_prob = array(NA,c(1000,length(idc_1),5))
    for(i in 1:length(idc_1)){
      for(j in 1:5){
        est_prob[,i,j] =  exp_softmax[,i,j]/(exp_softmax[,i,1] + exp_softmax[,i,2] +
                                               exp_softmax[,i,3] + exp_softmax[,i,4] +
                                               exp_softmax[,i,5])
      }
    }
    GG1[,tt,] = apply(est_prob,c(1,3),mean)
    message(tt)
  }

  saveRDS(GG, file = "AR_No_all_years_White.RDS")
  saveRDS(GG1, file = "AR_No_all_years_AsBl.RDS")  
}


#MT FtoF
{
  
  V1 = array(NA,c(1000,12,5))
  for(tt in 2:12){
    idc_1 = treat
    idc_1 = idc_1[which(MT[idc_1,tt] == 1  & GHQ_SFS[idc_1,tt-1] == 1 & idc_1 %in% T_ind[T_ind[,2]==tt,1])] 
    exp_softmax = array(NA,c(1000,length(idc_1),5))
    
    for(j in 1:(5)){
     if(length(idc_1)!=1){
        linear_pred = a_coefficients[,,j] %*% t(X_matrices[idc_1,tt,1:29]) + 
          a_coefficients_t[,,j] %*% t(X_matrices[idc_1,tt-1,30:33]) + 
          as.vector(d[,tt,j]) + g[,df_numeric[idc_1,tt],j] + 
          e[,idc_1,j] + as.vector(b0[,j])*tt  +
          as.vector(a0[,j]) +psi[,PSU_trans[idc_1,tt],j] + 
          k[,STR_trans[idc_1,tt],j]
        exp_softmax[,,j] = exp(linear_pred)
     }else{
        linear_pred = a_coefficients[,,j] %*% (X_matrices[idc_1,tt,1:29]) + 
          a_coefficients_t[,,j] %*% (X_matrices[idc_1,tt-1,30:33]) + 
          as.vector(d[,tt,j]) + g[,df_numeric[idc_1,tt],j] + 
          e[,idc_1,j] + as.vector(b0[,j])*tt  +
          as.vector(a0[,j]) +psi[,PSU_trans[idc_1,tt],j] + 
          k[,STR_trans[idc_1,tt],j]
        exp_softmax[,,j] = exp(linear_pred)
     }
    }
    est_prob = array(NA,c(1000,length(idc_1),5))
    for(i in 1:length(idc_1)){
      for(j in 1:5){
        est_prob[,i,j] =  exp_softmax[,i,j]/(exp_softmax[,i,1] + exp_softmax[,i,2] +
                                               exp_softmax[,i,3] + exp_softmax[,i,4] +
                                               exp_softmax[,i,5])
      }
    }
    V1[,tt,] = apply(est_prob,c(1,3),mean)
    message(tt)
  }
  
  }

  V2 = array(NA,c(1000,12,5))
  for(tt in 2:12){
    idc_1 = treat
    idc_1 = idc_1[which(MT[idc_1,tt] == 2  & GHQ_SFS[idc_1,tt-1] == 1 & idc_1 %in% T_ind[T_ind[,2]==tt,1])] 
    exp_softmax = array(NA,c(1000,length(idc_1),5))
    
    for(j in 1:(5)){
      if(length(idc_1)!=1){
        linear_pred = a_coefficients[,,j] %*% t(X_matrices[idc_1,tt,1:29]) + 
          a_coefficients_t[,,j] %*% t(X_matrices[idc_1,tt-1,30:33]) + 
          as.vector(d[,tt,j]) + g[,df_numeric[idc_1,tt],j] + 
          e[,idc_1,j] + as.vector(b0[,j])*tt  +
          as.vector(a0[,j]) +psi[,PSU_trans[idc_1,tt],j] + 
          k[,STR_trans[idc_1,tt],j]
        exp_softmax[,,j] = exp(linear_pred)
      }else{
        linear_pred = a_coefficients[,,j] %*% (X_matrices[idc_1,tt,1:29]) + 
          a_coefficients_t[,,j] %*% (X_matrices[idc_1,tt-1,30:33]) + 
          as.vector(d[,tt,j]) + g[,df_numeric[idc_1,tt],j] + 
          e[,idc_1,j] + as.vector(b0[,j])*tt  +
          as.vector(a0[,j]) +psi[,PSU_trans[idc_1,tt],j] + 
          k[,STR_trans[idc_1,tt],j]
        exp_softmax[,,j] = exp(linear_pred)
      }
    }
    est_prob = array(NA,c(1000,length(idc_1),5))
    for(i in 1:length(idc_1)){
      for(j in 1:5){
        est_prob[,i,j] =  exp_softmax[,i,j]/(exp_softmax[,i,1] + exp_softmax[,i,2] +
                                               exp_softmax[,i,3] + exp_softmax[,i,4] +
                                               exp_softmax[,i,5])
      }
    }
    V2[,tt,] = apply(est_prob,c(1,3),mean)
    message(tt)
  }
 
  saveRDS(V1, file = "AR_UC_all_years_Sin.RDS")
  saveRDS(V2, file = "AR_UC_all_years_Rel.RDS")
  
  V3 = array(NA,c(1000,12,5))
  for(tt in 2:12){
    idc_1 = contr
    idc_1 = idc_1[which(MT[idc_1,tt] == 1  & GHQ_SFS[idc_1,tt-1] == 1 & idc_1 %in% T_ind[T_ind[,2]==tt,1])] 
    exp_softmax = array(NA,c(1000,length(idc_1),5))
    
    for(j in 1:(5)){
     if(length(idc_1)!=1){
        linear_pred = a_coefficients[,,j] %*% t(X_matrices[idc_1,tt,1:29]) + 
          a_coefficients_t[,,j] %*% t(X_matrices[idc_1,tt-1,30:33]) + 
          as.vector(d[,tt,j]) + g[,df_numeric[idc_1,tt],j] + 
          e[,idc_1,j] + as.vector(b0[,j])*tt  +
          as.vector(a0[,j]) +psi[,PSU_trans[idc_1,tt],j] + 
          k[,STR_trans[idc_1,tt],j]
        exp_softmax[,,j] = exp(linear_pred)
     }else{
        linear_pred = a_coefficients[,,j] %*% (X_matrices[idc_1,tt,1:29]) + 
          a_coefficients_t[,,j] %*% (X_matrices[idc_1,tt-1,30:33]) + 
          as.vector(d[,tt,j]) + g[,df_numeric[idc_1,tt],j] + 
          e[,idc_1,j] + as.vector(b0[,j])*tt  +
          as.vector(a0[,j]) +psi[,PSU_trans[idc_1,tt],j] + 
          k[,STR_trans[idc_1,tt],j]
        exp_softmax[,,j] = exp(linear_pred)
     }
    }
    est_prob = array(NA,c(1000,length(idc_1),5))
    for(i in 1:length(idc_1)){
      for(j in 1:5){
        est_prob[,i,j] =  exp_softmax[,i,j]/(exp_softmax[,i,1] + exp_softmax[,i,2] +
                                               exp_softmax[,i,3] + exp_softmax[,i,4] +
                                               exp_softmax[,i,5])
      }
    }
    V3[,tt,] = apply(est_prob,c(1,3),mean)
    message(tt)
  }
  
 
  V4 = array(NA,c(1000,12,5))
  for(tt in 2:12){
    idc_1 = contr
    idc_1 = idc_1[which(MT[idc_1,tt] == 2  & GHQ_SFS[idc_1,tt-1] == 1 & idc_1 %in% T_ind[T_ind[,2]==tt,1])] 
    exp_softmax = array(NA,c(1000,length(idc_1),5))
    
    for(j in 1:(5)){
      if(length(idc_1)!=1){
        linear_pred = a_coefficients[,,j] %*% t(X_matrices[idc_1,tt,1:29]) + 
          a_coefficients_t[,,j] %*% t(X_matrices[idc_1,tt-1,30:33]) + 
          as.vector(d[,tt,j]) + g[,df_numeric[idc_1,tt],j] + 
          e[,idc_1,j] + as.vector(b0[,j])*tt  +
          as.vector(a0[,j]) +psi[,PSU_trans[idc_1,tt],j] + 
          k[,STR_trans[idc_1,tt],j]
        exp_softmax[,,j] = exp(linear_pred)
      }else{
        linear_pred = a_coefficients[,,j] %*% (X_matrices[idc_1,tt,1:29]) + 
          a_coefficients_t[,,j] %*% (X_matrices[idc_1,tt-1,30:33]) + 
          as.vector(d[,tt,j]) + g[,df_numeric[idc_1,tt],j] + 
          e[,idc_1,j] + as.vector(b0[,j])*tt  +
          as.vector(a0[,j]) +psi[,PSU_trans[idc_1,tt],j] + 
          k[,STR_trans[idc_1,tt],j]
        exp_softmax[,,j] = exp(linear_pred)
      }
    }
    est_prob = array(NA,c(1000,length(idc_1),5))
    for(i in 1:length(idc_1)){
      for(j in 1:5){
        est_prob[,i,j] =  exp_softmax[,i,j]/(exp_softmax[,i,1] + exp_softmax[,i,2] +
                                               exp_softmax[,i,3] + exp_softmax[,i,4] +
                                               exp_softmax[,i,5])
      }
    }
    V4[,tt,] = apply(est_prob,c(1,3),mean)
    message(tt)
  }
 
  saveRDS(V3, file = "AR_LW_all_years_Sin.RDS")
  saveRDS(V4, file = "AR_LW_all_years_Rel.RDS")
  
  V5 = array(NA,c(1000,12,5))
  for(tt in 2:12){
    idc_1 = contr_2
    idc_1 = idc_1[which(MT[idc_1,tt] == 1  & GHQ_SFS[idc_1,tt-1] == 1 & idc_1 %in% T_ind[T_ind[,2]==tt,1])] 
    exp_softmax = array(NA,c(1000,length(idc_1),5))
    
    for(j in 1:(5)){
      if(length(idc_1)!=1){
        linear_pred = a_coefficients[,,j] %*% t(X_matrices[idc_1,tt,1:29]) + 
          a_coefficients_t[,,j] %*% t(X_matrices[idc_1,tt-1,30:33]) + 
          as.vector(d[,tt,j]) + g[,df_numeric[idc_1,tt],j] + 
          e[,idc_1,j] + as.vector(b0[,j])*tt  +
          as.vector(a0[,j]) +psi[,PSU_trans[idc_1,tt],j] + 
          k[,STR_trans[idc_1,tt],j]
        exp_softmax[,,j] = exp(linear_pred)
      }else{
        linear_pred = a_coefficients[,,j] %*% (X_matrices[idc_1,tt,1:29]) + 
          a_coefficients_t[,,j] %*% (X_matrices[idc_1,tt-1,30:33]) + 
          as.vector(d[,tt,j]) + g[,df_numeric[idc_1,tt],j] + 
          e[,idc_1,j] + as.vector(b0[,j])*tt  +
          as.vector(a0[,j]) +psi[,PSU_trans[idc_1,tt],j] + 
          k[,STR_trans[idc_1,tt],j]
        exp_softmax[,,j] = exp(linear_pred)
      }
    }
    est_prob = array(NA,c(1000,length(idc_1),5))
    for(i in 1:length(idc_1)){
      for(j in 1:5){
        est_prob[,i,j] =  exp_softmax[,i,j]/(exp_softmax[,i,1] + exp_softmax[,i,2] +
                                               exp_softmax[,i,3] + exp_softmax[,i,4] +
                                               exp_softmax[,i,5])
      }
    }
    V5[,tt,] = apply(est_prob,c(1,3),mean)
    message(tt)
  }
  
  V6 = array(NA,c(1000,12,5))
  for(tt in 2:12){
    idc_1 = contr_2
    idc_1 = idc_1[which(MT[idc_1,tt] == 2  & GHQ_SFS[idc_1,tt-1] == 1 & idc_1 %in% T_ind[T_ind[,2]==tt,1])] 
    exp_softmax = array(NA,c(1000,length(idc_1),5))
    
    for(j in 1:(5)){
      if(length(idc_1)!=1){
        linear_pred = a_coefficients[,,j] %*% t(X_matrices[idc_1,tt,1:29]) + 
          a_coefficients_t[,,j] %*% t(X_matrices[idc_1,tt-1,30:33]) + 
          as.vector(d[,tt,j]) + g[,df_numeric[idc_1,tt],j] + 
          e[,idc_1,j] + as.vector(b0[,j])*tt  +
          as.vector(a0[,j]) +psi[,PSU_trans[idc_1,tt],j] + 
          k[,STR_trans[idc_1,tt],j]
        exp_softmax[,,j] = exp(linear_pred)
      }else{
        linear_pred = a_coefficients[,,j] %*% (X_matrices[idc_1,tt,1:29]) + 
          a_coefficients_t[,,j] %*% (X_matrices[idc_1,tt-1,30:33]) + 
          as.vector(d[,tt,j]) + g[,df_numeric[idc_1,tt],j] + 
          e[,idc_1,j] + as.vector(b0[,j])*tt  +
          as.vector(a0[,j]) +psi[,PSU_trans[idc_1,tt],j] + 
          k[,STR_trans[idc_1,tt],j]
        exp_softmax[,,j] = exp(linear_pred)
      }
    }
    est_prob = array(NA,c(1000,length(idc_1),5))
    for(i in 1:length(idc_1)){
      for(j in 1:5){
        est_prob[,i,j] =  exp_softmax[,i,j]/(exp_softmax[,i,1] + exp_softmax[,i,2] +
                                               exp_softmax[,i,3] + exp_softmax[,i,4] +
                                               exp_softmax[,i,5])
      }
    }
    V6[,tt,] = apply(est_prob,c(1,3),mean)
    message(tt)
  }

  saveRDS(V5, file = "AR_No_all_years_Sin.RDS")
  saveRDS(V6, file = "AR_No_all_years_Rel.RDS")  

