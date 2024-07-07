#### additional functions ------------------------------------------------------

#' Insert founders in pedigree data
#' 
#' @title Insert Founders in pedigree data  
#' @param ped the data frame with tree columns for Individual, female and male
#' @param founders can be add founders information 
#' @return will be return new pedigree file with founders added
#' 
insertPedM <- function (ped, founders = NULL){
  ped[, 1] <- as.character(ped[, 1])
  ped[, 2] <- as.character(ped[, 2])
  ped[, 3] <- as.character(ped[, 3])
  mmothers <- na.omit(ped[, 2][which(ped[, 2] %in% ped[, 1] == FALSE)])
  mfathers <- na.omit(ped[, 3][which(ped[, 3] %in% ped[, 1] == FALSE)])
  if (is.null(founders) == FALSE) {
    founders <- na.omit(founders[which(founders %in% ped[, 1] == FALSE)])
  }
  mparents <- unique(c(mmothers, mfathers, founders))
  nped <- ped[rep(1, length(mparents)), ]
  nped[, 1] <- mparents
  nped[, 2] <- 0
  nped[, 3] <- 0
  nped <- rbind(nped, ped)
  colnames(nped) <- colnames(ped)
  row.names(nped)  <- 1:nrow(nped)
  return(nped)
}

#' transform pedigree in numeric ID
#' 
#' @title renumering pedigree
#' @param ped the data frame with tree columns for Individual, female and male
#' @return Are returned a list with pedigree data frame, in new numeric scale and a second data frame  with IDs link from old to new pedigree.
#' 
renum <- function(ped){## from packager ggroups
  colnames(ped) = c("ID", "SIRE", "DAM")
  for (i in 1:3) ped[, i] = paste0("x", ped[, i])
  ped[ped == "x0"] = 0
  newped = ped
  xrf = data.frame(ID = c(), newID = c())
  curr.set = ped[ped$SIRE == 0 & ped$DAM == 0, ]$ID
  xrf = data.frame(ID = curr.set, newID = 1:length(curr.set))
  ped = ped[!ped$ID %in% curr.set, ]
  gen = 1
  while (nrow(ped) > 0) {
    curr.set = ped[!ped$SIRE %in% ped$ID & !ped$DAM %in% ped$ID, ]
    curr.set = curr.set[order(curr.set$DAM, curr.set$SIRE), ]$ID
    xrf = rbind(xrf, data.frame(ID = curr.set, newID = (xrf[nrow(xrf), 
    ]$newID + 1):(xrf[nrow(xrf), ]$newID + length(curr.set))))
    ped = ped[!ped$ID %in% curr.set, ]
    gen = gen + 1
  }
  newped[] = xrf$newID[match(unlist(newped), xrf$ID)]
  newped[is.na(newped)] = 0
  newped = newped[order(newped$ID), ]
  xrf$ID = substring(xrf$ID, 2)
  #message("Found ", gen, " generations")
  return(list(newped = newped, xrf = xrf))
}

#' Overlay matrix From Sommer packager 
#' @title Overlay matrix of design
#' @param name description
#' @details rom Sommer packager
#' @return will be return matrix with combined design matrix 
#' 
overlayM <- function (..., rlist = NULL, prefix = NULL, sparse = FALSE){
  init <- list(...)
  myTypes <- unlist(lapply(init, class))
  init0 <- init
  init <- lapply(init, as.character)
  namesInit <- as.character(substitute(list(...)))[-1L]
  dat <- as.data.frame(do.call(cbind, init))
  dat <- as.data.frame(dat)
  for (j in 1:length(myTypes)) {
    if (myTypes[j] == "factor") {
      levels(dat[, j]) <- c(levels(dat[, j]), setdiff(levels(init0[[j]]), 
                                                      levels(dat[, j])))
    }
  }
  if (is.null(dim(dat))) {
    stop("Please provide a data frame to the overlay function, not a vector.\\n", 
         call. = FALSE)
  }
  if (is.null(rlist)) {
    rlist <- as.list(rep(1, dim(dat)[2]))
  }
  ss1 <- colnames(dat)
  dat2 <- as.data.frame(dat[, ss1])
  head(dat2)
  colnames(dat2) <- ss1
  femlist <- list()
  S1list <- list()
  for (i in 1:length(ss1)) {
    femlist[[i]] <- ss1[i]
    dat2[, femlist[[i]]] <- as.factor(dat2[, femlist[[i]]])
    if (sparse) {
      S1 <- Matrix::sparse.model.matrix(as.formula(paste("~", 
                                                         femlist[[i]], "-1")), dat2)
    }
    else {
      S1 <- model.matrix(as.formula(paste("~", femlist[[i]], 
                                          "-1")), dat2)
    }
    colnames(S1) <- gsub(femlist[[i]], "", colnames(S1))
    S1list[[i]] <- S1
  }
  levo <- sort(unique(unlist(lapply(S1list, function(x) {
    colnames(x)
  }))))
  if (sparse) {
    S3 <- Matrix(0, nrow = dim(dat2)[1], ncol = length(levo))
  }
  else {
    S3 <- matrix(0, nrow = dim(dat2)[1], ncol = length(levo))
  }
  rownames(S3) <- rownames(dat2)
  colnames(S3) <- levo
  for (i in 1:length(S1list)) {
    if (i == 1) {
      S3[rownames(S1list[[i]]), colnames(S1list[[i]])] <- S1list[[i]] * 
        rlist[[i]]
    }
    else {
      S3[rownames(S1list[[i]]), colnames(S1list[[i]])] <- S3[rownames(S1list[[i]]), 
                                                             colnames(S1list[[i]])] + (S1list[[i]][rownames(S1list[[i]]), 
                                                                                                   colnames(S1list[[i]])] * rlist[[i]])
    }
  }
  if (!is.null(prefix)) {
    colnames(S3) <- paste(prefix, colnames(S3), sep = "")
  }
  attr(S3, "variables") <- namesInit
  return(S3)
}

#' Insert founders in pedigree data
#' 
#' @title 
#' @param name description
#' @param name description
#' @details SASA
#' @return description 
#' 
check_data_adj <- function(data=df_raw, traits_con = c('LDSR','FSEV','BSEV','YIELDCT','YIELDWT','VSEV','DSR')){
  data <- data |> 
    dplyr::mutate(observation_trait_val=as.numeric(observation_trait_val)) |> 
    dplyr::filter(!is.na(observation_trait_val))
  
  #### Adjust LDSR, BSEV, FSEV
  data[data$TRAIT %in% traits_con,'TRAIT'] <- 
    apply(data[data$TRAIT %in% traits_con,
               c('TRAIT','TRAIT_ABBREVIATION')],
          1,function(x){gsub('[ ;:]','_',paste(x,collapse = '_'))})
  
  
  data$TRAIT <- gsub('DSR_PATHOGEN','DSR',data$TRAIT)
  data$TRAIT <- gsub('__TISTY_','_',data$TRAIT)
  data$TRAIT <- gsub(' ','_',data$TRAIT)
  
  return(data)
}
#' Insert founders in pedigree data
#' 
#' @title 
#' @param name description
#' @param name description
#' @details SASA
#' @return description 
#' 

source_all_ID <- function(ID=ID,...){
  project <- "bcs-veg-rd" 
  sqlped <- paste("
SELECT  distinct  
CAST(fts_germplasm_id AS INT64) AS fts_germplasm_id,
  line_type,
    pedigree,
    origin,
    cross_name,
    germplasm_lineage_identifier
FROM
  `bcs-veg-mart.pipeline.veg_germplasms`
WHERE
  cross_name IN ('",paste(ID,collapse="','")  ,"')
  OR germplasm_lineage_identifier  IN ('",paste(ID,collapse="','")  ,"')
  OR pedigree IN ('",paste(ID,collapse="','")  ,"')
   OR fts_germplasm_id IN (",paste(na.omit(as.integer64(ID)),collapse=",")  ,")
                  ;",sep='')
  tb <- bigrquery::bq_project_query(project, sqlped)
  df_id <- bigrquery::bq_table_download(tb,bigint="integer64",page_size = 20000) |> 
    unique() 
  return(df_id)
}

#' Insert founders in pedigree data
#' 
#' @title 
#' @param name description
#' @param name description
#' @details SASA
#' @return description 
#' 

source_all_ID_ByOrigin_ID <- function(fts_germplasm_id,...){
  project <- "bcs-veg-rd" 
  sqlped <- paste("
SELECT 
DISTINCT
  CAST(vg.fts_germplasm_id AS INT64) AS fts_germplasm_id,
  vg.line_type,
  vg.pedigree,
  vg.origin,
  vg.cross_name,
  vg.germplasm_lineage_identifier,
  prod.primary_commercial_product_nm as PROD,
  CAST(pbc.male_parent_germplasm_id AS INT64) AS male_parent_germplasm_id,
  CAST(pbc.female_parent_germplasm_id AS INT64) AS female_parent_germplasm_id
  ##pbc.crosses_removed_nbr
FROM
  `bcs-veg-mart.pipeline.veg_germplasms` vg
  LEFT JOIN UNNEST(products) AS prod
  LEFT JOIN UNNEST(previous_5_binary_crosses) AS pbc
WHERE
  vg.origin IN (
  SELECT
  DISTINCT
    origin
  FROM `bcs-veg-mart.pipeline.veg_germplasms`
  WHERE
    fts_germplasm_id IN (",paste(unique(fts_germplasm_id),collapse=",")  ,"))
 AND (crosses_removed_nbr IS NULL OR  crosses_removed_nbr = 0);",sep='')
  
  tb <- bigrquery::bq_project_query(project, sqlped)
  df_id <- bigrquery::bq_table_download(tb,bigint="integer64",page_size = 20000) |> 
    unique() 
  ##df_id$fts_germplasm_id <- bit64::as.integer64(df_id$fts_germplasm_id)
  return(df_id)
}


#' Insert founders in pedigree data
#' 
#' @title 
#' @param name description
#' @param name description
#' @details SASA
#' @return description 
#' 

source_all_ID_ByOrigin <- function(origin,...){
  project <- "bcs-veg-rd" 
  sqlped <- paste("
SELECT  distinct  
CAST(fts_germplasm_id AS INT64) AS fts_germplasm_id,
  line_type,
    pedigree,
    origin,
    cross_name,
    germplasm_lineage_identifier,
    prod.primary_commercial_product_nm as PROD
FROM
  `bcs-veg-mart.pipeline.veg_germplasms`
  LEFT JOIN UNNEST(products) AS prod
WHERE
  origin IN ('",paste(unique(origin),collapse="','")  ,"')
  ;",sep='')
  tb <- bigrquery::bq_project_query(project, sqlped)
  df_id <- bigrquery::bq_table_download(tb,bigint="integer64",page_size = 20000) |> 
    unique() 
  return(df_id)
}




#' Insert founders in pedigree data
#' 
#' @title 
#' @param name description
#' @param name description
#' @details SASA
#' @return description 
#' 

source_ped <- function(cross_name,debug_folder,...){
  cross_name <- unique(cross_name)
  project <- "bcs-veg-rd" 
  sqlped <- paste("
SELECT  
distinct  
CAST(veoe.fts_germplasm_id AS INT64) AS fts_germplasm_id,
  veoe.line_type,
  veoe.pedigree,
  veoe.cross_name,
  veoe.germplasm_lineage_identifier,
  prod.primary_commercial_product_nm as PROD
FROM
  `bcs-veg-mart.pipeline.veg_germplasms` veoe
  left join (`bcs-veg-mart.pipeline.veg_germplasms` vg CROSS JOIN UNNEST (products) as prod) on veoe.fts_germplasm_id = vg.fts_germplasm_id
WHERE
  veoe.cross_name IN ('",paste(cross_name,collapse="','")  ,"')

UNION ALL
  SELECT  
distinct  
CAST(veoe.fts_germplasm_id AS INT64) AS fts_germplasm_id,
  veoe.line_type,
  veoe.pedigree,
  veoe.cross_name,
  veoe.germplasm_lineage_identifier,
  prod.primary_commercial_product_nm as PROD
FROM
  `bcs-veg-mart.pipeline.veg_germplasms` veoe
  left join (`bcs-veg-mart.pipeline.veg_germplasms` vg CROSS JOIN UNNEST (products) as prod) on veoe.fts_germplasm_id = vg.fts_germplasm_id
WHERE
  veoe.germplasm_lineage_identifier  IN ('",paste(cross_name,collapse="','")  ,"')

UNION ALL
  SELECT  
distinct  
CAST(veoe.fts_germplasm_id AS INT64) AS fts_germplasm_id,
  veoe.line_type,
  veoe.pedigree,
  veoe.cross_name,
  veoe.germplasm_lineage_identifier,
  prod.primary_commercial_product_nm as PROD
FROM
  `bcs-veg-mart.pipeline.veg_germplasms` veoe
  left join (`bcs-veg-mart.pipeline.veg_germplasms` vg CROSS JOIN UNNEST (products) as prod) on veoe.fts_germplasm_id = vg.fts_germplasm_id
WHERE
  veoe.pedigree IN ('",paste(cross_name,collapse="','")  ,"')
;",sep='')
  
  tb <- bigrquery::bq_project_query(project, sqlped)
  df_append <- bigrquery::bq_table_download(tb,bigint="integer64",page_size = 20000) |> 
    unique() |>  
    as.data.frame()
  not_found_cross_name <- df_append[apply(df_append,1,function(x){sum(is.na(x))})>0,]
  if(nrow(not_found_cross_name)>1){
    write.csv(not_found_cross_name,paste0(debug_folder,'\\out_InbredAppend_crossname_LID_notcomplete.csv'),row.names = FALSE)
  }
  return(df_append)
}

#' Insert founders in pedigree data
#' 
#' @title 
#' @param name description
#' @param name description
#' @details SASA
#' @return description 
#' 
source_gc_append_byLID <- function(fts_germplasm_id=fts_germplasm_id,traits='HETGP'){
  project <- "bcs-veg-rd" 
  sqlped <- paste("
  with ref as (
    SELECT 
    DISTINCT
    ##vg.fts_germplasm_id,
    ##gc.characteristic_val as LID
    vg.cross_name
    FROM
    `bcs-veg-mart.pipeline.veg_germplasms` vg,
    UNNEST(germplasm_characteristics) gc
    WHERE
    fts_germplasm_id IN (",paste(fts_germplasm_id,collapse=',')  ,")
    AND gc.characteristic_reference_cd = 'LID')
    
    SELECT DISTINCT
    CAST(vg.fts_germplasm_id AS INT64) AS fts_germplasm_id,
    vg.pedigree,
    vg.origin,
    vg.line_type,
    vg.germplasm_lineage_identifier,
    vg.cross_name,
    prod.primary_commercial_product_nm as PROD,
    gc.characteristic_reference_cd,
    gc.characteristic_val
    FROM
    `bcs-veg-mart.pipeline.veg_germplasms` vg
    LEFT JOIN UNNEST(germplasm_characteristics) gc
    LEFT JOIN UNNEST(products) AS prod
    WHERE
    ### fts_germplasm_id IN (select fts_germplasm_id FROM ref)
    cross_name IN (select cross_name FROM ref)
    ###AND gc.characteristic_reference_cd IN ('PMC','HETGP')
    AND gc.characteristic_reference_cd IN ('",paste(traits,collapse="','"),"')
    ;",sep='')
  
  tb <-  bigrquery::bq_project_query(project, sqlped)
  df_add <-  bigrquery::bq_table_download(tb,bigint="integer64",page_size = 20000) |> as.data.frame()
  
  df_add <- df_add |> 
    pivot_wider(id_cols =c('fts_germplasm_id',
                           'pedigree',
                           'origin',
                           'line_type',
                           'germplasm_lineage_identifier',
                           'cross_name',
                           'PROD'),
                names_from  = c('characteristic_reference_cd'),
                values_from =c('characteristic_val'),
                values_fn = max) |> 
    dplyr::mutate_all( ~replace(., lengths(.)==0, NA)) |> 
    as.data.frame()
  #head(df_add)
  return(df_add)
}

#' Insert founders in pedigree data
#' 
#' @title 
#' @param name description
#' @param name description
#' @details SASA
#' @return description 
#' 
source_gc_append <- function(fts_germplasm_id=fts_germplasm_id,
                             traits=c('HETGP','P50','S50','HBC','P1DATE','P2DATE',
                                      'PMC','PRODSTG','ENC')){
  project <- "bcs-veg-rd" 
  sqlped <- paste("SELECT DISTINCT
  CAST(fts_germplasm_id AS INT64) AS fts_germplasm_id,
  pedigree,
  origin,
  germplasm_lineage_identifier,
  line_type,
  cross_name,
  prod.primary_commercial_product_nm as PROD,
  characteristic_reference_cd,
  characteristic_val
 FROM
  `bcs-veg-mart.pipeline.veg_germplasms`
  LEFT JOIN UNNEST(germplasm_characteristics) gc
  LEFT JOIN UNNEST(products) AS prod
WHERE
  fts_germplasm_id IN (",paste(fts_germplasm_id,collapse=',')  ,")
 ;",sep='')
  
  tb <-  bigrquery::bq_project_query(project, sqlped)
  df_add <-  bigrquery::bq_table_download(tb,bigint="integer64",page_size = 20000) |> as.data.frame()
  
  df_add <- df_add |> 
    pivot_wider(id_cols =c('fts_germplasm_id',
                           'pedigree',
                           'origin',
                           'line_type',
                           'germplasm_lineage_identifier',
                           'cross_name',
                           'PROD'),
                names_from  = c('characteristic_reference_cd'),
                values_from = c('characteristic_val'),
                values_fn = max) |> 
    dplyr::mutate_all( ~replace(., lengths(.)==0, NA)) |> 
    as.data.frame()
  
  df_add <- df_add[,na.omit(base::match(c('fts_germplasm_id',
                                          'pedigree',
                                          'origin',
                                          'line_type',
                                          'germplasm_lineage_identifier',
                                          'cross_name',
                                          'PROD',traits), colnames(df_add)))]
  
  return(df_add)
}

#' Insert founders in pedigree data
#' 
#' @title 
#' @param name description
#' @param name description
#' @details SASA
#' @return description 
#' 
source_deep_ped2 <- function(fts_germplasm_id=fts_germplasm_id){
  project <- "bcs-veg-rd" 
  sqlped <- paste("
  SELECT DISTINCT
  CAST(pbc.child_germplasm_id AS INT64) AS child_germplasm_id,
  CAST(pbc.male_parent_germplasm_id AS INT64) AS male_parent_germplasm_id,
  CAST(pbc.female_parent_germplasm_id AS INT64) AS female_parent_germplasm_id
FROM
  `bcs-veg-mart.pipeline.veg_germplasms` vg
  LEFT JOIN UNNEST(previous_5_binary_crosses) AS pbc
WHERE
  vg.fts_germplasm_id IN (",paste(fts_germplasm_id,collapse=',')  ,");",sep='')
  
  tb <- bigrquery::bq_project_query(project, sqlped)
  df_ped <- bigrquery::bq_table_download(tb,bigint="integer64",
                                         page_size = 20000) |> 
    dplyr::filter(!is.na(child_germplasm_id)) |> 
    as.data.frame()
  return(df_ped)
}


#' Insert founders in pedigree data
#' 
#' @title 
#' @param name description
#' @param name description
#' @details SASA
#' @return description 
#' 
source_deep_ped <- function(fts_germplasm_id=fts_germplasm_id){
  project <- "bcs-veg-rd" 
  sqlped <- paste("
  SELECT
  DISTINCT
  CAST(vg.fts_germplasm_id AS INT64) AS fts_germplasm_id,
  vg.pedigree,
  vg.origin,
  vg.cross_name,
  vg.germplasm_lineage_identifier,
  prod.primary_commercial_product_nm as PROD,
  pbc.crosses_removed_nbr,
  CAST(pbc.child_germplasm_id AS INT64) AS  child_germplasm_id,
  CAST(pbc.male_parent_germplasm_id) AS male_parent_germplasm_id,
  CAST(pbc.female_parent_germplasm_id AS INT64) AS female_parent_germplasm_id,
FROM
  `bcs-veg-mart.pipeline.veg_germplasms` vg
  LEFT JOIN UNNEST(previous_5_binary_crosses) AS pbc
  LEFT JOIN UNNEST(products) AS prod
WHERE
  vg.fts_germplasm_id IN (",paste(fts_germplasm_id,collapse=',')  ,");",sep='')
  
  tb <- bigrquery::bq_project_query(project, sqlped)
  df_ped <- bigrquery::bq_table_download(tb,bigint="integer64",page_size = 20000)  |> 
    as.data.frame()
  return(df_ped)
}

#' Insert founders in pedigree data
#' 
#' @title 
#' @param name description
#' @param name description
#' @details SASA
#' @return description 
#' 

imp_fake_parentage2 <- function(df_ped=df_ped,data=df_raw){
  df_pedsel <- df_ped |> 
    dplyr::select(MIN_fts_germplasm_id,
                  MIN_male_parent_germplasm_id,
                  MIN_female_parent_germplasm_id) |>
    arrange(MIN_fts_germplasm_id) |> unique()
  
  Noparentage <- unique(data[!(data$MIN_fts_germplasm_id  %in% 
                                 df_pedsel$MIN_fts_germplasm_id),
                             c('MIN_fts_germplasm_id',
                               'MIN_male_parent_germplasm_id',
                               'MIN_female_parent_germplasm_id')]) |>
    as.data.frame()
  colnames(Noparentage) <- colnames(df_pedsel)
  
  con <- 1
  for(i in 1:nrow(Noparentage)){
    Noparentage$MIN_female_parent_germplasm_id[i] <- con
    Noparentage$MIN_male_parent_germplasm_id[i] <- con+1
    con <- con+2
  }
  return(Noparentage)
}

#' Insert founders in pedigree data
#' 
#' @title 
#' @param name description
#' @param name description
#' @details SASA
#' @return description 
#' 
converteM <- function(g){
  N <- nchar(g)
  g[N!=3] <- NA
  g <- gsub('\\[|\\,|\\]','',g)
  g[g=='**'] <- NA
  out <- cbind(substr(g,1,1),substr(g,2,2))
  out[out=='A'] <- 1
  out[out=='C'] <- 2
  out[out=='G'] <- 3
  out[out=='T'] <- 4
  out <- as.numeric(apply(out,1,function(x){paste0(x,collapse = '')}))
  return(out)
}


#' Insert founders in pedigree data
#' 
#' @title 
#' @param name description
#' @param name description
#' @details SASA
#' @return description 
#' 
converte <- function(g){
  out <- cbind(substr(g,1,1),substr(g,2,2))
  out[out=='A'] <- 1
  out[out=='C'] <- 2
  out[out=='G'] <- 3
  out[out=='T'] <- 4
  out <- as.numeric(apply(out,1,function(x){paste0(x,collapse = '')}))
  return(out)
}

#' Insert founders in pedigree data
#' 
#' @title 
#' @param name description
#' @param name description
#' @details SASA
#' @return description 
#' 
source_genotype_imp <- function(fts_germplasm_id){
  project <- "bcs-breeding-datasets" 
  sql <- paste0("SELECT
  CAST(MIDAS_GERMPLASM_ID as INT64) as MIDAS_GERMPLASM_ID,
  MARKER,
 SCORE
FROM
  `bcs-breeding-datasets.breeding_genomics.fpbeagleimputed_markers_corn`
WHERE
  MIDAS_GERMPLASM_ID IN (",paste0(fts_germplasm_id,collapse=","),")
;",sep='') 
  
  tb <- bq_project_query(project, sql)
  df_gen <- bq_table_download(tb,bigint="integer64",page_size = 20000) 
  df_gen <- df_gen |>  
    dplyr::mutate(SCORE=as.integer(SCORE+1)) |> 
    pivot_wider(id_cols = c('MIDAS_GERMPLASM_ID'), 
                names_from = MARKER, 
                values_from = SCORE,
                values_fill =NA,
                values_fn=max)  |> 
    as.data.frame() 
  colnames(df_gen) <- gsub('\\.', '_', colnames(df_gen))
  rownames(df_gen) <- df_gen[,1]
  df_gen <- as.matrix(df_gen[,-1])
  return(df_gen)
}

#' Insert founders in pedigree data
#' 
#' @title 
#' @param name description
#' @param name description
#' @details SASA
#' @return description 
#' 
source_genotype_raw <- function(fts_germplasm_id){
  project <- "product360-datasets"
  sql <- paste0("SELECT
  sub.id AS midas_germplasm_id,
  cal.marker_name AS marker_name,
  cal.allele_calls AS score
FROM
  `product360-datasets.genomics.MarkerCallSetsV2`,
  UNNEST(subjects) AS sub,
  UNNEST(calls) AS cal
WHERE
  sub.id IN ('",paste0(fts_germplasm_id,collapse="','"),"')
;",sep='') 
  
  tb <- bq_project_query(project, sql)
  df_gen <- bq_table_download(tb,bigint="integer64",page_size = 20000) 
  df_gen$score <- apply(plyr::ldply(df_gen$score),1,function(x) paste(x,collapse=','))
  df_gen$score1 <- converteM(df_gen$score)
  df_gen <- na.omit(df_gen)
  df_gen <- df_gen |>  
    pivot_wider(id_cols = c('midas_germplasm_id'), 
                names_from = marker_name, 
                values_from = score1,
                values_fill =NA,
                values_fn=max)  |> 
    as.data.frame() 
  colnames(df_gen) <- gsub('\\.', '_', colnames(df_gen))
  rownames(df_gen) <- df_gen[,1]
  df_gen <- df_gen[,apply(df_gen,2,function(x){sum(!is.na(x))}) >= nrow(df_gen)*0.80]
  df_gen <- (df_gen[,-1])
  for (i in 1:ncol(df_gen)){
    x <- data.frame(m=(as.character(df_gen[,i])))
    scorS <- table(df_gen[,i])
    if(nrow(scorS)>1){
      scorS <- data.frame(sort(scorS,decreasing = TRUE))
      scorS$homo <- substr(scorS$Var1,1,1) ==  substr(scorS$Var1,2,2)
      scorS$score <- NA
      scorS$score[scorS$homo == FALSE] <- 1
      scorS$score[scorS$homo == TRUE] <- c(2,0)
      scorS$Var1 <- (as.character(scorS$Var1))
      x$score <- NA
      for(k in 1:nrow(scorS)){
        x$score[x$m==scorS$Var1[k]] <- scorS$score[k]
      }
      df_gen[,i] <- x$score
    }else{
      df_gen[,i] <- 2
    }
  }
  df_gen <- as.matrix(df_gen[,-1])
  return(df_gen)
}

#' Insert founders in pedigree data
#' 
#' @title 
#' @param name description
#' @param name description
#' @details SASA
#' @return description 
#' 

makeHinv <- function(A=Am,Ai=Ai,G=Gm,tau=0.1,omega=0.1){
  ## based in function Hmatrix from AGHmatrix, reference: Martini et all 2018
  idA <- rownames(A)
  idAi <- rownames(Ai)
  idG <- rownames(G)
  if(sum(is.na(match(idG,idA)))>0){
    stop("G containg individual not include in Ainv")
  }
  if(sum(is.na(match(idA,idAi)))>0 | sum(is.na(match(idAi,idA)))>0){
    stop("Ai and A are was diferent individuals, please check")
  }
  
  idH <- unique(c(idG,idAi))
  Ai <- Ai[idH,idH]
  index <- is.na(match(idH,idG))
  Ai11 <- Ai[index,index]
  Ai12 <- Ai[index,!index]
  Ai21 <- Ai[!index,index]
  
  Ai11[Ai11 != 0] <- 0
  Ai12[Ai12 != 0] <- 0
  Ai21[Ai21 != 0] <- 0
  
  A22 <- A[!index,!index]
  A22inv <- solve(A22)
  
  G22 <- G[idH[!index],idH[!index]]
  G22inv = try(solve(G22),silent=TRUE)
  if(inherits(G22inv,"try-error")){
    cat(G22inv)
    stop("G22 not inverting with solve(), try a different/modified G matrix")
  }
  
  AG22 <- (tau*G22inv - omega*A22inv)
  
  Hinv <-  Ai+cbind(rbind(Ai11,Ai21),rbind(Ai12,AG22))
  Hinv <- Matrix::forceSymmetric(Hinv)
  return(Hinv)
}

#' Insert founders in pedigree data
#' 
#' @title 
#' @param name description
#' @param name description
#' @details SASA
#' @return description 
#' 
extract_trans_to_var <- function(fit,precision="Log precision for idz",...){
  INLA::inla.emarginal(function(x) x,
                       INLA::inla.tmarginal(function(x) 1/exp(x),
                                            fit$internal.marginals.hyperpar[[precision]]))
}

#' Insert founders in pedigree data
#' 
#' @title 
#' @param name description
#' @param name description
#' @details SASA
#' @return description 
makeZ <- function(data=df_model,female='femalef',male='malef',Cov_m=Cov_m){
  Z <- overlayM(data[,female],data[,male])
  Cov_mf <- colnames(Cov_m)[!(colnames(Cov_m) %in% colnames(Z))]
  AdM <- matrix(0,ncol=length(Cov_mf),nrow=nrow(Z))
  colnames(AdM) <- Cov_mf
  Z <- cbind(Z,AdM)
  Z <- Matrix::Matrix(Z[,colnames(Cov_m)],sparse = TRUE)
  return(Z)
}




#' TRUE names for use 
#' 
#' @title 
#' @param name description
#' @param name description
#' @details SASA
#' @return description 
#' 
#' 
# source_rename_fts <- function(crop='Corn'){
#   project <- "bcs-veg-rd" 
#   sqlped <- paste("SELECT DISTINCT
#   CAST(fts_germplasm_id AS INT64) AS fts_germplasm_id,
#   pedigree,
#   germplasm_lineage_identifier,
#   cross_name
# FROM
#   `bcs-veg-mart.pipeline.veg_germplasms`
# WHERE
#   crop = '",crop,"';",sep='')
#   tb <- bigrquery::bq_project_query(project, sqlped)
#   df_name <- bigrquery::bq_table_download(tb,bigint="integer64",page_size = 20000)  |> 
#     as.data.frame()
#   df_name$KEY_ID <- df_name$cross_name
#   df_name$KEY_ID[is.na(df_name$KEY_ID)] <- df_name$germplasm_lineage_identifier[is.na(df_name$KEY_ID)]
#   df_name$KEY_ID[is.na(df_name$KEY_ID)] <- df_name$pedigree[is.na(df_name$KEY_ID)]
#   df_name <- df_name |> dplyr::select(fts_germplasm_id,KEY_ID)
#   return(df_name)
# }



#' create key ID
#' 
#' @title 
#' @param name description
#' @param name description
#' @details SASA
#' @return description 
#' 
#' 
makeKEY_ID <- function(data=df_raw){
  if(sum(!(c('cross_name','germplasm_lineage_identifier','pedigree','PROD') %in% colnames(data)))>0){
    stop('Please provide a data frame with: cross_name, germplasm_lineage_identifier, pedigree, PROD')
  }
  ID_PROD <- data |> 
    dplyr::select(cross_name,PROD) |> 
    dplyr::mutate(cross_name = ifelse(cross_name == '',NA,cross_name) ,
                  PROD = ifelse(PROD == '',NA,PROD)) |> 
    unique() |> 
    na.omit()
  if(nrow(ID_PROD)>0){
    data <- data |>
      select(-PROD) |> 
      left_join(ID_PROD,
                relationship= 'many-to-one',
                by = join_by(cross_name),
                multiple='first')
  }
  data$KEY_ID <- as.character(data$PROD)
  data$KEY_ID[data$KEY_ID==''] <- NA
  data$KEY_ID[is.na(data$KEY_ID)] <- data$cross_name[is.na(data$KEY_ID)]
  data$KEY_ID[data$KEY_ID==''] <- NA
  data$KEY_ID[is.na(data$KEY_ID)] <- data$germplasm_lineage_identifier[is.na(data$KEY_ID)]
  data$KEY_ID[data$KEY_ID==''] <- NA
  data$KEY_ID[is.na(data$KEY_ID)] <- data$pedigree[is.na(data$KEY_ID)]
  data$KEY_ID[data$KEY_ID=='NA'] <- NA
  data <- data |> arrange(KEY_ID)
  return(data)
}


#' create key ID
#' 
#' @title 
#' @param name description
#' @param name description
#' @details SASA
#' @return description 
#' 
#' 
source_inv_byHBCPED <- function(owner_breeding_prog='6S'){
  project <- "bcs-veg-rd" 
  sqlped <- paste0("
  SELECT 
  characteristic_val as HBCPED, 
  count(inv.inv_barcode) as Ninv
  FROM
  `bcs-veg-mart.pipeline.veg_germplasms`,
  UNNEST(inventories) inv,
  UNNEST(germplasm_characteristics) gc
  WHERE 
  line_type = 'Hybrid' AND
  owner_breeding_prog_code IN ('",paste0(owner_breeding_prog,collapse="','"),"') AND 
  gc.characteristic_reference_cd = 'HBCPED'
  group by  characteristic_val;")
  tb <-  bigrquery::bq_project_query(project, sqlped)
  data <-  bigrquery::bq_table_download(tb,bigint="integer64",page_size = 20000) |> as.data.frame()
  return(data)
}



#' 
#' @title 
#' @param name description
#' @param name description
#' @details SASA
#' @return description 
#' 
#' 
#' 
#' 
#' 
st_mean <- function(x=x,sel_dir=1,min_range=1,max_range=10){
  if(sel_dir %in% c(-1,1)){
    x_std <- sel_dir*((x-mean(x,na.rm=T))/sd(x,na.rm=T))
    x_std[is.na(x_std)] <- 0
  }
  if(sel_dir %in% c(0)){
    mp <- mean(c(min_range,max_range))
    dif <- mp-min_range 
    xl <- c(min(x),
            min(c(mp-2*dif,min(x))),
            min(c(mp-1.5*dif,min(x))),
            mp-1*dif,
            mp-0.5*dif,
            mp-0*dif,
            mp+0.5*dif,
            mp+1*dif,
            min(c(mp+1.5*dif,max(x))),
            min(c(mp+2*dif,max(x))),
            max(x))
    yl <- -1.64+((2*1.64)*c(0.01,0.01,0.01,0.95,0.98,1.00,0.98,0.95,0.01,0.01,0.01))
    
    spl_func <- approxfun(xl,yl)
    x_std <- spl_func(x)
    x_std[is.na(x_std)] <- 0
  }
  return(x_std)
}


#' create key ID
#' 
#' @title 
#' @param name description
#' @param name description
#' @details SASA
#' @return description 
#' 
#' 
adv_plot <- function(data = df_sel,
                     x='PCTRC_predicted.value',
                     y='SCTPA_predicted.value',
                     z='LDSR_SPIRKU_predicted.value',
                     x_direct = "identity",#"reverse"
                     y_direct = "identity",#"reverse"
                     YLABS = 'SCTPA (t/ha)',
                     XLABS = 'PCTRC (%)',
                     ZLABS = 'SPIRKU',
                     prefix= 'g',
                     x_per_add=c(0.05,0.10),
                     y_per_add=c(0.05,0.10),
                     Nclass=5,
                     id_label='LIDm',
                     Stage='Stage',
                     x_add_margin = c(1,1),
                     y_add_margin = c(1,1),
                     ForceCheck=FALSE,
                     save_plot=F){
  
  if(is.null(z)){
    data <- data[,c(id_label,Stage,x,y)]
    colnames(data) <- c('id_label','Stage','x','y')
    z <- 1
    ZLABS <- NULL
  }else{
    data <- data[,c(id_label,Stage,x,y,z)]
    colnames(data) <- c('id_label','Stage','x','y','z')
  }
  check_mean <- c(mean(data$x[data$Stage %in% c('CHECK','CM')],na.rm=TRUE),
                  mean(data$y[data$Stage %in% c('CHECK','CM')],na.rm=TRUE))
  x_adj <- diff(range(data$x,na.rm = TRUE))*0.01
  y_adj <- diff(range(data$y,na.rm = TRUE))*0.015
  x_limits <- range(data$x,na.rm = T)
  y_limits <- range(data$y,na.rm = T)
  if(x_direct=="reverse"){
    x_per_add <- x_per_add[2:1]
    x_limits <- x_limits[2:1]
    x_add_margin <- x_add_margin[2:1]
    x_l1 <- check_mean[1]-abs(check_mean[1]*x_per_add[1])
    x_l2 <- check_mean[1]-abs(check_mean[1]*x_per_add[2])
  }else{
    x_l1 <- check_mean[1]+abs(check_mean[1]*x_per_add[1])
    x_l2 <- check_mean[1]+abs(check_mean[1]*x_per_add[2])
  }
  if(y_direct=="reverse"){
    y_per_add <- y_per_add[2:1]
    y_limits <- y_limits[2:1]
    y_add_margin <- y_add_margin[2:1]
    y_l1 <- check_mean[2]-abs(check_mean[2]*y_per_add[1])   
    y_l2 <- check_mean[2]-abs(check_mean[2]*y_per_add[2])
  }else{
    y_l1 <- check_mean[2]+abs(check_mean[2]*y_per_add[1])   
    y_l2 <- check_mean[2]+abs(check_mean[2]*y_per_add[2])
  }
  x_limits <- x_limits*x_add_margin
  y_limits <- y_limits*y_add_margin
  if(ForceCheck){
    data$id_label[!(data$Stage %in% c('CHECK','CM'))] <- ""
  }
  if(!(file.exists(file.path(getwd(),'out_plot')))){dir.create(file.path(getwd(),'out_plot'))}
  pp <- ggplot2::ggplot(data,aes(x=x,y=y,colour=Stage,label=id_label,size=z)) +
    ggplot2::geom_point()+
    ggplot2::scale_size_binned(name = ZLABS, n.breaks = Nclass,range = c(3, 7)) +
    ggplot2::guides(colour = guide_legend(override.aes = list(size=5),order = 1)) +
    ggplot2::xlab(XLABS)+
    ggplot2::ylab(YLABS)+
    ggplot2::scale_color_manual(breaks = c("CHECK","CM","SC2" ,"PS1", "PS2","PS3","PS4"),values=c("#F50057","#880E4F", "#01AAFF","#01BEFF", "#00314E","#32CD32","#006400")) +
    ggplot2::scale_x_continuous(limits = x_limits ,trans = x_direct) +
    ggplot2::scale_y_continuous(limits = y_limits ,trans = y_direct) +
    ggrepel::geom_text_repel() +
    ggplot2::geom_vline(xintercept = check_mean[1], color="gray40",lwd=1,linetype=2)+
    ggplot2::geom_hline(yintercept = check_mean[2], color="gray40",lwd=1,linetype=2)+
    ggplot2::geom_text(aes(x=check_mean[1],y=check_mean[2],label ='Check mean'),hjust = 0,nudge_x =  x_adj,nudge_y = y_adj,size=3,color="gray40") +
    ggplot2::geom_vline(xintercept = x_l1,color="gray60",lwd=1,linetype=3) +
    ggplot2::geom_hline(yintercept = y_l1,color="gray60",lwd=1,linetype=3) +
    ggplot2::geom_text(aes(x=x_l1, y=y_l1,label =paste('Check mean: ',round((x_per_add[1])*100,0),';',round((y_per_add[1])*100,0),'%')),
                       hjust = 0,nudge_x = x_adj,nudge_y = y_adj,size=3,color="gray60") + 
    ggplot2::geom_vline(xintercept = x_l2, color="gray80",lwd=1,linetype=4) +
    ggplot2::geom_hline(yintercept = y_l2, color="gray80",lwd=1,linetype=4) +
    ggplot2::geom_text(aes(x=x_l2,y=y_l2, label =paste('Check mean: ',round((x_per_add[2])*100,0),';',round((y_per_add[2])*100,0),'%')),
                       hjust = 0, nudge_x =  x_adj,nudge_y = y_adj,size=3,color="gray80") +
    ggplot2::theme_bw() +
    ggplot2::theme(text = element_text(size = 15),
                   legend.title = element_text(size= 15),
                   legend.text = element_text(size= 15),
                   plot.title = element_text(size = 15),
                   axis.title = element_text(size = 15),
                   axis.text.x = element_text(size = 15),
                   axis.text.y = element_text(size = 15))
  if(save_plot){
    ggplot2::ggsave(paste('out_plot/out_',prefix,'_',x,y,z,'.png',sep=''),w=12,h=8)
  }
  return(pp)
}


#' Make TIX 
#' 
#' @title 
#' @param name description
#' @param name description
#' @details SASA
#' @return description 
#' 
#' 
makeTIX <- function(dataP,traits){
  checkTraits <- data.frame(traits=unique(traits$trait),
                            include=unique(traits$trait) %in% colnames(dataP))
  checkTraits <- checkTraits |> dplyr::filter(include == FALSE)
  if(nrow(checkTraits)>0){
    warning(paste0('The trait: ',checkTraits$traits,' is no include in BLUP data avaliable, set 0 for all hybrid \n'))
  }
  df_fake <- as.data.frame(matrix(0,nrow = nrow(dataP),ncol = nrow(checkTraits)))
  colnames(df_fake) <- checkTraits$traits
  dataP <- cbind(dataP,df_fake)
  
  pdf(paste('post_Make_TIX/out_plotTIX.pdf',sep=''),w=12,h=8)
  for(it in unique(traits$TIX_name)){
    iTIX <- paste0(it,'_TIX')
    for(i in unique(traits$trait[traits$TIX_name==it])){
      if(i %in% colnames(dataP)){
        inamestan <- paste(it,i,'predS',sep='_')
        dataP <- dataP %>% 
          dplyr::mutate({{inamestan}} := st_mean(x=dataP[,i],
                                                 sel_dir=traits$trait_direction[traits$trait==i & traits$TIX_name==it],
                                                 min_range=traits$Minimal[traits$trait==i & traits$TIX_name==it],
                                                 max_range=traits$Maximum[traits$trait==i & traits$TIX_name==it]) *
                          traits$Weight[traits$trait==i & traits$TIX_name==it])
        
        pp <- ggplot(data.frame(x=dataP[,i],x_std=dataP[,inamestan]),aes(x=x,y=x_std)) + 
          geom_point() + 
          ggtitle(inamestan)+
          xlab('BLUP')+
          ylab('Standarized BLUP for TIX')+
          theme_bw()
        print(pp)
      }
    }
    x_temp <- dataP[,grep(paste0('^',it),colnames(dataP))]
    if(is.null(ncol(x_temp))){
      dataP <- dataP %>% 
        dplyr::mutate({{iTIX}} := x_temp)
    }else{
      dataP <- dataP %>% 
        dplyr::mutate({{iTIX}} := rowSums(x_temp))
    }
  }
  dev.off()
  return(dataP)
}


#' runTIX
#' 
#' @title 
#' @param name description
#' @param name description
#' @details SASA
#' @return description 
#'  
runTIX <- function(data,traits_TIX,traits_app,f_stages,out_file){
  data <- data |> dplyr::filter(PRODSTG %in% f_stages)
  dataP <- data |>
    pivot_wider(id_cols = c(KEY_ID,cross_name,pedigree,origin,
                            germplasm_lineage_identifier,
                            fts_germplasm_idH,PMC,PRODSTG),
                names_from = trait,
                values_from = EBV_GCA_SCA,
                values_fn = mean) |>
    dplyr::mutate_all( ~replace(., lengths(.)==0, NA)) |> 
    as.data.frame()
  
  ### Function to standardization and TIX 
  dataP <- makeTIX(dataP,traits_TIX)
  
  ## Append traits character 
  dataP <- dataP |> 
    left_join(traits_app,
              by = join_by(KEY_ID == KEY_ID))
  
  ### Save results 
  wb <- createWorkbook()
  addWorksheet(wb, "BLUP_wide")
  writeData(wb, 1, x =dataP, withFilter = TRUE)
  freezePane(wb, "BLUP_wide", firstActiveRow = 2, firstActiveCol = 2)
  saveWorkbook(wb, file = out_file, overwrite = TRUE)
}


#' rogers.wright 
#' 
#' @title 
#' @param name description
#' @param name description
#' @details SASA
#' @return description 
#'  
#' por Alexandre S. G. Coelho 12/06/2019
#' função para cálculo da distância de Rogers-Wright (1978)
rogers.wright <- function(g){
  tabela <- data.frame(locus=g@loc.fac, t(g@tab))
  d <- matrix(NA, nrow=ncol(tabela)-1, ncol=ncol(tabela)-1)
  pt <- lapply(1:length(levels(tabela$locus)), function(l)
    t(t(subset(tabela, locus==levels(locus)[l])[, -1])/colSums(subset(tabela, locus==levels(locus)[l])[, -1], na.rm=TRUE)))
  p <- matrix(do.call(rbind, pt), nrow=ncol(g@tab), ncol=nrow(g@tab))
  for (i in 1:ncol(p)){
    for (j in 1:ncol(p)){
      d[j, i] <- sqrt(sum((p[,i]-p[,j])^2, na.rm=TRUE)/(2*length(levels(tabela$locus))))
    }
  }
  D <- as.dist(d)
  labs <- row.names(g@tab)
  attr(D, 'Labels') <- labs
  return(D)
}


#'  
#' @title Export files from the TESS2.3 and STRUCTURE format to the ".geno" and ".lfmm" formats. 
#' @param name description
#' @param name description
#' @param file = a file in the STRUCTURE or TESS format
#' @param output.format = a character string. Allowed values are "geno" and "lfmm" 
#' @param TESS = TRUE if the TESS 2.3 format is used (geographic coordinates (Lon,Lat) are binded to the left of the genotypic matrix).
#' @param diploid = TRUE for diploids, = FALSE for haploids 
#' @param FORMAT = 1 for markers encoded using one row of data for each individual
#' @param FORMAT = 2 for markers encoded using two rows of data for each individual
#' @param extra.row = an integer indicating the number of extra rows in the input file (e.g., marker ids)
#' @param extra.col = an integer indicating the number of extra columns in the input file (e.g., individual ids, pop ids, phenotypes, etc)
#' @param geographic coordinates must be considered as extra columns if the flag TESS = FALSE   
#' @details Missing data are encoded as "-9" or any negative values
#' @return Output files in the ".geno" or ".lfmm" formats. Geographic coordinates in a separate file. 
#' por Alexandre S. G. Coelho 12/06/2019
#' função para cálculo da distância de Rogers-Wright (1978)
struct2geno <- function(file = NULL, 
                        output.format = "geno", 
                        TESS = FALSE, 
                        diploid = TRUE,
                        FORMAT = 1, 
                        extra.row = 0, 
                        extra.col = 0, 
                        output = "genotype.geno"){
  if (!diploid & FORMAT == 2) stop("FORMAT = 2 is for diploids only")
  dat = read.table(file)
  if (TESS == FALSE){ 
    if (extra.row > 0) dat = dat[-(1:extra.row),]
    if (extra.col > 0) dat = dat[,-(1:extra.col)]
    n = dim(dat)[1]
    L = dim(dat)[2]
    if (FORMAT == 1 & diploid == FALSE) {n.ind = n; n.loc = L}
    if (FORMAT == 1 & diploid == TRUE) {n.ind = n; n.loc = L/2}
    if (FORMAT == 2 & diploid == TRUE) {n.ind = n/2; n.loc = L}    
    cat("Input file in the STRUCTURE format. The genotypic matrix has", n.ind, "individuals and", n.loc,
        "markers.","\n")
    cat("The number of extra rows is", extra.row,"and the number of extra columns is",extra.col,".\n")
  }
  if (TESS == TRUE){ 
    if (extra.row > 0) dat = dat[-(1:extra.row),]
    if (extra.col > 0) dat = dat[,-(1:extra.col)]
    n = dim(dat)[1]
    L = dim(dat)[2]
    if (FORMAT == 1 & diploid == FALSE) {
      n.ind = n; n.loc = L-2;
      coord = dat[,1:2]
      dat = dat[,-(1:2)]
      write.table(file = "coordinates.coord", coord, quote = F, row.names = F, col.names = F)
    }
    if (FORMAT == 1 & diploid == TRUE) {
      n.ind = n; n.loc = L/2 - 1; 
      coord = dat[,1:2]
      dat = dat[,-(1:2)]
      write.table(file = "coordinates.coord", coord, quote = F, row.names = F, col.names = F)      
    }
    if (FORMAT == 2 & diploid == TRUE) {
      n.ind = n/2; n.loc = L - 2;
      coord = dat[seq(1,n,by = 2), 1:2]
      dat = dat[,-(1:2)]
      write.table(file = "coordinates.coord", coord, quote = F, row.names = F, col.names = F)       
    }    
    cat("Input file in the TESS format. The genotypic matrix has", n.ind, "individuals and", n.loc,
        "markers.","\n")
    cat("The number of extra rows is", extra.row, 
        "and the number of extra columns is",extra.col,".\n")
  }
  dat = as.matrix(dat)
  unique.dat = unique(as.numeric(dat))
  missing.dat = unique.dat[unique.dat < 0]
  if (length(missing.dat) == 0)  cat("The input file contains no missing genotypes","\n")
  if (length(missing.dat) == 1)  cat("Missing alleles are encoded as",missing.dat,"\n")
  if (length(missing.dat) > 1) stop("Multiple values for missing data","\n") 
  # Convert allelic data into absence/presence data at each locus
  # Results are stored in the "dat.binary" object
  L = dim(dat)[2]
  if (FORMAT == 1 & diploid == FALSE) {
    dat.binary = NULL
    for (j in 1:L){
      allele = sort(unique(dat[,j]))
      for (i in allele[allele >= 0]) dat.binary=cbind(dat.binary, dat[,j]== i )
      LL = dim(dat.binary)[2]
      ind = which(dat[,j] < 0)
      if (length(ind) != 0){dat.binary[ind, (LL - length(allele) + 2):LL] = -9}
    }}
  if (FORMAT == 1 & diploid == TRUE) {
    dat.2 = matrix(NA, ncol = L/2, nrow = 2*n)
    for (ii in 1:n){
      dat.2[2*ii-1,] = dat[ii, seq(1,L,by = 2)]
      dat.2[2*ii,] = dat[ii,seq(2,L,by = 2)]
    }
    L = dim(dat.2)[2]
    
    dat.binary = NULL
    
    for (j in 1:L){
      allele = sort(unique(dat.2[,j]))
      for (i in allele[allele >= 0]) dat.binary=cbind(dat.binary, dat.2[,j]==i)
      LL = dim(dat.binary)[2]
      ind = which(dat.2[,j] < 0)
      if (length(ind) != 0){dat.binary[ind, (LL - length(allele) + 2):LL] = -9}
    }}
  if (FORMAT == 2 & diploid == TRUE) {
    dat.binary = NULL
    for (j in 1:L){
      allele = sort(unique(dat[,j]))
      for (i in allele[allele >= 0]) dat.binary=cbind(dat.binary, dat[,j]==i)
      LL = dim(dat.binary)[2]
      ind = which(dat[,j] < 0)
      if (length(ind) != 0){dat.binary[ind, (LL - length(allele) + 2):LL] = -9}
    }}
  # Compute a genotype count for each allele (0,1,2 or 9 for a missing value)
  # The results are stored in 'genotype'  
  n = dim(dat.binary)[1]
  if (diploid == TRUE){
    n = n/2
    genotype = matrix(NA,nrow=n,ncol=dim(dat.binary)[2])
    for(i in 1:n){
      genotype[i,]= dat.binary[2*i-1,]+dat.binary[2*i,]
      genotype[i, (genotype[i,] < 0)] = 9
    }}
  if (FORMAT == 1 & diploid == FALSE){
    genotype = dat.binary
    for(i in 1:n){
      genotype[i, (genotype[i,] < 0)] = 9
    }}  
  cat("The output file is",output, "\n") 
  if (output.format == "geno"){
    # Export to the "geno" format
    write.table(file = output,t(genotype), row.names=F,col.names=F,quote=F, sep = "")
  } else {
    # Export to the "lfmm" format
    write.table(file = output, genotype, row.names=F,col.names=F,quote=F)
  }
  
}

#' rogers.wright 
#' 
#' @title 
#' @param name description
#' @param name description
#' @details SASA
#' @return description 
#'  
#' por Alexandre S. G. Coelho 12/06/2019
fst <- function(project,run = 1, K, ploidy = 2){
  library(LEA)
  ll = dim(G(project, K = K, run = run))[1]
  if (ploidy == 2) {freq = G(project, K = K, run = run)[seq(2,ll,by = 3),]/2 + G(project, K = K, run = run)[seq(3,ll,by = 3),] } 
  else {freq = G(project, K = K, run = run)[seq(2,ll,by = 2),]}
  q = apply(Q(project, K = K, run = run), MARGIN = 2, mean)
  H.s = apply(freq*(1-freq), MARGIN = 1, FUN = function(x) sum(q*x) )
  P.t = apply(freq, MARGIN = 1, FUN = function(x) sum(q*x) )
  return(1-H.s/P.t/(1-P.t))
}


#' Source genotype data by inventory 
#' 
#' @title 
#' @param name description
#' @param name description
#' @details SASA
#' @return description 
#' 

source_genotype_ByINV_old <- function(fts_germplasm_id){
  project <- "bcs-veg-rd" 
  sql <-  paste0("
WITH 
## Source all LID by fts_germplasm_id
DBLID AS (
  SELECT cast(inv.fts_inventory_id as STRING) AS fts_inventory_id, vg.fts_germplasm_id, vg.germplasm_lineage_identifier
          FROM  `bcs-veg-mart.pipeline.veg_germplasms` vg,
          UNNEST(inventories) inv
          where vg.germplasm_lineage_identifier IN (
            SELECT vg.germplasm_lineage_identifier
            FROM  `bcs-veg-mart.pipeline.veg_germplasms` vg
            WHERE
            vg.fts_germplasm_id IN (",paste0(fts_germplasm_id,collapse=","),")
            )),
## Source markers by all inventories from same LID
DB AS (
    SELECT name,
    sub.type AS type,
    CAST(sub.id AS STRING) AS fts_inventory_id,
    cal.marker_name AS marker_name,
    cal.allele_calls[OFFSET(0)] || cal.allele_calls[OFFSET(1)] AS score
    FROM
      `product360-datasets.genomics.MarkerCallSetsV2`,
      UNNEST(subjects) AS sub,
      UNNEST(calls) AS cal
    WHERE
      name IN (
      SELECT name FROM `product360-datasets.genomics.MarkerCallSetsV2` mc,
      UNNEST(subjects) AS sub
      WHERE
        sub.id IN (SELECT fts_inventory_id FROM DBLID))
      AND CHAR_LENGTH( cal.allele_calls[OFFSET(0)] || cal.allele_calls[OFFSET(1)]) = 2 ),
## pivot sample and inventory id
DB1 AS (
    SELECT cast(inventory as STRING) as inventory,sample,marker_name,score FROM DB
     PIVOT( MAX(fts_inventory_id) FOR type IN ('sample','inventory')))
## join markers data and LID information 
SELECT
    DBLID.germplasm_lineage_identifier,
    DBLID.fts_germplasm_id, 
    ## DB1.inventory,
    ## sample,
    marker_name,
    score FROM DB1
  LEFT JOIN DBLID ON DB1.inventory = DBLID.fts_inventory_id
  ;",sep='') 
  
  tb <- bq_project_query(project, sql)
  df_gen <- bq_table_download(tb,bigint="integer64",page_size = 20000) 
  
  df_gen$score1 <- converte(df_gen$score)
  df_gen <- df_gen |> dplyr::filter(!is.na(score1))
  
  fts_min <- df_gen |> 
    group_by(germplasm_lineage_identifier) |> 
    summarise(fts_min=min(fts_germplasm_id,na.rm=TRUE))
  
  df_gen <- df_gen |> 
    left_join(fts_min,
              by=join_by(germplasm_lineage_identifier == germplasm_lineage_identifier),
              relationship="many-to-one",
              multiple="first")
  
  df_gen <- df_gen |>  
    pivot_wider(id_cols = fts_min, 
                names_from = marker_name, 
                values_from = score1,
                values_fill =NA,
                values_fn=max)  |> 
    as.data.frame()
  
  colnames(df_gen) <- gsub('\\.', '_', colnames(df_gen))
  rownames(df_gen) <- df_gen[,1]
  
  df_gen <- df_gen[,apply(df_gen,2,function(x){sum(!is.na(x))}) >= nrow(df_gen)*0.80]
  df_gen <- df_gen[,-1]
  
  for (i in 1:ncol(df_gen)){
    x <- data.frame(m=(as.character(df_gen[,i])))
    scorS <- table(df_gen[,i])
    if(nrow(scorS)>1){
      scorS <- data.frame(sort(scorS,decreasing = TRUE))
      scorS$homo <- substr(scorS$Var1,1,1) ==  substr(scorS$Var1,2,2)
      scorS$score <- NA
      scorS$score[scorS$homo == FALSE] <- 1
      scorS$score[scorS$homo == TRUE] <- c(2,0)
      scorS$Var1 <- as.character(scorS$Var1)
      x$score <- NA
      for(k in 1:nrow(scorS)){
        x$score[x$m==scorS$Var1[k]] <- scorS$score[k]
      }
      df_gen[,i] <- x$score
    }else{
      df_gen[,i] <- 2
    }
  }
  df_gen <- as.matrix(df_gen[,-1])
  return(df_gen)
}


#' Source genotype data by inventory 
#' 
#' @title 
#' @param name description
#' @param name description
#' @details SASA
#' @return description 
#' 

source_genotype_ByINV <- function(fts_germplasm_id){
  project <- "bcs-veg-rd" 
  sql <-  paste0("
WITH 
## Source all LID by fts_germplasm_id
DBLID AS (
  SELECT cast(inv.fts_inventory_id as STRING) AS fts_inventory_id, vg.fts_germplasm_id, vg.germplasm_lineage_identifier
          FROM  `bcs-veg-mart.pipeline.veg_germplasms` vg,
          UNNEST(inventories) inv
          where vg.germplasm_lineage_identifier IN (
            SELECT vg.germplasm_lineage_identifier
            FROM  `bcs-veg-mart.pipeline.veg_germplasms` vg
            WHERE
            vg.fts_germplasm_id IN (",paste0(fts_germplasm_id,collapse=","),")
            )),
## Source markers by all inventories from same LID
DB AS (
    SELECT name,
    sub.type AS type,
    CAST(sub.id AS STRING) AS fts_inventory_id,
    cal.marker_name AS marker_name,
    cal.allele_calls[OFFSET(0)] || cal.allele_calls[OFFSET(1)] AS score
    FROM
      `product360-datasets.genomics.MarkerCallSetsV2`,
      UNNEST(subjects) AS sub,
      UNNEST(calls) AS cal
    WHERE
      name IN (
      SELECT name FROM `product360-datasets.genomics.MarkerCallSetsV2` mc,
      UNNEST(subjects) AS sub
      WHERE
        sub.id IN (SELECT fts_inventory_id FROM DBLID))
      AND CHAR_LENGTH( cal.allele_calls[OFFSET(0)] || cal.allele_calls[OFFSET(1)]) = 2),
## pivot sample and inventory id
DB1 AS (
    SELECT cast(inventory as STRING) as inventory,sample,marker_name,score FROM DB
     PIVOT( MAX(fts_inventory_id) FOR type IN ('sample','inventory')))
## join markers data and LID information 
SELECT
    DBLID.germplasm_lineage_identifier,
    DBLID.fts_germplasm_id, 
    ## DB1.inventory,
    ## sample,
    marker_name,
    score FROM DB1
  LEFT JOIN DBLID ON DB1.inventory = DBLID.fts_inventory_id
  ;",sep='') 
  
  tb <- bq_project_query(project, sql)
  df_gen <- bq_table_download(tb,bigint="integer64",page_size = 20000) 
  
  df_gen <- df_gen |> 
    group_by(germplasm_lineage_identifier) |> 
    dplyr::mutate(min_fts_germplasm_id = min(fts_germplasm_id,na.rm=TRUE)) |> 
    ungroup() |> 
    dplyr::mutate(score = gsub('[*]|[R]',NA,score)) |> 
    pivot_wider(id_cols = min_fts_germplasm_id, 
                names_from = marker_name, 
                values_from = score,
                values_fill =NA,
                values_fn=max) |> 
    as.data.frame()
  
  colnames(df_gen) <- gsub('\\.', '_', colnames(df_gen))
  rownames(df_gen) <- df_gen[,1]
  df_gen <- df_gen[,-1]
  return(df_gen)
}

#' Split traits of characters
#' 
#' @title 
#' @param name description
#' @param name description
#' @details SASA
#' @return description 
#' 
extract_string_traits <- function(data, traits_char){
  data <- data |> 
    dplyr::filter(TRAIT %in% traits_char) |> 
    makeKEY_ID() |> 
    dplyr::select(KEY_ID,observation_trait_val,TRAIT,field_name, season)
  if(nrow(data)>1){
    results_cmnts <- list()
    con <- 1
    for (i in unique(data$TRAIT)){
      data_i <- data |> dplyr::filter(TRAIT == i)
      for (j in unique(data_i$KEY_ID)){
        sel <- data_i$KEY_ID == j
        results_cmnts[[con]] <- data.frame(TRAIT = i,
                                           KEY_ID = j,
                                           VAR=paste(apply(data_i[sel,c('field_name','observation_trait_val')],1,function(x){paste(x,collapse = ':')}),collapse = " | "))
        con <- con+1
      }
    }
    results_cmnts <- plyr::ldply(results_cmnts) |> 
      pivot_wider(id_cols = KEY_ID,
                  names_from = TRAIT,
                  values_from = VAR)
    write.csv(results_cmnts,file='input_files/out_Traits_Char.csv',row.names = FALSE)
  }
  return()
}


#' Source genotype data by inventory 
#' 
#' @title 
#' @param name description
#' @param name description
#' @details SASA
#' @return description 
#' 

source_all_fts_byLID <- function(fts_germplasm_id){
  project <- "bcs-veg-rd" 
  sql <-  paste0("
  SELECT
  CAST(vg.fts_germplasm_id AS INT64) AS fts_germplasm_id, 
  vg.cross_name,
  vg.pedigree,
  vg.germplasm_lineage_identifier
          FROM  `bcs-veg-mart.pipeline.veg_germplasms` vg
          ### UNNEST(inventories) inv
          where vg.germplasm_lineage_identifier IN (
            SELECT vg.germplasm_lineage_identifier
            FROM  `bcs-veg-mart.pipeline.veg_germplasms` vg
            WHERE
            vg.fts_germplasm_id IN (",paste0(fts_germplasm_id,collapse=","),")
            )
  ;",sep='') 
  
  tb <- bq_project_query(project, sql)
  df_id <- bq_table_download(tb,bigint="integer64",page_size = 20000) 
  
  return(df_id)
}

#' Markers data quality 
#' 
#' @title 
#' @param name description
#' @param name description
#' @details SASA
#' @return description 
#' 
SNP_DQ <- function(data = df_gen,min_markers=50,per_min =0.70){
  print(paste0("Initial data: ",dim(data)[1]," Genotypes ",dim(data)[2]," Markers"))
  ## Filter by marker number and SNP data available
  MarkerNumber <- apply(data,1,function(x){sum(!is.na(x))})
  data <- data[MarkerNumber >= min_markers,]
  print(paste0("After remove by minimal markers: ",dim(data)[1]," Genotypes ",dim(data)[2]," Markers"))
  Nrow <- nrow(data)
  data <- data[,apply(data,2,function(x){sum(!is.na(x))}) >= Nrow * per_min] 
  print(paste0("After remove by missing data: ",dim(data)[1]," Genotypes ",dim(data)[2]," Markers"))
  ## Remove SNP with no variation
  #data <- data[,apply(data,2,function(x){var(na.omit(x))})>0]
  #print(paste0("After remove by no SNP variation: ",dim(data)[1]," Genotypes ",dim(data)[2]," Markers"))
  return(data)
}


#' predict hybrid markers
#' 
#' @title 
#' @param name description
#' @param name description
#' @details SASA
#' @return description 
#' 
pred_hyb <- function(
    M = as.matrix(Ma),
    ped = df_pedHyb,
    indiv = "ID",
    mother = "SIRE",
    father = "DAM"){
  ped <- ped[(ped[,mother] %in% row.names(M)) & (ped[,father] %in% row.names(M)),]
  Mpred <- matrix(NA,ncol=ncol(M),nrow=nrow(ped))
  colnames(Mpred) <- colnames(M)
  rownames(Mpred) <- rownames(ped)
  ROWnames <- rownames(M)
  for(i in 1:nrow(Mpred)){
    Mpred[i,] <- colMeans(M[ROWnames %in% ped[i,c(mother,father)],])
  }
  return(Mpred)
}


#' synthetic.crossXXX
#' 
#' @title 
#' @param name description
#' @param name description
#' @details SASA
#' @return description 
#' 
synthetic.crossXXX <- function (M = NULL, ped = NULL, indiv = NULL, mother = NULL, 
                                father = NULL, heterozygote.action = c("useNA", "exact", 
                                                                       "fail", "expected"), na.action = c("useNA", "expected"), 
                                message = TRUE) 
{
  na.action <- match.arg(na.action)
  heterozygote.action <- match.arg(heterozygote.action)
  check.data_(data_ = "M", class_ = "matrix")
  check.data_(data_ = "ped", class_ = "data.frame")
  check.args_(data_ = "ped", mandatory_ = TRUE, arg_ = "indiv", 
              class_ = "character")
  check.args_(data_ = "ped", mandatory_ = TRUE, arg_ = "mother", 
              class_ = "character")
  check.args_(data_ = "ped", mandatory_ = TRUE, arg_ = "father", 
              class_ = "character")
  heterozygote.action <- match.arg(heterozygote.action)
  na.action <- match.arg(na.action)
  check.logical_(arg_ = "message")
  if (heterozygote.action == "fail") {
    het.markers <- apply(X = M, MARGIN = 2, FUN = function(m) {
      any(m == 1, na.rm = TRUE)
    })
    if (any(het.markers)) {
      stop("Stop requested as some of the markers have have heterozygous states (1), e.g., ", 
           paste0(names(head(het.markers[het.markers], 5)), 
                  collapse = ", "), "...")
    }
  }
  if (!all(c(ped[[mother]], ped[[father]]) %in% rownames(M))) {
    stop("Some parents do not have genotypic information in matrix M.")
  }
  unique.comb <- paste(ped[[mother]], ped[[father]], sep = "_")
  unique.comb.rec <- paste(ped[[father]], ped[[mother]], sep = "_")
  recs <- unique.comb %in% unique.comb.rec
  dups <- duplicated(unique.comb)
  if (any(dups)) {
    message("A total of ", sum(dups), " duplicated rows were found in the supplied pedigree.")
    message("Removing 'ped' lines: ", which(dups), ".")
    ped <- ped[!dups, ]
  }
  if (any(recs)) {
    warning("Reciprocals found in the supplied pedigree. These were not removed.")
  }
  rm(dups, recs)
  initial.n.markers <- ncol(M)
  if (heterozygote.action == "exact") {
    if (message) {
      message(blue("\nExact method selected. Eliminating markers containing one or more heterozygotic read."))
    }
    het.markers <- apply(X = M, MARGIN = 2, FUN = function(m) {
      any(m == 1, na.rm = TRUE)
    })
    M <- M[, !het.markers, drop = FALSE]
    if (message) {
      message("Total number of dropped markers: ", initial.n.markers - 
                ncol(M))
      message("Total number of remaining markers: ", ncol(M))
    }
    if (ncol(M) == 0) {
      stop("No marker were left after removal of heterozygote reads.")
    }
  }
  marker.ids <- colnames(M)
  range.hybrids <- 1:nrow(ped)
  if (message) {
    message(blue("\nGenerating in hypothetical crosses genotypic information."))
  }
  get.off.gen <- function(cur.offspring = NULL) {
    M.mother <- M[ped[cur.offspring, mother], ]
    M.father <- M[ped[cur.offspring, father], ]
    M.offspring <- rbind(M.mother, M.father)
    if (heterozygote.action == "useNA") {
      M.offspring[M.offspring %in% 1] <- NaN
    }
    if (na.action == "useNA") {
      return(colMeans(M.offspring))
    }
    if (na.action == "expected") {
      evo <- NULL
      kid <- colMeans(M.offspring)
      missing.markers <- is.na(kid) & !is.nan(kid)
      if (any(missing.markers)) {
        freq2 <- colMeans(M, na.rm = T)/2
        evo <- sapply(X = which(missing.markers), function(m) {
          marker <- M.offspring[, m]
          f.pseudo.major <- freq2[m]
          f.pseudo.minor <- (1 - f.pseudo.major)
          if (sum(is.na(marker)) == 1) {
            par.gen.no.miss <- sum(marker, na.rm = TRUE)
            if (par.gen.no.miss == 0) {
              evo <- f.pseudo.minor^2 * 0 + 2 * f.pseudo.minor * 
                f.pseudo.major * 0.5 + f.pseudo.major^2 * 
                1
            }
            if (par.gen.no.miss == 1) {
              evo <- f.pseudo.minor^2 * 0.5 + 2 * f.pseudo.minor * 
                f.pseudo.major * 1 + f.pseudo.major^2 * 
                1.5
            }
            if (par.gen.no.miss == 2) {
              evo <- f.pseudo.minor^2 * 1 + 2 * f.pseudo.minor * 
                f.pseudo.major * 1.5 + f.pseudo.major^2 * 
                2
            }
          }
          if (sum(is.na(marker)) == 2) {
            f.pseudo.major <- freq2[m]
            evo <- f.pseudo.minor^2 * f.pseudo.minor^2 * 
              0 + 2 * (f.pseudo.minor^2 * 2 * f.pseudo.minor * 
                         f.pseudo.major) * 0.5 + 2 * (f.pseudo.minor^2 * 
                                                        f.pseudo.major^2) * 1 + 2 * f.pseudo.minor * 
              f.pseudo.major * 2 * f.pseudo.minor * f.pseudo.major * 
              1 + 2 * (2 * f.pseudo.minor * f.pseudo.major * 
                         f.pseudo.major^2) * 1.5 + f.pseudo.major^2 * 
              f.pseudo.major^2 * 2
          }
          return(evo)
        })
      }
      kid[which(missing.markers)] <- evo
      return(kid)
    }
  }
  M <- lapply(X = range.hybrids, FUN = get.off.gen)
  M <- do.call(rbind, M)
  M[is.nan(M)] <- NA
  rownames(M) <- ped[[indiv]]
  colnames(M) <- marker.ids
  return(M)
}


#' Source rawdata from CSW to Hybrid Analysis Pipeline
#' 
#' @title 
#' @param WHERE filter to apply in WHERE SQL query from CSW in table bcs-veg-mart.performance.veg_experiment_obs_entries
#' @details the columns mandatory for query is already fixed
#' @return data frame with raw data with no treatment from CSW
#' 
source_raw_data <- function(WHERE=NULL){
  project <- "bcs-veg-rd"
  sql <- paste("SELECT
  distinct
  veoe.trial_type,
  veoe.season,
  veoe.country, 
  veoe.field_name, 
  veoe.planting_year,
  veoe.planting_date,
  COALESCE(veoe.block_name,veoe.set_name) AS block_name,
  ### veoe.experiment_stage_cd,
  CAST(veoe.plot_id AS INT64) AS plot_id, 
  CAST(veoe.fts_germplasm_id AS INT64) AS fts_germplasm_id,
  veoe.pedigree,
  veoe.origin,
  veoe.germplasm_lineage_identifier, 
  veoe.cross_name,veoe.observation_id,
  prod.primary_commercial_product_nm as PROD, 
  veoe.line_type, veoe.observation_reference_cd,
  veoe.descriptor_abbreviation_cd,
  veoe.observation_trait_val,
  
  COALESCE(primary_question_cd,observation_reference_cd) AS TRAIT,
  COALESCE(secondary_answers_concat,descriptor_abbreviation_cd) AS TRAIT_ABBREVIATION,
  primary_question_cd,
  primary_question_cd_desc,
  secondary_answers_concat
FROM
  `bcs-veg-mart.performance.veg_experiment_obs_entries` veoe
  left join (`bcs-veg-mart.pipeline.veg_germplasms` vg CROSS JOIN UNNEST (products) as prod) on veoe.fts_germplasm_id = vg.fts_germplasm_id
WHERE
  trial_type IN ('Yield Trials','Research Trial')
  AND veoe.observation_trait_val IS NOT NULL
  AND veoe.is_plot_deactivated = false
  AND ",WHERE,";")
  tb <- bq_project_query(project, sql)
  df_raw <- bq_table_download(tb,bigint="integer64",page_size = 20000) |> as.data.frame()
  print(paste('DataFrame dimension: ',dim(df_raw)[1],' rows'))
  return(df_raw)
}



#' Source rawdata from CSW to Hybrid Analysis Pipeline
#' 
#' @title 
#' @param WHERE filter to apply in WHERE SQL query from CSW in table bcs-veg-mart.performance.veg_experiment_obs_entries
#' @details the columns mandatory for query is already fixed
#' @return data frame with raw data with no treatment from CSW
#' 

formatSeg <- function(x){
  return(sprintf("%02d:%02d:%02d",
                 x %% 86400 %/% 3600,
                 x %% 3600 %/% 60,  
                 x %% 60 %/% 1))
}

#' Source rawdata from CSW to Hybrid Analysis Pipeline
#' 
#' @title 
#' @param WHERE filter to apply in WHERE SQL query from CSW in table bcs-veg-mart.performance.veg_experiment_obs_entries
#' @details the columns mandatory for query is already fixed
#' @return data frame with raw data with no treatment from CSW
#' 
run_model <- function(
    df_raw = NULL,
    minGen = NULL,
    UseGenomeData = F, #only for Corn tested
    SNP_source = c('raw','imputed')[2],
    Inb_hyb_append = NULL,
    Use10genPed = NULL,
    Exclude_Hybrid_GCA = NULL,
    ##crop = NULL,
    file_out = NULL,
    debug_folder = 'out_DEBUG',
    traits_append = NULL,
    num_ParRun = 1,
    Zscore_filter = 4){
  
  # Adjust right functions 
  select <- dplyr::select
  filter <- dplyr::filter
  foreach <- foreach::foreach
  
  # Adjust folder to debug
  if(file.exists(file.path(getwd(),debug_folder))){
    unlink(file.path(getwd(),debug_folder),recursive=TRUE)
    dir.create(file.path(getwd(),debug_folder))
  }else{dir.create(file.path(getwd(),debug_folder))}
  
  # Adjust data issues expand tratis with sublevel 
  df_raw <- df_raw |> dplyr::filter(TRAIT != 'CMNTS')
  df_raw <- check_data_adj(data=df_raw)
  df_raw <- makeKEY_ID(data=df_raw)
  KEYID_backup <- unique(df_raw$KEY_ID)
  
  ## Summary data 
  df_rawS <- df_raw |> 
    group_by(TRAIT) |> 
    summarise(min=min(observation_trait_val,na.rm=TRUE),
              Q05=quantile(observation_trait_val,0.05),
              Q25=quantile(observation_trait_val,0.25),
              mean=mean(observation_trait_val,na.rm=TRUE),
              Q75=quantile(observation_trait_val,0.75),
              Q95=quantile(observation_trait_val,0.95),
              max=max(observation_trait_val,na.rm=TRUE),
              n_genotypes=dplyr::n_distinct(fts_germplasm_id),
              n=dplyr::n()) |> 
    dplyr::mutate(across(min:max,function(x) round(x,2))) |> 
    as.data.frame()
  print(df_rawS)
  
  ## Remove traits with low genotypes
  df_raw <- df_raw |> 
    dplyr::filter(TRAIT %in%
                    df_rawS$TRAIT[df_rawS$n_genotypes > minGen]) |> 
    droplevels()
  
  ## data quality plot 
  df_raw |> 
    ggplot(aes(x=observation_trait_val,fill=TRAIT))+
    geom_density()+
    guides(fill = FALSE)+
    facet_wrap(~TRAIT, scales='free')+
    theme_bw() 
  ggsave(paste0(debug_folder,'\\out_DQ.png'),w=15,h=15)
  
  ## Additional inbred and Hybrids load all fts_germplasm_ID for additional germplasm
  if(!is.null(Inb_hyb_append)){
    df_append <- source_ped(cross_name=Inb_hyb_append$germplasm_lineage_identifier,debug_folder = debug_folder)
    df_append <- makeKEY_ID(df_append)
    
    not_found_LID <- df_append[is.na(df_append$germplasm_lineage_identifier),]
    if(nrow(not_found_LID)>1){
      write.csv(not_found_LID,paste0(debug_folder,'\\out_InbredAppend_notFoundBy_LID.csv'),row.names = FALSE)
    }
    df_append <-  df_append |> 
      dplyr::select(fts_germplasm_id,line_type,KEY_ID) |> 
      unique() |> 
      as.data.frame()
    ## Join with data pedigree 
    ped_all <- rbind(df_append,df_raw[,c('fts_germplasm_id','line_type', 'KEY_ID')]) |>
      unique()
  }else{   ped_all <- df_raw[,c('fts_germplasm_id','line_type', 'KEY_ID')] |> unique()  }
  
  # Loading Pedigree 
  ## Source deep pedigree 5 generation or more 
  df_ped <- source_deep_ped2(fts_germplasm_id=unique(ped_all$fts_germplasm_id))
  if(nrow(df_ped) == 0){
    stop('No deep pedigre are returned in this analysis, check you selection')
  }
  if(Use10genPed){
    df_ped <- source_deep_ped2(
      unique(na.omit(c(df_ped$child_germplasm_id,
                       df_ped$male_parent_germplasm_id,
                       df_ped$female_parent_germplasm_id)))) |> 
      arrange(child_germplasm_id)
    print(dim(df_ped))
  }
  
  df_ped <- df_ped |> 
    unique() |> 
    arrange(child_germplasm_id,male_parent_germplasm_id,female_parent_germplasm_id)
  
  #df_ped |> View()
  ## Master ID file
  df_ID <- source_all_ID_ByOrigin_ID(fts_germplasm_id = 
                                       unique(c(df_ped$child_germplasm_id,
                                                df_ped$male_parent_germplasm_id,
                                                df_ped$female_parent_germplasm_id,
                                                df_raw$fts_germplasm_id))) 
  
  df_IDmissing <- df_raw[!(df_raw$fts_germplasm_id %in% df_ID$fts_germplasm_id),] |> 
    dplyr::select(c('fts_germplasm_id','line_type','pedigree','origin','cross_name','germplasm_lineage_identifier','PROD')) |> 
    unique() |> 
    dplyr::mutate(male_parent_germplasm_id = NA,
                  female_parent_germplasm_id = NA) |> 
    droplevels()
  
  df_ID <- rbind(df_ID,df_IDmissing) |> 
    makeKEY_ID() |> 
    dplyr::filter(!is.na(KEY_ID)) |> 
    group_by(KEY_ID) |> 
    dplyr::mutate(MIN_fts_germplasm_id = min(fts_germplasm_id),
                  rankID = dense_rank(fts_germplasm_id )) |> 
    ungroup() |> 
    arrange(KEY_ID,fts_germplasm_id) |> 
    as.data.frame()
  
  df_ID_male <- df_ID |> dplyr::select(fts_germplasm_id,MIN_fts_germplasm_id) |> dplyr::rename(MIN_male_parent_germplasm_id = MIN_fts_germplasm_id)
  df_ID_female <- df_ID |> dplyr::select(fts_germplasm_id,MIN_fts_germplasm_id) |> dplyr::rename(MIN_female_parent_germplasm_id = MIN_fts_germplasm_id)
  
  df_ID <-  df_ID |> 
    left_join(df_ID_male,by=join_by(male_parent_germplasm_id == fts_germplasm_id),
              multiple = 'first',relationship = 'many-to-one') |> 
    left_join(df_ID_female, by=join_by(female_parent_germplasm_id == fts_germplasm_id),
              multiple = 'first',relationship = 'many-to-one') |> 
    dplyr::mutate(MIN_male_parent_germplasm_id =
                    case_when(is.na(MIN_male_parent_germplasm_id) ~
                                male_parent_germplasm_id,.default = MIN_male_parent_germplasm_id),
                  MIN_female_parent_germplasm_id =
                    case_when(is.na(MIN_female_parent_germplasm_id) ~
                                female_parent_germplasm_id,.default = MIN_female_parent_germplasm_id)) |>  
    dplyr::select(-c(male_parent_germplasm_id, female_parent_germplasm_id))
  
  #df_ID |> View()
  
  df_ID_min <- df_ID |> dplyr::select(fts_germplasm_id,MIN_fts_germplasm_id)
  df_ID_min2 <- df_ID |> dplyr::filter(rankID == 1) |> dplyr::select(MIN_fts_germplasm_id,MIN_male_parent_germplasm_id, MIN_female_parent_germplasm_id) |> unique()
  
  df_ped <- df_ped |> 
    dplyr::select(child_germplasm_id) |> 
    left_join(df_ID_min,by=join_by(child_germplasm_id == fts_germplasm_id),
              multiple = 'first',relationship = 'many-to-one') |> 
    left_join(df_ID_min2,by=join_by(MIN_fts_germplasm_id == MIN_fts_germplasm_id),
              multiple = 'first',relationship = 'many-to-one') #|>  dplyr::rename(MIN_male_parent_germplasm_id = male_parent_germplasm_id,         MIN_female_parent_germplasm_id = female_parent_germplasm_id)
  
  ## needs improvement!!!!!!!!!!! 
  df_ped[is.na(df_ped)] <- 0
  df_ped <- df_ped |> dplyr::filter(!(MIN_fts_germplasm_id == MIN_female_parent_germplasm_id))
  df_ped <- df_ped |> dplyr::filter(!(MIN_fts_germplasm_id == MIN_male_parent_germplasm_id))
  
  ## Replace with MIN FTS in df_raw 
  #df_ped|> View()
  df_raw01 <- df_raw |> dplyr::filter(fts_germplasm_id %in% df_ID$MIN_fts_germplasm_id) |> droplevels() #fts_ID already is minimal
  df_raw02 <- df_raw |> dplyr::filter(!(fts_germplasm_id %in% df_ID$fts_germplasm_id)) |> droplevels()  #fts_ID no include in df_ID
  df_raw03 <- df_raw |> dplyr::filter(fts_germplasm_id %in% df_ID$fts_germplasm_id & 
                                        !(df_raw$fts_germplasm_id %in% df_ID$MIN_fts_germplasm_id)) |> droplevels() #fts_ID ins't minimal
  
  df_raw03 <- df_raw03 |> dplyr::select(-c(fts_germplasm_id,
                                           pedigree,
                                           origin,
                                           germplasm_lineage_identifier,
                                           cross_name,
                                           PROD))
  df_ID_min3 <- df_ID |> dplyr::filter(rankID == 1) |> dplyr::select(-c(MIN_male_parent_germplasm_id, MIN_female_parent_germplasm_id))
  
  df_raw03 <- df_raw03 |> 
    dplyr::select(-c(line_type)) |> 
    left_join(df_ID_min3, by=join_by(KEY_ID == KEY_ID ),
              multiple = 'first',relationship = 'many-to-one') |>
    dplyr::select(-c(fts_germplasm_id,rankID)) |> 
    dplyr::rename(fts_germplasm_id= MIN_fts_germplasm_id)
  df_raw03 <- df_raw03[,match(colnames(df_raw03),colnames(df_raw01))]
  
  df_raw <- rbind(df_raw01,df_raw02,df_raw03) |>
    dplyr::filter(!is.na(fts_germplasm_id)) |> 
    droplevels() |> 
    dplyr::rename(MIN_fts_germplasm_id = fts_germplasm_id)
  
  #df_raw |> View()
  ## joint parents ID with data 
  df_raw <- df_raw |> 
    left_join(unique(df_ped[,c('MIN_fts_germplasm_id',
                               'MIN_male_parent_germplasm_id',
                               'MIN_female_parent_germplasm_id')]),
              by= join_by(MIN_fts_germplasm_id == MIN_fts_germplasm_id),
              multiple = 'first',relationship = 'many-to-one')
  
  #df_raw |> View()
  df_pedsel <- df_ped |> 
    dplyr::select(MIN_fts_germplasm_id,
                  MIN_male_parent_germplasm_id,
                  MIN_female_parent_germplasm_id) |> 
    unique() 
  
  ## find no parentage (checks) 
  ### imputed fake parents
  if(anyNA(c(df_raw$MIN_male_parent_germplasm_id,df_raw$MIN_female_parent_germplasm_id))){
    Noparentage <- imp_fake_parentage2(df_ped=df_ped,data=df_raw)
    for(i in 1:nrow(Noparentage)){
      sel <- Noparentage$MIN_fts_germplasm_id[i]
      sel1 <- df_raw$MIN_fts_germplasm_id == sel
      sel2 <- Noparentage$MIN_fts_germplasm_id == sel
      df_raw[sel1, 'MIN_male_parent_germplasm_id'] <- Noparentage[sel2 ,'MIN_male_parent_germplasm_id']
      df_raw[sel1, 'MIN_female_parent_germplasm_id'] <- Noparentage[sel2,'MIN_female_parent_germplasm_id']
    }
    df_pedsel <- rbind(Noparentage,df_pedsel)
  }
  
  # make pedigree data 
  pedl <- insertPedM(df_pedsel)
  pedl <- pedl[!(pedl$MIN_fts_germplasm_id == pedl$MIN_male_parent_germplasm_id),]
  pedl[is.na(pedl)] <- 0
  pedln <- renum(pedl)
  pedln$newped[pedln$newped == 0] <- NA
  dim(pedln$newped)
  
  ## Merge id numeric 
  pedln$xrf$ID <- as.integer64(pedln$xrf$ID)
  df_raw <- df_raw |> 
    left_join(pedln$xrf,by=join_by(MIN_fts_germplasm_id == ID),relationship = 'many-to-one') |>
    dplyr::rename(fts_germplasm_idN=newID) |>
    left_join(pedln$xrf,by=join_by(MIN_female_parent_germplasm_id == ID),relationship = 'many-to-one') |>
    dplyr::rename(female_parent_germplasm_idN=newID) |>
    left_join(pedln$xrf,by=join_by(MIN_male_parent_germplasm_id==ID),relationship = 'many-to-one') |>
    dplyr::rename(male_parent_germplasm_idN=newID)
  #df_raw |> View()
  
  # Genotype matrix
  ## using genome and pedigree data 
  if(UseGenomeData){
    #### make G additive matrix 
    if(SNP_source == 'imputed'){
      df_gen <- source_genotype_imp(fts_germplasm_id=pedln$xrf$ID)
      df_gen <- df_gen[as.integer64(rownames(df_gen)) %in% df_ID$fts_germplasm_id,]
      rownames(df_gen) <- as.character(df_ID$MIN_fts_germplasm_id[match(as.integer64(rownames(df_gen)),
                                                                        df_ID$fts_germplasm_id)])
      Ma <- df_gen
    }
    if(SNP_source == 'raw'){
      df_gen <- source_genotype_ByINV(fts_germplasm_id=pedln$xrf$ID)
      df_gen <- df_gen[as.integer64(rownames(df_gen)) %in% df_ID$fts_germplasm_id,]
      rownames(df_gen) <- as.character(df_ID$MIN_fts_germplasm_id[match(as.integer64(rownames(df_gen)),
                                                                        df_ID$fts_germplasm_id)])
      
      df_gen <- SNP_DQ(df_gen,min_markers=50,per_min = 0.70) 
      Ma <- snpReady::raw.data(data = as.matrix(df_gen),
                               frame = "wide",
                               base = TRUE,
                               sweep.sample = 0.5,
                               call.rate = 0.85,
                               maf = 0.01,
                               imput = TRUE,
                               imput.type = "wright",
                               outfile = "012")$M.clean
    }
    row.names(Ma) <- pedln$xrf$newID[match(rownames(Ma),pedln$xrf$ID)]
    df_pedHyb <- pedln$newped[pedln$newped$ID %in% (pedln$xrf$newID[pedln$xrf$ID %in% df_ID$MIN_fts_germplasm_id[df_ID$line_type == 'Hybrid']]),]
    M_hyb_pred <- pred_hyb(M = as.matrix(Ma),ped = df_pedHyb,indiv = "ID",mother = "SIRE",father = "DAM")
    # Get genotype of crosses and use expected values.
    # M_Hyb <- synthetic.cross( M = as.matrix(Ma[,1:100]),ped = df_pedHyb,indiv = "ID",mother = "SIRE",father = "DAM",heterozygote.action = "expect",na.action = "expect",message = TRUE)
    Mall <- rbind(Ma,M_hyb_pred)
    if(ncol(Mall)<= 30){stop("No genotype data sourced, please check!!")}
    Gm <- Gmatrix(Mall, method="VanRaden", ploidy=2,missingValue=NA,integer = FALSE)
    Gm <- Matrix::forceSymmetric(Gm)
    diag(Gm) <- diag(Gm)+1e-03
    #colnames(Gm) <- rownames(Gm) <- pedln$xrf$newID[match(colnames(Gm),pedln$xrf$ID)]
    #image(Gm)
    
    #### make A Additive matrix 
    pedln_newped0 <- pedln$newped
    pedln_newped0[is.na(pedln_newped0)] <- 0
    Am <- AGHmatrix::Amatrix(pedln_newped0, ploidy=2)
    Am <- Matrix::forceSymmetric(Am)
    
    #### make A inverse Additive matrix 
    Ai <- MCMCglmm::inverseA(pedln$newped)$Ainv
    Ai <- Matrix::forceSymmetric(Ai)
    
    #### make inverse H Additive matrix 
    Cov_m <- makeHinv(A=Am,Ai=Ai,G=Gm,tau=(1-0.1),omega=0.1) #90% G and 10%A
    diag(Cov_m) <- diag(Cov_m) + 1e-03
    if(det(Cov_m)<= 0){stop("Cov_m are singular, det(Cov_m) <= 0, please check!!")}
    rm(Gm,Am,Ai)
  }else{
    ## using only pedigree data 
    # make A inverse Additive matrix 
    Cov_m <- MCMCglmm::inverseA(pedln$newped)$Ainv
    Cov_m <- Matrix::forceSymmetric(Cov_m)
    colnames(Cov_m) <- rownames(Cov_m) <- pedln$newped$ID
  }
  
  rinf <- apply(Cov_m,2,function(x){(sum(x!=0)-1)}) ## number of relationship in CovMatrix
  rinf <- data.frame(ID=as.integer(names(rinf)),rinf=rinf)
  
  png(paste0(debug_folder,'\\out_COV_Matrix.png'),w=1000,h=1000)
  print(image(Cov_m))
  dev.off()
  
  # Capture Germplasm Characteristic 
  unique_ids <- as.integer64(unique(c(pedln$xrf$ID,df_raw$MIN_fts_germplasm_id,df_IDmissing$fts_germplasm_id)))
  df_add <- source_gc_append(fts_germplasm_id=unique_ids,traits=unique(c('PRODSTG',traits_append)))
  df_add <- df_add |> 
    left_join(df_ID |> dplyr::select(fts_germplasm_id,KEY_ID),
              join_by(fts_germplasm_id == fts_germplasm_id),
              relationship = 'many-to-one',multiple = 'first')
  #df_add |> View()
  
  # Clear object not necessary 
  rm(df_ped,df_pedsel,Inb_hyb_append,ped_all,pedl,
     df_raw01,df_raw02,df_raw03,df_IDmissing,
     df_ID_female,df_ID_male,df_ID_min,df_ID_min2,df_ID_min3)
  
  # Run model 
  cl <- makePSOCKcluster(num_ParRun,
                         outfile = paste0(debug_folder,'\\out_debug.txt'),
                         revtunnel = TRUE)
  registerDoParallel(cl)
  result_list <- foreach::foreach(i=unique(df_raw$TRAIT),
                                  .inorder = FALSE,
                                  .verbose = TRUE,
                                  .packages =c('tidyverse','INLA'),
                                  .export = ls(globalenv())) %dopar% {
                                    tryCatch({
                                      ### Filter data 
                                      #i =  unique(df_raw$TRAIT)[2]
                                      T_i <- Sys.time()
                                      df_model <- df_raw |> 
                                        dplyr::filter(TRAIT == i) |> 
                                        dplyr::select(fts_germplasm_idN,
                                                      male_parent_germplasm_idN,
                                                      female_parent_germplasm_idN,
                                                      season,field_name,block_name,
                                                      observation_trait_val,
                                                      observation_id,
                                                      planting_year,
                                                      plot_id,
                                                      TRAIT) |> 
                                        dplyr::rename(genof = fts_germplasm_idN,
                                                      malef = male_parent_germplasm_idN,
                                                      femalef = female_parent_germplasm_idN) |> 
                                        dplyr::mutate(genof = as.character(genof),
                                                      season_field_name = paste(season,field_name,block_name,sep = '_')) 
                                      ### Model 
                                      Z <- makeZ(data=df_model,female='femalef',male='malef',Cov_m=Cov_m)
                                      df_model$idz <- 1:nrow(df_model)
                                      
                                      ### Standarized data  
                                      y_sd <- sd(df_model$observation_trait_val)
                                      y_mean <- mean(df_model$observation_trait_val)
                                      df_model$y_std <- (df_model$observation_trait_val-y_mean)/y_sd
                                      
                                      ### prior
                                      lgprior <- list(prec = list(prior = "loggamma", param = c(0.01, 0.01),initial = 1))
                                      
                                      ### Model 
                                      fit <- INLA::inla(y_std ~ 1 +
                                                          f(season_field_name,model='iid',constr=T, hyper = lgprior) +
                                                          f(idz,model='z',Z=Z,Cmatrix = Cov_m,constr=T,hyper = lgprior) +
                                                          f(genof,model='iid',constr=T, hyper = lgprior),
                                                        data=df_model,
                                                        verbose=T,
                                                        only.hyperparam = F,
                                                        control.compute=list(dic=F,config=TRUE),
                                                        control.predictor=list(compute = TRUE),
                                                        num.threads = 8)
                                      ##print(summary(fit))
                                      ### Data quality 
                                      df_model$residual <- (df_model$y_std - fit$summary.fitted.values$mean)
                                      df_model$Zscore <- round((df_model$residual-mean(df_model$residual))/sd(df_model$residual),2)
                                      plot_check <- df_model |> 
                                        dplyr::filter(abs(Zscore) > Zscore_filter) |> 
                                        dplyr::select(c('planting_year','season','field_name',
                                                        'block_name','plot_id',
                                                        'TRAIT','observation_id',
                                                        'observation_trait_val','Zscore'))
                                      
                                      ### Extract data model 
                                      int <- fit$summary.fixed[,1:3]
                                      int$mean <- int$mean + y_mean
                                      int$PEV <- int$sd^2 * y_sd
                                      
                                      inb <- fit$summary.random$idz[-(1:nrow(df_model)),1:3]
                                      inb$ID <- as.integer(colnames(Cov_m))
                                      inb$mean <- inb$mean * y_sd
                                      inb$PEV <- (inb$sd)^2 * y_sd
                                      
                                      hyb <- fit$summary.random$genof[,1:3]
                                      hyb$mean <- hyb$mean * y_sd
                                      hyb$PEV <- (hyb$sd)^2 * y_sd
                                      
                                      Va <- extract_trans_to_var(fit,precision="Log precision for idz") * y_sd 
                                      Vna <- extract_trans_to_var(fit,precision="Log precision for genof") * y_sd
                                      Venv <- extract_trans_to_var(fit,precision="Log precision for season_field_name") * y_sd 
                                      Ve <- extract_trans_to_var(fit,precision="Log precision for the Gaussian observations") * y_sd
                                      H2 <- (Va+Vna)/(Va+Vna+Ve) 
                                      h2 <- (Va)/(Va+Vna+Ve)
                                      
                                      ### Inbreds 
                                      df_x <- data.frame(ID = c(df_model$malef,df_model$femalef),
                                                         y = c(df_model$observation_trait_val,
                                                               df_model$observation_trait_val)) |> 
                                        group_by(ID) |> 
                                        summarise(pheno_mean=mean(y),
                                                  n=length(y))
                                      res_inb <- data.frame(trait=i,
                                                            ID=inb$ID,
                                                            u=inb$mean,
                                                            PEV=inb$PEV,
                                                            Intercept=int$mean,
                                                            Intercept_PEV=int$PEV)
                                      
                                      res_inb <- res_inb |> 
                                        left_join(rinf,by=join_by(ID==ID)) |> 
                                        relocate(rinf,.after=ID)
                                      res_inb$EBV_GCA <- rowSums(res_inb[,c('Intercept','u')],na.rm = TRUE)
                                      res_inb$Accuracy <- sqrt(1 - (res_inb$PEV/Va))
                                      res_inb$Reliability <- 1 - (res_inb$PEV/Va)
                                      res_inb$Heritability_ns <- h2
                                      res_inb <- merge(res_inb,df_x,all.x=TRUE,by='ID')
                                      #plot(res_inb$EBV_GCA,res_inb$pheno_mean,asp=1)
                                      
                                      ### Hybrids 
                                      df_x <- df_model |> 
                                        group_by(genof) |> 
                                        summarise(pheno_mean=mean(observation_trait_val),
                                                  n=table(genof))
                                      df_x$genof <- as.integer(df_x$genof)
                                      res_hyb <- data.frame(trait=i,
                                                            ID= hyb$ID,
                                                            u=hyb$mean,
                                                            PEV=hyb$PEV,
                                                            Intercept=int$mean,
                                                            Intercept_PEV=int$PEV)
                                      df_data <- unique(df_model[,c('femalef','malef','genof')])
                                      df_data$genof <- as.integer(df_data$genof)
                                      res_hyb <- merge(df_data,res_hyb,by.x='genof','ID',all.x = TRUE)
                                      res_hyb$femalef <- as.integer(res_hyb$femalef)
                                      res_inb$ID <- as.integer(res_inb$ID)
                                      res_hyb <- merge(res_hyb,res_inb[,c('ID','u','PEV')],
                                                       by.x='femalef',
                                                       by.y= 'ID',
                                                       all.x = TRUE,
                                                       suffixes = c("_h","_f"))
                                      res_hyb$malef <- as.integer(res_hyb$malef)
                                      res_inb$ID <- as.integer(res_inb$ID)
                                      res_hyb <- merge(res_hyb,res_inb[,c('ID','u','PEV')],
                                                       by.x='malef',
                                                       by.y='ID',
                                                       all.x = TRUE,
                                                       no.dups = TRUE) |> dplyr::rename(u_m=u,
                                                                                        PEV_m=PEV)
                                      res_hyb$EBV_GCA_SCA <- rowSums(res_hyb[,c('Intercept','u_h','u_m','u_f')],na.rm = TRUE)
                                      res_hyb$PEV_pre <- rowSums(res_hyb[,c('PEV_h','PEV_m','PEV_f')],na.rm = TRUE)
                                      res_hyb$Accuracy <- sqrt(1 - (res_hyb$PEV_pre/(2*Va+Vna)))
                                      res_hyb$Reliability <- 1 - (res_hyb$PEV_pre/(2*Va+Vna))
                                      res_hyb$Heritability_bs <- H2
                                      res_hyb$Va <- Va
                                      res_hyb$Vna <- Vna
                                      res_hyb$Ve <- Ve
                                      res_hyb$Venv <- Venv
                                      res_hyb <- merge(res_hyb,df_x,all.x=TRUE,by.x='genof',by.y='genof') 
                                      res_hyb <- res_hyb |> relocate(trait)
                                      #plot(res_hyb$EBV_GCA_SCA,res_hyb$pheno_mean,asp=1)
                                      Time <- difftime(Sys.time(), T_i,units = "secs")
                                      cat(paste('at:',Sys.time(),'## Done: ',i,'in :',round(Time,2),'secs\n'))
                                      out <- list(i = list(res_hyb=res_hyb,res_inb=res_inb,fit=fit,Time=Time,plot_check=plot_check))
                                      names(out) <- i
                                      return(out)
                                    },
                                    error =  function(e,...){
                                      out <- list(list(res_hyb=NA,res_inb=NA,fit=NA,Time=NA,plot_check=NA))
                                      names(out) <- i
                                      return(out)
                                    }
                                    )
                                  }
  parallel::stopCluster(cl)
  res_out_hybt <- res_out_inbt <- res_out_fit <- res_out_time <- res_out_plots_check <- list()
  con <- 1
  model_name <- c()
  for(i in 1:length(result_list)){
    res_out_hybt[[con]] <- result_list[[i]][[1]]$res_hyb  
    res_out_inbt[[con]] <- result_list[[i]][[1]]$res_inb
    res_out_fit[[con]] <- result_list[[i]][[1]]$fit
    res_out_time[[con]] <- result_list[[i]][[1]]$Time
    res_out_plots_check[[con]] <- result_list[[i]][[1]]$plot_check
    model_name <- c(model_name, names(result_list[[i]]))
    con <- con+1
  }
  
  res_out_hyb <- do.call(rbind,res_out_hybt)
  res_out_hyb <- res_out_hyb |> dplyr::filter(!is.na(trait))
  
  res_out_inb <- do.call(rbind,res_out_inbt)
  res_out_inb <- res_out_inb |> dplyr::filter(!is.na(trait))
  
  ## Data quality 
  res_out_plots_check <- do.call(rbind,res_out_plots_check)
  res_out_plots_check <-  res_out_plots_check |> dplyr::filter(!is.na(Zscore))
  
  if(nrow(res_out_plots_check)>1){
    ggplot(res_out_plots_check, aes(x=TRAIT, y=Zscore, fill=TRAIT)) +
      geom_point(size=2.5)+
      guides(fill='none')+
      ylab('Zscore (Standarized residuals)')+
      scale_x_discrete(guide = guide_axis(angle = 90)) +
      geom_hline(yintercept = c(-Zscore_filter,Zscore_filter),lty=2,color='red',lwd=1)+
      theme_bw()
    ggsave(paste0(debug_folder,'\\out_Zscore.png'),w=20,h=10)
  }
  
  # Export time consumed 
  res_out_Timep <- data.frame(Trait=model_name)
  for(i in 1:length(res_out_time)){
    x <- abs(as.numeric(res_out_time[i]))
    res_out_Timep$time_seg[i] <- x
    res_out_Timep$time_HHMMSS[i] <-  formatSeg(x)
  }
  
  # Model fitted plot diagnostic 
  names(res_out_fit) <- model_name
  pdf(paste0(debug_folder,'\\out_Diagplot.pdf'),h=7,w=12)
  for (i in 1:length(res_out_fit)){
    if(!is.na(res_out_fit[i])){
      plot(0,0,axes = F,main= model_name[i])
      try(plot(res_out_fit[[i]],
               single = T,
               plot.prior = T))
    }
  }
  dev.off()
  
  # Merge original id 
  res_out_inb <- res_out_inb |> 
    left_join(pedln$xrf,join_by(ID == newID), multiple = 'first', relationship = 'many-to-one')|> 
    dplyr::rename(fts_germplasm_id = ID.y)
  
  res_out_hyb <- res_out_hyb |> 
    left_join(pedln$xrf,join_by(genof == newID), multiple = 'first', relationship = 'many-to-one')|> 
    dplyr::rename(fts_germplasm_idH = ID) |> 
    left_join(pedln$xrf,join_by(femalef == newID), multiple = 'first', relationship = 'many-to-one')|> 
    dplyr::rename(fts_germplasm_idF = ID) |> 
    left_join(pedln$xrf,join_by(malef == newID), multiple = 'first', relationship = 'many-to-one')|> 
    dplyr::rename(fts_germplasm_idM = ID)
  
  # Merge results and append data 
  res_out_hyb <- res_out_hyb |> 
    left_join(df_add,join_by(fts_germplasm_idH == fts_germplasm_id),
              multiple = 'first', relationship = 'many-to-one') |> 
    relocate(trait,genof,pedigree,origin,germplasm_lineage_identifier)
  
  ## For fix checks with no df_add
  for(i in (unique(res_out_hyb[is.na(res_out_hyb$pedigree),'genof']))){
    df_filter <- res_out_hyb$genof %in% i
    res_out_hyb[df_filter,c('pedigree','origin','germplasm_lineage_identifier','cross_name','KEY_ID')] <- 
      df_raw[df_raw$MIN_fts_germplasm_id %in% res_out_hyb[df_filter,'fts_germplasm_idH'],
             c('pedigree','origin','germplasm_lineage_identifier','cross_name','KEY_ID')][1,]
  }
  res_out_hyb <- res_out_hyb  |> dplyr::select(-any_of(c('malef','femalef','genof')))
  
  res_out_hyb$PRODSTG <- as.character(res_out_hyb$PRODSTG)
  res_out_hyb[is.na(res_out_hyb$PRODSTG),'PRODSTG'] <- 'CHECK'
  
  res_out_inb <-  res_out_inb |> 
    left_join(df_add,join_by(fts_germplasm_id == fts_germplasm_id),
              multiple = 'first', relationship = 'many-to-one') |>
    dplyr::filter(!is.na(pedigree)) |> 
    relocate(trait,ID,pedigree,origin,germplasm_lineage_identifier) |> 
    dplyr::select(-any_of('ID'))
  
  if(Exclude_Hybrid_GCA){
    res_out_inb <-  res_out_inb |> dplyr::filter(line_type != 'Hybrid')
  }
  # Merge Summary 
  res_her <- merge(res_out_inb |> dplyr::select(trait,Heritability_ns) |> unique(),
                   res_out_hyb |> dplyr::select(trait,Heritability_bs,Va,Vna,Ve,Venv) |> unique()) |> 
    dplyr::mutate(p_Va=Va/(Va+Vna),p_Vna=Vna/(Va+Vna))
  
  df_rawS <- merge(df_rawS,res_her,by.x='TRAIT',by.y='trait',all.x=TRUE)
  
  ## Check if some KEYID is missing from hybrid results
  KEYID_backup <- KEYID_backup[!(KEYID_backup %in% res_out_hyb$KEY_ID)]
  
  # Plot Prediction 
  res_out_inb |> ggplot(aes(x=pheno_mean,y=EBV_GCA))+
    geom_point()+
    facet_wrap(~trait, scales='free')+
    geom_abline(intercept = 0, slope =1)+
    theme_bw()
  ggsave(paste0(debug_folder,'\\out_pred_inbred.png'),w=15,h=15)
  
  res_out_hyb |> ggplot(aes(x=pheno_mean,y=EBV_GCA_SCA))+
    geom_point()+
    facet_wrap(~trait, scales='free')+
    geom_abline(intercept = 0, slope =1)+
    theme_bw()
  ggsave(paste0(debug_folder,'\\out_pred_hybrid.png'),w=15,h=15)
  
  # Save results 
  wb <- createWorkbook()
  addWorksheet(wb, "SummaryRawData_Pred")
  addWorksheet(wb, "Results_Hybrid")
  addWorksheet(wb, "Results_Inbred")
  addWorksheet(wb, "Time_run_model")
  addWorksheet(wb, "Raw_data")
  addWorksheet(wb, "Missing_KEYID")
  addWorksheet(wb, "Zscore")
  writeData(wb, sheet = "SummaryRawData_Pred", x = df_rawS, withFilter = TRUE, na.string='')
  writeData(wb, sheet = "Results_Hybrid", x = res_out_hyb, withFilter = TRUE, na.string='')
  writeData(wb, sheet = "Results_Inbred", x = res_out_inb, withFilter = TRUE, na.string='')
  writeData(wb, sheet = "Time_run_model", x = res_out_Timep, withFilter = TRUE, na.string='')
  writeData(wb, sheet = "Raw_data", x = df_raw, withFilter = TRUE, na.string='')
  writeData(wb, sheet = "Missing_KEYID", x = KEYID_backup, withFilter = TRUE, na.string='')
  writeData(wb, sheet = "Zscore", x =  res_out_plots_check , withFilter = TRUE, na.string='')
  freezePane(wb, sheet = "SummaryRawData_Pred", firstRow = TRUE, firstCol = FALSE) 
  freezePane(wb, sheet = "Results_Hybrid", firstRow = TRUE, firstCol = FALSE) 
  freezePane(wb, sheet = "Results_Inbred", firstRow = TRUE, firstCol = FALSE) 
  freezePane(wb, sheet = "Time_run_model", firstRow = TRUE, firstCol = FALSE) 
  freezePane(wb, sheet = "Raw_data", firstRow = TRUE, firstCol = FALSE) 
  freezePane(wb, sheet = "Missing_KEYID", firstRow = TRUE, firstCol = FALSE) 
  freezePane(wb, sheet = "Zscore", firstRow = TRUE, firstCol = FALSE) 
  saveWorkbook(wb, file = paste0(file_out,".xlsx"), overwrite = TRUE)
  return(list(df_rawS,res_out_Timep))
}

#' Adjust prod stage related with velocity issue
#' 
#' @title 
#' @param 
#' @details 
#' @return 
#' 
adj_prodstg <- function(out_res_file='out_YT_Results_PAR.xlsx',
                        in_prod_file='input_files/in_PRODSTG.xlsx'){
  df_hyb <- read_excel(out_res_file,'Results_Hybrid')
  df_prodstg <- read_excel(in_prod_file)
  
  df_checks <- df_hyb |> 
    filter(PRODSTG == 'CHECK') |> 
    select(PRODSTG,KEY_ID,pedigree) |> 
    unique()
  
  df_prodstg <- rbind(df_prodstg,df_checks) |> 
    select(-pedigree)
  
  df_hyb <- df_hyb |> 
    left_join(df_prodstg,join_by(KEY_ID == KEY_ID)) |> 
    select(-PRODSTG.x) |> 
    mutate(PRODSTG.y = replace_na(PRODSTG.y,'DROP')) |> 
    rename(PRODSTG = PRODSTG.y)
  
  wb <- loadWorkbook(out_res_file)
  writeData(wb, sheet = "Results_Hybrid", x = df_hyb, withFilter = TRUE, na.string='')
  saveWorkbook(wb, file = out_res_file, overwrite = TRUE)
}


#' Adjust prod stage related with velocity issue
#' 
#' @title 
#' @param 
#' @details 
#' @return 
#' 
runHM <- function(
    data = read_excel('out_YT_Results_PAR.xlsx','Results_Inbred'),
    inbred = read_excel('input_files/in_Inbreds_to_HM.xlsx'),
    traits = read_excel('input_files/in_TIX_traits.xlsx'),
    CheckInventory = T,
    owner_breeding_prog = 'AV',
    out_res_folder = 'post_HM',
    GCT_append = c('P50', 'S50', 'HETGP','ENC')){
  
  traits <- traits |> filter(TIX_name == 'All')
  #TwoWayRep <- FALSE
  
  if(sum(traits$Weight) != 1){stop('The TIX in class All in input file do not sum 1, please check!')}
  
  female_KEY_ID = unique(na.omit(inbred$female_KEY_ID))
  male_KEY_ID = unique(na.omit(inbred$male_KEY_ID))
  length(female_KEY_ID)*length(male_KEY_ID)
  
  all_parents <- c(female_KEY_ID,male_KEY_ID) 
  all_parents <- all_parents[!(all_parents %in% data$KEY_ID)]
  if(length(all_parents)>1){
    write.csv(all_parents,paste0(out_res_folder,'/out_Parents_notFound_inGCA.csv'),row.names = FALSE)
  }
  # Hybrid Mining 
  # data and HM 
  dataP <- data |>
    filter(KEY_ID %in% c(female_KEY_ID,male_KEY_ID),
           trait %in% traits$trait) |>
    pivot_wider(id_cols = KEY_ID,
                names_from = trait,
                values_from = EBV_GCA,
                values_fn = mean)
  
  hybrid_pred <- expand.grid(female_KEY_ID=as.character(female_KEY_ID),
                             male_KEY_ID=as.character(male_KEY_ID)) |>
    filter(!(as.character(female_KEY_ID) == as.character(male_KEY_ID)))
  
  hybrid_pred$HBCPED <- apply(hybrid_pred,1,function(x){
    paste0(x[order(c(as.character(x)))],collapse = '+')})
  
  # if(TwoWayRep){
  #   hybrid_pred$HBCPED <- apply(hybrid_pred,1,function(x){paste0(sort(x),collapse = '+')})
  # }
  
  hybrid_pred <- hybrid_pred[!duplicated(hybrid_pred$HBCPED),]
  
  hybrid_pred <- hybrid_pred |>
    left_join(dataP, by = join_by(female_KEY_ID == KEY_ID)) |>
    left_join(dataP, by = join_by(male_KEY_ID == KEY_ID),
              suffix = c(".f", ".m"))
  
  ### Function to standardization 
  for(i in traits$trait){
    iname <- paste0(i,'_pred')
    inamestan <- paste0(i,'_predS')
    hybrid_pred <- hybrid_pred |>
      mutate({{iname}} := rowMeans(
        hybrid_pred[,grep(paste(paste('^',i,c('.f','.m'),sep=''),
                                collapse = '|'),colnames(hybrid_pred))]))
    hybrid_pred <- hybrid_pred |>
      mutate({{inamestan}} := st_mean(x=hybrid_pred[,iname],
                                      sel_dir=traits$trait_direction[traits$trait==i ],
                                      min_range=traits$Minimal[traits$trait==i ],
                                      max_range=traits$Maximum[traits$trait==i ]) *
               traits$Weight[traits$trait==i ])
  }
  
  hybrid_pred$TIX <- rowSums(hybrid_pred[,grep('_predS$',colnames(hybrid_pred))],na.rm = TRUE)
  
  ## Append Inbred data 
  dataAp <- data |>
    filter(KEY_ID %in% c(female_KEY_ID,male_KEY_ID)) |>
    select(all_of(c('KEY_ID','origin',GCT_append))) |>
    unique()
  
  hybrid_pred <- hybrid_pred |>
    left_join(dataAp, by = join_by(female_KEY_ID == KEY_ID)) |>
    left_join(dataAp, by = join_by(male_KEY_ID == KEY_ID),
              suffix = c(".f", ".m")) |>
    mutate(HETGP=paste0(HETGP.f,"/",HETGP.m),
           ASI50D=as.numeric(S50.f)-as.numeric(P50.m)) 
  
  hybrid_pred$HETGP_Alpha <- apply(hybrid_pred[,c('HETGP.f','HETGP.m')],
                                   1,
                                   function(x){ifelse(length(na.omit(x))==2,paste0(sort(x),collapse = '/'),NA)})
  
  ## Check is hybrid is already made 
  dataP <- data |> select(pedigree,KEY_ID) |> unique()
  
  hybrid_pred <- hybrid_pred |> 
    left_join(dataP,by=join_by(female_KEY_ID == KEY_ID)) |> 
    rename(pedigreeF=pedigree) |> 
    left_join(dataP,by=join_by(male_KEY_ID == KEY_ID)) |> 
    rename(pedigreeM=pedigree) |>
    mutate(pedigree=paste(pedigreeF,pedigreeM,sep='+'),
           pedigree_inv=paste(pedigreeM,pedigreeF,sep='+'))
  
  if(CheckInventory){
    find_cross <- source_inv_byHBCPED(owner_breeding_prog)
    hybrid_pred <- hybrid_pred |> 
      left_join(find_cross,by=join_by(HBCPED == HBCPED)) |> 
      dplyr::rename(Ninv_HBCPED= Ninv ) |> 
      mutate(Ninv_HBCPED = replace_na(Ninv_HBCPED,0),
             Status = if_else(Ninv_HBCPED > 0,'Already_make','No_maked'),
             TIX = if_else(TIX==0,NA,TIX))
  }else{
    hybrid_pred <- hybrid_pred |> 
      mutate(Ninv_HBCPED = 'NotChecked',
             Status = 'NotChecked',
             TIX = if_else(TIX==0,NA,TIX))
    
  }
  hybrid_pred <- hybrid_pred |>
    relocate(c(TIX,origin.f,origin.m,HETGP,HETGP_Alpha,ASI50D,Status,Ninv_HBCPED),.after = HBCPED)
  
  ### Save results 
  wb <- createWorkbook()
  addWorksheet(wb, "SummaryPrediction")
  writeData(wb, 1, x = hybrid_pred, withFilter = TRUE)
  freezePane(wb, "SummaryPrediction", firstActiveRow = 2, firstActiveCol = 5)
  saveWorkbook(wb, file = paste0(out_res_folder,"/out_HM_Prediction.xlsx"), overwrite = TRUE)
  
  # Inbred GCA 
  dataI <- data |>
    filter(trait %in% traits$trait) |>
    pivot_wider(id_cols = c(KEY_ID,pedigree),
                names_from = trait,
                values_from = EBV_GCA,
                values_fn = mean) |> 
    as.data.frame()
  
  ### Function to standardization 
  for(i in traits$trait){
    inamestan <- paste0(i,'_predS')
    dataI <- dataI |>
      mutate({{inamestan}} := st_mean(x=dataI[,i],
                                      sel_dir=traits$trait_direction[traits$trait==i ],
                                      min_range=traits$Minimal[traits$trait==i ],
                                      max_range=traits$Maximum[traits$trait==i ]) *
               traits$Weight[traits$trait==i ])
  }
  
  dataI$TIX <- rowSums(dataI[,grep('_predS$',colnames(dataI))])
  
  ## Append Inbred data 
  dataI <- dataI |>
    left_join(dataAp, by = join_by(KEY_ID == KEY_ID)) |>
    relocate(c(TIX),.after = pedigree)
  
  ### Save results -
  wb <- createWorkbook()
  addWorksheet(wb, "GCA_Inbred")
  writeData(wb, sheet = "GCA_Inbred", x = dataI, withFilter = TRUE)
  freezePane(wb, "GCA_Inbred", firstActiveRow = 2, firstActiveCol = 4)
  saveWorkbook(wb, file = paste0(out_res_folder,"/out_GCA_Inbred.xlsx"), overwrite = TRUE)
}


#' Adjust prod stage related with velocity issue
#' 
#' @title 
#' @param 
#' @details 
#' @return 
#' 
makePLOTS <- function(
    ## Load graphics parameter
  in_parp = "input_files/in_Graphics_Parameters.xlsx",
  dataPp = "post_Make_TIX/out_PivotData_BLUP_TIX_run.xlsx",
  ##head_file = readLines('Adv_MG_head.Rmd'),
  out_folder = "post_ADV_plots",
  selrows = 1:3){
  
  if(file.exists(file.path(getwd(),paste0(out_folder,'/out_plot')))){
    unlink(file.path(getwd(),paste0(out_folder,'/out_plot')),recursive=TRUE)
    dir.create(file.path(getwd(),paste0(out_folder,'/out_plot')))
  }else{
    dir.create(file.path(getwd(),paste0(out_folder,'/out_plot')))
  }
  dataP <-  readxl::read_excel(dataPp)
  if(!is.null(selrows)){
    in_par <- readxl::read_excel(in_parp)[selrows,]
  }else{
    in_par <- readxl::read_excel(in_parp)
  }
  
  head_file <-  paste0('
---
  title: "Adv Meeting graphics"
  output: 
    html_document:
    toc: TRUE
  toc_float: TRUE
---
  
```{css toc-content, echo = FALSE}
  #TOC {
  right: 50px;
  margin: 20px 0px 25px 0px;
  }
  .main-container {
  margin-left: 50px;
  }
```
```{r include=FALSE}
  pak <- c("openxlsx","readxl","tidyverse","ggrepel","DT","ggplot2")
  pakl <- unlist(lapply(pak, require, character.only = TRUE,quietly=TRUE))
  table(pakl)
  pak[!pakl]
  source("../fun/aux_functions.R")
```
  
  # Load data 
```{r echo=TRUE}
  dataP <- read_excel("../',dataPp,'")
```
  ')
  ## Check data 
  columnsSEL <- unique(c(in_par$x,in_par$y,in_par$z,in_par$id_label,in_par$Stage))
  
  if(sum(!(columnsSEL) %in% colnames(dataP)) > 0){
    stop(paste('Columns of data in Graphics parameter files is no include in BLUP Data:',
               paste(columnsSEL[!(columnsSEL %in% colnames(dataP))] ,collapse = ',')))
  }
  if(sum(duplicated(apply(in_par[,c('prefix','x','y','z')],1,paste0,collapse=''))) > 0){
    cat('WARNING: The combination of prefix, x, y, z, is not unique, the graphics export 
       will not cover all graphics requested, check Graphics parameters file')
  }
  
  ## Create string with graphics
  session <- in_par$session[1] 
  tmpG <- c("\n")
  
  for(i in 1:nrow(in_par)){
    if(i == 1){
      tmpG <- paste(tmpG,"# ",in_par$session[i],"\n",sep='')
    }
    if(i > 1 & session != in_par$session[i] ){
      tmpG <- paste(tmpG,"# ",in_par$session[i],"\n",sep='')
      session <- in_par$session[i] 
    }
    tmpG <- paste(tmpG,
                  "## ", in_par$x[i], " x ", in_par$y[i], " x ", in_par$z[i], "\n" ,
                  "\n",
                  "```{r echo=FALSE, fig.height=7, fig.width=12, warning=FALSE}\n",
                  "df_sel <- dataP |> dplyr::filter(",in_par$data_filter[i],")\n",
                  "adv_plot(data = ",in_par$data[i],",\n",
                  "x= '",in_par$x[i],"',\n",
                  "y= '",in_par$y[i],"',\n",
                  "z= '",in_par$z[i],"',\n",
                  "x_direct= '",in_par$x_direct[i],"',\n",
                  "y_direct= '",in_par$y_direct[i],"',\n",
                  "XLABS= '",in_par$XLABS[i],"',\n",
                  "YLABS= '",in_par$YLABS[i],"',\n",
                  "ZLABS= '",in_par$ZLABS[i],"',\n",
                  "x_per_add= ",in_par$x_per_add[i],",\n",
                  "y_per_add= ",in_par$y_per_add[i],",\n",
                  "prefix= '",in_par$prefix[i],"',\n",
                  "Nclass= ",in_par$Nclass[i],",\n",
                  "id_label= '",in_par$id_label[i],"',\n",
                  "save_plot= '",in_par$save_plot[i],"',\n",
                  "Stage= '",in_par$Stage[i],"')\n",
                  "```\n",
                  "\n",sep='')
  }
  #tmpG
  writeLines(paste(paste(head_file,collapse = '\n'), tmpG),paste0(out_folder,'/ADV_Meeting_graphics.Rmd'))
  rmarkdown::render(paste0(out_folder,"/ADV_Meeting_graphics.Rmd"))
}



# INBRED------------------------------------------------------------------------

#' Source rawdata from CSW to Hybrid Analysis Pipeline
#' 
#' @title 
#' @param WHERE filter to apply in WHERE SQL query from CSW in table bcs-veg-mart.performance.veg_experiment_obs_entries
#' @details the columns mandatory for query is already fixed
#' @return data frame with raw data with no treatment from CSW
#' 
source_raw_data_INB <- function(WHERE=NULL){
  project <- "bcs-veg-rd"
  sql <- paste("SELECT
  distinct
  veoe.trial_type, veoe.season, veoe.country, veoe.field_name, veoe.block_name,
  veoe.experiment_stage_cd, CAST(veoe.plot_id AS INT64) AS plot_id, CAST(veoe.fts_germplasm_id AS INT64) AS fts_germplasm_id,
  veoe.pedigree, veoe.origin, veoe.germplasm_lineage_identifier, veoe.cross_name,veoe.observation_id,
  prod.primary_commercial_product_nm as PROD, veoe.line_type, veoe.observation_reference_cd,
  veoe.descriptor_abbreviation_cd, veoe.observation_trait_val,
  
  COALESCE(primary_question_cd,observation_reference_cd) AS TRAIT,
  COALESCE(secondary_answers_concat,descriptor_abbreviation_cd) AS TRAIT_ABBREVIATION,
  primary_question_cd,
  primary_question_cd_desc,
  secondary_answers_concat
FROM
  `bcs-veg-mart.performance.veg_experiment_obs_entries` veoe
  left join (`bcs-veg-mart.pipeline.veg_germplasms` vg CROSS JOIN UNNEST (products) as prod) on veoe.fts_germplasm_id = vg.fts_germplasm_id
WHERE
  veoe.observation_trait_val IS NOT NULL
  AND planting_date IS NOT NULL
  AND veoe.is_plot_deactivated = false
  AND ",WHERE,";")
  tb <- bq_project_query(project, sql)
  df_raw <- bq_table_download(tb,bigint="integer64",page_size = 20000) |> as.data.frame()
  print(paste('DataFrame dimension: ',dim(df_raw)[1],' rows'))
  return(df_raw)
}



#' Source rawdata from CSW to Hybrid Analysis Pipeline
#' 
#' @title 
#' @param WHERE filter to apply in WHERE SQL query from CSW in table bcs-veg-mart.performance.veg_experiment_obs_entries
#' @details the columns mandatory for query is already fixed
#' @return data frame with raw data with no treatment from CSW
#' 
run_modelINB <- function(
    df_raw = NULL,
    minGen = NULL,
    Use10genPed = TRUE,
    file_out = NULL,
    debug_folder = 'out_DEBUG',
    traits_append = NULL,
    num_ParRun = 1){
  
  # Adjust right functions 
  select <- dplyr::select
  filter <- dplyr::filter
  foreach <- foreach::foreach
  
  # Adjust folder to debug
  if(file.exists(file.path(getwd(),debug_folder))){
    unlink(file.path(getwd(),debug_folder),recursive=TRUE)
    dir.create(file.path(getwd(),debug_folder))
  }else{dir.create(file.path(getwd(),debug_folder))}
  
  # Adjust data issues expand tratis with sublevel 
  df_raw <- df_raw |> dplyr::filter(TRAIT != 'CMNTS')
  df_raw <- check_data_adj(data=df_raw)
  df_raw <- makeKEY_ID(data=df_raw)
  KEYID_backup <- unique(df_raw$KEY_ID)
  
  ## Summary data 
  df_rawS <- df_raw |> 
    group_by(TRAIT) |> 
    summarise(min=min(observation_trait_val,na.rm=TRUE),
              Q05=quantile(observation_trait_val,0.05),
              Q25=quantile(observation_trait_val,0.25),
              mean=mean(observation_trait_val,na.rm=TRUE),
              Q75=quantile(observation_trait_val,0.75),
              Q95=quantile(observation_trait_val,0.95),
              max=max(observation_trait_val,na.rm=TRUE),
              n_genotypes=dplyr::n_distinct(fts_germplasm_id),
              n=dplyr::n()) |> 
    dplyr::mutate(across(min:max,function(x) round(x,2))) |> 
    as.data.frame()
  print(df_rawS)
  
  ## Remove traits with low genotypes
  df_raw <- df_raw |> 
    dplyr::filter(TRAIT %in%
                    df_rawS$TRAIT[df_rawS$n_genotypes > minGen]) |> 
    droplevels()
  
  ## data quality plot 
  df_raw |> 
    ggplot(aes(x=observation_trait_val,fill=TRAIT))+
    geom_density()+
    guides(fill = FALSE)+
    facet_wrap(~TRAIT, scales='free')+
    theme_bw() 
  ggsave(paste0(debug_folder,'\\out_DQ.png'),w=15,h=15)
  
  ped_all <- df_raw[,c('fts_germplasm_id','line_type', 'KEY_ID')] |> unique()
  
  # Loading Pedigree 
  ## Source deep pedigree 5 generation or more 
  df_ped <- source_deep_ped2(fts_germplasm_id=unique(ped_all$fts_germplasm_id))
  if(nrow(df_ped) == 0){
    stop('No deep pedigre are returned in this analysis, check you selection')
  }
  if(Use10genPed){
    df_ped <- source_deep_ped2(
      unique(na.omit(c(df_ped$child_germplasm_id,
                       df_ped$male_parent_germplasm_id,
                       df_ped$female_parent_germplasm_id)))) |> 
      arrange(child_germplasm_id)
    print(dim(df_ped))
  }
  
  df_ped <- df_ped |> 
    unique() |> 
    arrange(child_germplasm_id,male_parent_germplasm_id,female_parent_germplasm_id)
  
  #df_ped |> View()
  ## Master ID file
  df_ID <- source_all_ID_ByOrigin_ID(fts_germplasm_id = 
                                       unique(c(df_ped$child_germplasm_id,
                                                df_ped$male_parent_germplasm_id,
                                                df_ped$female_parent_germplasm_id,
                                                df_raw$fts_germplasm_id))) 
  
  df_IDmissing <- df_raw[!(df_raw$fts_germplasm_id %in% df_ID$fts_germplasm_id),] |> 
    dplyr::select(c('fts_germplasm_id','line_type','pedigree','origin','cross_name','germplasm_lineage_identifier','PROD')) |> 
    unique() |> 
    dplyr::mutate(male_parent_germplasm_id = NA,
                  female_parent_germplasm_id = NA) |> 
    droplevels()
  
  df_ID <- rbind(df_ID,df_IDmissing) |> 
    makeKEY_ID() |> 
    dplyr::filter(!is.na(KEY_ID)) |> 
    group_by(KEY_ID) |> 
    dplyr::mutate(MIN_fts_germplasm_id = min(fts_germplasm_id),
                  rankID = dense_rank(fts_germplasm_id )) |> 
    ungroup() |> 
    arrange(KEY_ID,fts_germplasm_id) |> 
    as.data.frame()
  
  df_ID_male <- df_ID |> dplyr::select(fts_germplasm_id,MIN_fts_germplasm_id) |> dplyr::rename(MIN_male_parent_germplasm_id = MIN_fts_germplasm_id)
  df_ID_female <- df_ID |> dplyr::select(fts_germplasm_id,MIN_fts_germplasm_id) |> dplyr::rename(MIN_female_parent_germplasm_id = MIN_fts_germplasm_id)
  
  df_ID <-  df_ID |> 
    left_join(df_ID_male,by=join_by(male_parent_germplasm_id == fts_germplasm_id),
              multiple = 'first',relationship = 'many-to-one') |> 
    left_join(df_ID_female, by=join_by(female_parent_germplasm_id == fts_germplasm_id),
              multiple = 'first',relationship = 'many-to-one') |> 
    dplyr::mutate(MIN_male_parent_germplasm_id =
                    case_when(is.na(MIN_male_parent_germplasm_id) ~
                                male_parent_germplasm_id,.default = MIN_male_parent_germplasm_id),
                  MIN_female_parent_germplasm_id =
                    case_when(is.na(MIN_female_parent_germplasm_id) ~
                                female_parent_germplasm_id,.default = MIN_female_parent_germplasm_id)) |>  
    dplyr::select(-c(male_parent_germplasm_id, female_parent_germplasm_id))
  
  #df_ID |> View()
  
  df_ID_min <- df_ID |> dplyr::select(fts_germplasm_id,MIN_fts_germplasm_id)
  df_ID_min2 <- df_ID |> dplyr::filter(rankID == 1) |> dplyr::select(MIN_fts_germplasm_id,MIN_male_parent_germplasm_id, MIN_female_parent_germplasm_id) |> unique()
  
  df_ped <- df_ped |> 
    dplyr::select(child_germplasm_id) |> 
    left_join(df_ID_min,by=join_by(child_germplasm_id == fts_germplasm_id),
              multiple = 'first',relationship = 'many-to-one') |> 
    left_join(df_ID_min2,by=join_by(MIN_fts_germplasm_id == MIN_fts_germplasm_id),
              multiple = 'first',relationship = 'many-to-one') #|>  dplyr::rename(MIN_male_parent_germplasm_id = male_parent_germplasm_id,         MIN_female_parent_germplasm_id = female_parent_germplasm_id)
  
  ## needs improvement!!!!!!!!!!! 
  df_ped[is.na(df_ped)] <- 0
  df_ped <- df_ped |> dplyr::filter(!(MIN_fts_germplasm_id == MIN_female_parent_germplasm_id))
  df_ped <- df_ped |> dplyr::filter(!(MIN_fts_germplasm_id == MIN_male_parent_germplasm_id))
  
  ## Replace with MIN FTS in df_raw 
  #df_ped|> View()
  df_raw01 <- df_raw |> dplyr::filter(fts_germplasm_id %in% df_ID$MIN_fts_germplasm_id) |> droplevels() #fts_ID already is minimal
  df_raw02 <- df_raw |> dplyr::filter(!(fts_germplasm_id %in% df_ID$fts_germplasm_id)) |> droplevels()  #fts_ID no include in df_ID
  df_raw03 <- df_raw |> dplyr::filter(fts_germplasm_id %in% df_ID$fts_germplasm_id & 
                                        !(df_raw$fts_germplasm_id %in% df_ID$MIN_fts_germplasm_id)) |> droplevels() #fts_ID ins't minimal
  
  df_raw03 <- df_raw03 |> dplyr::select(-c(fts_germplasm_id,
                                           pedigree,
                                           origin,
                                           germplasm_lineage_identifier,
                                           cross_name,
                                           PROD))
  df_ID_min3 <- df_ID |> dplyr::filter(rankID == 1) |> dplyr::select(-c(MIN_male_parent_germplasm_id, MIN_female_parent_germplasm_id))
  
  df_raw03 <- df_raw03 |> 
    dplyr::select(-c(line_type)) |> 
    left_join(df_ID_min3, by=join_by(KEY_ID == KEY_ID ),
              multiple = 'first',relationship = 'many-to-one') |>
    dplyr::select(-c(fts_germplasm_id,rankID)) |> 
    dplyr::rename(fts_germplasm_id= MIN_fts_germplasm_id)
  df_raw03 <- df_raw03[,match(colnames(df_raw03),colnames(df_raw01))]
  
  df_raw <- rbind(df_raw01,df_raw02,df_raw03) |>
    dplyr::filter(!is.na(fts_germplasm_id)) |> 
    droplevels() |> 
    dplyr::rename(MIN_fts_germplasm_id = fts_germplasm_id)
  
  #df_raw |> View()
  ## joint parents ID with data 
  df_raw <- df_raw |> 
    left_join(unique(df_ped[,c('MIN_fts_germplasm_id',
                               'MIN_male_parent_germplasm_id',
                               'MIN_female_parent_germplasm_id')]),
              by= join_by(MIN_fts_germplasm_id == MIN_fts_germplasm_id),
              multiple = 'first',relationship = 'many-to-one')
  
  #df_raw |> View()
  df_pedsel <- df_ped |> 
    dplyr::select(MIN_fts_germplasm_id,
                  MIN_male_parent_germplasm_id,
                  MIN_female_parent_germplasm_id) |> 
    unique() 
  
  ## find no parentage (checks) 
  ### imputed fake parents
  if(anyNA(c(df_raw$MIN_male_parent_germplasm_id,df_raw$MIN_female_parent_germplasm_id))){
    Noparentage <- imp_fake_parentage2(df_ped=df_ped,data=df_raw)
    for(i in 1:nrow(Noparentage)){
      sel <- Noparentage$MIN_fts_germplasm_id[i]
      sel1 <- df_raw$MIN_fts_germplasm_id == sel
      sel2 <- Noparentage$MIN_fts_germplasm_id == sel
      df_raw[sel1, 'MIN_male_parent_germplasm_id'] <- Noparentage[sel2 ,'MIN_male_parent_germplasm_id']
      df_raw[sel1, 'MIN_female_parent_germplasm_id'] <- Noparentage[sel2,'MIN_female_parent_germplasm_id']
    }
    df_pedsel <- rbind(Noparentage,df_pedsel)
  }
  
  # make pedigree data 
  pedl <- insertPedM(df_pedsel)
  pedl <- pedl[!(pedl$MIN_fts_germplasm_id == pedl$MIN_male_parent_germplasm_id),]
  pedl[is.na(pedl)] <- 0
  pedln <- renum(pedl)
  pedln$newped[pedln$newped == 0] <- NA
  dim(pedln$newped)
  
  ## Merge id numeric 
  pedln$xrf$ID <- as.integer64(pedln$xrf$ID)
  df_raw <- df_raw |> 
    left_join(pedln$xrf,by=join_by(MIN_fts_germplasm_id == ID),relationship = 'many-to-one') |>
    dplyr::rename(fts_germplasm_idN=newID) |>
    left_join(pedln$xrf,by=join_by(MIN_female_parent_germplasm_id == ID),relationship = 'many-to-one') |>
    dplyr::rename(female_parent_germplasm_idN=newID) |>
    left_join(pedln$xrf,by=join_by(MIN_male_parent_germplasm_id==ID),relationship = 'many-to-one') |>
    dplyr::rename(male_parent_germplasm_idN=newID)
  #df_raw |> View()
  # make A inverse Additive matrix 
  Cov_m <- MCMCglmm::inverseA(pedln$newped)$Ainv
  Cov_m <- Matrix::forceSymmetric(Cov_m)
  colnames(Cov_m) <- rownames(Cov_m) <- pedln$newped$ID
  
  rinf <- apply(Cov_m,2,function(x){(sum(x!=0)-1)}) ## number of relationship in CovMatrix
  rinf <- data.frame(ID=as.integer(names(rinf)),rinf=rinf)
  
  png(paste0(debug_folder,'\\out_COV_Matrix.png'),w=1000,h=1000)
  print(image(Cov_m))
  dev.off()
  
  # Capture Germplasm Characteristic 
  unique_ids <- as.integer64(unique(c(pedln$xrf$ID,df_raw$MIN_fts_germplasm_id,df_IDmissing$fts_germplasm_id)))
  df_add <- source_gc_append(fts_germplasm_id=unique_ids,traits=unique(c('PRODSTG',traits_append)))
  df_add <- df_add |> 
    left_join(df_ID |> dplyr::select(fts_germplasm_id,KEY_ID),
              join_by(fts_germplasm_id == fts_germplasm_id),
              relationship = 'many-to-one',multiple = 'first')
  #df_add |> View()
  
  # Clear object not necessary 
  rm(df_ped,df_pedsel,ped_all,pedl,
     df_raw01,df_raw02,df_raw03,df_IDmissing,
     df_ID_female,df_ID_male,df_ID_min,df_ID_min2,df_ID_min3)
  
  # Run model 
  cl <- makePSOCKcluster(num_ParRun,
                         outfile = paste0(debug_folder,'\\out_debug.txt'),
                         revtunnel = TRUE)
  registerDoParallel(cl)
  result_list <- foreach::foreach(i=unique(df_raw$TRAIT),
                                  .inorder = FALSE,
                                  .verbose = TRUE,
                                  .packages =c('tidyverse','INLA'),
                                  .export = ls(globalenv())) %dopar% {
                                    tryCatch({
                                      ### Filter data 
                                      T_i <- Sys.time()
                                      df_model <- df_raw |> 
                                        dplyr::filter(TRAIT == i) |> 
                                        dplyr::select(fts_germplasm_idN,
                                                      season,field_name,observation_trait_val) |> 
                                        dplyr::rename(genof = fts_germplasm_idN) |> 
                                        dplyr::mutate(genof = as.character(genof),
                                                      season_field_name = paste(season,field_name,sep = '_')) |> 
                                        dplyr::select(-c(season,field_name))
                                      ### Model 
                                      Z <- makeZadd(data=df_model,Cov_m=Cov_m)
                                      df_model$idz <- 1:nrow(df_model)
                                      
                                      ### Standarized data  
                                      y_sd <- sd(df_model$observation_trait_val)
                                      y_mean <- mean(df_model$observation_trait_val)
                                      df_model$y_std <- (df_model$observation_trait_val-y_mean)/y_sd
                                      
                                      ### prior
                                      lgprior <- list(prec = list(prior = "loggamma", param = c(5 ,1),
                                                                  initial = 1))
                                      
                                      ### Model 
                                      fit <- INLA::inla(y_std ~ 1 +
                                                          f(season_field_name,model='iid',constr=T, hyper = lgprior) +
                                                          f(idz,model='z',Z=Z,Cmatrix = Cov_m,constr=T,hyper = lgprior) 
                                                        #  f(genof,model='iid',constr=T, hyper = lgprior)
                                                        ,
                                                        data=df_model,
                                                        verbose=T,
                                                        only.hyperparam = F,
                                                        control.compute=list(dic=F,config=TRUE),
                                                        control.predictor=list(compute = TRUE),
                                                        num.threads = 8)
                                      ##print(summary(fit))
                                      ### Extract data model 
                                      int <- fit$summary.fixed[,1:3]
                                      int$mean <- int$mean + y_mean
                                      int$PEV <- int$sd^2 * y_sd
                                      
                                      inb <- fit$summary.random$idz[-(1:nrow(df_model)),1:3]
                                      inb$ID <- as.integer(colnames(Cov_m))
                                      inb$mean <- inb$mean * y_sd
                                      inb$PEV <- (inb$sd)^2 * y_sd
                                      Va <- extract_trans_to_var(fit,precision="Log precision for idz") * y_sd #* 4 
                                      Ve <- extract_trans_to_var(fit,precision="Log precision for the Gaussian observations") * y_sd
                                      h2 <- (Va)/(Va+Ve)
                                      
                                      ### Inbreds 
                                      df_x <- df_model |> 
                                        group_by(genof) |> 
                                        summarise(pheno_mean=mean(observation_trait_val),
                                                  pheno_sd =sd(observation_trait_val),
                                                  pheno_n=length(observation_trait_val))
                                      res_inb <- data.frame(trait=i,
                                                            ID=inb$ID,
                                                            u=inb$mean,
                                                            PEV=inb$PEV,
                                                            Intercept=int$mean,
                                                            Intercept_PEV=int$PEV)
                                      
                                      res_inb <- res_inb |> 
                                        left_join(rinf,by=join_by(ID==ID)) |> 
                                        relocate(rinf,.after=ID)
                                      res_inb$BLUP <- rowSums(res_inb[,c('Intercept','u')],na.rm = TRUE)
                                      res_inb$Accuracy <- sqrt(1 - (res_inb$PEV/Va))
                                      res_inb$Reliability <- 1 - (res_inb$PEV/Va)
                                      res_inb$Heritability_ns <- h2
                                      res_inb <- merge(res_inb,df_x,all.x=TRUE,by.x='ID',by.y='genof')
                                      Time <- difftime(Sys.time(), T_i,units = "secs")
                                      cat(paste('at:',Sys.time(),'## Done: ',i,'in :',round(Time,2),'secs\n'))
                                      out <- list(i = list(res_inb=res_inb,fit=fit,Time=Time))
                                      names(out) <- i
                                      return(out)
                                    },
                                    error =  function(e,...){
                                      out <- list(list(res_inb=NA,fit=NA,Time=NA))
                                      names(out) <- i
                                      return(out)
                                    }
                                    )
                                  }
  parallel::stopCluster(cl)
  res_out_inbt <- res_out_fit <- res_out_time <- list()
  con <- 1
  model_name <- c()
  for(i in 1:length(result_list)){
    res_out_inbt[[con]] <- result_list[[i]][[1]]$res_inb
    res_out_fit[[con]] <- result_list[[i]][[1]]$fit
    res_out_time[[con]] <- result_list[[i]][[1]]$Time
    model_name <- c(model_name, names(result_list[[i]]))
    con <- con+1
  }
  
  res_out_inb <- do.call(rbind,res_out_inbt)
  res_out_inb <- res_out_inb |> dplyr::filter(!is.na(trait))
  
  # Export time consumed 
  res_out_Timep <- data.frame(Trait=model_name)
  for(i in 1:length(res_out_time)){
    x <- abs(as.numeric(res_out_time[i]))
    res_out_Timep$time_seg[i] <- x
    res_out_Timep$time_HHMMSS[i] <-  formatSeg(x)
  }
  
  # Model fitted plot diagnostic 
  names(res_out_fit) <- model_name
  pdf(paste0(debug_folder,'\\out_Diagplot.pdf'),h=7,w=12)
  for (i in 1:length(res_out_fit)){
    if(!is.na(res_out_fit[i])){
      plot(0,0,axes = F,main= model_name[i])
      try(plot(res_out_fit[[i]],
               single = T,
               plot.prior = T))
    }
  }
  dev.off()
  
  # Merge original id 
  res_out_inb <- res_out_inb |> 
    left_join(pedln$xrf,join_by(ID == newID), multiple = 'first', relationship = 'many-to-one')|> 
    dplyr::rename(fts_germplasm_id = ID.y)
  
  res_out_inb <-  res_out_inb |> 
    left_join(df_add,join_by(fts_germplasm_id == fts_germplasm_id),
              multiple = 'first', relationship = 'many-to-one') |>
    dplyr::filter(!is.na(pedigree)) |> 
    relocate(trait,ID,pedigree,origin,germplasm_lineage_identifier) |> 
    dplyr::select(-any_of('ID'))
  
  res_out_inb <-  res_out_inb |> dplyr::filter(line_type != 'Hybrid')
  
  df_rawS <- merge(df_rawS,
                   res_out_inb |> dplyr::select(trait,Heritability_ns) |> unique()
                   ,by.x='TRAIT',by.y='trait',all.x=TRUE)
  
  # Plot Prediction 
  res_out_inb |> ggplot(aes(x=pheno_mean,y=BLUP))+
    geom_point()+
    facet_wrap(~trait, scales='free')+
    geom_abline(intercept = 0, slope =1)+
    theme_bw()
  ggsave(paste0(debug_folder,'\\out_pred_inbred.png'),w=15,h=15)
  
  # Save results 
  wb <- createWorkbook()
  addWorksheet(wb, "SummaryRawData_Pred")
  writeData(wb, 1, x = df_rawS, withFilter = TRUE, na.string='')
  addWorksheet(wb, "Results_Inbred")
  writeData(wb, 2, x = res_out_inb, withFilter = TRUE, na.string='')
  addWorksheet(wb, "Time_run_model")
  writeData(wb, 3, x = res_out_Timep, withFilter = TRUE, na.string='')
  addWorksheet(wb, "Raw_data")
  writeData(wb, 4, x = df_raw, withFilter = TRUE, na.string='')
  saveWorkbook(wb, file = paste0(file_out,".xlsx"), overwrite = TRUE)
  return(list('df_rawS'=df_rawS,
              'res_out_inb'=res_out_inb,
              'res_out_Timep' = res_out_Timep,
              'df_raw'=df_raw))
}

#' Source rawdata from CSW to Hybrid Analysis Pipeline
#' 
#' @title 
#' @param WHERE filter to apply in WHERE SQL query from CSW in table bcs-veg-mart.performance.veg_experiment_obs_entries
#' @details the columns mandatory for query is already fixed
#' @return data frame with raw data with no treatment from CSW
#' 

## Make augmented incidence matrix
makeZadd <- function(data=df_model,Cov_m=Cov_m){
  Z <- model.matrix(~genof-1,data)
  colnames(Z) <- gsub('genof','',colnames(Z))
  Cov_mf <- colnames(Cov_m)[!(colnames(Cov_m) %in% colnames(Z))]
  AdM <- matrix(0,ncol=length(Cov_mf),nrow=nrow(Z))
  colnames(AdM) <- Cov_mf
  Z <- cbind(Z,AdM)
  Z <- Matrix::Matrix(Z[,colnames(Cov_m)],sparse = TRUE)
  return(Z)
}