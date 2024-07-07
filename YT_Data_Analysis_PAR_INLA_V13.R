### Loading Yield test data and pedigree information to calculated GCA, SCA and 
### BLUP for hybrid performance 
### Eder Silva, 06/06/2024
# Loading library-------------------------------------------------------------
rm(list=ls())
pak <- c('bigrquery','lubridate','openxlsx','readxl','AGHmatrix',
         'bit64','foreach','doParallel','quantreg','INLA','tidyverse',
         'snpReady','ASRgenomics','MCMCglmm')
#install.packages(pak,dep=TRUE);
#install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
#if (!require("BiocManager", quietly = TRUE));install.packages("BiocManager");BiocManager::install("impute")
pakl <- unlist(lapply(pak, require, character.only = TRUE,quietly=TRUE));pak[!pakl]
bq_auth(email='eder.silva@bayer.com')

# active auxiliary functions ---------------------------------------------------
source('fun/aux_functions.R')

# Append GE e general parameters ---------------------------------------------
sourceNewData <- T
force_DQ <- T

# Load data YT Data ------------------------------------------------------------
if(sourceNewData){
  WHERE <- "
  veoe.crop IN ('Corn','Sweet Corn')
  AND veoe.country = 'Brazil'
  AND veoe.line_type IN ('Hybrid','HYBRID')
  AND veoe.planting_year >= 2021
  ##AND  COALESCE(veoe.block_name,veoe.set_name) IN ('24TSYP2_BR_L01','24TSYP2_BR_L02','24TSYP2_BR_L03','24TSYP2_BR_L04','24TSYP2_BR_L05','24TSYP2_BR_L06','24TSYP3_BR_L01','24TSYP3_BR_L02','24TSYP3_BR_L03','24TSYP3_BR_L04','24TSYP3_BR_L05','24TSYP3_BR_L06')
  AND COALESCE(primary_question_cd,observation_reference_cd) IN ('CMNTS',
  'STLP','RTLP','HVPOP','TDPP','KRP','FKCY','FY','POLSG','SLKEG','ASI50D','ERPOP',
  'CKP','PCK','PGW','GWP','POLSD','SLKED','EHT','DSR','PHT','HSC','RTLC',
  'STLC','TFIR','EARD','EARL','ENDC','KRN','FINAL','FLTXR','HSKER','SCCSA')
  "
  df_raw <- source_raw_data(WHERE)
  save(df_raw,file='input_files/in_Data.Rdata')
}else{
  load('input_files/in_Data.Rdata')
}

## Force remove observation ID -------------------------------------------------
if(force_DQ){
  observation_to_remove <- read_excel('input_files/in_DQ.xlsx')
  df_raw <- df_raw |> 
    filter(!(observation_id %in% observation_to_remove$observation_id))
}

## Split Comments and string traits --------------------------------------------
extract_string_traits(data = df_raw, traits_char = c('CMNTS','ENDC'))

#### Temporary fix traits for TSY 2024 ----------------------------------------- 
df_raw <- df_raw |> 
  mutate(observation_trait_val = as.numeric(observation_trait_val)) |> 
  mutate(TRAIT = case_when(TRAIT == 'PCK' ~ 'CKP',
                           TRAIT == 'PGW' ~ 'GWP',
                           TRAIT == 'SCCSA' ~ 'FKCY',
                           .default = TRAIT),
         observation_trait_val = case_when(TRAIT %in% c('ERPOP','HVPOP') & season =='Fall' ~ observation_trait_val / 1000,
                                           .default = observation_trait_val))
write.csv(df_raw,'df_raw.csv',row.names = F)
# SLKED
# SLKEG
# POLSD
# POLSG
# Duplucateds 


## Run model -------------------------------------------------------------------
##minGen = 5;UseGenomeData = F;SNP_source = c('raw','imputed')[2];Inb_hyb_append = read_excel('input_files/in_Append_Inbred.xlsx')[1:10,];Use10genPed = T;Exclude_Hybrid_GCA = T;file_out = 'out_YT_Results_PAR';debug_folder = 'out_DEBUG';traits_append = c('HETGP','P50','S50','HBC','LID','P1DATE','P2DATE','PMC','PRODSTG','ENC');num_ParRun = 1
outModel <- run_model(
  df_raw = df_raw,
  minGen = 30,
  UseGenomeData = F, 
  SNP_source = c('raw','imputed')[2],
  Inb_hyb_append = read_excel('input_files/in_Append_Inbred.xlsx')[1:3,],
  Use10genPed = T,
  Exclude_Hybrid_GCA = T,
  file_out = 'out_YT_Results_PAR',
  debug_folder = 'out_DEBUG',
  traits_append = c('HETGP','P50','S50','HBC','LID','P1DATE','P2DATE','PMC','PRODSTG','ENC'),
  num_ParRun = 1,
  Zscore_filter = 4)
outModel

## Adjust PRODSTAGE ------------------------------------------------------------
adj_prodstg(out_res_file='out_YT_Results_PAR.xlsx',
            in_prod_file='input_files/in_PRODSTG.xlsx')

## Run TIX ---------------------------------------------------------------------
runTIX(data = read_excel('out_YT_Results_PAR.xlsx','Results_Hybrid'),
       traits_TIX = read_excel('input_files/in_TIX_traits.xlsx'),
       traits_app = read.csv('input_files/out_Traits_Char.csv'),
       f_stages = c('SC2','PS1','PS2','PS3','PS4','CM','CHECK'),
       out_file = "post_Make_TIX/out_PivotData_BLUP_TIX_run.xlsx")

## Run advance graphics --------------------------------------------------------
makePLOTS(in_parp = "input_files/in_Graphics_Parameters.xlsx",
          dataPp = "post_Make_TIX/out_PivotData_BLUP_TIX_run.xlsx",
          out_folder = "post_ADV_plots",
          selrows = NULL)

## Run Hybrid mining -----------------------------------------------------------
runHM(data = read_excel('out_YT_Results_PAR.xlsx','Results_Inbred'),
      inbred = read_excel('input_files/in_Inbreds_to_HM.xlsx'),
      traits = read_excel('input_files/in_TIX_traits.xlsx'),
      CheckInventory = T,
      owner_breeding_prog = 'AV',
      out_res_folder = 'post_HM',
      GCT_append = c('P50','S50','HETGP','ENC'))

