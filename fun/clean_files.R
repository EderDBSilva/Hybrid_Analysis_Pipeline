
clear_files <- function(...){
  file.remove(list.files('out_DEBUG', include.dirs = F, full.names = T, recursive = T))
  file.remove('out_YT_Results_PAR.xlsx')
  file.remove('post_Make_TIX/out_PivotData_BLUP_TIX_run.xlsx')
  file.remove('post_Make_TIX/out_plotTIX.pdf')
  file.remove('post_HM/out_GCA_Inbred.xlsx')
  file.remove('post_HM/out_HM_Prediction.xlsx')
  file.remove('post_HM/out_Parents_notFound_inGCA.csv')
  file.remove('post_ADV_plots/ADV_Meeting_graphics.html')
  file.remove('post_ADV_plots/ADV_Meeting_graphics.Rmd')
  file.remove('out_DEBUG/out_Traits_Char.csv')
  file.remove('input_files/in_Data.Rdata')
  file.remove(list.files('post_ADV_plots/out_DEBUG/', include.dirs = F, full.names = T, recursive = T))
  file.remove(list.files('post_ADV_plots/out_plot/', include.dirs = F, full.names = T, recursive = T))
}

clear_files()
