#===================================================================================
# This file loads the look-up tables (luts) that are necessary to run analysis code
#===================================================================================
load_LUTs<-function(folder_path){
  # # Selectivity test groupings
  #     selectivity_groups<-
  #       read_csv(
  #         file = "luts/LUT_Selectivity_Test_Groups.csv",
  #         col_types = cols(
  #           .default = col_character()
  #         , Group_Number = col_double()
  #         )
  #       )
  #
  # # Abundance (JS) models
  #     model_files <-
  #       read_csv(
  #         file = "luts/LUT_abundance_JS_models.csv",
  #         col_types = cols(
  #           .default = col_character()
  #         , Model_Num = col_double()
  #         )
  #       )
  #
  # # Carcass conditions - numeric and alpha codes
  #     carc_cond<-
  #       read_csv(
  #         file = "luts/LUT_Carcass_Condition.csv",
  #         col_types = cols(
  #           .default = col_character()
  #         , CC_num  = col_double()
  #         )
  #     )

  # Carcass groupings options (i.e., options for stratifying carcass data; estimates of abundance are generated separately for each strata in J-S before being partitioned using biological data into other groupings of interest i.e., run, origin, and age)
      carc_group_options<-
        read_csv(
          file = glue("{folder_path}/LUT_carcass_grouping_options.csv"),
          col_types = cols(
            .default = col_character()
          , Grouping_Strategy_Num = col_double()
          , Num_of_strata = col_double()
          )
        )

  # Carcass strata_num and strata_names based on sex designation and grouping options
      carc_strata_by_sex_grouping<-
        read_csv(
          file = glue("{folder_path}/LUT_Carcass_Strata_by_Sex_for_Grouping_Options.csv"),
          col_types = cols(
              .default = col_character()
            , Grouping_Strategy_Num = col_double()
            , Strata_Num = col_double()
            )
       )

  # # Run_Sub classifications based on CWT recovery data
  #     run_sub_class<-
  #       read_csv(
  #         file = "luts/LUT_CWT_to_run_sub_classification.csv",
  #         col_types = cols(
  #           .default = col_character()
  #         )
  #       )
  #
  # # Run_Sub classifications based on CWT recovery data
  #     period_by_date_manual<-
  #       read_csv(
  #         file = "luts/LUT_periods_by_survey_date_manual.csv",
  #         col_types = cols(
  #           .default = col_double()
  #           , Survey_Date = col_date(format = "%d-%b-%y")
  #           , ui_project_name = col_character()
  #         )
  #       )

  return(list(  #selectivity_groups = selectivity_groups
              #, model_files = model_files
              #, carc_cond = carc_cond
              #,
              carc_group_options = carc_group_options
              , carc_strata_by_sex_grouping = carc_strata_by_sex_grouping
              #, run_sub_class = run_sub_class
              #, period_by_date_manual = period_by_date_manual
              )
         )

}
