#---------------------------------------------------------------------------------------------------------- -
# LOAD FUNCTIONS, PACKAGES, and look-up tables (LUTs)                                                    ----
#---------------------------------------------------------------------------------------------------------- -
# Load functions
sapply(FUN = source, paste(getwd(), "R/Kale_functions", list.files("R/Kale_functions"), sep="/"))

# Install/Load packages
package_list<-c("tidyverse", "glue")
install_or_load_pack(package_list)

# Load look-up tables
LUT<-load_LUTs(folder_path = "R/Kale_luts")
  LUT$carc_group_options
  LUT$carc_strata_by_sex_grouping

# Specify carcass grouping
ui_carcass_grouping<-c(3)
  LUT$carc_strata_by_sex_grouping |> filter(Grouping_Strategy_Num == ui_carcass_grouping)

#-----------------------------------------------------------------------------------------------------
# Import .csv datasets of interest
#-----------------------------------------------------------------------------------------------------
## Define filepath of datasets
folder_path <- "Data/Originals/NF_Lewis"
## Define columns to keep
columns_to_keep <- c("Return_Yr", "SPECIES", "Run", "TagDate", "TagReach",
                     "Tag1", "Tag2", "Sex", "FL", "Mark", "ScaleAge",
                     "RecapDate", "RecapTag1", "RecapTag2", "RecapReach", "SampleRate", "Comments")
# Import
combined_data <-
  import_csv_files(folder_path, columns_to_keep)

# View the structure of the final dataset
glimpse(combined_data$final_df)

#-----------------------------------------------------------------------------------------------------
# Review imported data
#-----------------------------------------------------------------------------------------------------
## Check for column mismatches (version 1)
column_mismatches <- map_df(list.files(folder_path, pattern = "\\.csv$", full.names = TRUE), ~{
  df <- read_and_standardize_csv(.x)
  tibble(file = basename(.x), columns = paste(colnames(df), collapse = ", "))
})
print(column_mismatches, width=Inf)


## Check for column mismatches (version 2)
column_sets <- map(combined_data$file_list, ~ colnames(read_and_standardize_csv(.x)))
column_differences <- map(column_sets, ~ setdiff(.x, column_sets[[1]])) # Identify the set of columns that are not in the first file

mismatch_summary <- tibble( # Combine results into a tibble for easier review
  file = basename(combined_data$file_list),
  extra_columns = map_chr(column_differences, ~ paste(.x, collapse = ", "))
)
print(mismatch_summary, n = Inf, width = 10000)

#-----------------------------------------------------------------------------------------------------
# Filter data set
#-----------------------------------------------------------------------------------------------------
##Species
combined_data$final_df |> distinct(SPECIES)
filt_species<-c("Chinook salmon")

dat_filt<-
  combined_data$final_df |>
  rename(Species = SPECIES) |>
  filter(is.na(Species)==FALSE & Species %in% filt_species)

#-----------------------------------------------------------------------------------------------------
# Format data set
#-----------------------------------------------------------------------------------------------------
dat_format_init<-
  dat_filt |>
  mutate(
      Location = "NF_Lewis"
    , Date_Maiden = dmy(TagDate)
    , Date_Recap  = dmy(RecapDate)
    , TagState = str_extract(TagReach, "\\d+")
    , RecapState = str_extract(RecapReach, "\\d+")
    , check = if_else((is.na(TagState)==FALSE & is.na(RecapState)==FALSE & TagState>RecapState), -1, 1)
  ) |>
  select(-TagDate, -TagReach, -RecapDate, -RecapReach) |>
  select(Location, Return_Yr, Species, Run, Sex, FL, Mark, ScaleAge, everything()) |>
  mutate(Tagged = if_else(is.na(Tag1)==FALSE | is.na(Tag2)==FALSE, 1, 0),                  # **tagged?** (where 1 = tagged; 0 = not tagged)
         Recap  = case_when(
           Tagged == 1 & (is.na(RecapTag1)==FALSE | is.na(RecapTag2)==FALSE) ~ 1, # **recapped?**  (where 1 = tagged & recapped, 0 = tagged but not recapped, NA = not tagged)
           Tagged == 1 & (is.na(RecapTag1)==TRUE & is.na(RecapTag2)==TRUE) ~ 0,
           Tagged == 0  ~ NA_real_
  ))

# Preview formatted data
dat_format_init |> count(Tagged, Recap)


#------------------------------------------------------------------------------------------------------------------------- -
# create master LUT of sample periods/dates  ----
#------------------------------------------------------------------------------------------------------------------------- -
# Run LUT data/period function
period_table<-
  create_master_period_tables(
      dat = dat_format_init
    , period_method = "week" # options: "week" (survey dates converted to periods based on survey week; default) or "manual_LUT" (user defines the periods for a set of distinct survey dates for a ui_project_name in "luts/LUT_periods_by_survey_date_manual.csv" )
    , date_capture_column = "Date_Maiden"
    , date_recapture_column = "Date_Recap"
    , year_column = "Return_Yr"
    )
period_table$master
period_table$periods_by_date

dat_format_periods<-
  dat_format_init |>
  rowid_to_column(var = "Index") |>
  left_join(period_table$periods_by_date, by = c("Date_Maiden" = "Date")) |>     # **period tagged**
  rename(Period_Cap = Period) |>
  mutate(Year_Period = paste(Return_Yr, Period_Cap, sep="_")) |>
  mutate(days_between_capture = if_else(Recap == 1, as.numeric(Date_Recap - Date_Maiden), NA_real_ ))

#--------------------------------------------------------------------------------------------------------------------------------------------------------
# supplemental section to account for jacks not being clearly identified by field crews (NOTE: protocols were updated to address this issue in 2022) ----
#--------------------------------------------------------------------------------------------------------------------------------------------------------
jack.cut.off<-59
Max_Rate_by_Period <- dat_format_periods |> mutate(SampleRate = replace_na(SampleRate, "1")) |> group_by(Year_Period) |> summarise(Max_Rate = as.numeric(max(SampleRate, na.rm = TRUE)))
dat_format_final<-
  dat_format_periods |>
  left_join(Max_Rate_by_Period, by = c("Year_Period")) |>
  #left_join(dat |> select(Index, Comments)) |>
  mutate(Sex_Temp = case_when(
    #mutate(Sex.NEW = case_when(
    is.na(Sex) == TRUE | Sex == "Adult" ~ "Adult",
    Sex == "Unknown" | Sex == "unknown" | Sex == "UNK" | Sex == "Unk" ~ "Adult",
    Sex == "Female" ~ "Female",
    Sex == "Jack" ~ "Jack",
    Sex == "Male" & SampleRate == 1 & Max_Rate > 1 ~ "Jack",
    Sex == "Male" & is.na(Comments) == FALSE & (Comments == "JACKS" | Comments == "Jacks" | Comments == "jacks" | Comments == "jack card") ~ "Jack",
    Sex == "Male" & SampleRate > 1 ~ "Male",
    Sex == "Male" & FL > jack.cut.off ~ "Male",
    Sex == "Male" & Max_Rate == 1 & FL <= jack.cut.off ~ "Jack",
    Sex == "Male" & Max_Rate == 1 & is.na(FL)==TRUE ~ "Male",
    Sex == "Female" & is.na(FL)==TRUE ~ "Female",
    Sex == "Jack" & Max_Rate == 1 & FL > jack.cut.off ~ "Male",
    TRUE ~ "Need.Another.Category"
  )
  )
dat_format_final |> group_by(Return_Yr) |> count(Sex_Temp) |> pivot_wider(names_from =  Sex_Temp, values_from = n)


#------------------------------------------------------------------------------------------------------------------------- -
# execute project specific data updates and generate final data set  ----
#------------------------------------------------------------------------------------------------------------------------- -
# Summarize carcass recoveries by sex proportions & capture period
  ui_sex_bio<-"tagged" # enter either "tagged" (only use sex calls from tagged fish) or "all" (use sex calls from any bio-sampled carcass)
  ui_sex_grouping<- "female_male" # enter either "jack_female_male" or "female_male" (choice depends on how fish were sub-sampled for sex)
  sex_props<-
      summ_sex_ratios(
          dat= dat_format_final
        , bio_data = ui_sex_bio
        , groups = ui_sex_grouping
        , periods = period_table$master)
  # NOTE: this summary table is used below to probabilistically assign a Sex(_final) designation to every carcass with a missing -- i.e., NA -- sex call; this step is needed to run the J-S as every carcass must have a sex call)
  sex_props

# create two LUTs based on carcass grouping options & ui_carcass_grouping
  (LUT_strata<-create_LUT_strata_sex(LUT_input = LUT$carc_strata_by_sex_grouping, ui_input = ui_carcass_grouping))

# **create** object "dat_sex_interpolated" which will be used to generate data objects used in abundance (J-S) model
  ui_seed_num = 298
  dat_sex_interpolated<-
    dat_format_final |> # Step 1 - create "Sex_final" column (see notes below)
    left_join(sex_props$sex_prop_period |> select(Year_Period, starts_with("Prop_")), by=c("Year_Period") ) |>
    rowwise() |>
    mutate(Prop_Sex_Jack = ifelse("Prop_Sex_Jack" %in% names(sex_props$sex_prop_period), Prop_Sex_Jack, 0)) |>
    rowwise() |>  # Step 1a - generate a "Sex_Interpol" for every carcass using the "sample_reprod" function; sample_reprod uses the period specific sex ratio
    mutate(Sex_Interpol = sample_reprod(prob = c(Prop_Sex_Male, Prop_Sex_Female, Prop_Sex_Jack), x = c("Male", "Female", "Jack"), size = 1, seed = ui_seed_num+Index),
           Sex_final = case_when( #Step 1b - assign "Sex_final" for every carcass...
             Sex_Temp == "Male" | Sex_Temp == "Female" | Sex_Temp == "Jack" ~ Sex_Temp, #...if Sex_Temp is designated as Jack, Female, or Male (i.e., designated in the field), assign Sex_final equal to Sex_Temp
             TRUE ~ Sex_Interpol #...if Sex_Temp is not designated as Jack, Female, or Male (i.e., sex not designated in the field), assign Sex_final equal to Sex_Interpol
           )
    )

# Generate final dataset
  dat_final<-
    dat_sex_interpolated  |>
    select(-starts_with("Prop_Sex")) |>
    ungroup() |>
    left_join(LUT_strata$by_sex , by=c("Sex_final")) |>
    filter(days_between_capture>0|is.na(Date_Recap)==TRUE) |>
    filter(is.na(Strata_Name)==FALSE) |>
    select(Location, Return_Yr, Period_Cap,	Year_Period, Species, Run, Sex_original = Sex, FL, Mark, ScaleAge, Tag1, Tag2, RecapTag1, RecapTag2, Date_Maiden,	Date_Recap, TagState,	RecapState, Sex_final, check) |>
    print(width=Inf)

#-----------------------------------------------------------------------------------------------------
# final data checks
#-----------------------------------------------------------------------------------------------------
# Total checks == -1
dat_final |>
  filter(check == -1) |>
  count()

# Break down of -1s by Year and Reach
dat_final |>
  filter(check == -1) |>
  group_by(Return_Yr, TagState, RecapState) |>
  summarise(n = n()) |>
  pivot_wider(names_from = RecapState, values_from = n) |>
  arrange(Return_Yr, TagState) |>
  print(n=Inf)

#-----------------------------------------------------------------------------------------------------
# export data
#-----------------------------------------------------------------------------------------------------
n_years<-dat_final |> distinct(Return_Yr) |> count() |> pull()
filepath_export<-"Data"
filename_export<-glue("NF_Lewis_combined-{min(range(dat_final$Return_Yr))}_{max(range(dat_final$Return_Yr))}-{Sys.Date()}.csv")
write.csv(x = dat_final, file = glue("{filepath_export}/{filename_export}"), row.names = FALSE)
