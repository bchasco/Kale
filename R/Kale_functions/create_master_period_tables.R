# Create a look-up table (LUT) of sample periods/dates
create_master_period_tables<-function(dat, period_method = "week", date_capture_column, date_recapture_column = NULL, year_column){

  # Identify the number of unique return years in the survey dataset
  unique_years<- dat %>%  distinct(.data[[year_column]]) %>% pull(.data[[year_column]])

  # Create "LUT_dates_by_week_init" depending on period_method
  if(period_method == "week"){
    # Create two blank tibbles to fill
    periods_by_date <- tibble(Date = as.Date(character()),Period = integer())
    periods_by_week <- tibble(Year = integer(), Period = integer(), week = integer(), n_surveys = numeric(), Date_Min = as.Date(character()), Date_Max = as.Date(character()), Date_Mid = as.Date(character()))
    #Loop through each unique year
    for(year in 1:length(unique_years)){
      dat_year<-dat |>
        filter(.data[[year_column]] == unique_years[year])

      Date_Captures<-
        dat_year |>
        distinct(Date = .data[[date_capture_column]])

      if(is.null(date_recapture_column)==FALSE){
        Date_Recaps<-
          dat_year |>
          distinct(Date = .data[[date_recapture_column]])
      }else{
        Date_Recaps<-tibble(Date = as.Date(character()))
      }

      dates_distinct<-
        rbind(Date_Captures, Date_Recaps) %>%
        drop_na |>
        distinct(Date) |>
        arrange(Date)

      LUT_dates_by_week_init<-
        dates_distinct|>
        mutate(week = as.integer(format(Date, "%U"))) |>
        group_by(week) |>
        summarise(
          n_surveys = as.numeric(n_distinct(Date)),
          Date_Min = min(Date),
          Date_Max = max(Date)
        ) |>
        ungroup() |>
        mutate(Year = as.integer(unique_years[year]), Date_Mid = if_else(Date_Max == Date_Min, Date_Min, Date_Min + round((Date_Max - Date_Min),0) - 1 )) |>
        arrange(Date_Mid) |>
        rowid_to_column(var = "Period") |>
        select(Year, Period, everything())

      LUT_periods_by_date<-
        dates_distinct %>%
        rowwise() %>%  # Ensures row-wise filtering for each Date
        mutate(Period = LUT_dates_by_week_init$Period[which(Date >= LUT_dates_by_week_init$Date_Min & Date <= LUT_dates_by_week_init$Date_Max)][1]) %>%
        ungroup()

    # Add year specific periods by date to "master" tables
      periods_by_week<-rbind(periods_by_week, LUT_dates_by_week_init)
      periods_by_date<-rbind(periods_by_date, LUT_periods_by_date)
    }
  }else{
    if(period_method == "manual_LUT"){
      LUT_dates_by_week_init<-
        LUT$period_by_date_manual |>
        filter(ui_return_year == ui_return_year, ui_project_name == ui_project_name) |>
        group_by(Period = Period_Num) |>
        summarise(n_surveys = n(), Date_Min = min(Survey_Date), Date_Max = max(Survey_Date)) |>
        ungroup() |>
        mutate(
          week = as.numeric(format(Date_Min, "%U"))
          , Date_Mid = if_else(Date_Max == Date_Min, Date_Min, Date_Min + round((Date_Max - Date_Min),0) - 1 )
        ) |>
        arrange(Date_Mid)

      periods_by_date<-
        LUT$period_by_date_manual |>
        select(Date = Survey_Date, Period = Period_Num)

    }else{
      print("poop")
    }
  }

  # Create "LUT_dates_by_week_final"
  LUT_dates_by_week_final<-
    periods_by_week |>
    mutate(
      Days_Btw_Periods = as.numeric(Date_Mid - lag(Date_Mid))
      , MidDate_Btw_Periods = as.Date((Date_Mid + round(lead(Days_Btw_Periods)/2, 0) ))
    ) |>
    mutate(Year_Period = paste(Year, Period, sep="_"), MidDate_Btw_Periods = if_else(is.na(MidDate_Btw_Periods)==TRUE, Date_Mid + median(Days_Btw_Periods, na.rm=TRUE)/2, MidDate_Btw_Periods)) |>
    select(Year, Period, Year_Period, week, n_surveys, Date_Min, Date_Mid, Date_Max,  Days_Btw_Periods, MidDate_Btw_Periods)

  # Return final LUT
  return(list(
    master = LUT_dates_by_week_final
    , periods_by_date = periods_by_date
  ))

}
