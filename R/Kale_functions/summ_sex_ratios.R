# summarize sex ratio data (i.e., proportions)
summ_sex_ratios<-function(dat, bio_data, groups, periods){
  sex_filter<-
    dat |>
    filter(case_when(
      bio_data == "tagged" ~ Tagged == 1,
      TRUE ~ Index>0
    )) |>
    filter(case_when(
          groups == "jack_female_male" ~ Sex_Temp == "Jack" | Sex_Temp == "Female" | Sex_Temp == "Male",
          groups == "female_male" ~ Sex_Temp == "Female" | Sex_Temp == "Male"
    ))

  missing_periods<-
    periods |>
    filter(!Year_Period %in% unique(sex_filter$Year_Period))|>
    pull(Year_Period)

  sex_prop_total<-
    sex_filter |>
    group_by(Sex_Temp) |>
    summarise(Count = n(), .groups = "drop") |>
    mutate(Prop_Sex = Count/sum(Count)) |>
    mutate(Year = "Total", Period = "Total", Year_Period = "Total") |>
    pivot_wider(
      id_cols = c(Year, Period,Year_Period),
      names_from = Sex_Temp,
      values_from = c(Prop_Sex, Count),
      values_fill = list(Prop_Sex = 0)
    )

  sex_prop_period<-
    sex_filter |>
    group_by(Return_Yr, Period_Cap, Sex_Temp) |>
    summarise(Count = n(), .groups = "drop") |>
    group_by(Return_Yr, Period_Cap) |>
    mutate(Prop_Sex = Count/sum(Count)) |>
    pivot_wider(
      id_cols = c(Return_Yr, Period_Cap),
      names_from = Sex_Temp,
      values_from =  c(Prop_Sex, Count),
      values_fill = list(Prop_Sex = 0)
    ) |>
    rename(Year = Return_Yr, Period = Period_Cap) |>
    mutate(Year_Period = paste(Year, Period, sep="_"))

  if(length(missing_periods)>0){
    temp_row<-c()
    for(n in 1:length(missing_periods)){
      temp_row<-
        sex_prop_total |>
        mutate(Year = str_split(missing_periods[n], "_", simplify = TRUE)[, 1]
                 , Period = str_split(missing_periods[n], "_", simplify = TRUE)[, 2]
               , Year_Period = missing_periods[n])
        #mutate(Return_Yr = "TOTAL", Period_Cap = as.integer(gsub("_.*", "", missing_periods[n])))

      sex_prop_period<-
      sex_prop_period |>
        mutate(Period = as.character(Period)) |>
        bind_rows(temp_row)
    }
  }

  sex_prop_period<-
    sex_prop_period |>
    mutate(Year = as.integer(Year), Period = as.integer(Period)) |>
    arrange(Year, Period) |>
    ungroup()

return(list(sex_prop_period=sex_prop_period, sex_prop_total=sex_prop_total, missing_periods = missing_periods))

}
