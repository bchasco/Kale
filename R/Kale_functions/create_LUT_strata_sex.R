# create look-up table (LUT_input) of carcass strata definitions based on LUT_input & ui_input
create_LUT_strata_sex<-function(LUT_input, ui_input){
  
  by_sex<-LUT_input |> filter(Grouping_Strategy_Num == ui_input) |> select(-Grouping_Strategy_Num)
  by_distinct<-by_sex|> distinct(Strata_Name, Strata_Num)

  return(list(by_sex=by_sex, by_distinct=by_distinct))
} 