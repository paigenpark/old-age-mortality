# debugging cohort data
cohort = 1895
age = 24

raw_deaths <- fra_deaths_hmd |>
  filter(Year == 1919 | Year == 1920, Age == age, Sex == "f")
raw_deaths

raw_expo <- fra_expo_cohort_hmd |>
  filter(Cohort == cohort, Age == age)
raw_expo

filtered_deaths <- cleaned_deaths_list[[1]] |>
  filter(Year == 1919 | Year == 1920, Age == age, Sex == "f")
filtered_deaths

calculated_cohort_deaths <- dfs_with_cohort[[1]] |>
  filter(Year == 1919 | Year == 1920, Age == age, Sex == "f")
print(calculated_cohort_deaths,width = Inf)

calculated_cohort_deaths_and_expo <- hmd_cohort_data |>
  filter(Cohort == cohort, Age == age, Sex == "f", country == "fra")
print(calculated_cohort_deaths_and_expo, width = Inf)

calculated_cohort_deaths_and_expo_mx <- cohort_val_with_mx |>
  filter(Cohort == cohort, Age == age, Sex == "f")
print(calculated_cohort_deaths_and_expo_mx, width = Inf)

hmd_mx <- fra_cohort_mx |>
  filter(Year == cohort, Age == age)
hmd_mx

test <- combined_mx_cohort |>
  filter(Cohort == cohort, Age == age, Sex == "f", country == "fra")
print(test, width = Inf)

calculated_cohort_younger <- combined_mx_younger_cohort |>
  filter(Cohort == cohort, Age == age, Sex == "f", country == "fra")
print(calculated_cohort_younger, width = Inf)

# something might be off with combine_mx_death_data_coh function - 
# Mxs look similar when calculated manually, and difference is much 
# larger in the combined data frame - need to debug 

rows_with_na <- rowSums(is.na(combined_mx)) > 0
combined_mx[rows_with_na, ]

rows_with_na <- rowSums(is.na(combined_mx_cohort)) > 0
combined_mx_cohort[rows_with_na, ]

rows_with_na <- rowSums(is.na(period_val)) > 0
period_val[rows_with_na, ]

rows_with_na <- rowSums(is.na(cohort_val)) > 0
cohort_val[rows_with_na, ]

period_val$Age[is.na(period_val$Deaths)]




