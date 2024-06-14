library(testthat)
library(here)
library(tidyverse)


# source file containing functions
path <- here("code")
source(paste(path, "model_fitting_functions.R", sep='/'))

# read in test data
data_path <- here("data") 
aus_deaths <- read.table(paste(data_path, "AUSdeath.txt", sep='/'), sep=',', header = TRUE)
aus_exposure <- read.table(paste(data_path, "AUSExposures_1x1.txt", sep='/'), sep='', 
                           header = TRUE, stringsAsFactors = FALSE, skip = 2)

aus_expo_long <- pivot_longer(aus_exposure, cols = c('Female', 'Male', 'Total'),
                              names_to = "Sex", values_to = "Exposure")

aus_expo_long <- aus_expo_long %>%
  mutate(Sex = ifelse(Sex == "Female", "f", 
                      ifelse(Sex == "Male", "m", 
                             ifelse(Sex == "Total", "b", Sex)))) %>%
  mutate(Age = ifelse(Age == "110+", "110", Age)) %>%
  mutate(Age = as.integer(Age)) %>%
  filter(!Sex == 'b')

aus_deaths <- aus_deaths |> 
  filter(LDB == 1) |>
  select(Year, Age, Sex, Deaths) |>
  filter(Age != "TOT" & Age != "UNK") |>
  mutate(Age = as.integer(Age))

combined_df <- aus_expo_long |>
  left_join(aus_deaths, by = c("Year", "Age", "Sex")) |>
  mutate(Deaths = as.integer(Deaths))

test_data <- combined_df |>
  filter(Age >= 80) |>
  replace_na(list(Deaths = 0))

# tests
test_that("log_lik_fn works correctly", {
  params <- c(0.1, 0.2)
  result <- log_lik_fn(params, test_data, kannisto_fn)
  
  expect_true(is.numeric(result))
  expect_false(is.na(result))
  expect_gt(result, -Inf)
})

test_that("get_initial_params works correctly", {
  initial_params <- get_initial_params(log_lik_fn = log_lik_fn, data = test_data,
                                      fn = kannisto_fn)
  
  expect_true(is.numeric(initial_params))
  expect_length(initial_params, 2)
  expect_named(initial_params, c("a", "b"))
})

test_that("kannisto_fn works correctly", {
  params <- c(0.1, 0.2)
  result <- kannisto_fn(params, test_data)
  
  expect_true(is.numeric(result))
  expect_length(result, nrow(test_data))
})

test_that("kannisto_get_grad works correctly", {
  params <- c(0.1, 0.2)
  result <- kannisto_get_grad(params, test_data, kannisto_fn)
  
  expect_true(is.numeric(result))
  expect_length(result, 2)
})

test_that("kannisto_get_mx works correctly", {
  a <- 0.1
  b <- 0.2
  result <- kannisto_get_mx(a, b)
  
  expect_true(is.matrix(result))
  expect_equal(ncol(result), 2)
  expect_equal(nrow(result), 31)
})

test_that("fit_kannisto works correctly", {
  result <- fit_kannisto(get_initial_params, log_lik_fn, kannisto_get_grad, kannisto_fn,
                         test_data)
  
  expect_true(is.matrix(result))
  expect_equal(ncol(result), 2)
  expect_gt(nrow(result), 0)
})



