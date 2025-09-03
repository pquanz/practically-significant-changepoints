library(dplyr)
library(tidyr)
library(zoo)

library(pscp)


load_data <- function(path) {
  tibble(read.csv(path)) %>%
    mutate(date = as.Date(date)) %>%
    arrange(date)
}

# set negative values to zero
# interpolate NA values and
# remove remaining (non-interpolated) NA values
clean_data <- function(data) {
  data %>%
    mutate(across(-date,~ ifelse(. < 0, 0, .))) %>%
    mutate(across(-date, ~na.approx(.x, x = date, na.rm = FALSE, rule = 2))) %>%
    select(where(~ !any(is.na(.x))))
}

data1 <- load_data("raw/pm10_data.csv")
data2 <- load_data("raw/pm25_data.csv")

data_pm <- full_join(data1, data2[, 2:ncol(data2)], by = "date") %>%
  complete(date = seq(min(date), max(date), by = "day"))

data_weekly <- clean_data(data_pm %>% slice(seq(1, n(), by = 7)))


# only the numeric values: drop date and 'X' columns
df <- data.matrix(data_weekly %>% select(-c(1:2)))

# store parameters in tibble
parameters <- tibble(
  n = nrow(df),
  p = ncol(df),
  khat = estimate_changepoint(df),
  khat_date = format(data_weekly[khat, ]$date, "%B %d, %Y")
)



# visually selecting m1 and m2
# future::plan(future::multisession, workers = 10)
# del_f1 <- get_delta_f1(df, parameters$khat, parallelize = TRUE)
# plot(del_f1, xlab = "m", ylab = "ΔF1", type = "l")

# del_f2 <- get_delta_f2(df, parameters$khat, parallelize = TRUE)
# plot(del_f2, xlab = "m", ylab = "ΔF2", type = "l")


parameters$m1 <- 15
parameters$m2 <- 18


print(parameters)




# running the tests (short), use cp_test, since we chose m's ourselves
m <- unlist(parameters |> select(m1, m2))

test_n <- cp_test_normalized(df, m = m, parallelize = TRUE)
test_s <- cp_test_sparsity_adj(df, m = m, parallelize = TRUE)


print("95% Confidence intervals for normalized l2-norm")
print(sqrt(ci_sq_norm(test_n)))
print(sqrt(ci_sq_norm_twoside(test_n)))


print("95% Confidence intervals for sparsity adjusted l2-norm")
print(sqrt(ci_sq_norm(test_s)))
print(sqrt(ci_sq_norm_twoside(test_s)))
