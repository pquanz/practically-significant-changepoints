# install.packages("remotes")
# install_remotes::github("openaq/openaq-r@*release")

# install.packages("progressr")


library(openaq)
library(dplyr)
library(tidyr)
library(progressr)


set_api_key(Sys.getenv("OPENAQ_API_KEY"))

# store separate by PM10 and PM25
pollutant_ids <- c(1, 2)

all_sensor_ids <- read.csv("raw/sensors.csv")$id
message(sprintf("Found a total of %s sensors", length(all_sensor_ids)))



# Fetching sensor data for single id
fetch_all_sensor_measurements <- function(id) {
  page <- 1
  out_list <- list()

  repeat {
    df <- list_sensor_measurements(
      id, data = "days", limit = 1000, page = page
    )
    df <- tibble::as_tibble(df)

    Sys.sleep(runif(1, min = 0.5, max = 1.5))
    if (nrow(df) == 0) break
    out_list[[page]] <- df

    # if fewer than 1000 rows, this was the last page
    if (nrow(df) < 1000) break
    page <- page + 1
  }

  bind_rows(out_list)
}



# store measurements in tibble
measures_df <- tibble()


message("Starting download of sensor measurement data.")
handlers("txtprogressbar")
with_progress({   # progress bar
  p <- progressor(steps = length(all_sensor_ids))

  for (i in seq_along(all_sensor_ids)) {
    p(message = sprintf("i=%d", i))
    id <- allall_sensor_ids[i]
    tryCatch({
      res <- fetch_all_sensor_measurements(id) %>%
        transmute(
          value = value,
          pollutant = parameter_name,
          date = as.Date(datetime_from),
          sensor_id = id
        )
      measures_df <- bind_rows(measures_df, res)
    },
    error = function(e) {
      warning(sprintf("Data of sensor %s not stored: (%s)", id, e$message))
    })
  }
})


message(
  "Download of sensor data complete. Now splitting and storing by pollutant."
)


for (poll in unique(measures_df$pollutant)) {
  poll_wide <- measures_df %>%
    filter(pollutant == poll) %>%     # filter by pollutant
    mutate(date = as.Date(date)) %>%  # correctly format dates
    group_by(date, sensor_id) %>%     # remove duplicates and sort by date
    summarise(
      value = mean(value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(date, sensor_id) %>%
    pivot_wider(
      names_from = sensor_id,
      values_from = value,
      names_prefix = "id_"
    ) %>%               # stack columns next to each other (order by date)
    complete(
      date = seq(min(date), max(date), by = "day")
    ) %>%               # fill in missing dates and sort again
    arrange(date)

  # save data
  write.csv(poll_wide, paste0("raw/", poll, "_data", ".csv"))
}
