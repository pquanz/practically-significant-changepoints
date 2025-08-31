# install.packages("remotes")
# install_remotes::github("openaq/openaq-r@*release")

library(openaq)
library(dplyr)


set_api_key(Sys.getenv("OPENAQ_API_KEY"))



# fetch only PM10 and PM25
pollutant_ids <- c(1, 2)



fetch_all_locations <- function(bbox) {
  page <- 1
  pieces <- list()
  limit <- 1000

  repeat {
    loc_list <- tryCatch(
      tibble::as_tibble(
        list_locations(bbox = bbox, limit = limit, page = page)
      ),
      error = function(e) {
        warning(sprintf("Locations page %d failed: %s", page, e$message))
        NULL
      }
    )

    if (is.null(loc_list) || nrow(loc_list) == 0) break
    pieces[[length(pieces) + 1]] <- loc_list

    if (nrow(loc_list) < limit) break  # last page
    page <- page + 1
  }
  bind_rows(pieces)
}



# fully southwest europe (almost up to paris, not full france)
# Uncomment following lines for fetching locations:
locations <- fetch_all_locations(
  bbox = c(xmin = -9.75, ymin = 36.00, xmax = 5.91, ymax = 48.53)
)
write.csv(locations, "raw/locations.csv")

# otherwise read in fetched data
# locations <- tibble(read.csv("raw/locations.csv"))


message(sprintf("Found a total of %s locations", length(locations$id)))
message("Getting sensor data of all locations.")

all_sensors <- tibble()
for (i in seq_along(locations$id)) {
  if (i %% 50 == 0) {
    message(sprintf("Progress: %-3d / %s", i, length(locations$id)))
  }
  sensors <- tibble(list_location_sensors(locations_id = locations$id[i])) %>%
    filter(parameters_id %in% pollutant_ids)
  all_sensors <- bind_rows(all_sensors, sensors)

  Sys.sleep(runif(1, min = 0.5, max = 1.5))
}

write.csv(all_sensors, "raw/sensors.csv")
