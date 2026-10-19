### Data sources ###
## Location of primary schools: https://database.espon.eu/indicator/885/
## Temperature observations (ERA5 post-processed daily statistics on single levels from 1940 to present): https://cds.climate.copernicus.eu/datasets/derived-era5-single-levels-daily-statistics

library(sf)
library(terra)
library(openxlsx)

## Import and prepare data
schools <- st_read("Schools/POI_SchoolsPrimary_2021.shp")
tmax <- rast("Temperature extremes Europe june 2026.nc")

schools_ll <- st_transform(schools, "EPSG:4326")

## Attribute observed temperature to schools
# For every primary school, extract the maximum temperature observed during June 2026
val <- extract(tmax, vect(schools_ll))
# Verify temperature attribution
stopifnot(nrow(val) == nrow(schools_ll))
dim(val)        # Expected return: 169672  31
head(val)       # Expected column headers: t2m_1 ... t2m_30

mat <- as.matrix(val[, grep("^t2m_", names(val))])
dim(mat)        # doit afficher : 169672  30

## Identify schools having encountered temperatures >30°C
Days     <- rowSums(mat > 303.15, na.rm = TRUE)
# Exclude schools beyond range
BeyondRange <- rowSums(!is.na(mat)) == 0
schools_ll$HotDays <- NA
schools_ll$HotDays[val$ID] <- Days
schools_ll$exposed <- NA
schools_ll$exposed[val$ID] <- ifelse(BeyondRange, NA, Days > 0)

# Count, for every country, the number of schools having encountered temperatures >30°C
sum(schools_ll$exposed, na.rm = TRUE)
sum(BeyondRange)                                    # Sum of the schools beyond range
sort(table(schools_ll$cc[schools_ll$exposed == TRUE]), decreasing = TRUE)

## Data visualization
tmp <- data.frame(country = schools_ll$cc,
 exposed = schools_ll$exposed == TRUE)
agg <- aggregate(exposed ~ country, data = tmp,
 FUN = function(x) c(total = length(x), SchoolsExposed = sum(x)))
tab_country <- do.call(data.frame, agg)
tab_country$part_pct <- round(100 * tab_country$exposed.SchoolsExposed /
 tab_country$exposed.total, 1)
tab_country <- tab_country[order(-tab_country$exposed.SchoolsExposed), ]
names(tab_country) <- c("Country", "Number of schools", "Schools exposed to temperatures > 30 °C", "Schools exposed to temperatures > 30 °C")
tab_country
write.xlsx(tab_country, 'Schools temperature exposure.xlsx')