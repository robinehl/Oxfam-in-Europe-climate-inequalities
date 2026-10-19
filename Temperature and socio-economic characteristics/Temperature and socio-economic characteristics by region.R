### Data sources ###
## Income NUTS-2 Eurostat: https://ec.europa.eu/eurostat/databrowser/view/ilc_mdsd18/default/table?category=livcon.ilc.ilc_md.ilc_mdsd
## NUTS-2 boundaries: https://ec.europa.eu/eurostat/web/gisco/geodata/statistical-units/territorial-units-statistics
## Temperature observations (ERA5 post-processed daily statistics on single levels from 1940 to present) https://cds.climate.copernicus.eu/datasets/derived-era5-single-levels-daily-statistics


library(sf)
library(terra)
library(readxl)
library(ggplot2)

## Import and prepare data
nuts <- st_read("NUTS-2 boundaries.gpkg")
deprivation <- read_excel("Income NUTS-2 Eurostat.xlsx", sheet = "Severe deprivation")
revenue <- read_excel("Income NUTS-2 Eurostat.xlsx", sheet = "Disposable income")
tmax <- rast("Temperature observations.nc")

# Filter NUTS shapes to keep only NUTS-2 information
nuts2 <- nuts[nuts$LEVL_CODE == 2, ]
nuts2 <- st_transform(nuts2, "EPSG:4326")
nuts2$idx <- seq_len(nrow(nuts2))

# Attribute data on deprivation to the NUTS-2 regions
names(deprivation) <- c("NUTS_ID", "label", "privation")
names(revenue)      <- c("NUTS_ID", "label", "revenue")
deprivation$NUTS_ID <- trimws(deprivation$NUTS_ID)
revenue$NUTS_ID      <- trimws(revenue$NUTS_ID)
revenue$revenue <- suppressWarnings(as.numeric(revenue$revenue))
revenue[!is.na(revenue$NUTS_ID) & is.na(revenue$revenue), ]
# For certain regions, data is coded using the NUTS 2021 schema instead of the NUTS 2024 schema. These attributions are updated when the identifiant changed.
recoded <- c(
 "NL31" = "NL35",   # Utrecht
 "NL33" = "NL36",   # Zuid-Holland
 "NL42" = "NL23",   # Flevoland
 "SI01" = "SI03",   # Vzhodna Slovenija
 "SI02" = "SI04",   # Zahodna Slovenija
 "PT18" = "PT1C"    # Alentejo
)

deprivation$NUTS_ID_2024 <- ifelse(
 deprivation$NUTS_ID %in% names(recoded),
 unname(recoded[deprivation$NUTS_ID]), deprivation$NUTS_ID)
revenue$NUTS_ID_2024 <- ifelse(
 revenue$NUTS_ID %in% names(recoded),
 unname(recoded[revenue$NUTS_ID]), revenue$NUTS_ID)

# However, some regions are in 2021 shapes without having a direct successor. These regions are excluded from further analysis, since socio-economic data is thus unavailable from Eurostat sources
without_successor <- c("HU10", "PT16", "PT17", "HR04",
                    "NO01", "NO03", "NO04", "NO05")

sum(deprivation$NUTS_ID %in% without_successor)
sum(revenue$NUTS_ID %in% without_successor)

deprivation <- deprivation[!deprivation$NUTS_ID %in% without_successor, ]
revenue <- revenue[!grepl("ZZ$", revenue$NUTS_ID) &
                      !revenue$NUTS_ID %in% without_successor, ]
setdiff(deprivation$NUTS_ID_2024, nuts2$NUTS_ID)
setdiff(revenue$NUTS_ID_2024, nuts2$NUTS_ID)

# Attribution of the maximum temperature to every region. Therefore, rasterize maximum temperature data
zones <- rasterize(vect(nuts2), tmax_eu, field = "idx")
zres  <- zonal(tmax_eu, zones, fun = "mean", na.rm = TRUE)

# Identify regions excluded from rasterization. This can be the case in two scenarios: region too small or beyond range.
missing <- setdiff(nuts2$idx, zres[[1]])
nuts2[missing, c("NUTS_ID", "NAME_FREN")]

# For regions two small for the raster, attribution of the next known max temp value (using the center of the NUTS-2 region
if (length(missing) > 0) {
  pts       <- st_centroid(nuts2[missing, ])
  vpts      <- extract(tmax_eu, vect(pts))
  daily_sup <- as.data.frame(vpts[, -1])
  daily_sup$idx <- missing
} else {
  daily_sup <- NULL
}

# Calculate daily values for every region during the three summer months of 2026
daily <- as.data.frame(zres[, -1])
daily$idx <- zres[[1]]
daily <- rbind(daily, daily_sup)

# Using the daily values, calculate maximum temperature for the summer of 2026 per region and the medium temperature
day <- grep("^t2m_", names(daily))
m <- as.matrix(daily[, day])
tab_temp <- data.frame(
  NUTS_ID         = nuts2$NUTS_ID[daily$idx],
  region          = nuts2$NAME_FREN[daily$idx],
  temp_max_summer    = round(apply(m, 1, max, na.rm = TRUE) - 273.15, 1),
  temp_median_summer = round(apply(m, 1, median, na.rm = TRUE) - 273.15, 1)
)
stopifnot(nrow(tab_temp) == nrow(nuts2))

# Some regions still do not have any temperature values assigned. It is assumed that these regions are beyond the range of the analysis, eg oversea territories
tab_temp[!is.finite(tab_temp$temp_max_summer), c("NUTS_ID", "region")]

## Crossing information on temperature and socio-economic indicators for every NUTS-2 region included in the analysis

tab_final <- merge(tab_temp,
 deprivation[, c("NUTS_ID_2024", "privation")],
 by.x = "NUTS_ID", by.y = "NUTS_ID_2024", all.x = TRUE)
tab_final <- merge(tab_final,
 revenue[, c("NUTS_ID_2024", "revenue")],
 by.x = "NUTS_ID", by.y = "NUTS_ID_2024", all.x = TRUE)

tab_final <- tab_final[order(-tab_final$temp_max_summer),
 c("NUTS_ID", "region", "temp_max_summer",
 "temp_median_summer", "privation", "revenue")]
rownames(tab_final) <- NULL

head(tab_final, 15)
summary(tab_final[, c("temp_max_summer", "temp_median_summer", "privation", "revenue")])
nrow(tab_final)

# Commented line to write the results in a table
## write.csv2(tab_final, "Temperature deprivation summer 2026 analysis.csv", row.names = FALSE)

## Regression analysis of temperature and socio-economic characteristics

# Limit the regression analysis to NUTS-2 regions with temperature information (and continue working with “dat” variable”
dat <- tab_final[!is.na(tab_final$privation) &
 is.finite(tab_final$temp_max_summer), ]
nrow(dat)

# Regression analysis share of deprivation and maximum summer temperature
modele <- lm(privation ~ temp_max_summer, data = dat)
summary(modele)
confint(modele)

# Regression analysis share of deprivation and median summer temperature
summary(lm(privation ~ temp_median_summer, data = dat))

## Create scatter plots

ggplot(dat, aes(x = temp_max_summer, y = privation)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "firebrick") +
  labs(title = "Severe social and material privation against maximum temperature during summer 2026",
       subtitle = paste0("European NUTS-2 regions, n = ", nrow(dat)),
       x = "Maximum temperature during summer 2026 (°C)",
       y = "Share of severe social and material deprivation (%)") +
  theme_minimal()