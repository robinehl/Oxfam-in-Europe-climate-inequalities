### Data sources
### Income NUTS-2 Eurostat: https://ec.europa.eu/eurostat/databrowser/view/ilc_mdsd18/default/table?category=livcon.ilc.ilc_md.ilc_mdsd
### NUTS-2 boundaries: https://ec.europa.eu/eurostat/web/gisco/geodata/statistical-units/territorial-units-statistics
### Temperature observations (ERA5 post-processed daily statistics on single levels from 1940 to present): https://cds.climate.copernicus.eu/datasets/derived-era5-single-levels-daily-statistics


library(sf)
library(terra)
library(readxl)
library(ggplot2)

## Import and prepare data

nuts <- st_read("NUTS-2 boundaries.gpkg")
deprivation <- read_excel("Income NUTS-2 Eurostat.xlsx", sheet = "Severe deprivation")
revenue     <- read_excel("Income NUTS-2 Eurostat.xlsx", sheet = "Disposable income")

# Temperature data are obtained in three seperate files
fichiers <- file.path("2026 Maximum temperature", c(
 "2026 Maximum Temperature June.nc",
 "2026 Maximum Temperature July.nc",
 "2026 Maximum Temperature August.nc"))
stopifnot(all(file.exists(fichiers)))

tmax <- rast(fichiers)
print(tmax)
tps <- time(tmax)
stopifnot(nlyr(tmax) == 92)

# Prepare the geographic information of the temperature data. Stored with between 0 to 180, they are centered on Greenwich with a range of -180 to 180
tmax <- rotate(tmax)

# Impose chronological order of temperature files
tmax <- tmax[[order(tps)]]

# Rename date Indicator of temperature files
names(tmax) <- format(time(tmax), "%Y-%m-%d")

# Temperature conversion from °K to °C
tmax <- tmax - 273.15

# Use same geographic projection of the temperature data than for the shapes of the NUTS-2 regions
nuts2 <- nuts[nuts$LEVL_CODE == 2, ]
nuts2 <- st_transform(nuts2, terra::crs(tmax))

# Harmonize column names of the NUTS-2 shapes and the regional socio-economic indicators
names(deprivation) <- c("NUTS_ID", "label", "privation")
names(revenue)      <- c("NUTS_ID", "label", "revenue")
deprivation$NUTS_ID <- trimws(deprivation$NUTS_ID)
revenue$NUTS_ID     <- trimws(revenue$NUTS_ID)
revenue$revenue <- suppressWarnings(as.numeric(revenue$revenue))
revenue[!is.na(revenue$NUTS_ID) & is.na(revenue$revenue), ]

# Reassign some regions whose code changed because of the 2024 NUTS labelling compared to the previous version from 2021. Some regions remain unassigned
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

without_successor <- c("HU10", "PT16", "PT17", "HR04",
                       "NO01", "NO03", "NO04", "NO05")
deprivation <- deprivation[!deprivation$NUTS_ID %in% without_successor, ]
revenue <- revenue[!grepl("ZZ$", revenue$NUTS_ID) &
                       !revenue$NUTS_ID %in% without_successor, ]

# Since the data from Eurostat stores information on regions modified because of NUTS-2 reforms in two seperate lines, keep only the line which should be the most recent (without a mention to an older NUTS version)
dedupe_2024 <- function(df, valeur) {
    df <- df[order(df$NUTS_ID != df$NUTS_ID_2024, is.na(df[[valeur]])), ]
    df[!duplicated(df$NUTS_ID_2024), ]
}
deprivation <- dedupe_2024(deprivation, "privation")
revenue     <- dedupe_2024(revenue, "revenue")

## For every region and for every summer day, calculate the median of the maximum temperatures observed in a region "temp_median_daily"
mat_median <- exact_extract(tmax, nuts2, "median", stack_apply = TRUE)
rownames(mat_median) <- nuts2$NUTS_ID

# Since no temperature value is attributed for some regions, the next possible value is allocated using the point at the center of the region
Without_Value <- rowSums(is.na(mat_median) | mat_median <= 0,
                       na.rm = TRUE) == nlyr(tmax)
if (any(Without_Value)) {
    pts <- st_centroid(nuts2[Without_Value, ])
    v   <- terra::extract(tmax, vect(pts))
    mat_median[Without_Value, ] <- as.matrix(v[, -1])
    print(nuts2$NUTS_ID[Without_Value])
}

# For every region, calculate the median of the daily median value of the maximum températures
tab_temp <- data.frame(
    NUTS_ID = nuts2$NUTS_ID,
    region  = nuts2$NAME_LATN,
    temp_median_summer = round(apply(mat_median, 1, median, na.rm = TRUE), 1)
)

# Exclude regions still without value
tab_temp[!is.finite(tab_temp$temp_median_summer), c("NUTS_ID", "region")]

## Temperature indicature is calculated for every region. Now, cross these values with the socio-economic indicators
tab_final <- merge(tab_temp,
                   deprivation[, c("NUTS_ID_2024", "privation")],
                   by.x = "NUTS_ID", by.y = "NUTS_ID_2024", all.x = TRUE)
tab_final <- merge(tab_final,
                   revenue[, c("NUTS_ID_2024", "revenue")],
                   by.x = "NUTS_ID", by.y = "NUTS_ID_2024", all.x = TRUE)

# Limit analysis to EU-27

pays_ue27 <- c("AT", "BE", "BG", "HR", "CY", "CZ", "DK", "EE", "FI",
               "FR", "DE", "EL", "HU", "IE", "IT", "LV", "LT", "LU",
               "MT", "NL", "PL", "PT", "RO", "SK", "SI", "ES", "SE")
tab_final <- tab_final[substr(tab_final$NUTS_ID, 1, 2) %in% pays_ue27, ]

tab_final <- tab_final[order(-tab_final$temp_median_summer),
  c("NUTS_ID", "region", "temp_median_summer", "privation", "revenue")]
rownames(tab_final) <- NULL

head(tab_final, 15)
summary(tab_final[, c("temp_median_summer", "privation", "revenue")])
nrow(tab_final)

#Linear regression analysis
dat <- tab_final[!is.na(tab_final$privation) &
                     !is.na(tab_final$temp_median_summer) &
                     tab_final$temp_median_summer > 0, ]
nrow(dat)

modele <- lm(privation ~ temp_median_summer, data = dat)
summary(modele)
confint(modele)

modele <- lm(revenue ~ temp_median_summer, data = dat)
summary(modele)
confint(modele)


# Correlation test
cor.test(dat$privation, dat$temp_median_summer, method = "spearman")

# Scatter plot for temperature & severe social and economic deprivation
ggplot(dat, aes(x = temp_median_summer, y = privation)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "firebrick") +
    labs(title = "Severe social and material deprivation against daily maximum temperature during summer 2026",
         subtitle = paste0("EU NUTS-2 regions, n = ", nrow(dat)),
         x = "Heat exposure indicator (regional median of median daily maximum temperature in °C)",
         y = "Severe social and economic deprivation (%)") +
    theme_minimal()

# Scatter plot for temperature & disposable income
ggplot(dat, aes(x = temp_median_summer, y = revenue)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "firebrick") +
    labs(title = "Income against daily maximum temperature during summer 2026",
         subtitle = paste0("EU NUTS-2 regions, n = ", nrow(dat)),
         x = "Heat exposure indicator",
         y = "Annual disposable income (€)") +
    theme_minimal()