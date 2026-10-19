### Data sources ###
## Hospital location: https://ec.europa.eu/eurostat/web/gisco/geodata/basic-services#Healthcare
##  Potential flood prone area extent: https://sdi.eea.europa.eu/data/28c36420-c31b-440e-80c5-8064696f3517

library(sf)
library(dplyr)

## Import data
hosp <- st_read("Hospital location.geojson", quiet = TRUE)
st_layers(dsn = "Potential_flood-prone_area_extent.gdb") # List available layers of the data on flood extent
flood <- st_read(
 dsn   = "Potential_flood-prone_area_extent.gdb",
 layer = "Potential_flood_prone_area_extent_wgs84",
 quiet = TRUE
)

## Prepare spatial analysis by correcting corrupted data or, where unavoidable, removing corrupted data
hosp <- hosp[!st_is_empty(hosp), ] # Remove two hospitals without information on their location
# Correct data on flood-prone areas which has non-closed polygones
flood <- st_make_valid(flood)

hosp <- st_transform(hosp, st_crs(flood)) # Uniform projection of the geographic data

Area_Analyzed <- st_as_sfc(st_bbox(flood)) # Exclude hospitals from the analysis that are located beyond Europe, as information on the flood-prone areas is only available for Europe
hosp$hors_Area_Analyzed <-
  lengths(st_intersects(hosp, Area_Analyzed)) == 0
table(hosp$hors_Area_Analyzed)

## Spatial analysis to calculate the share of hospitals in flood-prone areas. Therefore, check whether the location of the hospital is located within a food-prone area. Note that every hospital is only represented by one data point to indicate a specific location, even though the area of the hospital may be larger
hosp$inond <- lengths(st_intersects(hosp, flood)) > 0

# Calculate global share of hospitals in flood-prone areas
mean(hosp$inond) * 100
mean(hosp$inond[!hosp$hors_Area_Analyzed]) * 100 # Correct the mean by taking into account the fact that some hospitals are not in Europe

# 7f. Résultats par pays ---------------------------------------
res_pays <- hosp |>
  st_drop_geometry() |>
  filter(!hors_Area_Analyzed) |>
  group_by(cntr_id) |>
  summarise(
    n_hopitaux     = n(),
    n_inond        = sum(inond),
    part_inond_pct = round(100 * n_inond / n_hopitaux, 1)
  ) |>
  arrange(desc(part_inond_pct))
res_pays

# 7g. Aperçu cartographique du résultat -------------------------
plot(st_geometry(flood), col = "lightblue", border = NA)
plot(st_geometry(hosp[hosp$inond, ]),
     pch = 16, cex = 0.6, col = "red", add = TRUE)