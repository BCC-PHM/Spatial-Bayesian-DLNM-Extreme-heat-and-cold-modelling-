setwd("C:/Users/TMPACGAG/OneDrive - Birmingham City Council/Documents/R projects/PHM/Extreme heat and cold")


library(tidyverse)
library(readr)
library(sf)

# ------------------------------------------------------------
# 1. Read LSOA boundary files
# ------------------------------------------------------------

LSOA_21_map = read_sf(
  "data/external/boundaries/boundaries-lsoa-2021-birmingham/boundaries-lsoa-2021-birmingham.shp"
)

LSOA_11_map = read_sf(
  "data/external/boundaries/LSOA_2011_map/Birmingham LSOA shape.shp"
)

# ------------------------------------------------------------
# 2. Reproject BOTH layers to British National Grid (EPSG:27700)
#    (required for valid area calculations)
# ------------------------------------------------------------

LSOA_11_map_bng = st_transform(LSOA_11_map, 27700)
LSOA_21_map_bng = st_transform(LSOA_21_map, 27700)

# ------------------------------------------------------------
# 3. Compute geometric intersections
#    (returns polygon fragments, lines, and points)
# ------------------------------------------------------------

intersections = st_intersection(
  LSOA_11_map_bng %>% select(lsoa11 = LSOA11CD),
  LSOA_21_map_bng %>% select(lsoa21 = LSOA21CD)
)

# ------------------------------------------------------------
# 4. Keep polygonal overlaps only (area-based intersections)
# ------------------------------------------------------------

intersections_poly = intersections %>%
  filter(st_dimension(geometry) == 2)

# ------------------------------------------------------------
# 5. Collapse geometry fragments to LSOA11–LSOA21 pairs
#    and remove negligible boundary slivers
# ------------------------------------------------------------

weights_pair = intersections_poly %>%
  mutate(area_intersection = as.numeric(st_area(geometry))) %>%
  filter(area_intersection > 1) %>%        # remove numerical artefacts
  st_drop_geometry() %>%
  group_by(lsoa11, lsoa21) %>%
  summarise(
    area_intersection = sum(area_intersection),
    .groups = "drop"
  )

# ------------------------------------------------------------
# 6. Compute total area of each LSOA11 (denominator)
# ------------------------------------------------------------

lsoa11_area = LSOA_11_map_bng %>%
  mutate(area_lsoa11 = as.numeric(st_area(geometry))) %>%
  st_drop_geometry() %>%
  select(lsoa11 = LSOA11CD, area_lsoa11)

# ------------------------------------------------------------
# 7. Calculate area-based weights
# ------------------------------------------------------------

weights_pair = weights_pair %>%
  left_join(lsoa11_area, by = "lsoa11") %>%
  mutate(weight = area_intersection / area_lsoa11)

# ------------------------------------------------------------
# 8. Renormalise weights within each LSOA11
#    (guards against floating-point error)
# ------------------------------------------------------------

weights_pair = weights_pair %>%
  group_by(lsoa11) %>%
  mutate(weight = weight / sum(weight)) %>%
  ungroup()

# ------------------------------------------------------------
# 9. Final validation (optional but recommended)
# ------------------------------------------------------------

# Each LSOA11 should sum to 1
# weights_pair %>%
#   group_by(lsoa11) %>%
#   summarise(w = sum(weight)) %>%
#   filter(abs(w - 1) > 1e-12)


write.csv(weights_pair, "data/processed/lsoa11_to_lsoa21_weights.csv")

