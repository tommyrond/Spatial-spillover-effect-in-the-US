## Libraries, dataset --------------------

library(tidyverse)
library(arrow)
library(sf)
library(tmap)
library(spdep)
library(spatialreg)

set.seed(1)

setwd(...) # If needed

# The dataset is read
df <- read_parquet("election_df.parquet")

sf_data <- df %>%
  mutate(geometry = st_as_sfc(geometry)) %>%
  st_as_sf()

sf_data <- st_set_crs(sf_data, 4269)

sf_data <- sf_data %>%
  mutate(Median_income = Median_income/1000)

summary(sf_data)

print(sf_data)

ggplot(data = sf_data) +
  geom_sf(aes(fill=Republican_votes)) +
  theme_minimal()

# The neighbours by contiguity of each county are identified
neighbours <- poly2nb(sf_data)
no_neighbours <- sapply(neighbours, function(x) x[1] == 0)
indexes_no_neighbors <- which(no_neighbours)
# Alaska and some other islands clearly have no neighbours
weights_matrix <- nb2mat(neighbours, style = "B", zero.policy = TRUE)
nei_listw <- nb2listw(neighbours,style="W",zero.policy = T)

## Elhorst's algorithm -------------------

# Linear model
OLSmodel <- lm(Republican_votes ~ Median_income + Non_employed + Low_educated + Immigrants_23 + Delta,
               data=sf_data)
summary(OLSmodel)

# Lagrange multiplier test
natOLSlmTests <- lm.RStests(OLSmodel, nei_listw, 
                            test=c("LMerr", "LMlag", "RLMerr", "RLMlag"))
summary(natOLSlmTests)

# Spatial Durbin Model
SDM <- lagsarlm(Republican_votes ~ Median_income + Non_employed + Low_educated + Immigrants_23 + Delta,
                data=sf_data, Durbin=T, listw = nei_listw, zero.policy=T)
summary(SDM)

# Spatial Autoregressive Model
SAR <- lagsarlm(Republican_votes ~ Median_income + Non_employed + Low_educated + Immigrants_23 + Delta,
                data=sf_data, listw = nei_listw, zero.policy=T)
summary(SAR)

# Spatial Error Model
SEM <- errorsarlm(Republican_votes ~ Median_income + Non_employed + Low_educated + Immigrants_23 + Delta,
                  data=sf_data, listw = nei_listw, zero.policy=T)
summary(SEM)

# Likelihood-ratio test on the restrictions of the SDM parameters
anova(SDM, SAR)
anova(SDM, SEM)

## Impact measures ----------------------

impSDM <- impacts(SDM, listw=nei_listw, R=100)
summary(impSDM, zstats=TRUE, short=TRUE)

## Spillover effects -----------------------

# Global spillover effect
moran.test(sf_data$Republican_votes, nei_listw, randomisation=F)
moran.test(sf_data$Republican_votes, nei_listw, randomisation=T)
moran.mc(sf_data$Republican_votes, nei_listw, nsim=999)

# Moran's I test in the OLS residuals
lm.morantest(OLSmodel,nei_listw,resfun=rstudent)

# Moran scatterplot
mplot <- moran.plot(sf_data$Republican_votes, nei_listw, main="Moran scatterplot")

# Local Moran's I index
lmI <- localmoran(sf_data$Republican_votes, nei_listw)
head(lmI)
sf_data$lmI <- lmI[,1]

# Choropleth map of the influential counties
sf_data$locmpv <- p.adjust(lmI[, "Pr(z != E(Ii))"], "bonferroni")

sf_data$quadrant <- attr(lmI, "quadr")$mean
sf_data$quadrant[sf_data$locmpv>0.05] <- NA
tm_shape(sf_data, bbox=c(xmin = -130, ymin = 23, xmax = -67, ymax = 50)) + tm_polygons("quadrant") +
  tm_layout(legend.outside = TRUE)

# Choropleth map of the Local Moran's I value
centroids <- st_centroid(sf_data)
coords <- st_coordinates(centroids)

# Remove non interesting counties
filtered_sf_data <- sf_data[coords[,1] > -140,]

tm_shape(filtered_sf_data) + 
  tm_polygons("lmI", title = "Local Moran's I values") + 
  tm_layout(legend.outside = TRUE)
