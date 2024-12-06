library(gstat)
library(sf)
library(terra)
library(ggplot2)
library(ggrepel)
library(ggspatial)
library(tidyterra)

# Boundaries of the soil samples
site_bound <- vect("data/shp/site_bound.shp")

# grid for interpolation
xrange <- seq(ext(site_bound)[1], ext(site_bound)[2], 1)
yrange <- seq(ext(site_bound)[3], ext(site_bound)[4], 1)
zrange <- seq(10,90,20)

site_grid <- expand.grid(xrange, yrange, zrange)
colnames(site_grid) <- c("x", "y", "z")
site_grid <- st_as_sf(site_grid, coords=c("x", "y", "z"), crs=st_crs(site_bound))

# Raw soil measurements
raw_data <- read.csv("data/caldeirao_data.csv")

data_points <- st_as_sf(raw_data, coords=c("x", "y", "Depth"), crs=st_crs(site_bound))

# Select only relevant columns
data <- raw_data[,c(7:20)]
names <- names(data)

# Some columns are integer
data <- as.data.frame(sapply(data, as.numeric))

# Define pretic areas (criteria from World Reference Base for Soil Resources)
pretic <- function(x) {x$C >= 6 & (x$Ca + x$Mg) >= 1 & x$P >= 100}

# Symbols on plot will be based on pretic x non-pretic, depth and ceramic weight
groups <- as.factor(as.numeric(pretic(raw_data)))
depths <- as.factor(raw_data$Depth)
sizes <- raw_data$PotteryWeight * 100

# For reproducibility
set.seed(100)

# PCA
res.pca <- prcomp(data, center=TRUE, scale=TRUE)
print(summary(res.pca))

fig2 <- ggplot(res.pca$x, aes(x=PC1, y=PC2)) +
        geom_jitter(aes(color=groups, shape=depths, size=sizes)) +
        geom_segment(data=res.pca$rotation, aes(x=0,y=0,xend=(PC1*10), yend=(PC2*10)), arrow=arrow(length=unit(1/2, "picas"))) +
        annotate("text", x=(res.pca$rotation[,1]*12), y=(res.pca$rotation[,2]*12), label=rownames(res.pca$rotation)) +
        scale_shape_discrete(name="Depth (cm)", labels=c("0-20", "20-40", "40-60", "60-80", "80-100")) +
        scale_color_discrete(name="Pretic horizon", labels=c("No", "Yes")) + 
        scale_size(name="Pottery weight (g)") +
        xlab("Dimension 1") +
        ylab("Dimension 2") +
        theme_minimal()

png("figs/Figure2.png", res=300, width=2000, height=1600)
plot(fig2)
dev.off()

# Store Dimension 1 on data frame for interpolation
data_points$PC1 <- res.pca$x[,1]

# Kriging
vgm <- variogram(PC1~1, data_points)
plot(vgm)

# Range ~ 200, sill ~ 8, nugget ~ 2
fit <- fit.variogram(vgm, vgm(8, c("Exp", "Sph", "Gau"), 200, 2))
plot(vgm, fit)

# 3D Kriging
krig <- krige(PC1~1, data_points, site_grid, fit)

krig_stack <- c()

# Extract raster for each depth
for (i in seq(10,90,20))
{
    idx <- which(st_coordinates(krig)[,3] == i)
    sel <- krig[idx,]
    krig_sp <- vect(sel)
    raster_template <- rast(ext(krig_sp), res=1)
    krig_stack <- c(krig_stack, rasterize(krig_sp, raster_template, "var1.pred"))
}

krig_stack <- rast(krig_stack)
crs(krig_stack) <- crs(site_bound)
krig_stack <- mask(krig_stack, site_bound)
names(krig_stack) <- c("0-20cm", "20-40cm", "40-60cm", "60-80cm", "80-100cm")

data_points$lyr <- factor(raw_data$Depth, levels=seq(10,90,20), labels=names(krig_stack))
data_points$PotteryWeight[data_points$PotteryWeight == 0] <- NA
dates <- raw_data[raw_data$C14Age != "",]
dates$lyr <- factor(dates$Depth, levels=seq(10,90,20), labels=names(krig_stack))

zmin <- min(values(krig_stack), na.rm=T)
zmax <- max(values(krig_stack), na.rm=T)

fig3 <- ggplot() + geom_spatraster(data=krig_stack) +
    scale_fill_viridis_c(option="turbo", limits=c(zmin,zmax), na.value="transparent", name="Dimension 1", direction=-1) +
    geom_sf(data=data_points, aes(size=PotteryWeight*100), pch=1) +
    scale_size(name="Pottery weight (g)") +
    geom_point(data=dates, aes(x=x, y=y, shape="Dated sample")) +
    geom_text_repel(data=dates, aes(x=x, y=y, label=C14Age), size=2.5, box.padding=0.5) +
    scale_shape_manual(name="", values=c(3)) +
    geom_sf(data=site_bound, fill=NA) +
    annotation_scale(location="bl", width_hint=0.4, style="ticks") +
    annotation_north_arrow(location="bl", which_north="true", style=north_arrow_minimal,
                           height=unit(1, "cm"), pad_y=unit(0.5, "cm"), pad_x=unit(0.0, "cm")) +
    facet_wrap(~lyr) +
    theme_void()

png("figs/Figure3.png", res=300, width=2000, height=1600)
plot(fig3)
dev.off()
