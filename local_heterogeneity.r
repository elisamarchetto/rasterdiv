library(raster)
library(rasterdiv)
library(RStoolbox)
library(ggplot2)
library(scico)

# set working directory
# sentinel-2, instrument:MSI,31/07/2020
# image elabrated with SNAP: bands 2,3,4,8
# majella <- stack("subset_0_of_S2A_MSIL1C_20200731T095041_N0209_R079_T33TVG_20200731T102505.tif")

# RGB images
ggRGB(majella, r = 4, g = 2, b = 3,stretch = "lin") + ggtitle("False colors(NIR,GREEN,RED), Sentinel-2 image (10m)") + theme_light() + theme(plot.title = element_text(size = rel(0.9)))
ggRGB(majella, r = 3, g = 2, b = 1,stretch = "lin") + ggtitle("Natural colors(RED,GREEN,BLUE), Sentinel-2 image (10m)") + theme_light() + theme(plot.title = element_text(size = rel(0.9)))

# First crop on the original sentinel image
extent <- c(413685,441135,4631415,4686315)
majella_crop <- crop(majella, extent)


# NDVI computation, NDVI rescaled to 8-bit
ndvi_mj <- spectralIndices(majella_crop, red = "subset_0_of_S2A_MSIL1C_20200731T095041_N0209_R079_T33TVG_20200731T102505.3", nir = "subset_0_of_S2A_MSIL1C_20200731T095041_N0209_R079_T33TVG_20200731T102505.4", indices = "NDVI")
ndvi_mjr <- aggregate(ndvi_mj, fact=2)
ndvi_mjrs <- stretch(ndvi_mjr, minv=0, maxv=255)
storage.mode(ndvi_mjrs[]) = "integer"

# RGB images
ggRGB(majella_crop, r = 4, g = 2, b = 3,stretch = "lin") + ggtitle("False colors Sentinel-2 image (10m)") + theme_light() + theme(plot.title = element_text(size = rel(0.9))) 
ggRGB(majella_crop, r = 3, g = 2, b = 1,stretch = "lin") + ggtitle("Natural colors Sentinel-2 image (10m)") + theme_light() + theme(plot.title = element_text(size = rel(0.9))) 

# NDVI visualisation
ggR(ndvi_mjrs, geom_raster = TRUE) + ggtitle("NDVI") + theme_light() + theme(plot.title.position ='plot', plot.title = element_text(hjust = 0.5, size = rel(1))) 

# computation indices first crop with rasterdiv functions
sha <- Shannon(ndvi_mjrs, window=9, rasterOut=TRUE, np=3,na.tolerance=0.9, cluster.type="SOCK", debugging=FALSE)
bepar <- BergerParker(ndvi_mjrs, window=9, rasterOut=TRUE, np=3, na.tolerance=0.9, cluster.type="SOCK", debugging=FALSE)
pielou <- Pielou(ndvi_mjrs, window=9, rasterOut=TRUE, np=3, na.tolerance=0.9, cluster.type="SOCK", debugging=FALSE)
renyi5 <- Renyi(ndvi_mjrs, window=9, alpha=5, base=exp(1), rasterOut=TRUE, np=3, na.tolerance=0.9, cluster.type="SOCK", debugging=FALSE)
renyi50 <- Renyi(ndvi_mjrs, window=9, alpha=50, base=exp(1), rasterOut=TRUE, np=3, na.tolerance=0.9, cluster.type="SOCK", debugging=FALSE)
rao <- Rao(ndvi_mjrs, dist_m="euclidean", window=9, rasterOut = TRUE, mode="classic",lambda=0, shannon=FALSE, rescale=FALSE, na.tolerance=0.9, simplify=3, np=3, cluster.type="SOCK", debugging=FALSE)
parao3 <- paRao(ndvi_mjrs, dist_m="euclidean", window=9, alpha=3, method="classic", rasterOut=TRUE, lambda=0, na.tolerance=0.9, rescale=FALSE, diag=TRUE,simplify=1, np=3,cluster.type="SOCK", debugging=FALSE)
parao20 <- paRao(ndvi_mjrs, dist_m="euclidean", window=9, alpha=20, method="classic", rasterOut=TRUE, lambda=0, na.tolerance=0.9, rescale=FALSE, diag=TRUE,simplify=1, np=3,cluster.type="SOCK", debugging=FALSE)
parao100 <- paRao(ndvi_mjrs, dist_m="euclidean", window=9, alpha=100, method="classic", rasterOut=TRUE, lambda=0, na.tolerance=0.9, rescale=FALSE, diag=TRUE,simplify=1, np=3,cluster.type="SOCK", debugging=FALSE)

# indices values visualisation 
ggR(sha, geom_raster = TRUE) + scale_fill_scico(palette = "hawaii") + ggtitle("Shannon's entropy") + theme_light() + theme(plot.title.position ='plot', plot.title = element_text(hjust = 0.5))
ggR(bepar, geom_raster = TRUE) + scale_fill_scico(palette = "hawaii") + ggtitle("Berger-Parker index") + theme_light() + theme(plot.title.position ='plot', plot.title = element_text(hjust = 0.5))  
ggR(pielou, geom_raster = TRUE) + scale_fill_scico(palette = "hawaii") + ggtitle("Pielou's evenness index") + theme_light() + theme(plot.title.position ='plot', plot.title = element_text(hjust = 0.5)) 
ggR(renyi5, geom_raster = TRUE) +scale_fill_scico(palette = "hawaii") + ggtitle("Rényi's index alpha=5") + theme_light() + theme(plot.title.position ='plot', plot.title = element_text(hjust = 0.5))  
ggR(renyi50, geom_raster = TRUE) + scale_fill_scico(palette = "hawaii") + ggtitle("Rényi's index alpha=50") + theme_light() + theme(plot.title.position ='plot', plot.title = element_text(hjust = 0.5))  
ggR(rao, geom_raster = TRUE) + scale_fill_scico(palette = "hawaii") + ggtitle("Rao's Q heterogeneity index") + theme_light() + theme(plot.title.position ='plot', plot.title = element_text(hjust = 0.5))  
ggR(parao3, geom_raster = TRUE) + scale_fill_scico(palette = "hawaii") + ggtitle("Parametric Rao's Q alpha=3") + theme_light() + theme(plot.title.position ='plot', plot.title = element_text(hjust = 0.5))  
ggR(parao20, geom_raster = TRUE) + scale_fill_scico(palette = "hawaii") + ggtitle("Parametric Rao's Q alpha=20") + theme_light() + theme(plot.title.position ='plot', plot.title = element_text(hjust = 0.5))  
ggR(parao100, geom_raster = TRUE) + scale_fill_scico(palette = "hawaii") + ggtitle("Parametric Rao's Q alpha=100") + theme_light() + theme(plot.title.position ='plot', plot.title = element_text(hjust = 0.5))   

# second crop 
ext <- c(422000,434000,4655000,4676000)
m_amaro <- crop(majella_crop,ext)

# NDVI computation, NDVI rescaled to 8-bit
ndvi_m_a <- spectralIndices(m_amaro, red = "subset_0_of_S2A_MSIL1C_20200731T095041_N0209_R079_T33TVG_20200731T102505.3", nir = "subset_0_of_S2A_MSIL1C_20200731T095041_N0209_R079_T33TVG_20200731T102505.4", indices = "NDVI")
ndvi_m_as <- stretch(ndvi_m_a, minv=0, maxv=255)
storage.mode(ndvi_m_as[]) = "integer" 

# RGB images
ggRGB(m_amaro, r = 4, g = 2, b = 3,stretch = "lin") + ggtitle("False colors Sentinel-2 image (10m)") + theme_light() + theme(plot.title = element_text(size = rel(1))) 
ggRGB(m_amaro, r = 3, g = 2, b = 1,stretch = "lin") + ggtitle("Natural colors Sentinel-2 image (10m)") + theme_light() + theme(plot.title = element_text(size = rel(1))) 

# NDVI visualisation  
ggR(ndvi_m_as, geom_raster = TRUE) + ggtitle("NDVI") + theme_light() + theme(plot.title.position ='plot', plot.title = element_text(hjust = 0.5, size = rel(1))) 


# computation indices second crop with rasterdiv functions
renyi1 <- Renyi(ndvi_m_as, window=9, alpha=1, base=exp(1), rasterOut=TRUE, np=3, na.tolerance=0.9, cluster.type="SOCK", debugging=FALSE)
renyi70 <- Renyi(ndvi_m_as, window=9, alpha=70, base=exp(1), rasterOut=TRUE, np=3, na.tolerance=0.9, cluster.type="SOCK", debugging=FALSE)
rao_m_a <- Rao(ndvi_m_as, dist_m="euclidean", window=9, rasterOut = TRUE, mode="classic",lambda=0, shannon=FALSE, rescale=FALSE, na.tolerance=0.9, simplify=3, np=3, cluster.type="SOCK", debugging=FALSE)

# indices values visualisation 
ggR(renyi1, geom_raster = TRUE) + scale_fill_scico(palette = "hawaii") + ggtitle("Rényi's index alpha=1") + theme_light() + theme(plot.title.position ='plot', plot.title = element_text(hjust = 0.5))
ggR(renyi70, geom_raster = TRUE) + scale_fill_scico(palette = "hawaii") + ggtitle("Rényi's index alpha=70") + theme_light() + theme(plot.title.position ='plot', plot.title = element_text(hjust = 0.5))
ggR(rao_m_a, geom_raster = TRUE) + scale_fill_scico(palette = "hawaii") + ggtitle("Rao's Q heterogeneity index") + theme_light() + theme(plot.title.position ='plot', plot.title = element_text(hjust = 0.5))  

