library(tidyverse)
library(magrittr)
library(plyr)
#library(dplyr)
library(lubridate)
library(ggpmisc)
#df1 <- read.csv('/Net/Groups/BGI/people/relghawi/Lightning_Hybrid/FluxNet_Paper1/Flux_P_UP.csv',header=TRUE, row.names="Datetime")
library(ncdf4) # package for netcdf manipulation
# # package for raster manipulation
library(terra)
library(rgdal) # package for geospatial analysis
library(ggplot2) # package for plotting
library(tibble)
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(ggprism)
library(RcppRoll)
library(bigleaf)

set.seed(75431)


df<- read.csv("/Net/Groups/BGI/people/relghawi/Julia_hyb/Clima-Hyb2/Experiments/predictions3.csv") 

df <- slice_sample(df,n=200,replace = FALSE )

my.formula <- y ~ x

p<-df %>%  ggplot() + geom_abline(slope=1,inter=0)+ 
       geom_smooth(data = df, aes(y=gs,x =α),method = "lm", se=TRUE,alpha=0.2,span = 1.0,n=1000, formula = my.formula) + geom_point(data = df, aes(y=gs,x =α),alpha=0.05, stroke=0)+

   labs(y=expression("Leaf diffusive conductance to water H₂O [mol m⁻² s⁻¹]"), x=expression("Pred Leaf diffusive conductance to water H₂O [mol m⁻² s⁻¹]"))+  scale_color_viridis_c()+
#    geom_density_2d(data = df, aes(y=gs,x =α),alpha=0.3,bins=6)+
    theme(panel.background = element_blank(),axis.line = element_line(colour = "black"),
    axis.text.x=element_text(size=4),axis.text.y=element_text(size=4),legend.title = element_text(size=10),
    axis.title = element_text(size=8),plot.title =element_text(size=15),plot.margin=unit(c(0.3,0.3,0.3,0.3), "cm")) 
p2<- p + scale_x_continuous(guide = "prism_minor",breaks= seq(0,5e-3,1e-3), 
                  limits = c(0,5e-3), expand = c(0,0), minor_breaks = seq(0,5e-3,1e-4))+
                  scale_y_continuous(guide = "prism_minor",breaks= seq(0,5e-3,1e-3), 
                  limits = c(0,5e-3), expand = c(0,0), minor_breaks = seq(0,5e-3,1e-4)) + geom_abline(slope=1,inter=0)

p2<- p2+ stat_poly_eq(formula = my.formula,data =df, aes(y=gs,x =α,label = paste(..adj.rr.label.., sep = "~~~")),rr.digits=3,size=3, parse = TRUE,label.x = "right", label.y = "bottom")                  
ggsave(plot = p2, width =5, height = 5, dpi = 500, filename = "/Net/Groups/BGI/scratch/relghawi/paper_2/GPP_siteID_ICON.png")


p<-df %>%  ggplot() + geom_abline(slope=1,inter=0)+ 
       geom_smooth(data = df, aes(y=y,x =ŷ),method = "lm", se=TRUE,alpha=0.2,span = 1.0,n=1000, formula = my.formula) + geom_point(data = df, aes(y=y,x =ŷ),alpha=0.05, stroke=0)+

   labs(y=expression("Transpiration"), x=expression("Pred Transpiration"))+  scale_color_viridis_c()+
#    geom_density_2d(data = df, aes(y=gs,x =α),alpha=0.3,bins=6)+
    theme(panel.background = element_blank(),axis.line = element_line(colour = "black"),
    axis.text.x=element_text(size=4),axis.text.y=element_text(size=4),legend.title = element_text(size=10),
    axis.title = element_text(size=8),plot.title =element_text(size=15),plot.margin=unit(c(0.3,0.3,0.3,0.3), "cm")) 
p2<- p + scale_x_continuous(guide = "prism_minor",breaks= seq(0,0.015,1e-3), 
                  limits = c(0,0.015), expand = c(0,0), minor_breaks = seq(0,0.015,1e-4))+
                  scale_y_continuous(guide = "prism_minor",breaks= seq(0,0.015,1e-3), 
                  limits = c(0,0.015), expand = c(0,0), minor_breaks = seq(0,0.015,1e-4)) + geom_abline(slope=1,inter=0)

p2<- p2+ stat_poly_eq(formula = my.formula,data =df, aes(y=y,x =ŷ,label = paste(..adj.rr.label.., sep = "~~~")),rr.digits=3,size=3, parse = TRUE,label.x = "right", label.y = "bottom")                  
ggsave(plot = p2, width =5, height = 5, dpi = 500, filename = "/Net/Groups/BGI/scratch/relghawi/paper_2/GPP_siteID_ICON.png")

df<- read.csv("/Net/Groups/BGI/people/relghawi/Julia_hyb/Clima-Hyb2/Experiments/full.csv")

p<-df %>%  ggplot() +  
       geom_smooth(data = df, aes(y=g_lw_un,x =RAD),method = "lm", se=TRUE,alpha=0.2,span = 1.0,n=1000, formula = my.formula) + geom_point(data = df, aes(y=g_lw_un,x =RAD),alpha=0.05, stroke=0)+

   labs(y=expression("Leaf diffusive conductance"), x=expression("Radiation"))+  scale_color_viridis_c()+
#    geom_density_2d(data = df, aes(y=gs,x =α),alpha=0.3,bins=6)+
    theme(panel.background = element_blank(),axis.line = element_line(colour = "black"),
    axis.text.x=element_text(size=4),axis.text.y=element_text(size=4),legend.title = element_text(size=10),
    axis.title = element_text(size=8),plot.title =element_text(size=15),plot.margin=unit(c(0.3,0.3,0.3,0.3), "cm")) 
ggsave(plot = p, width =5, height = 5, dpi = 500, filename = "/Net/Groups/BGI/scratch/relghawi/paper_2/GPP_siteID_ICON.png")

p<-df %>%  ggplot() +  
       geom_smooth(data = df, aes(y=g_lw_un,x =T_AIR),method = "lm", se=TRUE,alpha=0.2,span = 1.0,n=1000, formula = my.formula) + geom_point(data = df, aes(y=g_lw_un,x =T_AIR),alpha=0.05, stroke=0)+

   labs(y=expression("Leaf diffusive conductance"), x=expression("Air Temperature"))+  scale_color_viridis_c()+
#    geom_density_2d(data = df, aes(y=gs,x =α),alpha=0.3,bins=6)+
    theme(panel.background = element_blank(),axis.line = element_line(colour = "black"),
    axis.text.x=element_text(size=4),axis.text.y=element_text(size=4),legend.title = element_text(size=10),
    axis.title = element_text(size=8),plot.title =element_text(size=15),plot.margin=unit(c(0.3,0.3,0.3,0.3), "cm")) 
ggsave(plot = p, width =5, height = 5, dpi = 500, filename = "/Net/Groups/BGI/scratch/relghawi/paper_2/GPP_siteID_ICON.png")


p<-df %>%  ggplot() +  
       geom_smooth(data = df, aes(y=g_lw_un,x =vpd),method = "lm", se=TRUE,alpha=0.2,span = 1.0,n=1000, formula = my.formula) + geom_point(data = df, aes(y=g_lw_un,x =vpd),alpha=0.05, stroke=0)+

   labs(y=expression("Leaf diffusive conductance"), x=expression("VPD"))+  scale_color_viridis_c()+
#    geom_density_2d(data = df, aes(y=gs,x =α),alpha=0.3,bins=6)+
    theme(panel.background = element_blank(),axis.line = element_line(colour = "black"),
    axis.text.x=element_text(size=4),axis.text.y=element_text(size=4),legend.title = element_text(size=10),
    axis.title = element_text(size=8),plot.title =element_text(size=15),plot.margin=unit(c(0.3,0.3,0.3,0.3), "cm")) 
ggsave(plot = p, width =5, height = 5, dpi = 500, filename = "/Net/Groups/BGI/scratch/relghawi/paper_2/GPP_siteID_ICON.png")
