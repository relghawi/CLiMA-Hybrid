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

#df <- slice_sample(df,n=200,replace = FALSE )
df <- df %>% filter(gs> 0.001)

my.formula <- y ~ x

p<-df %>%  ggplot() + geom_abline(slope=1,inter=0)+ 
       geom_smooth(data = df, aes(y=gs,x =α),method = "lm", se=TRUE,alpha=0.2,span = 1.0,n=1000, formula = my.formula) + 
       geom_point(data = df, aes(y=gs,x =α),alpha=0.5, stroke=0)+

   labs(y=expression("Leaf diffusive conductance to water H₂O [mol m⁻² s⁻¹]"), x=expression("Pred Leaf diffusive conductance to water H₂O [mol m⁻² s⁻¹]"))+  scale_color_viridis_c()+
#    geom_density_2d(data = df, aes(y=gs,x =α),alpha=0.3,bins=6)+
    theme(panel.background = element_blank(),axis.line = element_line(colour = "black"),
    axis.text.x=element_text(size=8),axis.text.y=element_text(size=8),legend.title = element_text(size=10),
    axis.title = element_text(size=8),plot.title =element_text(size=15),plot.margin=unit(c(0.3,0.3,0.3,0.3), "cm")) 
p2<- p + scale_x_continuous(guide = "prism_minor",breaks= seq(0,5e-3,1e-3), 
                  limits = c(0,5e-3), expand = c(0,0), minor_breaks = seq(0,5e-3,1e-4))+
                  scale_y_continuous(guide = "prism_minor",breaks= seq(0,5e-3,1e-3), 
                  limits = c(0,5e-3), expand = c(0,0), minor_breaks = seq(0,5e-3,1e-4)) + geom_abline(slope=1,inter=0)

p2<- p2+ stat_poly_eq(formula = my.formula,data =df, aes(y=gs,x =α,label = paste(..adj.rr.label.., sep = "~~~")),rr.digits=3,size=3, parse = TRUE,label.x = "right", label.y = "bottom")                  
ggsave(plot = p2, width =5, height = 3, dpi = 500, filename = "/Net/Groups/BGI/scratch/relghawi/paper_2/GPP_siteID_ICON.png")


p<-df %>%  ggplot() + geom_abline(slope=1,inter=0)+ 
       geom_smooth(data = df, aes(y=y,x =ŷ),method = "lm", se=TRUE,alpha=0.2,span = 1.0,n=1000, formula = my.formula) + 
       geom_point(data = df, aes(y=y,x =ŷ),alpha=0.5, stroke=0)+

   labs(y=expression("Leaf Transpiration"), x=expression("Pred Leaf Transpiration"))+  scale_color_viridis_c()+
#    geom_density_2d(data = df, aes(y=gs,x =α),alpha=0.3,bins=6)+
    theme(panel.background = element_blank(),axis.line = element_line(colour = "black"),
    axis.text.x=element_text(size=8),axis.text.y=element_text(size=8),legend.title = element_text(size=10),
    axis.title = element_text(size=10),plot.title =element_text(size=15),plot.margin=unit(c(0.3,0.3,0.3,0.3), "cm")) 
p2<- p + scale_x_continuous(guide = "prism_minor",breaks= seq(0,0.015,5e-3), 
                  limits = c(0,0.015), expand = c(0,0), minor_breaks = seq(0,0.015,1e-3))+
                  scale_y_continuous(guide = "prism_minor",breaks= seq(0,0.015,5e-3), 
                  limits = c(0,0.015), expand = c(0,0), minor_breaks = seq(0,0.015,1e-3)) + geom_abline(slope=1,inter=0)

p2<- p2+ stat_poly_eq(formula = my.formula,data =df, aes(y=y,x =ŷ,label = paste(..adj.rr.label.., sep = "~~~")),rr.digits=3,size=3, parse = TRUE,label.x = "right", label.y = "bottom")                  
ggsave(plot = p2, width =5, height = 3, dpi = 500, filename = "/Net/Groups/BGI/scratch/relghawi/paper_2/GPP_siteID_ICON.png")

### hybrid-clima
df1<- read.csv("/Net/Groups/BGI/people/relghawi/Julia_hyb/Clima-Hyb2/ClimaLand_examples/debug.outputgsleaf.csv")
df2<- read.csv("/Net/Groups/BGI/people/relghawi/Julia_hyb/Clima-Hyb2/ClimaLand_examples/clima_ds_full2.csv")

df <- merge(df1, df2)

#df <- slice_sample(df,n=200,replace = FALSE )
df <- df %>% filter(T_VEG_un> 0.0001)

my.formula <- y ~ x

p<-df %>%  ggplot() + geom_abline(slope=1,inter=0)+ 
       geom_smooth(data = df, aes(y=g_lw_un,x =g_lw_un_pred),method = "lm", se=TRUE,alpha=0.2,span = 1.0,n=1000, formula = my.formula) + 
       geom_point(data = df, aes(y=g_lw_un,x =g_lw_un_pred),alpha=0.5, stroke=0)+

   labs(y=expression("Leaf diffusive conductance to water H₂O [mol m⁻² s⁻¹]"), x=expression("Pred Leaf diffusive conductance to water H₂O [mol m⁻² s⁻¹]"))+  scale_color_viridis_c()+
#    geom_density_2d(data = df, aes(y=gs,x =α),alpha=0.3,bins=6)+
    theme(panel.background = element_blank(),axis.line = element_line(colour = "black"),
    axis.text.x=element_text(size=8),axis.text.y=element_text(size=8),legend.title = element_text(size=10),
    axis.title = element_text(size=8),plot.title =element_text(size=15),plot.margin=unit(c(0.3,0.3,0.3,0.3), "cm")) 

p2<- p+ stat_poly_eq(formula = my.formula,data =df, aes(y=g_lw_un,x =g_lw_un_pred,label = paste(..adj.rr.label.., sep = "~~~")),rr.digits=3,size=3, parse = TRUE,label.x = "right", label.y = "bottom")                  
ggsave(plot = p2, width =5, height = 3, dpi = 500, filename = "/Net/Groups/BGI/scratch/relghawi/paper_2/GPP_siteID_ICON.png")


p<-df %>%  ggplot() + geom_abline(slope=1,inter=0)+ 
       geom_smooth(data = df, aes(y=T_VEG_un,x =T_VEG_un_pred),method = "lm", se=TRUE,alpha=0.2,span = 1.0,n=1000, formula = my.formula) + 
       geom_point(data = df, aes(y=T_VEG_un,x =T_VEG_un_pred),alpha=0.5, stroke=0)+

   labs(y=expression("Leaf Transpiration"), x=expression("Pred Leaf Transpiration"))+  scale_color_viridis_c()+
#    geom_density_2d(data = df, aes(y=gs,x =α),alpha=0.3,bins=6)+
    theme(panel.background = element_blank(),axis.line = element_line(colour = "black"),
    axis.text.x=element_text(size=8),axis.text.y=element_text(size=8),legend.title = element_text(size=10),
    axis.title = element_text(size=10),plot.title =element_text(size=15),plot.margin=unit(c(0.3,0.3,0.3,0.3), "cm")) 
p2<- p + scale_x_continuous(guide = "prism_minor",breaks= seq(0,3e-3,5e-4), 
                  limits = c(0,3e-3), expand = c(0,0))+
                  scale_y_continuous(guide = "prism_minor",breaks= seq(0,3e-3,5e-4), 
                  limits = c(0,3e-3), expand = c(0,0)) + geom_abline(slope=1,inter=0)

p2<- p2+ stat_poly_eq(formula = my.formula,data =df, aes(y=T_VEG_un,x =T_VEG_un_pred,label = paste(..adj.rr.label.., sep = "~~~")),rr.digits=3,size=3, parse = TRUE,label.x = "right", label.y = "bottom")                  
ggsave(plot = p2, width =5, height = 3, dpi = 500, filename = "/Net/Groups/BGI/scratch/relghawi/paper_2/GPP_siteID_ICON.png")

### hybrid-clima gc
df1<- read.csv("/Net/Groups/BGI/people/relghawi/Julia_hyb/Clima-Hyb2/ClimaLand_examples/debug.outputgc.csv")
df2<- read.csv("/Net/Groups/BGI/people/relghawi/Julia_hyb/Clima-Hyb2/ClimaLand_examples/clima_ds_full2.csv")

df <- merge(df1, df2)

#df <- slice_sample(df,n=200,replace = FALSE )
df <- df %>% filter(T_VEG_un> 0.0001)

my.formula <- y ~ x

p<-df %>%  ggplot() + geom_abline(slope=1,inter=0)+ 
       geom_smooth(data = df, aes(y=g_lw,x =g_lw_pred),method = "lm", se=TRUE,alpha=0.2,span = 1.0,n=1000, formula = my.formula) + 
       geom_point(data = df, aes(y=g_lw,x =g_lw_pred),alpha=0.5, stroke=0)+

   labs(y=expression("Canopy conductance"), x=expression("Pred Canopy conductance"))+  scale_color_viridis_c()+
#    geom_density_2d(data = df, aes(y=gs,x =α),alpha=0.3,bins=6)+
    theme(panel.background = element_blank(),axis.line = element_line(colour = "black"),
    axis.text.x=element_text(size=8),axis.text.y=element_text(size=8),legend.title = element_text(size=10),
    axis.title = element_text(size=8),plot.title =element_text(size=15),plot.margin=unit(c(0.3,0.3,0.3,0.3), "cm")) 

p2<- p+ stat_poly_eq(formula = my.formula,data =df, aes(y=g_lw,x =g_lw_pred,label = paste(..adj.rr.label.., sep = "~~~")),rr.digits=3,size=3, parse = TRUE,label.x = "right", label.y = "bottom")                  
ggsave(plot = p2, width =5, height = 3, dpi = 500, filename = "/Net/Groups/BGI/scratch/relghawi/paper_2/GPP_siteID_ICON.png")


p<-df %>%  ggplot() + geom_abline(slope=1,inter=0)+ 
       geom_smooth(data = df, aes(y=T_VEG,x =T_VEG_pred),method = "lm", se=TRUE,alpha=0.2,span = 1.0,n=1000, formula = my.formula) + 
       geom_point(data = df, aes(y=T_VEG,x =T_VEG_pred),alpha=0.5, stroke=0)+

   labs(y=expression("Canopy Transpiration"), x=expression("Pred Canopy Transpiration"))+  scale_color_viridis_c()+
#    geom_density_2d(data = df, aes(y=gs,x =α),alpha=0.3,bins=6)+
    theme(panel.background = element_blank(),axis.line = element_line(colour = "black"),
    axis.text.x=element_text(size=8),axis.text.y=element_text(size=8),legend.title = element_text(size=10),
    axis.title = element_text(size=10),plot.title =element_text(size=15),plot.margin=unit(c(0.3,0.3,0.3,0.3), "cm")) 
p2<- p + scale_x_continuous(guide = "prism_minor",breaks= seq(0,6e-5,1e-5), 
                  limits = c(0,6e-5), expand = c(0,0), minor_breaks = seq(0,6e-5,5e-6))+
                  scale_y_continuous(guide = "prism_minor",breaks= seq(0,6e-5,1e-5), 
                  limits = c(0,6e-5), expand = c(0,0), minor_breaks = seq(0,6e-5,5e-6)) + geom_abline(slope=1,inter=0)

p2<- p2+ stat_poly_eq(formula = my.formula,data =df, aes(y=T_VEG,x =T_VEG_pred,label = paste(..adj.rr.label.., sep = "~~~")),rr.digits=3,size=3, parse = TRUE,label.x = "right", label.y = "bottom")                  
ggsave(plot = p2, width =5, height = 3, dpi = 500, filename = "/Net/Groups/BGI/scratch/relghawi/paper_2/GPP_siteID_ICON.png")


### oz site

df<- read.csv("/Net/Groups/BGI/people/relghawi/Julia_hyb/Clima-Hyb2/Experiments/predictions3_oz.csv")

#df <- slice_sample(df,n=200,replace = FALSE )
df <- df %>% filter(gs> 0.004)
df <- df %>% filter(α> 0.004)

my.formula <- y ~ x


p<-df %>%  ggplot() + geom_abline(slope=1,inter=0)+ 
       geom_smooth(data = df, aes(y=gs,x =α),method = "lm", se=TRUE,alpha=0.2,span = 1.0,n=1000, formula = my.formula) + 
       geom_point(data = df, aes(y=gs,x =α),alpha=0.5, stroke=0)+

   labs(y=expression("Leaf diffusive conductance to water H₂O [mol m⁻² s⁻¹]"), x=expression("Pred Leaf diffusive conductance to water H₂O [mol m⁻² s⁻¹]"))+  scale_color_viridis_c()+
#    geom_density_2d(data = df, aes(y=gs,x =α),alpha=0.3,bins=6)+
    theme(panel.background = element_blank(),axis.line = element_line(colour = "black"),
    axis.text.x=element_text(size=8),axis.text.y=element_text(size=8),legend.title = element_text(size=10),
    axis.title = element_text(size=8),plot.title =element_text(size=15),plot.margin=unit(c(0.3,0.3,0.3,0.3), "cm")) 
p2<- p + scale_x_continuous(guide = "prism_minor",breaks= seq(0,1e-2,1e-3), 
                  limits = c(0,1e-2), expand = c(0,0), minor_breaks = seq(0,1e-2,1e-4))+
                  scale_y_continuous(guide = "prism_minor",breaks= seq(0,1e-2,1e-3), 
                  limits = c(0,1e-2), expand = c(0,0), minor_breaks = seq(0,1e-2,1e-4)) + geom_abline(slope=1,inter=0)

p2<- p2+ stat_poly_eq(formula = my.formula,data =df, aes(y=gs,x =α,label = paste(..adj.rr.label.., sep = "~~~")),rr.digits=3,size=3, parse = TRUE,label.x = "right", label.y = "bottom")                  
ggsave(plot = p2, width =5, height = 3, dpi = 500, filename = "/Net/Groups/BGI/scratch/relghawi/paper_2/GPP_siteID_ICON.png")


p<-df %>%  ggplot() + geom_abline(slope=1,inter=0)+ 
       geom_smooth(data = df, aes(y=y,x =ŷ),method = "lm", se=TRUE,alpha=0.2,span = 1.0,n=1000, formula = my.formula) + 
       geom_point(data = df, aes(y=y,x =ŷ),alpha=0.5, stroke=0)+

   labs(y=expression("Leaf Transpiration"), x=expression("Pred Leaf Transpiration"))+  scale_color_viridis_c()+
#    geom_density_2d(data = df, aes(y=gs,x =α),alpha=0.3,bins=6)+
    theme(panel.background = element_blank(),axis.line = element_line(colour = "black"),
    axis.text.x=element_text(size=8),axis.text.y=element_text(size=8),legend.title = element_text(size=10),
    axis.title = element_text(size=10),plot.title =element_text(size=15),plot.margin=unit(c(0.3,0.3,0.3,0.3), "cm")) 

p2<- p+ stat_poly_eq(formula = my.formula,data =df, aes(y=y,x =ŷ,label = paste(..adj.rr.label.., sep = "~~~")),rr.digits=3,size=3, parse = TRUE,label.x = "right", label.y = "bottom")                  
ggsave(plot = p2, width =5, height = 3, dpi = 500, filename = "/Net/Groups/BGI/scratch/relghawi/paper_2/GPP_siteID_ICON.png")

### hybrid-clima oz site
# df1<- read.csv("/Net/Groups/BGI/people/relghawi/Julia_hyb/Clima-Hyb2/ClimaLand_examples/debug.output_Hyb_oz.csv")
# df2<- read.csv("/Net/Groups/BGI/people/relghawi/Julia_hyb/Clima-Hyb2/ClimaLand_examples/oz_full.csv")

df<- read.csv("/Net/Groups/BGI/people/relghawi/Julia_hyb/Clima-Hyb2/ClimaLand_examples/Hyb_oz_merge.csv")

# df <- merge(df1, df2)

#df <- slice_sample(df,n=200,replace = FALSE )
df <- df %>% filter(T_VEG_un_pred> 0.0001)
df <- df %>% filter(T_VEG_un> 0.0001)

my.formula <- y ~ x

p<-df %>%  ggplot() + geom_abline(slope=1,inter=0)+ 
       geom_smooth(data = df, aes(y=g_lw_un,x =g_lw_un_pred),method = "lm", se=TRUE,alpha=0.2,span = 1.0,n=1000, formula = my.formula) + 
       geom_point(data = df, aes(y=g_lw_un,x =g_lw_un_pred),alpha=0.5, stroke=0)+

   labs(y=expression("Leaf diffusive conductance to water H₂O [mol m⁻² s⁻¹]"), x=expression("Pred Leaf diffusive conductance to water H₂O [mol m⁻² s⁻¹]"))+  scale_color_viridis_c()+
#    geom_density_2d(data = df, aes(y=gs,x =α),alpha=0.3,bins=6)+
    theme(panel.background = element_blank(),axis.line = element_line(colour = "black"),
    axis.text.x=element_text(size=8),axis.text.y=element_text(size=8),legend.title = element_text(size=10),
    axis.title = element_text(size=8),plot.title =element_text(size=15),plot.margin=unit(c(0.3,0.3,0.3,0.3), "cm")) 

p2<- p+ stat_poly_eq(formula = my.formula,data =df, aes(y=g_lw_un,x =g_lw_un_pred,label = paste(..adj.rr.label.., sep = "~~~")),rr.digits=3,size=3, parse = TRUE,label.x = "right", label.y = "bottom")                  
ggsave(plot = p2, width =5, height = 3, dpi = 500, filename = "/Net/Groups/BGI/scratch/relghawi/paper_2/GPP_siteID_ICON.png")


p<-df %>%  ggplot() + geom_abline(slope=1,inter=0)+ 
       geom_smooth(data = df, aes(y=T_VEG_un,x =T_VEG_un_pred),method = "lm", se=TRUE,alpha=0.2,span = 1.0,n=1000, formula = my.formula) + 
       geom_point(data = df, aes(y=T_VEG_un,x =T_VEG_un_pred),alpha=0.5, stroke=0)+

   labs(y=expression("Leaf Transpiration"), x=expression("Pred Leaf Transpiration"))+  scale_color_viridis_c()+
#    geom_density_2d(data = df, aes(y=gs,x =α),alpha=0.3,bins=6)+
    theme(panel.background = element_blank(),axis.line = element_line(colour = "black"),
    axis.text.x=element_text(size=8),axis.text.y=element_text(size=8),legend.title = element_text(size=10),
    axis.title = element_text(size=10),plot.title =element_text(size=15),plot.margin=unit(c(0.3,0.3,0.3,0.3), "cm")) 
p2<- p + scale_x_continuous(guide = "prism_minor",breaks= seq(0,1e-2,1e-3), 
                  limits = c(0,6e-2), expand = c(0,0), minor_breaks = seq(0,6e-2,1e-4))+
                  scale_y_continuous(guide = "prism_minor",breaks= seq(0,6e-2,1e-3), 
                  limits = c(0,6e-2), expand = c(0,0), minor_breaks = seq(0,6e-2,1e-4)) + geom_abline(slope=1,inter=0)

p2<- p+ stat_poly_eq(formula = my.formula,data =df, aes(y=T_VEG_un,x =T_VEG_un_pred,label = paste(..adj.rr.label.., sep = "~~~")),rr.digits=3,size=3, parse = TRUE,label.x = "right", label.y = "bottom")                  
ggsave(plot = p2, width =5, height = 3, dpi = 500, filename = "/Net/Groups/BGI/scratch/relghawi/paper_2/GPP_siteID_ICON.png")



###### oz gc canopy level#####
df<- read.csv("/Net/Groups/BGI/people/relghawi/Julia_hyb/Clima-Hyb2/Experiments/predictions3_gc.csv")

df <- df %>% filter(gs> 0.01)

my.formula <- y ~ x


p<-df %>%  ggplot() + geom_abline(slope=1,inter=0)+ 
       geom_smooth(data = df, aes(y=gs,x =α),method = "lm", se=TRUE,alpha=0.2,span = 1.0,n=1000, formula = my.formula) +
        geom_point(data = df, aes(y=gs,x =α),alpha=0.5, stroke=0)+

   labs(y=expression("Canopy conductance"), x=expression("Pred Canopy Conductance"))+  scale_color_viridis_c()+
#    geom_density_2d(data = df, aes(y=gs,x =α),alpha=0.3,bins=6)+
    theme(panel.background = element_blank(),axis.line = element_line(colour = "black"),
    axis.text.x=element_text(size=8),axis.text.y=element_text(size=8),legend.title = element_text(size=10),
    axis.title = element_text(size=10),plot.title =element_text(size=15),plot.margin=unit(c(0.3,0.3,0.3,0.3), "cm")) 
p2<- p + scale_x_continuous(guide = "prism_minor",breaks= seq(0,0.07,1e-2), 
                  limits = c(0,0.07), expand = c(0,0), minor_breaks = seq(0,0.07,5e-3))+
                  scale_y_continuous(guide = "prism_minor",breaks= seq(0,0.07,1e-2), 
                  limits = c(0,0.07), expand = c(0,0), minor_breaks = seq(0,0.07,5e-3))+ geom_abline(slope=1,inter=0)

p2<- p2+ stat_poly_eq(formula = my.formula,data =df, aes(y=gs,x =α,label = paste(..adj.rr.label.., sep = "~~~")),rr.digits=3,size=3, parse = TRUE,label.x = "right", label.y = "bottom")                  
ggsave(plot = p2, width =5, height =3, dpi = 500, filename = "/Net/Groups/BGI/scratch/relghawi/paper_2/GPP_siteID_ICON.png")


p<-df %>%  ggplot() + geom_abline(slope=1,inter=0)+ 
       geom_smooth(data = df, aes(y=y,x =ŷ),method = "lm", se=TRUE,alpha=0.2,span = 1.0,n=1000, formula = my.formula) +
        geom_point(data = df, aes(y=y,x =ŷ),alpha=0.5, stroke=0)+

   labs(y=expression("Canopy Transpiration"), x=expression("Pred Canopy Transpiration"))+  scale_color_viridis_c()+
#    geom_density_2d(data = df, aes(y=gs,x =α),alpha=0.3,bins=6)+
    theme(panel.background = element_blank(),axis.line = element_line(colour = "black"),
    axis.text.x=element_text(size=8),axis.text.y=element_text(size=8),legend.title = element_text(size=10),
    axis.title = element_text(size=10),plot.title =element_text(size=15),plot.margin=unit(c(0.3,0.3,0.3,0.3), "cm")) 

p2<- p+ stat_poly_eq(formula = my.formula,data =df, aes(y=y,x =ŷ,label = paste(..adj.rr.label.., sep = "~~~")),rr.digits=3,size=3, parse = TRUE,label.x = "right", label.y = "bottom")                  
ggsave(plot = p2, width =5, height = 3, dpi = 500, filename = "/Net/Groups/BGI/scratch/relghawi/paper_2/GPP_siteID_ICON.png")


### hybrid clima gc oz
df<- read.csv("/Net/Groups/BGI/people/relghawi/Julia_hyb/Clima-Hyb2/ClimaLand_examples/Hyb_oz_gc_merge.csv")

df <- df %>% filter(T_VEG> 0.0001)
df <- df %>% filter(T_VEG_pred> 0.0001)

my.formula <- y ~ x


p<-df %>%  ggplot() + geom_abline(slope=1,inter=0)+ 
       geom_smooth(data = df, aes(y=g_lw,x =g_lw_pred),method = "lm", se=TRUE,alpha=0.2,span = 1.0,n=1000, formula = my.formula) + 
       geom_point(data = df, aes(y=g_lw,x =g_lw_pred),alpha=0.5, stroke=0)+

   labs(y=expression("Canopy conductance"), x=expression("Pred Canopy conductance"))+  scale_color_viridis_c()+
#    geom_density_2d(data = df, aes(y=gs,x =α),alpha=0.3,bins=6)+
    theme(panel.background = element_blank(),axis.line = element_line(colour = "black"),
    axis.text.x=element_text(size=8),axis.text.y=element_text(size=8),legend.title = element_text(size=10),
    axis.title = element_text(size=8),plot.title =element_text(size=15),plot.margin=unit(c(0.3,0.3,0.3,0.3), "cm")) 

p2<- p+ stat_poly_eq(formula = my.formula,data =df, aes(y=g_lw,x =g_lw_pred,label = paste(..adj.rr.label.., sep = "~~~")),rr.digits=3,size=3, parse = TRUE,label.x = "right", label.y = "bottom")                  
ggsave(plot = p2, width =5, height = 3, dpi = 500, filename = "/Net/Groups/BGI/scratch/relghawi/paper_2/GPP_siteID_ICON.png")


p<-df %>%  ggplot() + geom_abline(slope=1,inter=0)+ 
       geom_smooth(data = df, aes(y=T_VEG,x =T_VEG_pred),method = "lm", se=TRUE,alpha=0.2,span = 1.0,n=1000, formula = my.formula) + 
       geom_point(data = df, aes(y=T_VEG,x =T_VEG_pred),alpha=0.5, stroke=0)+

   labs(y=expression("Canopy Transpiration"), x=expression("Pred Canopy Transpiration"))+  scale_color_viridis_c()+
#    geom_density_2d(data = df, aes(y=gs,x =α),alpha=0.3,bins=6)+
    theme(panel.background = element_blank(),axis.line = element_line(colour = "black"),
    axis.text.x=element_text(size=8),axis.text.y=element_text(size=8),legend.title = element_text(size=10),
    axis.title = element_text(size=10),plot.title =element_text(size=15),plot.margin=unit(c(0.3,0.3,0.3,0.3), "cm")) 
p2<- p + scale_x_continuous(guide = "prism_minor",breaks= seq(0,1e-3,1e-4), 
                  limits = c(0,1e-3), expand = c(0,0), minor_breaks = seq(0,1e-3,5e-5))+
                  scale_y_continuous(guide = "prism_minor",breaks= seq(0,1e-3,1e-4), 
                  limits = c(0,1e-3), expand = c(0,0), minor_breaks = seq(0,1e-3,5e-5))+ geom_abline(slope=1,inter=0)

p2<- p2+ stat_poly_eq(formula = my.formula,data =df, aes(y=T_VEG,x =T_VEG_pred,label = paste(..adj.rr.label.., sep = "~~~")),rr.digits=3,size=3, parse = TRUE,label.x = "right", label.y = "bottom")                  
ggsave(plot = p2, width =5, height = 3, dpi = 500, filename = "/Net/Groups/BGI/scratch/relghawi/paper_2/GPP_siteID_ICON.png")



df<- read.csv("/Net/Groups/BGI/people/relghawi/Julia_hyb/Clima-Hyb2/Experiments/full.csv")

p<-df %>%  ggplot() +  
       geom_smooth(data = df, aes(y=g_lw_un,x =RAD),method = "lm", se=TRUE,alpha=0.2,span = 1.0,n=1000, formula = my.formula) + geom_point(data = df, aes(y=g_lw_un,x =RAD),alpha=0.05, stroke=0)+

   labs(y=expression("Leaf diffusive conductance"), x=expression("Radiation"))+  scale_color_viridis_c()+
#    geom_density_2d(data = df, aes(y=gs,x =α),alpha=0.3,bins=6)+
    theme(panel.background = element_blank(),axis.line = element_line(colour = "black"),
    axis.text.x=element_text(size=8),axis.text.y=element_text(size=8),legend.title = element_text(size=10),
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
