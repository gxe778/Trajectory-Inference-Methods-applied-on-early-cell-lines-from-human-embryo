#HEE Scoup

#load HEE
load("Hee.Rda")

install.packages("SCORPIUS")

library(SCORPIUS)



expression<-HEE@assayData$exprs
group<-HEE@phenoData@data$Sample

DF<-paste(group,expression)
DF<-data.frame(DF)
space <- reduce_dimensionality(expression, "spearman")
draw_trajectory_plot(space, group, contour = TRUE)

traj <- infer_trajectory(space)
draw_trajectory_plot(space, group, traj$path, contour = TRUE)
