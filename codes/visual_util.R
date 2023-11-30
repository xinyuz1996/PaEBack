# Visual utility

#### Result visualization ################################################################
library(ggplot2)
library(ggpubr)
library(stringr)

# n=100;N=100;model_name="1";l=3

#' kplot: Visualization for optimla k: -----
#' @model_name = 1 - 4
kplot <- function(n, N, model_name="1", l=3, ymin = 0, ymax=0.5){
  filename = paste("./simu/n",n,"_N",N,"/n", n,"_", model_name, "_ar5_l",l,"_N",N,".RDS", sep="");filename
  final_result <- readRDS(filename)
  k_opt <- unlist( sapply(final_result, '[', 2)  )
  k = as.numeric(names(table(k_opt)))
  frequency = table(k_opt)/length(k_opt)
  data <- data.frame(k,frequency)
  peak = which.max(table(k_opt)/1000)
  if(is.null(ymin)){ymin = min(frequency)*0.95}
  if(is.null(ymax)){ymax = NA}
  # Plot
  koptgg <- ggplot(data, aes(x=k, y=frequency)) +
    geom_segment( aes(x=k, xend=k, y=ymin, yend=frequency), color="grey") +
    geom_point( color="orange", size=.3) +
    theme_light() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.border = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    xlab("optimal k") +  ylim(ymin, ymax) + 
    ylab("frequency") +
    ggtitle(paste("l=",l, "",sep="")) +
    annotate(geom="text", x=k[peak] + n/20, y=frequency[peak], label=c(paste("kopt =", k[peak]), sep=""), vjust=1.5, size=4) +
    annotate(geom="point", x=k[peak], y=frequency[peak], size=1.5, shape=19, color="red")
  print(koptgg)
  return(koptgg)
}

kplot(n=n, N=N, model_name="1", l=5, ymin=0, ymax=0.13)

#' kplot: Visualization for optimla k: -----
#' @model_name = 1 - 4
kplotcombineL <- function(n=1000, N=1000, model_name="1", main_folder = "./figure"){
  model_fullname = switch(model_name,
                          "1" = "1. Oracle Model",
                          "2" = "2. Lasso Model",
                          "3" = "3. Fixed Elastic",
                          "4" = "3. Tuned Elastic")
  g1 = kplot(n, N, model_name, 3) + xlab("") + ylab("")
  g2 = kplot(n, N, model_name, 5) + xlab("") + ylab("")
  g3 = kplot(n, N, model_name, 7) + xlab("") + ylab("")
  g = ggarrange(g1, g2, g3,
                ncol = 3, nrow = 1, 
                legend = "top",
                label.x = 0,
                label.y = 0,
                font.label = list(size = 15, face = "bold"),
                common.legend = T, 
                align = "v")
  yvar =  c("Frequency")
  xvar = paste("Selected k for Model ",model_fullname," (n=",n,", N=", N,")", sep="");xvar
  g4 = annotate_figure(g, left = text_grob(yvar, just ="centre", size = 12, rot = 90, hjust = .5, vjust = 2),
                       #top = text_grob(c("Frequency of kopt"), size=15),
                       bottom = text_grob(xvar, just ="centre", size = 12, rot = 0, hjust = .35, vjust = -1))
  width = 10
  height = 3.5
  figurepathsep = paste("./simu/n",n,"_N",N,"/kopt_",model_name,"_", "n",str_count(n,"0"),"N",str_count(N,"0"),".png", sep="");
  cat(figurepathsep)
  figurepathtog = paste(main_folder, "/kopt_",model_name,"_", "n",str_count(n,"0"),"N",str_count(N,"0"),".png", sep="")
  cat(figurepathtog)
  ggsave( filename = figurepathtog, plot = g4, width = width,  height = height,  units = c("in"), dpi = 500)
  return(g4)
}

kplotcombineL(n, N, model_name=2, main_folder = "./simu/n100_N100")
n=100;N=100
# for(m in c(1:4)){
#   kplotcombineL(n, N, model_name=m)
# }
# 
# n=100;N=1000
# for(m in c(1:4)){
#   kplotcombineL(n, N, model_name=m)
# }
# 
# n=100;N=100
# for(m in c(1:4)){
#   kplotcombineL(n, N, model_name=m)
# }





# Result visualization -----

#' Look at the distribution of RMSE(k) versus k for the s^th simulation
#' @param  s The index of which simulation to be visualized.
plotRMSEforASimu <- function(s){
  RMSEs_s <- as.vector(final_result[[s]]$RMSE)
  ks_s <- as.vector(final_result[[s]]$k_list) # <=> as.vector('['(final_result[[s]],3))
  plot(ks_s, RMSEs_s, type="l", main=paste("Simu", s), xlab="k", ylab="RMSE")
}

RMSEplot <- function(version_name, folder_name="AR5_Elas"){
  filename = paste("./simu/",folder_name,"/",version_name,".RDS", sep="")
  final_result <- readRDS(filename)
  str(final_result)
  pdfname <- paste("./simu/", folder_name, "/", version_name, ".pdf", sep="");pdfname
  pdf(pdfname)
  par(mfcol=c(2,1))
  kplot(version_name=version_name, folder_name=folder_name)
  indexes = as.list(sort(sample.int(1000, 16*4) ) )
  par(mfrow=c(4,4), mgp = c(2, 1, 0), mar=c(3, 3, 3, 2))
  lapply(indexes, plotRMSEforASimu)
  dev.off()
}
# kplot("ar5_l3", folder_name=folder_name)
# folder_name = fname="elas.5"
# version_name = "ar5_l3"
# RMSEplot("ar5_l3", folder_name=fname)
# RMSEplot("ar5_l5", folder_name=fname)
# RMSEplot("ar5_l7", folder_name=fname)
# dev.off()