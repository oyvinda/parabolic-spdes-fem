# Libraries
library(ggplot2)

## Loading
  # Load simulation for gamma = 0.25
  load('../results/Des12_2D_statistics_g025.RData')
  yG025 = ySolAll[[1]][[1]][[2]]
  
  # Load simulation for gamma = 0.75
  load('../results/Des12_2D_statistics_g075.RData')
  yG075 = ySolAll[[1]][[1]][[2]]

## Set up mapping
  # Grid
  xRes = seq(0, 1, length.out = 300)
  yRes = seq(0, 1, length.out = 300)
  xs = rep(xRes, 300)
  ys = rep(yRes, each = 300)
  
  # Matrix to map coefficients to grid
  Amap = fm_evaluator_mesh_2d(mesh = meshList[[2]], as.matrix(cbind(xs, ys)))$A
  Amap = Amap[,meshList[[2]]$int]

## Make figures for gamma = 025
  tIdx = 398
  df1 = data.frame(x = xs,
                   y = ys,
                   u = as.vector(Amap%*%yG025[,tIdx]),
                   time = paste("t =", round((tIdx-1)/(ncol(yG025)-1), 3)))
  tIdx = tIdx + 10
  df2 = data.frame(x = xs,
                   y = ys,
                   u = as.vector(Amap%*%yG025[,tIdx]),
                   time = paste("t =", round((tIdx-1)/(ncol(yG025)-1), 3)))
  tIdx = tIdx + 20
  df3 = data.frame(x = xs,
                   y = ys,
                   u = as.vector(Amap%*%yG025[,tIdx]),
                   time = paste("t =", round((tIdx-1)/(ncol(yG025)-1), 3)))
  df = rbind(df1, df2, df3)
  fig = ggplot(data = df, aes(x, y)) + 
    geom_raster(aes(fill = u)) +
    scale_fill_viridis_c(name = "u(t)", limits = c(-max(abs(df$u))-1e-8,max(abs(df$u))+1e-8)) + 
    coord_fixed() + 
    theme(axis.text.x = element_text(size = 16), 
          axis.title.x =element_text(size = 16),
          axis.text.y = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          legend.text = element_text(size=16),
          legend.title = element_text(size=16),
          strip.text = element_text(size = 16),
          panel.spacing = unit(2, "lines"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    facet_wrap(~time)
  fig
  ggsave(filename = "../figures/Example2D_NonStat_g025.png",
         plot = fig,
         width = 28,
         heigh = 10,
         units = "cm")

## Make figures for gamma = 075
  tIdx = 398
  df1 = data.frame(x = xs,
                   y = ys,
                   u = as.vector(Amap%*%yG075[,tIdx]),
                   time = paste("t =", round((tIdx-1)/(ncol(yG025)-1), 3)))
  tIdx = tIdx + 10
  df2 = data.frame(x = xs,
                   y = ys,
                   u = as.vector(Amap%*%yG075[,tIdx]),
                   time = paste("t =", round((tIdx-1)/(ncol(yG025)-1), 3)))
  tIdx = tIdx + 20
  df3 = data.frame(x = xs,
                   y = ys,
                   u = as.vector(Amap%*%yG075[,tIdx]),
                   time = paste("t =", round((tIdx-1)/(ncol(yG025)-1), 3)))
  df = rbind(df1, df2, df3)
  fig = ggplot(data = df, aes(x, y)) + 
    geom_raster(aes(fill = u)) +
    scale_fill_viridis_c(name = "u(t)", limits = c(-max(abs(df$u))-1e-8,max(abs(df$u))+1e-8)) + 
    coord_fixed() + 
    theme(axis.text.x = element_text(size = 16), 
          axis.title.x =element_text(size = 16),
          axis.text.y = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          legend.text = element_text(size=16),
          legend.title = element_text(size=16),
          strip.text = element_text(size = 16),
          panel.spacing = unit(2, "lines"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    facet_wrap(~time)
  fig
  ggsave(filename = "../figures/Example2D_NonStat_g075.png",
         plot = fig,
         width = 28,
         heigh = 10,
         units = "cm")

  
  
  
  
  
  
  # 
  # 
  # print("solved and saved for gamma = 0.0. Computing other gammas")
  # xRes = seq(0, 1, length.out = 300)
  # yRes = seq(0, 1, length.out = 300)
  # xs = rep(xRes, 300)
  # ys = rep(yRes, each = 300)
  # 
  # kRes = 2
  # tRes = 1
  # Atmp = fm_evaluator_mesh_2d(mesh = meshList[[kRes]], as.matrix(cbind(xs, ys)))$A
  # Atmp = Atmp[,meshList[[kRes]]$int]
  # tSave = 10
  # 
  # library(fields)
  # for(i in 1:(350+40)){
  #   print(i)
  #   zs = Atmp%*%yG075[,i]
  #   image.plot(x = xRes, y = yRes, z = matrix(zs, ncol = 300), zlim = c(-2,2),
  #              main = paste("t =", (i-1)/2^(tSave)))
  #   Sys.sleep(0.05)
  # }
  # 
  # 