  # # # # # # # # # # # # # # # # # # # # # # # # # #
  # # # #                                     # # # #
  # # # #  plot_3D_Backtransform_morphospace  # # # #
  # # # #   Written by N. Chatar June 2021    # # # #
  # # # #                                     # # # #
  # # # # # # # # # # # # # # # # # # # # # # # # # #

plot_3D_bt_shape <- function(GPA, grid_density = c(4,4), axes = c(1,2), color_points= "#696969") 
{

  ## Packages of interest
  packages = c("plot3D", "geomorph")

  ## Load or install&load the packages
    package.check <- lapply(packages,
      FUN = function(x) 
        {
        if (!require(x, character.only = TRUE)) 
          {
          install.packages(x, dependencies = TRUE)
          library(x, character.only = TRUE)
          }
        })
    
  
  # # # # # # # # # # # # # # # # # # # # # # # #
  # # # #                                 # # # #
  # # # #            ARGUMENTs            # # # #
  # # # #                                 # # # #
  # # # # # # # # # # # # # # # # # # # # # # # #
  
  # # #  GPA # # # 
  
  # An object of class gpagen 
  
  # # #  grid_density  # # #

  # a vector c(x,y) defining the number/density of backtransformed shape to be plotted 

  # # #  color_points  # # #

  # Color to be used to plot the landmarks, default = #696969 (Dim gray) 
    
  # # #  Axes # # # 
    
  # a vector c(x,y) defining the the PCA axes to be plotted, default = PC1 & PC2
    
    
  # # # Check 2D or 3D. If 2D, creates an extra column for the z coordinates full of zeros # # #
    
  if (dim(GPA$coords)[2] == 2)
    {
      gpa_array <- abind(GPA$coords, array(0, replace(dim(GPA$coords), 2, 1)), along = 2)
      colnames(gpa_array)[3] <- "Z"
    }
    
  else
      
    {
      gpa_array <- GPA$coords
    }
  
  # # # PCA # # # 
  
  PCA <- gm.prcomp(gpa_array)
  
  eigenvalues <- PCA$d
  scores <- PCA$x

  df_pca <-cbind(as.data.frame(scores[,axes[1]:axes[2]]))
  df_pca <- as.matrix(df_pca)
  
  # Compute the extent of the two first axis in the PCA
  extent_PC1= max(df_pca[,1]) - (min(df_pca[,1]))
  extent_PC2= max(df_pca[,2]) - (min(df_pca[,2]))
  
  # First loop to compute the widest extant in x and y, this value will be then used to shift the shapes
  
  extant_shape_x <- 0
  extant_shape_y <- 0
  
  for(i in 1:grid_density[1])
    
  {
    coord1 <- min(df_pca[,1]) + ((1/grid_density[1])/2*extent_PC1) + (i-1)*((1/grid_density[1])*extent_PC1)
    
    for(j in 1:grid_density[2])
    {
      coord2 <- max(df_pca[,2]) - ((1/grid_density[2])/2*extent_PC2) - (j-1)*((1/grid_density[2])*extent_PC2)
      
      preds <- shape.predictor(gpa_array, x= df_pca, pred = c(coord1, coord2))
      bt_shape <- preds$pred
      
      extant_shape_x_current <- max(bt_shape[,1]) - min(bt_shape[,1])
      extant_shape_y_current <- max(bt_shape[,2]) - min(bt_shape[,2])
      
      if (extant_shape_x_current >  extant_shape_x)
        {
        extant_shape_x <- extant_shape_x_current
        }
      
      if (extant_shape_y_current >  extant_shape_y)
        {
          extant_shape_y <- extant_shape_y_current
        }

    }
    
  }
  
  # Now that we have the max extent in x and y we can compute the shift between the shapes. It will be
  # 10% higher than the extant we just computed just to make sure shapes are well separated. 
  shift_x <- extant_shape_x + 0.1*extant_shape_x
  shift_y <- extant_shape_y + 0.1*extant_shape_y
  
  # Second loop to plot the shapes 
  # coord1 and coord2 are computed based on the grid density and the extant of the first two PCs
  # shape.predictor is then used to obtain the theoritical shape corresponding to coord1 and coord2
  
  for(i in 1:grid_density[1])
    
  {
    coord2 <- max(df_pca[,2]) - ((1/grid_density[2])/2*extent_PC2) - (i-1)*((1/grid_density[2])*extent_PC2)
    
    for(j in 1:grid_density[2])
    {
      
      coord1 <- min(df_pca[,1]) + ((1/grid_density[1])/2*extent_PC1) + (j-1)*((1/grid_density[1])*extent_PC1)
      
      preds <- shape.predictor(GPA$coords, x= df_pca, pred = c(coord1, coord2))
      bt_shape <- preds$pred
      
      if (length(bt_shape[1,]) == 2)
      {
        bt_shape <- cbind(bt_shape, rep(0, length(bt_shape[,1])))
        colnames(bt_shape)[3] <- "z"
      }
      
      # To plot the shape that has just been computed limits of the plot have to be measured,
      # those are given by the maximua and minima on each axis: x, y, z
      
      limit_max_plot <- max(max(bt_shape[,1]), max(bt_shape[,2]), max(bt_shape[,3]))
      limit_min_plot <- min(min(bt_shape[,1]), min(bt_shape[,2]), min(bt_shape[,3]))
      limit_plot <- c(limit_min_plot,limit_max_plot)
      
      # Now there is 4 possible scenarios to plot the shapes that require slightly different syntax in each case
      # i == 1 && j == 1
        # That's the first time we entered the loop and thus the first shape to be plotted.
        # In this case add=FALSE (that will be the only one) and there is no shift to be added. 
      
      # i == 1 && j > 1
      # We are drawing the first line of shape but the very first one is already done.
      # In this case add=TRUE and we need a shift on the x axis. 
      
      # i > 1 && j == 1
      # We jumped to the other lines but are still on the first column.
      # In this case add=TRUE and we need a shift on the y axis.
      
      # else, meaning i > 1 && j > 1
      # We draw all the other shapes, those that are not in the first line or column. 
      # In this case add=TRUE and we need a shift on both the x AND the y axis.
      
      if (i == 1 && j == 1) 
      {
        plot3d(x=bt_shape[,1], y=bt_shape[,2], z=bt_shape[,3], size = 5, add=FALSE, box=FALSE, axes=FALSE, xlim=limit_plot, ylim=limit_plot, zlim=limit_plot, col = color_points)
      } 
      else if(i == 1 && j > 1)
      {
        plot3d(x=bt_shape[,1] + shift_x*(j-1), y=bt_shape[,2], z=bt_shape[,3], size = 5, add=TRUE, box=FALSE, axes=FALSE, xlim=limit_plot + shift_x*(j-1), ylim=limit_plot, zlim=limit_plot, col = color_points)
      } 
      
      else if(i > 1 && j == 1)
      {
        plot3d(x=bt_shape[,1], y=bt_shape[,2] - shift_y*(i-1), z=bt_shape[,3], size = 5, add=TRUE, box=FALSE, axes=FALSE, xlim=limit_plot, ylim=limit_plot - shift_y*(i-1), zlim=limit_plot, col = color_points)
      }
      else
      {
        plot3d(x=bt_shape[,1] + shift_x*(j-1), y=bt_shape[,2] - shift_y*(i-1), z=bt_shape[,3], size = 5, add=TRUE, box=FALSE, axes=FALSE, xlim=limit_plot + shift_x*(j-1), ylim=limit_plot - shift_y*(i-1), zlim=limit_plot, col = color_points)
      }
      
    }
    
  }
}
