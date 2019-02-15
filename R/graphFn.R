# graphFn.R
# Chamith Fonseka
# 30 Dec 2015
#
# Objective: Commonly used functions for graphing CyTOF data with ggplot2
#########################################################################################

## LIBRARIES ##
require(ggplot2)
require(RColorBrewer)
require(scales)
require(reshape2)
require(ggrepel)
require(ggthemes)
require(gridExtra)
require(grid)
require(viridis)
require(ggbeeswarm)

grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
    
    plots <- list(...)
    position <- match.arg(position)
    g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
    legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
    lheight <- sum(legend$height)
    lwidth <- sum(legend$width)
    gl <- lapply(plots, function(x) x + theme(legend.position="none"))
    gl <- c(gl, ncol = ncol, nrow = nrow)
    
    combined <- switch(position,
                       "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                              legend,
                                              ncol = 1,
                                              heights = unit.c(unit(1, "npc") - lheight, lheight)),
                       "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                             legend,
                                             ncol = 2,
                                             widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
    grid.newpage()
    grid.draw(combined)
    
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

asinh_breaks <- function(x) {
  br <- function(r) {
    lmin <- round(log10(r[1]))
    lmax <- round(log10(r[2]))
    lbreaks <- seq(lmin, lmax, by = 1)
    breaks <- 10 ^ lbreaks
  }
  p.rng <- range(x[x > 0], na.rm = TRUE)
  breaks <- br(p.rng)
  if (min(x) <= 0) {breaks <- c(0, breaks)}
  if (sum(x < 0) > 1) { #more negative values that expected from expanding scale that includes zero
    n.rng <- -range(x[x < 0], na.rm = TRUE)
    breaks <- c(breaks, -br(n.rng))
  }
  return(sort(breaks))
}
# test_that("asinh_breaks make sense", {
#   expect_equal(asinh_breaks(c(-0.05, 0, 1, 101)), c(0, 1, 10, 100))
#   expect_equal(asinh_breaks(c(-0.11, -0.05, 0, 1, 101)), c(-0.1, 0, 1, 10, 100))
#   expect_equal(asinh_breaks(c(0, 10, 1001)), c(0, 10, 100, 1000))
#   expect_equal(asinh_breaks(c(0, 0.05, 0.07, 0.1, 0.2)), c(0, 0.1))
#   expect_equal(asinh_breaks(c(0.01, 0.02)), c(0.01))
# })

asinh_trans <- function() {
  trans_new("asinh",
            transform = asinh,
            inverse   = sinh,
            breaks = asinh_breaks)
}

asinh.x <- function(marker) {
  max <- max(marker)
  order.of.mag <- ceiling(log10(max))
  break.intervals <- c(-1, 0, 10^seq(0, order.of.mag))
  return(scale_x_continuous(trans = asinh_trans(), breaks = c(-1, 0, 10^seq(0, order.of.mag)), limits = c(0, 10^order.of.mag), labels = expression(-1, 0, 10^0, 10^1, 10^2, 10^3, 10^4, 10^5, 10^6)[1:(order.of.mag + 3)]))
}

asinh.y <- function(marker) {
  max <- max(marker)
  order.of.mag <- ceiling(log10(max))
  break.intervals <- c(-1, 0, 10^seq(0, order.of.mag))
  return(scale_y_continuous(trans = asinh_trans(), breaks = c(-1, 0, 10^seq(0, order.of.mag)), limits = c(0, 10^order.of.mag), labels = expression(-1, 0, 10^0, 10^1, 10^2, 10^3, 10^4, 10^5, 10^6)[1:(order.of.mag + 3)]))
}

SetLegend <- function(alpha = 1, size = 1) {
  return(guides(colour = guide_legend(override.aes = list(alpha = alpha, size = size)), fill = guide_legend(override.aes = list(alpha = alpha, size = size))))
}

# Generic theme settings for ggplot
GenericTheme <- theme(legend.text = element_text(size = 18),
                      legend.text.align = 0,
                      legend.title = element_text(size = 18, face = 'bold'),
                      legend.key = element_rect(fill = NA),
                      panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(), 
                      panel.background = element_blank(), 
                      axis.line.x = element_line(size = 0.5, color = 'black'), 
                      axis.line.y = element_line(size = 0.5, color = 'black'), 
                      strip.text = element_text(size = 24, face = 'bold'), 
                      axis.text = element_text(size = 18, colour = 'black'), 
                      axis.title = element_text(size = 24, face = 'bold'),
                      plot.title = element_text(size = 32, face = 'bold'),
                      plot.subtitle = element_text(size = 28))

# Theme setting for creating SNE plots
sneTheme <- theme(axis.text = element_blank(), 
                  axis.title = element_blank(), 
                  axis.ticks = element_blank(), 
                  axis.line = element_blank(), 
                  plot.title = element_text(size = 28, hjust = 0.5), 
                  panel.border = element_rect(fill = NA, size = 0.5))

# Theme settings for publication, larger font size, etc
PubTheme <- theme(legend.text = element_text(size = 24),
                  legend.text.align = 0,
                  legend.title = element_text(size = 24, face = 'bold'),
                  legend.key=element_rect(fill = NA),
                  panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(), 
                  panel.background = element_blank(), 
                  axis.line.x = element_line(size = 0.5, color = 'black'), 
                  axis.line.y = element_line(size = 0.5, color = 'black'), 
                  strip.text = element_text(size = 32, face = 'bold'), 
                  axis.text = element_text(size = 28, colour = 'black'), 
                  axis.title = element_text(size = 32, face = 'bold'),
                  plot.title = element_text(size = 44, face = 'bold'),
                  plot.subtitle = element_text(size = 36))

my.cols <- rev(brewer.pal(11, "RdYlBu"))
flowjo.cols <- c("#0000FF", "#55FFF6", "#4FFF13", "#FCFF00", "#ED5200", "#EA0000")
# # Set colorscale for plotting marker expression
# stain.cols <- rev(c(rainbow(12, end = 4/6)[1:3],
#                     brewer.pal(11, "RdYlBu")[5],
#                     brewer.pal(11, "RdYlGn")[7:8],
#                     colorRampPalette(c(brewer.pal(11, "RdYlGn")[9], brewer.pal(11, "RdYlBu")[9]))(4),
#                     brewer.pal(11, "RdYlBu")[10:11]))

# Better marker expression coloring - less hue variation
stain.cols <- c("#313695", "#008C48", "#FFB900FF", "#FF0000FF")
status.tissue.cols <- c("OA Arthro" = "#6A3D9A", "RA Arthro" = "#FFD8B2", "RA Biopsy" = "#FF7F00")
case.colors <- rev(ggthemes_data$ptol$qualitative[[2]])


markerColorbar <- function(marker) {
  # Calculate order of magnitude by using the 99.99% value to prevent a few, 
  # real outlier signals from completely changing the colorscale
  max <- quantile(marker, probs = 0.9999)
  order.of.mag <- ceiling(log10(max))
  break.intervals <- c(0, 10^seq(0, order.of.mag))
  colorscale <- scale_colour_gradientn(trans = asinh_trans(), breaks = break.intervals, limits = c(0, 10^order.of.mag), colors = stain.cols, 
                                       labels = expression(0, 10^0, 10^1, 10^2, 10^3, 10^4, 10^5, 10^6)[1:(order.of.mag + 2)], 
                                       guide = guide_colorbar(barheight = 8, nbin = 1000)) 
  return(colorscale)
}


# Fancy formatting for scientifc numbers
fancy_scientific <- function(l, digit_precision = 3, hide_significand = FALSE,
                             return_char = FALSE) { 
  # turn in to character string in scientific notation 
  l <- format(l, digits = digit_precision, scientific = TRUE)
  # Fix 0 case
  l <- gsub("0e\\+00","0",l)
  # quote the part before the exponent to keep all the digits 
  l <- gsub("^(.*)e", "'\\1'e", l) 
  # remove + after exponent, if exists. E.g.: (3x10^+2 -> 3x10^2) 
  l <- gsub("e\\+","e",l)
  # turn the 'e+' into plotmath format 
  l <- gsub("e", "%*%10^", l) 
  # convert 1x10^ or 1.000x10^ -> 10^
  if (hide_significand == TRUE) {
    l <- gsub("\\'1[\\.0]*\\'\\%\\*\\%", "", l)
  }
  # return as an expression or as a character vector if specified
  if (return_char == TRUE) {
    return(l)
  } else {
    return(parse(text = l))
  }
} 




# Get density of points in 2 dimensions.
# @param x A numeric vector.
# @param y A numeric vector.
# @param n Create a square n by n grid to compute density.
# @return The density within each square.
get_density <- function(x, y, n = 100) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}



