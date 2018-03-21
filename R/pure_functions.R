
status <- function(...) message(Sys.time(), "\t", paste(...))

which_top <- function(x, n = 10) {
  which(x >= sort(x, decreasing = TRUE)[n])
}

outlier_sum <- function(x) {
  # Subtract the median, divide by median absolute deviation.
  x <- (x - median(x)) / mad(x)
  # Find the 25th and 75th quantiles.
  xq <- quantile(x, c(0.25, 0.75)) 
  # Find the interquartile range.
  iqr <- xq[2] - xq[1]
  # Compute a statistic for the positive and negative outliers.
  w1 <- abs(sum(x[x > xq[2] + iqr]))
  w2 <- abs(sum(x[x < xq[1] - iqr]))
  # The outlier-sum statistic is the maximum of these values.
  return(max(w1, w2))
}

counts_to_log2cpm <- function(x) log2(x / sum(x) * 1e6 + 1)

minmax <- function(data, min, max) {
  data[data > max] <- max
  data[data < min] <- min
  return(data)
}

lrt_test <- function(x, y, xmin = 1) {
  lrtX     <- bimodLikData(x)
  lrtY     <- bimodLikData(y)
  lrtZ     <- bimodLikData(c(x,y))
  lrt_diff <- 2 * (lrtX + lrtY - lrtZ)
  return(pchisq(lrt_diff, 3, lower.tail = F))
}

bimodLikData <- function(x, xmin = 0) {
  # Small values.
  x1   <- x[x <= xmin]
  # Big values.
  x2   <- x[x > xmin]
  # Proportion of big values.
  xal  <- minmax(length(x2) / length(x), min = 1e-5, max = (1 - 1e-5))
  likA <- length(x1) * log(1 - xal)
  if (length(x2) < 2) {
    x2_sd <- 1
  } else {
    x2_sd <- sd(x2)
  }
  likB <- length(x2) * log(xal) + sum(dnorm(x2, mean(x2), x2_sd, log = TRUE))
  return(likA + likB)
}

extreme_idx <- function(x, lo = 0.01, hi = 0.01) {
  x < quantile(x, lo) | x > quantile(x, 1 - hi)
}

extreme_n <- function(x, n = 5) {
  rank_x <- rank(x)
  which(rank_x < (n + 0.5) / 2 | rank_x > (length(x) - (n + 0.5) / 2))
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}

fancy_scientific <- function(l) {
  l <- format(l, scientific = TRUE)
  l <- gsub("\\+", "", gsub(".+e(.+)$", "10^\\1", l))
  parse(text=l)
}

get_density <- function(x, y, n = 100, by = NULL) {
  if (is.null(by)) {
    dens <- MASS::kde2d(x = x, y = y, n = n)
    ix   <- findInterval(x, dens$x)
    iy   <- findInterval(y, dens$y)
    ii   <- cbind(ix, iy)
    return(dens$z[ii])
  }
  by <- factor(by)
  zlist <- lapply(levels(by), function(i) {
    idx  <- which(by == i)
    dens <- MASS::kde2d(x = x[idx], y = y[idx], n = n)
    ix   <- findInterval(x[idx], dens$x)
    iy   <- findInterval(y[idx], dens$y)
    ii   <- cbind(ix, iy)
    return(dens$z[ii])
  })
  return(unlist(zlist))
}

scale_rows <- function(...) t(scale(t(...)))

quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

quantile_factor <- function(x, n = 7) {
  retval_breaks <- quantile(x, seq(0, 1, length.out = n))
  retval_cuts <- cut(x, retval_breaks, labels = FALSE, include.lowest = TRUE)
  retval <- factor(signif(unname(retval_breaks[retval_cuts + 1]), 2))
  levels(retval) <- paste(
    "<",
    sitools::f2si(as.numeric(as.character(levels(retval))))
  )
  return(retval)
}

optimize_png <- function(filename) {
  opt_filename <- sprintf(
    "%s-fs8.png", substr(filename, 1, nchar(filename) - 4)
  )
  command <- sprintf(
    "pngquant --ext -fs8.png -- %s && mv -f %s %s",
    filename, opt_filename, filename
  )
  system(command)
}

ggsave_optimize_png <- function(filename, ...) {
  ggplot2::ggsave(filename, ...)
  n <- nchar(filename)
  if (substr(filename, n - 3, n) == ".png") {
    optimize_png(filename)
  }
}

#' Plot the extreme values in one column of a matrix, labeled by row names.
#' @param mat A matrix with row names.
#' @param col A column in the matrix.
#' @param n The number of extreme values to plot.
#' @param x,y,title,subtitle Additional arguments passed to \code{ggplot2::labs}.
#' @return A ggplot2 object.
plot_extreme_values <- function(
  mat, col = 1, n = 20, x = NULL, y = NULL, title = NULL, subtitle = NULL
) {
  dat       <- data.frame(value = mat[,col])
  dat$label <- rownames(mat)
  dat       <- dat[order(dat$value),]
  dat       <- rbind(head(dat, n / 2), tail(dat, n / 2))
  ggplot2::ggplot(
    dat,
    ggplot2::aes(reorder(label, value), value, group = label)
  ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x    = reorder(label, value),
        xend = reorder(label, value),
        y    = 0,
        yend = value
      )
    ) +
    ggplot2::geom_point(size = 2) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = x, y = y, title = title, subtitle = subtitle)
}

plot_pca_loading <- function(pca, pc = 1, n = 20) {
  matpx <- data.frame(PC = pca$x[,pc])
  matpx$gene <- rownames(pca$x)
  pca$variance <- pca$sdev ^ 2 / sum(pca$sdev ^ 2)
  matpx <- matpx[order(matpx$PC),]
  matpx <- rbind(head(matpx, n / 2), tail(matpx, n / 2))
  ggplot(matpx, aes(reorder(gene, PC), PC, group = gene)) +
    geom_segment(
      aes(
        x    = reorder(gene, PC),
        xend = reorder(gene, PC),
        y    = 0,
        yend = PC
      )
    ) +
    geom_point(size = 2) +
    coord_flip() +
    labs(
      x = NULL,
      y = NULL,
      title = sprintf(
        "PC%s (%.1f%%)", pc, 100 * pca$variance[pc]
      )
    )
}

