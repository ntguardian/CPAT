################################################################################
# Plotting.R
################################################################################
# 2018-09-06
# Curtis Miller
################################################################################
# Functions for creating plots used in papers
################################################################################

################################################################################
# DISTRIBUTION CONVERGENCE PLOTS
################################################################################

#' Create Tikz Plot Demonstrating Rényi-Type Statistic's Convergence in
#' Distribution
#'
#' Create a Tikz file containing a plot demonstrating that the Rényi-type
#' statistic converges in distribution. Optionally, create a PDF as well.
#'
#' @param obj The list containing the simulations
#' @param dist The identifier of the data-generating process that generated the
#'             datasets on which the Rényi-type statistic was computed
#' @param trim The identifier of the trimming parameter of the Rényi-type
#'             statistic
#' @param size The sample size of the simulated data sets
#' @param title The title of the plot
#' @param width The width of the plot
#' @param height The height of the plot
#' @param filename The name of the output file (without extensions; \code{.tex}
#'                 and maybe \code{.pdf} files will be created); if \code{NULL},
#'                 the name will automatically be determined (of the form
#'                 \code{dist_conv_dist_nsize_trim})
#' @param verbose Print updates about progress (via \code{link[base]{cat}})
#' @param makePDF Automatically compile the resulting \code{.tex} file
#' @examples
#' \dontrun{
#' ZnSimulations <- list(
#'   "norm" = list(
#'     "log" = list(
#'       "n500" = c(3.206, 1.32, 0.776, 1.262, 0.676)
#'     )
#'   )
#' )
#' dist_conv_plot_tikz(ZnSimulations, "norm", "log", 500)
#' }
dist_conv_plot_tikz <- function(obj, dist, trim, size, title = "", width = 4,
                                height = 3, filename = NULL, makePDF = TRUE,
                                verbose = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE) |
    !requireNamespace("tikzDevice", quietly = TRUE)) {
    stop("ggplot2 and tikzDevice are needed to use this function")
  }

  qplot <- ggplot2::qplot
  aes <- ggplot2::aes
  geom_line <- ggplot2::geom_line
  xlab <- ggplot2::xlab
  ylab <- ggplot2::ylab
  ggtitle <- ggplot2::ggtitle
  coord_cartesian <- ggplot2::coord_cartesian
  theme_bw <- ggplot2::theme_bw
  unit <- ggplot2::unit
  theme <- ggplot2::theme
  tikz <- tikzDevice::tikz

  if (is.null(filename)) {
    filename <- paste0("dist_conv_", dist, "_n", size, "_", trim)
  }
  filename <- gsub(".", "", filename, fixed = TRUE)
  filename.tex <- paste0(filename, ".tex")
  filename.pdf <- paste0(filename, ".pdf")
  
  x <- seq(0,4, length = 1000)
  
  p <- qplot(obj[[dist]][[trim]][[paste0("n", size)]], geom = "density") +
    geom_line(data = data.frame("xval" = x, "yval" = dZn(x)), aes(x = xval,
                                                                  y = yval),
              linetype = "dashed", size = 1.2) +
    # xlab("$x$") +
    # ylab("Density") +
    xlab("") +
    ylab("") +
    # ggtitle(title) +
    ggtitle("") +
    coord_cartesian(xlim = c(0, 4)) +
    theme_bw() +
    theme(legend.position="none",
      plot.margin = unit(c(.125, .125, .125, .125), "in")) #, 
      # plot.title = element_text(size = rel(.6)),
      # axis.title = element_text(size = rel(.6)),
      # axis.text = element_text(size = rel(.45)))
  if (makePDF) {
    tikz(filename.tex, width = width, height = height, standAlone = TRUE)
    print(p)
    dev.off()
    
    tools::texi2pdf(, clean = TRUE, quiet = verbose, clean = TRUE, quiet = verbose)
  }
  if (verbose) {cat("Creating", filename.tex, "\n")}
  tikz(filename.tex, width = width, height = height)
  print(p)
  dev.off()
}

#' Power Curve Plot
#'
#' Create a Tikz plot of the power curves of simulated statistics.
#'
#' @param data A \code{data.frame} containing the data to plot
#' @param d Label for data-generating process used to simulate the data on which
#'          the statistics were computed
#' @param t Label for the trimming parameter of the Rényi-type statistic
#' @param c Label for the process that computes the location of change points
#' @param N The sample size of the simulated data sets on which the statistics
#'          were computed
#' @param statlines A character vector where the names of the entries are the
#'                  labels of the statistics in the \code{stat} column of
#'                  \code{data} and the entries define the line types used by
#'                  the \code{values} entry of
#'                  \code{\link[ggplot2]{scale_linetype_manual}}
#' @param title The title of the plot
#' @param legend_pos A string to be passed to \code{link[ggplot2]{theme}} (the
#'                   \code{legend.position} argument) identifying where to place
#'                   the legend
#' @param width The width of the plot
#' @param height The height of the plot
#' @param verbose Print progress information
#' @param filename The name of the file to save output (without stems; files
#'                 with this string appended with \code{.tex} and maybe
#'                 \code{.pdf} will be created); if \code{NULL}, the name will
#'                 be automatically determined
#' @param verbose Print updates about progress (via \code{link[base]{cat}})
#' @param makePDF Automatically compile the resulting \code{.tex} file
#' @examples
#' \dontrun{
#' pdat <- data.frame(power = c(0.8926, 0.8714, 0.8296, 0.7936),
#'                    stat  = c("de",   "de",   "de",   "de"),
#'                    dist  = c("norm", "norm", "norm", "norm"),
#'                    kn    = c("log",  "log",  "log",  "log"),
#'                    n     = c(50,     50,     50,     50),
#'                    cpt   = c("c4rt", "c4rt", "c4rt", "c4rt"),
#'                    delta = c(-2.0,   -1.9,   -1.8,   -1.7))
#' power_plot_tikz_by_n(pdat, "norm", "log", "c4rt", 50, c("de" = "solid"))
#' }
power_plot_tikz_by_n <- function(data, d, t, c, N, statlines, title = "",
                                 legend_pos = "none", width = 4.5, height = 3.5,
                                 filename = NULL, verbose = FALSE,
                                 makePDF = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE) |
    !requireNamespace("tikzDevice", quietly = TRUE) |
    !requireNamespace("dplyr", quietly = TRUE)) {
    stop("ggplot2, dplyr, and tikzDevice are needed to use this function")
  }

  ggplot <- ggplot2::ggplot
  aes <- ggplot2::aes
  scale_y_continuous <- ggplot::scale_y_continuous
  scale_linetype_manual <- ggplot::scale_linetype_manual
  guides <- ggplot2::guides
  guide_legend <- ggplot2::guide_legend
  geom_hline <- ggplot2::geom_hline
  geom_line <- ggplot2::geom_line
  unit <- ggplot2::unit
  xlab <- ggplot2::xlab
  ylab <- ggplot2::ylab
  ggtitle <- ggplot2::ggtitle
  coord_cartesian <- ggplot2::coord_cartesian
  theme_bw <- ggplot2::theme_bw
  theme <- ggplot2::theme
  `%>%` <- dplyr::`%>%`
  filter <- dplyr::filter
  tikz <- tikzDevice::tikz

  if (is.null(filename)) {
    filename <- paste("power_plot", d, "n" %s0% N, t, c, sep = "_")
  }
  filename <- gsub(".", "", filename, fixed = TRUE)
  filename.tex <- paste0(filename, ".tex")
  filename.pdf <- paste0(filename, ".pdf")
  
  p <- ggplot(data %>% filter(dist == d & kn == t & cpt == c & n == N & 
      stat %in% names(statlines)),
              aes(x = delta, y = power, linetype = as.factor(stat),
                group = as.factor(stat))) +
    geom_line() +
    scale_y_continuous(breaks = c(0.05, .25, .5, .75, 1.0),
                       labels = c("$\\alpha$", "$0.25$", "$0.50$", "$0.75$",
                         "$1.00$"),
                       limits = c(0, 1)) +
    scale_linetype_manual(values = statlines) +
    # xlab("$\\Delta$") +
    xlab("") +
    # ylab("Power") +
    ylab("") +
    # ggtitle(title) +
    ggtitle("") +
    guides(linetype = guide_legend("$T$")) +
    geom_hline(yintercept = 0.05, color = "black", linetype = "dotted") +
    theme_bw() +
    theme(legend.position = legend_pos, legend.spacing = unit(0, "in"),
      plot.margin = unit(c(.125, .125, .125, .125), "in")) #,
      # strip.background = element_blank(),
      # panel.border = element_rect(color = "black"))
      # plot.title = element_text(size = rel(.8)),
      # axis.title = element_text(size = rel(.8)),
      # axis.text = element_text(size = rel(.7)),
      # legend.title = element_text(size = rel(.75)),
      # legend.text = element_text(size = rel(.7)),
      # plot.margin = unit(c(.125, .125, .125, .125), "in"))
  
  if (makePDF) {
    tikz(filename.tex, width = width, height = height, standAlone = TRUE)
    print(p)
    dev.off()

    tools::texi2pdf(, clean = TRUE, quiet = verbose, clean = TRUE, quiet = verbose)
  }
  if (verbose) {cat("Creating", filename.tex, "\n")}
  tikz(filename.tex, width = width, height = height)
  print(p)
  dev.off()
}

#' Power Curve Plot (By Statistic)
#'
#' Create a Tikz plot of the power curves of a statistic, with each sample size
#' having its own curve.
#'
#' @param s The statistic for which to plot a power curve
#' @inheritParams power_plot_tikz_by_n
#' @examples
#' \dontrun{
#' pdat <- data.frame(power = c(0.8926, 0.8714, 0.8296, 0.7936),
#'                    stat  = c("de",   "de",   "de",   "de"),
#'                    dist  = c("norm", "norm", "norm", "norm"),
#'                    kn    = c("log",  "log",  "log",  "log"),
#'                    n     = c(50,     50,     50,     50),
#'                    cpt   = c("c4rt", "c4rt", "c4rt", "c4rt"),
#'                    delta = c(-2.0,   -1.9,   -1.8,   -1.7))
#' power_plot_tikz(pdat, "norm", "log", "c4rt", "de")
#' }
power_plot_tikz <- function(data, d, t, c, s, title = "", legend_pos = "none",
                            width = 4.5, height = 3.5, filename = NULL,
                            verbose = FALSE, makePDF = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE) |
    !requireNamespace("tikzDevice", quietly = TRUE) |
    !requireNamespace("dplyr", quietly = TRUE)) {
    stop("ggplot2, dplyr, and tikzDevice are needed to use this function")
  }

  ggplot <- ggplot2::ggplot
  aes <- ggplot2::aes
  scale_y_continuous <- ggplot::scale_y_continuous
  guides <- ggplot2::guides
  guide_legend <- ggplot2::guide_legend
  geom_hline <- ggplot2::geom_hline
  geom_line <- ggplot2::geom_line
  unit <- ggplot2::unit
  xlab <- ggplot2::xlab
  ylab <- ggplot2::ylab
  ggtitle <- ggplot2::ggtitle
  coord_cartesian <- ggplot2::coord_cartesian
  theme_bw <- ggplot2::theme_bw
  theme <- ggplot2::theme
  `%>%` <- dplyr::`%>%`
  filter <- dplyr::filter
  tikz <- tikzDevice::tikz

  if (is.null(filename)) {
    filename <- paste("power_plot", d, s, t, c, sep = "_")
  }
  filename <- gsub(".", "", filename, fixed = TRUE)
  filename.tex <- paste0(filename, ".tex")
  filename.pdf <- paste0(filename, ".pdf")
  
  p <- ggplot(data %>% filter(dist == d & kn == t & cpt == c & stat == s),
            aes(x = delta, y = power, linetype = as.factor(n),
              group = as.factor(n))) +
    geom_line() +
    scale_y_continuous(breaks = c(0.05, .25, .5, .75, 1.0),
                       labels = c("$\\alpha$", "$0.25$", "$0.50$", "$0.75$",
                                  "$1.00$"),
                       limits = c(0, 1)) +
    # xlab("$\\Delta$") +
    xlab("") +
    # ylab("Power") +
    ylab("") +
    # ggtitle(title) +
    ggtitle("") +
    guides(linetype = guide_legend("$T$")) +
    geom_hline(yintercept = 0.05, color = "red", linetype = "dashed") +
    theme_bw() +
    theme(legend.position = legend_pos, legend.spacing = unit(0, "in"),
      plot.margin = unit(c(.125, .125, .125, .125), "in")) #,
      # strip.background = element_blank(),
      # panel.border = element_rect(color = "black"))
      # plot.title = element_text(size = rel(.8)),
      # axis.title = element_text(size = rel(.8)),
      # axis.text = element_text(size = rel(.7)),
      # legend.title = element_text(size = rel(.75)),
      # legend.text = element_text(size = rel(.7)),
      # plot.margin = unit(c(.125, .125, .125, .125), "in"))
  
  if (makePDF) {
    tikz(filename.tex, width = width, height = height, standAlone = TRUE)
    print(p)
    dev.off()
    
    tools::texi2pdf(, clean = TRUE, quiet = verbose, clean = TRUE, quiet = verbose)
  }
  if (verbose) {cat("Creating", filename.tex, "\n")}
  tikz(filename.tex, width = width, height = height)
  print(p)
  dev.off()
}

#' Long-Run Variance Estimation Simulations Plot
#'
#' Create a Tikz plot of the estimated distribution of LRV estimators
#'
#' @param data A \code{data.frame} containing the data to plot
#' @param n The sample size of simulated data sets for which to plot an
#'          estimated distribution
#' @param phi The autocorrelation parameter of the simulated data sets to plot;
#'            if \code{NULL}, the data is assumed to have been generated with a
#'            GARCH(1,1) process
#' @param ker_name The name of the kernel function used in the LRV estimator
#' @param true_lrv The value of the true long-run variance
#' @param xrange The limits of the horizontal axis of the plot
#' @param width The width of the plot
#' @param height The height of the plot
#' @param filename The name of the file to save output (without stems; files
#'                 with this string appended with \code{.tex} and maybe
#'                 \code{.pdf} will be created); if \code{NULL}, a file name
#'                 will automatically be chosen (of the form
#'                 \code{lrv_est_plot_ker_name_phi})
#' @param verbose Print updates about progress (via \code{link[base]{cat}})
#' @param makePDF Automatically compile the resulting \code{.tex} file
#' @examples
#' \dontrun{
#' plotlist <- data.frame(val = c(0.649, 0.965, 0.905),
#'                        n   = c(50,    50,    50),
#'                        phi = c(0,     0,     0))
#' lrv_plot_tikz(plotlist, 50, ker_name = "bartlett", true_lrv = 1)
#' }
lrv_plot_tikz <- function(data, n, ker_name, true_lrv, phi = NULL,
                          xrange = NULL, width = 4.5, height = 3.5,
                          filename = NULL, verbose = FALSE, makePDF = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE) |
    !requireNamespace("tikzDevice", quietly = TRUE) |
    !requireNamespace("dplyr", quietly = TRUE)) {
    stop("ggplot2, dplyr, and tikzDevice are needed to use this function")
  }

  ggplot <- ggplot2::ggplot
  aes <- ggplot2::aes
  geom_vline <- ggplot2::geom_vline
  geom_density <- ggplot2::geom_density
  unit <- ggplot2::unit
  xlab <- ggplot2::xlab
  ylab <- ggplot2::ylab
  ggtitle <- ggplot2::ggtitle
  theme_bw <- ggplot2::theme_bw
  theme <- ggplot2::theme
  `%>%` <- dplyr::`%>%`
  filter <- dplyr::filter
  tikz <- tikzDevice::tikz

  if (is.null(phi)) {
    phi <- "garch"
    data$phi <- "garch"
  }
  
  if (is.null(filename)) {
    filename <- paste("lrv_est_plot", ker_name, ifelse(phi == "garch", phi,
        "ar" %s0% "_" %s0% phi), n, sep = "_")
  }
  filename <- gsub(".", "", filename, fixed = TRUE)
  filename.tex <- paste0(filename, ".tex")
  filename.pdf <- paste0(filename, ".pdf")
  
  print("n =" %s% n %s% "phi =" %s% phi)
  n0 <- n
  phi0 <- phi
  
  p <- ggplot(data %>% filter(n == n0 & phi == phi0),
              aes(x = val)) +
    geom_density() +
    geom_vline(aes(xintercept = true_lrv)) +
    # xlab("Estimated Long Run Variance") +
    xlab("") +
    # ylab("Density") +
    ylab("") +
    # ggtitle(title) +
    ggtitle("") +
    xlim(xrange) +
    theme_bw() +
    theme(plot.margin = unit(c(.125, .125, .125, .125), "in"))
  
  print(p)
    
  if (verbose) {cat("Creating", filename.tex, "\n")}
  if (makePDF) {
    tikz(filename.tex, width = width, height = height, standAlone = TRUE)
    print(p)
    dev.off()

    tools::texi2pdf(filename.tex, clean = TRUE, quiet = verbose)
  }
  tikz(filename.tex, width = width, height = height)
  print(p)
  dev.off()
}
