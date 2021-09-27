#' Plot permutation in 2D
#'
#' Plot permutation in 2D
#'
#' @inheritParams plot_var_2D
#' @inheritParams plot2D
#' @inheritParams get_bootstrap
#' @param perm A permutation object (see \code{\link[RGCCA]{rgcca_permutation}})
#' @param bars A character giving the variability among "points" (all the points
#'  are shown), "sd" (standard deviations bars) , "stderr"(bars of standard
#'  deviation divided by sqrt(n)) or "quantile" (for the 0.05-0.95 quantiles bars)
#' @param type A character giving the type of the index to look at (among 'crit' for
#'  the RGCCA criterion and 'zstat' for the pseudo Z-score)
#' @param correction A Character (either 'none' or 'bonferroni') indicating which
#' method is applied for multiple testing correction. Only useful
#' when "type" is set to "zstat" (default is "none").
#' @param display_all A boolean indicating is all parameter combinations have to
#' be displayed (default is FALSE).
#' @param show_legend A boolean indicating if the legend is displayed (default is
#' TRUE).
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 geom_vline
#' @importFrom stats aggregate
plot_permut_2D <- function(
    perm,
    type = "crit",
    cex = 1,
    title = NULL,
    cex_main = 14 * cex,
    cex_sub = 12 * cex,
    cex_point = 3 * cex,
    cex_lab = 10 * cex,
    bars = "points",
    colors = c("red", "grey"),
    correction = "none",
    display_all = FALSE,
    show_legend = TRUE
) {

    xend <- yend <- NULL
    stopifnot(is(perm, "permutation"))
    match.arg(type, c("crit", "zstat"))
    match.arg(bars, c("points","sd","stderr","quantile"))
    match.arg(correction, c("none", "bonferroni"))
    if ((correction %in% c("bonferroni")) && !(type %in% c("zstat")))
        warning("A correction is only applied when 'type' is equal to 'zstat'.
                Parameter 'correction' is thus ignored.")
    for (i in c("cex", "cex_main", "cex_sub", "cex_point", "cex_lab"))
        check_integer(i, get(i))
    check_colors(colors)
    if (length(colors) < 2)
        colors <- rep(colors, 2)
    crit_title = ifelse(perm$call$method %in% c("sgcca","spls"),
                        "SGCCA criterion",
                        "RGCCA criterion")

    switch(type,
           "zstat" =  y_title <- "Z-score",
           "crit"  = y_title <-  crit_title
    )

    check_ncol(list(perm$zstat), 1)

    y    = unlist(perm[type])
    N    = nrow(perm$penalties)

    df <- setNames(data.frame(y,
                              rep("Other parameter set", N),
                              apply(perm$penalties, 1, function(x)
                                  paste0(round(x, 3), collapse = "/"))),
                   c(type, "Non_permuted", "label"))
    idx_order = sort(df[[type]], decreasing = F, index.return = T)$ix
    df        = df[idx_order, ]
    best      = which.max(unlist(perm["zstat"])[idx_order])

    axis <- function(margin){
        element_text(
            face = "italic",
            size = cex_lab * 0.75,
            margin = margin
        )
    }

    if (is.null(title))
        title <- paste0(
            "Permutation scores (",perm$call$n_perms, " runs) \n Best parameters : ",
            df$label[best]
        )
    else
        title <- paste0(title, collapse = " ")

    df$label              = factor(df$label, levels = rev(df$label), ordered = T)
    df$Non_permuted[best] = "Best parameter set"
    breaks                = rev(levels(df$label))
    labels                = as.expression(breaks)
    if (!display_all){
        jitter         = floor(N/50)
        breaks[max(best-jitter, 1):min(best+jitter, N)] = ""
        breaks[best]   = labels[[best]]
    }
    labels[[best]] = bquote(underline(bold(.(labels[[best]]))))

    p <- ggplot(data = df, mapping = aes_string(x = type, y = "label", xmin = 0)) +
        theme_classic() +
        geom_point(aes(shape = Non_permuted, size = Non_permuted)) +
        labs(
            title = title,
            y     = "Combinations",
            x     = y_title
        ) +
        theme_perso(cex, cex_main, cex_sub) +
        theme(axis.text         = element_text(size = 10, face = "bold"),
              axis.title.x      = axis(margin(0, 20, 0, 0)),
              axis.title.y      = axis(margin(20, 0, 0, 0)),
              axis.line         = element_line(size = 0.5),
              axis.ticks        = element_line(size = 0.5),
              axis.ticks.length = unit(2, "mm")) +
        scale_shape_manual("Non Permuted", values = c("Best parameter set" = 17, "Other parameter set" = 2)) +
        scale_size_manual("Non Permuted", values = c("Best parameter set" = 4, "Other parameter set" = 1)) +
        scale_y_discrete(labels = labels, breaks = breaks, guide = guide_axis(check.overlap = T))
    if (type == "zstat"){
        xintercept = qt(p = c(0.975, 0.995, 0.9995), df = perm$call$n_perms -1)
        if (correction == "bonferroni") xintercept = xintercept*length(df$label)
        p <- p + geom_vline(size = 0.5, color = colors[2], linetype = "dashed",
                            xintercept = xintercept)
    } else {
        dft = data.frame(combinations = rep(df$label, NCOL(perm$permcrit)),
                         permcrit     = c(perm$permcrit[idx_order, ]),
                         Permuted     = rep(df$Non_permuted, NCOL(perm$permcrit)))
        if (bars == "points"){
            p$layers = c(geom_boxplot(data = dft,
                                      aes(x      = permcrit,
                                          y      = combinations,
                                          colour = Permuted),
                                      size = 0.8),
                         p$layers)
        }else{
            switch(bars,
                   "sd" = {
                       tab = aggregate(permcrit ~ combinations + Permuted, dft,
                                       function(x) c(mean   = mean(x),
                                                     xstart = mean(x) - sd(x),
                                                     xend   = mean(x) + sd(x)))
                   },
                   "stderr" = {
                       tab = aggregate(permcrit ~ combinations + Permuted, dft,
                                       function(x) c(mean   = mean(x),
                                                     xstart = mean(x) - sd(x)/sqrt(length(x)),
                                                     xend   = mean(x) + sd(x)/sqrt(length(x))))
                   },
                   "quantile" = {
                       tab = aggregate(permcrit ~ combinations + Permuted, dft,
                                       function(x) c(mean   = mean(x),
                                                     xstart = quantile(x,0.05)[[1]],
                                                     xend   = quantile(x,0.95)[[1]]))
                   })
            print(tab)
            p$layers = c(geom_point(data = tab, aes(y      = as.numeric(combinations),
                                                    x      = permcrit[, "mean"],
                                                    colour = Permuted)),
                         geom_segment(data = tab, aes(y    = as.numeric(combinations),
                                                      x    = permcrit[, "xstart"],
                                                      yend = as.numeric(combinations),
                                                      xend = permcrit[, "xend"],
                                                      colour = Permuted), size=0.5),
                         p$layers)
        }
        p <- p + scale_colour_manual("Permuted", values = c("Best parameter set"  = "black",
                                                            "Other parameter set" = "grey"))
    }

    p <- p +
        theme(plot.title  = element_text(vjust=5),
              plot.margin = margin(5, 0, 0, 0, "mm"),
              legend.position = c(0.7, 0.8))
    if (!show_legend){
        p <- p + theme(legend.position = "none")
    }

    attributes(p)$penalties <- perm$penalties[idx_order, ]

    return(p)
}
