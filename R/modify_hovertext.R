# p: a ggplot function
# hovertext : a boolean for the use of hovertext (if TRUE) as the attribute
# to parse the onMouseOver text ("text" attribute, if FALSE)
# type: class of graphic among regular, var1D, perm
# p_perm: permutation object
modify_hovertext <- function(p, hovertext = TRUE, type = "regular",
                             perm = NULL) {
  attr <- ifelse(hovertext, "hovertext", "text")
  # identify the order / id of the traces which corresponds to x- and y-axis
  # (should be before the splitting function)

  if (type == "perm") {
    k <- 2
  } else if (type == "var1D") {
    k <- 1
  } else {
    k <- seq(2)
  }

  if (type == "perm") {
    traces <- c(1, 2, 5)
  } else if (type == "cv") {
    traces <- 2
    n <- 3
  } else if (type %in% c("regular", "boot1D")) {
    traces <- which(lapply(
      p$x$data,
      function(x) x$mode == "lines"
    ) == TRUE)

    # length of groups of points without traces and circle points
    if (length(traces) > 1) {
      n <- min(traces) - 1
    } else {
      n <- length(p$x$data)
    }
  } else {
    n <- length(p$x$data)
  }

  if (type != "perm") {
    n <- seq(n)
  } else {
    n <- c(3, 4)
  }

  for (i in n) {
    # For each lines of each group in the legend
    for (j in seq(length(p$x$data[[i]][attr][[1]]))) {
      if (!is.null(p$x$data[[i]][attr][[1]][j])) {
        # Distinguish each duplicate by splitting with "<br>" and separe
        # them in key/value by splitting with ": " (like a dictionnary)
        l_text <- lapply(
          as.list(strsplit(p$x$data[[i]][attr][[1]][j], "<br />")[[1]]),
          function(x) strsplit(x, ": ")[[1]]
        )

        l_text <- unlist(lapply(l_text, function(x) {
          if (!is.na(x[1])) {
            if ((type == "boot" && !is.character2(x[2])) ||
              x[1] %in% paste0("df[, ", k, "]")) {
              round(as.numeric(x[2]), 3)
            } else if (type == "boot1D" && x[1] != "order") {
              if (!is.character2(x[2])) {
                round(as.numeric(x[2]), 3)
              } else {
                x[2]
              }
            } else if (type != "var1D" && x[1] == "resp") {
              x[2]
            } else if (type == "cv" && x[1] == "mean") {
              round(as.numeric(x[2]), 3)
            }
          } else {
            NA
          }
        }))

        if (type == "regular" || grepl("boot", type)) {
          # do not print names because text = FALSE plots have not names
          name <- ifelse(hovertext,
            paste0("name: ", p$x$data[[i]]$text[j], "<br />"),
            ""
          )

          parse <- function(x) gsub(" ?</?.> ?|\n", "", x)
          if (type == "boot") {
            x_lab <- parse(p$x$layout$xaxis$title$text)
            y_lab <- parse(p$x$layout$yaxis$title$text)
          } else if (type == "boot1D") {
            x_lab <- gsub(
              "\\(\\d* bootstraps\\)", "",
              parse(p$x$layout$title$text)
            )
            y_lab <- parse(p$x$layout$annotations[[1]]$text)
            if (length(y_lab) < 1) {
              y_lab <- parse(p$x$data[[length(p$x$data)]]$marker$colorbar$title)
            }
          } else {
            x_lab <- "x"
            y_lab <- "y"
          }
          # Overwrite the onMouseOver text with the (x, y) coordinates and
          # the response if exists
          p$x$data[[i]][attr][[1]][j] <-
            paste0(
              name, x_lab, ": ", l_text[1], "<br />", y_lab, ": ", l_text[2],
              ifelse(
                length(l_text) == 3,
                paste0("<br />response: ", l_text[3]),
                ""
              )
            )
        } else {
          if (type == "perm" || (type == "cv" && i != 2)) {
            if (type == "cv" && i == 3) {
              comb <- p$x$data[[i]]$x[1]
            } else {
              if (i == 4) {
                comb <- p$x$data[[i]]$x[1]
              } else {
                comb <- j
              }
            }
            res <- as.list(round(perm$penalties[comb, ], 3))
            if (type == "perm") {
              label <- gsub("</?i>", "", p$x$layout$yaxis$title$text)
            } else {
              label <- "RMSE"
            }
            l_text <- paste0(
              label,
              ": ",
              l_text
            )
            for (b in seq(length(names(res)))) {
              l_text <- paste0(
                l_text, "<br />",
                paste0(names(res)[b], ": ", res[b])
              )
            }
          }
          suppressWarnings(p$x$data[[i]][attr][[1]][j] <- l_text)
        }
      }
    }
  }

  if (type %in% c("perm", "cv")) {
    if (!is.null(p$x$layout$xaxis$title)) {
      p$x$layout$xaxis$title$text <- paste0(c(
        rep("&nbsp;", 20),
        rep("\n&nbsp;", 3),
        p$x$layout$xaxis$title$text,
        rep("&nbsp;", 20)
      ),
      collapse = ""
      )
    }
    if (!is.null(p$x$layout$yaxis$title)) {
      p$x$layout$yaxis$title$text <- paste0(c(
        rep("&nbsp;", 20),
        p$x$layout$yaxis$title$text,
        rep("&nbsp;", 20),
        rep("\n&nbsp;", 3)
      ),
      collapse = ""
      )
    }
    if (!is.null(p$x$layout$xaxis$tickfont)) {
      p$x$layout$xaxis$tickfont$size <- 9
    }
    if (!is.null(p$x$layout$yaxis$tickfont)) {
      p$x$layout$yaxis$tickfont$size <- 9
    }
  }

  if (type == "boot1D") {
    for (i in traces) {
      p$x$data[[i]]$error_x$width <- NULL
    }
  }

  if (type %in% c("boot1D", "cv", "perm", "var1D")) {
    p <- plotly::config(p, scrollZoom = FALSE)
    p <- plotly::layout(p,
      xaxis = list(fixedrange = TRUE),
      yaxis = list(fixedrange = TRUE)
    )
  }

  # Remove the x- and y- axis onOverMouse
  if (
    type %in% c("regular", "cv", "perm") && (length(traces) > 1) ||
      type == "boot1D"
  ) {
    p <- plotly::style(p, hoverinfo = "none", traces = traces)
  }

  if (type %in% c("perm", "cv")) {
    (plotly::layout(p, margin = list(l = 30, r = 10, b = 20, t = 70)))
  } else {
    p
  }
}
