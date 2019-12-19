#' Save a ggplot object
#'
#' Save a ggplot in various output formats
#'
#' @param f A character giving the name of a file
#' @param p A ggplot object
#' @examples
#' library('ggplot2')
#' df = as.data.frame(matrix(runif(20), 10, 2))
#' p = ggplot(df, aes(df[, 1], df[, 2]))
#' #save_plot('Rplot.png', p)
#' @export
save_plot <- function(f, p) {

    # get suffixe of filename
    format <- unlist(strsplit(f, ".", fixed = "TRUE"))
    format <- format[length(format)]

    # dynamic loading of formattion depending of the extension
    if (format == "dat")
        formatFunc <- pdf
    else
        formatFunc <- get(format)

    # save
    if (format %in% c("pdf", "dat"))
        formatFunc(f, width = 10, height = 8)
    else
        formatFunc(
            f,
            width = 10,
            height = 8,
            units = "in",
            res = 200
        )

    if (is.function(p))
        p()
    else
        plot(p)

    invisible(dev.off())
}
