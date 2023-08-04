#' @title Color utilities
#'
#' @description
#' From given values, generate colors for plotting.
#'
#' @param x Vector of values to evaluate.
#'
#' @return
#' \code{get_palette} returns a vector of 99 hex color values. \code{colvec}
#'     returns a hex color vector of same length as \code{x}, corresponding to
#'     values of \code{x}.
#'
#' @examples
#' # generate data
#' x <- y <- seq(-3, 3, length = 30)
#'
#' # color points by numeric values
#' plot(x, y, col=colvec(x), pch=16)
#'
#' # color points by factor values
#' trmt <- gl(3, 10, labels = c('Q', 'R', 'S'))
#' plot(x, y, col=colvec(rev(trmt)), pch=16)
#'
#' @export
#' @rdname utils_color
`get_palette` <- function() {
  pal <- c('#414487E6','#404688E6','#3F4889E6','#3E4989E6','#3E4C8AE6',
           '#3D4E8AE6','#3C508BE6','#3B528BE6','#3A548CE6','#39558CE6',
           '#38588CE6','#375A8CE6','#365C8DE6','#355E8DE6','#35608DE6',
           '#34618DE6','#33638DE6','#32658EE6','#31678EE6','#30698EE6',
           '#306A8EE6','#2F6C8EE6','#2E6E8EE6','#2D708EE6','#2C718EE6',
           '#2C738EE6','#2B748EE6','#2A768EE6','#2A788EE6','#297A8EE6',
           '#287C8EE6','#287D8EE6','#277F8EE6','#26818EE6','#26828EE6',
           '#25848EE6','#24868EE6','#24878EE6','#23898EE6','#228B8DE6',
           '#228D8DE6','#218F8DE6','#21908CE6','#20928CE6','#20938CE6',
           '#1F958BE6','#1F978BE6','#1F998AE6','#1F9A8AE6','#1E9C89E6',
           '#1F9E89E6','#1FA088E6','#1FA187E6','#20A386E6','#20A486E6',
           '#21A685E6','#22A884E6','#24AA83E6','#25AC82E6','#26AD81E6',
           '#28AE80E6','#2AB07FE6','#2DB27DE6','#2FB47CE6','#32B67AE6',
           '#34B679E6','#37B878E6','#3ABA76E6','#3DBC74E6','#40BD72E6',
           '#43BF71E6','#47C06FE6','#4AC16DE6','#4EC36BE6','#52C569E6',
           '#55C668E6','#59C864E6','#5DC863E6','#60CA60E6','#65CB5EE6',
           '#68CD5BE6','#6DCD59E6','#71CF57E6','#75D054E6','#7AD151E6',
           '#7FD34EE6','#83D44CE6','#87D549E6','#8CD646E6','#90D743E6',
           '#95D840E6','#9AD93CE6','#9FDA3AE6','#A3DA37E6','#A8DB34E6',
           '#ADDC30E6','#B2DD2DE6','#B7DE2AE6','#BBDF27E6')
  return(pal)
}
#' @export
#' @rdname utils_color
`colvec` <- function(x) {
  pal <- get_palette()
  return(pal[cut(as.numeric(x), breaks=length(pal), include.lowest=TRUE)])
}




