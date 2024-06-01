get_mic <- function(x,
                    ids,
                    id_col,
                    ab_col,
                    as_mic = TRUE,
                    prefer_high_mic = TRUE,
                    simplify = TRUE) {
  pre_sort_ids <- data.frame(id = x[[id_col]])
  x <- x[order(x[[id_col]], x[[ab_col]], decreasing = prefer_high_mic),]
  x <- x[!duplicated(x[[id_col]]),]

  output <- data.frame(ids)
  names(output) <- id_col
  output <- merge(output, x,
        sort = FALSE,
        all.x = TRUE)
  if (as_mic) {
    output[[ab_col]] <- AMR::as.mic(output[[ab_col]])
  }
  if (simplify) {
    return(output[[ab_col]])
  } else {
    return(output)
  }
}
