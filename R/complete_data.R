complete_data <- function(data) {
  return(data[!(is.na(data$y1) | is.na(data$y2)), ])
}