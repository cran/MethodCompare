aggregate_data <- function(data) {
  return(data[!duplicated(data$id), ])
}