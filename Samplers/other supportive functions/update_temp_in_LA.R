update_temp_in_LA <- function(LA, new_temp) {
  LA$log_post <- new_temp*LA$llh + LA$lmp
  return(LA)
}