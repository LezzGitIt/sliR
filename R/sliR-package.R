#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom rlang .data :=
## usethis namespace: end
NULL

# NSE column names used inside dplyr verbs. Declared so R CMD check does not flag them as undefined globals.
utils::globalVariables(c(
  "Append", "Mass", "Temp", "Append_log", "Mass_log",
  "b_sli_avg", "sli", "slope", ".log_app", ".log_mass"
))
