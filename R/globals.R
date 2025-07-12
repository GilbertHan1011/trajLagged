# Global variables to avoid R CMD check notes
# These are used in dplyr operations and are not true global variables

utils::globalVariables(c(
  "time",
  "time_bin", 
  "reference",
  "target"
)) 