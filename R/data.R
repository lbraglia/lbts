#' fmri data
#'
#' subjects are at rest
#' 
#' 24 subjects
#' 70 brain region for each subject series
#' 2 set of test: test=1, retest=2
#' 
#' 404 misllisecond time series
#' 
#' @export
fmri <- function(pt=18, brain_region=64, test_retest = 1){
  unname(FMRI[brain_region,, pt, test_retest])
}
