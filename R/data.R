#'
#'
#' ip,input are simulated data, simulating reads count data from a 200-site MeRIP-Seq experiment with 30 sample conditions..
#' ip :
#'    Three biclusters are implanted in the data, the sizes of which are 80×13,8×8 and 80×12, respectively.
#'    The number of reads for each site under each condition in the IP matrix is generated following (2) in the paper EBBM.
#'    The parameters p of the binomial distribution, which the three biclusters and background follow, are set to 0.65, 0.94, 0.75, and 0.5 respectively,
#'    to make the whole distribution statistic of the simulated data similar to the real data.The three biclusters did not overlap in sites.
#'    In conditions, there are 5 overlaps between the first two biclusters, and 3 overlaps between the last two biclusters
#'
#'
#'
#' @examples
#'   head(ip)
"ip"

#' input
#' input :
#'    Three biclusters are implanted in the data, the sizes of which are 80×13,8×8 and 80×12, respectively.
#'    The number of reads for each site under each condition in the IP matrix is generated following (2) in the paper EBBM.
#'    The parameters p of the binomial distribution, which the three biclusters and background follow, are set to 0.65, 0.94, 0.75, and 0.5 respectively,
#'    to make the whole distribution statistic of the simulated data similar to the real data.The three biclusters did not overlap in sites.
#'    In conditions, there are 5 overlaps between the first two biclusters, and 3 overlaps between the last two biclusters
#'
#'
#'
#'
#' @examples
#'   head(input)
"input"
