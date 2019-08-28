#' Data set of Liver Toxicity Data, from \emph{mixOmics}
#'
#' This data set contains the expression measure of \eqn{3116} genes
#' and \eqn{10} clinical measurements for \eqn{64} subjects (rats) that
#'  were exposed to non-toxic, moderately toxic or severely toxic doses
#'   of acetaminophen in a controlled experiment.
#'
#'   The data come from a liver toxicity study (Bushel \emph{et al.}, 2007) in
#'   which 64 male rats of the inbred strain Fisher 344 were exposed to
#'   non-toxic (50 or 150 mg/kg), moderately toxic (1500 mg/kg) or severely
#'   toxic (2000 mg/kg) doses of acetaminophen (paracetamol) in a controlled
#'   experiment. Necropsies were performed at 6, 18, 24 and 48 hours after
#'   exposure and the mRNA from the liver was extracted. Ten clinical chemistry
#'    measurements of variables containing markers for liver injury are available
#'     for each subject and the serum enzymes levels are measured numerically.
#'      The data were further normalized and pre-processed by Bushel \emph{et al.} (2007).
#'
#'@usage
#'data(liverToxicity)
#'
#' @format A list containing the following components:
#' \describe{
#'   \item{gene}{data frame with \eqn{64} rows and \eqn{3116} columns. The
#'   expression measure of 3116 genes for the \eqn{64} subjects (rats).}
#'   \item{clinic}{weight of the diamond, in carats.}
#'   \item{treatment}{data frame with \eqn{64} rows and \eqn{4} columns,
#'   containing the treatment information on the \eqn{64} subjects, such as doses of
#'   acetaminophen
#'    and times of necropsies.}
#'  \item{gene.ID}{data frame with \eqn{3116} rows and \eqn{2} columns,
#'    containing geneBank IDs and gene titles of the annotated genes.}
#' }
#'
#' @source The liver toxicity dataset has been downloaded from the \strong{mixOmics} package.
#' \url{http://mixomics.org/methods/pls-da/}.
#'
#' @references{
#'   \insertRef{bushel2007simultaneous}{ddsPLS}
#' }
#'
#' @importFrom Rdpack reprompt
#'
"liverToxicity"

#' Data set of three species of Penicillium fungi, from \emph{sparseLDA}
#'
#' The data set \code{penicilliumYES} has \emph{36} rows and \emph{3754} columns. The variables
#' are 1st order statistics from multi-spectral images of three species of Penicillium fungi:
#'  \emph{Melanoconidium}, \emph{Polonicum}, and \emph{Venetum}. These are the data used in
#'  the Clemmemsen et al "Sparse Discriminant Analysis" paper.
#'
#'@usage
#'data(penicilliumYES)
#'
#' @format This data set contains the following matrices:
#' \describe{
#'   \item{X}{A matrix with 36 columns and 3754 rows. The training and test data. The first 12
#'    rows are \emph{P. Melanoconidium} species, rows 13-24 are \emph{P. Polonicum} species,
#'     and the last 12 rows are \emph{P. Venetum} species. The samples are ordered so that
#'     each pair of three is from the same isolate.}
#'   \item{Y}{A matrix of dummy variables for the training data.}
#'   \item{Z}{Z matrix of probabilities for the subcalsses of the training data.}
#' }
#'
#' @details
#' The X matrix is not normalized.
#'
#' @source
#' \url{http://www.imm.dtu.dk/~lhc}.
#'
#' @references{
#'   \insertRef{clemmensen2007method}{ddsPLS}
#' }
#'
"penicilliumYES"
