#' Bythotrephes Functional Response Data
#'
#' Functional response dataset for \emph{Bythotrephes} spp. (water fleas)
#' preying on prey items of different sizes.
#'
#' @docType data
#' @name bythotrephes
#' @format A data frame with the following columns:
#' \describe{
#'   \item{n_prey_initial}{An integer. The initial density of prey.}
#'   \item{n_prey_eaten}{An integer. The number of prey eaten.}
#'   \item{n_prey_remaining}{An integer. The number of prey left alive.}
#'   \item{size}{A factor with levels \code{'small'}, \code{'medium'}, and
#'   \code{'large'}. The size of prey items.}
#' }
#' @details \emph{Bythotrephes} spp. (water fleas) preying on prey items of
#' different sizes. Prey were not replaced during the experiment, and the total experimental time was 12 hours. Provides an example dataset for type-III and flexible exponent models.
#' @source Daniel Barrios-O'Neill's Frair GitHub:
#' https://github.com/dpritchard/frair
#' @examples
#' data(bythotrephes)
"bythotrephes"
