.api_header <- function(x) x@api_header

#' @name oncoKB-class
#'
#' @title A class for representing the oncoKB API protocol
#'
#' @description The `oncoKB` class is a representation of the `oncoKB`
#'     API protocol that directly inherits from the `Service` class in the
#'     `AnVIL` package. For more information, see the
#'     \link[AnVIL:Service]{AnVIL} package.
#'
#' @slot api_header `named character()` vector passed on to the `.headers`
#'   argument in `add_headers()`
#'
#' @details This class takes the static API as provided at
#'     \url{https://www.oncoKB.org/api/api-docs} and creates an R object
#'     with the help from underlying infrastructure (i.e.,
#'     \link[rapiclient:rapiclient-package]{rapiclient} and
#'     \link[AnVIL:Service]{AnVIL}) to give the user a unified representation
#'     of the API specification provided by the oncoKB group. Users are not
#'     expected to interact with this class other than to use it as input
#'     to the functionality provided by the rest of the package.
#'
#' @importFrom methods new
#'
#' @seealso  \link{oncoKB}, \linkS4class{Service}
#'
#' @return An `oncoKB` class instance.
#'
#' @md
#'
#' @examples
#' oncoKB()
#'
#' @exportClass oncoKB
.oncoKB <- setClass(
    "oncoKB",
    contains = "Service",
    slots = c(api_header = "character")
)

#' @describeIn oncoKB-class
#'
#' @importFrom AnVIL operations
#' @importFrom methods callNextMethod
#'
#' @param x A \linkS4class{Service} instance or API representation as
#'     given by the \link{oncoKB} function.
#'
#' @inheritParams AnVIL::operations
#'
#' @export
setMethod(
    "operations", "oncoKB",
    function(x, ..., .deprecated = FALSE)
{
    callNextMethod(
        x, .headers = .api_header(x), ..., .deprecated = .deprecated
    )
})
