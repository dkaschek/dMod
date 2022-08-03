#' Get Symbols and Numeric constants from a character
#'
#' @param char Character vector (e.g. equation)
#' @param exclude Character vector, the symbols to be excluded from the return value
#' 
#' @export
#' 
#' @examples getElements(c("A*AB+B^2"))
#' 
getElements <- function (char, exclude = NULL) 
{
  if (is.null(char)) 
    return(NULL)
  char <- char[char != "0"]
  out <- parse(text = char, keep.source = TRUE)
  out <- utils::getParseData(out)
  names <- out$text[out$token == "SYMBOL" | out$token == "NUM_CONST"]
  if (!is.null(exclude)) 
    names <- names[!names %in% exclude]
  return(names)
}