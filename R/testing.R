#' Extract an example from a unit test file.
#' 
#' @param test File name of the test file, source.
#' @param testPath Path to test files. Defaults to "inst/tests".
#' @param examplePath Path under which the example is stored. Defaults to 
#'   "inst/examples".
#' @details From \option{test}, an example is extracted and written to a file of
#'   the same name, but with removed "test-". The file is saved under 
#'   \option{examplePath}. The special tag #-! can be used to hide code in the 
#'   test which is enabled again in the example.
#'   
#' @author Wolfgang Mader, \email{Wolfgang.Mader@@fdm.uni-freiburg.de}
#'   
#' @export
exmpextr <- function(test, testPath = NULL, examplePath = NULL) {
  
  # Construct paths and file names
  # Input: test file
  if (is.null(testPath)) {
    testPath = file.path("inst", "tests")
  }
  testFile = file.path(testPath, test)
  
  # Output: example file
  if (is.null(examplePath)) {
    examplePath = file.path("inst", "examples")
  }
  exampleFile = file.path(examplePath, substr(test, 6, nchar(test)))
  if (!file.exists(examplePath)) {
    stop(paste("Folder", file.path(getwd(), examplePath), "does not exist. Can not write example file"))
  }
  
  
  # Read test file
  if (file.exists(testFile)) {
    code <- readLines(testFile, warn = FALSE)
  } else {
    stop(paste("File", file.path(getwd(), testFile), "not found."))
  }
  
  
  # Data format
  tag <- "#-!"
  tagStart <- paste0(tag, "Start example code")
  tagEnd <- paste0(tag, "End example code")
  
  
  # Find and check example blocks
  blocksStart <- grep(tagStart, code, fixed = TRUE)
  blocksEnd <- grep(tagEnd, code, fixed = TRUE)
  
  if (length(blocksStart) != length(blocksEnd)) {
    stop("Start and end tags do not match.")
  }
  
  blocks <- cbind(blocksStart, blocksEnd)
  colnames(blocks) <- c("start", "end")
  
  # Check if all Start..End tags are in proper order.
  if (nrow(blocks) == 0) {
    warning(paste("No example found in", file.path(getwd(), testFile)))
    return(1)
  } else if (nrow(blocks) == 1) {
    if (blocks[1, 1] > blocks[1, 2]) {
      stop("Start and end tags are not interleaved.")
    }
  } else {
    for (i in 2:nrow(blocks)) {
      if (!(blocks[i - 1, 1] < blocks[i - 1, 2] &&
            blocks[i - 1, 2] < blocks[i, 1] &&
            blocks[i, 1] < blocks[i, 2])) {
        stop("Start and end tags are not interleaved.")
      }
    }
  }
  
  
  # Parse each example block
  examples <- mapply(function(s, e) {
    # Get the current example block
    blockCode <- code[(s + 1):(e - 1)]
    
    # Inspect lines in order to remove our tag.
    # We do not guarantee, that 'tag' is the first non-whitespace character.
    # Thus, the first occurence of 'tag' is removed irrespective of its position
    # in l. This seems to be a corner case which is not worth introducing an
    # extra if statement.
    blockParsed <- do.call(rbind, lapply(blockCode, function(l) {
      return(sub(tag, "", l))
    }))
  }, s = blocks[, "start"], blocks[, "end"])
  

  # Write examples to file.
  outFile <- file(exampleFile, open = "wt")
  
  if (!is.list(examples)) {
    writeLines(examples, outFile)
  } else {
    lapply(examples, function(block) {
      writeLines(block, outFile)
      cat("\n\n", file = outFile)
    })
  }
  
  close(outFile)

  return(0)
}



#' Extract example from unit tests.
#'
#' @details Calls exmpextr for all unit test files.
#' @param testPath character, path to test folder 
#' @param examplePath character path to example folder
#'
#' @author Wolfgang Mader, \email{Wolfgang.Mader@@fdm.uni-freiburg.de}
#' 
#' @export
extractExamples <- function(testPath = NULL, examplePath = NULL) {
  
  # Construct paths
  # Input: test path
  if (is.null(testPath)) {
    testPath = file.path("inst", "tests")
  }

  # Output: example path
  if (is.null(examplePath)) {
    examplePath = file.path("inst", "examples")
  }

  testList <- dir(testPath, pattern = "test-")
  if (length(testList) == 0 ) {
    warning(paste("No unit tests found under", file.path(getwd(), testPath)))
  }
  invisible(lapply(testList, function(test) {
    exmpextr(test, testPath, examplePath)
  }))
}



#' Open a unit test template.
#' 
#' @param test Name of the unit test used as file name with "test-" prefixed.
#' @param testPath Unit test folder. Defaults to "inst/tests".
#'
#' @author Wolfgang Mader, \email{Wolfgang.Mader@@fdm.uni-freiburg.de}
#' 
#' @export
createExample <- function(test, testPath = NULL) {
  
  # Construct path and file names
  if (is.null(testPath)) {
    testPath = file.path("inst", "tests")
  }
  if (!file.exists(testPath)) {
    stop(paste("Folder", file.path(getwd(), testPath), "does not exist. Can not copy unit test template."))
  }
  testFile = file.path(testPath, paste0("test-", test, ".R"))
  
  dModPath <- path.package("dMod")
  file.copy(file.path(dModPath, "templates", "unitTestTemplate.R"), testFile, overwrite = FALSE)
  file.edit(testFile)
}