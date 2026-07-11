#' Read a YAML-like Configuration File
#'
#' Reads a simple key-value configuration file and returns a list of parameters.
#' This provides a lightweight alternative to the yaml package when formatting
#' is simple.
#'
#' @param file Character. Path to configuration file.
#' @param format Character. Config format. \code{"simple"} for key=value or
#'   key: value lines; \code{"auto"} auto-detects.
#'
#' @return A list of configuration parameters.
#'
#' @details
#' The simple format supports:
#' - Comments starting with \code{#}
#' - Key-value pairs separated by \code{=} or \code{:}
#' - Section headers as \code{[section]} for organizational grouping
#' - Blank lines are ignored
#'
#' @export
read_config <- function(file, format = c("auto", "simple")) {
  format <- match.arg(format)
  if (!file.exists(file)) stop("Configuration file not found: ", file)

  lines <- readLines(file, warn = FALSE)
  lines <- trimws(lines)
  # Remove comments (unless they're inside values — not supported in simple format)
  lines <- gsub("\\s*#.*$", "", lines)
  lines <- lines[lines != ""]

  config <- list()
  section <- NULL

  for (line in lines) {
    # Section header
    if (grepl("^\\[", line)) {
      section <- gsub("\\[|\\]", "", line)
      section <- trimws(section)
      next
    }

    # Parse key=value or key: value
    parts <- regexec("^([^:=]+)[:=](.+)$", line)
    m <- regmatches(line, parts)[[1]]
    if (length(m) == 3) {
      key <- trimws(m[2])
      value <- trimws(m[3])

      # Try to coerce value types
      value <- .coerce_config_value(value)

      # Nest under section if applicable
      if (!is.null(section)) {
        if (is.null(config[[section]])) config[[section]] <- list()
        config[[section]][[key]] <- value
      } else {
        config[[key]] <- value
      }
    }
  }

  config
}

.coerce_config_value <- function(value) {
  # Handle quoted strings
  if (grepl('^".*"$', value) || grepl("^'.*'$", value)) {
    return(substr(value, 2, nchar(value) - 1))
  }
  # Numeric
  if (grepl("^[0-9.eE+-]+$", value)) {
    return(as.numeric(value))
  }
  # Integer
  if (grepl("^[0-9]+$", value)) {
    return(as.integer(value))
  }
  # Logical
  if (tolower(value) %in% c("true", "yes")) return(TRUE)
  if (tolower(value) %in% c("false", "no")) return(FALSE)
  # NULL
  if (tolower(value) == "null" || tolower(value) == "none") return(NULL)
  # Default as string
  value
}

#' Validate Configuration Parameters
#'
#' Checks that a configuration list contains all required fields and that
#' values are of the expected types.
#'
#' @param config List. Configuration from \code{\link{read_config}}.
#' @param required_fields Character vector. Names of required top-level fields.
#' @param param_schema Optional list. Named list specifying expected types
#'   (e.g., \code{list(threshold = "numeric", method = "character")}).
#'
#' @return Invisibly returns \code{TRUE} if all checks pass; throws an error
#'   otherwise.
#'
#' @export
validate_config <- function(config, required_fields = NULL, param_schema = NULL) {
  # Check required fields
  if (!is.null(required_fields)) {
    missing <- setdiff(required_fields, names(config))
    if (length(missing) > 0) {
      stop("Missing required config fields: ",
           paste(missing, collapse = ", "))
    }
  }

  # Validate types
  if (!is.null(param_schema)) {
    for (name in names(param_schema)) {
      if (name %in% names(config)) {
        expected_type <- param_schema[[name]]
        actual_value <- config[[name]]
        if (expected_type == "numeric" && !is.numeric(actual_value)) {
          stop(sprintf("Config field '%s' should be numeric, got %s",
                       name, class(actual_value)[1]))
        }
        if (expected_type == "character" && !is.character(actual_value)) {
          stop(sprintf("Config field '%s' should be character, got %s",
                       name, class(actual_value)[1]))
        }
        if (expected_type == "logical" && !is.logical(actual_value)) {
          stop(sprintf("Config field '%s' should be logical, got %s",
                       name, class(actual_value)[1]))
        }
      }
    }
  }

  invisible(TRUE)
}
