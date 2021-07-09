#' @export
#' @rdname geese-common
names.flock <- function(x) {
  names.geese(x)
}


#' @export
#' @rdname geese-common
print.geese <- function(x, ...) {

  print_geese(x)

}

#' @export
#' @rdname geese-common
print.flock <- function(x, ...) {

  print_geese(x)

}
