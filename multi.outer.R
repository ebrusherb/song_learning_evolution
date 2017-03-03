list_args <- Vectorize( function(a,b) c( as.list(a), as.list(b) ), 
                        SIMPLIFY = FALSE)

make_args_mtx <- function( alist ) {
  Reduce(function(x, y) outer(x, y, list_args), alist)
}

multi.outer <- function(f, ... ) {
  args <- make_args_mtx(list(...))
  apply(args, 1:length(dim(args)), function(a) do.call(f, a[[1]] ) )
}

