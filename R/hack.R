###  hack at automaticaly generating method for IRanges and Granges
### to avoid get global .method, use a closure.
setOldClass("uneval")
.geoms.ggbio <- paste0("geom_", .ggbio.geom)
.stats.ggbio <- paste0("stat_", .ggbio.stat)
.geoms.ggplot <- paste0("geom_", .ggplot.geom)
.stats.ggplot <- paste0("stat_", .ggplot.stat)

.layouts <-  c("layout_circle", "layout_karyogram")

.gr.name.ggbio <- c(.geoms.ggbio, .stats.ggbio, .layouts)
.gr.name.ggbio <- setdiff(.gr.name.ggbio, c(.geoms.ggplot, .stats.ggplot))
.gr.name.ggplot <- c(.geoms.ggplot, .stats.ggplot)


for(method in .gr.name.ggbio){
  ## for IRanges
  ifun <- function(method){
    .method <- method
    if(hasMethod(.method, "GRanges") && !hasMethod(.method, "IRanges")){
      setMethod(.method, "IRanges", function(data, ...){
        .fun <- selectMethod(.method, signature = "GRanges")
        df <- values(data)
        values(data) <- NULL
        gr <- GRanges("chr_non", data)
        values(gr) <- df
        .fun(gr, ...)
      })
    }
  }
  ifun(method)

  ## for GRangesList

  gfun <- function(method){
    .method <- method
    if(hasMethod(.method, "GRanges") && !hasMethod(.method, "GRangesList")){
      setMethod(.method, "GRangesList", function(data, ...){
        .fun <- selectMethod(.method, signature = "GRanges")
        gr <- biovizBase:::flatGrl(data)
        .fun(gr, ...)
      })
    }
  }
  gfun(method)

  ## hacking for ggplot2-like API without using proto
  ## is data is missing, return a call and parse the data
  mfun <- function(method){
    .method <- method
      setMethod(.method, "missing", function(data,...){
          mc <- match.call()
          mc[-1L] <- list(...)
        return(mc)
      })
  }
  mfun(method)

  ufun <- function(method){
    .method <- method
    setMethod(.method, "uneval", function(data, ...){
      lst <- as.list(match.call())
      idx <- names(lst) != "data"
      aes.u <- unname(lst[!idx])
      res <- lst[idx]
      res <- c(res, aes.u)
      return(as.call(res))
    })
  }
  ufun(method)
}



for(method in .gr.name.ggplot){
  ## for IRanges
  ifun <- function(method){
    .method <- method
    if(hasMethod(.method, "GRanges")) {
      setMethod(.method, "IRanges", function(data, ...){
        .fun <- selectMethod(.method, signature = "GRanges")
        df <- values(data)
        values(data) <- NULL
        gr <- GRanges("chr_non", data)
        values(gr) <- df
        .fun(gr, ...)
      })
    }
  }
  ifun(method)

  ## for GRangesList

  gfun <- function(method){
    .method <- method
    if(hasMethod(.method, "GRanges")) {
      setMethod(.method, "GRangesList", function(data, ...){
        .fun <- selectMethod(.method, signature = "GRanges")
        gr <- biovizBase::flatGrl(data)
        .fun(gr, ...)
      })
    }
  }
  gfun(method)

  ## hacking for ggplot2-like API without using proto
  mfun <- function(method){
      .method <- method

      setMethod(.method, "missing", function(data, ...){
          method0 <- getFromNamespace(method, "ggplot2")
          args <- list(...)
          args.aes <- parseArgsForAes(args)
          args.non <- parseArgsForNonAes(args)
          args.non <- remove_args(args.non, "nbin")
          args <- c(args.non, list(args.aes))
          tm <- try({res <- do.call(method0, args)}, silent = TRUE)
          if(inherits(tm, "try-error")){
              res <-  match.call()
          }else{
              mc <- match.call()
              attr(res, "call") <- TRUE
              attr(res, "mc") <- mc
          }
        return(res)
      })
  }
  mfun(method)

  ufun <- function(method){
      .method <- method
      setMethod(.method, "uneval", function(data, ...){
          method0 <- getFromNamespace(method, "ggplot2")
          args <- list(...)
          args.non <- remove_args(args, "facets")
          args.aes <- data
          args <- c(args.non, list(args.aes))
          tm <- try({res <- do.call(method0, args)}, silent = TRUE)
          if(inherits(tm, "try-error")){
              res <- match.call()
          }else{
              mc <- match.call()
              attr(res, "call") <- TRUE
              attr(res, "mc") <- mc
          }
          return(res)
      })
  }
  ufun(method)
}

