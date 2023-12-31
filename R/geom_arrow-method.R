## FIXME: the group.selfish doesn't work
setGeneric("geom_arrow", function(data, ...) standardGeneric("geom_arrow"))

setMethod("geom_arrow", "GRanges", function(data, ...,
                                            xlab, ylab, main, 
                                            angle = 30,
                                            length = unit(0.12, "cm"),
                                            type = "open", 
                                            stat = c("stepping", "identity"),
                                            facets = NULL, arrow.rate = 0.03,
                                            group.selfish = TRUE){



  ## remove width = 1
  idx <- width(data) > 1
  data <- data[idx]
  stat <- match.arg(stat)
  ## shape <- match.arg(shape)
  arrow.r <- max(1L, round(width(range(ranges(data))) * arrow.rate, 0))

  args <- list(...)
  args$facets <- facets
  
  args.aes <- parseArgsForAes(args)
  args.non <- remove_args(parseArgsForNonAes(args), c("fill", "facets"))


  facet <- build_facet(data, args)

  if(length(data)){
    ## small arrow
    if(stat == "stepping"){
      if(!"stepping" %in% colnames(values(data))){
        grl <- splitByFacets(data, facets)
        res <- endoapply(grl, make_addStepping, args.aes, group.selfish)
        data <- unlist(res)
      }
      df <- mold(data)
      lst <- apply(df, 1, function(x){
        x <- as.data.frame(t(x))
        x.s <- as.numeric(as.character(x$start))
        x.e <- as.numeric(as.character(x$end))
        N <- (x.e - x.s) %/% arrow.r
        N <- ifelse(N <= 2, 2, N )
        res <- approx(c(x.s, x.e),
                      rep(as.numeric(as.character(x$stepping)), 2),n = N)
        res.df <- do.call(rbind,lapply(1:N, function(i){
          x
        }))
        res.df$temp.x <- res$x
        .res <- res.df[-N,]
        .res$temp.x2 <- res.df[-1, "temp.x"]
        .res
      })
      res <- do.call(rbind,lst)
      res$stepping <- as.numeric(res$stepping)
      args.aes$x <- as.name("temp.x")
      args.aes$xend <- as.name("temp.x2")
      args.aes$y <- args.aes$yend <- as.name("stepping")

      ## need to split to two direction/maybe three?
      p <- by2(res, res$strand, function(x){
        s <- unique(as.character(x$strand))
        aes.temp <- do.call(aes, args.aes)
        aes.temp <- remove_args(aes.temp, "fill")
        if (identical(s, "+")) {
            args.non$arrow <- arrow(length = length, ends = "last",
                                    type = type, angle = angle)
        } else if (identical(s, "-")) {
            args.non$arrow <- arrow(length = length, ends = "first",
                                    type = type, angle = angle)
        }
        p <- do.call(ggplot2::geom_segment, c(list(data = x), list(aes.temp), args.non))
        p
      })
    }
    if(stat == "identity"){
    if(!"y" %in% names(args.aes)){
      if(!all(c("x","xend", "y", "yend") %in% names(args.aes))){
        stop("aes(x =, xend= , y =, yend= ) is required for stat 'identity',
              you could also specify aes(y =) only as alternative")
      }
    }else{
      .y <- args.aes$y
      args.aes$x <- as.name("start")
      args.aes$xend <- as.name("end")
      args.aes$y <-  args.aes$yend <- .y
    }
      
      df <- mold(data)

    lst <- lapply(split(df, seq_len(nrow(df))), function(x){
        x.s <- x$start
        x.e <- x$end
        N <- (x.e - x.s) %/% arrow.r
        N <- ifelse(N <= 2, 2, N )
        res <- approx(c(x.s, x.e ),
                      rep(0, 2),n = N)
        res.df <- x[rep(1L, N),]
        res.df$start <- res$x
        .res <- res.df[-N,]
        .res$end <- res.df[-1L, "start"]
        .res
      })
      res <- do.call(rbind,lst)

      p <- by2(res, res$strand, function(x){
        s <- unique(as.character(x$strand))        
        p <- switch(s,
                    "+" = {
                      args.non$arrow <- arrow(length = length, ends = "last",
                                              type = type, angle = angle)
                      aes.temp <- do.call(aes, args.aes)
                      do.call(ggplot2::geom_segment, c(list(data = x), list(aes.temp), args.non))
                      
                    },
                    "-" = {
                      args.non$arrow <- arrow(length = length, ends = "first",
                                              type = type, angle = angle)
                      aes.temp <- do.call(aes, args.aes)
                      do.call(ggplot2::geom_segment, c(list(data = x), list(aes.temp), args.non))
                    },
                    "*" = {
                      aes.temp <- do.call(aes, args.aes)                    
                      do.call(ggplot2::geom_segment, c(list(data = x), list(aes.temp), args.non))
                    })
        p
      })
    }
  }else{
    p <- NULL
  }
    p <- c(list(p) , list(facet))

  labels <- Labels(xlab, ylab, main, fallback = c(x = ""))
  p <- c(p, labels)
  p

})


