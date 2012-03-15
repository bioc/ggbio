setGeneric("geom_rect", function(data, ...) standardGeneric("geom_rect"))
setMethod("geom_rect", "data.frame", function(data, ...){
  args <- as.list(match.call(call = sys.call(sys.parent(2)))[-1])
  do.call(ggplot2::geom_rect, args)
  ## ggplot2::geom_rect(data, ...)
})
## alignment should be convenient toggle with chevron...
setMethod("geom_rect", "GRanges", function(data,...,
                                           xlab, ylab, main,
                                           facets = NULL,
                                           stat = c("stepping", "identity"),
                                           rect.height = 0.4,
                                           group.selfish = TRUE){

  args <- as.list(match.call(call = sys.call(sys.parent(2)))[-1])
  args.aes <- parseArgsForAes(args)
  args.non <- parseArgsForNonAes(args)
  args.facets <- subsetArgsByFormals(args, facet_grid, facet_wrap)
  args.non <- args.non[!names(args.non) %in% c("data", "facets", "rect.height", "geom", "stat",
                                               "xlab", "ylab", "main")]
  facet <- .buildFacetsFromArgs(data, args.facets)
  
  stat <- match.arg(stat)

  rect.height <- force(rect.height)
  
  if(stat == "stepping"){
    ## if(rect.height <= 0 | rect.height >= 0.5)
    ##   stop("rect.height must be a value in (0,0.5)")
    
    grl <- splitByFacets(data, facets)
    res <- endoapply(grl,
                     function(dt){
                       if("group" %in% names(args.aes))
                         dt <- addStepping(dt, group.name = as.character(args.aes$group),
                                            group.selfish = group.selfish)
                       else
                         dt <- addStepping(dt)
                     })
    res <- unlist(res)
    df <- fortify(model = res)

    args.aes <- args.aes[!(names(args.aes) %in% c("xmin", "xmax", "ymin", "ymax", "data"))]
    args.non <- args.non[!(names(args.non) %in% c("xmin", "xmax", "ymax", "ymax", "data"))]
    if("group" %in% names(args.aes))
      gpn <- as.character(args.aes$group)
    else
      gpn <- "stepping"
    args.aes <- args.aes[names(args.aes) != "group"]
    args.aes <- c(args.aes, list(xmin = substitute(start),
                                 xmax = substitute(end),
                                 ymin = substitute(stepping - rect.height),
                                 ymax = substitute(stepping + rect.height)))

    args.aes <- args.aes[names(args.aes) != "size"]
    aes.res <- do.call(aes, args.aes)
    args.res <- c(list(data = df), list(aes.res),
                  args.non)
    p <- list(do.call(ggplot2::geom_rect,args.res))
    p <- .changeStrandColor(p, args.aes)
    .df.lvs <- unique(df$stepping)
    .df.sub <- df[, c("stepping", gpn)]
    .df.sub <- .df.sub[!duplicated(.df.sub$stepping),]
    ## FIXME:
    if(gpn != "stepping" & group.selfish){
      p <- c(p , list(scale_y_continuous(breaks = .df.sub$stepping,
                                         labels = as.character(.df.sub[, gpn]))))
    } else{
      p <- c(p, list(scale_y_continuous(breaks = NULL)))
    }
  }
  
  if(stat == "identity"){
    if(!"y" %in% names(args.aes)){
      if(!all(c("ymin","ymax", "x", "xmax") %in% names(args.aes))){
        stop("aes(xmin =, xmax= , ymin =, ymax= ) is required for stat 'identity',
              you could also specify aes(y =) only as alternative")
      }
    }else{
      .y <- args.aes$y
      args.aes$xmin <- as.name("start")
      args.aes$xmax <- as.name("end")
      args.aes$ymin <- substitute(y + rect.height, list(y = .y, rect.height = rect.height))
      args.aes$ymax <- substitute(y - rect.height , list(y = .y, rect.height = rect.height))
    }
    df <- fortify(data)
    args.aes <- args.aes[names(args.aes) != "group"]
    args.aes <- args.aes[names(args.aes) != "size"]
    
    aes.res <- do.call(aes, args.aes)
    args.res <- c(list(data = df), list(aes.res),
                  args.non)
    p <- list(do.call(ggplot2::geom_rect,args.res))
    p <- .changeStrandColor(p, args.aes)
  }
  p <- c(list(p) , list(facet))

  if(!missing(xlab))
    p <- c(p, list(ggplot2::xlab(xlab)))
  else
    p <- c(p, list(ggplot2::xlab("Genomic Coordinates")))
  if(!missing(ylab))
    p <- c(p, list(ggplot2::ylab(ylab)))
  if(!missing(main))
    p <- c(p, list(opts(title = main)))
  
  p
})
