## reference:http://groups.google.com/group/ggplot2/browse_thread/thread/72403c6997b79c3b?pli=1
## "link" and "ribbon" requie a special data structure
## 1. we could implement it as GRangesList with no direction specification
## 2. we could use a to.gr as element meta data which is a granges with direction,
## this is also a GRanges object which is general(implemented this first)
setGeneric("layout_circle", function(data,...) standardGeneric("layout_circle"))
setMethod("layout_circle",  "GRanges",
          function(data, ..., geom = c("point", "line", "link", "ribbon","rect", "bar",
                                       "segment", "rectangle", "hist", "scale", 
                                "ideogram", "text"), linked.to,
                          radius = 10, trackWidth = 5, 
                          space.skip = 0.015, direction = c("clockwise", "anticlockwise"),
                          link.fun = function(x, y, n = 30) bezier(x, y, evaluation = n),
                   rect.inter.n = 5, rank, 
                   scale.n = 60, scale.unit = NULL, scale.type = c("M", "B", "sci")){

  args <- as.list(match.call(call = sys.call(sys.parent(2)))[-1])
  ## args <- force(args)
  args.aes <- parseArgsForAes(args)
  args.non <- parseArgsForNonAes(args)
  args.non <- args.non[!names(args.non) %in% c("geom", "data", "rank")]

  dots <- as.list(match.call(call = sys.call(sys.parent(2)))[-1])
  scale.type <- match.arg(scale.type)
  geom <- match.arg(geom)
  if(geom == "rect")
    geom <- "rectangle"

  ## rank
  if(!missing(rank)){
    radius <- radius + rank * trackWidth
  }
  ## idoegram parse seqlengths
  if(geom == "ideogram"){
    ## data.back <- data
    data <- getIdeoGR(data)
    res <- rectInter(data, y = as.character(args.aes$y),
                    space.skip = space.skip, trackWidth = trackWidth, radius = radius,
                    direction = direction, n = rect.inter.n)
    
    df <- as.data.frame(res)
    idx <- order(df$.biovizBase.group, df$.int.id)
    df <- df[idx, ]

    ## args.aes.text <- args.aes
    ## args.aes.text$y <- as.name(".biovizBase.y")
    ## args.aes.text$x <- as.name(".biovizBase.x")
    
    args.aes <- args.aes[names(args.aes) != "label"]
    args.aes$y <- as.name(".biovizBase.y")
    args.aes$x <- as.name(".biovizBase.x")
    

    
    args.aes$group <- as.name(".biovizBase.group")

    
    if("fill" %in% names(args.aes)){
      if(!"color" %in% names(args.aes)){
        args.aes$color <- args.aes$fill
      }
    }
    aes <- do.call("aes", args.aes)
    if(!"color" %in% names(args.aes)){
      col <- I("black")
      args.non$color <- col
    }

    args.tot <- c(list(data = df, aes), args.non)
    res <- do.call(geom_polygon, args.tot)
    p <- list(res)
  }

  if(geom == "text"){

    if("label" %in% names(args.aes)){
      lbs <- as.character(args.aes$label)
      if(!lbs %in% c(colnames(values(data)),"start", "end", "seqnames","width"))
        stop("label must be one of column names")
    }else{
      stop("missing label argument in aes()")
    }
    obj <- gr2newLinear(data, space.skip)    
    obj <- gr2circle(obj, y = as.character(args.aes$y), radius= radius, width = trackWidth,
                     direction = direction)
    ## compute angle
    if("angle" %in% names(args.aes)){
      ags <- eval(args.aes$angle, data)
      ags <-  - values(obj)$.biovizBase.angle * 180 / pi + ags
      values(obj)$.processed.angle <- ags
      args.aes$angle <- as.name(".processed.angle")      
    }else{
      ags <-  - values(obj)$.biovizBase.angle * 180/pi 
      values(obj)$.processed.angle <- ags
      args.aes$angle <- as.name(".processed.angle")      
    }

    if("angle" %in% names(dots)){
      ags <-  - values(obj)$.biovizBase.angle * 180 / pi +
        as.numeric(paste(as.character(dots$angle), collapse = ""))
      values(obj)$.processed.angle <- ags
      args.aes$angle <- as.name(".processed.angle")      
    }

    df <- as.data.frame(obj)
    args.aes$y <- as.name(".biovizBase.y")
    args.aes$x <- as.name(".biovizBase.x")
    aes <- do.call("aes", args.aes)
    args.tot <- c(list(data = df, aes), args.non)
    res <- do.call(geom_text, args.tot)
    p <- list(res)
    
  }

  if(geom == "point"){
    obj <- gr2newLinear(data, space.skip)
    if(!"y" %in% names(args.aes)){
      .y <- 1
      warning("y is missing in aes(), use equal y")
    }else{
      .y <- as.character(args.aes$y)
    }
    obj <- gr2circle(obj, y = .y, radius= radius, width = trackWidth,
                     direction = direction)
    df <- as.data.frame(obj)
    args.aes$y <- as.name(".biovizBase.y")
    args.aes$x <- as.name(".biovizBase.x")
    aes <- do.call("aes", args.aes)
    args.tot <- c(list(data = df, aes), args.non)
    res <- do.call(geom_point, args.tot)
    p <- list(res)
  }
  
  if(geom == "line"){
    if(!"y" %in% names(args.aes))
      stop("y is missing in aes()")
    obj <- gr2newLinear(data, space.skip)    
    obj <- gr2circle(obj, y = as.character(args.aes$y), radius= radius, width = trackWidth,
                     direction = direction)
    df <- as.data.frame(obj)
    args.aes$y <- as.name(".biovizBase.y")
    args.aes$x <- as.name(".biovizBase.x")
    args.aes$group <- as.name("seqnames")
    aes <- do.call("aes", args.aes)
    args.tot <- c(list(data = df, aes), args.non)
    res <- do.call(geom_path, args.tot)
    p <- list(res)
  }
  
  if(geom == "segment"){
    ## TODO
    res <- segInter(data, y = as.character(args.aes$y),
                    space.skip = space.skip, trackWidth = trackWidth, radius = radius,
                      direction = direction)
    df <- as.data.frame(res)
    args.aes$y <- as.name(".biovizBase.y")
    args.aes$x <- as.name(".biovizBase.x")
    args.aes$group <- as.name(".biovizBase.group")    
    aes <- do.call("aes", args.aes)
    args.tot <- c(list(data = df, aes), args.non)
    res <- do.call(geom_path, args.tot)
    p <- list(res)

  }
  
  if(geom == "scale"){
    ## like ideogram
    res <- getIdeoGR(data)
    res <- getScale(res, scale.unit, n = scale.n, scale.type)
    values(res)$.biovizBase.group <- seq_len(length(res))
    res0 <- res
    values(res0)$scale.y <- 0
    values(res0)$.biovizBase.group <- seq_len(length(res0))
    res <- c(res, res0)
    res <- gr2newLinear(res, space.skip)    
    res <- gr2circle(res, y = "scale.y", radius= radius, width = trackWidth,
                     direction = direction)
    df <- as.data.frame(res)
    idx <- order(df$.biovizBase.group)
    df <- df[idx, ]
    N <- nrow(df)
    res <- df[seq(1, N-1, by = 2),]
    res[,c(".biovizBase.xend", ".biovizBase.yend")] <-
      df[seq(2, N, by = 2), c(".biovizBase.x", ".biovizBase.y")]
    args.aes$y <- as.name(".biovizBase.y")
    args.aes$x <- as.name(".biovizBase.x")
    args.aes$yend <- as.name(".biovizBase.yend")
    args.aes$xend <- as.name(".biovizBase.xend")
    ## aes <- do.call("aes", args.aes)
    args.aes.text <- args.aes[!names(args.aes) %in% c("xend", "yend")]
    if("angle" %in% names(args.aes)){
      ags <- eval(args.aes$angle, data)
      ags <- 90 - res$.biovizBase.angle * 180 / pi + ags
      res$.processed.angle <- ags
      args.aes.text$angle <- as.name(".processed.angle")      
    }else{
      ags <- 90 - res$.biovizBase.angle * 180/pi 
      res$.processed.angle <- ags
      args.aes.text$angle <- as.name(".processed.angle")      
    }
    args.aes.text$label <- as.name("text.major")
    if(!"hjust" %in% c(names(args.non), names(args.aes.text)))
      args.non$hjust <- 0
    if(!"size" %in% c(names(args.non), names(args.aes.text)))    
      args.non$size <- 3
    
    aes <- do.call("aes", args.aes)
    aes.text <- do.call("aes", c(args.aes.text))
    args.tot <- c(list(data = res), list(aes.text), args.non)
    res.text <- do.call(geom_text, args.tot)
    res.seg <- do.call(ggplot2::geom_segment,c(list(data = res), list(aes)))
    p <- c(list(res.text), list(res.seg))
  }

  if(geom == "rectangle"){
    res <- rectInter(data, y = as.character(args.aes$y),
                    space.skip = space.skip, trackWidth = trackWidth, radius = radius,
                    direction = direction, n = rect.inter.n)
    df <- as.data.frame(res)
    idx <- order(df$.biovizBase.group, df$.int.id)
    df <- df[idx, ]
    args.aes.p <- args.aes
    args.aes.p$y <- as.name(".biovizBase.y")
    args.aes.p$x <- as.name(".biovizBase.x")
    args.aes.p$group <- as.name(".biovizBase.group")

    if("fill" %in% names(args.aes.p)){
      if(!"color" %in% names(args.aes.p)){
        args.aes.p$color <- args.aes.p$fill
      }
    }
    aes.p <- do.call("aes", args.aes.p)
    if(!"color" %in% names(args.aes.p)){
      col <- I("black")
      args.non$color <- col
    }
    args.tot <- c(list(data = df, aes.p), args.non)
    res <- do.call(geom_polygon, args.tot)
    p <- list(res)
  }
  
  if(geom == "bar"){
    res <- barInter(data, y = as.character(args.aes$y),
                    space.skip = space.skip, trackWidth = trackWidth, radius = radius,
                    direction = direction)
    df <- as.data.frame(res)
    idx <- order(df$.biovizBase.group)
    df <- df[idx, ]
    N <- nrow(df)
    res <- df[seq(1, N-1, by = 2),]
    res[,c(".biovizBase.xend", ".biovizBase.yend")] <-
      df[seq(2, N, by = 2), c(".biovizBase.x", ".biovizBase.y")]
    args.aes$y <- as.name(".biovizBase.y")
    args.aes$x <- as.name(".biovizBase.x")
    args.aes$yend <- as.name(".biovizBase.yend")
    args.aes$xend <- as.name(".biovizBase.xend")
    aes <- do.call("aes", args.aes)

    args.tot <- c(list(data = df, aes), args.non)
    res <- do.call(geom_segment, args.tot)
    p <- list(res)
  }
  
  if(geom == "link"){
    res <- linkInter(data, space.skip = space.skip, linked.to = linked.to,
                     link.fun = link.fun, trackWidth = trackWidth, radius = radius,
                     direction = direction)
    args.aes$y <- as.name("y")
    args.aes$x <- as.name("x")
    args.aes$group <- as.name(".biovizBase.group")
    aes <- do.call("aes", args.aes)
    args.tot <- c(list(data = res, aes), args.non )
    res <- do.call(geom_path, args.tot)
    p <- list(res)

  }
  if(geom == "ribbon"){
    stop("geom(ribbon) is not implemented yet")
  }
  p <- c(p, list(opts(aspect.ratio = 1), theme_null()))
  p 
})



