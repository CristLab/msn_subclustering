#This is the plot_constellation function from Allen Brain's scrattch.hicat with minor changes
#The fill argument is included in the geom_mark_hull function. The scale_fill_identity function is added on the next line
#The function outputs PNG and EPS files rather than PDF

plot_constellation_fill <- function(knn.cl.df, 
                               cl.center.df, 
                               out.dir, 
                               node.label="cluster_id", 
                               exxageration=2, 
                               curved = TRUE, 
                               plot.parts=FALSE, 
                               plot.hull = NULL, 
                               plot.height=25, 
                               plot.width=25, 
                               node.dodge=FALSE, 
                               label.size=5, 
                               max_size=10, 
                               label_repel=FALSE,
                               return.list = FALSE) { 
  
  library(gridExtra)
  library(sna)
  library(Hmisc)
  library(reshape2)
  #library(ggalt)
  library(ggforce)
  library(dplyr)
  
  st=format(Sys.time(), "%Y%m%d_%H%M%S_")
  
  
  if(!file.exists(out.dir)){
    dir.create(out.dir)
  }
  ###==== Cluster nodes will represent both cluster.size (width of point) and edges within cluster (stroke of point)
  
  # select rows that have edges within cluster
  knn.cl.same <- knn.cl.df[knn.cl.df$cl.from == knn.cl.df$cl.to,] 
  
  #append fraction of edges within to cl.center.umap for plotting of fraction as node linewidth
  cl.center.df$edge.frac.within <- knn.cl.same$frac[match(cl.center.df$cl, knn.cl.same$cl.from)] 
  
  
  ###==== plot nodes
  labels <- cl.center.df[[node.label]] 
  
  p.nodes <-   ggplot() +     
                  geom_point(data=cl.center.df,
                             shape=19, 
                             aes(x=x, 
                                 y=y, 
                                 size=cluster_size, 
                                 color=alpha(cluster_color, 0.8))) +
                  scale_size_area(trans="sqrt",
                                  max_size=max_size,
                                  breaks = c(100,1000,10000,100000)) +
                  scale_color_identity() +  
                  geom_text(data=cl.center.df,
                            aes(x=x, 
                                y=y, 
                                label=labels),
                            size = label.size)
  

  #+ theme_void()
  #p.nodes
  if (plot.parts == TRUE) {
    ggsave(file.path(out.dir,paste0(st,"nodes.org.pos.pdf")), p.nodes, width = plot.width, height = plot.height, units="cm",useDingbats=FALSE) }
  
  
  ###==== extract node size/stroke width to replot later without scaling
  g <- ggplot_build(p.nodes)
  dots <-g[["data"]][[1]] #dataframe with geom_point size, color, coords
  
  nodes <- left_join(cl.center.df, dots, by=c("x","y")) %>% ungroup()

  
  ###==== if node.dodge==TRUE new xy coords are calculated for overlapping nodes.
  
    if (node.dodge==TRUE){
  
  #<><><># make update here to convert units by scale. check geom_mark_hull code for oneliner
      
  # dodge nodes starting at center of plot moving outward 
      
      nodes$r<- (nodes$size/10)/2
      
      
      x.list <- c(mean(nodes$x), nodes$x )
      y.list <- c(mean(nodes$y), nodes$y)
      dist.test <- as.matrix(dist(cbind(x.list, y.list)))
      nodes$distance <- dist.test[2:nrow(dist.test), 1]
      nodes <- nodes[order(nodes$distance),]
       
      
  for (d1 in 1:(nrow(nodes)-1)) {
    j <- d1+1
      for (d2 in j:nrow(nodes)) {
      print(paste(d1,d2))
      
        distSq <- sqrt(((nodes$x[d1]-nodes$x[d2])*(nodes$x[d1]-nodes$x[d2]))+((nodes$y[d1]-nodes$y[d2])*(nodes$y[d1]-nodes$y[d2])))
      
        radSumSq <- (nodes$r[d1] *1.25)+ (nodes$r[d2]*1.25) # overlapping radius + a little bit extra
       
        if (distSq < radSumSq) {
          print(paste(d1,d2))
          
          subdfk <- nodes[c(d1,d2),]
          subdfk.mod <- subdfk
          subdfd1 <- subdfk[1,]
          subdfd2  <- subdfk[2,]
          angsk <- seq(0,2*pi,length.out=nrow(subdfd2)+1)
          subdfd2$x <- subdfd2$x+cos(angsk[-length(angsk)])*(subdfd1$r+subdfd2$r+0.5)#/2
          subdfd2$y <- subdfd2$y+sin(angsk[-length(angsk)])*(subdfd1$r+subdfd2$r+0.5)#/2
          subdfk.mod[2,] <- subdfd2
          nodes[c(d1,d2),] <- subdfk.mod
        }
      }
    }
      
      
      for (d1 in 1:(nrow(nodes)-1)) {
        j <- d1+1
        for (d2 in j:nrow(nodes)) {
          print(paste(d1,d2))
          
          distSq <- sqrt(((nodes$x[d1]-nodes$x[d2])*(nodes$x[d1]-nodes$x[d2]))+((nodes$y[d1]-nodes$y[d2])*(nodes$y[d1]-nodes$y[d2])))
          
          radSumSq <- (nodes$r[d1] *1.25)+ (nodes$r[d2]*1.25) # overlapping radius + a little bit extra
          
          if (distSq < radSumSq) {
            print(paste(d1,d2))
            
            subdfk <- nodes[c(d1,d2),]
            subdfk.mod <- subdfk
            subdfd1 <- subdfk[1,]
            subdfd2  <- subdfk[2,]
            angsk <- seq(0,2*pi,length.out=nrow(subdfd2)+1)
            subdfd2$x <- subdfd2$x+cos(angsk[-length(angsk)])*(subdfd1$r+subdfd2$r+0.5)#/2
            subdfd2$y <- subdfd2$y+sin(angsk[-length(angsk)])*(subdfd1$r+subdfd2$r+0.5)#/2
            subdfk.mod[2,] <- subdfd2
            nodes[c(d1,d2),] <- subdfk.mod
          }
        }
      }
      
  }
  
  nodes <- nodes[order(nodes$cluster_id),]
  
  
  
  ## when printing lines to pdf the line width increases slightly. This causes the edge to extend beyond the node. Prevent this by converting from R pixels to points. 
  conv.factor <- ggplot2::.pt*72.27/96
  
  
  ## line width of edge can be scaled to node point size 
  nodes$node.width <- nodes$size 
  
  
  if (plot.parts == TRUE) { 
    if (node.dodge == TRUE) {
    write.csv(nodes, file=file.path(out.dir,paste0(st,"nodes.dodge.csv"))) }
    else {
      write.csv(nodes, file=file.path(out.dir,paste0(st,"nodes.csv")))
    }
  }
  
  ###==== prepare data for plotting of edges between nodes
  
  ##filter out all edges that are <5% of total for that cluster
  #knn.cl <- knn.cl.df[knn.cl.df$frac >0.05,] #1337 lines
  knn.cl <- knn.cl.df
  ##from knn.cl data frame remove all entries within cluster edges.
  knn.cl.d <- knn.cl[!(knn.cl$cl.from == knn.cl$cl.to),] 
  nodes$cl=as.numeric(as.character(nodes$cl))
  knn.cl.d$cl.from <- as.numeric(as.character(knn.cl.d$cl.from))
  knn.cl.d$cl.to <- as.numeric(as.character(knn.cl.d$cl.to))
  
  knn.cl.d <- left_join(knn.cl.d, select(nodes, cl, node.width), by=c("cl.from"="cl"))
  colnames(knn.cl.d)[colnames(knn.cl.d)=="node.width"]<- "node.pt.from"
  knn.cl.d$node.pt.to <- ""
  knn.cl.d$Freq.to <- ""
  knn.cl.d$frac.to <- ""
  
  
  #bidirectional 
  knn.cl.bid <- NULL
  for (i in 1:nrow(knn.cl.d)) {
    
    line <- subset(knn.cl.d[i,])
    r <- subset(knn.cl.d[i:nrow(knn.cl.d),])
    r <- r[(line$cl.from == r$cl.to & line$cl.to == r$cl.from ),] 
    
    if (dim(r)[1] != 0) {
      line$Freq.to <- r$Freq
      line$node.pt.to <- r$node.pt.from
      line$frac.to <- r$frac
      knn.cl.bid <- rbind(knn.cl.bid, line)
    }
    #print(i)
  }
  
  #unidirectional
  knn.cl.uni <- NULL
  for (i in 1:nrow(knn.cl.d)) {
    
    line <- subset(knn.cl.d[i,])
    r <- knn.cl.d[(line$cl.from == knn.cl.d$cl.to & line$cl.to == knn.cl.d$cl.from ),] 
    
    if (dim(r)[1] == 0) {
      knn.cl.uni <- rbind(knn.cl.uni, line)
    }
   #print(i)
  }
  
  
  #min frac value = 0.01
  knn.cl.uni$node.pt.to <- nodes$node.width[match(knn.cl.uni$cl.to, nodes$cl)]
  knn.cl.uni$Freq.to <- 1
  knn.cl.uni$frac.to <- 0.01
  knn.cl.lines <- rbind(knn.cl.bid, knn.cl.uni)
  
  
  ###==== create line segments
  
  line.segments <- knn.cl.lines %>% select(cl.from, cl.to)
  nodes$cl <- as.numeric((as.character(nodes$cl)))
  line.segments <- left_join(line.segments,select(nodes, x, y, cl), by=c("cl.from"="cl"))
  line.segments <- left_join(line.segments,select(nodes, x, y, cl), by=c("cl.to"="cl"))
  colnames(line.segments) <- c("cl.from", "cl.to", "x.from", "y.from", "x.to", "y.to")
  
  line.segments <- data.frame(line.segments,
                              freq.from = knn.cl.lines$Freq,
                              freq.to = knn.cl.lines$Freq.to,
                              frac.from = knn.cl.lines$frac,
                              frac.to =  knn.cl.lines$frac.to,
                              node.pt.from =  knn.cl.lines$node.pt.from,
                              node.pt.to = knn.cl.lines$node.pt.to)
  
  
  ##from points to native coords
  line.segments$node.size.from <- line.segments$node.pt.from/10
  line.segments$node.size.to <- line.segments$node.pt.to/10
  
  
  line.segments$line.width.from <- line.segments$node.size.from*line.segments$frac.from
  line.segments$line.width.to <- line.segments$node.size.to*line.segments$frac.to
  
  ##max fraction to max point size 
  line.segments$line.width.from<- (line.segments$frac.from/max(line.segments$frac.from, line.segments$frac.to))*line.segments$node.size.from
  
  line.segments$line.width.to<- (line.segments$frac.to/max(line.segments$frac.from, line.segments$frac.to))*line.segments$node.size.to
  
  
  ###=== create edges, exaggerated width
  
  line.segments$ex.line.from <-line.segments$line.width.from #true to frac
  line.segments$ex.line.to <-line.segments$line.width.to #true to frac
  
  line.segments$ex.line.from <- pmin((line.segments$line.width.from*exxageration),line.segments$node.size.from) #exxagerated width
  line.segments$ex.line.to <- pmin((line.segments$line.width.to*exxageration),line.segments$node.size.to) #exxagerated width
  

  line.segments <- na.omit(line.segments)
  
  print("calculating edges")
  
  allEdges <- lapply(1:nrow(line.segments), edgeMaker, len = 500, curved = curved, line.segments=line.segments)
  allEdges <- do.call(rbind, allEdges)  # a fine-grained path with bend
  
  
  groups <- unique(allEdges$Group)
  
  poly.Edges <- data.frame(x=numeric(), y=numeric(), Group=character(),stringsAsFactors=FALSE)
  imax <- as.numeric(length(groups))
  
  for(i in 1:imax) { 
    #svMisc::progress(i)
    #svMisc::progress(i, progress.bar=TRUE)
    select.group <- groups[i]
    #print(select.group)
    select.edge <- allEdges[allEdges$Group %in% select.group,]
    
    x <- select.edge$x
    y <- select.edge$y
    w <- select.edge$fraction
    
    N <- length(x)
    leftx <- numeric(N)
    lefty <- numeric(N)
    rightx <- numeric(N)
    righty <- numeric(N)
    
    ## Start point
    perps <- perpStart(x[1:2], y[1:2], w[1]/2)
    leftx[1] <- perps[1, 1]
    lefty[1] <- perps[1, 2]
    rightx[1] <- perps[2, 1]
    righty[1] <- perps[2, 2]
    
    ### mid points
    for (ii in 2:(N - 1)) {
      seq <- (ii - 1):(ii + 1)
      perps <- perpMid(as.numeric(x[seq]), as.numeric(y[seq]), w[ii]/2)
      leftx[ii] <- perps[1, 1]
      lefty[ii] <- perps[1, 2]
      rightx[ii] <- perps[2, 1]
      righty[ii] <- perps[2, 2]
    }
    ## Last control point
    perps <- perpEnd(x[(N-1):N], y[(N-1):N], w[N]/2)
    leftx[N] <- perps[1, 1]
    lefty[N] <- perps[1, 2]
    rightx[N] <- perps[2, 1]
    righty[N] <- perps[2, 2]
    
    lineleft <- data.frame(x=leftx, y=lefty)
    lineright <- data.frame(x=rightx, y=righty)
    lineright <- lineright[nrow(lineright):1, ]
    lines.lr <- rbind(lineleft, lineright)
    lines.lr$Group <- select.group
    
    poly.Edges <- rbind(poly.Edges,lines.lr)
    
    Sys.sleep(0.01)
    cat("\r", i, "of", imax)
    
  }
  
  if (plot.parts == TRUE) {
    write.csv(poly.Edges, file=file.path(out.dir,paste0(st,"poly.edges.csv"))) }
 
  
  #############################
  ##                         ##
  ##        plotting         ##
  ##                         ##
  #############################
  
  labels <- nodes[[node.label]] 
  
   
  ####plot edges
  p.edges <- ggplot(poly.Edges, aes(group=Group))
  p.edges <- p.edges +geom_polygon(aes(x=x, y=y), alpha=0.2) + theme_void()
  #p.edges
  
  if (!is.null(plot.hull)) {
  #### plot all layers
    plot.all <-  ggplot()+
          geom_polygon(data=poly.Edges, 
                        alpha=0.2, 
                        aes(x=x, y=y, group=Group))+ 
          geom_point(data=nodes,
                     alpha=0.8, 
                     shape=19,
                     aes(x=x, 
                         y=y, 
                         size=cluster_size, 
                         color=cluster_color)) +
          scale_size_area(trans="sqrt",
                          max_size=max_size,
                          breaks = c(100,1000,10000,100000)) +
          scale_color_identity()  + 
          theme_void()+ 
          geom_mark_hull(data=nodes,
                         concavity = 8,
                         radius = unit(5,"mm"),
                         aes(filter = nodes$clade_id %in% plot.hull,x, y, 
                             color=nodes$clade_color,fill=nodes$clade_color)) +
          scale_fill_identity() +
		  theme(legend.position = "none")
    
      if(label_repel ==TRUE){
        plot.all <- plot.all +
          ggrepel::geom_text_repel(data=nodes,
                          aes(x=x, 
                              y=y, 
                              label=labels),
                          size = label.size,
                          min.segment.length = Inf) 
      }
      else{
        plot.all <- plot.all +
            geom_text(data=nodes,
                      aes(x=x, 
                          y=y, 
                          label=labels),
                      size = label.size)  }

  #plot.all
    } else {
    #### plot all layers
    plot.all <-  ggplot()+
          geom_polygon(data=poly.Edges, 
                       alpha=0.2, 
                       aes(x=x, y=y, group=Group))+ 
          geom_point(data=nodes,
                     alpha=0.8, 
                     shape=19,
                     aes(x=x, 
                         y=y, 
                         size=cluster_size, 
                         color=cluster_color)) +
          scale_size_area(trans="sqrt",
                          max_size=max_size,
                          breaks = c(100,1000,10000,100000)) +
          scale_color_identity() 
      if(label_repel ==TRUE){
        plot.all <- plot.all +
                  geom_text_repel(data=nodes,
                        aes(x=x, 
                            y=y, 
                            label=labels),
                        size = label.size,
                        min.segment.length = Inf) + 
                  theme_void() +
                  theme(legend.position="none")
      }
        else{
          plot.all <- plot.all +
                  geom_text(data=nodes,
                            aes(x=x, 
                                y=y, 
                                label=labels),
                            size = label.size) + 
                  theme_void() +
                  theme(legend.position="none") }
    #plot.all
  }
  
  
  
  segment.color = NA
  
  if (plot.parts == TRUE) {
    ggsave(file.path(out.dir,paste0(st,"comb.constellation.pdf")), plot.all, width = plot.width, height = plot.height, units="cm",useDingbats=FALSE) }
  
  
  
  #############################
  ##                         ##
  ##      plot legends       ##
  ##                         ##
  #############################
  
  
  ### plot node size legend (1)
  plot.dot.legend <- ggplot()+
              geom_polygon(data=poly.Edges, 
                           alpha=0.2, 
                           aes(x=x, y=y, group=Group))+ 
              geom_point(data=nodes,
                         alpha=0.8, 
                         shape=19,
                         aes(x=x, 
                             y=y, 
                             size=cluster_size, 
                             color=cluster_color)) +
              scale_size_area(trans="sqrt",
                              max_size=max_size,
                              breaks = c(100,1000,10000,100000)) +
              scale_color_identity() + 
              geom_text(data=nodes,
                        aes(x=x, 
                            y=y, 
                            label=labels),
                        size = label.size)+
              theme_void()
  dot.size.legend <- cowplot::get_legend(plot.dot.legend)
  
  ### plot cluster legend (3)
  cl.center.df$cluster.label <-  cl.center.df$cluster_label
  cl.center.df$cluster.label <- as.factor(cl.center.df$cluster.label)
  label.col <- setNames(cl.center.df$cluster_color, cl.center.df$cluster.label)
  cl.center.df$cluster.label <- as.factor(cl.center.df$cluster.label)
  leg.col.nr <- min((ceiling(length(cl.center.df$cluster_id)/20)),5)
   
  cl.center <- ggplot(cl.center.df, 
                      aes(x=cluster_id, y=cluster_size)) + 
              geom_point(aes(color=cluster.label))+
              scale_color_manual(values=as.vector(label.col[levels(cl.center.df$cluster.label)]))+  
              guides(color = guide_legend(override.aes = list(size = 8), ncol=leg.col.nr))
  
  cl.center.legend <- cowplot::get_legend(cl.center)  
  #plot(cl.center.legend)
  
  
  ###plot legend line width (2)
  width.1 <- max(line.segments$frac.from,line.segments$frac.to) 
  width.05 <- width.1/2
  width.025 <- width.1/4
  
  
  edge.width.data <- tibble(node.width = c(1,1,1), x=c(2,2,2), y=c(5,3.5,2), line.width=c(1,0.5,0.25), fraction=c(100, 50, 25),frac.ex=c(width.1, width.05, width.025))
  edge.width.data$fraction.ex <- round((edge.width.data$frac.ex*100), digits = 0)

  poly.positions <- data.frame(id=rep(c(1,2,3), each = 4), x=c(1,1,2,2,1,1,2,2,1,1,2,2), y=c(4.9,5.1,5.5,4.5,3.4,3.6,3.75,3.25,1.9,2.1,2.125,1.875)) 
  
if (exxageration !=1) {
 edge.width.legend <- ggplot()  +  
        geom_polygon(data=poly.positions, aes(x=x,y=y, group=id), fill="grey60")+
        geom_circle(data=edge.width.data, aes(x0=x, y0=y, r=node.width/2), fill="grey80", color="grey80", alpha=0.4)+ 
        scale_x_continuous(limits=c(0,3)) + 
        theme_void() +
        coord_fixed() + 
        geom_text(data=edge.width.data, aes(x= 2.7, y=y, label=fraction.ex, hjust=0, vjust=0.5)) + 
        annotate("text", x = 2, y = 6, label = "Fraction of edges \n to node") } 
    else { edge.width.legend <- ggplot()  +  
        geom_polygon(data=poly.positions, aes(x=x,y=y, group=id), fill="grey60")+
        geom_circle(data=edge.width.data, aes(x0=x, y0=y, r=node.width/2), fill="grey80", color="grey80", alpha=0.4)+ 
        scale_x_continuous(limits=c(0,3)) + 
        theme_void() +coord_fixed() + 
        geom_text(data=edge.width.data, aes(x= 2.7, y=y, label=fraction, hjust=0, vjust=0.5)) + 
        annotate("text", x = 2, y = 6, label = "Fraction of edges \n to node")}
    

 #############################
 ##                         ##
 ##      save elements      ##
 ##                         ##
 #############################
 
 
  
  layout_legend <- rbind(c(1,3,3,3,3),c(2,3,3,3,3))  
  
  if (plot.parts == TRUE) {
    ggsave(file.path(out.dir,paste0(st,"comb.LEGEND.pdf")),gridExtra::marrangeGrob(list(dot.size.legend,edge.width.legend,cl.center.legend),nrow = 3, ncol=6, layout_matrix=layout_legend),height=20,width=20,useDingbats=FALSE)  }
  
  
  g2 <- gridExtra::arrangeGrob(grobs=list(dot.size.legend,edge.width.legend,cl.center.legend), layout_matrix=layout_legend)
  
  
  ggsave(file.path(out.dir,paste0(st,"constellation.pdf")),marrangeGrob(list(plot.all,g2),nrow = 1, ncol=1),width = plot.width, height = plot.height, units="cm",useDingbats=FALSE)
  ggsave(file.path(out.dir,paste0(st,"constellation.png")),plot.all, width = plot.width, height = plot.height, units="cm")
  ggsave(file.path(out.dir,paste0(st,"constellation.eps")),marrangeGrob(list(plot.all,g2),nrow = 1, ncol=1),width = plot.width, height = plot.height, units="cm")
  
  if(return.list==TRUE){
  return(list(constellation = plot.all, edges.df = poly.Edges, nodes.df = nodes))}
}






