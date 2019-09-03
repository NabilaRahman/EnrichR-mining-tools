function( enrichData=mysummary, desiredGenes = desiredGenes) { 
  
  genes<-unique(unlist(strsplit(as.character(desiredGenes), split = "; ")))
  
  
  for (i in 1:nrow(enrichData)) {
    b<-enrichData[i,"Genes"]
    b<-unique(unlist(strsplit(as.character(b), split = ";")))
    bcount<- length(b)
    desiredGenes <- paste(subset(b, b %in% genes), collapse="; ")
    remainder <- paste(subset(b, !(b %in% genes)), collapse="; ")
    desiredGeneCount <- length(subset(b, b %in% genes))
    enrichData[i,"desiredGenes"] <- desiredGenes
    enrichData[i,"otherGenes"] <- remainder
    enrichData[i,"desiredGeneCount"] <- desiredGeneCount
    enrichData[i,"Remainder"] <- bcount - desiredGeneCount

    
  }
  return ( enrichData )
}


# PlotStackedBars ---------------------------------------------------------
plotStackedBars <- function(data=tf.Down.Gene.Up
                            , refcol = "Term"
                            , title = "C"
                            , color = "steelblue"
                            , transp = 0.4 
                            , x.lab = element_blank()
                            , y.lab = element_blank()
                            , barMax = 15
                            , type="desGenes"
                            , overwrite = T
) {
  require(ggplot2)
  require(reshape)
  plotData<- data.frame( Term=reorder(data[,refcol], -data$Adjusted.P.value), "Desired Genes"=data$desiredGeneCount, Others=data$Remainder,  Pval=-log10(data$Adjusted.P.value) )
  if( nrow(plotData) >= barMax ) {
    plotData <- plotData[c(1:barMax),]
  }
  plotData<-melt(plotData, id=c("Term","Pval") )
  require(dplyr)
  plotData <- plotData %>%
    group_by(Term) %>%
    arrange(Term, desc(variable)) %>%
    mutate(lab_ypos = cumsum(value) - 0.5 * value) 
  
  
  plotEnrichR <- ggplot(data=plotData, aes(x=Term, y=value, fill = variable  ) ) +
    geom_bar(stat = "identity", alpha = transp, width=0.7,color=color, ) + 
    #  geom_bar(fill="steelblue", alpha = 0.6) +
    geom_text(aes(y = lab_ypos, label = value, group =variable), size=3.5 ) +
    coord_flip() +
    theme_classic() +
    theme(
      panel.grid.major = element_blank()
      , panel.grid.minor = element_blank()
      , plot.title = element_text( hjust= 0 )
      , panel.background = element_blank()
      , axis.ticks.y = element_blank()
    ) +
    labs(
      title = title
      , x = x.lab
      , y = y.lab
    ) +
    scale_y_continuous(expand = c(0,0)) +   ## Removes gap  
    scale_fill_manual(values = c(  color, "grey" )) +
    theme(legend.position="bottom", legend.direction="horizontal",
          legend.title = element_blank()) + 
    guides(shape = guide_legend(override.aes = list(size = 0.5)))
  
  filename = paste("barH.des",type,project,"svg",sep=".")
  
  if( !file.exists(filename) | overwrite == T ) {
    ggsave(
      file=filename 
      , plot = plotEnrichR
      , device = "svg"
      , dpi = 300
      , width = 4.35
    )
  }
  
  return (plotEnrichR)
}
