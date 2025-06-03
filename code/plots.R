plot_data_matrix_heatmap <- function(pwm,
                                     data.l,
                                     chip.df,
                                     rank,
                                     cov_name,
                                     data_name,
                                     chip_name,
                                     zMax_data = c(1, 3),
                                     zMax_chip = c(200, 1),
                                     title,
                                     narrow_panel_width = 1,
                                     wide_panel_width = 3.5){

  for ( i in c(1:length(data.l))){
    data.l[[i]][data.l[[i]] > zMax_data[[i]]] = zMax_data[[i]]
  }

  pred.df = chip.df
  for ( i in c(1:ncol(chip.df))){
    chip.df[,i][chip.df[,i] > zMax_chip[[i]]] = zMax_chip[[i]]
  }

  myColors <- colorRampPalette(c(rgb(1,1,1),rgb(0,0,1),rgb(0,0,0.5)))(1000);

  colNumber = 1+length(data.l)+ncol(chip.df)

  layout(matrix(c(1:colNumber), nrow = 1, ncol = colNumber, byrow = F),
         widths=c(rep(narrow_panel_width, 1), rep(wide_panel_width,length(data.l)), rep(narrow_panel_width,ncol(chip.df))),
         heights = rep(1,colNumber))
  cex.text = 0.7
  padj.x = -2.1
  padj.y = 1.8
  par(oma = c(2,2,2,2))

  # cat('plotting pwm \n')
  par(mai=c(0.1,0.2,0.5,0.2))
  image(t(as.matrix(pwm[rank])), col = myColors, axes=FALSE)
  box()
  mtext("PWM\nscore",3,line=0.5,cex =cex.text)

  # cat('plotting data \n')
  for ( i in c(1:length(data.l))){

    par(mai=c(0.1,0.2,0.5,0.2))
    data = data.l[[i]]
    x = c(1:ncol(data))
    y = c(1:nrow(data))
    image(x, y , z = t(data[rank,]), col = myColors, axes=FALSE, zlim = c(0,zMax_data[[i]]))
    box()
    if(ncol(data) < 200){
      axis(1,at=c(1, 51, ncol(data)-50, ncol(data)),labels=c(-50, "","" ,50),padj=padj.x,cex.axis= cex.text,tck=-0.05, tick = T, cex = cex.text)}
    else{
      axis(1,at=c(1, 101, ncol(data)-100, ncol(data)),labels=c(-100, "","" ,100),padj=padj.x,cex.axis= cex.text,tck=-0.05, tick = T, cex = cex.text)
    }
    mtext(paste(data_name[i]),3,line=0.5,cex =cex.text)
    mtext('Dist. to motif (bp)',1,line=1, cex = cex.text-0.1)
  }

  # cat('plotting chip \n')
  for ( i in 1: ncol(chip.df)){
    par(mai=c(0.1,0.2,0.5,0.2))
    image(t(as.matrix(chip.df[rank,i])), col = myColors, axes=FALSE)
    box()
    mtext(chip_name[i],3,line=0.5,cex = cex.text)
  }

  mtext(title, side=3, outer=TRUE, line= -0.5, cex = cex.text + 0.1)

}

plot_pos_neg_profiles <- function(pos_profile, neg_profile,
                                  title = "",
                                  flank = 100,
                                  ylab = "Average counts",
                                  ylim = c(0, max(c(pos_profile, neg_profile), na.rm= TRUE)+0.1),
                                  type = "l"){
  plot(pos_profile, type = type, col = "red", ylab = ylab, ylim = ylim,
       xlab = "Dist. to motif (bp)", xaxt="n", main = title)
  lines(neg_profile, type = type, col = "blue")
  axis(1,at=c(1, flank+1, length(pos_profile)-flank, length(pos_profile)),
       labels=c(-flank, "","" ,flank),
       padj=-0.5,cex.axis= 0.8, tck=-0.02, tick = T)
  legend("topleft", legend = c("positive sites", "negative sites"),
         col = c("red", "blue"), bty = "n", lty = 1)
}

colors_25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)
