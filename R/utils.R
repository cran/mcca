print.mcca.ccp=function(x,...){
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Overall Correct Classification Probability:\n",x$measure,"\n\n")
  cat("Category-specific Correct Classification Probability:\n")
  print(x$table)
}


print.mcca.pdi=function(x,...){
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Overall Polytomous Discrimination Index:\n",x$measure,"\n\n")
  cat("Category-specific Polytomous Discrimination Index:\n")
  print(x$table)
}

print.mcca.rsq=function(x,...){
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Overall R-squared value:\n",x$measure,"\n\n")
  cat("Category-specific R-squared value:\n")
  print(x$table)
}

print.mcca.hum=function(x,...){
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Overall HUM value:\n",x$measure,"\n\n")
}

plot.mcca.hum=function(x,labs=levels(x$y),coords=1:3,nticks=5,filename='fig.png',cex=0.7,...){

  n=100

  y_raw=x$y
  y=as.numeric(y_raw)

  #true positive rate for three classes
  tpP<-function(c1,c2,x) sum(x[,1]>c1 & y==1)/sum(y==1)
  tpN<-function(c1,c2,x) sum(!(x[,1]>c1) & x[,2]>c2 & y==2)/sum(y==2)
  #tpN<-function(c1,c2,x) sum(x[,2]>c2 & y==2)/sum(y==2)
  tpZ<-function(c1,c2,x) sum(!(x[,2]>c2) &! (x[,1]>c1) & y==3)/sum(y==3)

  #The outer product
  qP0<-seq(0,1,length.out=n)
  qN0<-seq(0,1,length.out=n)

  X<-outer(qP0,qN0,Vectorize(tpP,vectorize.args=c("c1","c2")),x$pm)
  Y<-outer(qP0,qN0,Vectorize(tpN,vectorize.args=c("c1","c2")),x$pm)
  Z<-outer(qP0,qN0,Vectorize(tpZ,vectorize.args=c("c1","c2")),x$pm)
  X=1-X
  Z=1-Z

  #visulization
  rgl::rgl.surface(X,Z,Y,coords = coords,color=rainbow(10)[cut(Z, breaks = 10)],
                   back = "fill",front = "fill")
  rgl::rgl.bbox(xlen=0, ylen=0, zlen=0)
  rgl::axes3d(c('x','y','z'),color='white',nticks=nticks,family = "serif",cex = cex)
  rgl::title3d('','',labs[coords][1],labs[coords][2],labs[coords][3],color='white',
               family = "serif",cex = cex)
  rgl::view3d( theta = 210, phi = 10)
  rgl::snapshot3d(filename, fmt = "png", top = TRUE)
  #rgl.viewpoint( theta = 1, phi = 15, fov = 60, zoom = 0, interactive = TRUE )
  #axes <- rbind(c(0.5, 0, 0), c(0, 0.5, 0),
  #              c(0, 0, 0.5))
  #rgl::rgl.texts(axes, text = levels(y_raw)[coords], color = "white",
  #               adj = c(-1, 0), size = 4)
}
