
precision_plot <- function(object){
        data_old <- object$Ref
        data_new <- object$New
        min_y <- min(min(data_old$sig_resid_y2,na.rm=TRUE),min(data_new$sig_resid_y1_corrected,na.rm=TRUE))
        max_y <- max(max(data_old$sig_resid_y2,na.rm=TRUE),max(data_new$sig_resid_y1_corrected,na.rm=TRUE))
        min_y <- floor(min_y)
        range = max_y - min_y
        max_y <- ceiling(max_y+range*0.2)
        par(mar=c(3.5,3.5,2,2)+0.1)
        plot(data_old$y2_hat,data_old$sig_resid_y2,xlab="",
             ylab="",axes = F,
             cex=0.8,ylim=c(min_y,max_y))
        title(main="Precision plot",cex.main=0.9)
        ### Add the y axis
        axis(2,col="black",las=1)
        mtext("Standard deviation of measurement errors",side = 2,line=2,cex=0.9)
        box(col="black")
        ### Add the x axis
        axis(1)
        mtext("BLUP of x",side=1,col="black",line=2,cex=0.9)

        points(data_new$y2_hat,data_new$sig_resid_y1_corrected,
               cex=0.8,pch=19,col="blue")
       legend("topleft",legend=c("Reference method (y2)","New method (corrected y1)"),
               pch=c(1,19),col=c("black","blue"),xpd=T,horiz = F,
               box.lwd=0,cex=0.7)
}
