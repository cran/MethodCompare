
bias_plot <- function(object){
        ### extract the objects from the output
        Bias <- object$Bias
        data_old <- object$Ref
        data_new <- object$New
        Models <- object$Models
        subtitle <- paste("Differential bias:",round(Bias[1,1],3),"; ",
                          "Proportional bias:",round(Bias[2,1],3),sep="")
        min_y <- min(min(data_old$y2,na.rm=TRUE),min(data_new$y1,na.rm=TRUE))
        max_y <- max(max(data_old$y2,na.rm=TRUE),max(data_new$y1,na.rm=TRUE))
        min_y <- floor(min_y)
        range = max_y-min_y
        max_y <- ceiling(max_y+range*0.2)
        par(mar=c(3.5,3.5,3,4)+0.1)
        ### plot scatter for y1 and y2 with respect to BLUP of x
        plot(data_old$y2_hat,data_old$y2,axes=F,xlab="",ylab="",
             col="grey",ylim = c(min_y,max_y),cex=0.5)
        title(main="Bias",cex.main=0.9)
        points(data_new$y2_hat,data_new$y1,col="blue",pch=19,cex=0.5)
        ### Add the fitted line from model 2 and model 4,
        abline(Models[[2]]$coefficients,lwd=2)
        abline(Models[[4]]$coefficients,lwd=2,lty=2,col="blue")
        ### Add the subtitle
        mtext(subtitle,side=3,cex=0.8)
        ### Add the left y axis
        axis(2,col="black",las=1)
        mtext("y1 and y2",side = 2,line=2)
        box(col="black")

        ### Add second plot: bias axis
        par(new=TRUE)
        ## Plot the bias plot and put axis scale on right
        plot(data_new$y2_hat,data_new$Bias,xlab="",ylab="",axes=FALSE,
             col="red",lty=1,type="l")
        abline(h=0,col="black",lwd=2)
        ## Add the right y axis and label
        mtext("Bias",side=4,col="black",line=2.5)
        axis(4,col="black",col.axis="black",las=1)
        ### Draw the x axis and add the label
        axis(1)
        mtext("BLUP of x",side=1,col="black",line=2)
        ### Add the legend
        lm_bias <- lm(data_new$Bias ~data_new$y2_hat)
        if (coef(lm_bias)[2] <0){
        legend("top",legend=c("Reference method (y2)","New method (y1)","Bias"),
               pch = c(1,19),col = c("black","blue","red"),pt.cex=c(1,1,0.0001),
               y.intersp = 0.7,yjust=0.2,lty=c(1,2,1),bty = "n",cex=0.8)}
        else {legend("topleft",legend=c("Reference method (y2)","New method (y1)","Bias"),
                     pch = c(1,19),col = c("black","blue","red"),pt.cex=c(1,1,0.0001),
                     y.intersp = 0.7,yjust=0.2,lty=c(1,2,1),bty = "n",cex=0.8)}
}
