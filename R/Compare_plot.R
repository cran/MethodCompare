#' Plot used to visualize the recalibration of the new method after estimating the bias
#'
#' This function allows the visualization of the bias-corrected values (i.e.
#' recalibrated values, variable y1_corr) of the new measurement method.
#'
#' @param object an object retunred by a call to \link{measure_compare}
#' @author  Mingkai Peng
#' @export
#' @examples
#' ### load the data
#' data(data1)
#' ### analysis
#' measure_model <- measure_compare(data1)
#' ### compare plot
#' compare_plot(measure_model)
#'

compare_plot <- function(object){
        data_old <- object$Ref
        data_new <- object$New
        Models <- object$Models
        max <- max(data_old$y2,data_new$y1,data_new$y1_corrected,na.rm=TRUE)
        min <- min(data_old$y2,data_new$y1,data_new$y1_corrected,na.rm=TRUE)
        range <- max-min
        par(mar=c(3.5,3.5,2,2)+0.1)
        plot(data_old$y2_hat,data_old$y2,pch=1,cex=0.5,col="grey",axes = F,
             xlab="",ylab="",ylim=c(min-range*0.1,max+range*0.2))
        title(main="Comparison of the methods",cex.main=0.9)
        xlab="Average:(y1+y2)/2"
        ### Add the y axis
        axis(2,col="black",las=1)
        mtext("Measurement method",side = 2,line=2)
        box(col="black")
        ### Add the x axis
        axis(1)
        mtext("BLUP of x",side=1,col="black",line=2)
        points(data_new$y2_hat,data_new$y1,pch=19,col="blue",cex=0.5)
        points(data_new$y2_hat,data_new$y1_corrected,pch=18,col="red",cex=0.5)
        #abline(Models[[2]]$coefficients,lwd=2)
        abline(c(0,1),lwd=2)
        abline(Models[[4]]$coefficients,lwd=2,lty=1,col="blue")
        abline(Models[[6]]$coefficients,lwd=2,lty=2,col="red")
        legend("topleft",legend=c("Reference method (y2)","New method (y1)",
                                  "New method(corrected)"),
               pch=c(1,19,18),lty = c(1,1,2),col=c("black","blue","red"),
               y.intersp = 0.7,yjust=0.2,bty = "n",cex=0.8)
}
