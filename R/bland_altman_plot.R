#' Bland and Atlman's limits of agreement plot
#'
#' This function produces Bland and Altman's Limits of Agreement plot
#' (LoA) when there are repeated measurements with possibly heteroscedastic
#' measurement errors.
#'
#'
#' @param data a dataframe contains the object identification number (id),
#' the measurement values from new measurement method (y1) and those from
#' the reference standard (y2)
#' @param new specify the variable name or location for the new measurement method
#' @param Ref specify the variable name or location for the reference  measuerment method
#' @param ID specify the variable name for location for the subject identification number (id)
#' @param fill logical. if \code{TRUE} use the avarage value for new methods to
#' fill out the missing value (only useful for drawing a plot with all the
#' measurements by the reference standard)
#'
#' @author  Mingkai Peng
#' @details  This functions computes the limits of agreement (LoA) when there are
#' repeated measurements and possibly the measurement error are heteroscedastic
#' @export
#' @importFrom stats lm
#' @examples
#' ### Load the data
#' data(data1)
#' ### Bland and Altman's plot
#' bland_altman_plot(data1)



bland_altman_plot <- function(data,new="y1",Ref="y2",ID="id",fill=TRUE){
        data_sub <-data[,c(ID,new,Ref)]
        colnames(data_sub) <- c("id","y1","y2")
        #### calculate the difference and average
        data_sub$Diff_M <- data_sub$y1-data_sub$y2
        data_sub$AVG_M <- (data_sub$y1+data_sub$y2)/2
        #### Calculate the harmonic means
        y1 <- as.numeric(data.frame(table(na.omit(data_sub)$id))$Freq)
        y2 <- as.numeric(data.frame(table(data_sub$id))$Freq)
        coef_1<- 1/mean(1/y1)
        coef_1 <- 1-1/coef_1
        coef_2<- 1/mean(1/y2)
        coef_2 <- 1- 1/coef_2
        #### With-object variance
        model_y2 <- lme(y2~1,data=data_sub,random = ~1|id)
        VIF<- model_y2$sigma^2*coef_2
        if (max(y1) > 1){
                model_y1 <- lme(y1~1,data=data_sub,random = ~1|id,na.action = na.exclude)
                VIF<- model_y2$sigma^2*coef_2+model_y1$sigma^2*coef_1
        }
        ####    Model on difference based on average
        data_sub_1 <- na.omit(data_sub)
        Model_1<- lm(Diff_M~AVG_M,data=data_sub_1)
        data_sub_1$fitted <- fitted(Model_1)
        data_sub_1$resid_fitted <- residuals(Model_1)
        data_sub_1$resid_fitted_abs <- abs(data_sub_1$resid_fitted)
        #####   Regression on absoulte residuals based on averge
        Model_2 <- lm(resid_fitted_abs~AVG_M,data=data_sub_1)
        data_sub_1$sig2_abs_res=fitted(Model_2)*sqrt(pi/2)
        data_sub_1$sig2_abs_res <- data_sub_1$sig2_abs_res^2
        data_sub_1$upper <- data_sub_1$fitted+1.96*sqrt(data_sub_1$sig2_abs_res+VIF)
        data_sub_1$lower <- data_sub_1$fitted-1.96*sqrt(data_sub_1$sig2_abs_res+VIF)
        ####
        if (fill) {
                mean_y1 <- aggregate(y1~id,data=data_sub,mean)
                colnames(mean_y1)[2] <- "y1.mean"
                data_sub <-merge(data_sub,mean_y1,by="id")
                data_sub$y1 <- ifelse(is.na(data_sub$y1),data_sub$y1.mean,data_sub$y1)
                data_sub$Diff_M <- data_sub$y1-data_sub$y2
                data_sub$AVG_M <- (data_sub$y1+data_sub$y2)/2
        }
        ### model for plot
        Model_3 <- lm(upper~AVG_M,data=data_sub_1)
        Model_4 <- lm(lower~AVG_M,data=data_sub_1)
        max <- max(abs(data_sub$Diff_M),na.rm = T)
        max <- max + max*0.2
        ##### final plot
        par(mar=c(3.5,3.5,2,2)+0.1)
        plot(data_sub$AVG_M,data_sub$Diff_M,xlab = "",ylab = "",axes = F,
             col="grey",ylim=c(-max,max))
        title(main="Bland and Altman' limit of agreement (LoA)",cex.main=0.9)
        xlab="Average:(y1+y2)/2"
        ### Add the y axis
        axis(2,col="black",las=1)
        mtext("Difference:y1-y2",side = 2,line=2)
        box(col="black")
        ### Add the x axis
        axis(1)
        mtext("Average:(y1+y2)/2",side=1,col="black",line=2)


        abline(Model_1$coefficients,col="blue",lwd=2)
        abline(h=0,col="black",lwd=2)
        abline(Model_3$coefficients,col="blue",lty=2,lwd=2)
        abline(Model_4$coefficients,col="blue",lty=2,lwd=2)
        legend("topleft",legend=c("Regression line","95% LoA"),col = c("blue","blue"),
               y.intersp = 0.7,yjust=0.2,lty=c(2,2),bty = "n",cex=0.8)
}


