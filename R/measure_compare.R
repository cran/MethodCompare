
measure_compare <- function(data,new="y1",Ref="y2",ID="id"){
  #### Check whether required package is installed.
  data_sub <-(data[,c(ID,new,Ref)])
  colnames(data_sub) <- c("id","y1","y2")
  #### Calculate average value for reference method and categorize it into 10
  #### or 5 categories
  Mean_y2 <- aggregate(y2~id,data=data_sub,mean)
  N_cat <- ifelse(dim(Mean_y2)[1]>=100,11,6)
  Mean_y2$cat_y2_mean <- cut(Mean_y2$y2,
                             breaks=quantile(Mean_y2$y2, probs=seq(0,1,length.out = N_cat)),
                             include.lowest=TRUE)
  colnames(Mean_y2)[2] <- "y2_mean"
  data_sub <- merge(data_sub,Mean_y2,by="id")
  ### create a dataset for reference method by excluding missing values
  data_sub_y2 <- data_sub[ , c("y2")]
  data_sub_nomiss_y2 <- data_sub[complete.cases(data_sub_y2), ]
  ### Model 1: mixed model for BLUP estimate of x
  model_1 <- lme(y2 ~ 1, data = data_sub_nomiss_y2, random = ~ 1|id,
                 weights=varIdent(form = ~1|cat_y2_mean), na.action = na.exclude)
    data_sub_nomiss_y2$y2_hat <- fitted(model_1)
  ### create a dataset with y2_hat
  data_sub_y2_hat <- data_sub_nomiss_y2[ , c("id","y2_hat")]
  data_sub_y2_hat <- aggregate(y2_hat~id,data=data_sub_y2_hat,mean)
  data_sub <- merge(data_sub,data_sub_y2_hat,by="id")
  ### create a dataset for new method by excluding missing value of new method y1
  data_sub_y1 <- data_sub[ , c("y1")]
  data_sub_nomiss_y1 <- data_sub[complete.cases(data_sub_y1), ]
  data_newmethod <- data_sub_nomiss_y1
  ########         Models on the reference method
  #### Model 2: regression of y2 based on BLUP of x
  model_2 <- lm(y2~y2_hat,data=data_sub, na.action = na.exclude)
  data_sub$fitted_y2 <- data_sub$y2_hat
  data_sub$resid_y2 <- data_sub$y2-data_sub$y2_hat
  data_sub$resid_y2_abs <- abs(data_sub$resid_y2)
  #### Model 3: Estimation of variance function for y2
  model_3 <- lm(resid_y2_abs~y2_hat,data=data_sub, na.action = na.exclude)
  ### Smooth standard deviation estimate
  data_sub$sig_resid_y2 <- fitted(model_3)*sqrt(pi/2)
  #######         Models on the new method - original measured value
  #### Model 4: Regression of y1 based on BLUP of x
  model_4 <- lm(y1~y2_hat,data=data_newmethod)
  #### Differential and proportional bias of new method
  Bias <- cbind(model_4$coefficients,confint(model_4))
  rownames(Bias) <- c("Differential bias","Proportional bias")
  colnames(Bias)[1] <- "Estimate"
  #### Residuals and corrected y1
  data_newmethod$resid_y1 <- residuals(model_4)
  data_newmethod$resid_y1_abs <- abs(data_newmethod$resid_y1)
  data_newmethod$y1_corrected <- (data_newmethod$y1-Bias[1])/Bias[2]
  data_newmethod$Bias <- Bias[1]+data_newmethod$y2_hat*(Bias[2]-1)
  #### Model 5: Variance estimation for y1 based on BLUP of x
  model_5<- lm(resid_y1_abs ~ y2_hat,data=data_newmethod, na.action = na.exclude)
  data_newmethod$fitted_y1 <- fitted(model_5)
  ### Smooth standard deviation estimate
  data_newmethod$sig_resid_y1 <- fitted(model_5)*sqrt(pi/2)
  ########        Models on new method- corrected value
  #### Model 6: regression on corrected y1 using BLUP_y2
  model_6 <- lm(y1_corrected~y2_hat,data=data_newmethod, na.action = na.exclude)
  data_newmethod$y1_corrected_fitted <- fitted(model_6)
  data_newmethod$resid_y1_corrected <- residuals(model_6)
  data_newmethod$resid_y1_corrected_abs <- abs(data_newmethod$resid_y1_corrected)
  #### Model 7: variance function estimation for corrected y1
  model_7 <- lm(resid_y1_corrected_abs ~ y2_hat,data=data_newmethod, na.action = na.exclude)
  data_newmethod$sig_resid_y1_corrected <- fitted(model_7)*sqrt(pi/2)
  ######  Final output
  ## 1. Bias (differential and proportional) from New method
  ## 2. Output values including y2_hat, sig_resid_y2, sig_resid_y1_corr, Bias,
  ## 3. list of all models: mixed model
  Model_list <-  list(model_1,model_2, model_3,model_4,model_5,
                      model_6,model_7)
  data_sub <- data_sub[,c("id","y1","y2","y2_hat","sig_resid_y2")]
  data_newmethod <- data_newmethod[,c("id","y1","y2_hat","Bias","y1_corrected",
                                      "sig_resid_y1_corrected")]
  output <- list(Bias=Bias,Models=Model_list,Ref=data_sub,New=data_newmethod)
  return(output)
}
