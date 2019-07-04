

library(MASS)
library(dplyr)
library(ggpubr)
library(pROC)


## repeat Chi-Squared Test for Lipase and Biulibrubin elevation vs. Cmin Nilotinib 
## Ref: Eur J Clin Pharmaocl (2013) 69:813-823

## Grade 3/4 Lipase
df = data.frame(withAE=c(17,16,18,27), withoutAE=c((122-17),(114-16), (118-18), (122-27) ))
row.names(df) = c("Q1", "Q2", "Q3", "Q4") # Quantiles of Cmin Q1 < 429 Q2 429-615 Q3 615-8550 Q4 > 850 ng/mL

chisq.test(df)



##AllGrade Lipase
df = data.frame(withAE=c(38,39,61,63), withoutAE=c((116-38),(115-39), (121-61), (122-63) ))
row.names(df) = c("Q1", "Q2", "Q3", "Q4")


chisq.test(df)


## Grade 3/4 Bilirubin
df = data.frame(withAE=c(7,12,12,17), withoutAE=c((123-7),(123-12), (123-12), (124-17) ))
row.names(df) = c("Q1", "Q2", "Q3", "Q4")


chisq.test(df)



##AllGrade Bilirubin
df = data.frame(withAE=c(55,97,89,110), withoutAE=c((107-55),(133-97), (123-89), (126-110) ))
row.names(df) = c("Q1", "Q2", "Q3", "Q4")


chisq.test(df)

#### -------------------------------------------------------------------------------------------------


## Wilcoxon rank sum Test cthrough Dabrafenib/Trametinib and dose reduction event
## Ref: Clinica Chimica Acta 472 (2017) 26-29
set.seed(8)

TotalPatients = 27
Pat_without_reduction = 19
mean_cthroug_no_red = 33.5 #ng/mL
sd_cthrough_no_red = 18.5 #ng/mL

Pat_with_reduction = 8
mean_cthroug_with_red = 118.6 #ng/mL
sd_cthrough_with_red = 84.7 #ng/mL

## Simulate 19 Patients with typical c_through in the no reduction group
## Use log normal distribution (information in paper does not support normal dist)
c_trough_no_red <-log( rlnorm(n=Pat_without_reduction, 
                         mean=mean_cthroug_no_red, 
                         sd=sd_cthrough_no_red))

## Simulate 8 Patients with typical c_through in the reduction group
## Use log normal distribution (information in paper does not support normal dist)
c_trough_with_red <-log( rlnorm(n=Pat_with_reduction, 
                           mean=mean_cthroug_with_red, 
                           sd=sd_cthrough_with_red))

df <- data.frame(group=c(rep("no_red", Pat_without_reduction), 
                         rep("red", Pat_with_reduction)),
                 c_through=c(c_trough_no_red, c_trough_with_red))
df

group_by(df, group) %>% summarise(count=n(),
                                  median = median(c_through, na.rm=TRUE),
                                  IQR = IQR(c_through, na.rm=TRUE)
                                  )

ggboxplot(df, x="group", y="c_through")

res <- wilcox.test(c_trough_no_red, c_trough_with_red)
res

### ------ Create a ROC curve with pROC

### Explore the data
#c_through threshold for prediction

thresh = 0 #ng/mL for discrimination

####
##-----Test differenct thresholds to generate ROC curve and Youden Index, highest Youden Index => Threshold for discrimination fo patients
####

roc_data = NULL

TP = NULL
FP = NULL
YI= NULL
TH = NULL

for (thresh in 1:500){
## These patients exceed threshold -> dose reduction
assigned_dose_red <- subset.data.frame(df, c_through > thresh)

## These patients do not exceed the threshold and have no risk for dose reduction
assigned_no_dose_red <- subset.data.frame(df, c_through <= thresh)

assigned_dose_red

true_positive <- subset.data.frame(assigned_dose_red, group == "red")
false_positive <- subset.data.frame(assigned_dose_red, group == "no_red")

true_negative <- subset.data.frame(assigned_no_dose_red, group == "no_red")
false_negative <- subset.data.frame(assigned_no_dose_red, group == "red")

TP = c(TP, length(true_positive$group)/Pat_with_reduction)
FP = c(FP, length(false_positive$group)/Pat_without_reduction)

YI = c(YI, (length(true_positive$group)/(length(true_positive$group)+length(false_negative$group))+length(true_negative$group)/(length(true_negative$group)+length(false_positive$group)) )-1 )
TH = c(TH, thresh)

}
roc_data = data.frame(TP, FP, YI, TH)

max_YI <- subset(roc_data, roc_data$YI == max(roc_data$YI))

threshold <- subset(max_YI, max_YI$TH == min(max_YI$TH))$TH

threshold ## c_through threshold in ng/mL to discriminate patients with risk for dose reduction from the rest

ggplot(data=roc_data) + geom_point(aes(x=FP, y=TP)) + geom_abline(aes(intercept=0,slope=1))



### Create the ROC curve with package pROC

df <- cbind(df, class= as.numeric(df$c_through >thresh))



pred.z.01 <- prediction(df$c_through, df$group)
perf.z.01 <- performance(pred.z.01, "tpr", "fpr")

plot.new()
plot(perf.z.01, col="green")
abline(0,1,col="grey")
auc.z.01 <- performance(pred.z.01, "auc")
legend("bottomright", 
       paste(round(as.numeric(auc.z.01@y.values),
                   digits= 2)),
       col=c("green"),
       pch = c(3))

