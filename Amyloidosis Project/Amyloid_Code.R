library(dplyr)
setwd("C:/Users/annik/OneDrive/Desktop/BU Fall 2024/MA675/Amyloid Project")

###
am_neg <- read.csv("Amyloid_negative.csv")
am_pos <- read.csv("Amyloid_positive.csv")
## assigning column names
am_pos$X.1 <- NULL
am_pos$X.2 <- NULL
am_pos$X.3 <- NULL
a <- which(colnames(am_neg) != colnames(am_pos))
colnames(am_neg)[a]
colnames(am_pos)[a]
###


am_neg$amyloid_status <- 0
am_pos$amyloid_status <- 1
amyloid <- rbind(am_neg, am_pos)
amyloid <- amyloid[1:189 ,]
amyloid$amyloid_status[c(150, 151, 152)] <- 0

amyloid$amyloid_status <- as.factor(amyloid$amyloid_status)
plot(amyloid$amyloid_status, amyloid$Points, xlab= "Amyloidosis Status", 
     ylab="Nomogram Score")
qnt <- quantile(amyloid$Age , seq(0,1,.25))

amyloid$agegroup <- cut(amyloid$Age,unique(qnt),include.lowest=TRUE)
plot(amyloid$agegroup, amyloid$amyloid_status , xlab="Age Categories", ylab="Amyloid Status")

### Make Grade ordered numeric
### Fit GLM of amyloid as a model of points, simple linear regression
## report true pos and false pos rate, specificity and sensitivity
## Amyloid as a model of risk factors used for nomogram report the same
## Then code severity ("Grade") as ordered numeric and add to model
## report new findings

## Ask michael about 0 in grade column
## read pages 1-5 up until recoverability

##################################################################
## Null model
amy_null <- glm(amyloid_status ~ 1 , data=amyloid, family = binomial)
amy_null

### Restructure for factor variables used as risk scores
amyloid$Race <- factor(amyloid$Race, levels = c("White", "AA"), labels = c("white", "black"))
amyloid$Monoclonal.Gammopathy <- as.factor(amyloid$Monoclonal.Gammopathy)
amyloid$Rheumatoid.Arthritis <- as.factor(amyloid$Rheumatoid.Arthritis)
amyloid$Afib <- as.factor(amyloid$Afib)
amyloid$Degenerative.Spine.Disease <- as.factor(amyloid$Degenerative.Spine.Disease)
amyloid$Bilateral. <- as.factor(amyloid$Bilateral.)
###################################################################
### First Model Amyloid ~ Nomogram
library(pROC)
library(ggplot2)
library(lmtest)

amy_nomo <- glm(amyloid_status ~ Points, data = amyloid, family=binomial)
summary(amy_nomo)
amyloid$predicted_prob <- predict(amy_nomo, type = "response")


# ROC curve
roc_curve <- roc(amyloid$amyloid_status, amyloid$predicted_prob)

# Extract AUC
auc_value <- auc(roc_curve)

# Plot the ROC curve
ggroc <- ggplot(data.frame(tpr = roc_curve$sensitivities, 
                           fpr = 1 - roc_curve$specificities), aes(x = fpr, y = tpr)) +
  geom_line(color = "blue", size = 1.2) +
  geom_abline(linetype = "dashed") +
  labs(title = paste0("ROC Curve For Nomogram Simple Model (AUC = ", round(auc_value, 2), ")"),
       x = "False Positive Rate (1 - Specificity)", 
       y = "True Positive Rate (Sensitivity)") +
  theme_minimal()

print(ggroc)

tpr <- roc_curve$sensitivities
fpr <- 1 - roc_curve$specificities    
fnr <- 1 - roc_curve$sensitivities    
tnr <- roc_curve$specificities        

roc_metrics <- data.frame(Threshold = roc_curve$thresholds,
                          TPR = tpr, 
                          FPR = fpr, 
                          FNR = fnr, 
                          TNR = tnr)
roc_metrics

AIC(amy_nomo)
#################################################################
#### Second Model Amyloid ~ Risk Factors

risk_factor_model <- glm(amyloid_status~ Afib + Degenerative.Spine.Disease+ Bilateral. +
                           Monoclonal.Gammopathy + Rheumatoid.Arthritis + Race + Age ,
                         data = amyloid, family = binomial)
summary(risk_factor_model)
#plot(risk_factor_model)
amyloid$predicted_risk <- predict(risk_factor_model, type="response")

lrtest(amy_nomo, amy_null)

roc_curve1 <- roc(amyloid$amyloid_status, amyloid$predicted_risk)

# Extract AUC
auc_value1 <- auc(roc_curve1)

# Plot the ROC curve
ggroc1 <- ggplot(data.frame(tpr = roc_curve1$sensitivities, 
                           fpr = 1 - roc_curve1$specificities), aes(x = fpr, y = tpr)) +
  geom_line(color = "blue", size = 1.2) +
  geom_abline(linetype = "dashed") +
  labs(title = paste0("ROC Curve For Risk Factor Model (AUC = ", round(auc_value1, 2), ")"),
       x = "False Positive Rate (1 - Specificity)", 
       y = "True Positive Rate (Sensitivity)") +
  theme_minimal()

print(ggroc1)

tpr1 <- roc_curve1$sensitivities        
fpr1 <- 1 - roc_curve1$specificities    
fnr1 <- 1 - roc_curve1$sensitivities    
tnr1 <- roc_curve1$specificities       


roc_metrics1 <- data.frame(Threshold = roc_curve1$thresholds,
                          TPR = tpr1, 
                          FPR = fpr1, 
                          FNR = fnr1, 
                          TNR = tnr1)
roc_metrics1


#########################################################################
### Model 3 Using Severity as well
amyloid$Grade <- factor(amyloid$Grade, 
                   levels = c(0, "mild", "moderate", "severe"), 
                   labels = c("none", "mild", "moderate", "severe"))

grade_model <- glm(amyloid_status~ Afib + Degenerative.Spine.Disease+ Bilateral. +
                           Monoclonal.Gammopathy + Rheumatoid.Arthritis + Race + Age
                   + Grade, data = amyloid, family = binomial)
summary(grade_model)

amyloid$predicted_grade <- predict(grade_model, type="response")

roc_curve2 <- roc(amyloid$amyloid_status, amyloid$predicted_grade)

# Extract AUC
auc_value2 <- auc(roc_curve2)

# Plot the ROC curve
ggroc2 <- ggplot(data.frame(tpr = roc_curve2$sensitivities, 
                            fpr = 1 - roc_curve2$specificities), aes(x = fpr, y = tpr)) +
  geom_line(color = "blue", size = 1.2) +
  geom_abline(linetype = "dashed") +
  labs(title = paste0("ROC Curve For Risk Factor And Severity (AUC = ", round(auc_value2, 2), ")"),
       x = "False Positive Rate (1 - Specificity)", 
       y = "True Positive Rate (Sensitivity)") +
  theme_minimal()

print(ggroc2)

tpr2 <- roc_curve2$sensitivities        
fpr2 <- 1 - roc_curve2$specificities    
fnr2 <- 1 - roc_curve2$sensitivities    
tnr2 <- roc_curve2$specificities       


roc_metrics2 <- data.frame(Threshold = roc_curve2$thresholds,
                           TPR = tpr2, 
                           FPR = fpr2, 
                           FNR = fnr2, 
                           TNR = tnr2)
roc_metrics2

###########################################################################

## Extra things to add for presentation tomorrow
###AIC for all three models

AIC_df <- data.frame(
  Model = c("Simple Amyloidosis", "Risk Factor Model", "Risk Factors and Severity"),
  AIC = c(AIC(amy_nomo), AIC(risk_factor_model), AIC(grade_model))
)
AIC_df
#View(AIC_df)
########################################################################
######### Binned plots
arm::binnedplot(fitted(amy_nomo), resid(amy_nomo), xlab = "Fitted Values", ylab = "Residuals", main = "Binned Plot Simple Amyloidosis Model")
arm::binnedplot(fitted(risk_factor_model), resid(risk_factor_model), xlab = "Fitted Values", ylab = "Residuals",
                main="Binned Plot Risk Factor Model")
arm::binnedplot(fitted(grade_model), resid(grade_model), xlab = "Fitted Values", ylab = "Residuals",
                main = "Binned Plot Risk Factor And Severity Model")

#########################################################################
######## Table for severity and amyloidosis status

table(amyloid$amyloid_status, amyloid$Grade)

######################################################################
#######################################################################

## AGE COMPARISONS

all_data <- read.csv("Amyloid All Patients.csv")
View(all_data)

sampled_ages <- all_data[1:189,1]
unsampled_ages <- all_data[190:341, 1]
summary(sampled_ages)
summary(unsampled_ages)
par(mfrow=c(1,2))
boxplot(sampled$Age, main="Boxplot of Sampled Age",ylim=c(20, 100), col="darkblue")
boxplot(unsampled$Age, main="Boxplot of Unsampled Age", ylim=c(20, 100), col="darkorange")


##################################################################
################################################################
## points comparisons

which(colnames(all_data)  == "Points")

sampled <- all_data[1:189, c(1, 31)]
unsampled <- all_data[190:341, c(1, 31)]

summary(sampled$Points)
summary(unsampled$Points)

### 
par(mfrow=c(1,2))
boxplot(sampled$Points, main="Boxplot of Sampled Points",ylim=c(25, 185), col="darkblue")
boxplot(unsampled$Points, main="Boxplot of Unsampled Points", ylim=c(25, 185), col="darkorange")


boxplot(sampled$Points, main="Comparison of Sampled & Unsampled Points", ylim=c(25, 185), col="lightblue")
boxplot(unsampled$Points, add=TRUE, col=rgb(1, 0, 0, 0.5))  # Semi-transparent red

#######################################################################
## Model
all_data <- all_data[1:341, ]
all_data <- all_data %>%
  mutate(indicator = if_else(row_number() %in% c(1:189), 1, 0))

indicatormod <- glm(indicator ~ Points, data=all_data, family=binomial) 

summary(indicatormod)
arm::binnedplot(fitted(indicatormod), resid(indicatormod), xlab = "Fitted Values", ylab = "Residuals", main = "Binned Plot Simple Indicator Model")

all_data$predicted_indicator <- predict(indicatormod, type="response")

roc_curve_ind <- roc(all_data$indicator, all_data$predicted_indicator)

# Extract AUC
auc_value_ind <- auc(roc_curve_ind)

# Plot the ROC curve
ggroc_ind <- ggplot(data.frame(tpr = roc_curve_ind$sensitivities, 
                            fpr = 1 - roc_curve_ind$specificities), aes(x = fpr, y = tpr)) +
  geom_line(color = "blue", size = 1.2) +
  geom_abline(linetype = "dashed") +
  labs(title = paste0("ROC Curve For Risk Factor And Sample Status (AUC = ", round(auc_value_ind, 2), ")"),
       x = "False Positive Rate (1 - Specificity)", 
       y = "True Positive Rate (Sensitivity)") +
  theme_minimal()

print(ggroc_ind)



#### Model with weights
all_data$weight <- ifelse(all_data$indicator == 1, 1 / all_data$predicted_indicator, NA)
amyloid <- amyloid %>%
  mutate(weight = all_data$weight[1:189])


amyloid$severity <- ifelse(amyloid$Grade == "severe", 1 ,0)
fit_weighted <- glm(amyloid_status ~ Points + severity,
                    
                    data = amyloid,
                    family = binomial,
                    weights = weight)
summary(fit_weighted)
