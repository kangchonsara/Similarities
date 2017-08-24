#############################################
# Analysis of VE Zost et al. (2017) 
#
#
# This script includes the code to generate 
# the figures and perform statistical
# analysis of the serological data.
#
#
# SECTIONS
# - Data import & preparation (must execute)
# - (Figures) Pre- & post-vaccination trends
# - (Statistics) Correlations by age
# - (Figures) Fold change by vaccine group
# - (Figures) Vaccination history & boost
#
#############################################


# WT = glycosylated (T160)
# 160K = unglycosylated

# FZ = Fluzone, egg 
# FB = Flublok, recombinant 
# FCV = Flucelvax, MDCK 

require(dplyr)
require(ggplot2)

############################################
####### Data import & preparation ##########
############################################

# Import FRNT data (with two columns of Vx history)
titers_raw <- read.table("2017-04-19_FRNT_for_import.csv",header=TRUE,sep=",",na.strings="Unknown")
H3_exp <- read.table("H3_exp.csv",sep=",")
colnames(H3_exp) <- c("YOB","pH3")
T160K_exp <- read.table("160K_exp.csv",sep=",",header=FALSE)
colnames(T160K_exp) <- c("YOB","p160K")

# Merge data and add basic statistics
titers <- titers_raw
titers$YOB <- 2017 - titers$Age

geo_mean <- function(t1,t2,t3) {
	tmp_df <- data.frame(t1,t2,t3)
	tmp_df <- log(tmp_df)
	tmp_df <- rowMeans(tmp_df)
	tmp_df <- exp(tmp_df)
	return(tmp_df)
}
titers$GMT_Pre_WT <- geo_mean(titers$Pre_WT_T1,titers$Pre_WT_T2,titers$Pre_WT_T3)
titers$GMT_Pre_160K <- geo_mean(titers$Pre_160K_T1,titers$Pre_160K_T2,titers$Pre_160K_T3)
titers$GMT_Post_WT <- geo_mean(titers$Post_WT_T1,titers$Post_WT_T2,titers$Post_WT_T3)
titers$GMT_Post_160K <- geo_mean(titers$Post_160K_T1,titers$Post_160K_T2,titers$Post_160K_T3)
titers$GM_FC_WT <- titers$GMT_Post_WT/titers$GMT_Pre_WT
titers$GM_FC_160K <- titers$GMT_Post_160K/titers$GMT_Pre_160K
titers$Pre_ratio_WT_160K <- titers$GMT_Pre_WT/titers$GMT_Pre_160K


# Define ordinal variables based on vaccination history
# and add probability of primary exposure to H3
titers$Vx_hx <- NA
titers$H3_exp <- NA
titers$T160K_exp <- NA
for (p in 1:dim(titers)[1]) {
	if ( is.na(titers$Vx_2014[p])==FALSE & is.na(titers$Vx_2015[p])==FALSE ) {
		if ( (titers$Vx_2014[p] == 'Inactivated') & (titers$Vx_2015[p] =='Inactivated') ) {
			titers$Vx_hx[p] <- 2
		} else if ( (titers$Vx_2014[p] == 'Inactivated') & (titers$Vx_2015[p] =='None') ) {
				titers$Vx_hx[p] <- 1
		} else if ( (titers$Vx_2014[p] == 'None') & (titers$Vx_2015[p]=='Inactivated') ) {
			titers$Vx_hx[p] <- 1
		} else if ( (titers$Vx_2014[p] == 'None') & (titers$Vx_2015[p] =='None') ) {
			titers$Vx_hx[p] <- 0
		}
	}
	titers$H3_exp[p] <- subset(H3_exp,H3_exp$YOB==titers$YOB[p])$pH3
	titers$T160K_exp[p] <- subset(T160K_exp,T160K_exp$YOB==titers$YOB[p])$p160K
}

#### Calculate people with any unknown status
any_unknown <- titers %>% filter(is.na(Vx_2014) | is.na(Vx_2015))
nrow(any_unknown)

#### Define partitions
titers.FB <- subset(titers,titers$Group=="FB")
titers.FCV <- subset(titers,titers$Group=="FCV")
titers.FZ <- subset(titers,titers$Group=="FZ")
titers.FB_FZ <- titers %>% filter(Group=="FB" | Group=="FZ")
titers.FB_FCV <- titers %>% filter(Group=="FB" | Group=="FCV")
No_vax <- titers %>% filter(Vx_2014=="None" & Vx_2015=="None")
Only_2014 <- titers %>% filter(Vx_2014=="Inactivated" & Vx_2015=="None")
Only_2015 <- titers %>% filter(Vx_2014=="None" & Vx_2015=="Inactivated")
Both <- titers %>% filter(Vx_2014=="Inactivated" & Vx_2015=="Inactivated")
titers.post1979 <- titers %>% filter(YOB >= 1979)
titers.pre1979 <- titers %>% filter(YOB < 1979)

both2015 <- data.frame(rbind(Both,Only_2015))
none2014 <- data.frame(rbind(Only_2014,No_vax))
nonboth_df <- data.frame(rbind(No_vax,Only_2014,Only_2015))
both_df <- data.frame(Both)
any_df <- data.frame(rbind(Both,Only_2014,Only_2015))
none_df <- data.frame(No_vax)

FZ_col <- rgb(237, 163, 28,100,maxColorValue=255)
FCV_col <- rgb(41,137,192,100,maxColorValue=255)
FB_col <- rgb(18,153,20,100,maxColorValue=255)
WT_col <- rgb(0,0,1,1/2)
K_col <- rgb(1,0,0,1/2)

############################################
# (Figures) Pre- and post-vaccination trends 
############################################

pdf("Pre_post_vacc_by_YOB.pdf",width=6,height=6)
par(mfrow=c(2,2), mar=c(4.1, 4, 0.9, 0.8) + 0.1)
plot(jitter(titers$YOB),jitter(titers$GMT_Pre_WT,2),log="y",xlab="Birth year",ylab="Pre-vaccination WT GMT",xlim=c(1967,2000),ylim=c(7,10245),pch=16,col=WT_col,axes=FALSE)
axis(1)
axis(2, at = c(10,80,640,5120))

plot(jitter(titers$YOB),jitter(titers$GMT_Pre_160K,2),log="y",xlab="Birth year",ylab="Pre-vaccination 160K GMT",xlim=c(1967,2000),ylim=c(7,10245),pch=16,col=K_col,axes=FALSE)
axis(1)
axis(2, at = c(10,80,640,5120))

plot(jitter(titers$YOB),jitter(titers$GMT_Post_WT,2),log="y",xlab="Birth year",ylab="Post-vaccination WT GMT",xlim=c(1967,2000),ylim=c(7,10245),pch=16,col=WT_col,axes=FALSE)
axis(1)
axis(2, at = c(10,80,640,5120))

plot(jitter(titers$YOB),jitter(titers$GMT_Post_160K,2),log="y",xlab="Birth year",ylab="Post-vaccination 160K GMT",xlim=c(1967,2000),ylim=c(7,10245),pch=16,col=K_col,axes=FALSE)
axis(1)
axis(2, at = c(10,80,640,5120))
dev.off()


############################################
#### (Statistics) Correlations by age ######
############################################

cor.test(titers$Age,titers$GMT_Pre_WT,method='spearman')
# Spearman's rank correlation rho

# data:  titers$Age and titers$GMT_Pre_WT
# S = 65015, p-value = 0.2563
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
       # rho 
# -0.1375183 

# Warning message:
# In cor.test.default(titers$Age, titers$GMT_Pre_WT, method = "spearman") :
  # Cannot compute exact p-value with ties


cor.test(titers$Age,titers$GMT_Pre_160K,method='spearman')
	# Spearman's rank correlation rho

# data:  titers$Age and titers$GMT_Pre_160K
# S = 75621, p-value = 0.006371
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
     # rho 
# -0.32309 

# Warning message:
# In cor.test.default(titers$Age, titers$GMT_Pre_160K, method = "spearman") :
  # Cannot compute exact p-value with ties 

cor.test(titers$Age,titers$GMT_Post_WT,method='spearman')
	# Spearman's rank correlation rho

# data:  titers$Age and titers$GMT_Post_WT
# S = 68533, p-value = 0.09849
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
      # rho 
# -0.199081 

# Warning message:
# In cor.test.default(titers$Age, titers$GMT_Post_WT, method = "spearman") :
  # Cannot compute exact p-value with ties

cor.test(titers$Age,titers$GMT_Post_160K,method='spearman')
# Spearman's rank correlation rho

# data:  titers$Age and titers$GMT_Post_160K
# S = 79850, p-value = 0.0006653
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
       # rho 
# -0.3970838 

# Warning message:
# In cor.test.default(titers$Age, titers$GMT_Post_160K, method = "spearman") :
  # Cannot compute exact p-value with ties
  
  
############################################
## (Figures) Fold change by vaccine group ##
############################################

pdf("FC_group.pdf",width=6,height=3.2)
#postscript("FC_group_easy_axes.ps")
par(mfrow=c(1,2),mar=c(4, 4, 1, 1) + 0.1)
plot(titers$Group,log(titers$GM_FC_WT,2),xlab="Vaccine type",ylab=expression("WT fold change"),col=alpha('blue',0.5),ylim=c(-1,7),axes=FALSE)
axis(1,at=c(1,2,3),labels=c("FB","FCV","FZ"))
axis(2,at=c(0,2,4,6),labels=c(1,4,16,64))
plot(titers$Group,log(titers$GM_FC_160K,2),xlab="Vaccine type",ylab=expression("160K fold change"),col=alpha('red',0.5),ylim=c(-1,7),axes=FALSE)
axis(1,at=c(1,2,3),labels=c("FB","FCV","FZ"))
axis(2,at=c(0,2,4,6),labels=c(1,4,16,64))
dev.off()


############################################
## (Figures) Vaccination history and boost #
############################################
pdf("Pre_FC_Vx_hx.pdf",width=6,height=6)
par (mfrow=c(2,2))

# 160K pre-vacc titer by vaccination history
plot(jitter(titers$Vx_hx),titers$GMT_Pre_160K,xlab="Vaccination history",ylab="Pre-vaccination 160K titer",log="y",ylim=c(7,10245),axes=FALSE,col=alpha('red',0.5),pch=16)
axis(1,at=c(0,1,2),labels=c("None","One","Both"))
axis(2, at = c(10,80,640,5120))

# 160K change by vaccination history
plot(jitter(titers$Vx_hx),titers$GM_FC_160K,xlab="Vaccination history",ylab="160K fold change",ylim=c(0.4,81),axes=FALSE,log="y",col=alpha('red',0.5),pch=16)
axis(1,at=c(0,1,2),labels=c("None","One","Both"))
axis(2)

# WT pre-vacc by vaccination history
plot(jitter(titers$Vx_hx),titers$GMT_Pre_WT,xlab="Vaccination history",ylab="Pre-vaccination WT titer",log="y",ylim=c(7,10245),axes=FALSE,col=alpha('blue',0.5),pch=16)
axis(1,at=c(0,1,2),labels=c("None","One","Both"))
axis(2, at = c(10,80,640,5120))

# WT change by vaccination history
plot(jitter(titers$Vx_hx),titers$GM_FC_WT,xlab="Vaccination history",ylab="WT fold change",log="y",ylim=c(0.4,81),axes=FALSE,col=alpha('blue',0.5),pch=16)
axis(1,at=c(0,1,2),labels=c("None","One","Both"))
axis(2)
dev.off()  

############################################
########## Regression models ###############
############################################

############################################
## Prevacc titers to K160

# Age effects
fit <- lm( log(GMT_Pre_160K,2) ~ Age + factor(Vx_hx), data=titers)
summary(fit); logLik(fit); extractAIC(fit)

# Call:
# lm(formula = log(GMT_Pre_160K, 2) ~ Age + factor(Vx_hx), data = titers)

# Residuals:
    # Min      1Q  Median      3Q     Max 
# -3.9844 -1.0637  0.2044  1.1688  3.2474 

# Coefficients:
               # Estimate Std. Error t value Pr(>|t|)    
# (Intercept)     8.10997    0.88653   9.148 7.57e-13 ***
# Age            -0.07295    0.02448  -2.980   0.0042 ** 
# factor(Vx_hx)1  0.69526    0.71296   0.975   0.3335    
# factor(Vx_hx)2  0.65545    0.64459   1.017   0.3134    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 1.726 on 58 degrees of freedom
  # (8 observations deleted due to missingness)
# Multiple R-squared:  0.1424,	Adjusted R-squared:  0.09801 
# F-statistic: 3.209 on 3 and 58 DF,  p-value: 0.02952

# 'log Lik.' -119.739 (df=5)
# [1]  4.00000 71.52954


# K160 effects
fit <- lm( log(GMT_Pre_160K,2) ~ T160K_exp + factor(Vx_hx), data=titers)
summary(fit); logLik(fit); extractAIC(fit)

# Call:
# lm(formula = log(GMT_Pre_160K, 2) ~ T160K_exp + factor(Vx_hx), 
    # data = titers)

# Residuals:
    # Min      1Q  Median      3Q     Max 
# -3.5892 -0.8808  0.0517  1.0969  3.3825 

# Coefficients:
               # Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      4.3926     0.8458   5.193 2.78e-06 ***
# T160K_exp        1.9661     0.7085   2.775  0.00742 ** 
# factor(Vx_hx)1   0.4940     0.7181   0.688  0.49430    
# factor(Vx_hx)2   0.5838     0.6492   0.899  0.37227    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 1.741 on 58 degrees of freedom
  # (8 observations deleted due to missingness)
# Multiple R-squared:  0.127,	Adjusted R-squared:  0.08179 
# F-statistic: 2.811 on 3 and 58 DF,  p-value: 0.04725

# 'log Lik.' -120.2914 (df=5)
# [1]  4.00000 72.63443

############################################
## Prevacc titers to WT

# Age effects
fit <- lm( log(GMT_Pre_WT,2) ~ Age + factor(Vx_hx), data=titers)
summary(fit); logLik(fit); extractAIC(fit)
# Call:
# lm(formula = log(GMT_Pre_WT, 2) ~ Age + factor(Vx_hx), data = titers)

# Residuals:
     # Min       1Q   Median       3Q      Max 
# -1.92020 -1.08809 -0.07095  0.83215  2.53211 

# Coefficients:
               # Estimate Std. Error t value Pr(>|t|)    
# (Intercept)     4.96254    0.61633   8.052 5.02e-11 ***
# Age            -0.02056    0.01702  -1.208    0.232    
# factor(Vx_hx)1  0.69079    0.49567   1.394    0.169    
# factor(Vx_hx)2  0.17487    0.44813   0.390    0.698    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 1.2 on 58 degrees of freedom
  # (8 observations deleted due to missingness)
# Multiple R-squared:  0.06539,	Adjusted R-squared:  0.01705 
# F-statistic: 1.353 on 3 and 58 DF,  p-value: 0.2663

# 'log Lik.' -97.20013 (df=5)
# [1]  4.00000 26.45189


# K160 effects
fit <- lm( log(GMT_Pre_WT,2) ~ T160K_exp + factor(Vx_hx), data=titers)
summary(fit); logLik(fit); extractAIC(fit)

# Call:
# lm(formula = log(GMT_Pre_WT, 2) ~ T160K_exp + factor(Vx_hx), 
    # data = titers)

# Residuals:
     # Min       1Q   Median       3Q      Max 
# -1.80732 -1.16085  0.00002  0.82933  2.71024 

# Coefficients:
               # Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      3.6870     0.5760   6.401 2.95e-08 ***
# T160K_exp        0.8165     0.4826   1.692    0.096 .  
# factor(Vx_hx)1   0.6258     0.4891   1.280    0.206    
# factor(Vx_hx)2   0.1621     0.4422   0.367    0.715    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 1.186 on 58 degrees of freedom
  # (8 observations deleted due to missingness)
# Multiple R-squared:  0.08695,	Adjusted R-squared:  0.03972 
# F-statistic: 1.841 on 3 and 58 DF,  p-value: 0.1498

# 'log Lik.' -96.4768 (df=5)
# [1]  4.00000 25.00523

############################################
## Postvacc titers to K160

# Age effects
fit <- lm( log(GMT_Post_160K,2) ~ Age + factor(Vx_hx) + Group + log(GMT_Pre_160K,2), data=titers)
summary(fit); logLik(fit); extractAIC(fit)

# Call:
# lm(formula = log(GMT_Post_160K, 2) ~ Age + factor(Vx_hx) + Group + 
    # log(GMT_Pre_160K, 2), data = titers)

# Residuals:
     # Min       1Q   Median       3Q      Max 
# -2.09618 -0.78038 -0.08893  0.77508  2.39687 

# Coefficients:
                     # Estimate Std. Error t value Pr(>|t|)    
# (Intercept)           7.74734    1.00506   7.708 2.61e-10 ***
# Age                  -0.03798    0.01804  -2.105 0.039838 *  
# factor(Vx_hx)1       -1.38345    0.48948  -2.826 0.006551 ** 
# factor(Vx_hx)2       -1.85282    0.45119  -4.107 0.000134 ***
# GroupFCV             -0.36871    0.36462  -1.011 0.316334    
# GroupFZ              -0.08527    0.38386  -0.222 0.825020    
# log(GMT_Pre_160K, 2)  0.46420    0.09288   4.998 6.24e-06 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 1.172 on 55 degrees of freedom
  # (8 observations deleted due to missingness)
# Multiple R-squared:  0.5352,	Adjusted R-squared:  0.4844 
# F-statistic: 10.55 on 6 and 55 DF,  p-value: 9.078e-08

# 'log Lik.' -94.12552 (df=8)
# [1]  7.00000 26.30267


# K160 effects
fit <- lm( log(GMT_Post_160K,2) ~ T160K_exp + factor(Vx_hx) + Group + log(GMT_Pre_160K,2), data=titers)
summary(fit); logLik(fit); extractAIC(fit)

# Call:
# lm(formula = log(GMT_Post_160K, 2) ~ T160K_exp + factor(Vx_hx) + 
    # Group + log(GMT_Pre_160K, 2), data = titers)

# Residuals:
     # Min       1Q   Median       3Q      Max 
# -1.98130 -0.63666 -0.05465  0.59496  2.46341 

# Coefficients:
                     # Estimate Std. Error t value Pr(>|t|)    
# (Intercept)           5.77175    0.72138   8.001 8.69e-11 ***
# T160K_exp             1.18691    0.53950   2.200  0.03203 *  
# factor(Vx_hx)1       -1.49163    0.48488  -3.076  0.00326 ** 
# factor(Vx_hx)2       -1.92834    0.44808  -4.304 6.96e-05 ***
# GroupFCV             -0.45796    0.36877  -1.242  0.21956    
# GroupFZ               0.02737    0.38904   0.070  0.94416    
# log(GMT_Pre_160K, 2)  0.45201    0.09423   4.797 1.27e-05 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 1.168 on 55 degrees of freedom
  # (8 observations deleted due to missingness)
# Multiple R-squared:  0.5383,	Adjusted R-squared:  0.488 
# F-statistic: 10.69 on 6 and 55 DF,  p-value: 7.607e-08

# 'log Lik.' -93.91377 (df=8)
# [1]  7.00000 25.87916

############################################
## Postvacc titers to WT

# Age effects
fit <- lm( log(GMT_Post_WT,2) ~ Age + factor(Vx_hx) + Group + log(GMT_Pre_WT,2), data=titers)
summary(fit); logLik(fit); extractAIC(fit)

# Call:
# lm(formula = log(GMT_Post_WT, 2) ~ Age + factor(Vx_hx) + Group + 
    # log(GMT_Pre_WT, 2), data = titers)

# Residuals:
    # Min      1Q  Median      3Q     Max 
# -3.2809 -0.7337 -0.2436  0.8351  2.7716 

# Coefficients:
                   # Estimate Std. Error t value Pr(>|t|)    
# (Intercept)         4.94816    1.02500   4.827 1.14e-05 ***
# Age                -0.01837    0.01897  -0.968  0.33719    
# factor(Vx_hx)1     -1.66731    0.55470  -3.006  0.00399 ** 
# factor(Vx_hx)2     -2.20743    0.50569  -4.365 5.65e-05 ***
# GroupFCV           -1.04215    0.40471  -2.575  0.01274 *  
# GroupFZ            -0.89441    0.43076  -2.076  0.04255 *  
# log(GMT_Pre_WT, 2)  0.84208    0.14638   5.753 4.02e-07 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 1.318 on 55 degrees of freedom
  # (8 observations deleted due to missingness)
# Multiple R-squared:  0.5547,	Adjusted R-squared:  0.5062 
# F-statistic: 11.42 on 6 and 55 DF,  p-value: 2.974e-08

# 'log Lik.' -101.3966 (df=8)
# [1]  7.00000 40.84475


# K160 effects
fit <- lm( log(GMT_Post_WT,2) ~ T160K_exp + factor(Vx_hx) + Group + log(GMT_Pre_WT,2), data=titers)
summary(fit); logLik(fit); extractAIC(fit)

# Call:
# lm(formula = log(GMT_Post_WT, 2) ~ T160K_exp + factor(Vx_hx) + 
    # Group + log(GMT_Pre_WT, 2), data = titers)

# Residuals:
    # Min      1Q  Median      3Q     Max 
# -3.2706 -0.8431 -0.2298  0.7079  2.9355 

# Coefficients:
                   # Estimate Std. Error t value Pr(>|t|)    
# (Intercept)          4.0998     0.8654   4.737 1.57e-05 ***
# T160K_exp            0.4408     0.5739   0.768  0.44581    
# factor(Vx_hx)1      -1.7134     0.5541  -3.092  0.00312 ** 
# factor(Vx_hx)2      -2.2431     0.5065  -4.429 4.55e-05 ***
# GroupFCV            -1.0723     0.4094  -2.619  0.01137 *  
# GroupFZ             -0.8548     0.4378  -1.953  0.05596 .  
# log(GMT_Pre_WT, 2)   0.8349     0.1502   5.559 8.21e-07 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 1.322 on 55 degrees of freedom
  # (8 observations deleted due to missingness)
# Multiple R-squared:  0.5519,	Adjusted R-squared:  0.5031 
# F-statistic: 11.29 on 6 and 55 DF,  p-value: 3.497e-08

# 'log Lik.' -101.5898 (df=8)
# [1]  7.00000 41.23129


