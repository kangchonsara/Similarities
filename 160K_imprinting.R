########################################################################
# Analysis of H3 imprinting 
# Zost et al. (2017)
#
########################################################################

########################################################################
## Parameters
########################################################################

H3_AR <- 0.25 # mean H3 attack rate
maxY <- 12 # oldest age by which first infection with H3 must occur

########################################################################
# Fraction of H3 (of all flu A) per year, normalized to obtain intensity
# (and then detrended)
########################################################################

H3_frac <- read.table("H3_subtype_frac.csv",header=TRUE,sep=",")
H3_int <- data.frame(year = H3_frac$year, int = H3_frac$H3/mean(H3_frac$H3))

int_fit <- lm(H3_int$int ~ H3_int$year)
# > summary(int_fit)

# Call:
# lm(formula = H3_int$int ~ H3_int$year)

# Residuals:
    # Min      1Q  Median      3Q     Max 
# -0.6285 -0.1348  0.0662  0.1847  0.3875 

# Coefficients:
             # Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 32.055186   6.162570   5.202 5.85e-06 ***
# H3_int$year -0.015613   0.003098  -5.039 9.88e-06 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 0.2521 on 41 degrees of freedom
# Multiple R-squared:  0.3825,	Adjusted R-squared:  0.3674 
# F-statistic:  25.4 on 1 and 41 DF,  p-value: 9.879e-06

for (y in 1:length(H3_int$year)) {
	H3_int$adj[y] <- H3_int$int[y]+0.015613*(y-1)
}

########################################################################
# Frequency of 160 alleles by year for H3 isolates
########################################################################

f160 <- read.table("160_freq.csv",header=TRUE,sep=",")
f160_aa <- cbind(f160$T,f160$K,f160$A,f160$R,f160$Others)
colnames(f160_aa) <- c('T','K','A','R','other')
T_col <- rgb(0.5,0,0.5,0.75)
K_col <- rgb(1,0,0,0.75)
A_col <- rgb(0,0.5,0,0.75)
R_col <- rgb(0.5,0.5,0.5,0.75)
oth_col <- rgb(0,0,0)
col_map <- c(T_col,K_col,A_col,R_col,oth_col)


pdf("Frequencies_by_year_barplot.pdf",width=6,height=3)
par(mar=c(4.1, 4.1, 0.9, 4.1), xpd=TRUE)
pl <- barplot(as.matrix(t(f160_aa)),axes=F,ylab="Allele frequency",ylim=c(0,1),col=col_map)
axis(1,labels=f160$year,at=pl)
axis(2)
legend(x="right",legend=colnames(f160_aa),bty="n",inset=c(-0.15,0),pch=15,col=col_map)
dev.off()

########################################################################
# Reconstruction by birth year
########################################################################

imp160 <- data.frame(yob=seq(from=1968,to=1998))
imp160$T <- 0
imp160$K <- 0
imp160$A <- 0
imp160$R <- 0
imp160$other <- 0
for (by in 1:length(imp160$yob)) {
	yob <- imp160$yob[by]
	prob_uninf <- 1
	for (y in 1:maxY) {
		thisY <- yob + y # assuming first year of risk is year after birth
		thisY_int <- subset(H3_int,H3_int$year==thisY)$adj
		AR_this_year <- H3_AR*thisY_int
		prob_first_inf <- prob_unif*AR_this_year
		imp160$T[by] <- imp160$T[by] + prob_first_inf*subset(f160,f160$year==thisY)$T
		imp160$K[by] <- imp160$K[by] + prob_first_inf*subset(f160,f160$year==thisY)$K
		imp160$A[by] <- imp160$A[by] + prob_first_inf*subset(f160,f160$year==thisY)$A
		imp160$R[by] <- imp160$R[by] + prob_first_inf*subset(f160,f160$year==thisY)$R
		imp160$other[by] <- imp160$other[by] + prob_first_inf*subset(f160,f160$year==thisY)$Others		
		prob_uninf <- prob_uninf*(1-AR_this_year) 
	}	
}

imp160_aa <- cbind(imp160$T,imp160$K, imp160$A, imp160$R, imp160$other)
imp160_aa <- imp160_aa/rowSums(imp160_aa) # kluge for asymptotic infection risk
colnames(imp160_aa) <- c('T','K','A','R','other')

pdf("160_imprinting_barplot.pdf",width=6,height=3)
par(mar=c(4.1, 4.1, 0.9, 4.1), xpd=TRUE)
pl <- barplot(as.matrix(t(imp160_aa)),axes=F,ylab="Probability",xlab="Birth year",ylim=c(0,1),col=col_map)
axis(1,labels=imp160$yob,at=pl)
axis(2)
legend(x="right",legend=colnames(imp160_aa),bty="n",inset=c(-0.15,0),pch=15,col=col_map)
dev.off()

