################################################################################
# File name: cjh12X1.R                                                         #
# Last modified by: Ralph Hurtado                                              #
# Original file by: Christopher J. House                                       #
#                                                                              #
# Purpose: This file calculates diagnostics for a linear, weighted least       #
#          squares model.                                                      #
# Inputs: none                                                                 #
# Outputs: graphical and non-graphical tests of goodness of fit.               #
################################################################################

## LINEAR with Weighting
dfw=NULL
y=NULL
x=c(5,10,25,50,100,250,500,1000)
dfw=data.frame(x)
dfw
##      x
## 1    5
## 2   10
## 3   25
## 4   50
## 5  100
## 6  250
## 7  500
## 8 1000
k=1

beta <- 9.81E-03
sgsq0 <- 6.06E-02 / 1000  # constant variance term
Vmat <- rnorm(8,0,x) # diagonal terms of V matrix, with 1/x weighting
eps <- sgsq0 * Vmat
Vmat <- abs(Vmat)

#y=9.81E-03*x-6.06E-02*x/1000*rnorm(8,0,1)
y <- beta * x - eps

dfw<-cbind(dfw,y)
while(k<30) {
        extra=NULL
        k=k+1
        #y=9.81E-03*x-6.06E-02*x/1000*rnorm(8,0,1)
        y <- beta * x - sgsq0 * rnorm(8,0,x)
        extra=data.frame(x)
        extra=cbind(extra,y)
        dfw=rbind(dfw,extra)
}
#dfw
par(mfrow = c(3, 2))
plot(dfw[,1],dfw[,2],col=1,pch=1,xlab="Concentration",ylab="Area",main="Linear with Weighting")

#-------------------------------------------------------------------------------------------

w <- rep(1/x, k) # create weight vector
fw <- lm(dfw$y ~ dfw$x, weights = w) # fit the weighted model

abline(fw, col = 'red') # plot regression line

fq <- lm(dfw$y ~ dfw$x^2 + dfw$x) # fit an unweighted quadratic model

fqw <- lm(dfw$y ~ dfw$x^2 + dfw$x, weights = w^2) # fit a weighted quadratic model

a <- anova(fw, fq, fqw) # compare the models with an ANOVA table

# --------------------------------------------------------------------------------

# From here down is for the weighted linear model

wRes <- fw$residuals * sqrt(w) # compute weighted residuals
dfw <- cbind(dfw, wRes)

wy <- dfw$y * sqrt(w) # weight the y values
dfw <- cbind(dfw, wy)

# Compute variance of y values at each cal point
yvar <- aggregate(dfw$y ~ dfw$x, dfw, var)
dfxvary <- data.frame(x, yvar)
# Compute variance of weighted y values at each cal point
wyvar <- aggregate(dfw$wy ~ dfw$x, dfw, var)

# Plot Var(y values) vs. unweighted x
plot(dfxvary$x, dfxvary$dfw.y, main = "Weighted 1/x Model", 
                xlab = "x", ylab = "Variance of Y")

# Plot Var(weighted y values) vs. x
plot(wyvar$'dfw$x', wyvar$`dfw$wy`, main = "Weighted 1/x Model", 
                xlab = "x", ylab = "Var. of Weighted Y Values")

# Plot Var(weighted y values) vs. weighted x
plot(wyvar$'dfw$x' * sqrt(w[1:length(x)]), wyvar$`dfw$wy`, main = "Weighted 1/x Model", 
     xlab = "Weighted x Values", ylab = "Var. of Weighted Y Values")

# Plot residuals vs weighted x values
plot(x = dfw$x * sqrt(w), y = wRes, main = "Weighted 1/x Model", xlab = "Weighted x",
     ylab = "Weighted Residuals")

# Plot residuals vs weighted fitted values        
plot(x = fw$fitted.values * sqrt(w), y = wRes, main = "Weighted 1/x Model", 
     xlab = "Weighted Fitted Values",
     ylab = "Weighted Residuals")
                