#Math 156 Final Project Script 

#summary
mydata <- read.csv("modAGdata.csv") #read in data
colnames(mydata)
head(mydata) 
nrow(mydata)

#------
#BMI & Gamma distribution 
#Bootstrapping to approximate sampling distribution -- Gamma distribution?
  #shape = (mean/sd)^2
  #rate = mean/sd^2
BMI2 <- mydata$BMI
par(mfrow=c(1,1))
hist(BMI2,  freq = FALSE, main = "Histogram of BMI")
mean1 <- mean(BMI2); mean1
sd1 <- sd(BMI2); sd1
shape1 <- (mean1/sd1)^2; shape1
rate1 <- mean1/sd1^2; rate1
#overlay gamma distribution
curve(dgamma(x, shape1, rate1), col = "red", add = TRUE)  

#bootstrapping to estimate sampling distribution 
n <- length(BMI2)
N <- 10^4
BMI2.mean <- numeric(N)
for (i in 1:N)
{x <- sample(BMI2, n, replace = TRUE)
BMI2.mean[i] <- mean(x)
}
hist(BMI2.mean, freq = FALSE, main = "Bootstrap distribution of means")
abline(v = mean(BMI2), col = "blue", lty = 2) #observed means

mean2 <- mean(BMI2.mean); mean2
sd2 <- sd(BMI2.mean); sd2
shape2 <- (mean2/sd2)^2; shape2
rate2 <- mean2/sd2^2; rate2
#overlay gamma distribution
curve(dgamma(x, shape2, rate2), col = "red", add = TRUE)

#bootstrap distribution of means approximately gamma 

#____________

#regressions
#logistic regression - BMI and Bulimia (0 or 1)
par(mfrow=c(1,1))
BMI <- mydata$BMI
bulimia <- as.integer(mydata$Bulimia)
plot(BMI, bulimia, main = "Bulimia vs. BMI")

library(stats4)
MLL<- function(alpha, beta)
  -sum( log( exp(alpha+beta*BMI)/(1+exp(alpha+beta*bulimia)) )*bulimia+ log(1/(1+exp(alpha+beta*BMI)))*(1-bulimia) )
results<-mle(MLL,start = list(alpha = -0.1, beta = -0.02))
results@coef
curve(exp(results@coef[1]+results@coef[2]*x)/ (1+exp(results@coef[1]+results@coef[2]*x)), col = "blue", xlab = "BMI", ylab = "Log-odds ratio of Bulimia", main = "Log-odds ratio of Bulimia vs. BMI")

#a one unit increase in BMI will result in a 0.021 decrease in the log-odds ratio of success:failure
#in other words, if we increase BMI, the odds of Bulimia versus no bulimia will decrease
#resulting in bulimia being less likely than it was before the increase

#linear regression - Age and BMI 
age <- mydata$Age
plot(age, BMI, main = "BMI vs. Age")
max(age) #there's no way someone is 146 years old, remove incorrect outliers
mydata2 <- mydata[which(mydata$Age < 100), ]
age2 <- mydata2$Age; max(age2)
BMI2 <- mydata2$BMI
plot(age2, BMI2, main = "BMI vs. Age: w/o outliers")

lm <- lm(BMI2~age2, data = mydata2);lm
a <- lm$coefficients[1]; a
b <- lm$coefficients[2]; b

#We can add this regression line to the plot of the data
abline(a, b, col = "red")
#a one year increase in age increases BMI by approximately .09 units

#we can split the observations into predicted values and residuals
PredictBMI = a + b * age2
points(age2,PredictBMI, col = "green",add = TRUE, pch = 20)  #these lie on the regression line
#The residual is the observation minus the prediction
ResidBMI = BMI2 - PredictBMI

#We can plot the residuals separately
par(mfrow=c(1,1))
plot(age2, BMI2)      #for comparison
abline(a, b, col = "red")
points(age2,PredictBMI, col = "green",add = TRUE, pch = 20, main = "Predicted BMI")

plot(age2,ResidBMI, main = "Residuals")    #by construction, the residuals have a mean of zero
abline(h=0, col = "blue")

#Weak positive relationship between age and BMI which is what we might expect
#probably want to look at time series for individual patients for further affirmation
#and in order to rule out time specific factors (do older people have higher BMIs than younger people for reason other than age?)

#_______________

#Student T distribution - difference in Bacteroidaceae between males and females
par(mfrow=c(1,1))
males <- subset(mydata, select = Bacteroidaceae, subset = (Sex == "FALSE"), drop = TRUE) #subset males in Bacteroidaceae
females <- subset(mydata, select = Bacteroidaceae, subset = (Sex == "TRUE"), drop = TRUE) #subset females Bacteroidaceae
Mm<-mean(males); Mf<-mean(females); Mm; Mf #is this difference significant?
nm <- length(males); nf <- length(females); nm; nf #number of males and females
Sm<- sd(males); Sf<- sd(females); Sm; Sf #calculate sample standard deviations

#calculate sample variance - this is the appropriate student denominator 
SE <- sqrt(Sm^2/nm + Sf^2/nf); SE 
#Welch's approximation- reasonable value for degrees of freedom
ndf <- (Sm^2/nm + Sf^2/nf)^2/(Sm^4/(nm^2*(nm-1)) + Sf^4/(nf^2*(nf-1))); ndf 

#Assume that the difference in means divided by the standard error is Student t
#Get the quantiles for the t distribution with this many degrees of freedom
tqnt<- qt(c(.025, .975),ndf); tqnt

#Multiply by the sample standard error to create a 95% confidence interval
interval<- (Mm-Mf)+ SE*tqnt; interval #large but does not include zero

#automated test - same results
t.test(males, females)
interval 

#Use Student t to compute P value 
#calculate the probability that such a large rescaled mean could arise by chance
curve(dt(x, df  = ndf), from = -5, to =5)
abline(v = (Mm-Mf)/SE, col = "red" )
pval <-pt((Mm-Mf)/SE, df = ndf, lower.tail = FALSE); pval
2*pval    #double for two-sided

#p-value and confidence interval are inconsistent- what's the problem?
#population is not normal 
hist(mydata$Bacteroidaceae) #by inspection 
library(dplyr)
sample <- sample_n(mydata, 1000)
shapiro.test(sample$Bacteroidaceae) #Shapiro Wilkinson test for normality
#not consistent with normal distribution which is a condition for the Welch's two sample t test
#lesson - don't assume normal distributions in real world data

#Instead --> Bootstrapping difference of means
N <- 10^4
sex.diff.mean <- numeric(N)
for (i in 1:N)
{
  male.sample <- sample(males, nm, replace = TRUE)
  female.sample <- sample(females, nf, replace = TRUE)
  sex.diff.mean[i] <- mean(male.sample) - mean(female.sample)
}
hist(sex.diff.mean,
     main = "Bootstrap distribution of difference in means")
abline(v = mean(males) - mean(females), col = "blue", lty = 2)

#numeric summaries 
mean(males) - mean(females) #observed diff in means
mean(sex.diff.mean) #bootstrapped diff in means
sd(sex.diff.mean)
quantile(sex.diff.mean, c(0.025, 0.975)) 
mean(sex.diff.mean) - (mean(males) - mean(females)) #bias

#confidence interval display 
counter <- 0
plot(x =c(-2500,1000), y = c(1,105), type = "n", xlab = "Mean difference between sex", ylab = "Frequency", main = "Random Confidence Intervals") #blank plot
observed.mean <- mean(males) - mean(females)
mean(x); var(x)
for (i in 1:1000) {
  x <- sample(sex.diff.mean, 30, replace = TRUE)  #random sample of diff in mean
  L <- quantile(x, 0.025) #usually less than the true mean
  U <- quantile(x, 0.975) #usually greater than the true mean
  if (L < observed.mean && U > observed.mean) counter <- counter + 1 #count +1 if we were correct
  if(i <= 100) segments(L, i, U, i)
}
abline (v = observed.mean, col = "red") #vertical line at true mean

#seems that females have higher counts of this type of bacteria than males 



