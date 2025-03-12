library(tableone)
library(Matching)
library(MatchIt)
library(ipw)
library(survey)

data(lalonde)
age <- as.numeric(lalonde$age)
educ <- as.numeric(lalonde$educ)
hispan <- as.numeric(lalonde$race=='hispan')
black <- as.numeric(lalonde$race=='black')
married <- as.numeric(lalonde$married)
nodegree <- as.numeric(lalonde$nodegree)
re74 <- as.numeric(lalonde$re74)
re75 <- as.numeric(lalonde$re75)
treat <- as.numeric(lalonde$treat)
re78 <- as.numeric(lalonde$re78)
treat <- as.numeric(lalonde$treat)
mydata <- cbind(age,
                educ,
                black,
                hispan,
                married,
                nodegree,
                re74,
                re75,
                treat,
                re78)
mydata <- data.frame(mydata)
myvars <- c("age","educ","black","hispan",
            "married", "nodegree","re74",
            "re75")
psmodel <- glm(treat~age+educ+black+hispan+married+nodegree+re74+re75,
                family=binomial(),
                data=mydata)

# value of propensity score
ps <- psmodel$fitted.values

# creat weights
weight <- ifelse(treat == 1, 1 / ps, 1 / (1 - ps))

print(summary(weight))

# Apply weights to data
weighteddata <- svydesign(ids = ~1, data = mydata, weights = ~weight)

weightedtable <- svyCreateTableOne(
    vars = myvars, strat = "treat",
    data = weighteddata, test = FALSE
)

print(weightedtable, smd = TRUE)

# Fit MSM
msm <- (svyglm(re78 ~ treat,
    design = svydesign(~1, weights = ~weight, data = mydata)
))
print(coef(msm))
print(confint(msm))

# Analysis with truncation
weightmodel <- ipwpoint(
    exposure = treat, family = "binomial",
    link = "logit",
    denominator = ~age+educ+black+hispan+married+nodegree+re74+re75,
    data = mydata,
    trunc = .01
)

print(summary(weightmodel$weights.trunc))

msm <- (svyglm(re78 ~ treat,
    design = svydesign(~1, weights = ~weightmodel$weights.trunc, data = mydata)
))
print(coef(msm))
print(confint(msm))