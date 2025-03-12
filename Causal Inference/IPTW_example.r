library(tableone)
library(ipw)
library(sandwich)
library(survey)

rhc <- read.csv("rhc.csv")
# View(rhc)

arf <- as.numeric(rhc$cat1 == "ARF")
chf <- as.numeric(rhc$cat1 == "CHF")
cirr <- as.numeric(rhc$cat1 == "Cirrhosis")
colcan <- as.numeric(rhc$cat1 == "Colon Cancer")
coma <- as.numeric(rhc$cat1 == "Coma")
copd <- as.numeric(rhc$cat1 == "COPD")
lungcan <- as.numeric(rhc$cat1 == "Lung Cancer")
mosf <- as.numeric(rhc$cat1 == "MOSF w/Malignancy")
sepsis <- as.numeric(rhc$cat1 == "MOSF w/Sepsis")
female <- as.numeric(rhc$sex == "Female")
died <- as.numeric(rhc$death == "Yes")
age <- rhc$age
treatment <- as.numeric(rhc$swang1 == "RHC")
meanbp1 <- rhc$meanbp1
aps <- rhc$aps1

# new dataset
xvars <- c(
    "arf", "chf", "cirr", "colcan", "coma",
    "lungcan", "mosf", "sepsis", "age",
    "female", "meanbp1"
)

mydata <- cbind(
    arf, chf, cirr, colcan, coma,
    lungcan, mosf, sepsis, age,
    female, meanbp1, aps, treatment, died
)
mydata <- data.frame(mydata)

# propensity score model
psmodel <- glm(
    treatment ~ age + female + meanbp1 + arf + chf +
        cirr + colcan + coma + lungcan + mosf + sepsis,
    family = binomial(link = "logit")
)

# value of propensity score
ps <- predict(psmodel, type = "response")

print(summary(psmodel))

# creat weights
weight <- ifelse(treatment == 1, 1 / ps, 1 / (1 - ps))

# Apply weights to data
weighteddata <- svydesign(ids = ~1, data = mydata, weights = ~weight)

weightedtable <- svyCreateTableOne(
    vars = xvars, strat = "treatment",
    data = weighteddata, test = FALSE
)

print(weightedtable, smd = TRUE)

# Manual weighted mean
print(mean(weight[treatment == 1] * age[treatment == 1]) / mean(weight[treatment == 1]))

# Marginal structural models
glm.obj <- glm(died ~ treatment, weights = weight, family = binomial(link = "identity"))
print(summary(glm.obj))

betaiptw <- coef(glm.obj)

# To account for weighting, use sandwich variance
SE <- sqrt(diag(vcovHC(glm.obj, type = "HC0")))

# get point estimate and CI for relative risk ( need to exponentiate)
# causal relative risk
causalrr <- (betaiptw[2])
lcl <- (betaiptw[2] - 1.96 * SE[2])
ucl <- (betaiptw[2] + 1.96 * SE[2])

print(c(lcl, causalrr, ucl))

weightmodel <- ipwpoint(
    exposure = treatment, family = "binomial",
    link = "logit",
    denominator = ~ age + female + meanbp1 + arf + chf +
        cirr + colcan + coma + lungcan + mosf + sepsis,
    data = mydata
)

print(summary(weightmodel$weights.trun))
# ipwplot(weights = weightmodel$ipw.weights, logscale = FALSE,
# main = "weights",xlim = c(0,22))

# Fit MSM
msm <- (svyglm(died ~ treatment,
    design = svydesign(~1, weights = ~weight, data = mydata)
))
print(coef(msm))
print(confint(msm))

# truncating weights
truncweight <- replace(weight, weight > 10, 10)

# get causal risk difference
glm.obj <- glm(died ~ treatment,
    weights = truncweight,
    family = binomial(link = "identity")
)

weightmodel <- ipwpoint(
    exposure = treatment, family = "binomial",
    link = "logit",
    denominator = ~ age + female + meanbp1 + arf + chf +
        cirr + colcan + coma + lungcan + mosf + sepsis,
    data = mydata,
    trunc = .01
)

print(summary(weightmodel$ipw.weights))
# ipwplot(weights = weightmodel$ipw.weights, logscale = FALSE,
# main = "weights",xlim = c(0,22))

# Fit MSM
msm <- (svyglm(died ~ treatment, design = svydesign(~1, weights = ~weight, data = mydata,trunc=0.01)))
print(coef(msm))
print(confint(msm))
