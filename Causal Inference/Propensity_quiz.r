library(tableone)
library(Matching)
library(MatchIt)

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
table1 <- CreateTableOne(data = mydata, 
                    vars = myvars,
                    strata = "treat",
                    test = FALSE)
print(table1,smd=TRUE)
print(mean(mydata$re78[mydata$treat==1])-mean(mydata$re78[mydata$treat==0]))

psmodel <- glm(treat~age+educ+black+hispan+married+nodegree+re74+re75,
                family=binomial(),
                data=mydata)

print(summary(psmodel))
pscore<-psmodel$fitted.values

print(min(pscore))
print(max(pscore))

set.seed(931139)

psmatch <- Match(Tr=mydata$treat,M=1,X=pscore,replace=FALSE,caliper = .1)
matched <- mydata[unlist(psmatch[c("index.treated","index.control")]),]
matchedtab1 <- CreateTableOne(vars=myvars,
                              data=matched,
                              test=FALSE,
                              strata="treat")
print(matchedtab1,smd=TRUE)

print(mean(matched$re78[matched$treat==1])-mean(matched$re78[matched$treat==0]))

y_trt <- matched$re78[matched$treat==1]
y_con <- matched$re78[matched$treat==0]

diffy<-y_trt-y_con

print(t.test(diffy))