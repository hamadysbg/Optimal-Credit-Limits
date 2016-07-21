#
#
#

require(quantreg)
require(boot)
require(rpart)

# Read Sample Data

sanofi=read.csv("Sanofi_CLA_06152016.csv") # Sanofi Sample
summary(sanofi[,-1])
dim(sanofi)

sanofi = sanofi[order(sanofi$totempl),]  # Order by Size Proxy


vaxs=read.csv("vaxserve_CLA_06152016.csv") # Vaxserve Sample
summary(vaxs[,-1])
dim(vaxs)

vaxs = vaxs[order(vaxs$totempl),]  # Order by Size Proxy


# Aggregate Industry Segments

levels(sanofi$industry_segment)[table(sanofi$industry_segment) < 100] = "OTHER"
levels(vaxs$industry_segment)[table(vaxs$industry_segment) < 100] = "OTHER"

# EDA

with(sanofi,hist(base_high_credit,breaks="FD"))
with(sanofi,summary(curr_high_credit / base_high_credit))
with(subset(sanofi,curr_high_credit < 5 * base_high_credit),hist(curr_high_credit / base_high_credit, breaks = "scott"))

with(sanofi,table(risk_cat,risk_group))

with(sanofi,table(risk_cat,industry_segment))
with(sanofi,plot(table(industry_segment,risk_cat==1),las=3))

with(sanofi,tapply(cons_cred_limit,risk_cat,min))
with(sanofi,tapply(cons_cred_limit,risk_cat,max))


boxplot(totempl~risk_cat,sanofi,log="y",col=4)
boxplot(totempl~risk_cat,vaxs,log="y",col=2)

plot(I(curr_high_credit+1) ~ tegroup,sanofi,log="y" , main = "SANOFI")
plot(I(curr_high_credit+1) ~ industry_segment,sanofi,log="y", las=2, cex.axis=0.7 , main = "SANOFI")
plot(I(curr_high_credit+1) ~ tegroup,vaxs,log="y" , main = "VAXSERVE")
plot(I(curr_high_credit+1) ~ industry_segment,vaxs,log="y", las=2, cex.axis=0.7 , main = "VAXSERVE")

q=quantile(sanofi$ccs, seq(0,1,0.1), type = 1)
names(q) = paste("R",seq(1,length(q)),sep="")
cbind(Lower=q[-length(q)],Upper=q[-1])

# Statistical Analysis

# 90th Percentile Regression
s.fit = rq(base_high_credit ~ ccs + fss + factor(risk_cat) + totempl, tau = 0.90, data = sanofi)
summary(s.fit)
# fss and ccs and risk_cat not significant

s1.fit = rq(curr_high_credit ~ ccs + fss + base_high_credit, tau = 0.90, data = sanofi)
summary(s1.fit)
s1b.fit = rq(curr_high_credit ~ fss + base_high_credit, tau = 0.90, data = sanofi)
summary(s1b.fit)
anova(s1.fit,s1b.fit)

s2.fit = rq(curr_high_credit ~ ccs + fss + totempl, tau = 0.90, data = sanofi)
summary(s2.fit)
s2b.fit = rq(curr_high_credit ~ fss + totempl, tau = 0.90, data = sanofi)
summary(s2b.fit)
anova(s2.fit,s2b.fit)

AIC(s1b.fit)
AIC(s2b.fit)


# Apply Log Transformation
sb.fit = rq(log(curr_high_credit+1) ~ log(base_high_credit+1) + log(totempl), tau = 0.75, data = sanofi)
summary(sb.fit)
sc.fit = rq(log(curr_high_credit+1) ~ log(base_high_credit+1) , tau = 0.75, data = sanofi)
summary(sc.fit)
anova(sb.fit,sc.fit)
AIC(sb.fit)
AIC(sc.fit)

# Relation Between curr_high_credit and base_high_credit

plot(curr_high_credit ~ base_high_credit, sanofi, log="xy")

rpart(log(curr_high_credit+1) ~ . , sanofi, control = list(xref=10))

summary(lm(log(curr_high_credit+1) ~ log(base_high_credit) , sanofi))

###############
# FINAL MODEL #
###############

# Multiple Percentile Specifications
s2b.fit = rq(log(curr_high_credit+1) ~ log(totempl) + sqrt(totempl) + industry_segment, tau = c(0.50, 0.75, 0.90), data = sanofi)
summary(s2b.fit)
anova(s2b.fit)

v2b.fit = rq(log(curr_high_credit+1) ~ log(totempl) + sqrt(totempl) + industry_segment, tau = c(0.50, 0.75, 0.90), data = vaxs)
summary(v2b.fit)
anova(v2b.fit,joint=FALSE)

# Test Fit
sanofi$pred   = exp(predict(s2b.fit))-1
vaxs$pred     = exp(predict(v2b.fit))-1

tau = c(0.50, 0.75, 0.90)

with(sanofi,tapply(curr_high_credit,cut(totempl,c(0,4,9,23,Inf)),quantile,tau))

plot(I(curr_high_credit+1)~totempl,sanofi,log="xy")
with(sanofi,lines(lowess(curr_high_credit),col=2))


sanofi$tegroup = with(sanofi,cut(totempl,quantile(totempl,seq(0,1,0.1), type = 1), include.lowest = T) )
vaxs$tegroup   = with(vaxs ,cut(totempl,quantile(totempl,seq(0,1,0.2), type = 1), include.lowest = T) )

# SANOFI RESULTS
with(sanofi,plot(tapply(curr_high_credit,tegroup,quantile,0.50),type="b",log="y"))
with(sanofi,lines(tapply(pred[,1],tegroup,median),col=2))

with(sanofi,plot(tapply(curr_high_credit,tegroup,quantile,0.75),type="b",log="y"))
with(sanofi,lines(tapply(pred[,2],tegroup,median),col=2))

with(sanofi,plot(tapply(curr_high_credit,tegroup,quantile,0.90),type="b",log="y"))
with(sanofi,lines(tapply(pred[,3],tegroup,median),col=2))


plot(sanofi$curr_high_credit+1,log="y",ylab = "Actual")
points(sanofi$pred[,3],col=2)

# VAXSERVE RESULTS
with(vaxs,plot(tapply(curr_high_credit,tegroup,quantile,0.50),type="b"))
with(vaxs,lines(tapply(pred[,1],tegroup,median),col=2))

with(vaxs,plot(tapply(curr_high_credit[order(totempl)],tegroup,quantile,0.75),type="b"))
with(vaxs,lines(tapply(pred[,2],tegroup,median),col=2))

with(vaxs,plot(tapply(curr_high_credit[order(totempl)],tegroup,quantile,0.90),type="b"))
with(vaxs,lines(tapply(pred[,3],tegroup,median),col=2))


plot(vaxs$curr_high_credit+1,log="y",ylab = "Actual")
points(vaxs$pred[,3],col=2)


# How many clients are below Model Credit Limit? 
with(sanofi,table(curr_high_credit <= pred[,3]))


# How many clients are below Model Credit Limit? 
with(vaxs,table(curr_high_credit <= pred[,3]))


############################
# Export Example to GAMS

sanofi = sanofi[order(sanofi$ccs),]  # Order by CCS Pct
set.seed(2016)

ex=read.csv("example_bad_v_ccs.csv")

sanofi$bad = ex$bad # Use Example Bad Values
mean(sanofi$bad)

sanofi$pd  = round(with(sanofi, inv.logit(0.05 * (12 - ccs)) ),5) # Use Example PD Values
mean(sanofi$pd)
with(sanofi,tapply(pd,bad,mean))

quantile(sanofi$pd,seq(0,1,0.1))
sanofi$pdgroup = with(sanofi,cut(pd,quantile(sanofi$pd,seq(0,1,0.1)),include.lowest = T))
with(sanofi,table(pdgroup,tegroup))

#sanofi$lgd = round(rbeta(1000,2,5),5)
sanofi$lgd = with(sanofi, 0.95 - 0.08 * as.numeric(tegroup) )
hist(sanofi$lgd)

#write.csv(cbind(name=paste("C",1:nrow(sanofi),sep=""),sanofi[,c("totempl","base_high_credit","ccs","pd","lgd")]),
#          "toGams.csv",row.names = F)
write.csv(cbind(name=paste("C",1:nrow(sanofi),sep=""),sanofi[,c("totempl","base_high_credit","ccs","pd","lgd")]),
          "toGams2.csv",row.names = F)

summary(sanofi$lgd)
summary(sanofi$pd)

with(sanofi,tapply(base_high_credit*pd*lgd,list(pdgroup,tegroup),mean))
with(sanofi,tapply(pd*lgd,list(pdgroup,tegroup),mean))
