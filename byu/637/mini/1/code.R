library(lme4)       # glmer(), mixed effects (subject specific)
library(geepack)    # geeglm(), population average
data(respiratory)

head(respiratory)
respiratory$center = as.factor(respiratory$center)
respiratory$treat = ifelse(respiratory$treat == "P", "Placebo", "Active")
respiratory$sex = ifelse(respiratory$sex == "M", "Male", "Female")
respiratory$baseline = as.factor(respiratory$baseline)
respiratory$baseline = ifelse(respiratory$baseline == "0", "Poor", "Good")
respiratory$unique.id = respiratory$id
respiratory$unique.id[225:444] = respiratory$unique.id[225:444] + 56
head(respiratory)

### population average model
pa.mod = geeglm(outcome ~ treat + center + age + sex + baseline,
    data = respiratory, id = unique.id, family = binomial, corstr = "ar1")

### subject-specific model (uncorrelated observations?)
ss.mod = glmer(outcome ~ treat + center + age + sex + baseline + 
    (1|unique.id), data = respiratory, family = binomial)


# get random effects
b = ranef(ss.mod)[[1]][,1]

# "should" look normal
#plot(density(b))

# compare two models' coefficients
summary(pa.mod)
summary(ss.mod)
