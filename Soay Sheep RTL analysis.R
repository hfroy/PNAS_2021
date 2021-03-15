rm(list=ls(all=TRUE))
library(glmmTMB)
library(MCMCglmm)
library(MasterBayes)

RTL <- read.table('SoaySheepRTL.txt', sep='\t', header=T)
nrow(RTL) # 3641
Ped <- read.table('PrunedPed.txt', sep='\t', header=T)
nrow(Ped) # 2411


##############################
### AGE FUNCTION
##############################

### THRESHOLD FUNCTION
threshold1 <- function(Age, T1)
{ Age.1 <- Age.2 <- Age
Age.1[Age.1 > T1] <- T1
Age.2[Age.2 <= T1] <- 0
Age.2[Age.2 > T1] <- Age.2[Age.2 > T1] - T1
cbind(Age1.1 = Age.1, Age1.2 = Age.2) }

temp <- RTL
A1 <- glmmTMB (RTL ~ Adult + threshold1(Age, 2) + (1|SampleYear) + (1|ID) + (1|qPCRPlate) + (1|qPCRRow), data=temp, REML=F)
A2 <- glmmTMB (RTL ~ Adult + threshold1(Age, 3) + (1|SampleYear) + (1|ID) + (1|qPCRPlate) + (1|qPCRRow), data=temp, REML=F)
A3 <- glmmTMB (RTL ~ Adult + threshold1(Age, 4) + (1|SampleYear) + (1|ID) + (1|qPCRPlate) + (1|qPCRRow), data=temp, REML=F)
A4 <- glmmTMB (RTL ~ Adult + threshold1(Age, 5) + (1|SampleYear) + (1|ID) + (1|qPCRPlate) + (1|qPCRRow), data=temp, REML=F)
A5 <- glmmTMB (RTL ~ Adult + threshold1(Age, 6) + (1|SampleYear) + (1|ID) + (1|qPCRPlate) + (1|qPCRRow), data=temp, REML=F)
A6 <- glmmTMB (RTL ~ Adult + threshold1(Age, 7) + (1|SampleYear) + (1|ID) + (1|qPCRPlate) + (1|qPCRRow), data=temp, REML=F)
A7 <- glmmTMB (RTL ~ Adult + threshold1(Age, 8) + (1|SampleYear) + (1|ID) + (1|qPCRPlate) + (1|qPCRRow), data=temp, REML=F)
A8 <- glmmTMB (RTL ~ Adult + threshold1(Age, 9) + (1|SampleYear) + (1|ID) + (1|qPCRPlate) + (1|qPCRRow), data=temp, REML=F)
A9 <- glmmTMB (RTL ~ Adult + threshold1(Age, 10) + (1|SampleYear) + (1|ID) + (1|qPCRPlate) + (1|qPCRRow), data=temp, REML=F)
A10 <- glmmTMB (RTL ~ Adult + threshold1(Age, 11) + (1|SampleYear) + (1|ID) + (1|qPCRPlate) + (1|qPCRRow), data=temp, REML=F)
A11 <- glmmTMB (RTL ~ Adult + poly(Age, 2) + (1|SampleYear) + (1|ID) + (1|qPCRPlate) + (1|qPCRRow), data=temp, REML=F)
A12 <- glmmTMB (RTL ~ Adult + poly(Age, 3) + (1|SampleYear) + (1|ID) + (1|qPCRPlate) + (1|qPCRRow), data=temp, REML=F)
A13 <- glmmTMB (RTL ~ Adult + Age + (1|SampleYear) + (1|ID) + (1|qPCRPlate) + (1|qPCRRow), data=temp, REML=F)
A14 <- glmmTMB (RTL ~ Adult + (1|SampleYear) + (1|ID) + (1|qPCRPlate) + (1|qPCRRow), data=temp, REML=F)
A15 <- glmmTMB (RTL ~ Age + (1|SampleYear) + (1|ID) + (1|qPCRPlate) + (1|qPCRRow), data=temp, REML=F)
A16 <- glmmTMB (RTL ~ 1 + (1|SampleYear) + (1|ID) + (1|qPCRPlate) + (1|qPCRRow), data=temp, REML=F)

AICs <- AIC(A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12, A13, A14, A15, A16)
new <- AICs
new$dAIC <- new$AIC - min(new$AIC)
new <- new[order(new$AIC),]

### BEST MODEL
priorRTL <- list(G=list(G1=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*100), 
                        G2=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*100),
                        G3=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*100),
                        G4=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*100)),
                 R=list(V=diag(1), nu=0.002))

RTL_Age <- MCMCglmm(RTL ~ Adult + Age, 
                    random = ~ ID + SampleYear + qPCRPlate + qPCRRow, 
                    family = "gaussian", data = temp, prior = priorRTL, 
                    nitt = 110000, burnin = 10000, thin = 50)


##############################
### SEX EFFECTS
##############################

temp <- RTL
A1 <- glmmTMB (RTL ~ Adult + Age + (1|SampleYear) + (1|ID) + (1|qPCRPlate) + (1|qPCRRow), data=temp, REML=F)
A1a <- glmmTMB (RTL ~ Adult + Age + Sex + (1|SampleYear) + (1|ID) + (1|qPCRPlate) + (1|qPCRRow), data=temp, REML=F)
A1b <- glmmTMB (RTL ~ Adult * Sex + Age + (1|SampleYear) + (1|ID) + (1|qPCRPlate) + (1|qPCRRow), data=temp, REML=F)
A1c <- glmmTMB (RTL ~ Adult + Age * Sex + (1|SampleYear) + (1|ID) + (1|qPCRPlate) + (1|qPCRRow), data=temp, REML=F)

anova(A1, A1a)
anova(A1, A1b)
anova(A1, A1c)


##############################
### HERITABILITY
##############################

temp <- subset(RTL, !is.na(MumID))
temp$ANIMAL <- as.factor(temp$ID)
temp$MOTHER <- as.factor(temp$MumID)
Ainv <- inverseA(Ped)$Ainv

### UNIVARIATE
priorRTL <- list(G=list(G1=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*100), 
                        G2=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*100),
                        G3=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*100),
                        G4=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*100),
                        G5=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*100),
                        G6=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*100)),
                 R=list(V=diag(1), nu=0.002))

RTLAni <- MCMCglmm(RTL ~ Adult + Sex + Age, 
                      random = ~ ANIMAL + ID + MOTHER + SampleYear + qPCRPlate + qPCRRow, 
                      ginverse=list(ANIMAL=Ainv), family = "gaussian", 
                      data = temp, prior = priorRTL, 
                      nitt = 110000, burnin = 10000, thin = 50)

summary(RTLAni)

### BIVARIATE
temp <- subset(RTL, !is.na(MumID))
temp$ANIMAL <- as.factor(temp$ID)
temp$MOTHER <- as.factor(temp$MumID)
Ainv <- inverseA(Ped)$Ainv

priorRTL <- list(G=list(G1=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*100), 
                        G2=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*100),
                        G3=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*100),
                        G4=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*100),
                        G5=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*100),
                        G6=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*100)),
                 R=list(V=diag(2), nu=2.002))

BiAni <- MCMCglmm(RTL ~ Adult - 1 + at.level(Adult, 1):(Sex) +
                    at.level(Adult, 2):(Sex + Age), 
                  random = ~ us(Adult):ANIMAL + us(Adult):SampleYear +
                    us(at.level(Adult, 1)):MOTHER + us(at.level(Adult, 2)):ID +
                    qPCRPlate + qPCRRow, 
                  rcov = ~us(Adult):units, 
                  ginverse=list(ANIMAL=Ainv), family = "gaussian", 
                  data = temp, prior = priorRTL, 
                  nitt = 210000, burnin = 10000, thin = 200)

summary(BiAni)


##############################
### MULTIVARIATE 
##############################

temp <- subset(RTL, !is.na(Survival)  & !is.na(Weight) & !is.na(MumID))
temp$RTL <- scale(temp$RTL)
temp$Weight <- scale(temp$Weight)

### PHENOTYPIC
priorRTL <- list(G=list(G1=list(V=diag(3), nu=3, alpha.mu=c(0,0,0), alpha.V=diag(3)*100), 
                        G2=list(V=diag(3), nu=3, alpha.mu=c(0,0,0), alpha.V=diag(3)*100),
                        G3=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*100),
                        G4=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*100),
                        G5=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*100),
                        G6=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*100),
                        G7=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*100),
                        G8=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*100)), 
                 R=list(V=diag(3), nu=3.002, fix=3))

Tri_pheno <- MCMCglmm(cbind(RTL, Weight, Survival) ~ trait - 1  + 
                        at.level(trait, "RTL"):(Adult) +
                        at.level(trait, "Weight"):(Adult * Sex + Sex * poly(Age, 3)) +
                        at.level(trait, "Survival"):(Adult * Sex + poly(Age, 2)), 
                      random = ~ us(trait):ID + us(trait):SampleYear + 
                        us(at.level(trait, "RTL")):qPCRPlate + us(at.level(trait, "RTL")):qPCRRow +
                        us(at.level(trait, "Weight")):MumID + us(at.level(trait, "Weight")):BirthYear +
                        us(at.level(trait, "Survival")):MumID + us(at.level(trait, "Survival")):BirthYear,
                      rcov = ~us(trait):units, family = c("gaussian", "gaussian", "threshold"),
                      trunc = T, data = temp, prior = priorRTL, 
                      nitt = 530000, burnin = 30000, thin = 500)

summary(Tri_pheno)

### GENETIC
temp <- subset(RTL, !is.na(Survival)  & !is.na(Weight) & !is.na(MumID))
temp$RTL <- scale(temp$RTL)
temp$Weight <- scale(temp$Weight)
temp$ANIMAL <- as.factor(temp$ID)
temp$MOTHER <- as.factor(temp$MumID)
Ainv <- inverseA(Ped)$Ainv

priorRTL <- list(G=list(G1=list(V=diag(3), nu=3, alpha.mu=c(0,0,0), alpha.V=diag(3)*100), 
                        G2=list(V=diag(3), nu=3, alpha.mu=c(0,0,0), alpha.V=diag(3)*100),
                        G3=list(V=diag(3), nu=3, alpha.mu=c(0,0,0), alpha.V=diag(3)*100),
                        G4=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*100),
                        G5=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*100),
                        G6=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*100),
                        G7=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*100),
                        G8=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*100),
                        G9=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*100)), 
                 R=list(V=diag(3), nu=3.002, fix=3))

Tri_ani <- MCMCglmm(cbind(RTL, Weight, Survival) ~ trait - 1 + 
                      at.level(trait, "RTL"):(Adult) + 
                      at.level(trait, "Weight"):(Adult * Sex + Sex * poly(Age, 3)) +
                      at.level(trait, "Survival"):(Adult * Sex + poly(Age, 2)), 
                    random = ~ us(trait):ANIMAL + us(trait):ID + us(trait):SampleYear + 
                      us(at.level(trait, "RTL")):qPCRPlate + us(at.level(trait, "RTL")):qPCRRow +
                      us(at.level(trait, "Weight")):MumID + us(at.level(trait, "Weight")):BirthYear +
                      us(at.level(trait, "Survival")):MumID + us(at.level(trait, "Survival")):BirthYear,
                    rcov = ~us(trait):units, ginverse=list(ANIMAL=Ainv),
                    family = c("gaussian", "gaussian", "threshold"),
                    trunc = T, data = temp, prior = priorRTL, 
                    nitt = 530000, burnin = 30000, thin = 500)

summary(Tri_ani)

