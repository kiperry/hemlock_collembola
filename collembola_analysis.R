############################################################################
# Non-target impacts of insecticides used to treat elongate hemlock scale
# on soil arthropods
#
# Cleveland Metroparks
# North Chagrin Reservation
# AB Williams Woods
#
# Analyses for abundance and family richness of Collembola
#
# Alyssa Mills & Kayla I Perry
#
# 31 October 2025
############################################################################

# import data
dat <- read.csv("collembola_counts.csv")
str(dat)
dat$Treatment <- as.factor(dat$Treatment)
dat$Block <- as.factor(dat$Block)
dat$Tree_Rep <- as.factor(dat$Tree_Rep)
str(dat)

colSums(dat[,8:14])

boxplot(Onychiuridea ~ Treatment, data = dat)
boxplot(Isotomidea ~ Treatment, data = dat)
boxplot(Entomobryidae ~ Treatment, data = dat)
boxplot(Sminthuridae ~ Treatment, data = dat)
boxplot(Tomoceridae ~ Treatment, data = dat)
boxplot(Hypogastruridae ~ Treatment, data = dat)
boxplot(Neanuridae ~ Treatment, data = dat)

# load packages
library(tidyverse)
library(lme4)
library(lmerTest)
library(blmeco)
library(emmeans)
library(car)
library(DHARMa)
library(hillR)
library(RColorBrewer)

# separate data for analyses
con <- dat[which(dat$Treatment == "Control"),]
maple <- dat[which(dat$Treatment == "Sugar_Maple"),]
TR22 <- dat[which(dat$Treatment == "TR_2022"),]
TR23 <- dat[which(dat$Treatment == "TR_2023"),]
TRwk <- dat[which(dat$Treatment == "TR_2024_1wk"),]
TRpr <- dat[which(dat$Treatment == "TR_2024_pre"),]
TRmo <- dat[which(dat$Treatment == "TR_2024_1mo"),]

# across year treatment comparison: control, 2022, 2023, 2024 (week or month)
# within year treatment comparison: pre, week, month
# species comparison: maple, control

trmt.w.year <- rbind(TRpr, TRwk, TRmo) %>% droplevels()
str(trmt.w.year)

trmt.w.year$abund <- rowSums(trmt.w.year[,8:14])
trmt.w.year$rich <- hill_taxa(trmt.w.year[,8:14], q = 0, MARGIN = 1)
str(trmt.w.year)

mapl <- rbind(maple, con) %>% droplevels()
str(mapl)

mapl$abund <- rowSums(mapl[,8:14])
mapl$rich <- hill_taxa(mapl[,8:14], q = 0, MARGIN = 1)

#
trmt.a.year <- rbind(con, TR22, TR23, TRmo) %>% droplevels()
str(trmt.a.year)

trmt.a.year$abund <- rowSums(trmt.a.year[,8:14])
trmt.a.year$rich <- hill_taxa(trmt.a.year[,8:14], q = 0, MARGIN = 1)
str(trmt.a.year)


### models: within year
dotchart(trmt.w.year$abund, group = trmt.w.year$Treatment)
hist(trmt.w.year$abund)
boxplot(trmt.w.year$abund ~ trmt.w.year$Treatment)

dotchart(trmt.w.year$rich, group = trmt.w.year$Treatment)
hist(trmt.w.year$rich)
boxplot(trmt.w.year$rich ~ trmt.w.year$Treatment)

colSums(trmt.w.year[,8:14])

dotchart(trmt.w.year$Onychiuridea, group = trmt.w.year$Treatment)
hist(trmt.w.year$Onychiuridea)
boxplot(trmt.w.year$Onychiuridea ~ trmt.w.year$Treatment)

dotchart(trmt.w.year$Isotomidea, group = trmt.w.year$Treatment)
hist(trmt.w.year$Isotomidea)
boxplot(trmt.w.year$Isotomidea ~ trmt.w.year$Treatment)

boxplot(trmt.w.year$Entomobryidae ~ trmt.w.year$Treatment)
boxplot(trmt.w.year$Hypogastruridae ~ trmt.w.year$Treatment)

#
abund.w.year <- glmer(abund ~ Treatment + (Tree_No|Block), 
                      family = poisson, data = trmt.w.year)
summary(abund.w.year)
Anova(abund.w.year, type = "III")
emmeans(abund.w.year, pairwise ~ Treatment)
testDispersion(abund.w.year)
res.mod.1 <- simulateResiduals(abund.w.year)
plot(abund.w.year)
testCategorical(abund.w.year, catPred = trmt.w.year$Treatment)
testZeroInflation(abund.w.year)

#
rich.w.year <- glmer(rich ~ Treatment + (Tree_No|Block), 
                      family = poisson, data = trmt.w.year)
summary(rich.w.year)
Anova(rich.w.year, type = "III")
emmeans(rich.w.year, pairwise ~ Treatment)
testDispersion(rich.w.year)
res.mod.1 <- simulateResiduals(rich.w.year)
plot(rich.w.year)
testCategorical(rich.w.year, catPred = trmt.w.year$Treatment)
testZeroInflation(rich.w.year)


# change the names of the variables and reorder them for the figures
levels(trmt.w.year$Treatment)
levels(trmt.w.year$Treatment)[levels(trmt.w.year$Treatment)=="TR_2024_1mo"] <- "Post-Month"
levels(trmt.w.year$Treatment)[levels(trmt.w.year$Treatment)=="TR_2024_1wk"] <- "Post-Week"
levels(trmt.w.year$Treatment)[levels(trmt.w.year$Treatment)=="TR_2024_pre"] <- "Pre-treatment"
levels(trmt.w.year$Treatment)
trmt.w.year$Treatment <- factor(trmt.w.year$Treatment, levels = c("Pre-treatment", "Post-Week", "Post-Month"))

par(mfrow=c(1,2))
par(mar=c(5,7,4,2))

boxplot(abund ~ Treatment, data = trmt.w.year, col = c("seagreen4", "seagreen2", "goldenrod2"),
        ylim = c(0,8), ylab = "Abundance", xlab = "", cex.lab = 1.6, cex.axis = 1.5)
stripchart(abund ~ Treatment, data = trmt.w.year, pch = 19, cex = 2, add = TRUE, 
           vertical = TRUE, method = "jitter", jitter = 0.2)
text(0.5,7.4, "A", pos = 3, font = 1, cex = 1.5)

boxplot(rich ~ Treatment, data = trmt.w.year, col = c("seagreen4", "seagreen2", "goldenrod2"),
        ylim = c(0,6), ylab = "Family Richness", xlab = "", cex.lab = 1.6, cex.axis = 1.5)
stripchart(rich ~ Treatment, data = trmt.w.year, pch = 19, cex = 2, add = TRUE, 
           vertical = TRUE, method = "jitter", jitter = 0.2)
text(0.5,5.5, "B", pos = 3, font = 1, cex = 1.5)


###############################
### models: across years
dotchart(trmt.a.year$abund, group = trmt.a.year$Treatment)
hist(trmt.a.year$abund)
boxplot(trmt.a.year$abund ~ trmt.a.year$Treatment)

dotchart(trmt.a.year$rich, group = trmt.a.year$Treatment)
hist(trmt.a.year$rich)
boxplot(trmt.a.year$rich ~ trmt.a.year$Treatment)

colSums(trmt.a.year[,8:14])

dotchart(trmt.a.year$Onychiuridea, group = trmt.a.year$Treatment)
hist(trmt.a.year$Onychiuridea)
boxplot(trmt.a.year$Onychiuridea ~ trmt.a.year$Treatment)

dotchart(trmt.a.year$Isotomidea, group = trmt.a.year$Treatment)
hist(trmt.a.year$Isotomidea)
boxplot(trmt.a.year$Isotomidea ~ trmt.a.year$Treatment)

boxplot(trmt.a.year$Entomobryidae ~ trmt.a.year$Treatment)
boxplot(trmt.a.year$Hypogastruridae ~ trmt.a.year$Treatment)

#
abund.a.year <- glmer.nb(abund ~ Treatment + (1|Block), 
                      family = poisson, data = trmt.a.year)
summary(abund.a.year)
Anova(abund.a.year, type = "III")
emmeans(abund.a.year, pairwise ~ Treatment)
testDispersion(abund.a.year)
res.mod.1 <- simulateResiduals(abund.a.year)
plot(abund.a.year)
testCategorical(abund.a.year, catPred = trmt.a.year$Treatment)
testZeroInflation(abund.a.year)

#
rich.a.year <- glmer(rich ~ Treatment + (1|Block), 
                     family = poisson, data = trmt.a.year)
summary(rich.a.year)
Anova(rich.a.year, type = "III")
emmeans(rich.a.year, pairwise ~ Treatment)
testDispersion(rich.a.year)
res.mod.1 <- simulateResiduals(rich.a.year)
plot(rich.a.year)
testCategorical(rich.a.year, catPred = trmt.a.year$Treatment)
testZeroInflation(rich.a.year)

# change the names of the variables and reorder them for figure
levels(trmt.a.year$Treatment)
levels(trmt.a.year$Treatment)[levels(trmt.a.year$Treatment)=="TR_2024_1mo"] <- "2024"
levels(trmt.a.year$Treatment)[levels(trmt.a.year$Treatment)=="TR_2022"] <- "2022"
levels(trmt.a.year$Treatment)[levels(trmt.a.year$Treatment)=="TR_2023"] <- "2023"
levels(trmt.a.year$Treatment)[levels(trmt.a.year$Treatment)=="Control"] <- "Untreated"
levels(trmt.a.year$Treatment)
trmt.a.year$Treatment <- factor(trmt.a.year$Treatment,
                                levels = c("Untreated", "2022", "2023", "2024"))


boxplot(abund ~ Treatment, data = trmt.a.year, col = c("seagreen4", "seagreen2", "goldenrod2", "lightgoldenrod2"),
        ylim = c(0,8), ylab = "Abundance", xlab = "", cex.lab = 1.6, cex.axis = 1.5)
stripchart(abund ~ Treatment, data = trmt.a.year, pch = 19, cex = 2, add = TRUE, 
           vertical = TRUE, method = "jitter", jitter = 0.2)
text(0.5,7.4, "A", pos = 3, font = 1, cex = 1.5)

boxplot(rich ~ Treatment, data = trmt.a.year, col = c("seagreen4", "seagreen2", "goldenrod2", "lightgoldenrod2"),
        ylim = c(0,6), ylab = "Family Richness", xlab = "", cex.lab = 1.6, cex.axis = 1.5)
stripchart(rich ~ Treatment, data = trmt.a.year, pch = 19, cex = 2, add = TRUE, 
           vertical = TRUE, method = "jitter", jitter = 0.2)
text(0.5,5.5, "B", pos = 3, font = 1, cex = 1.5)



# maple vs hemlock
mapl <- rbind(maple, con) %>% droplevels()
str(mapl)
mapl$abund <- rowSums(mapl[,8:14])
mapl$rich <- hill_taxa(mapl[,8:14], q = 0, MARGIN = 1)
boxplot(mapl$abund ~ mapl$Treatment, pch = 19)
boxplot(mapl$rich ~ mapl$Treatment, pch = 19)

colSums(mapl[,8:14])

abund.mapl <- glmer(abund ~ Treatment + (1|Block), 
                         family = poisson, data = mapl)
summary(abund.mapl)
Anova(abund.mapl, type = "III")
emmeans(abund.mapl, pairwise ~ Treatment)
testDispersion(abund.mapl)
res.mod.1 <- simulateResiduals(abund.mapl)
plot(abund.mapl)
testCategorical(abund.mapl, catPred = mapl$Treatment)
testZeroInflation(abund.mapl)

rich.mapl <- glmer(rich ~ Treatment + (1|Block), 
                    family = poisson, data = mapl)
summary(rich.mapl)
Anova(rich.mapl, type = "III")
emmeans(rich.mapl, pairwise ~ Treatment)
testDispersion(rich.mapl)
res.mod.1 <- simulateResiduals(rich.mapl)
plot(rich.mapl)
testCategorical(rich.mapl, catPred = mapl$Treatment)
testZeroInflation(rich.mapl)

# change the names of the variables and reorder them for figure
levels(mapl$Treatment)
levels(mapl$Treatment)[levels(mapl$Treatment)=="Control"] <- "Hemlock"
levels(mapl$Treatment)[levels(mapl$Treatment)=="Sugar_Maple"] <- "Sugar Maple"
levels(mapl$Treatment)
mapl$Treatment <- factor(mapl$Treatment,
                                levels = c("Hemlock", "Sugar Maple"))


png("Figures/Collembola_hemlock_maple.png", width = 2300, height = 1000, pointsize = 30)

par(mfrow=c(1,2))
par(mar=c(5,7,4,2))

boxplot(abund ~ Treatment, data = mapl, col = c("#807DBA", "#DADAEB"),
        ylim = c(0,40), ylab = "Abundance", xlab = "", cex.lab = 1.6, cex.axis = 1.5)
stripchart(abund ~ Treatment, data = mapl, pch = 19, cex = 1.8, add = TRUE, 
           vertical = TRUE, method = "jitter", jitter = 0.2)
text(0.5,36.5, "A", pos = 3, font = 1, cex = 1.5)
text(1,9, "a", pos = 3, font = 1, cex = 1.2)
text(2,17, "b", pos = 3, font = 1, cex = 1.2)

boxplot(rich ~ Treatment, data = mapl, col = c("#807DBA", "#DADAEB"),
        ylim = c(0,6), ylab = "Family Richness", xlab = "", cex.lab = 1.6, cex.axis = 1.5)
stripchart(rich ~ Treatment, data = mapl, pch = 19, cex = 1.8, add = TRUE, 
           vertical = TRUE, method = "jitter", jitter = 0.2)
text(0.5,5.5, "B", pos = 3, font = 1, cex = 1.5)
text(1,2.2, "a", pos = 3, font = 1, cex = 1.2)
text(2,4.2, "b", pos = 3, font = 1, cex = 1.2)

dev.off()




################
# Create a panel figure

brewer.pal(9, "Purples")
brewer.pal(9, "Greens")
brewer.pal(9, "Blues")

png("Figures/Collembola_treatment_panel.png", width = 2300, height = 2000, pointsize = 30)

par(mfrow=c(2,2))
par(mar=c(5,7,4,2))

# abundance
boxplot(abund ~ Treatment, data = trmt.w.year, col = c("#C6DBEF", "#6BAED6", "#2171B5"),
        ylim = c(0,8), ylab = "Abundance", xlab = "", cex.lab = 1.6, cex.axis = 1.5)
stripchart(abund ~ Treatment, data = trmt.w.year, pch = 19, cex = 2, add = TRUE, 
           vertical = TRUE, method = "jitter", jitter = 0.2)
text(0.5,7.4, "A", pos = 3, font = 1, cex = 1.5)

boxplot(abund ~ Treatment, data = trmt.a.year, col = c("#E5F5E0", "#C7E9C0", "#A1D99B", "#238B45"),
        ylim = c(0,8), ylab = "Abundance", xlab = "", cex.lab = 1.6, cex.axis = 1.5)
stripchart(abund ~ Treatment, data = trmt.a.year, pch = 19, cex = 2, add = TRUE, 
           vertical = TRUE, method = "jitter", jitter = 0.2)
text(0.5,7.4, "C", pos = 3, font = 1, cex = 1.5)


# richness
boxplot(rich ~ Treatment, data = trmt.w.year, col = c("#C6DBEF", "#6BAED6", "#2171B5"),
        ylim = c(0,6), ylab = "Family Richness", xlab = "", cex.lab = 1.6, cex.axis = 1.5)
stripchart(rich ~ Treatment, data = trmt.w.year, pch = 19, cex = 2, add = TRUE, 
           vertical = TRUE, method = "jitter", jitter = 0.2)
text(0.5,5.5, "B", pos = 3, font = 1, cex = 1.5)

boxplot(rich ~ Treatment, data = trmt.a.year, col = c("#E5F5E0", "#C7E9C0", "#A1D99B", "#238B45"),
        ylim = c(0,6), ylab = "Family Richness", xlab = "", cex.lab = 1.6, cex.axis = 1.5)
stripchart(rich ~ Treatment, data = trmt.a.year, pch = 19, cex = 2, add = TRUE, 
           vertical = TRUE, method = "jitter", jitter = 0.2)
text(0.5,5.5, "D", pos = 3, font = 1, cex = 1.5)

dev.off()
