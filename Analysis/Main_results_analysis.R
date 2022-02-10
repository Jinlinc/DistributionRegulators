######### Distribution manuscript data analysis #########
######### Jinlin Chen #########
######### Last update: 2022 Feb #########

### Content
## 0. general preparation
## 1. distribution
## 2. thermal performance curve
## 3. regression: thermal traits - distribution 
## 4. short-term competition and equilirium states
## 5. long-term competition

## 0. preparation
## packages
library(tidyverse) # plot and data wrangling 
library(lme4) # mix effect model
library(ggrepel) # work with ggplot to add non-overlapping text labels to points
library(cowplot) # plot multiple ggplot into one graph. Compare with library(ggpubr), this package can alige the plots better
library(rstan); options(mc.cores = parallel::detectCores()); rstan_options(auto_write = TRUE) # to fit tpc by Bayesian statistics
library(bayesplot) # use its function "ppc_dens_overlay" to do fitting diagnostics 
library(ape) # to read phylogenetic tree and do relavant stuff
library(brms) # phylogenetic regression


## 1. distribution
########################## data formating #############################
distribution.dat <- read.csv("Data/pupaeSamplingCore.csv")
dist.dat.core <- distribution.dat %>%
  filter(is.na(Host) != TRUE) %>%
  filter(is.na(Site) != TRUE) %>%
  filter(Host != "Drosophila sp. 1" & Host != "sulfurigaster complex" & Host != "serrata" & Host != "immigrans") %>%
  mutate(elevation = as.numeric(substr(Site,2,4))) %>%
  mutate(t_sp = paste0(Host, "_", Transect)) %>%
  mutate(elev.s = ifelse(Site == "K070" | Site == "P070", "Low",
                         ifelse(Site == "K390" | Site == "P350", "Medium", "High"))) %>%
  mutate(hIndex = ifelse(elev.s == "Low", 0,
                         ifelse(elev.s == "Medium", 0.5, 1)))
dist.dat.core$Host <- factor(dist.dat.core$Host, levels = c("bunnanda", "pandora", "bipectinata", "pseudoananassae", "sulfurigaster", "rubida", "birchii", "palidifrons", "pseudotakahashii"))
dist.dat.core$Site <- factor(dist.dat.core$Site, levels = c("K070", "P070", "P350", "K390", "K730", "P880"))
dist.dat.core$elev.s <- factor(dist.dat.core$elev.s, levels = c("High", "Medium", "Low"))

########################## distribution visualization #############################
## The following code generates Figure 1 in the main text
## bar plot: bar length represents the proportion of records of a species among all records of this species 
hcc <- c("Low" = "#FFCC33", "Medium" = "#009900", "High" = "#0066CC")
dist.dat.core %>%
  group_by(Host, elev.s, Transect) %>%
  summarise(count = n()) -> host.summary.trans
host.summary.trans %>%
  ggplot(aes(x = Host, y = count, group = elev.s)) + 
  geom_col(aes(fill = elev.s), position = "fill") + 
  scale_fill_manual(values = hcc, name = "Elevation") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  geom_text(aes(label = count), position = position_fill(vjust = 0.5), alpha = 0.6, size = 6) +
  ylab("Proportion") + xlab("Drosophila species") +
  facet_grid(Transect ~ .) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15), axis.text.y = element_text(hjust = 1, size = 15)) + 
  theme(axis.title.y = element_text(size = 18), axis.title.x = element_text(size = 18)) + 
  theme(legend.text=element_text(size=18), legend.title = element_text(size=18), legend.position = "top") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave("distributionByTransect_direct.png", width = 6, height = 7)

########################## hIndex #############################
## hIndex: low- 0.0, mid - 0.5, high - 1. Then calculate an average ##
## these two transects differ in the exact elevation of the sampling site. So it can't be directly integrated (Paluma high is higher than Kimura high, but this is only due to sampling rather than reality)
## better to use hIndex than direct height.
dist.dat.core %>%  
  group_by(Host) %>%
  summarise(h.mean = mean(elevation), hIndex.mean = mean(hIndex)) -> h.index

########################## presence/absence - elevation regression  #############################
## construct presence/absence dataset for each species ("1" means the entry is the focal species, "0" means the entry is not the focal spcies)
## generalized linear model: y(0/1) ~ elevation * transect, family  = binomial 
## "bunnanda", "pandora", "pseudotakahashii": It was unable to include the interactive term and/or the transect term for these species who have very skewed distribution (only a few samples are present at one elevation in either or both mountain range). Fitting with those terms generated warning message: “glm.fit: fitted probabilities numerically 0 or 1 occurred”. 
## "bunnanda" still generate warning message even only including one fix effect - elevation, because there are too few samples on both transects. 
## store the coefficients and p values in glm.r dataframe
sp.list <- c("bunnanda", "pandora", "bipectinata", "pseudoananassae", "rubida", "sulfurigaster", "birchii", "palidifrons", "pseudotakahashii")
glm.r <- data.frame(sp = NA, coef.elev = NA, p.elev = NA, coef.trans = NA, p.trans = NA, coef.int = NA, p.int = NA)
for (i in 1:length(sp.list)){
  sp.name <- sp.list[i]
  temp.sp <- dist.dat.core%>%
    mutate(temp.sp = ifelse(Host == sp.name, 1, 0)) %>%
    mutate(elev.scale = elevation/1000)
  if (sp.name=="bunnanda"|sp.name=="pandora"|sp.name=="pseudotakahashii") {
    temp.glm <- glm(formula = temp.sp ~ elev.scale, family = binomial(link = "logit"), data = temp.sp)
    coef.elev <- coef(temp.glm)["elev.scale"]
    coef.trans <- "NA"
    coef.int <- "NA"
    p.elev <- coef(summary(temp.glm))[,'Pr(>|z|)']["elev.scale"]
    p.trans <- "NA"
    p.int <- "NA"
    } 
  else{
    temp.glm <- glm(formula = temp.sp ~ elev.scale * Transect, family = binomial(link = "logit"), data = temp.sp)
    coef.elev <- coef(temp.glm)["elev.scale"]
    coef.trans <- coef(temp.glm)["TransectPaluma"]
    coef.int <- coef(temp.glm)["elev.scale:TransectPaluma"]
    p.elev <- coef(summary(temp.glm))[,'Pr(>|z|)']["elev.scale"]
    p.trans <- coef(summary(temp.glm))[,'Pr(>|z|)']["TransectPaluma"]
    p.int <- coef(summary(temp.glm))[,'Pr(>|z|)']["elev.scale:TransectPaluma"]
  }
  glm.r[i,] <- c(sp.name, coef.elev, p.elev, coef.trans, p.trans, coef.int, p.int)
}


########################## hIndex - regression correlation #############################
## evaluate their correlations
dist.index <- cbind(glm.r, h.index)
## Supplemenratry Table 2
write.csv(dist.index, file = "distribution index.csv")
## Results - 1. Field distribution: values of hIndex and regression coefficient of elevation is highly correlated
rank.cor <- cor.test(dist.index$hIndex.mean, as.numeric(dist.index$coef.elev), method = "spearman") # spearman v.s. kendall
rank.cor 


## 2. thermal performance curve
########################## data formating #############################
tpc_data <- read.csv("Data/TPC_data_without na.csv")
## combine fecundity in day1-2 and day7-8, and root transform the offspring counts.
## create an matrix to hold the infomation of species and temperature for model fitting
tpc_data %>%
  filter(period != "recovery") %>%
  mutate(adultT = adultF + adultM) %>%
  group_by(sp, round, rep, temp, temp.corr) %>%
  summarize(RS = sum(adultT)) %>%
  ungroup() %>%
  select(RS, temp, temp.corr, sp) %>%
  mutate(dailyRS = RS/8, dailyRS.sqrt = sqrt(RS/8))-> tpc_data_RS
swfun.sp2 <- function(x){
  switch (x,
          "BIP" = "bipectinata", "bir" = "birchii", "BUN" = "bunnanda",
          "MEL" = "melanogaster", "PAL" = "palidifrons", "PAN" = "pandora",
          "PSA" = "pseudoananassae", "SIM" = "simulans", "SUL" = "sulfurigaster"
  )
}
tpc_data_RS$spName <- sapply(as.character(tpc_data_RS$sp), swfun.sp2)
tpc_data_RS$temp <- as.factor(tpc_data_RS$temp)
tpc_data_RS$sp <- as.factor(tpc_data_RS$sp)
tpc_data_RS$spName <- as.factor(tpc_data_RS$spName)
SP <- model.matrix(~sp, data = tpc_data_RS, 
                   contrasts.arg = list(sp=contrasts(tpc_data_RS$sp, contrasts = FALSE)))
SP <- SP[,-1]
TEMPS <- model.matrix(~temp, data = tpc_data_RS, 
                      contrasts.arg = list(temp=contrasts(tpc_data_RS$temp, contrasts = FALSE)))
TEMPS <- TEMPS[,-1]
## overview of the data
tpc_data_RS %>%
  ggplot(aes(x = temp.corr, y = dailyRS.sqrt)) + geom_point() + 
  facet_grid(.~sp)

########################## model fitting #############################
## stan assuming normal distribution with sqrt transformed input, varing sd among temperature treatment
# use guassian distribution to fit sqrt daily fecundity 
# use briere 2 function to describe the relationships between sqrt daily fecundity and temperature
# the 9 values(as for 9 species) of each parameter (a, b, RTmin, RTmax) follow a normal distribution
# the variation of fecundity is different for different temperature but the same for every species
# no bound on b also works - only one parameter Rhat is 1.01, others are 1
# the fitting takes about 10-20 mins to run
fit_tpc_sqrt_varingSD_noBound <- stan(file = "StanModels/tpc_sqrt_varingSD.stan",
                              data = list(N = nrow(tpc_data_RS),
                                          Nsp = 9,
                                          Ntemps = 7,
                                          Temp = tpc_data_RS$temp.corr,
                                          sp = SP,
                                          temps = TEMPS,
                                          y = tpc_data_RS$dailyRS.sqrt),
                              chains = 4,
                              seed = 1,
                              iter = 3000,
                              control = list(adapt_delta = 0.9, max_treedepth = 12))
# save the stanfit for future use (both are ways to save fitted stan model)
save(fit_tpc_sqrt_varingSD_noBound, file = 'StanFits/fit_tpc_sqrt_varingSD_noBound')
# if only want to examine the fitted results, you may directly read the fitted model run by Jinlin before
load("StanFits/fit_tpc_sqrt_varingSD_noBound")

## model diagnostics
# 1. print out fitted values and examine the RHat
print(fit_tpc_sqrt_varingSD_noBound, pars = c("a","b","RTmax","RTmin", "sigma_y", "mu_a", "mu_b", "mu_RTmin", "mu_RTmax", "sigma_a","sigma_b","sigma_RTmin","sigma_RTmax"))
# 2. distribution of est and obs - Supplementary Figure 2A
yrep <- rstan::extract(fit_tpc_sqrt_varingSD_noBound, pars = c("y_new"), permuted = TRUE)
yrep <- as.data.frame(yrep); yrep <- as.matrix(yrep)
yobs <- tpc_data_RS$dailyRS.sqrt
ppc_dens_overlay(yobs, yrep[1:100, ])
ggsave("tpc_diagnostics.png", width = 6, height = 4)
# 3. residual plots - Supplementary Figure 2B
yest <- rstan::extract(fit_tpc_sqrt_varingSD_noBound, pars = c("theta"), permuted = TRUE)
yest <- as.data.frame(yest); yest <- as.matrix(yest)
sigma_y.est <- rstan::extract(fit_tpc_sqrt_varingSD_noBound, pars = c("sigma_y_temps"), permuted = TRUE)
sigma_y.est <- as.data.frame(sigma_y.est); sigma_y.est <- as.matrix(sigma_y.est)
y.res <- yobs - colMeans(yest)
sigma_y.est.1 <- colMeans(sigma_y.est)
png(filename = "tpc_residuals.png", width = 1000, height = 400)
par(mfrow=c(1,2))
plot(x = colMeans(yest), y = y.res, ylab="residuals", xlab="predicted value"); abline(lm(y.res ~ colMeans(yest)), col = "red")
plot(x = colMeans(yest), y = y.res/sigma_y.est.1, ylab="standardized residuals", xlab="predicted value"); abline(lm(y.res/sigma_y.est.1 ~ colMeans(yest)), col = "red")
dev.off()
par(mfrow=c(1,1))


########################## tpc parameters and plotting #############################
## extract the model parameters
## construct predicted curves 
## plot the predicted thermal performance curves of nine species together
sp.list <- c("bipectinata", "birchii", "bunnanda", "melanogaster", "palidifrons", "pandora", "pseudoananassae", "simulans", "sulfurigaster")
pars <- rstan::extract(fit_tpc_sqrt_varingSD_noBound, pars = c("a","b","RTmax","RTmin", "RTopt"), permuted = TRUE)
a.med <- apply(pars$a, 2, FUN = median) 
b.med <- apply(pars$b, 2, FUN = median)
RTmin.med <- apply(pars$RTmin, 2, FUN = median)
RTmax.med <- apply(pars$RTmax, 2, FUN = median)
RTopt.med <- apply(pars$RTopt, 2, FUN = median)
pars.est <- data.frame(sp = rep(NA,9), a = NA, b = NA, RTmin = NA, RTmax = NA, RTopt = NA)
pars.est <- pars.est %>%
  mutate(sp = sp.list, a = a.med, b = b.med, RTmin = RTmin.med, RTmax = RTmax.med, RTopt=RTopt.med)
# ci
a.lower <- apply(pars$a, 2, quantile, probs=c(0.05))
b.lower <- apply(pars$b, 2, quantile, probs=c(0.05))
RTmax.lower <- apply(pars$RTmax, 2, quantile, probs=c(0.05))
RTmin.lower <- apply(pars$RTmin, 2, quantile, probs=c(0.05))
RTopt.lower <- apply(pars$RTopt, 2, quantile, probs=c(0.05))
a.upper <- apply(pars$a, 2, quantile, probs=c(0.95))
b.upper <- apply(pars$b, 2, quantile, probs=c(0.95))
RTmax.upper <- apply(pars$RTmax, 2, quantile, probs=c(0.95))
RTmin.upper <- apply(pars$RTmin, 2, quantile, probs=c(0.95))
RTopt.upper <- apply(pars$RTopt, 2, quantile, probs=c(0.95))
pars.est_ci <- data.frame(sp.list, a.med, a.lower, a.upper, b.med, b.lower, b.upper, 
                          RTmin.med, RTmin.lower, RTmin.upper, RTmax.med, RTmax.lower, RTmax.upper, 
                          RTopt.med, RTopt.lower, RTopt.upper)
## Table 1: TPC parameters
write.csv(pars.est_ci, "tpc parameters.csv")

## calculate predicted fecundity across temperature
temp <- seq(7.5, 33, 0.05); n.temp <- length(temp)
tpc.pred <- data.frame(spName=NA, temp=rep(temp,9), est=NA)
for (sp in 1:9){
  a = pars.est$a[sp]
  b = pars.est$b[sp]
  Tmin = pars.est$RTmin[sp]
  Tmax = pars.est$RTmax[sp]
  for (t in 1:n.temp){
    tpc.pred$spName[((sp-1)*n.temp + t)] <- as.character(sp.list[sp])
    tpc.pred$est[((sp-1)*n.temp + t)] <- ifelse(temp[t]<= Tmin, 0, 
                                            ifelse(temp[t] >= Tmax, 0, a*temp[t]*(temp[t] - Tmin)*(Tmax - temp[t])^(1/b)))
  }
}
tpc.pred <- tpc.pred %>%
  mutate(dailyRS = est^2)

## Figure 2: plotting tpc
cc.tpc <- c("melanogaster" = "grey",
         "simulans" = "#663300", "pandora" = "#FF9933", "bipectinata" = "#FFCC99", "bunnanda" = "#CC6600",
         "pseudoananassae"="#BDE51F", "sulfurigaster" = "#5FCB0B",
         "palidifrons" = "#004C99", "birchii" = "#025139")
tpc.pred$spName <- factor(tpc.pred$spName, levels = c("palidifrons", "birchii", "sulfurigaster", "pseudoananassae", "bipectinata", "pandora", "bunnanda", "simulans", "melanogaster"))
tpc.pred %>%
  ggplot(aes(x = temp, y = dailyRS, color = spName)) + geom_path(size = 1) + 
  scale_color_manual(values = cc.tpc, name = "Drosophila species") + 
  xlim(7,33) + xlab("Temperature") + ylab("Fecundity (per female per day)") + ylim(0,25) + 
  theme(axis.title.y = element_text(size = 18), axis.title.x = element_text(size = 18)) +
  theme(axis.text.y = element_text(size = 15), axis.text.x = element_text(size = 15)) +
  theme(legend.text=element_text(size=18), legend.title = element_text(size=18)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
ggsave("tpc.png", width = 10, height = 6)

## Supplementary Figure 2B
tpc_data_RS %>%
  ggplot(aes(x = temp.corr, y = dailyRS)) + geom_point() + 
  geom_path(data = tpc.pred, aes(x = temp, y = dailyRS, color = "blue")) + scale_color_manual(values = "blue") +  
  xlim(7,33) + xlab("Temperature") + ylab("Fecundity (per female per day)") + ylim(0,25) + 
  facet_grid(. ~ spName) + 
  theme(axis.title.y = element_text(size = 18), axis.title.x = element_text(size = 18)) +
  theme(axis.text.y = element_text(size = 15), axis.text.x = element_text(angle = 45, size = 15)) +
  theme(legend.position = "none") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave("tpc_summary.png", width = 10, height = 5)

## 3. regression: thermal traits - distribution 
########################## summurizing thermal traits #############################
swfun.hIndex2 <- function(x){
  switch (x,
          "bipectinata" = 0.1176, "birchii" = 0.6400, "bunnanda" = 0,
          "melanogaster" = NA, "palidifrons" = 0.7629, "pandora" = 0.0714,
          "pseudoananassae" = 0.4048, "simulans" = NA, "sulfurigaster" = 0.4286
  )
}
## use the pars.est_ci from the previous tpc analysis for RTmin, RTmax and RTopt
# need to add a t.tree column to do the phylogenetic analysis
pars.est_ci <- read.csv(file.choose())
RT.traits2 <- pars.est_ci %>% mutate(t.tree = sp.list) %>% filter(sp.list != "melanogaster" & sp.list != "simulans")
RT.traits2$hIndex <- sapply(as.character(RT.traits2$sp.list), swfun.hIndex2)

## get exact fecundity traits from tpc raw data (tpc_data)
# overview of data
tpc_data %>%
  ggplot(aes(x = temp, y = adultF+adultM)) + geom_point() + facet_grid(period~sp, switch = "both")
# ggsave("summary.jpg", width = 10)
# start extracting
tpc_data$spName <- sapply(as.character(tpc_data$sp), swfun.sp2)
## combined (day 1-2 + day 7-8) daily fecundity ##
tpc_data %>%
  mutate(adultT = adultF + adultM) %>%
  filter(period != "recovery") %>%
  group_by(spName, round, rep, temp, temp.corr) %>%
  summarize(adultT.com = sum(adultT)) %>% 
  mutate(dailyRS = adultT.com/8) -> fec.data # RS is the offspring (male+female) number for 4 days, dailyFec is only measured by female 
fec.data$spName <- factor(fec.data$spName, levels = c("bunnanda", "pandora", "bipectinata", "pseudoananassae", "sulfurigaster", "birchii", "palidifrons", "simulans", "melanogaster"))
fec.data$hIndex <- sapply(as.character(fec.data$spName), swfun.hIndex2)
fec.data$round <- as.factor(fec.data$round)
## recovered (day 9-12) daily fecundity ##
rec.data <- tpc_data %>%
  filter(period == "recovery") %>% 
  mutate(RS = adultF + adultM, dailyRS.rec = (adultF + adultM)/8) %>%  # I haven't correct for the dead flies before recovery
  select(spName, round, rep, temp, temp.corr, RS, dailyRS.rec)
rec.data$spName <- factor(rec.data$spName, levels = c("bunnanda", "pandora", "bipectinata", "pseudoananassae", "sulfurigaster", "birchii", "palidifrons", "simulans", "melanogaster"))
rec.data$hIndex <- sapply(as.character(rec.data$spName), swfun.hIndex2)
rec.data$round <- as.factor(rec.data$round)

## get knockdown tolerance traits from knockdown experiment
## hot knockdown ##
hot.data <- read.csv("Data/hot_tolerance_data.csv")
hot.data$spName <- sapply(as.character(hot.data$species), swfun.sp2)
hot.data$hIndex <- sapply(as.character(hot.data$spName), swfun.hIndex2)
hot.data$spName <- factor(hot.data$spName, levels = c("bunnanda", "pandora", "bipectinata", "pseudoananassae", "sulfurigaster", "birchii", "palidifrons", "simulans", "melanogaster"))
hot.data$kd.t <- as.numeric(hot.data$kd.t)
hot.data$hIndex <- as.numeric(hot.data$hIndex)
## cold knockdown ##
cold.data <- read.csv("Data/cold_tolerance_data.csv")
cold.data$spName <- sapply(as.character(cold.data$species), swfun.sp2)
cold.data$hIndex <- sapply(as.character(cold.data$spName), swfun.hIndex2)
cold.data$spName <- factor(cold.data$spName, levels = c("bunnanda", "pandora", "bipectinata", "pseudoananassae", "sulfurigaster", "birchii", "palidifrons", "simulans", "melanogaster"))
cold.data$kd.t <- as.numeric(cold.data$kd.t)
cold.data$hIndex <- as.numeric(cold.data$hIndex)


########################## regression and plot #############################
# loading .tre file from Filton
t <- read.tree("Data/All taxa_Bayesian.con.tre")
# prune the tree to just my species
sp.list.fullname <- c("D.bunnanda", "D.pandora", "D.bipectinata", "D.pseudoananassae", "D.rubida", "D.sulfurigaster", "D.birchii", "D.palidifrons", "D.pseudotakahashii", "D.melanogaster", "D.simulans")
t$tip.label[match(sp.list.fullname, t$tip.label)] # first check the species overlap
# palidifrons and pandora are unavailable
tip.avail <- c("D.bunnanda", "D.bipectinata", "D.pseudoananassae", "D.rubida", "D.sulfurigaster", "D.birchii", "D.pseudotakahashii", "D.melanogaster", "D.simulans")
t.avail <- keep.tip(t, t$tip.label[match(tip.avail, t$tip.label)])
write.tree(t.avail)# this is how the phylogenetic structure is written in text format. Tips and nodes are represented by thier name, the branth distance from this tip/node to the nearest node is indicated by the value after ":"
# Therefore, I can copy the above style to incorporate D.palidifrons and D.pandora into the tree
# It's known that D.palidiforns is sister species with D. sulfurigaster, and D.pandora with D.bipectinata.
# For simplicity, assume that the newly added tip is halfway along the edge. Then the tree can be written as below:
t.all <- read.tree(text = "(((((D.melanogaster:0.023028,D.simulans:0.016768)Node1:0.069284,D.pseudotakahashii:0.070281)Node2:0.061871,
                   (D.bunnanda:0.05867,D.birchii:0.061723)Node3:0.079749)Node4:0.018618,
                   (D.pseudoananassae:0.022019,(D.bipectinata:0.0091415,D.pandora:0.0091415)Node5:0.0091415)Node6:0.130325)Node7:0.105535,
                   (D.rubida:0.104915,(D.sulfurigaster:0.055463,D.palidifrons:0.055463)Node8:0.055463)Node9:0.071313)Node10;")
## Supplementary Figure 3A - phylogeny of all species
plot(t.all)
# the species (tips) that are examnined in the tpc experiment are below:
tip.tpc <- c("D.bunnanda", "D.pandora", "D.bipectinata", "D.pseudoananassae", "D.sulfurigaster", "D.birchii", "D.palidifrons")
# prune my.tree accordingly:
t.tree <- keep.tip(t.all, t.all$tip.label[match(tip.tpc, t.all$tip.label)])
# ultrametric tree - needed to calculate the inverse matrix
is.ultrametric(t.tree) # not!
# convert to ultrametric tree - 1
# t.tree <- phytools::nnls.t.tree(cophenetic(t.tree),t.tree, rooted=TRUE, trace=0)
# because I can't install phytools in my current R version (it is already up to date), I found the code to manually compute the ultrametric tree by extend method.
force.ultrametric.extend<-function(tree){
  h<-diag(vcv(tree))
  d<-max(h)-h
  ii<-sapply(1:Ntip(tree),function(x,y) which(y==x),
             y=tree$edge[,2])
  tree$edge.length[ii]<-tree$edge.length[ii]+d
  return(tree)
}
# convert to ultrametric tree - 2
t.tree <- force.ultrametric.extend(t.tree)
## Supplementary Figure 3B - uphylogeny of the seven species used in thermal traits regression
plot(t.tree)
# to match the tip name here with the species name in other datasets
new_tip <- c("bunnanda", "pandora", "bipectinata", "pseudoananassae", "sulfurigaster", "birchii", "palidifrons")
t.tree$tip.label <- new_tip[match(t.tree$tip.label,tip.tpc)]
## create the between-level covariance matrix based on the phylo dataset. VERY USEFUL!!
inv.phylo <- MCMCglmm::inverseA(t.tree, nodes = "TIPS", scale = TRUE) # inverseA(): Inverse Relatedness Matrix and Phylogenetic Covariance Matrix
A <- solve(inv.phylo$Ainv)
rownames(A) <- rownames(inv.phylo$Ainv)


# RTopt #
fit.RTopt <- brm(
  RTopt.med ~ hIndex + (1|gr(t.tree, cov = A)),
  data = RT.traits2, 
  family = gaussian(),
  data2 = list(A = A))
# Result-2.Thermal perfomance curves
# "The temperature for optimal reproductive performance, RTopt, did not correlate with their distribution patterns (Coefficient = 0.10, 95% credible interval -2.81 – 3.05). "
summary(fit.RTopt)
plot(fit.RTopt, pars = c("hIndex"))

# RTmax #
fit.RTmax <- brm(
  RTmax.med ~ hIndex + (1|gr(t.tree, cov = A)),
  data = RT.traits2, 
  family = gaussian(),
  data2 = list(A = A))
# Result-4.Heat tolerance
# "species whose distribution were biased towards lowland consistently had higher RTmax (Figure 3e. Coefficient = -3.01, 95% credible interval = - 5.08 – - 0.72). "
summary(fit.RTmax)
plot(fit.RTmax, pars = c("hIndex"))

# RTmin #
fit.RTmin <- brm(
  RTmin.med ~ hIndex + (1|gr(t.tree, cov = A)),
  data = RT.traits2, 
  family = gaussian(),
  data2 = list(A = A))
# Result-3.Cold tolerance
# "Values of RTmin were not correlated with the species distribution patterns (Figure 3a. Coefficient = -0.35, 95% credible interval = -4.69 – 3.49). "
summary(fit.RTmin)
plot(fit.RTmin, pars = c("hIndex"))

# fecundity at 29C #
fec.data %>%
  filter(temp == 29) %>%
  filter(spName != "melanogaster" & spName != "simulans") %>%
  ungroup() %>%
  select(adultT.com, dailyRS, hIndex, round, spName) %>%
  mutate(t.tree = spName)-> test.data
# no phylogenetic correction
test.data %>%
  glmer.nb(formula = adultT.com ~ hIndex + round + (1 | spName)) %>%
  summary()
# with phylogenetic correction by brms
# lowland species are significantly producing more in 29C.
fit.29C <- brm(
  adultT.com ~ hIndex + round + (1|spName) + (1|gr(t.tree, cov = A)),
  data = test.data, 
  family = negbinomial(link = "log"),
  data2 = list(A = A))
# Result-4.Heat tolerance
# "Reproductive performance at 29°C also decreased with hIndex (Figure 3f. Coefficient = -5.70, 95% credible interval -8.79 – -2.54). "
summary(fit.29C) 
plot(fit.29C, pars = c("hIndex"))

# fecundity at 17C #
fec.data %>%
  filter(temp == 17) %>%
  filter(spName != "melanogaster" & spName != "simulans") %>%
  ungroup() %>%
  select(adultT.com, dailyRS, hIndex, round, spName) %>%
  mutate(t.tree = spName)-> test.data
# no phylogenetic correction
test.data %>%
  glmer.nb(formula = adultT.com ~ hIndex + round + (1 | spName)) %>%
  summary()
# with phylogenetic correction by brms
# no relationship between distribtion and fecundity at 17C
fit.17C <- brm(
  adultT.com ~ hIndex + round + (1|spName) + (1|gr(t.tree, cov = A)),
  data = test.data, 
  family = negbinomial(link = "log"),
  data2 = list(A = A))
# Result-3.Cold tolerance
# "Similarly, upland-biased species did not show higher fecundity at the low temperature, 17°C (Figure 3b. Coefficient = -0.21, 95% credible interval -3.96 – 3.42). "
summary(fit.17C)  
plot(fit.17C, pars = c("hIndex"))

# fecundity recovered from 29C # 
rec.data %>%
  filter(temp == 29) %>%
  filter(spName != "melanogaster" & spName != "simulans") %>%
  ungroup() %>%
  select(RS, dailyRS.rec, hIndex, round, spName) %>%
  mutate(t.tree = spName, rec = ifelse(RS > 0, 1, 0))-> test.data
# with phylogenetic correction by brms - zero_inflated_negative binomial(RS) - no convergence
# with phylogenetic correction by brms - zero_inflated_binomial(recovery/non-recovery) - poor convergence but acceptable
# upland species are less likely to recover their reproduction, but not significant
fit.rec29C <- brm(
  rec ~ hIndex + round + (1|spName) + (1|gr(t.tree, cov = A)),
  data = test.data, 
  family = zero_inflated_binomial(link = "logit", link_zi = "logit"),
  data2 = list(A = A),
  control = list(adapt_delta = 0.90))
# Result-4.Heat tolerance - quantitative results not included
summary(fit.rec29C) 

# fecundity recovered from 14C # 
rec.data %>%
  filter(temp == 14) %>%
  filter(spName != "melanogaster" & spName != "simulans") %>%
  ungroup() %>%
  select(RS, dailyRS.rec, hIndex, round, spName) %>%
  mutate(t.tree = spName)-> test.data
# with phylogenetic correction by brms
# upland species have higher recovered fecundity but not significant
fit.rec14C <- brm(
  RS ~ hIndex + round + (1|spName) + (1|gr(t.tree, cov = A)),
  data = test.data, 
  family = negbinomial(link = "log"),
  data2 = list(A = A),
  control = list(adapt_delta = 0.90))
# Results-3.cold tolerance
# "This recovered fecundity showed a minor but non-significant increase among upland species (Figure 3c. Coefficient = 0.32, 95% credible interval -0.62 – 1.09). "
summary(fit.rec14C) 

# knockdown in 40C #
# Knockdown time at lethal high temperature (40°C) was lower among upland species (Figure 3h. Male: coefficient = -7.83 (-14.67 – -1.20); female: coefficient = -3.12 (-10.43 – 3.77)), 
hot.data2F <- hot.data %>%
  filter(spName != "melanogaster" & spName != "simulans") %>%
  filter(gender == "F") %>%
  mutate(t.tree = spName)
hot.data2M <- hot.data %>%
  filter(spName != "melanogaster" & spName != "simulans") %>%
  filter(gender == "M") %>%
  mutate(t.tree = spName)
# no phylogenetic correction
hot.data2F %>%
  lmer(formula = kd.t ~ hIndex + round + position + (1|spName)) %>% summary()
hot.data2M %>%
  lmer(formula = kd.t ~ hIndex + round + position + (1|spName)) %>% summary()
# with phylogenetic correction
# upland species has lower heat tolerance for both male and female. Female is not significant but male is significant
fit.kdt.F <- brm(
  kd.t ~ hIndex + round + position + (1|spName) + (1|gr(t.tree, cov = A)),
  data = hot.data2F, 
  family = gaussian(),
  data2 = list(A = A))
summary(fit.kdt.F)
fit.kdt.M <- brm(
  kd.t ~ hIndex + round + position + (1|spName) + (1|gr(t.tree, cov = A)),
  data = hot.data2M, 
  family = gaussian(),
  data2 = list(A = A),
  control = list(adapt_delta = 0.90))
summary(fit.kdt.M)

# knockdown and recover by 5C #
# "It took longer for upland species to regain mobility after the chill coma (Figure 3d. Male: coefficient = 14.14 (-8.43 – 36.8); female: coefficient = 9.44 (-1.77 – 22.26)), "
cold.data2F <- cold.data %>%
  filter(spName != "melanogaster" & spName != "simulans") %>%
  filter(gender == "F") %>%
  mutate(t.tree = spName)
cold.data2M <- cold.data %>%
  filter(spName != "melanogaster" & spName != "simulans") %>%
  filter(gender == "M") %>%
  mutate(t.tree = spName)
# no phylogenetic correction
cold.data2F %>%
  lmer(formula = kd.t ~ hIndex + round + position + (1|spName)) %>% summary()
cold.data2M %>%
  lmer(formula = kd.t ~ hIndex + round + position + (1|spName)) %>% summary()
cold.data2F %>%
  lmer(formula = rc.t ~ hIndex + round + position + (1|spName)) %>% summary()
cold.data2M %>%
  lmer(formula = rc.t ~ hIndex + round + position + (1|spName)) %>% summary()
# no relationship for knockdown time, so only analyze recovery time with the phylogenetic correction
# upland species has lower heat tolerance for both male and female. Female is not significant but male is significant
fit.rct.F <- brm(
  rc.t ~ hIndex + round + position + (1|spName) + (1|gr(t.tree, cov = A)),
  data = cold.data2F, 
  family = gaussian(),
  data2 = list(A = A))
summary(fit.rct.F)
fit.rct.M <- brm(
  rc.t ~ hIndex + round + position + (1|spName) + (1|gr(t.tree, cov = A)),
  data = cold.data2M, 
  family = gaussian(),
  data2 = list(A = A),
  control = list(adapt_delta = 0.90))
summary(fit.rct.M)
# "When exposed to acute sublethal low temperature (5°C), all seven tropical Drosophila species showed similarly weak resistance compared to D. simulans and D. melanogaster "
cold.kd.aov <- aov(kd.t ~ species + round + position + gender, data = cold.data) 
summary(cold.kd.aov)
aov.r <- TukeyHSD(cold.kd.aov)
write.csv(as.data.frame(aov.r$species), file = "coldKnockdownAOV.csv")


# print out all the regression results
sink("regression_results.txt")
print("========= RTopt ==========")
print(summary(fit.RTopt))
print("========= RTmax ==========")
print(summary(fit.RTmax))
print("========= RTmin ==========")
print(summary(fit.RTmin))
print("========= 29C ==========")
print(summary(fit.29C))
print("========= 17C ==========")
print(summary(fit.17C))
print("========= rec29C ==========")
print(summary(fit.rec29C))
print("========= rec14C ==========")
print(summary(fit.rec14C))
print("========= heat kd F ==========")
print(summary(fit.kdt.F))
print("========= heat kd M ==========")
print(summary(fit.kdt.M))
print("========= cold rc F ==========")
print(summary(fit.rct.F))
print("========= cold rc M ==========")
print(summary(fit.rct.M))
sink()

## Figure 3 - Thermal traits ~ hIndex plot (4*2)
# define the levels for the factor(spName) for all the dataset needed
pars.est_ci$sp.list <- factor(pars.est_ci$sp.list,levels = c("bunnanda", "pandora", "bipectinata", "pseudoananassae", "sulfurigaster", "birchii", "palidifrons", "simulans", "melanogaster"))
fec.data$spName <- factor(fec.data$spName, levels = c("bunnanda", "pandora", "bipectinata", "pseudoananassae", "sulfurigaster", "birchii", "palidifrons", "simulans", "melanogaster"))
rec.data$spName <- factor(rec.data$spName, levels = c("bunnanda", "pandora", "bipectinata", "pseudoananassae", "sulfurigaster", "birchii", "palidifrons", "simulans", "melanogaster"))
hot.data$spName <- factor(hot.data$spName, levels = c("bunnanda", "pandora", "bipectinata", "pseudoananassae", "sulfurigaster", "birchii", "palidifrons", "simulans", "melanogaster"))
cold.data$spName <- factor(cold.data$spName, levels = c("bunnanda", "pandora", "bipectinata", "pseudoananassae", "sulfurigaster", "birchii", "palidifrons", "simulans", "melanogaster"))
# individual plots
p1.1 <- pars.est_ci %>%
  filter(sp.list != "melanogaster" & sp.list != "simulans") %>%
  ggplot(aes(x = sp.list, y = RTmin.med)) + geom_point(size = 2) + 
  geom_errorbar(aes(ymin=RTmin.lower, ymax=RTmin.upper), width=.2) + 
  ylim(12, 18) + ylab("RTmin(°C)") + 
  theme(axis.title.y = element_text(size = 10), axis.text.y = element_text(size = 10)) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p1.2 <- pars.est_ci %>%
  filter(sp.list != "melanogaster" & sp.list != "simulans") %>%
  ggplot(aes(x = sp.list, y = RTmax.med)) + geom_point(size = 2) + 
  geom_errorbar(aes(ymin=RTmax.lower, ymax=RTmax.upper), width=.2) + 
  ylim(28, 32.5) + ylab("RTmax(°C)") + 
  theme(axis.title.y = element_text(size = 10), axis.text.y = element_text(size = 10)) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p2.1 <- fec.data %>%
  filter(spName != "melanogaster" & spName != "simulans") %>%
  filter(temp == 17) %>%
  ggplot(aes(x = spName, y = dailyRS)) + geom_boxplot() +
  ylab("Fecundity in 17°C") + theme(legend.position = "none") +
  theme(axis.title.y = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p2.2 <- fec.data %>%
  filter(spName != "melanogaster" & spName != "simulans") %>%
  filter(temp == 29) %>%
  ggplot(aes(x = spName, y = dailyRS)) + geom_boxplot() +
   ylab("Fecundity in 29°C") + theme(legend.position = "none") +
  theme(axis.title.y = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
p3.1 <- rec.data %>%
  filter(spName != "melanogaster" & spName != "simulans") %>%
  filter(temp == 14) %>%
  ggplot(aes(x = spName, y = dailyRS.rec)) + geom_boxplot() +
  ylab("Resumed fecundity from 14°C") + theme(legend.position = "none") +
  theme(axis.title.y = element_text(size = 9), axis.text.y = element_text(size = 10)) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
p3.2 <- rec.data %>%
  filter(spName != "melanogaster" & spName != "simulans") %>%
  filter(temp == 29) %>%
  ggplot(aes(x = spName, y = dailyRS.rec)) + geom_boxplot() +
  ylab("Resumed fecundity from 29°C") + theme(legend.position = "none") +
  theme(axis.title.y = element_text(size = 9), axis.text.y = element_text(size = 10)) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p4.1 <- cold.data %>%
  filter(spName != "melanogaster" & spName != "simulans") %>%
  ggplot(aes(x = spName, y = rc.t)) + geom_boxplot(aes(fill = gender)) + 
  scale_fill_manual(values=c("#FFFFFF", "#999999")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15)) +
  ylab("Coldshock recovery (mins)") + ylim(0,61) + xlab("Species") +
  theme(legend.position = "none") +
  theme(axis.title.y = element_text(size = 9), axis.text.y = element_text(size = 10), 
        axis.title.x = element_text(size = 10), axis.text.x = element_text(size = 10)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p4.2 <- hot.data %>%
  filter(spName != "melanogaster" & spName != "simulans") %>%
  ggplot(aes(x = spName, y = kd.t)) + geom_boxplot(aes(fill = gender)) + 
  scale_fill_manual(values=c("#FFFFFF", "#999999")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15)) +
  ylab("Heatshock knockdown (mins)") + ylim(0,35) + xlab("Species") +
  theme(legend.position = "none") +
  theme(axis.title.y = element_text(size = 9), axis.text.y = element_text(size = 10), 
        axis.title.x = element_text(size = 10), axis.text.x = element_text(size = 10)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
# put plots together (figure 3)
plot_grid(p1.1, p1.2, p2.1, p2.2, p3.1, p3.2,
          p4.1 + theme(axis.title.x=element_blank(), axis.text.x=element_blank()), 
          p4.2 + theme(axis.title.x=element_blank(), axis.text.x=element_blank()),
          ncol=2, nrow = 4, align="v")
ggsave("all_traits.png", width = 6, height = 7)
plot_grid(p4.1, p4.2, 
          ncol = 2, nrow = 1)
ggsave("physiology_traits.png", width = 6, height = 3)

## Supplementary Figure 6.A - RTopt ~ hIndex
pars.est_ci %>%
  filter(sp.list != "melanogaster" & sp.list != "simulans") %>%
  ggplot(aes(x = sp.list, y = RTopt.med)) + geom_point(size = 2) + 
  geom_errorbar(aes(ymin=RTopt.lower, ymax=RTopt.upper), width=.2) + 
  ylim(22, 27) + ylab("RTopt(°C)") + xlab("Drosophila species") + 
  theme(axis.title.y = element_text(size = 10), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10), axis.text.x = element_text(size = 10, angle = 45, hjust = 1)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave("RTopt.png", width = 4, height = 3)

## trade-off between RTmax and RTmin
pars <- rstan::extract(fit_tpc_sqrt_varingSD_noBound, pars = c("a","b","RTmax","RTmin", "RTopt"), permuted = TRUE)
# extracting pars
swfun.sp <- function(x){
  switch (x,
          "A" = "bipectinata", "B" = "birchii", "C" = "bunnanda",
          "D" = "melanogaster", "E" = "palidifrons", "F" = "pandora",
          "G" = "pseudoananassae", "H" = "simulans", "I" = "sulfurigaster"
  )
}
RTmaxs <- as.data.frame(as.table(pars[["RTmax"]]))
RTmins <- as.data.frame(as.table(pars[["RTmin"]]))
RTmaxs$spName <- sapply(as.character(RTmaxs$Var2), swfun.sp)
RTmaxs %>%
  mutate(RTmax = Freq, RTmin = RTmins$Freq) %>%
  dplyr::select(spName, RTmax, RTmin) -> RT.traits
RT.traits$spName <- factor(RT.traits$spName, levels = c("palidifrons", "birchii", "sulfurigaster", "pseudoananassae", "bipectinata", "pandora", "bunnanda", "simulans", "melanogaster"))
## Result-2.thermal performance curves
## "There was no general trade-off between cold tolerance (RTmin) versus heat tolerance (RTmax) that correspond to their distribution types (Spearman’s rank correlation rho = -0.6, p value = 0.10). "
cor.test(pars.est_ci$RTmax.med, pars.est_ci$RTmin.med, method = "spearman")
## Supplementary Figure 6.B 
RT.traits %>% 
  ggplot(aes(x = RTmin, y = RTmax, color = spName)) + 
  geom_point(size = 0.6) + 
  scale_color_manual(values = cc.tpc, name = "Drosophila species") + 
  theme(axis.title.y = element_text(size = 18), axis.title.x = element_text(size = 18)) +
  theme(axis.text.y = element_text(size = 15), axis.text.x = element_text(size = 15)) +
  theme(legend.text=element_text(size=12), legend.title = element_text(size=15)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
ggsave("tradeOff.png", width = 7, height = 4)
  


## 4. short-term competition and equilirium states
########################## preparation #############################
# species list - the order of the species must not be changed in this list!!
species <- c('BIP', 'PAL', 'PAN', 'PST', 'SUL') 
# pair indicated by their species index
pairs <- rbind(c(1,2), c(1,3), c(1,4), c(1,5), c(2,5), c(2,3))
# array of Aij (inter- or intra-specific competition terms) for all 6 pairs
locateA17 <- cbind(c(1,2,3,4,5), c(6,7,8,NA,9), c(10,11,12,NA,NA), c(13,NA,NA,14,NA), c(15,16,NA,NA,17))
## colorcode
cols = c('PST'  = 'skyblue', 
         'PAL' = 'darkblue',
         'PAN' = 'red',
         'BIP'= 'gold',  
         'SUL' = 'purple', 
         "NA" = "black")


########################## fit upland data #############################
## read in data
dat <- read.csv('Data/sixPairs_cold.csv')

## prepare the data for model fitting
dat <- dat %>%
  mutate(r.ind = NA, Aii.ind = NA, Aij.ind = NA, blockN = NA) 
# "ind"" indicates which growth rate/alpha should be used for this entry
# the following code assign the "ind"s to each entry
for (k in 1:length(dat$tubeID)){
  dat$r.ind[k] = dat$i.ind[k]
  dat$Aii.ind[k] = locateA17[dat$i.ind[k], dat$i.ind[k]]
  dat$Aij.ind[k] = ifelse(dat$type[k] == "intra", 0, locateA17[dat$j.ind[k], dat$i.ind[k]])
  temp <- tail(strsplit(dat$block[k],"")[[1]], n = 1)
  dat$blockN[k] = ifelse(dat$block[k] == "SEP_ coldA", 4,  # for PAN and PAL, those blocks done during SEP are block 4 and 5
                         ifelse(dat$block[k] == "SEP_ coldB", 5,
                                ifelse(temp == "A", 1, ifelse(temp == "B", 2, 3))))
}

## write model.stan file
# the function creat the stan file. This stan file works on one treatment at one time, it analyzes the 6 pairs together. There are 3 blocks for five pairs done in Dec and 2 blocks for the PAN_PAL pair done in Sep.
source('NegBFunction_1treatment_6Pairs_5block.R')
Write_Stan_NegB_1layer_6pairs_5blocks(PreparedData1 = dat,
                                      '6pairs_cold',
                                      AlphaLowerBound = TRUE, AlphaSD= 1) # the .stan doc is stored in StanModels file

## fit the competition model by calling stan
Num_Draws <- 2000
fit_cold <- stan(file = 'StanModels/BuiltModel_NegB_6pairs_cold.stan',
                 data = list(N = nrow(dat), 
                             y = c(dat$obs_count)) , 
                 chains = 4, seed = 1, iter = Num_Draws)
save(fit_cold, file = 'StanFits/fit_6pairs_cold')

## Diagnostics of fitting cold treatment data
# check values and convergence: all but one Rhat = 1.00, satisfying convergence
print(fit_cold, pars = c('mur', 'A', 'phi', "r1", "r2","r3","r4","r5","sigmar")) 
## Supplementary figure 5A - distribution of est and obs
yrep <- rstan::extract(fit_cold, pars = c("y_sim"), permuted = TRUE)
yrep <- as.data.frame(yrep); yrep <- as.matrix(yrep)
yobs <- dat$obs_count  
ppc_dens_overlay(yobs, yrep[1:100, ]) + ggtitle("cold treatment")
ggsave("competition_cold_diagnostic.png", width = 6, height = 4)
## standardized residual plots - not included
# theoretically it is feasible to calculate each expected y and the residuals. 
# However, the way to calculate expected y in this script is not one function or matrix multiplication, but each expected y is expressed explicitely
# So it is very tedious to get all the expected y here. Therefore, bypass this diagnostics


########################## fit lowland data #############################
## read in data
dat <- read.csv('Data/sixPairs_hot.csv')

## prepare the data for model fitting
dat <- dat %>%
  mutate(r.ind = NA, Aii.ind = NA, Aij.ind = NA, blockN = NA) 
# "ind"" indicates which growth rate/alpha should be used for this entry
# the following code assign the "ind"s to each entry
for (k in 1:length(dat$tubeID)){
  dat$r.ind[k] = dat$i.ind[k]
  dat$Aii.ind[k] = locateA17[dat$i.ind[k], dat$i.ind[k]]
  dat$Aij.ind[k] = ifelse(dat$type[k] == "intra", 0, locateA17[dat$j.ind[k], dat$i.ind[k]])
  temp <- tail(strsplit(dat$block[k],"")[[1]], n = 1)
  dat$blockN[k] = ifelse(dat$block[k] == "SEP_ hotA", 4,  # for PAN and PAL, those blocks were done during SEP are block 4 and 5
                         ifelse(dat$block[k] == "SEP_ hot", 5,
                                ifelse(temp == "A", 1, ifelse(temp == "B", 2, 3))))
}

## write model.stan file
# the function creat the stan file. This stan file works on one treatment at one time, it analyzes the 6 pairs together. There are 3 blocks for five pairs done in Dec and 2 blocks for the PAN_PAL pair done in Sep.
source('NegBFunction_1treatment_6Pairs_5block.R')
Write_Stan_NegB_1layer_6pairs_5blocks(PreparedData1 = dat,
                                      '6pairs_hot',
                                      AlphaLowerBound = TRUE, AlphaSD= 1)

## fit the competition model by calling stan
Num_Draws <- 2000
fit_hot <- stan(file = 'StanModels/BuiltModel_NegB_6pairs_hot.stan',
                data = list(N = nrow(dat), 
                            y = c(dat$obs_count)) , 
                chains = 4, seed = 1, iter = Num_Draws, 
                control = list(adapt_delta = 0.85))
save(fit_hot, file = 'StanFits/fit_6pairs_hot')

## Diagnostics of fitting cold treatment data
# check values and convergence: the maximum Rhat is 1.09, and most of Rhat are 1.00-1.03. 
# the model didn't converge well. Changing adapt_delta to 0.90 make the convergence worse
# the hot treatment data has some species having very small offspring count (zero or close to zero), which may be the reason for the algorithm to converge.
# I have no other way to improve convergence, so I will just use this fitted result. 
print(fit_hot, pars = c('mur', 'A', 'phi', "r1", "r2","r3","r4","r5","sigmar")) 
## Supplementary figure 5B - distribution of est and obs
yrep <- rstan::extract(fit_hot, pars = c("y_sim"), permuted = TRUE)
yrep <- as.data.frame(yrep); yrep <- as.matrix(yrep)
yobs <- dat$obs_count  
ppc_dens_overlay(yobs, yrep[1:100, ]) + ggtitle("hot treatment")
ggsave("competition_hot_diagnostic.png", width = 6, height = 4)
## standardized residual plots - not included for the same reason


########################## visualizing inter-specific competition  #############################
# optional code to conduct the follwoing analysis without fitting the model again - use the model I fitted before
load("StanFits/fit_6pairs_hot")
load("StanFits/fit_6pairs_cold")

## functions to predict the effect of spB to spA
prediction <- function(spA.ind, spB.ind, stanfit){
  # construct fitted values
  other.den <- seq(0, 16, 0.1)
  # get model parameters from stanfit
  draw_dat <- rstan::extract(stanfit)
  mur <- draw_dat$mur
  A <- draw_dat$A
  # ci of the fittled line
  drawN <- 200 # using the first 200 draws (should be enough)
  mur.spA <- mur[1:drawN, spA.ind]
  mur.spB <- mur[1:drawN, spB.ind]
  A.AA <- A[1:drawN, locateA17[spA.ind, spA.ind]]
  A.BA <- A[1:drawN, locateA17[spB.ind, spA.ind]]
  A.AB <- A[1:drawN, locateA17[spA.ind, spB.ind]]
  A.BB <- A[1:drawN, locateA17[spB.ind, spB.ind]]
  pred.spA <- matrix(data = NA, nrow = drawN, ncol = length(other.den))
  pred.spB <- matrix(data = NA, nrow = drawN, ncol = length(other.den))
  for (i in 1:drawN){
    pred.spA[i,] = (mur.spA[i]*8)/((1+A.AA[i]*8 + A.BA[i]*other.den)^1)
    pred.spB[i,] = (mur.spB[i]*8)/((1+A.BB[i]*8 + A.AB[i]*other.den)^1)
  }
  pred.spA.lower  <- apply(pred.spA,  2, quantile, probs=c(0.05))
  pred.spA.upper  <- apply(pred.spA,  2, quantile, probs=c(0.95))
  pred.spB.lower  <- apply(pred.spB,  2, quantile, probs=c(0.05))
  pred.spB.upper  <- apply(pred.spB,  2, quantile, probs=c(0.95))
  # est of the fitted line
  mur.spA <- median(mur[, spA.ind])
  mur.spB <- median(mur[, spB.ind])
  b.spA <- 1
  b.spB <- 1
  A.AA <- median(A[, locateA17[spA.ind, spA.ind]])
  A.BA <- median(A[, locateA17[spB.ind, spA.ind]])
  A.AB <- median(A[, locateA17[spA.ind, spB.ind]])
  A.BB <- median(A[, locateA17[spB.ind, spB.ind]])
  fitted.spA <- (mur.spA*8)/((1+A.AA*8 + A.BA*other.den)^b.spA)
  fitted.spB <- (mur.spB*8)/((1+A.BB*8 + A.AB*other.den)^b.spB)
  pred.A <- data.frame(FocalSp = species[spA.ind], FocalSp.ind = spA.ind,
                       OtherSp = species[spB.ind], OtherSp.ind = spB.ind, OtherSp_Den = other.den, 
                       median.pred = fitted.spA, upper.pred = pred.spA.upper, lower.pred = pred.spA.lower)
  pred.B <- data.frame(FocalSp = species[spB.ind], FocalSp.ind = spB.ind,
                       OtherSp = species[spA.ind], OtherSp.ind = spA.ind, OtherSp_Den = other.den, 
                       median.pred = fitted.spB, upper.pred = pred.spB.upper, lower.pred = pred.spB.lower)
  pred <- rbind(pred.A, pred.B)
  return(pred)
}

# dataframe "allPredxxx" record the predicted effect of competition
# cold treatment
allPred.cold <- data.frame(FocalSp = NA, FocalSp.ind = NA, OtherSp = NA, OtherSp.ind = NA, 
                           OtherSp_Den = NA, 
                           median.pred = NA, upper.pred = NA, lower.pred = NA)
for (i in 1:length(pairs[,1])){
  spA.ind <- pairs[i,1]
  spB.ind <- pairs[i,2]
  temp <- prediction(spA.ind, spB.ind, fit_cold)
  allPred.cold <- rbind(allPred.cold, temp)
}
allPred.cold <- allPred.cold[-1,]

# hot treatment
allPred.hot <- data.frame(FocalSp = NA, FocalSp.ind = NA, OtherSp = NA, OtherSp.ind = NA, 
                          OtherSp_Den = NA, 
                          median.pred = NA, upper.pred = NA, lower.pred = NA)
for (i in 1:length(pairs[,1])){
  spA.ind <- pairs[i,1]
  spB.ind <- pairs[i,2]
  temp <- prediction(spA.ind, spB.ind, fit_hot)
  allPred.hot <- rbind(allPred.hot, temp)
}
allPred.hot <- allPred.hot[-1,]

# combine the above two dataset together
allPred.hot <- allPred.hot %>%
  mutate(site = "24/28.5°C (lowland)", group = paste(FocalSp, "_", OtherSp))
allPred.cold <- allPred.cold %>%
  mutate(site = "21/23°C (upland)", group = paste(FocalSp, "_", OtherSp))
allall <- rbind(allPred.cold, allPred.hot)
allall$FocalSp <- factor(allall$FocalSp, levels = c("PST", "PAL", "SUL", "BIP", "PAN"))
alltext <- allall %>% filter(OtherSp_Den == 16) %>% mutate(OtherSp_Den = 16)

# Figure 4 - effect of shorterm interspecific competition 
allall %>% 
  ggplot(aes(x = OtherSp_Den, group = group)) +
  geom_line(aes(x = OtherSp_Den, y = median.pred/4, col = FocalSp)) + 
  geom_ribbon(aes(ymin=lower.pred/4, ymax=upper.pred/4, fill = FocalSp), alpha = 0.15) +
  geom_text_repel(aes(x = OtherSp_Den, y = median.pred/4, label = group), 
                  data = alltext, size = 3, 
                  min.segment.length = 0, 
                  hjust = 0, nudge_x = 1, direction = "y") + 
  facet_grid(. ~ site) + theme(strip.text.x = element_text(size = 15)) + 
  scale_color_manual('Focal species', values=cols) + 
  scale_fill_manual('Focal species', values=cols) + 
  xlab('Number of individuals of competing species') + xlim(0, 22) + 
  ylab('Fitted fecundity of focal species') + ylim(0, 22) + 
  theme(legend.position = "top") +
  theme(axis.title.y = element_text(size = 15), axis.title.x = element_text(size = 15)) +
  theme(axis.text.y = element_text(size = 15), axis.text.x = element_text(size = 15)) +
  theme(legend.text=element_text(size=15), legend.title = element_text(size=15)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave("figure4-shortterm competition_ci.png", width = 8, height = 7)


########################## model parameter output  #############################
## hot treatment
draw_dat <- rstan::extract(fit_hot, pars= c("mur", "A"))
i <- c(1,1,1,1,2,2,2,3,3,4,5,5) # all possibility of focal species
j <- c(2,3,4,5,1,3,5,1,2,1,1,2) # all possibility of competing species
par.table <- data.frame(FocalSp = NA, R0.est = NA, R0.ci90 = NA, a.est = NA, a.ci90 = NA, Competitor = NA, beta.est = NA, beta.ci90 = NA)
for (k in 1:length(i)){
  FocalSp  <- species[i[k]]
  Competitor <- species[j[k]]
  R0.est <- round(median(draw_dat$mur[,i[k]]), 2)
  R0.ci90 <- paste0(round(quantile(draw_dat$mur[,i[k]], probs = c(0.05)), 2), "-", 
                    round(quantile(draw_dat$mur[,i[k]], probs = c(0.95)), 2))
  a.est <- round(median(draw_dat$A[, locateA17[i[k], i[k]]]), 2)
  a.ci90 <- paste0(round(quantile(draw_dat$A[, locateA17[i[k], i[k]]], probs = c(0.05)), 2), "-", 
                   round(quantile(draw_dat$A[, locateA17[i[k], i[k]]], probs = c(0.95)), 2))
  beta.all <- draw_dat$A[, locateA17[j[k], i[k]]]/draw_dat$A[, locateA17[i[k], i[k]]]
  beta.est <- round(median(beta.all), 2)
  beta.lower  <- round(quantile(beta.all, probs = c(0.05)), 2)
  beta.upper  <- round(quantile(beta.all, probs = c(0.95)), 2)
  beta.ci90 <- paste0(beta.lower, "-", beta.upper)
  temp <- cbind(FocalSp, R0.est, R0.ci90, a.est, a.ci90, Competitor, beta.est, beta.ci90)
  par.table <- rbind(par.table, temp)
}
par.table <- par.table[-1, ]
par.table.hot <- par.table %>% mutate(temperature = "Lowland")
write.csv(par.table, "competition model parameter_hot.csv")

# cold treatment
draw_dat <- rstan::extract(fit_cold, pars= c("mur", "A"))
i <- c(1,1,1,1,2,2,2,3,3,4,5,5) # all possibility of focal species
j <- c(2,3,4,5,1,3,5,1,2,1,1,2) # all possibility of competing species
par.table <- data.frame(FocalSp = NA, R0.est = NA, R0.ci90 = NA, a.est = NA, a.ci90 = NA, Competitor = NA, beta.est = NA, beta.ci90 = NA)
for (k in 1:length(i)){
  FocalSp  <- species[i[k]]
  Competitor <- species[j[k]]
  R0.est <- round(median(draw_dat$mur[,i[k]]), 2)
  R0.ci90 <- paste0(round(quantile(draw_dat$mur[,i[k]], probs = c(0.05)), 2), "-", 
                    round(quantile(draw_dat$mur[,i[k]], probs = c(0.95)), 2))
  a.est <- round(median(draw_dat$A[, locateA17[i[k], i[k]]]), 2)
  a.ci90 <- paste0(round(quantile(draw_dat$A[, locateA17[i[k], i[k]]], probs = c(0.05)), 2), "-", 
                   round(quantile(draw_dat$A[, locateA17[i[k], i[k]]], probs = c(0.95)), 2))
  beta.all <- draw_dat$A[, locateA17[j[k], i[k]]]/draw_dat$A[, locateA17[i[k], i[k]]]
  beta.est <- round(median(beta.all), 2)
  beta.lower  <- round(quantile(beta.all, probs = c(0.05)), 2)
  beta.upper  <- round(quantile(beta.all, probs = c(0.95)), 2)
  beta.ci90 <- paste0(beta.lower, "-", beta.upper)
  temp <- cbind(FocalSp, R0.est, R0.ci90, a.est, a.ci90, Competitor, beta.est, beta.ci90)
  par.table <- rbind(par.table, temp)
}
par.table <- par.table[-1, ]
par.table.cold <- par.table %>% mutate(temperature = "Upland")

## Table 2 - single-generation competition model parameters
par.table.all <- rbind(par.table.cold, par.table.hot)
write.csv(par.table, "competition model parameter.csv")


########################## infer coexistence status  #############################
## preparation: function to plot isocline plot
plot.coexis.HC <- function(pair.ind, pars){
  p <- pars[pair.ind, ]
  temp <- rbind(c("isolineX", 0, as.numeric(p$isoXonY)), c("isolineX", as.numeric(p$isoXonX), 0),
                c("isolineY", 0, as.numeric(p$isoYonY)), c("isolineY", as.numeric(p$isoYonX), 0))
  colnames(temp) <- c("type", "X", "Y")
  temp <- as.data.frame(temp)
  temp$X <- as.numeric(temp$X)
  temp$Y <- as.numeric(temp$Y)
  p <- temp %>%
    ggplot(aes(x = X, y = Y, group = type, col = type)) + geom_line(aes(group = type), size = 1) + 
    theme_bw() + theme(aspect.ratio=1) + 
    xlab(paste0("X(", p$spA, ")")) + ylab(paste0("Y(", p$spB, ")")) + 
    xlim(0, max(temp$X)*1.1) + ylim(0, max(temp$X)*1.1) + 
    scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) + 
    theme(axis.line = element_line(size = 0.5, colour = "black", linetype=1)) + 
    theme(axis.ticks = element_line(size = 0.5, color="black") , axis.ticks.length = unit(.5, "cm")) + 
    theme(legend.title = element_blank(), legend.position = c(0.7,0.8), legend.background = element_blank())
  return(p)
}

## cold treatment
# calculate all the parameters needed for infering equilibrium-based coexistence
draw_dat <- rstan::extract(fit_cold, pars= c("mur", "A"))
coexist.cold_HC <- data.frame(Pair = NA, spA = NA, spB = NA, 
                              R0.A = NA, R0.B = NA, 
                              a.A = NA, a.B = NA, 
                              beta.BA = NA, beta.AB = NA)
for (i in 1:6){
  spA.ind <- pairs[i,1]
  spB.ind <- pairs[i,2]
  spA <- species[spA.ind]
  spB <- species[spB.ind]
  Pair <- paste0(spA, "-", spB)
  R0.A <- median(draw_dat$mur[,spA.ind])
  R0.B <- median(draw_dat$mur[,spB.ind])
  a.A <- median(draw_dat$A[, locateA17[spA.ind, spA.ind]])
  a.B <- median(draw_dat$A[, locateA17[spB.ind, spB.ind]])
  beta.BA <- median(draw_dat$A[, locateA17[spB.ind, spA.ind]]/draw_dat$A[, locateA17[spA.ind, spA.ind]])
  beta.AB <- median(draw_dat$A[, locateA17[spA.ind, spB.ind]]/draw_dat$A[, locateA17[spB.ind, spB.ind]])
  temp <- cbind(Pair, spA, spB, R0.A , R0.B, a.A, a.B, beta.BA, beta.AB)
  coexist.cold_HC <- rbind(coexist.cold_HC, temp)
}
coexist.cold_HC <- coexist.cold_HC[-1,]
coexist.cold_HC <- coexist.cold_HC %>%
  mutate(R0.A = as.numeric(R0.A), R0.B = as.numeric(R0.B), a.A = as.numeric(a.A), 
         a.B = as.numeric(a.B), beta.BA = as.numeric(beta.BA), beta.AB = as.numeric(beta.AB))

# calculated derived values for ploting the isocline
coexist.cold_HC <- coexist.cold_HC %>%
  mutate(theta.A = 1/R0.A, theta.B = 1/R0.B, gamma.A = theta.A*a.A, gamma.B = theta.B*a.B) %>%
  mutate(stable.Co = ( (1-theta.A)/beta.BA > (1-theta.B) ) & ( (1-theta.B)/beta.AB > (1-theta.A) )  )
coexist.cold_HC <- coexist.cold_HC %>%
  mutate(isoXonY = (1-theta.A)/(gamma.A * beta.BA), isoXonX = (1-theta.A)/gamma.B, 
         isoYonX = (1-theta.B)/(gamma.B * beta.AB), isoYonY = (1-theta.B)/gamma.A)

# plot the isocline for the six pairs
p1 <- plot.coexis.HC(1, coexist.cold_HC)
p2 <- plot.coexis.HC(2, coexist.cold_HC)
p3 <- plot.coexis.HC(3, coexist.cold_HC)
p4 <- plot.coexis.HC(4, coexist.cold_HC)
p5 <- plot.coexis.HC(5, coexist.cold_HC)
p6 <- plot.coexis.HC(6, coexist.cold_HC)
## Table 2 last column - infer equilibrium state
plot_grid(p1, p2, p3, p4, p5, p6, ncol = 3, nrow = 2)
ggsave("coexistence_HassellComins_cold.png", width = 10)


## hot treatment
draw_dat <- rstan::extract(fit_hot, pars= c("mur", "A"))
coexist.hot_HC <- data.frame(Pair = NA, spA = NA, spB = NA, 
                             R0.A = NA, R0.B = NA, 
                             a.A = NA, a.B = NA, 
                             beta.BA = NA, beta.AB = NA)
for (i in 1:6){
  spA.ind <- pairs[i,1]
  spB.ind <- pairs[i,2]
  spA <- species[spA.ind]
  spB <- species[spB.ind]
  Pair <- paste0(spA, "-", spB)
  R0.A <- median(draw_dat$mur[,spA.ind])
  R0.B <- median(draw_dat$mur[,spB.ind])
  a.A <- median(draw_dat$A[, locateA17[spA.ind, spA.ind]])
  a.B <- median(draw_dat$A[, locateA17[spB.ind, spB.ind]])
  beta.BA <- median(draw_dat$A[, locateA17[spB.ind, spA.ind]]/draw_dat$A[, locateA17[spA.ind, spA.ind]])
  beta.AB <- median(draw_dat$A[, locateA17[spA.ind, spB.ind]]/draw_dat$A[, locateA17[spB.ind, spB.ind]])
  temp <- cbind(Pair, spA, spB, R0.A , R0.B, a.A, a.B, beta.BA, beta.AB)
  coexist.hot_HC <- rbind(coexist.hot_HC, temp)
}
coexist.hot_HC <- coexist.hot_HC[-1,]
coexist.hot_HC <- coexist.hot_HC %>%
  mutate(R0.A = as.numeric(R0.A), R0.B = as.numeric(R0.B), a.A = as.numeric(a.A), 
         a.B = as.numeric(a.B), beta.BA = as.numeric(beta.BA), beta.AB = as.numeric(beta.AB))
coexist.hot_HC <- coexist.hot_HC %>%
  mutate(theta.A = 1/R0.A, theta.B = 1/R0.B, gamma.A = theta.A*a.A, gamma.B = theta.B*a.B) %>%
  mutate(stable.Co = ( (1-theta.A)/beta.BA > (1-theta.B) ) & ( (1-theta.B)/beta.AB > (1-theta.A) )  )
coexist.hot_HC <- coexist.hot_HC %>%
  mutate(isoXonY = (1-theta.A)/(gamma.A * beta.BA), isoXonX = (1-theta.A)/gamma.B, 
         isoYonX = (1-theta.B)/(gamma.B * beta.AB), isoYonY = (1-theta.B)/gamma.A)

p1 <- plot.coexis.HC(1, coexist.hot_HC)
p2 <- plot.coexis.HC(2, coexist.hot_HC)
p3 <- plot.coexis.HC(3, coexist.hot_HC)
p4 <- plot.coexis.HC(4, coexist.hot_HC)
p5 <- plot.coexis.HC(5, coexist.hot_HC)
p6 <- plot.coexis.HC(6, coexist.hot_HC)

## Table 2 last column - infer equilibrium state
# for any pairs that contains PAL or PST, as their R0 is smaller than 1, their population will eventually die out.
plot_grid(p2, p4, ncol = 1, nrow = 2)
ggsave("coexistence_HassellComins_hot.png", width = 3.3)



## 5. long-term competition
########################## data formating  #############################
## preparing data
longdat <- read.csv("Data/long_term_competition_PAN-PAL.csv")
longdat <- longdat %>%
  mutate(compo = ifelse(intra, "single", "mix"), gr = paste0(species, "_", compo))
longdat$compo <- factor(longdat$compo, c("single", "mix"))
longdat$gr <- factor(longdat$gr, c("pandora_single", "pallidifrons_single", "pandora_mix", "pallidifrons_mix"))

longdat.g <- longdat %>%
  group_by(vialID, treatment, block, intra, species, rep, compo, gr) %>%
  summarize(totalN = sum(count))

longdat.g <- longdat.g %>%
  mutate(gr2 = paste0(treatment, "_", compo))
longdat.g$gr2 <- factor(longdat.g$gr2, c("HOT_single","HOT_mix","COLD_single","COLD_mix"))
longdat.g$compo <- factor(longdat.g$compo, c("single", "mix"))


########################## visualization: reaction norm  #############################
# Figure 5 - reaction norms plot: plot both species and competition treatments together
longdat.g %>%
  ggplot(aes(x = treatment, y = totalN, group = gr2, color = species)) + 
  geom_jitter(aes(shape = compo), width = 0.08, height = 0, alpha = 0.6, size = 3) +
  stat_summary(fun = mean, geom = "line", aes(group = gr, linetype = compo, color = species), size = 1) + 
  scale_color_manual("Species", values = c('pallidifrons' = 'darkblue', 'pandora' = 'red')) +
  scale_shape_manual("Competition", values = c('single' = 1, 'mix' = 17)) + 
  scale_linetype_manual("Competition", values = c("single" = "dotted", "mix" = "solid") ) + 
  xlab("Temperature") + ylab("Population size") +
  theme(legend.position = c(0.8, 0.3)) + 
  theme(axis.title.y = element_text(size = 18), axis.title.x = element_text(size = 18)) +
  theme(axis.text.y = element_text(size = 15), axis.text.x = element_text(size = 15)) +
  theme(legend.text=element_text(size=12), legend.title = element_text(size=12)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave("figure5_longterm competition.png", width = 4, height = 7)

########################## examining the three way interation  #############################
## statistical test
# tried negbinomial(), zero_inflated_negbinomial(), zero_inflated_lpoisson(). Convergence is the best with zero_inflated_negbinomial(). 
longdat.g %>%
  brm(formula = totalN ~ 1 + treatment * species * compo + (1|vialID), 
      family = zero_inflated_negbinomial(link = "log"), 
      control = list(adapt_delta = 0.95, max_treedepth = 20)) -> fit2
save(fit, file = 'StanFits/fit_NegB_PANPAL_longterm')
load('StanFits/fit_NegB_PANPAL_longterm')

## Supplementarty figure 4C - diagnostics
# convergence - acceptable ESS
print(fit)
# est vs. obs
yrep <- posterior_predict(fit)
dim(yrep)
ppc_dens_overlay(longdat.g$totalN, yrep[1:100, ])
ggsave("longterm_model_diagnosis.png")

## the effect of temperature to population size with or without competitors
## Figure 5b - the effect of hot temperature
# extract every posterior that is related to temperature (other will only contribute to a varying intercept, but we just focus of the effect size (slope) of temperature)
b_HOT <- posterior_samples(fit, pars = c("treatmentHOT"))
l <- length(b_HOT[,1])
temp.effect <- data.frame(sp = NA, competition = NA, effect = NA)
temp.effect <- rbind(temp.effect, cbind(sp = rep("PAL", l), competition = rep("Single", l), 
                                        effect = b_HOT[,1]))
temp.effect <- rbind(temp.effect, cbind(sp = rep("PAL", l), competition = rep("Mix", l), 
                                        effect = b_HOT[,1] + b_HOT[,3]))
temp.effect <- rbind(temp.effect, cbind(sp = rep("PAN", l), competition = rep("Single", l), 
                                        effect = b_HOT[,1] + b_HOT[,2]))
temp.effect <- rbind(temp.effect, cbind(sp = rep("PAN", l), competition = rep("Mix", l), 
                                        effect = b_HOT[,1] + b_HOT[,2] + b_HOT[,3] + b_HOT[,4]))

temp.effect <- temp.effect[-1, ]
temp.effect$effect <- as.numeric(temp.effect$effect)

p.pal <- temp.effect %>%
  filter(sp == "PAL") %>%
  ggplot(aes(x = competition, y = effect)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") + 
  geom_violin(trim = FALSE, fill = "gray") + coord_flip() +
  ylim(-20, 5) + ylab("Effect of high temperature") + 
  xlab("") + ggtitle("D. pallidifrons") +
  theme(axis.title.y = element_text(size = 13), axis.title.x = element_text(size = 13)) +
  theme(axis.text.y = element_text(size = 13, angle = 45), axis.text.x = element_text(size = 13)) +
  theme(plot.title = element_text(size=18)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

p.pan <- temp.effect %>%
  filter(sp == "PAN") %>%
  ggplot(aes(x = competition, y = effect)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") + 
  geom_violin(trim = FALSE, fill = "gray") + coord_flip() +
  ylim(-1, 2) + ylab("Effect of high temperature") + 
  xlab("") + ggtitle("D. pandora") + 
  theme(axis.title.y = element_text(size = 13), axis.title.x = element_text(size = 13)) +
  theme(axis.text.y = element_text(size = 13, angle = 45), axis.text.x = element_text(size = 13)) +
  theme(plot.title = element_text(size=18)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

plot_grid(p.pal, p.pan, ncol = 1, nrow = 2)
ggsave("effect_temperature_longterm.png", width = 3.5, height = 7)

## the effect of competition for different species in hot or cold temperature:
## Figure 5c - the effect of competition
# extract every posterior that is related to competition (single vs. mix)
b_MIX <- posterior_samples(fit, pars = c("compomix"))
l <- length(b_MIX[,1])
mix.effect <- data.frame(sp = NA, temperature = NA, effect = NA)
mix.effect <- rbind(mix.effect, cbind(sp = rep("PAL", l), temperature = rep("Cold", l), 
                                      effect = b_MIX[,1]))
mix.effect <- rbind(mix.effect, cbind(sp = rep("PAL", l), temperature = rep("Hot", l), 
                                      effect = b_MIX[,1] + b_MIX[,2]))
mix.effect <- rbind(mix.effect, cbind(sp = rep("PAN", l), temperature = rep("Cold", l), 
                                      effect = b_MIX[,1] + b_MIX[,3]))
mix.effect <- rbind(mix.effect, cbind(sp = rep("PAN", l), temperature = rep("Hot", l), 
                                      effect = b_MIX[,1] + b_MIX[,2] + b_MIX[,3] + b_MIX[,4]))

mix.effect <- mix.effect[-1, ]
mix.effect$effect <- as.numeric(mix.effect$effect)

p.pal2 <- mix.effect %>%
  filter(sp == "PAL") %>%
  ggplot(aes(x = temperature, y = effect)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") + 
  geom_violin(trim = FALSE, fill = "gray") + coord_flip() +
  ylab("Effect of competition") + ylim(-10, 1) + 
  xlab("") + ggtitle("D. pallidifrons") + 
  theme(axis.title.y = element_text(size = 13), axis.title.x = element_text(size = 13)) +
  theme(axis.text.y = element_text(size = 13, angle = 45), axis.text.x = element_text(size = 13)) +
  theme(plot.title = element_text(size=18)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

p.pan2 <- mix.effect %>%
  filter(sp == "PAN") %>%
  ggplot(aes(x = temperature, y = effect)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") + 
  geom_violin(trim = FALSE, fill = "gray") + coord_flip() +
  ylab("Effect of competition") + ylim(-2.2, 0.5) + 
  xlab("") + ggtitle("D. pandora") + 
  theme(axis.title.y = element_text(size = 13), axis.title.x = element_text(size = 13)) +
  theme(axis.text.y = element_text(size = 13, angle = 45), axis.text.x = element_text(size = 13)) +
  theme(plot.title = element_text(size=18)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

plot_grid(p.pal2, p.pan2, ncol = 1, nrow = 2)
ggsave("effect_competition_longterm.png", width = 3.5, height = 7)


