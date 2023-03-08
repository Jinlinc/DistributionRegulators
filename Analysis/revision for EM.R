######### Distribution manuscript additional analysis for revision #########
######### Jinlin Chen #########
######### Last update: 2023 Jan #########

## note: these codes need to run together with the main_results_analysis to set the global environment ##

## the coolest temperature ##
climate <- read.csv(file.choose())
summary(climate)
climate %>% 
  filter(Site2016=="P780"| Site2016=="K730") -> climate.highland
quantile(climate.highland$Celsius,0.01)
summary(climate.highland$Celsius)

## CTmin analysis ## 
ctmin <- read.csv(file.choose())
plot(ctmin$Ctmin.in.Kellermann ~ ctmin$hIndex)
cor(ctmin$Ctmin.in.Kellermann, ctmin$CTmin.in.Overgaard, method = c("pearson"))
cor(ctmin$Ctmin.in.Kellermann, ctmin$RTmin.in.ours, method = c("pearson"))
cor(ctmin$CTmin.in.Overgaard, ctmin$RTmin.in.ours, method = c("pearson"))
cor(ctmin$Ctmin.in.Kellermann, ctmin$hIndex, method = c("pearson"))
cor(ctmin$CTmin.in.Overgaard, ctmin$hIndex, method = c("pearson"))
lm(ctmin$CTmin.in.Overgaard ~ ctmin$hIndex) %>% summary()


## distribution analysis ##
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
                         ifelse(elev.s == "Medium", 0.5, 1))) %>%
  mutate(Host = ifelse(Host == "palidifrons", "pallidifrons", Host))  # correct the mis-spealing of the species name
dist.dat.core$Host <- factor(dist.dat.core$Host, levels = c("bunnanda", "pandora", "bipectinata", "pseudoananassae", "sulfurigaster", "rubida", "birchii", "pallidifrons", "pseudotakahashii"))
dist.dat.core$Site <- factor(dist.dat.core$Site, levels = c("K070", "P070", "P350", "K390", "K730", "P880"))
dist.dat.core$elev.s <- factor(dist.dat.core$elev.s, levels = c("High", "Medium", "Low"))

# ranks of weighted centre elevation in two transects
dist.dat.core %>%  
  filter(Transect == "Paluma") %>%
  group_by(Host) %>%
  summarise(h.mean = mean(elevation), hIndex.mean = mean(hIndex)) -> h.index.Paluma
# rank
#[1] bunnanda         bipectinata      pandora          pseudoananassae  sulfurigaster   
#[6] rubida           birchii          pallidifrons     pseudotakahashii

dist.dat.core %>%  
  filter(Transect == "Kirrama") %>%
  group_by(Host) %>%
  summarise(h.mean = mean(elevation), hIndex.mean = mean(hIndex)) -> h.index.Kirrama
# rank
#[1] bunnanda         pandora          bipectinata      sulfurigaster     rubida
#[6] pseudoananassae           birchii          pallidifrons     pseudotakahashii

pars.est_ci_2 <- read.csv(file.choose())
pars.est_ci_2 <- pars.est_ci_2 %>% filter(t.tree != "melanogaster" & t.tree != "simulans")  
pars.est_ci_2$hIndex.P <- as.numeric(pars.est_ci_2$hIndex.P)
pars.est_ci_2$hIndex.K <- as.numeric(pars.est_ci_2$hIndex.K)

ggplot(data = pars.est_ci_2, aes(x = hIndex.P, y = RTmax.med)) + geom_point()
lm(data = pars.est_ci_2, RTmax.med ~ hIndex.P) %>% summary 

ggplot(data = pars.est_ci_2, aes(x = hIndex.K, y = RTmax.med)) + geom_point()
lm(data = pars.est_ci_2, RTmax.med ~ hIndex.K) %>% summary 


## alternatives to hIndex - do they change the regression conclusion? ##
# The previous way of calculating hIndex across two transects need modification.
# I should not calculate the mean directly from both transect (as the transect differ in their sample size, one distribution pattern will be over represent)
# I should calculate the mean for each of transect, then average them.
dist.dat.core <- distribution.dat %>%
  filter(is.na(Host) != TRUE) %>%
  filter(is.na(Site) != TRUE) %>%
  filter(Host != "Drosophila sp. 1" & Host != "sulfurigaster complex" & Host != "serrata" & Host != "immigrans") %>%
  mutate(elevation = as.numeric(substr(Site,2,4))) %>%
  mutate(t_sp = paste0(Host, "_", Transect)) %>%
  mutate(elev.s = ifelse(Site == "K070" | Site == "P070", "Low",
                         ifelse(Site == "K390" | Site == "P350", "Medium", "High"))) %>%
  mutate(hIndex = ifelse(elev.s == "Low", 0,
                         ifelse(elev.s == "High", 1, 
                                ifelse(Site == "P350", 350/(880-70), 390/(730-70))))) %>%
  mutate(Host = ifelse(Host == "palidifrons", "pallidifrons", Host))  # correct the mis-spealing of the species name

dist.dat.core %>%  
  group_by(Host, Transect) %>%
  summarise(h.mean = mean(elevation), hIndex.mean = mean(hIndex)) %>%
  ungroup(Host, Transect) %>%
  group_by(Host) %>%
  summarise(h.mean = mean(h.mean), hIndex.mean = mean(hIndex.mean)) -> h.index.3

## regression by the equidistant hIndex
swfun.hIndex2 <- function(x){
  switch (x,
          "bipectinata" = 0.163, "birchii" = 0.651, "bunnanda" = 0,
          "melanogaster" = NA, "pallidifrons" = 0.771, "pandora" = 0.039,
          "pseudoananassae" = 0.341, "simulans" = NA, "sulfurigaster" = 0.419
  )
}
# the rest regression is using the main code

## regression by the actual km elevation
swfun.hIndex3 <- function(x){
  switch (x,
          "bipectinata" = 0.169, "birchii" = 0.498, "bunnanda" = 0.070,
          "melanogaster" = NA, "pallidifrons" = 0.617, "pandora" = 0.095,
          "pseudoananassae" = 0.280, "simulans" = NA, "sulfurigaster" = 0.367
  )
}
# ... main code
RT.traits2$hMean <- sapply(as.character(RT.traits2$sp.list), swfun.hIndex3)
# ... main code to prepare data
fit.RTmax.hMean <- brm(
  RTmax.med ~ hMean + (1|gr(t.tree, cov = A)),
  data = RT.traits2, 
  family = gaussian(),
  data2 = list(A = A),
  control = list(adapt_delta = 0.95))
# Result-4.Heat tolerance
# "species whose distribution were biased towards lowland consistently had higher RTmax (Figure 3e. Coefficient = -3.06, 95% credible interval = - 5.30 – - 0.88). "
summary(fit.RTmax.hMean)


## figure 3 modification ##
cc.dist.type <- c("pandora" = "#FF9933", "bipectinata" = "#FF9933", "bunnanda" = "#FF9933",
                  "pseudoananassae"="#025139", "sulfurigaster" = "#025139", "birchii" = "#025139",
                  "pallidifrons" = "#004C99")

## RTmin ~ hIndex plot ##
p1.1.revision <- RT.traits2 %>%
  filter(sp.list != "melanogaster" & sp.list != "simulans") %>%
  ggplot(aes(x = hIndex, y = RTmin.med, color = sp.list)) + geom_point(size = 2) + 
  geom_errorbar(aes(ymin=RTmin.lower, ymax=RTmin.upper), width=.2) + 
  scale_color_manual(values = cc.dist.type) + theme(legend.position = "none") + 
  ylim(12, 18) + ylab("RTmin(°C)") + xlab("hIndex") + 
  theme(axis.title.y = element_text(size = 10), axis.text.y = element_text(size = 10)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave("RTmin_regression_revision.png", width = 4, height = 3)


## figure 2 modification ##
# plot Topt ~ species ordered by hIndex again
cc.tpc <- c("melanogaster" = "grey",
            "simulans" = "#663300", "pandora" = "#FF9933", "bipectinata" = "#FFCC99", "bunnanda" = "#CC6600",
            "pseudoananassae"="#BDE51F", "sulfurigaster" = "#5FCB0B",
            "pallidifrons" = "#004C99", "birchii" = "#025139")
pars.est_ci %>%
  filter(sp.list != "melanogaster" & sp.list != "simulans") %>%
  ggplot(aes(x = sp.list, y = RTopt.med, color = sp.list)) + geom_point(size = 2) + 
  geom_errorbar(aes(ymin=RTopt.lower, ymax=RTopt.upper), width=.2) + 
  scale_color_manual(values = cc.tpc) + theme(legend.position = "none") + 
  ylim(22, 27) + ylab("RTopt(°C)") + xlab("Drosophila species") + 
  theme(axis.title.y = element_text(size = 10), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10, face = "italic"), axis.text.x = element_text(size = 10, angle = 45, hjust = 1, face = "italic")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave("RTopt_color.png", width = 3, height = 3)


## Conservative test: RTmax ~ distribution rank ##
pars.est_ci <- read.csv(file.choose())

pars.est_ci %>%
  filter(t.tree != "melanogaster" & t.tree != "simulans") %>%
  mutate(distRank = as.numeric(distRank)) %>%
  lm(formula = RTmax.med ~ distRank) %>% summary()
pars.est_ci %>%
  filter(t.tree != "melanogaster" & t.tree != "simulans") %>%
  mutate(distRank = as.numeric(distRank)) %>%
  lm(formula = RTmin.med ~ distRank) %>% summary()


## thermal performance curves with ci ##
temperature <- seq(10, 33, 0.1)
prediction.tpc <- function(sp.ind, pars) {
  drawN <- 200
  
  # ci
  pred.sp <- matrix(data = NA, nrow = drawN, ncol = length(temperature))
  for (i in 1:drawN){
    a.current <- pars$a[i, sp.ind]
    b.current <- pars$b[i, sp.ind]
    RTmax.current <- pars$RTmax[i, sp.ind]
    RTmin.current <- pars$RTmin[i, sp.ind]
    for (j in 1:length(temperature)){
      t <- temperature[j]
      pred.sp[i, j] <- ifelse(t <= RTmin.current, 0, 
                                  ifelse(t >= RTmax.current, 0, 
                                         a.current*t*(t - RTmin.current)*(RTmax.current - t)^(1/b.current)))
      pred.sp[i, j] <- pred.sp[i, j]^2
    }
  }
  pred.sp.lower <- apply(pred.sp, 2, quantile, probs=c(0.05))
  pred.sp.upper <- apply(pred.sp, 2, quantile, probs=c(0.95))
  
  # est of the fitted line
  a.med <- median(pars$a[, sp.ind])
  b.med <- median(pars$b[, sp.ind])
  RTmax.med <- median(pars$RTmax[, sp.ind])
  RTmin.med <- median(pars$RTmin[, sp.ind])
  est.sp <- rep(NA, length(temperature))
  for (j in 1:length(temperature)){
    t <- temperature[j]
    est.sp[j] <- ifelse(t <= RTmin.med, 0, 
                            ifelse(t >= RTmax.med, 0, 
                                   a.med*t*(t - RTmin.med)*(RTmax.med - t)^(1/b.med)))
    est.sp[j] <- est.sp[j]^2
  }
  pred.sp <- data.frame(sp = sp.list[sp.ind], temp = temperature, est = est.sp, 
                        upper.pred = pred.sp.upper, lower.pred = pred.sp.lower)
  return(pred.sp)
}

pred.all <- data.frame(sp = NA, temp = NA, est = NA, upper.pred = NA, lower.pred = NA)
for (sp in 1:9){
  pred.spx <- prediction.tpc(sp, pars)
  pred.all <- rbind(pred.all, pred.spx)
}
pred.all <- pred.all[-1,]
pred.all <- pred.all %>% mutate(spName = sp)
View(pred.all)

# Supplementary Figure 2C
tpc_data_RS %>%
  ggplot(aes(x = temp.corr, y = dailyRS)) + geom_point() + 
  geom_ribbon(data = pred.all, aes(x = temp, y = est, ymin=lower.pred, ymax=upper.pred), alpha = 0.6) +
  geom_path(data = pred.all, aes(x = temp, y = est, color = "blue")) + scale_color_manual(values = "blue") +  
  xlim(10,33) + xlab("Temperature (°C)") + ylab("Productivity (per female parent, per day)") + ylim(0,32) + 
  facet_grid(. ~ spName) + 
  theme(axis.title.y = element_text(size = 18), axis.title.x = element_text(size = 18)) +
  theme(axis.text.y = element_text(size = 12), axis.text.x = element_text(angle = 45, size = 12)) +
  theme(strip.text.x = element_text(size = 9, face = "italic")) + 
  theme(legend.position = "none") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave("tpc_summary_ci.png", width = 12, height = 5)


## figure 5 ##
load('StanFits/fit_NegB_PANPAL_longterm')
longdat.g %>%
  ggplot(aes(x = temperature, y = totalN, group = compo, color = species)) + 
  geom_point(aes(shape = compo), alpha = 0.6, size = 3) +
  stat_summary(fun = mean, geom = "line", aes(group = gr, linetype = compo, color = species), size = 1) + 
  scale_color_manual("Species", values = c('pallidifrons' = 'darkblue', 'pandora' = 'red')) +
  scale_shape_manual("Competition", values = c('monoculture' = 1, 'mixed' = 17)) + 
  scale_linetype_manual("Competition", values = c("monoculture" = "dotted", "mixed" = "solid") ) + 
  xlab("Temperature regime") + ylab("Population size") +
  facet_grid(.~ species) + 
  theme(legend.position = c(0.15, 0.8)) + 
  theme(strip.text.x = element_text(size = 15)) + 
  theme(axis.title.y = element_text(size = 18), axis.title.x = element_text(size = 18)) +
  theme(axis.text.y = element_text(size = 15), axis.text.x = element_text(size = 15)) +
  theme(legend.text=element_text(size=12), legend.title = element_text(size=12)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave("figure5_longterm competition_simple.png", width = 8, height = 6)
  




