## load packages
library(ape)
library(wesanderson)
library(brms)
library(readxl)
library(paleotree)
library(phytools)
library(mice)
library(rstan)
library(MCMCglmm)
library(bayesplot)
library(rethinking)


# set working directory
setwd("C:/Model S5 Ancestral state estimation of all marsupials")

########################################################## PREPARE DATA###########################################################
d_mas<- as.data.frame(read_excel ("S5.xlsx"))

## prepare phylogeny
mas4<-read.nexus("outputS5.nex")

### remove species without any data(normally in this excel all sp has info)
d_mas2<- d_mas[!is.na(d_mas$social_state),]

#### check species names match and trim trees to dataset
setdiff(d_mas2$Genus_species, mas4[[1]]$tip.label) # character(0)

##### prune tree to only retain species in the dataset, for each of the 100 trees in the sample
mas.prun<- mas4
for(i in 1:100){
  mas.prun[[i]]<- drop.tip(mas.prun[[i]], setdiff(mas.prun[[i]]$tip.label,d_mas2$Genus_species))
}

###### check again whether ultrametric
ultra<- rep(NA, 100)
for(i in 1:100){
  ultra[i]<- is.ultrametric(mas.prun[[i]])
}
sum(ultra) # all false

# try if it works nonetheless
# convert to covariance matrix (see https://cran.r-project.org/web/packages/brms/vignettes/brms_phylogenetics.html)
A<- list() 
for(i in 1:100){
  A[[i]]<- vcv.phylo(mas.prun[[i]])
}


###################################################################################################################################
############################################# FIT PHYLOGENTIC MODELS AND CORRECT MISSPLEING #######################################
###################################################################################################################################
d_mas2$phylo<- d_mas2$Genus_species #rename the Genus_species column (that is our phylogeny)
d_mas2$social_state<- relevel(as.factor(d_mas2$social_state), ref="S") #S as reference so it's not in the following prior


#Center predictors
d_mas2$Habitat_hetero<-(d_mas2$Habitat_hetero)-1
d_mas2$Nbr_studies<-(d_mas2$Nbr_studies)-1
d_mas2$logBM.foss<- (log(d_mas2[,"body_mass(g)"])-log(37.75)/sd(log(d_mas2[,"body_mass(g)"]), na.rm=TRUE)) # 37.75 is the estimated body mass of marsupial ancestor (Djarthia murgonensis)


#model4
bform<- bf(social_state~1 + Nbr_studies+ Habitat_hetero + logBM.foss +(1|phylo)+(1|Genus_species))

get_prior(bform, data=d_mas2, family=categorical)
prior=c(prior(normal(0,10), class=Intercept),
        prior(cauchy(0,2), class=sd, dpar="muP"), 
        prior(cauchy(0,2), class=sd, dpar="muG"),
        prior(cauchy(0,2), class=sd, dpar="muSP"),
        prior(cauchy(0,2), class=sd, dpar="muPG"),
        prior(cauchy(0,2), class=sd, dpar="muSPG"),
        prior(cauchy(0,2), class=sd, dpar="muVG"),
        prior(normal(0,1), class=b))

m3_mas<- brm(bform, data=d_mas2, family=categorical, cov_ranef = list(phylo = A[[1]]), prior = prior, 
             chains = 2, cores = 2, iter = 2000, warmup = 1000, thin=5, control = list(adapt_delta = 0.99, max_treedepth = 15))

summary(m3_mas)

#save result(1 repeat)
save(m3_mas,file="m3_mas.RData")
load("m3_mas.RData")

post_m3_mas<- posterior_samples(m3_mas)

str(post_m3_mas[,c(1:36)])

post_m3_all_mas<- post_m3_mas[,c(1:36)]



#########################################repeat 100 times###############################################
for(i in 2:100){
  m3_i_mas<- brm(bform, data=d_mas2, family=categorical, cov_ranef = list(phylo = A[[i]]), prior = prior,
                 chains = 2, cores = 2, iter = 2000, warmup = 1000, thin=5, control = list(adapt_delta = 0.99, max_treedepth = 15))
  post_m3_mas<-posterior_samples(m3_i_mas)
  post_m3_all_mas<- rbind(post_m3_all_mas, post_m3_mas[,c(1:36)])
  save(post_m3_all_mas, file="post_m3_all_mas.robj")
}



save(post_m3_all_mas, file="post_m3_all_mas.RData")
save(post_m3_mas, file="post_m3_mas.RData")
save(m3_i_mas, file="m3_i_mas.RData") 
load("post_m3_all_mas.RData")
load("post_m3_mas.RData")
load("m3_i_mas.RData")

nrow(post_m3_all_mas) # 40 000 --> ok

summary(m3_i_mas)



# phylogenetic signal: sum of all SO variances / sum of all variance components + distribution-specific variance
summary(post_m3_all_mas)
VarPhy.P<-post_m3_all_mas$sd_phylo__muP_Intercept 
VarPhy.G<-post_m3_all_mas$sd_phylo__muG_Intercept
VarPhy.SP<-post_m3_all_mas$sd_phylo__muSP_Intercept 
VarPhy.PG<-post_m3_all_mas$sd_phylo__muPG_Intercept 
VarPhy.SPG<-post_m3_all_mas$sd_phylo__muSPG_Intercept
VarPhy.VG<-post_m3_all_mas$sd_phylo__muVG_Intercept


VarSpec.P<-post_m3_all_mas$sd_Genus_species__muP_Intercept                            
VarSpec.G<-post_m3_all_mas$sd_Genus_species__muG_Intercept                             
VarSpec.SP<-post_m3_all_mas$sd_Genus_species__muSP_Intercept                             
VarSpec.PG<-post_m3_all_mas$sd_Genus_species__muPG_Intercept                            
VarSpec.SPG<-post_m3_all_mas$sd_Genus_species__muSPG_Intercept                            
VarSpec.VG<-post_m3_all_mas$sd_Genus_species__muVG_Intercept    

VarDistro<- pi^2 / 3

lambda<- (VarPhy.P+VarPhy.G+VarPhy.SP+VarPhy.PG+VarPhy.SPG+VarPhy.VG)/
  (VarPhy.P+VarPhy.G+VarPhy.SP+VarPhy.PG+VarPhy.SPG+VarPhy.VG
   +VarSpec.P+VarSpec.G+VarSpec.SP+VarSpec.PG+VarSpec.SPG+VarSpec.VG
   +VarDistro)
mean(lambda); median(lambda); HPDI(lambda, prob=0.95) 
#[1] 0.1780213
#[1] 0.168798
#|0.95      0.95| 
#  0.06470337 0.31212483

{    K <- 7 #number of social states
  ns <- nrow(post_m3_all_mas)
  n <- 1
  
  softmax2 <- function(x) {
    x <- max(x) - x
    exp(-x)/sum(exp(-x))
  }
  
  p <- list()
  
  for ( i in 1:n ) {
    p[[i]] <- sapply( 1:K , function(k) {
      if ( k < K ) {
        ptemp <- post_m3_all_mas[,k] 
      } else {
        ptemp <- rep(0,ns)
      }
      return(ptemp)
    })
    ## The values are converted to probabilities using the softmax function
    ## which ensures that the predicted values across categories sum to
    ## 100% probabilities.
    for ( s in 1:ns ) p[[i]][s,] <- softmax2( p[[i]][s,] )
  }
}

#probabilities estimations at root:
p_mean.brms_mas <- sapply( 1:length(p) , function(i) apply(p[[i]],2,mean) )
p_HPDI.brms_mas <- sapply( 1:length(p) , function(i) apply(p[[i]],2,HPDI, prob=0.95) )
pred_probs_mod5_mas.brms<- cbind(p_mean.brms_mas, p_HPDI.brms_mas[c(1,3,5,7,9,11,13),], p_HPDI.brms_mas[c(2,4,6,8,10,12,14),])
rownames(pred_probs_mod5_mas.brms)<- c("G","P","PG","SP","SPG","VG","S")
colnames(pred_probs_mod5_mas.brms)<- c("mean", "lwr95", "upr95")
round(pred_probs_mod5_mas.brms,2)
#mean lwr95 upr95
#G   0.09     0  0.53
#P   0.22     0  0.90
#PG  0.25     0  0.91
#SP  0.03     0  0.15
#SPG 0.05     0  0.27
#VG  0.02     0  0.06
#S   0.35     0  0.88



######################################################################################################
#plot 1 - social organization probabilities at root
par(mar = c(3, 6, 4, 1), xpd=NA)
plot(c(0,0.14, 0.28,0.43,0.57,0.71,1)~c(1,2,3,4,5,6,7), col="white", xaxt="n", ylab="Probability at root", xlab="", main="",cex=0.05)
axis(1, at=c(1,2,3,4,5,6,7), labels=c("S", "P", "G","SP","PG","SPG","VG"),cex=0.05)


# S on line 2
points(1, pred_probs_mod5_mas.brms[7,1], pch=16,  col="darkorange", cex=2)
arrows(1, pred_probs_mod5_mas.brms[7,2], 1, pred_probs_mod5_mas.brms[7,3], length=0, angle=90, col="darkorange", lwd=2)

# P on line 3
points(2, pred_probs_mod5_mas.brms[2,1], pch=16, col="darkorchid", cex=2)
arrows(2, pred_probs_mod5_mas.brms[2,2], 2, pred_probs_mod5_mas.brms[2,3], length=0, angle=90, col="darkorchid", lwd=2)

#G on line 1
points(3, pred_probs_mod5_mas.brms[1,1], pch=16, col="cornflowerblue", cex=2)
arrows(3, pred_probs_mod5_mas.brms[1,2], 3, pred_probs_mod5_mas.brms[1,3], length=0, angle=90, col="cornflowerblue", lwd=2)

# SP on line 2
points(4, pred_probs_mod5_mas.brms[4,1], pch=16,  col="brown", cex=2)
arrows(4, pred_probs_mod5_mas.brms[4,2], 4, pred_probs_mod5_mas.brms[4,3], length=0, angle=90, col="brown", lwd=2)

#PG on line 1
points(5, pred_probs_mod5_mas.brms[3,1], pch=16, col="blue4", cex=2)
arrows(5, pred_probs_mod5_mas.brms[3,2], 5, pred_probs_mod5_mas.brms[3,3], length=0, angle=90, col="blue4", lwd=2)

# SPG on line 2
points(6, pred_probs_mod5_mas.brms[5,1], pch=16,  col="black", cex=2)
arrows(6, pred_probs_mod5_mas.brms[5,2], 6, pred_probs_mod5_mas.brms[5,3], length=0, angle=90, col="black", lwd=2)

# VG on line 3
points(7, pred_probs_mod5_mas.brms[6,1], pch=16, col="burlywood", cex=2)
arrows(7, pred_probs_mod5_mas.brms[6,2], 7, pred_probs_mod5_mas.brms[6,3], length=0, angle=90, col="burlywood", lwd=2)



