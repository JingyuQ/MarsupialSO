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
library(brms)
library(bayesplot)
library(rethinking)

# set working directory
setwd("C:/Model S7 Ancestral state estimation on the species level with main social organisation of all marsupials")

########################################################## PREPARE DATA###########################################################
d_mas<- as.data.frame(read_excel ("S7.xlsx"))

## prepare phylogeny
mas4<-read.nexus("outputS7.nex")

### remove species without any data
d_mas2<- d_mas[!is.na(d_mas$Main_SO),]


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
d_mas2$Main_SO<- relevel(as.factor(d_mas2$Main_SO), ref="G") #G as reference so it's not in the following prior

#center predictor
d_mas2$Nbr_studies<-(d_mas2$Nbr_studies)-1

#model2
bform<- bf(Main_SO~1+Nbr_studies+(1|phylo)+(1|Genus_species))

get_prior(bform, data=d_mas2, family=categorical)
prior=c(prior(normal(0,10), class=Intercept),
        prior(cauchy(0,2), class=sd, dpar="muS"), 
        prior(cauchy(0,2), class=sd, dpar="muP"),
        prior(normal(0,1), class=b))

m2_mas<- brm(bform, data=d_mas2, family=categorical, cov_ranef = list(phylo = A[[1]]), prior = prior, 
                       chains = 2, cores = 2, iter = 2000, warmup = 1000, thin=5, control = list(adapt_delta = 0.99, max_treedepth = 15))
summary(m2_mas)

#save result(1 repeat)
save(m2_mas,file="m2_mas.RData")
load("m2_mas.RData")


post_m2_mas<- posterior_samples(m2_mas)

str(post_m2_mas[,c(1:20)])

post_m2_all_mas<- post_m2_mas[,c(1:20)]



#############repeat 100 times#####################
for(i in 2:100){
  m2_i_mas<- brm(bform, data=d_mas2, family=categorical, cov_ranef = list(phylo = A[[i]]), prior = prior,
                          chains = 2, cores = 2, iter = 2000, warmup = 1000, thin=5, control = list(adapt_delta = 0.99, max_treedepth = 15))
  post_m2_mas<-posterior_samples(m2_i_mas)
  post_m2_all_mas<- rbind(post_m2_all_mas, post_m2_mas[,c(1:20)])
}


save(post_m2_all_mas, file="post_m2_all_mas.RData")
save(post_m2_mas, file="post_m2_mas.RData")
save(m2_i_mas, file="m2_i_mas.RData")
load("post_m2_all_mas.RData")
load("post_m2_mas.RData")
load("m2_i_mas.RData")

nrow(post_m2_all_mas) # 40 000 --> ok

summary(m2_i_mas)
summary(post_m2_all_mas)

VarPhy.S<-post_m2_all_mas$sd_phylo__muS_Intercept#no references(in this case Group-living)
VarPhy.P<-post_m2_all_mas$sd_phylo__muP_Intercept


VarSpec.S<-post_m2_all_mas$sd_Genus_species__muS_Intercept                            
VarSpec.P<-post_m2_all_mas$sd_Genus_species__muP_Intercept                            

VarDistro<- pi^2/3

lambda<- (VarPhy.S+VarPhy.P)/(VarPhy.S+VarPhy.P+VarSpec.S+VarSpec.P+VarDistro)
mean(lambda); median(lambda); HPDI(lambda, prob=0.95)



{    K <- 3 #number of social states used
  ns <- nrow(post_m2_all_mas)
  n <- 1
  
  softmax2 <- function(x) {
    x <- max(x) - x
    exp(-x)/sum(exp(-x))
  }
  
  p <- list()
  
  for ( i in 1:n ) {
    p[[i]] <- sapply( 1:K , function(k) {
      if ( k < K ) {
        ptemp <- post_m2_all_mas[,k] 
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

p_mean.brms_mas <- sapply( 1:length(p) , function(i) apply(p[[i]],2,mean) )
p_HPDI.brms_mas <- sapply( 1:length(p) , function(i) apply(p[[i]],2,HPDI, prob=0.95) )
pred_probs_mod2_mas.brms<- cbind(p_mean.brms_mas, p_HPDI.brms_mas[c(1,3,5),], p_HPDI.brms_mas[c(2,4,6),])
rownames(pred_probs_mod2_mas.brms)<- c("P","S","G")
colnames(pred_probs_mod2_mas.brms)<- c("mean", "lwr95", "upr95")
round(pred_probs_mod2_mas.brms,2)
#mean lwr95 upr95
#P 0.24     0  1.00
#S 0.61     0  1.00
#G 0.16     0  0.97


