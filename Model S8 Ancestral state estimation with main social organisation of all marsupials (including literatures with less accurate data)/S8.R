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
setwd("C:/Model S8 Ancestral state estimation with main social organisation of all marsupials (including literatures with less accurate data)")

########################################################## PREPARE DATA###########################################################
d_mas<- as.data.frame(read_excel ("S8.xlsx"))

## prepare phylogeny
mas4<-read.nexus("outputS8.nex")

### remove species without any data(normally in this excel all sp has info)
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
d_mas2$Main_SO<- relevel(as.factor(d_mas2$Main_SO), ref="S") #S as reference so it's not in the following prior


#Center predictors
d_mas2$Habitat_hetero<-(d_mas2$Habitat_hetero)-1    #this model is to estimate ancestral state. I excluded habitat so we can have more species
d_mas2$Nbr_studies<-(d_mas2$Nbr_studies)-1
d_mas2$logBM.foss<- (log(d_mas2[,"body_mass(g)"])-log(37.75)/sd(log(d_mas2[,"body_mass(g)"]), na.rm=TRUE)) # 4.9 is the estimated body mass of marsupial ancestor (Djarthia murgonensis)



#model4
bform<- bf(Main_SO~1 + Nbr_studies + Habitat_hetero + logBM.foss +(1|phylo)+(1|Genus_species))

get_prior(bform, data=d_mas2, family=categorical)
prior=c(prior(normal(0,10), class=Intercept),
        prior(cauchy(0,2), class=sd, dpar="muP"), 
        prior(cauchy(0,2), class=sd, dpar="muG"),

        prior(normal(0,1), class=b))

m3_mas<- brm(bform, data=d_mas2, family=categorical, cov_ranef = list(phylo = A[[1]]), prior = prior, 
             chains = 2, cores = 2, iter = 2000, warmup = 1000, thin=5, control = list(adapt_delta = 0.99, max_treedepth = 15))
summary(m3_mas)

#save result(1 repeat)
save(m3_mas,file="m3_mas.RData")
load("m3_mas.RData")

post_m3_mas<- posterior_samples(m3_mas)

str(post_m3_mas[,c(1:12)])

post_m3_all_mas<- post_m3_mas[,c(1:12)]



#########################################repeat 100 times###############################################
for(i in 2:100){
  m3_i_mas<- brm(bform, data=d_mas2, family=categorical, cov_ranef = list(phylo = A[[i]]), prior = prior,
                 chains = 2, cores = 2, iter = 2000, warmup = 1000, thin=5, control = list(adapt_delta = 0.99, max_treedepth = 15))
  post_m3_mas<-posterior_samples(m3_i_mas)
  post_m3_all_mas<- rbind(post_m3_all_mas, post_m3_mas[,c(1:12)])
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
VarPhy.P<-post_m3_all_mas$sd_phylo__muP_Intercept #no references
VarPhy.G<-post_m3_all_mas$sd_phylo__muG_Intercept


VarSpec.P<-post_m3_all_mas$sd_Genus_species__muP_Intercept                            
VarSpec.G<-post_m3_all_mas$sd_Genus_species__muG_Intercept                             


VarDistro<- pi^2 / 3

lambda<- (VarPhy.P+VarPhy.G)/
  (VarPhy.P+VarPhy.G+VarSpec.P+VarSpec.G+VarDistro)
mean(lambda); median(lambda); HPDI(lambda, prob=0.95) 
#[1] 0.1452618
#[1] 0.1346708
#|0.95      0.95| 
#  0.04463299 0.26619006 

{    K <- 3 #number of social states used
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
pred_probs_mod5_mas.brms<- cbind(p_mean.brms_mas, p_HPDI.brms_mas[c(1,3,5),], p_HPDI.brms_mas[c(2,4,6),])
rownames(pred_probs_mod5_mas.brms)<- c("G","P","S")
colnames(pred_probs_mod5_mas.brms)<- c("mean", "lwr95", "upr95")
round(pred_probs_mod5_mas.brms,2)
#  mean lwr95 upr95
#G 0.12  0.00  0.82
#P 0.29  0.00  0.97
#S 0.59  0.01  1.00

