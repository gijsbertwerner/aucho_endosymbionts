# Match and Mismatch between Dietary Switches and Microbial Partners in Plant Sap-Feeding Insects
#Louis Bell-Roberts, Angela Douglas, Gijsbert Werner

#Script (C) Gijsbert Werner, University of Oxford, 2019

# Packages ----------------------------------------------------------------

#Load packages required for these analyses
library(ape)
library(phytools)
library(stringr)
library(diversitree)
library(RColorBrewer)
library(phangorn)
library(corHMM)
library(qpcR)
library(parallel)
library(ggplot2)

# Data --------------------------------------------------------------------

#Load the data-files

#Phylogenies read-in
aucho_subfam_dated<-read.tree("./Data/Dating/r8s/dist/aucho_till_subfamily_boosterweb_dated.phy")

#Quick filecheck
aucho_subfam_dated
is.binary(aucho_subfam_dated)
is.ultrametric(aucho_subfam_dated,tol = 0.0001) #Is actually ultrametric, rounding errors. 
plot(aucho_subfam_dated)

#Data file
dat_aucho<-read.csv("./Data/Master_db_standardised_unique_with_taxonomy.csv",as.is=T,strip.white = T,row.names = NULL)
head(dat_aucho)
sapply(dat_aucho, class)

# Main Data Formatting ----------------------------------------------------

#Create variables needed for downstream analysis
head(dat_aucho)

#Xylem binary
table(dat_aucho$diet,useNA = "ifany")
dat_aucho$diet_xylem<-ifelse(dat_aucho$diet=="Xylem",1,0)
table(dat_aucho$diet_xylem,useNA = "ifany")

#Any symbiont
table(dat_aucho$primary.endosymbiont,dat_aucho$companion.symbiont,useNA = "ifany")
dat_aucho$any_symbiont<-ifelse(dat_aucho$primary.endosymbiont==1|dat_aucho$companion.symbiont==1,1,0)
table(dat_aucho$any_symbiont,useNA = "ifany")
dat_aucho %>% dplyr::select(primary.endosymbiont,companion.symbiont,any_symbiont) #Yes, this has propagated correctly. Good

#Any beta-symbiont
table(dat_aucho$companion.symbiont.taxonomic.group,useNA = "ifany")
dat_aucho$any_beta_symbiont<-ifelse(dat_aucho$companion.symbiont.taxonomic.group %in% c("beta_proteobacteria","beta_proteobacteria_plus_gamma_proteobacteria"),1,0)
dat_aucho$any_beta_symbiont<-ifelse(is.na(dat_aucho$companion.symbiont.taxonomic.group),NA,dat_aucho$any_beta_symbiont)
table(dat_aucho$any_beta_symbiont,useNA = "ifany")

#Any gamma-symbiont
table(dat_aucho$companion.symbiont.taxonomic.group,useNA = "ifany")
dat_aucho$any_gamma_symbiont<-ifelse(dat_aucho$companion.symbiont.taxonomic.group %in% c("gamma_proteobacteria","beta_proteobacteria_plus_gamma_proteobacteria"),1,0)
dat_aucho$any_gamma_symbiont<-ifelse(is.na(dat_aucho$companion.symbiont.taxonomic.group),NA,dat_aucho$any_gamma_symbiont)
table(dat_aucho$any_gamma_symbiont,useNA = "ifany")

#Any YSL
table(dat_aucho$companion.symbiont.taxonomic.group,useNA = "ifany")
dat_aucho$any_YLS_symbiont<-ifelse(dat_aucho$companion.symbiont.taxonomic.group %in% c("YLS"),1,0)
dat_aucho$any_YLS_symbiont<-ifelse(is.na(dat_aucho$companion.symbiont.taxonomic.group),NA,dat_aucho$any_YLS_symbiont)
table(dat_aucho$any_YLS_symbiont,useNA = "ifany")

#Remove the genera for which we have no Sulcia data. #See if we keep this. 
table(dat_aucho$primary.endosymbiont,useNA = "ifany")
dat_aucho <-dat_aucho %>% dplyr::filter(!is.na(primary.endosymbiont))
table(dat_aucho$primary.endosymbiont,useNA = "ifany")

#Subset tree and dataset to their overlap. 
small_aucho_subfam_dated<-drop.tip(phy = aucho_subfam_dated,tip = setdiff(aucho_subfam_dated$tip.label,dat_aucho$genus))
small_aucho_subfam_dated
#Same tip label ordering (this is important later on for tip label colours)
small_aucho_subfam_dated$tip.label==small_aucho_subfam_dated$tip.label 
#Reduce dataset
reduced_dat_aucho<-dat_aucho %>% dplyr::filter(genus %in% small_aucho_subfam_dated$tip.label)
reduced_dat_aucho_small_aucho_subfam_order<-reduced_dat_aucho[match(small_aucho_subfam_dated$tip.label,reduced_dat_aucho$genus),]

#Print the analysed genus list #Potentially still reverse to achieve match with figure order. 
write.csv(rev(small_aucho_subfam_dated$tip.label),"./Output/Genera_Analysed.csv")

# Data Descriptives -------------------------------------------------------
head(reduced_dat_aucho_small_aucho_subfam_order)

#Let's have a look at our variables of interest

#The primary symbionts
table(reduced_dat_aucho_small_aucho_subfam_order$primary.endosymbiont,useNA = "ifany")
table(reduced_dat_aucho_small_aucho_subfam_order$primary.endosymbiont.identity,useNA = "ifany")
table(reduced_dat_aucho_small_aucho_subfam_order$primary.endosymbiont.taxonomic.group,useNA = "ifany")
#The primary endosymbiont, both as a binary (presence/absence) variable, all. are Bacteroidetes (with the exeption of the outgroup)

#The secondary symbionts
table(reduced_dat_aucho_small_aucho_subfam_order$companion.symbiont,useNA = "ifany")
table(reduced_dat_aucho_small_aucho_subfam_order$companion.symbiont.identity,useNA = "ifany")
table(reduced_dat_aucho_small_aucho_subfam_order$companion.symbiont.taxonomic.group,useNA = "ifany")
#Presence absence again, identity (too much detail). The taxonomic group level-idenity is less detailed than Louis' report Figure 2, though. Talk to him. 

#Diet
table(reduced_dat_aucho_small_aucho_subfam_order$diet,useNA = "ifany")


# Data Analysis - ASRs -------------------------------------------------------

set.seed(01865)

######Diet
head(reduced_dat_aucho_small_aucho_subfam_order)
table(reduced_dat_aucho_small_aucho_subfam_order$diet)

#Data formatting. We need two columns, tip label and trait. We'll set the outgroup as NA, since we don't want to introduce two more states
dat_ASR_diet<-reduced_dat_aucho_small_aucho_subfam_order %>% dplyr::select(genus,diet)
dat_ASR_diet$diet[dat_ASR_diet$diet %in% c("Moss","Predatory")] <-"Parenchyma&Phloem&Xylem"
table(dat_ASR_diet$diet,useNA = "ifany")

#Run analyses
ASR_diet_ER_yang_subfam_dated<-rayDISC(phy = small_aucho_subfam_dated,data = dat_ASR_diet,ntraits = 1,
                                       model="ER",node.states = "marginal",root.p="yang",
                                       verbose = T,lewis.asc.bias = T)
ASR_diet_SYM_yang_subfam_dated<-rayDISC(phy = small_aucho_subfam_dated,data = dat_ASR_diet,ntraits = 1,
                                        model="SYM",node.states = "marginal",root.p="yang",
                                        verbose = T,lewis.asc.bias = T)
ASR_diet_ARD_yang_subfam_dated<-rayDISC(phy = small_aucho_subfam_dated,data = dat_ASR_diet,ntraits = 1,
                                        model="ARD",node.states = "marginal",root.p="yang",
                                        verbose = T,lewis.asc.bias = T)

#Which is the best model? 
akaike.weights(c(ASR_diet_ER_yang_subfam_dated$AICc,
                 ASR_diet_SYM_yang_subfam_dated$AICc,
                 ASR_diet_ARD_yang_subfam_dated$AICc))

#ARD by far the best
ASR_diet_ARD_yang_subfam_dated

#Save all model ran
save(ASR_diet_ER_yang_subfam_dated,file="./Output/ASR_diet_ER_yang_subfam_dated")
save(ASR_diet_SYM_yang_subfam_dated,file="./Output/ASR_diet_SYM_yang_subfam_dated")
save(ASR_diet_ARD_yang_subfam_dated,file="./Output/ASR_diet_ARD_yang_subfam_dated")

#What is the ancestral state at the aucho root? 
getMRCA(ASR_diet_ARD_yang_subfam_dated$phy,tip = c("Ledra","Cedusa"))
plot.phylo(ASR_diet_ARD_yang_subfam_dated$phy,show.tip.label = F)
nodelabels(cex=0.5)
round(head(ASR_diet_ARD_yang_subfam_dated$states,12),3)
head(ASR_diet_ARD_yang_subfam_dated$phy$edge,12)

pdf("./Output/FigS2_HighDefinition.pdf",width = 12,height = 12)
plotRECON(ASR_diet_ARD_yang_subfam_dated$phy,ASR_diet_ARD_yang_subfam_dated$states,
          piecolors = c("#fdae61","#d9ef8b","#1a9641"),
          label.offset=3.85,cex=0.45,pie.cex=0.325)
tiplabels(pch=23,offset=1.5,cex=1.25,
          bg= c("#fdae61","black","#d9ef8b","#1a9641")[as.factor(dat_ASR_diet$diet)])
polygon(x=c(0,50,50,0),y=c(0,0,145,145),col="#80808025",border="#80808025")
polygon(x=c(100,150,150,100),y=c(0,0,145,145),col="#80808025",border="#80808025")
polygon(x=c(200,250,250,200),y=c(0,0,145,145),col="#80808025",border="#80808025")
par(new=T)
plotRECON(ASR_diet_ARD_yang_subfam_dated$phy,ASR_diet_ARD_yang_subfam_dated$states,
          piecolors = c("#fdae61","#d9ef8b","#1a9641"),
          label.offset=3.85,cex=0.45,pie.cex=0.325)
tiplabels(pch=23,offset=1.5,cex=1.25,
          bg= c("#fdae61","black","#d9ef8b","#1a9641")[as.factor(dat_ASR_diet$diet)])
dev.off()

######Primary Symbiosis
head(reduced_dat_aucho_small_aucho_subfam_order)
table(reduced_dat_aucho_small_aucho_subfam_order$primary.endosymbiont,useNA = "ifany")

#Data formatting. We need two columns, tip label and trait. 
#We'll set the outgroup as NA with respect to Sulcia/primary status, since we don't want to introduce additional states to account for their symbiont status. 
dat_ASR_primary.endosymbiont<-reduced_dat_aucho_small_aucho_subfam_order %>% dplyr::select(genus,primary.endosymbiont)
dat_ASR_primary.endosymbiont$primary.endosymbiont[is.na(dat_ASR_primary.endosymbiont$primary.endosymbiont)] <-"0&1"
dat_ASR_primary.endosymbiont_corhmm<-dat_ASR_primary.endosymbiont %>% dplyr::filter(primary.endosymbiont!="0&1") #Remove the NAs (necessary for corHMM analysis to work)
table(dat_ASR_primary.endosymbiont_corhmm$primary.endosymbiont,useNA = "ifany")
nrow(dat_ASR_primary.endosymbiont)

corenum<-detectCores() #This determines the number of cores avaialbe on your machine to run the subsequent (very time-consuming) analyses

corHMM_primary.endosymbiont_rate1_yang_subfam_dated<-corHMM(phy = drop.tip(small_aucho_subfam_dated,
                                                                           tip = small_aucho_subfam_dated$tip.label[!small_aucho_subfam_dated$tip.label %in% dat_ASR_primary.endosymbiont_corhmm$genus]),
                                                            data = dat_ASR_primary.endosymbiont_corhmm,
                                                            rate.cat = 1,
                                                            n.core = corenum - 1,nstarts=100,optim.method = "subplex",
                                                            node.states = "marginal",root.p=c(0,1),
                                                            diagn = T)
corHMM_primary.endosymbiont_rate1_yang_subfam_dated
corHMM_primary.endosymbiont_rate2_yang_subfam_dated<-corHMM(phy = drop.tip(small_aucho_subfam_dated,
                                                                           tip = small_aucho_subfam_dated$tip.label[!small_aucho_subfam_dated$tip.label %in% dat_ASR_primary.endosymbiont_corhmm$genus]),
                                                            data = dat_ASR_primary.endosymbiont_corhmm,
                                                            rate.cat = 2,
                                                            n.core = corenum - 1,nstarts=100,optim.method = "subplex",
                                                            node.states = "marginal",root.p=c(0,1),
                                                            diagn = T)
corHMM_primary.endosymbiont_rate2_yang_subfam_dated
corHMM_primary.endosymbiont_rate3_yang_subfam_dated<-corHMM(phy = drop.tip(small_aucho_subfam_dated,
                                                                           tip = small_aucho_subfam_dated$tip.label[!small_aucho_subfam_dated$tip.label %in% dat_ASR_primary.endosymbiont_corhmm$genus]),
                                                            data = dat_ASR_primary.endosymbiont_corhmm,
                                                            rate.cat = 3,
                                                            n.core = corenum - 1,nstarts=100,optim.method = "subplex",
                                                            node.states = "marginal",root.p=c(0,1),
                                                            diagn = T)
corHMM_primary.endosymbiont_rate3_yang_subfam_dated
corHMM_primary.endosymbiont_rate4_yang_subfam_dated<-corHMM(phy = drop.tip(small_aucho_subfam_dated,
                                                                           tip = small_aucho_subfam_dated$tip.label[!small_aucho_subfam_dated$tip.label %in% dat_ASR_primary.endosymbiont_corhmm$genus]),
                                                            data = dat_ASR_primary.endosymbiont_corhmm,
                                                            rate.cat = 4,
                                                            n.core = corenum - 1,nstarts=100,optim.method = "subplex",
                                                            node.states = "marginal",root.p=c(0,1),
                                                            diagn = T)
corHMM_primary.endosymbiont_rate4_yang_subfam_dated
corHMM_primary.endosymbiont_rate5_yang_subfam_dated<-corHMM(phy = drop.tip(small_aucho_subfam_dated,
                                                                           tip = small_aucho_subfam_dated$tip.label[!small_aucho_subfam_dated$tip.label %in% dat_ASR_primary.endosymbiont_corhmm$genus]),
                                                            data = dat_ASR_primary.endosymbiont_corhmm,
                                                            rate.cat = 5,
                                                            n.core = corenum - 1,nstarts=100,optim.method = "subplex",
                                                            node.states = "marginal",root.p=c(0,1),
                                                            diagn = T)
corHMM_primary.endosymbiont_rate5_yang_subfam_dated

#Which is best? 
akaike.weights(c(corHMM_primary.endosymbiont_rate1_yang_subfam_dated$AICc,
                 corHMM_primary.endosymbiont_rate2_yang_subfam_dated$AICc,
                 corHMM_primary.endosymbiont_rate3_yang_subfam_dated$AICc,
                 corHMM_primary.endosymbiont_rate4_yang_subfam_dated$AICc,
                 corHMM_primary.endosymbiont_rate5_yang_subfam_dated$AICc))

#By far the two rate class model is best
corHMM_primary.endosymbiont_rate2_yang_subfam_dated
round(corHMM_primary.endosymbiont_rate2_yang_subfam_dated$solution*100,2)
#Relative Gain/Loss rates
corHMM_primary.endosymbiont_rate2_yang_subfam_dated$solution[2,1]/corHMM_primary.endosymbiont_rate2_yang_subfam_dated$solution[1,2]
corHMM_primary.endosymbiont_rate2_yang_subfam_dated$solution[4,3]/corHMM_primary.endosymbiont_rate2_yang_subfam_dated$solution[3,4]

#Save all models to output
save(corHMM_primary.endosymbiont_rate1_yang_subfam_dated,file="./Output/corHMM_primary.endosymbiont_rate1_yang_subfam_dated")
save(corHMM_primary.endosymbiont_rate2_yang_subfam_dated,file="./Output/corHMM_primary.endosymbiont_rate2_yang_subfam_dated")
save(corHMM_primary.endosymbiont_rate3_yang_subfam_dated,file="./Output/corHMM_primary.endosymbiont_rate3_yang_subfam_dated")
save(corHMM_primary.endosymbiont_rate4_yang_subfam_dated,file="./Output/corHMM_primary.endosymbiont_rate4_yang_subfam_dated")
save(corHMM_primary.endosymbiont_rate5_yang_subfam_dated,file="./Output/corHMM_primary.endosymbiont_rate5_yang_subfam_dated")

#For plotting
dat_ASR_primary.endosymbiont_corhmm<-
  dat_ASR_primary.endosymbiont_corhmm[match(corHMM_primary.endosymbiont_rate2_yang_subfam_dated$phy$tip.label,dat_ASR_primary.endosymbiont_corhmm$genus),]
dat_ASR_primary.endosymbiont_corhmm$genus== corHMM_primary.endosymbiont_rate2_yang_subfam_dated$phy$tip.label
nrow(dat_ASR_primary.endosymbiont_corhmm)

pdf("./Output/FigS4_HighDefinition.pdf",width = 20,height = 20)
plotRECON(corHMM_primary.endosymbiont_rate2_yang_subfam_dated$phy,corHMM_primary.endosymbiont_rate2_yang_subfam_dated$states,
          piecolors = c("#33a02c","#1f78b4","#b2df8a","#a6cee3"),
          label.offset=3.85,cex=0.55,pie.cex=0.325)
tiplabels(pch=23,offset=1.5,cex=1.25,
          bg=brewer.pal(n=3,"Greys")[c(2,3)][as.factor(dat_ASR_primary.endosymbiont_corhmm$primary.endosymbiont)])
polygon(x=c(0,50,50,0),y=c(0,0,145,145),col="#80808025",border="#80808025")
polygon(x=c(100,150,150,100),y=c(0,0,145,145),col="#80808025",border="#80808025")
polygon(x=c(200,250,250,200),y=c(0,0,145,145),col="#80808025",border="#80808025")
par(new=T)
plotRECON(corHMM_primary.endosymbiont_rate2_yang_subfam_dated$phy,corHMM_primary.endosymbiont_rate2_yang_subfam_dated$states,
          piecolors = c("#33a02c","#1f78b4","#b2df8a","#a6cee3"),
          label.offset=3.85,cex=0.55,pie.cex=0.325)
tiplabels(pch=23,offset=1.5,cex=1.25,
          bg=brewer.pal(n=3,"Greys")[c(2,3)][as.factor(dat_ASR_primary.endosymbiont_corhmm$primary.endosymbiont)])
dev.off()

####Calculate the number of loss events

#Extract the states per node
states_recon_binarised<-data.frame(
  abs_perc=corHMM_primary.endosymbiont_rate2_yang_subfam_dated$states[,1]+corHMM_primary.endosymbiont_rate2_yang_subfam_dated$states[,3],
  pres_perc=corHMM_primary.endosymbiont_rate2_yang_subfam_dated$states[,2]+corHMM_primary.endosymbiont_rate2_yang_subfam_dated$states[,4]
)
#Assign a vector to the internal nodes, with the right edge numbers
states_recon_binarised$vec<-
  (143):(143+nrow(states_recon_binarised)-1)
#Extract the edges
edges<-as.data.frame(corHMM_primary.endosymbiont_rate2_yang_subfam_dated$phy$edge)
names(edges)<-c("parent","daughter")
#Merge the two, mathching 'parent nodes' to the node vector
states_binarised_parent_daughter<-merge(x=edges,y = states_recon_binarised,
                                        by.x="parent",by.y="vec")
#And then the 'daughter nodes' 
states_binarised_parent_daughter<-merge(x=states_binarised_parent_daughter,y=states_recon_binarised,
                                        by.x="daughter",by.y = "vec",all.x = T)
states_binarised_parent_daughter
names(states_binarised_parent_daughter)<-c("parent","daughter",
                                           "par_abs_perc","par_pres_perc","dau_abs_perc","dau_pres_perc")

##Now calculate the losses, per node
states_binarised_parent_daughter$primary_loss<-
  states_binarised_parent_daughter$par_pres_perc-states_binarised_parent_daughter$dau_pres_perc
states_binarised_parent_daughter$primary_gain<-
  states_binarised_parent_daughter$par_abs_perc-states_binarised_parent_daughter$dau_abs_perc
#Sum ove the tree
sum_primary_loss<-sum(states_binarised_parent_daughter$primary_loss[states_binarised_parent_daughter$primary_loss>0],
                      na.rm=T)
sum_primary_gain<-sum(states_binarised_parent_daughter$primary_gain[states_binarised_parent_daughter$primary_gain>0],
                      na.rm=T)
#Print them
sum_primary_loss
sum_primary_gain

######Secondary Symbiosis
#Categorical
head(reduced_dat_aucho_small_aucho_subfam_order)
table(reduced_dat_aucho_small_aucho_subfam_order$companion.symbiont.taxonomic.group,useNA = "ifany")

#Data formatting. We need two columns, tip label and trait. 
dat_ASR_companion.symbiont.taxonomic.group<-reduced_dat_aucho_small_aucho_subfam_order %>% dplyr::select(genus,companion.symbiont.taxonomic.group)
#The follwoing evaluates all the NAs as potentially having each of the potential states.
dat_ASR_companion.symbiont.taxonomic.group$companion.symbiont.taxonomic.group[is.na(dat_ASR_companion.symbiont.taxonomic.group$companion.symbiont.taxonomic.group)]<-"Absent&alfa_proteobacteria&beta_proteobacteria&gamma_proteobacteria&YLS"
#A few cases are coded both as beta and gamma, analyse them accordinly. 
dat_ASR_companion.symbiont.taxonomic.group$companion.symbiont.taxonomic.group[dat_ASR_companion.symbiont.taxonomic.group$companion.symbiont.taxonomic.group=="beta_proteobacteria_plus_gamma_proteobacteria"] <-"beta_proteobacteria&gamma_proteobacteria"

table(dat_ASR_companion.symbiont.taxonomic.group$companion.symbiont.taxonomic.group,useNA = "ifany")
#Looks good 

#Run analyses
ASR_companion.symbiont.taxonomic.group_ER_yang_subfam_dated<-rayDISC(phy = small_aucho_subfam_dated,data = dat_ASR_companion.symbiont.taxonomic.group,ntraits = 1,
                                                                     model="ER",node.states = "marginal",root.p="yang",
                                                                     verbose = T,lewis.asc.bias = T)
ASR_companion.symbiont.taxonomic.group_SYM_yang_subfam_dated<-rayDISC(phy = small_aucho_subfam_dated,data = dat_ASR_companion.symbiont.taxonomic.group,ntraits = 1,
                                                                      model="SYM",node.states = "marginal",root.p="yang",
                                                                      verbose = T,lewis.asc.bias = T)
ASR_companion.symbiont.taxonomic.group_ARD_yang_subfam_dated<-rayDISC(phy = small_aucho_subfam_dated,data = dat_ASR_companion.symbiont.taxonomic.group,ntraits = 1,
                                                                      model="ARD",node.states = "marginal",root.p="yang",
                                                                      verbose = T,lewis.asc.bias = T)

#Which is the best model? 
akaike.weights(c(ASR_companion.symbiont.taxonomic.group_ER_yang_subfam_dated$AICc,
                 ASR_companion.symbiont.taxonomic.group_SYM_yang_subfam_dated$AICc,
                 ASR_companion.symbiont.taxonomic.group_ARD_yang_subfam_dated$AICc))

#SYM by far the best
ASR_companion.symbiont.taxonomic.group_SYM_yang_subfam_dated

#Save them all
save(ASR_companion.symbiont.taxonomic.group_ER_yang_subfam_dated,file="./Output/ASR_companion.symbiont.taxonomic.group_ER_yang_subfam_dated")
save(ASR_companion.symbiont.taxonomic.group_SYM_yang_subfam_dated,file="./Output/ASR_companion.symbiont.taxonomic.group_SYM_yang_subfam_dated")
save(ASR_companion.symbiont.taxonomic.group_ARD_yang_subfam_dated,file="./Output/ASR_companion.symbiont.taxonomic.group_ARD_yang_subfam_dated")

#What is the auco ancestral state?
getMRCA(ASR_companion.symbiont.taxonomic.group_SYM_yang_subfam_dated$phy,tip = c("Ledra","Cedusa"))
head(ASR_companion.symbiont.taxonomic.group_SYM_yang_subfam_dated$phy$edge,12)
round(head(ASR_companion.symbiont.taxonomic.group_SYM_yang_subfam_dated$states,12),3)

pdf("./Output/FigS5_HighDefinition.pdf",width = 12,height = 12)
plotRECON(ASR_companion.symbiont.taxonomic.group_SYM_yang_subfam_dated$phy,
          ASR_companion.symbiont.taxonomic.group_SYM_yang_subfam_dated$states,
          piecolors = c(brewer.pal(n=7,"Set2")[1:3],brewer.pal(n=7,"Set2")[5],brewer.pal(n=7,"Set2")[7]),
          label.offset=3.85,cex=0.45,pie.cex=0.325)
tiplabels(pch=23,offset=1.5,cex=1.25,
          bg= c(brewer.pal(n=7,"Set2")[1:5],brewer.pal(n=7,"Set2")[7])[as.factor(reduced_dat_aucho_small_aucho_subfam_order$companion.symbiont.taxonomic.group)])
polygon(x=c(0,50,50,0),y=c(0,0,145,145),col="#80808025",border="#80808025")
polygon(x=c(100,150,150,100),y=c(0,0,145,145),col="#80808025",border="#80808025")
polygon(x=c(200,250,250,200),y=c(0,0,145,145),col="#80808025",border="#80808025")
par(new=T)
plotRECON(ASR_companion.symbiont.taxonomic.group_SYM_yang_subfam_dated$phy,
          ASR_companion.symbiont.taxonomic.group_SYM_yang_subfam_dated$states,
          piecolors = c(brewer.pal(n=7,"Set2")[1:3],brewer.pal(n=7,"Set2")[5],brewer.pal(n=7,"Set2")[7]),
          label.offset=3.85,cex=0.45,pie.cex=0.325)
tiplabels(pch=23,offset=1.5,cex=1.25,
          bg= c(brewer.pal(n=7,"Set2")[1:5],brewer.pal(n=7,"Set2")[7])[as.factor(reduced_dat_aucho_small_aucho_subfam_order$companion.symbiont.taxonomic.group)])
dev.off()

# Data Analysis - Correlated Models of Evolution  -------------------------------------------------------

###Any primary symbiont vs diet

#Data prep
head(reduced_dat_aucho_small_aucho_subfam_order)
dat_fitPagel_diet_xylem_primary.endosymbiont<-reduced_dat_aucho_small_aucho_subfam_order %>% 
  dplyr::filter(!(is.na(diet_xylem)|is.na(primary.endosymbiont))) %>%
  dplyr::select(genus,primary.endosymbiont,diet_xylem)
dat_fitPagel_diet_xylem_primary.endosymbiont$diet_xylem_primary.endosymbiont_joined<-
  paste0(dat_fitPagel_diet_xylem_primary.endosymbiont$primary.endosymbiont,
         dat_fitPagel_diet_xylem_primary.endosymbiont$diet_xylem)
nrow(dat_fitPagel_diet_xylem_primary.endosymbiont)
head(dat_fitPagel_diet_xylem_primary.endosymbiont)

vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont<-
  dat_fitPagel_diet_xylem_primary.endosymbiont$primary.endosymbiont
vec_fitPagel_diet_xylem_primary.endosymbiont_diet_xylem<-
  dat_fitPagel_diet_xylem_primary.endosymbiont$diet_xylem
names(vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont)<-dat_fitPagel_diet_xylem_primary.endosymbiont$genus
names(vec_fitPagel_diet_xylem_primary.endosymbiont_diet_xylem)<-dat_fitPagel_diet_xylem_primary.endosymbiont$genus
table(vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont)
table(vec_fitPagel_diet_xylem_primary.endosymbiont_diet_xylem)
head(vec_fitPagel_diet_xylem_primary.endosymbiont_diet_xylem,25)
head(vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont,25)

#Run
corevol_diet_xylem_primary.endosymbiont_dated_subfam<-
  fitPagel(tree = drop.tip(small_aucho_subfam_dated,
                           tip = small_aucho_subfam_dated$tip.label[!small_aucho_subfam_dated$tip.label %in% dat_fitPagel_diet_xylem_primary.endosymbiont$genus]),
           x=vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont,
           y=vec_fitPagel_diet_xylem_primary.endosymbiont_diet_xylem)
corevol_diet_xylem_primary.endosymbiont_dated_subfam
save(corevol_diet_xylem_primary.endosymbiont_dated_subfam,file="./Output/corevol_diet_xylem_primary.endosymbiont_dated_subfam")

#Reconstruct
bestrates_corevol_diet_xylem_primary.endosymbiont_dated_subfam<-
  as.vector(corevol_diet_xylem_primary.endosymbiont_dated_subfam$dependent.Q)[c(2,3,5,8,9,12,14,15)]

recon_corevol_diet_xylem_primary.endosymbiont_dated_subfam<-
  ancRECON(phy = drop.tip(small_aucho_subfam_dated,
                          tip = small_aucho_subfam_dated$tip.label[!small_aucho_subfam_dated$tip.label %in% dat_fitPagel_diet_xylem_primary.endosymbiont$genus]),
           p=bestrates_corevol_diet_xylem_primary.endosymbiont_dated_subfam,
           data = dat_fitPagel_diet_xylem_primary.endosymbiont,
           method = "marginal",ntraits = 2,model = "ARD",
           hrm = F,rate.cat = NULL,
           root.p = "yang")

#Plot
round(corevol_diet_xylem_primary.endosymbiont_dated_subfam$dependent.Q*100,2)

###Any beta symbiont vs diet

#Data prep
head(reduced_dat_aucho_small_aucho_subfam_order)
dat_fitPagel_diet_xylem_any_beta_symbiont<-reduced_dat_aucho_small_aucho_subfam_order %>% 
  dplyr::filter(!(is.na(diet_xylem)|is.na(any_beta_symbiont))) %>%
  dplyr::select(genus,any_beta_symbiont,diet_xylem)
dat_fitPagel_diet_xylem_any_beta_symbiont$diet_xylem_any_beta_symbiont_joined<-
  paste0(dat_fitPagel_diet_xylem_any_beta_symbiont$any_beta_symbiont,
         dat_fitPagel_diet_xylem_any_beta_symbiont$diet_xylem)
nrow(dat_fitPagel_diet_xylem_any_beta_symbiont)
head(dat_fitPagel_diet_xylem_any_beta_symbiont)

vec_fitPagel_diet_xylem_any_beta_symbiont_any_beta_symbiont<-
  dat_fitPagel_diet_xylem_any_beta_symbiont$any_beta_symbiont
vec_fitPagel_diet_xylem_any_beta_symbiont_diet_xylem<-
  dat_fitPagel_diet_xylem_any_beta_symbiont$diet_xylem
names(vec_fitPagel_diet_xylem_any_beta_symbiont_any_beta_symbiont)<-dat_fitPagel_diet_xylem_any_beta_symbiont$genus
names(vec_fitPagel_diet_xylem_any_beta_symbiont_diet_xylem)<-dat_fitPagel_diet_xylem_any_beta_symbiont$genus
table(vec_fitPagel_diet_xylem_any_beta_symbiont_any_beta_symbiont)
table(vec_fitPagel_diet_xylem_any_beta_symbiont_diet_xylem)
head(vec_fitPagel_diet_xylem_any_beta_symbiont_diet_xylem)

#Run
corevol_diet_xylem_any_beta_symbiont_dated_subfam<-
  fitPagel(tree = drop.tip(small_aucho_subfam_dated,
                           tip = small_aucho_subfam_dated$tip.label[!small_aucho_subfam_dated$tip.label %in% dat_fitPagel_diet_xylem_any_beta_symbiont$genus]),
           x=vec_fitPagel_diet_xylem_any_beta_symbiont_any_beta_symbiont,
           y=vec_fitPagel_diet_xylem_any_beta_symbiont_diet_xylem)
corevol_diet_xylem_any_beta_symbiont_dated_subfam
save(corevol_diet_xylem_any_beta_symbiont_dated_subfam,file="./Output/corevol_diet_xylem_any_beta_symbiont_dated_subfam")


###Any gamma symbiont vs diet

#Data prep
head(reduced_dat_aucho_small_aucho_subfam_order)
dat_fitPagel_diet_xylem_any_gamma_symbiont<-reduced_dat_aucho_small_aucho_subfam_order %>% 
  dplyr::filter(!(is.na(diet_xylem)|is.na(any_gamma_symbiont))) %>%
  dplyr::select(genus,any_gamma_symbiont,diet_xylem)
dat_fitPagel_diet_xylem_any_gamma_symbiont$diet_xylem_any_gamma_symbiont_joined<-
  paste0(dat_fitPagel_diet_xylem_any_gamma_symbiont$any_gamma_symbiont,
         dat_fitPagel_diet_xylem_any_gamma_symbiont$diet_xylem)
nrow(dat_fitPagel_diet_xylem_any_gamma_symbiont)
head(dat_fitPagel_diet_xylem_any_gamma_symbiont)

vec_fitPagel_diet_xylem_any_gamma_symbiont_any_gamma_symbiont<-
  dat_fitPagel_diet_xylem_any_gamma_symbiont$any_gamma_symbiont
vec_fitPagel_diet_xylem_any_gamma_symbiont_diet_xylem<-
  dat_fitPagel_diet_xylem_any_gamma_symbiont$diet_xylem
names(vec_fitPagel_diet_xylem_any_gamma_symbiont_any_gamma_symbiont)<-dat_fitPagel_diet_xylem_any_gamma_symbiont$genus
names(vec_fitPagel_diet_xylem_any_gamma_symbiont_diet_xylem)<-dat_fitPagel_diet_xylem_any_gamma_symbiont$genus
table(vec_fitPagel_diet_xylem_any_gamma_symbiont_any_gamma_symbiont)
table(vec_fitPagel_diet_xylem_any_gamma_symbiont_diet_xylem)
head(vec_fitPagel_diet_xylem_any_gamma_symbiont_diet_xylem)

#Run
corevol_diet_xylem_any_gamma_symbiont_dated_subfam<-
  fitPagel(tree = drop.tip(small_aucho_subfam_dated,
                           tip = small_aucho_subfam_dated$tip.label[!small_aucho_subfam_dated$tip.label %in% dat_fitPagel_diet_xylem_any_gamma_symbiont$genus]),
           x=vec_fitPagel_diet_xylem_any_gamma_symbiont_any_gamma_symbiont,
           y=vec_fitPagel_diet_xylem_any_gamma_symbiont_diet_xylem)
corevol_diet_xylem_any_gamma_symbiont_dated_subfam
save(corevol_diet_xylem_any_gamma_symbiont_dated_subfam,file="./Output/corevol_diet_xylem_any_gamma_symbiont_dated_subfam")

###Any YLS symbiont vs diet

#Data prep
head(reduced_dat_aucho_small_aucho_subfam_order)
dat_fitPagel_diet_xylem_any_YLS_symbiont<-reduced_dat_aucho_small_aucho_subfam_order %>% 
  dplyr::filter(!(is.na(diet_xylem)|is.na(any_YLS_symbiont))) %>%
  dplyr::select(genus,any_YLS_symbiont,diet_xylem)
dat_fitPagel_diet_xylem_any_YLS_symbiont$diet_xylem_any_YLS_symbiont_joined<-
  paste0(dat_fitPagel_diet_xylem_any_YLS_symbiont$any_YLS_symbiont,
         dat_fitPagel_diet_xylem_any_YLS_symbiont$diet_xylem)
nrow(dat_fitPagel_diet_xylem_any_YLS_symbiont)
head(dat_fitPagel_diet_xylem_any_YLS_symbiont)

vec_fitPagel_diet_xylem_any_YLS_symbiont_any_YLS_symbiont<-
  dat_fitPagel_diet_xylem_any_YLS_symbiont$any_YLS_symbiont
vec_fitPagel_diet_xylem_any_YLS_symbiont_diet_xylem<-
  dat_fitPagel_diet_xylem_any_YLS_symbiont$diet_xylem
names(vec_fitPagel_diet_xylem_any_YLS_symbiont_any_YLS_symbiont)<-dat_fitPagel_diet_xylem_any_YLS_symbiont$genus
names(vec_fitPagel_diet_xylem_any_YLS_symbiont_diet_xylem)<-dat_fitPagel_diet_xylem_any_YLS_symbiont$genus
table(vec_fitPagel_diet_xylem_any_YLS_symbiont_any_YLS_symbiont)
table(vec_fitPagel_diet_xylem_any_YLS_symbiont_diet_xylem)
head(vec_fitPagel_diet_xylem_any_YLS_symbiont_diet_xylem)

#Run
corevol_diet_xylem_any_YLS_symbiont_dated_subfam<-
  fitPagel(tree = drop.tip(small_aucho_subfam_dated,
                           tip = small_aucho_subfam_dated$tip.label[!small_aucho_subfam_dated$tip.label %in% dat_fitPagel_diet_xylem_any_YLS_symbiont$genus]),
           x=vec_fitPagel_diet_xylem_any_YLS_symbiont_any_YLS_symbiont,
           y=vec_fitPagel_diet_xylem_any_YLS_symbiont_diet_xylem)
corevol_diet_xylem_any_YLS_symbiont_dated_subfam
save(corevol_diet_xylem_any_YLS_symbiont_dated_subfam,file="./Output/corevol_diet_xylem_any_YLS_symbiont_dated_subfam")



# Main text Figure --------------------------------------------------------

#Create a data frame to plot the trait data
plot_dat_bin<-reduced_dat_aucho_small_aucho_subfam_order %>% dplyr::select(primary.endosymbiont,diet,companion.symbiont.taxonomic.group)
nrow(plot_dat_bin)
row.names(plot_dat_bin)<-reduced_dat_aucho_small_aucho_subfam_order$genus
plot_dat_bin$diet<-as.numeric(as.factor(plot_dat_bin$diet))-1
plot_dat_bin$companion.symbiont.taxonomic.group<-as.numeric(as.factor(plot_dat_bin$companion.symbiont.taxonomic.group))-1
table(reduced_dat_aucho_small_aucho_subfam_order$companion.symbiont.taxonomic.group)
table(plot_dat_bin$companion.symbiont.taxonomic.group)
table(plot_dat_bin$diet)
head(plot_dat_bin)

#Create a branch colour vector based on companion symbionts status
head(ASR_companion.symbiont.taxonomic.group_SYM_yang_subfam_dated$phy$edge)
tail(ASR_companion.symbiont.taxonomic.group_SYM_yang_subfam_dated$phy$edge)
max(ASR_companion.symbiont.taxonomic.group_SYM_yang_subfam_dated$phy$edge)

#Get the ancetral node state from the diet reconstruction
companionnodelabs<-ASR_companion.symbiont.taxonomic.group_SYM_yang_subfam_dated$phy$node.label
companionnodelabs
names(companionnodelabs)<-143:276
companionnodelabs

colnodevec_companion<-companionnodelabs[
  match(ASR_companion.symbiont.taxonomic.group_SYM_yang_subfam_dated$phy$edge[,1],names(companionnodelabs))]
table(colnodevec_companion)
colnodevec_companion<-brewer.pal(n=7,"Set2")[c(1,3,5,7)][as.factor(colnodevec_companion)]
colnodevec_companion

#Main text superfigure
pdf("./Output/MainTextFigure.pdf",width = 20,height = 20)
trait.plot(tree = small_aucho_subfam_dated,dat = plot_dat_bin,
           cols = list(primary.endosymbiont=brewer.pal(n=3,"Greys")[c(2,3)],
                       diet = c("white","#fdae61","#d9ef8b","#1a9641"),
                       companion.symbiont.taxonomic.group=c(brewer.pal(n=7,"Set2")[1:5],brewer.pal(n=7,"Set2")[7])),
           type="p",legend=F,w=1/40,edge.width =2,
           cex.lab = 0.01,tip.color="white",
           edge.color=colnodevec_companion)
nodelabels(pie = recon_corevol_diet_xylem_primary.endosymbiont_dated_subfam$lik.anc.states,
           piecol = c("#a65628","#377eb8","#984ea3","#ff7f00"),cex=0.3)
add.scale.bar()
dev.off()


############################################################################################################################
######Sensitivity Analyses  -------------------------------------------------------
############################################################################################################################

#####1. Phylogenetic uncertainty -------------------------------------------------

#########
#1.1 Undated subfamily tree boostrap
aucho_subfam_undated_boots<-read.tree("./Data/RAxML_bootstrap.out_in_lbrem_outdefined_constraint_till_subfamily_minus_Formo_Xyphon_k")
aucho_subfam_undated_boots

#Subset to the genera where we have correlated diet vs primary data:
aucho_subfam_undated_boots_reduced<-list()
for(i in 1:length(aucho_subfam_undated_boots)){
  aucho_subfam_undated_boots_reduced[[i]]<-drop.tip(aucho_subfam_undated_boots[[i]],
                                                    tip = aucho_subfam_undated_boots[[i]]$tip.label[!aucho_subfam_undated_boots[[i]]$tip.label %in% names(vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont)])
}

#Repeat the key diet vs primary symbiont on each of them. 
length(aucho_subfam_undated_boots_reduced)
corevol_diet_xylem_primary.endosymbiont_undated_subfam_boots_reduced<-list()
for(i in 1:length(aucho_subfam_undated_boots_reduced)){
  corevol_diet_xylem_primary.endosymbiont_undated_subfam_boots_reduced[[i]]<-
    fitPagel(tree = aucho_subfam_undated_boots_reduced[[i]],
             x=vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont,
             y=vec_fitPagel_diet_xylem_primary.endosymbiont_diet_xylem)
}  

save(corevol_diet_xylem_primary.endosymbiont_undated_subfam_boots_reduced,
     file="./Output/corevol_diet_xylem_primary.endosymbiont_undated_subfam_boots_reduced")

##Analyse the differences
#Correlated model
p_vec<-NULL
aic_dif_vec<-NULL
corevol_diet_xylem_primary.endosymbiont_undated_subfam_boots_reduced[[1]]
for(i in 1:length(corevol_diet_xylem_primary.endosymbiont_undated_subfam_boots_reduced)){
  p_vec[i]<-corevol_diet_xylem_primary.endosymbiont_undated_subfam_boots_reduced[[i]]$P
  aic_dif_vec[i]<-corevol_diet_xylem_primary.endosymbiont_undated_subfam_boots_reduced[[i]]$independent.AIC-
    corevol_diet_xylem_primary.endosymbiont_undated_subfam_boots_reduced[[i]]$dependent.AIC
}

#Calculate the differences for the best model
best_mod_p<-corevol_diet_xylem_primary.endosymbiont_dated_subfam$P
best_mod_aic_dif<-corevol_diet_xylem_primary.endosymbiont_dated_subfam$independent.AIC-
  corevol_diet_xylem_primary.endosymbiont_dated_subfam$dependent.AIC
best_mod_p
best_mod_aic_dif

aic_dif_density_plot<-
  ggplot()+
  geom_density(aes(x = aic_dif_vec))+
  geom_vline(aes(xintercept=best_mod_aic_dif),color="red")+
  geom_vline(aes(xintercept=median(aic_dif_vec)),color="blue")+
  xlim(10,25)+
  xlab("delta AIC")+
  theme_bw()+
  theme(axis.title.x = element_text(size = rel(1.5)),axis.title.y = element_text(size = rel(1.5)))+
  theme(axis.text.x = element_text(size = rel(1.75)),axis.text.y = element_text(size = rel(1.75)))
aic_dif_density_plot

p_density_plot<-
  ggplot()+
  geom_density(aes(x = p_vec))+
  geom_vline(aes(xintercept=best_mod_p),color="red")+
  geom_vline(aes(xintercept=median(p_vec)),color="blue")+
  xlim(0,0.0004)+
  xlab("LRT p-value")+
  theme_bw()+
  theme(axis.title.x = element_text(size = rel(1.5)),axis.title.y = element_text(size = rel(1.5)))+
  theme(axis.text.x = element_text(size = rel(1.75)),axis.text.y = element_text(size = rel(1.75)))
p_density_plot

library(gridExtra)
png("./Output/FigureS6.png",width = 800)
grid.arrange(p_density_plot,aic_dif_density_plot,ncol=2)
dev.off()

#########
#2. Data uncertainty  -----------------------------------------------------

#Check the two relevant data vectors
vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont
vec_fitPagel_diet_xylem_primary.endosymbiont_diet_xylem
small_aucho_subfam_dated_corevol_diet_xylem_primary.endosymbiont<-drop.tip(small_aucho_subfam_dated,
                                                                           tip = small_aucho_subfam_dated$tip.label[!small_aucho_subfam_dated$tip.label %in% dat_fitPagel_diet_xylem_primary.endosymbiont$genus])

#Range of different fp anf fs
correct_prob_negatives<-0.95
correct_prob_positives<-0.95
corevol_diet_xylem_primary.endosymbiont_dated_subfam_05fn_05fp<-list()
for(i in 1:100){
  absence_mutate<-1-rbinom(n=length(which(vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont==0)),
                           size=1,prob = correct_prob_negatives)
  presence_mutate<-rbinom(n=length(which(vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont==1)),
                          size=1,prob = correct_prob_positives)
  vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont_mutated<- vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont
  vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont_mutated[vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont_mutated==0]<-absence_mutate
  corevol_diet_xylem_primary.endosymbiont_dated_subfam_05fn_05fp[[i]]<-
    fitPagel(tree = small_aucho_subfam_dated_corevol_diet_xylem_primary.endosymbiont,
             x = vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont_mutated,
             y = vec_fitPagel_diet_xylem_primary.endosymbiont_diet_xylem)
}

save(corevol_diet_xylem_primary.endosymbiont_dated_subfam_05fn_05fp,
     file="./Output/corevol_diet_xylem_primary.endosymbiont_dated_subfam_05fn_05fp")

correct_prob_negatives<-0.95
correct_prob_positives<-0.85
corevol_diet_xylem_primary.endosymbiont_dated_subfam_05fn_15fp<-list()
for(i in 1:100){
  absence_mutate<-1-rbinom(n=length(which(vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont==0)),
                           size=1,prob = correct_prob_negatives)
  presence_mutate<-rbinom(n=length(which(vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont==1)),
                          size=1,prob = correct_prob_positives)
  vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont_mutated<- vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont
  vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont_mutated[vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont_mutated==0]<-absence_mutate
  corevol_diet_xylem_primary.endosymbiont_dated_subfam_05fn_15fp[[i]]<-
    fitPagel(tree = small_aucho_subfam_dated_corevol_diet_xylem_primary.endosymbiont,
             x = vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont_mutated,
             y = vec_fitPagel_diet_xylem_primary.endosymbiont_diet_xylem)
}

save(corevol_diet_xylem_primary.endosymbiont_dated_subfam_05fn_15fp,
     file="./Output/corevol_diet_xylem_primary.endosymbiont_dated_subfam_05fn_15fp")


correct_prob_negatives<-0.95
correct_prob_positives<-0.75
corevol_diet_xylem_primary.endosymbiont_dated_subfam_05fn_25fp<-list()
for(i in 1:100){
  absence_mutate<-1-rbinom(n=length(which(vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont==0)),
                           size=1,prob = correct_prob_negatives)
  presence_mutate<-rbinom(n=length(which(vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont==1)),
                          size=1,prob = correct_prob_positives)
  vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont_mutated<- vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont
  vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont_mutated[vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont_mutated==0]<-absence_mutate
  corevol_diet_xylem_primary.endosymbiont_dated_subfam_05fn_25fp[[i]]<-
    fitPagel(tree = small_aucho_subfam_dated_corevol_diet_xylem_primary.endosymbiont,
             x = vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont_mutated,
             y = vec_fitPagel_diet_xylem_primary.endosymbiont_diet_xylem)
}

save(corevol_diet_xylem_primary.endosymbiont_dated_subfam_05fn_25fp,
     file="./Output/corevol_diet_xylem_primary.endosymbiont_dated_subfam_05fn_25fp")

correct_prob_negatives<-0.85
correct_prob_positives<-0.95
corevol_diet_xylem_primary.endosymbiont_dated_subfam_15fn_05fp<-list()
for(i in 1:100){
  absence_mutate<-1-rbinom(n=length(which(vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont==0)),
                           size=1,prob = correct_prob_negatives)
  presence_mutate<-rbinom(n=length(which(vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont==1)),
                          size=1,prob = correct_prob_positives)
  vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont_mutated<- vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont
  vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont_mutated[vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont_mutated==0]<-absence_mutate
  corevol_diet_xylem_primary.endosymbiont_dated_subfam_15fn_05fp[[i]]<-
    fitPagel(tree = small_aucho_subfam_dated_corevol_diet_xylem_primary.endosymbiont,
             x = vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont_mutated,
             y = vec_fitPagel_diet_xylem_primary.endosymbiont_diet_xylem)
}

save(corevol_diet_xylem_primary.endosymbiont_dated_subfam_15fn_05fp,
     file="./Output/corevol_diet_xylem_primary.endosymbiont_dated_subfam_15fn_05fp")

correct_prob_negatives<-0.85
correct_prob_positives<-0.85
corevol_diet_xylem_primary.endosymbiont_dated_subfam_15fn_15fp<-list()
for(i in 1:100){
  absence_mutate<-1-rbinom(n=length(which(vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont==0)),
                           size=1,prob = correct_prob_negatives)
  presence_mutate<-rbinom(n=length(which(vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont==1)),
                          size=1,prob = correct_prob_positives)
  vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont_mutated<- vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont
  vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont_mutated[vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont_mutated==0]<-absence_mutate
  corevol_diet_xylem_primary.endosymbiont_dated_subfam_15fn_15fp[[i]]<-
    fitPagel(tree = small_aucho_subfam_dated_corevol_diet_xylem_primary.endosymbiont,
             x = vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont_mutated,
             y = vec_fitPagel_diet_xylem_primary.endosymbiont_diet_xylem)
}

save(corevol_diet_xylem_primary.endosymbiont_dated_subfam_15fn_15fp,
     file="./Output/corevol_diet_xylem_primary.endosymbiont_dated_subfam_15fn_15fp")

correct_prob_negatives<-0.85
correct_prob_positives<-0.75
corevol_diet_xylem_primary.endosymbiont_dated_subfam_15fn_25fp<-list()
for(i in 1:100){
  absence_mutate<-1-rbinom(n=length(which(vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont==0)),
                           size=1,prob = correct_prob_negatives)
  presence_mutate<-rbinom(n=length(which(vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont==1)),
                          size=1,prob = correct_prob_positives)
  vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont_mutated<- vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont
  vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont_mutated[vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont_mutated==0]<-absence_mutate
  corevol_diet_xylem_primary.endosymbiont_dated_subfam_15fn_25fp[[i]]<-
    fitPagel(tree = small_aucho_subfam_dated_corevol_diet_xylem_primary.endosymbiont,
             x = vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont_mutated,
             y = vec_fitPagel_diet_xylem_primary.endosymbiont_diet_xylem)
}

save(corevol_diet_xylem_primary.endosymbiont_dated_subfam_15fn_25fp,
     file="./Output/corevol_diet_xylem_primary.endosymbiont_dated_subfam_15fn_25fp")

correct_prob_negatives<-0.75
correct_prob_positives<-0.95
corevol_diet_xylem_primary.endosymbiont_dated_subfam_25fn_05fp<-list()
for(i in 1:100){
  absence_mutate<-1-rbinom(n=length(which(vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont==0)),
                           size=1,prob = correct_prob_negatives)
  presence_mutate<-rbinom(n=length(which(vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont==1)),
                          size=1,prob = correct_prob_positives)
  vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont_mutated<- vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont
  vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont_mutated[vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont_mutated==0]<-absence_mutate
  corevol_diet_xylem_primary.endosymbiont_dated_subfam_25fn_05p[[i]]<-
    fitPagel(tree = small_aucho_subfam_dated_corevol_diet_xylem_primary.endosymbiont,
             x = vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont_mutated,
             y = vec_fitPagel_diet_xylem_primary.endosymbiont_diet_xylem)
}

save(corevol_diet_xylem_primary.endosymbiont_dated_subfam_25fn_05fp,
     file="./Output/corevol_diet_xylem_primary.endosymbiont_dated_subfam_25fn_05fp")

correct_prob_negatives<-0.75
correct_prob_positives<-0.85
corevol_diet_xylem_primary.endosymbiont_dated_subfam_25fn_15p<-list()
for(i in 1:100){
  absence_mutate<-1-rbinom(n=length(which(vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont==0)),
                           size=1,prob = correct_prob_negatives)
  presence_mutate<-rbinom(n=length(which(vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont==1)),
                          size=1,prob = correct_prob_positives)
  vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont_mutated<- vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont
  vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont_mutated[vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont_mutated==0]<-absence_mutate
  corevol_diet_xylem_primary.endosymbiont_dated_subfam_25fn_15p[[i]]<-
    fitPagel(tree = small_aucho_subfam_dated_corevol_diet_xylem_primary.endosymbiont,
             x = vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont_mutated,
             y = vec_fitPagel_diet_xylem_primary.endosymbiont_diet_xylem)
}

save(corevol_diet_xylem_primary.endosymbiont_dated_subfam_25fn_15fp,
     file="./Output/corevol_diet_xylem_primary.endosymbiont_dated_subfam_25fn_15fp")

correct_prob_negatives<-0.75
correct_prob_positives<-0.75
corevol_diet_xylem_primary.endosymbiont_dated_subfam_25fn_25p<-list()
for(i in 1:100){
  absence_mutate<-1-rbinom(n=length(which(vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont==0)),
                           size=1,prob = correct_prob_negatives)
  presence_mutate<-rbinom(n=length(which(vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont==1)),
                          size=1,prob = correct_prob_positives)
  vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont_mutated<- vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont
  vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont_mutated[vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont_mutated==0]<-absence_mutate
  corevol_diet_xylem_primary.endosymbiont_dated_subfam_25fn_25p[[i]]<-
    fitPagel(tree = small_aucho_subfam_dated_corevol_diet_xylem_primary.endosymbiont,
             x = vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont_mutated,
             y = vec_fitPagel_diet_xylem_primary.endosymbiont_diet_xylem)
}

save(corevol_diet_xylem_primary.endosymbiont_dated_subfam_25fn_25fp,
     file="./Output/corevol_diet_xylem_primary.endosymbiont_dated_subfam_25fn_25fp")
save.image()

###Analyse & visualise data uncertainty
#Extract the aic/p-value, and visualise the results per reconstruction (9 in total, each replicated 100 times)

p_vec_05fn_05fp<-NULL
aic_dif_vec_05fn_05fp<-NULL
for(i in 1:length(corevol_diet_xylem_primary.endosymbiont_dated_subfam_05fn_05fp)){
  p_vec_05fn_05fp[i]<-corevol_diet_xylem_primary.endosymbiont_dated_subfam_05fn_05fp[[i]]$P
  aic_dif_vec_05fn_05fp[i]<-corevol_diet_xylem_primary.endosymbiont_dated_subfam_05fn_05fp[[i]]$independent.AIC-
    corevol_diet_xylem_primary.endosymbiont_dated_subfam_05fn_05fp[[i]]$dependent.AIC
}

aic_dif_density_plot_05fn_05fp<-
  ggplot()+
  geom_density(aes(x = aic_dif_vec_05fn_05fp))+
  geom_vline(aes(xintercept=best_mod_aic_dif),color="red")+
  geom_vline(aes(xintercept=median(aic_dif_vec_05fn_05fp)),color="blue")+
  xlim(5,30)+
  xlab("delta AIC")+
  theme_bw()+
  theme(axis.title.x = element_text(size = rel(1.5)),axis.title.y = element_text(size = rel(1.5)))+
  theme(axis.text.x = element_text(size = rel(1.75)),axis.text.y = element_text(size = rel(1.75)))+
  theme(plot.title = element_text(size = rel(1.5)))+
  ggtitle("A. 5% false positive, 5% false negative")
aic_dif_density_plot_05fn_05fp

p_vec_05fn_15fp<-NULL
aic_dif_vec_05fn_15fp<-NULL
for(i in 1:length(corevol_diet_xylem_primary.endosymbiont_dated_subfam_05fn_15fp)){
  p_vec_05fn_15fp[i]<-corevol_diet_xylem_primary.endosymbiont_dated_subfam_05fn_15fp[[i]]$P
  aic_dif_vec_05fn_15fp[i]<-corevol_diet_xylem_primary.endosymbiont_dated_subfam_05fn_15fp[[i]]$independent.AIC-
    corevol_diet_xylem_primary.endosymbiont_dated_subfam_05fn_15fp[[i]]$dependent.AIC
}

aic_dif_density_plot_05fn_15fp<-
  ggplot()+
  geom_density(aes(x = aic_dif_vec_05fn_15fp))+
  geom_vline(aes(xintercept=best_mod_aic_dif),color="red")+
  geom_vline(aes(xintercept=median(aic_dif_vec_05fn_15fp)),color="blue")+
  xlim(5,30)+
  xlab("delta AIC")+
  theme_bw()+
  theme(axis.title.x = element_text(size = rel(1.5)),axis.title.y = element_text(size = rel(1.5)))+
  theme(axis.text.x = element_text(size = rel(1.75)),axis.text.y = element_text(size = rel(1.75)))+
  theme(plot.title = element_text(size = rel(1.5)))+
  ggtitle("B. 5% false positive, 15% false negative")
aic_dif_density_plot_05fn_15fp

p_vec_05fn_25fp<-NULL
aic_dif_vec_05fn_25fp<-NULL
for(i in 1:length(corevol_diet_xylem_primary.endosymbiont_dated_subfam_05fn_25fp)){
  p_vec_05fn_25fp[i]<-corevol_diet_xylem_primary.endosymbiont_dated_subfam_05fn_25fp[[i]]$P
  aic_dif_vec_05fn_25fp[i]<-corevol_diet_xylem_primary.endosymbiont_dated_subfam_05fn_25fp[[i]]$independent.AIC-
    corevol_diet_xylem_primary.endosymbiont_dated_subfam_05fn_25fp[[i]]$dependent.AIC
}

aic_dif_density_plot_05fn_25fp<-
  ggplot()+
  geom_density(aes(x = aic_dif_vec_05fn_25fp))+
  geom_vline(aes(xintercept=best_mod_aic_dif),color="red")+
  geom_vline(aes(xintercept=median(aic_dif_vec_05fn_25fp)),color="blue")+
  xlim(5,30)+
  xlab("delta AIC")+
  theme_bw()+  
  theme(axis.title.x = element_text(size = rel(1.5)),axis.title.y = element_text(size = rel(1.5)))+
  theme(axis.text.x = element_text(size = rel(1.75)),axis.text.y = element_text(size = rel(1.75)))+
  theme(plot.title = element_text(size = rel(1.5)))+
  ggtitle("C. 5% false positive, 25% false negative")
aic_dif_density_plot_05fn_25fp

p_vec_15fn_05fp<-NULL
aic_dif_vec_15fn_05fp<-NULL
for(i in 1:length(corevol_diet_xylem_primary.endosymbiont_dated_subfam_15fn_05fp)){
  p_vec_15fn_05fp[i]<-corevol_diet_xylem_primary.endosymbiont_dated_subfam_15fn_05fp[[i]]$P
  aic_dif_vec_15fn_05fp[i]<-corevol_diet_xylem_primary.endosymbiont_dated_subfam_15fn_05fp[[i]]$independent.AIC-
    corevol_diet_xylem_primary.endosymbiont_dated_subfam_15fn_05fp[[i]]$dependent.AIC
}

aic_dif_density_plot_15fn_05fp<-
  ggplot()+
  geom_density(aes(x = aic_dif_vec_15fn_05fp))+
  geom_vline(aes(xintercept=best_mod_aic_dif),color="red")+
  geom_vline(aes(xintercept=median(aic_dif_vec_15fn_05fp)),color="blue")+
  xlim(5,30)+
  xlab("delta AIC")+
  theme_bw()+
  theme(axis.title.x = element_text(size = rel(1.5)),axis.title.y = element_text(size = rel(1.5)))+
  theme(axis.text.x = element_text(size = rel(1.75)),axis.text.y = element_text(size = rel(1.75)))+
  theme(plot.title = element_text(size = rel(1.5)))+
  ggtitle("D. 15% false positive, 5% false negative")
aic_dif_density_plot_15fn_05fp

p_vec_15fn_15fp<-NULL
aic_dif_vec_15fn_15fp<-NULL
for(i in 1:length(corevol_diet_xylem_primary.endosymbiont_dated_subfam_15fn_15fp)){
  p_vec_15fn_15fp[i]<-corevol_diet_xylem_primary.endosymbiont_dated_subfam_15fn_15fp[[i]]$P
  aic_dif_vec_15fn_15fp[i]<-corevol_diet_xylem_primary.endosymbiont_dated_subfam_15fn_15fp[[i]]$independent.AIC-
    corevol_diet_xylem_primary.endosymbiont_dated_subfam_15fn_15fp[[i]]$dependent.AIC
}

aic_dif_density_plot_15fn_15fp<-
  ggplot()+
  geom_density(aes(x = aic_dif_vec_15fn_15fp))+
  geom_vline(aes(xintercept=best_mod_aic_dif),color="red")+
  geom_vline(aes(xintercept=median(aic_dif_vec_15fn_15fp)),color="blue")+
  xlim(5,30)+
  xlab("delta AIC")+
  theme_bw()+
  theme(axis.title.x = element_text(size = rel(1.5)),axis.title.y = element_text(size = rel(1.5)))+
  theme(axis.text.x = element_text(size = rel(1.75)),axis.text.y = element_text(size = rel(1.75)))+
  theme(plot.title = element_text(size = rel(1.5)))+
  ggtitle("E. 15% false positive, 15% false negative")
aic_dif_density_plot_15fn_15fp

p_vec_15fn_25fp<-NULL
aic_dif_vec_15fn_25fp<-NULL
for(i in 1:length(corevol_diet_xylem_primary.endosymbiont_dated_subfam_15fn_25fp)){
  p_vec_15fn_25fp[i]<-corevol_diet_xylem_primary.endosymbiont_dated_subfam_15fn_25fp[[i]]$P
  aic_dif_vec_15fn_25fp[i]<-corevol_diet_xylem_primary.endosymbiont_dated_subfam_15fn_25fp[[i]]$independent.AIC-
    corevol_diet_xylem_primary.endosymbiont_dated_subfam_15fn_25fp[[i]]$dependent.AIC
}

aic_dif_density_plot_15fn_25fp<-
  ggplot()+
  geom_density(aes(x = aic_dif_vec_15fn_25fp))+
  geom_vline(aes(xintercept=best_mod_aic_dif),color="red")+
  geom_vline(aes(xintercept=median(aic_dif_vec_15fn_25fp)),color="blue")+
  xlim(5,30)+
  xlab("delta AIC")+
  theme_bw()+
  theme(axis.title.x = element_text(size = rel(1.5)),axis.title.y = element_text(size = rel(1.5)))+
  theme(axis.text.x = element_text(size = rel(1.75)),axis.text.y = element_text(size = rel(1.75)))+
  theme(plot.title = element_text(size = rel(1.5)))+
  ggtitle("F. 15% false positive, 25% false negative")
aic_dif_density_plot_15fn_25fp

p_vec_25fn_05fp<-NULL
aic_dif_vec_25fn_05fp<-NULL
for(i in 1:length(corevol_diet_xylem_primary.endosymbiont_dated_subfam_25fn_05fp)){
  p_vec_25fn_05fp[i]<-corevol_diet_xylem_primary.endosymbiont_dated_subfam_25fn_05fp[[i]]$P
  aic_dif_vec_25fn_05fp[i]<-corevol_diet_xylem_primary.endosymbiont_dated_subfam_25fn_05fp[[i]]$independent.AIC-
    corevol_diet_xylem_primary.endosymbiont_dated_subfam_25fn_05fp[[i]]$dependent.AIC
}

aic_dif_density_plot_25fn_05fp<-
  ggplot()+
  geom_density(aes(x = aic_dif_vec_25fn_05fp))+
  geom_vline(aes(xintercept=best_mod_aic_dif),color="red")+
  geom_vline(aes(xintercept=median(aic_dif_vec_25fn_05fp)),color="blue")+
  xlim(5,30)+
  xlab("delta AIC")+
  theme_bw()+
  theme(axis.title.x = element_text(size = rel(1.5)),axis.title.y = element_text(size = rel(1.5)))+
  theme(axis.text.x = element_text(size = rel(1.75)),axis.text.y = element_text(size = rel(1.75)))+
  theme(plot.title = element_text(size = rel(1.5)))+
  ggtitle("G. 25% false positive, 5% false negative")
aic_dif_density_plot_25fn_05fp

p_vec_25fn_15fp<-NULL
aic_dif_vec_25fn_15fp<-NULL
for(i in 1:length(corevol_diet_xylem_primary.endosymbiont_dated_subfam_25fn_15fp)){
  p_vec_25fn_15fp[i]<-corevol_diet_xylem_primary.endosymbiont_dated_subfam_25fn_15fp[[i]]$P
  aic_dif_vec_25fn_15fp[i]<-corevol_diet_xylem_primary.endosymbiont_dated_subfam_25fn_15fp[[i]]$independent.AIC-
    corevol_diet_xylem_primary.endosymbiont_dated_subfam_25fn_15fp[[i]]$dependent.AIC
}

aic_dif_density_plot_25fn_15fp<-
  ggplot()+
  geom_density(aes(x = aic_dif_vec_25fn_15fp))+
  geom_vline(aes(xintercept=best_mod_aic_dif),color="red")+
  geom_vline(aes(xintercept=median(aic_dif_vec_25fn_15fp)),color="blue")+
  xlim(5,30)+
  xlab("delta AIC")+
  theme_bw()+
  theme(axis.title.x = element_text(size = rel(1.5)),axis.title.y = element_text(size = rel(1.5)))+
  theme(axis.text.x = element_text(size = rel(1.75)),axis.text.y = element_text(size = rel(1.75)))+
  theme(plot.title = element_text(size = rel(1.5)))+
  ggtitle("H. 25% false positive, 15% false negative")
aic_dif_density_plot_25fn_15fp

p_vec_25fn_25fp<-NULL
aic_dif_vec_25fn_25fp<-NULL
for(i in 1:length(corevol_diet_xylem_primary.endosymbiont_dated_subfam_25fn_25fp)){
  p_vec_25fn_25fp[i]<-corevol_diet_xylem_primary.endosymbiont_dated_subfam_25fn_25fp[[i]]$P
  aic_dif_vec_25fn_25fp[i]<-corevol_diet_xylem_primary.endosymbiont_dated_subfam_25fn_25fp[[i]]$independent.AIC-
    corevol_diet_xylem_primary.endosymbiont_dated_subfam_25fn_25fp[[i]]$dependent.AIC
}

aic_dif_density_plot_25fn_25fp<-
  ggplot()+
  geom_density(aes(x = aic_dif_vec_25fn_25fp))+
  geom_vline(aes(xintercept=best_mod_aic_dif),color="red")+
  geom_vline(aes(xintercept=median(aic_dif_vec_25fn_25fp)),color="blue")+
  xlim(5,30)+
  xlab("delta AIC")+
  theme_bw()+
  theme(axis.title.x = element_text(size = rel(1.5)),axis.title.y = element_text(size = rel(1.5)))+
  theme(axis.text.x = element_text(size = rel(1.75)),axis.text.y = element_text(size = rel(1.75)))+
  theme(plot.title = element_text(size = rel(1.5)))+
  ggtitle("I. 25% false positive, 25% false negative")
aic_dif_density_plot_25fn_25fp


###Create figures combining all the above in a single figure
library(gridExtra)
png("./Output/FigureS7.png",width = 1400,height = 1800)
grid.arrange(aic_dif_density_plot_05fn_05fp,aic_dif_density_plot_05fn_15fp,aic_dif_density_plot_05fn_25fp,
             aic_dif_density_plot_15fn_05fp,aic_dif_density_plot_15fn_15fp,aic_dif_density_plot_15fn_25fp,
             aic_dif_density_plot_25fn_05fp,aic_dif_density_plot_25fn_15fp,aic_dif_density_plot_25fn_25fp,
             ncol=3)
dev.off()

#########
#3. Genus Sampling  -----------------------------------------------------

#Do they have the same order?
names(vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont)==names(vec_fitPagel_diet_xylem_primary.endosymbiont_diet_xylem)

corevol_diet_xylem_primary.endosymbiont_dated_genus_sampling_10loss<-list()
for(i in 1:100){
  loss_num<-sample(x=1:length(vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont),
                   size = round(length(vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont)*0.10))
  vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont_reduced<-vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont[-loss_num]
  vec_fitPagel_diet_xylem_primary.endosymbiont_diet_xylem_reduced<-vec_fitPagel_diet_xylem_primary.endosymbiont_diet_xylem[-loss_num]
  small_aucho_subfam_dated_corevol_diet_xylem_primary.endosymbiont_reduced<-
    drop.tip(phy = small_aucho_subfam_dated_corevol_diet_xylem_primary.endosymbiont,
             tip = small_aucho_subfam_dated_corevol_diet_xylem_primary.endosymbiont$tip.label[!small_aucho_subfam_dated_corevol_diet_xylem_primary.endosymbiont$tip.label %in% 
                                                                                                names(vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont_reduced)])
  corevol_diet_xylem_primary.endosymbiont_dated_genus_sampling_10loss[[i]]<-
    fitPagel(tree = small_aucho_subfam_dated_corevol_diet_xylem_primary.endosymbiont_reduced,
             x = vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont_reduced,
             y = vec_fitPagel_diet_xylem_primary.endosymbiont_diet_xylem_reduced)
}

save(corevol_diet_xylem_primary.endosymbiont_dated_genus_sampling_10loss,
     file="./Output/corevol_diet_xylem_primary.endosymbiont_dated_genus_sampling_10loss")

corevol_diet_xylem_primary.endosymbiont_dated_genus_sampling_20loss<-list()
for(i in 1:100){
  loss_num<-sample(x=1:length(vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont),
                   size = round(length(vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont)*0.20))
  vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont_reduced<-vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont[-loss_num]
  vec_fitPagel_diet_xylem_primary.endosymbiont_diet_xylem_reduced<-vec_fitPagel_diet_xylem_primary.endosymbiont_diet_xylem[-loss_num]
  small_aucho_subfam_dated_corevol_diet_xylem_primary.endosymbiont_reduced<-
    drop.tip(phy = small_aucho_subfam_dated_corevol_diet_xylem_primary.endosymbiont,
             tip = small_aucho_subfam_dated_corevol_diet_xylem_primary.endosymbiont$tip.label[!small_aucho_subfam_dated_corevol_diet_xylem_primary.endosymbiont$tip.label %in% 
                                                                                                names(vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont_reduced)])
  corevol_diet_xylem_primary.endosymbiont_dated_genus_sampling_20loss[[i]]<-
    fitPagel(tree = small_aucho_subfam_dated_corevol_diet_xylem_primary.endosymbiont_reduced,
             x = vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont_reduced,
             y = vec_fitPagel_diet_xylem_primary.endosymbiont_diet_xylem_reduced)
}

save(corevol_diet_xylem_primary.endosymbiont_dated_genus_sampling_20loss,
     file="./Output/corevol_diet_xylem_primary.endosymbiont_dated_genus_sampling_20loss")

corevol_diet_xylem_primary.endosymbiont_dated_genus_sampling_30loss<-list()
for(i in 1:100){
  loss_num<-sample(x=1:length(vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont),
                   size = round(length(vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont)*0.30))
  vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont_reduced<-vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont[-loss_num]
  vec_fitPagel_diet_xylem_primary.endosymbiont_diet_xylem_reduced<-vec_fitPagel_diet_xylem_primary.endosymbiont_diet_xylem[-loss_num]
  small_aucho_subfam_dated_corevol_diet_xylem_primary.endosymbiont_reduced<-
    drop.tip(phy = small_aucho_subfam_dated_corevol_diet_xylem_primary.endosymbiont,
             tip = small_aucho_subfam_dated_corevol_diet_xylem_primary.endosymbiont$tip.label[!small_aucho_subfam_dated_corevol_diet_xylem_primary.endosymbiont$tip.label %in% 
                                                                                                names(vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont_reduced)])
  corevol_diet_xylem_primary.endosymbiont_dated_genus_sampling_30loss[[i]]<-
    fitPagel(tree = small_aucho_subfam_dated_corevol_diet_xylem_primary.endosymbiont_reduced,
             x = vec_fitPagel_diet_xylem_primary.endosymbiont_primary.endosymbiont_reduced,
             y = vec_fitPagel_diet_xylem_primary.endosymbiont_diet_xylem_reduced)
}

save(corevol_diet_xylem_primary.endosymbiont_dated_genus_sampling_30loss,
     file="./Output/corevol_diet_xylem_primary.endosymbiont_dated_genus_sampling_30loss")

save.image()

###Analyse and visualise the differences for the three genus sampling reconstructions (each replicated 100 times)
p_vec_sampling_10loss<-NULL
aic_dif_vec_sampling_10loss<-NULL
for(i in 1:length(corevol_diet_xylem_primary.endosymbiont_dated_genus_sampling_10loss)){
  p_vec_sampling_10loss[i]<-corevol_diet_xylem_primary.endosymbiont_dated_genus_sampling_10loss[[i]]$P
  aic_dif_vec_sampling_10loss[i]<-corevol_diet_xylem_primary.endosymbiont_dated_genus_sampling_10loss[[i]]$independent.AIC-
    corevol_diet_xylem_primary.endosymbiont_dated_genus_sampling_10loss[[i]]$dependent.AIC
}

aic_dif_density_plot_sampling_10loss<-
  ggplot()+
  geom_density(aes(x = aic_dif_vec_sampling_10loss))+
  geom_vline(aes(xintercept=best_mod_aic_dif),color="red")+
  geom_vline(aes(xintercept=median(aic_dif_vec_sampling_10loss)),color="blue")+
  xlim(0,25)+
  xlab("delta AIC")+
  theme_bw()+
  ggtitle("A. Genus loss 10%")
aic_dif_density_plot_sampling_10loss

p_vec_sampling_20loss<-NULL
aic_dif_vec_sampling_20loss<-NULL
for(i in 1:length(corevol_diet_xylem_primary.endosymbiont_dated_genus_sampling_20loss)){
  p_vec_sampling_20loss[i]<-corevol_diet_xylem_primary.endosymbiont_dated_genus_sampling_20loss[[i]]$P
  aic_dif_vec_sampling_20loss[i]<-corevol_diet_xylem_primary.endosymbiont_dated_genus_sampling_20loss[[i]]$independent.AIC-
    corevol_diet_xylem_primary.endosymbiont_dated_genus_sampling_20loss[[i]]$dependent.AIC
}

aic_dif_density_plot_sampling_20loss<-
  ggplot()+
  geom_density(aes(x = aic_dif_vec_sampling_20loss))+
  geom_vline(aes(xintercept=best_mod_aic_dif),color="red")+
  geom_vline(aes(xintercept=median(aic_dif_vec_sampling_20loss)),color="blue")+
  xlim(0,25)+
  xlab("delta AIC")+
  theme_bw()+
  ggtitle("B. Genus loss 20%")
aic_dif_density_plot_sampling_20loss

p_vec_sampling_30loss<-NULL
aic_dif_vec_sampling_30loss<-NULL
for(i in 1:length(corevol_diet_xylem_primary.endosymbiont_dated_genus_sampling_30loss)){
  p_vec_sampling_30loss[i]<-corevol_diet_xylem_primary.endosymbiont_dated_genus_sampling_30loss[[i]]$P
  aic_dif_vec_sampling_30loss[i]<-corevol_diet_xylem_primary.endosymbiont_dated_genus_sampling_30loss[[i]]$independent.AIC-
    corevol_diet_xylem_primary.endosymbiont_dated_genus_sampling_30loss[[i]]$dependent.AIC
}

aic_dif_density_plot_sampling_30loss<-
  ggplot()+
  geom_density(aes(x = aic_dif_vec_sampling_30loss))+
  geom_vline(aes(xintercept=best_mod_aic_dif),color="red")+
  geom_vline(aes(xintercept=median(aic_dif_vec_sampling_30loss)),color="blue")+
  xlim(0,25)+
  xlab("delta AIC")+
  theme_bw()+
  ggtitle("C. Genus loss 30%")
aic_dif_density_plot_sampling_30loss

png("./Output/FigureS8.png",width = 1200)
grid.arrange(aic_dif_density_plot_sampling_10loss,aic_dif_density_plot_sampling_20loss,aic_dif_density_plot_sampling_30loss,
             ncol=3)
dev.off()
