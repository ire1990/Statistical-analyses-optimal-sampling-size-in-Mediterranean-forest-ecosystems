###############################################################
##                                                           ##
##      Sampling forest soils to describe fungal diversity   ##
##   and composition. Which is the optimal sampling size in  ##
##     Mediterranean pure and mixed pine-oak forests?        ##                                  
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
##                                                           ##
##              Adamo et al. 2020  Journal Name              ##
##                                                           ##
###############################################################

#Last update: 2020
#Citation: XXcode citing X doi.
#This is the code used in Adamo et al. to perform the statistical analyses 
#that can be used in similar works.

###############################################################
#                        TABLE OF CONTENTS                    #
#                                                             #
#   Line  31: Paper briefing and aims                         #
#   Line  51: Required libraries                              #
#   Line  61: Simulating and Importing the data               #
#   Line  88: Formatting the data                             #                                                 #
#   Line  116: iNEXT extrapolated curves (aim 1)              #                                                     #
#   Line  171: Variance of Bray-Curtis matrix (aim 2)         # 
#   Line  248: NMDS of differences between sample pools (aim 2) #                                                  #
#   Line 376: Using the beta-indices (tbi, aim 3)             #
#                                                             #                                       #
###############################################################

#Briefing:####

#By using high throughput-DNA sequencing data (Illumina MiSeq), we identified the minimum number of pooled samples needed to obtain a reliable description of fungal communities
#in terms of diversity and composition in three different Mediterranean forests' stands.. 
#Twenty soil samples were randomly taken and kept separated in each of the three plots per forest type. After keeping one random sample which was not pooled, we obtained 5 composite samples 
#consisting of pools of 3, 6, 10, 15 and 20 samples. We then sequenced the ITS2 using Illumina MiSeq.
#We further tested:
#1. The effect of sample pooling on fungal diversity, the iNEXT function was used to build 
#rarefactions curves pooling together the individual samples. 

#2. The variance of Bray-Curtis matrix between the number of sample pools for each forest type was compared 
#using the betadisper function which is analogue to a Levene's test.

#3. Beta-diversity patterns and whether the core of most abundant fungal species 
#is maintained between sites, we evaluated for each pool the species (or abundances-per-species) 
#losses (B) and species gains (C) using the beta-indices (tbi function, Legendre, 2019, Ecology and Evolution). 
#We used one-sample pool per each forest (sample 1) as a reference, and we compared pools with 
#increasing number of samples (sample 3, 6, 10, 15 and 20) to identify species losses and gains.

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Required libraries####

library(vegan)
library(lattice)
library(ggplot2)
library(ggpubr)
library(adespatial)
library(iNEXT)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Data simulation####

#The original data is still used in different unpublished works, thus we simulate a similar dataset here.

spp_data<-read.table("OTU_table.txt", header=T)# similar dataset to the original and you can perform your own simulation
m_env<-read.table("attribute_matrix.txt", header=T)

sps <- 2800# we can tell the number of species we want in our dataset
datos  <- spp_data[,1:sps]#we take the number of species from the dataset we have previously imported, which corresponds here to 2800 species

vec <- c(as.matrix(datos))
vec2 <- round(runif(length(vec[vec>0]), min=0, max=100),0)#here we simulate the data changing the abundance of the species but keeping the same richness

vec[vec>0] <- vec2

df <- data.frame(matrix(vec, nrow(datos), sps))
names(df) <- paste0("sp", c(1:sps))
df                   
row.names(df) <- m_env$SAMPLE# this is the final simulated data
#read data

m_env<-read.table("attribute_matrix.txt", header=T)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Formatting the data####

#we transform the data per forest type for downstream analyses
mixed1 <- df[1:18,]
mixed1.hell <- decostand(mixed1, 'hell')

Mixed_att <- factor(m_env[1:18, -c(1:3)])
Mixed_att

Pinus_att <- factor(m_env[36:53,-c(1:3)])
Pinus_att
Quercus_att <- factor(m_env[19:35, -c(1:3)])
Quercus_att

Pinus <-df[36:53,]
P_env <- factor(m_env[36:53,-c(1:3)])
Pinus.hell <- decostand(Pinus, 'hell')

Quercus <- df[19:35,]

Q_env <- factor(m_env[19:35, -c(1:3)])
Q_env
Quercus.hell <- decostand(Quercus, 'hell')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# iNEXT extrapolated curves between forest types (aim 1) ####

# We use the function iNEXT to assess the differences in richness between the number of soil sample pools per each forest type
# therefore, we will set the order q in the iNEXT function to zero
# sub corresponds to the composite samples, and the number indicate the numbers of subsamples that represent the given composite sample.

#Mixed stands, we create one matrix for each composite sample number and we will do so for each forest type
Sub_1<- as.matrix(t(df[c(1,7,13), ])) 
Sub_3<- as.matrix(t(df[c(2,8,14), ]))
Sub_6 <- as.matrix(t(df[c(3,9,15),]))
Sub_10 <-  as.matrix(t(df[c(4,10,16), ]))
Sub_15 <-  as.matrix(t(df[ c(5,11,17), ]))
Sub_20 <-  as.matrix(t(df[c(6,12,18), ]))
mixed1_hill <-  list(SU.1 = Sub_1, SU.3 = Sub_3, SU.6 =Sub_6 , SU_10 = Sub_10 , SU_15 = Sub_15, SU_20 =Sub_20  )
typemixed1 <- lapply(mixed1_hill, as.abucount)
curve_mixed1 <- iNEXT((typemixed1 ), q=c(0), datatype = "abundance", endpoint = 100000 )#we want to test the differences in richness between composite samples 
# richness so we will only add in the function q=c(0), which is the hill number 0 that stands for richness
curve_mixed1$DataInfo
mixed1_plot <- ggiNEXT(curve_mixed1 , type=1, color="site") ;mixed1_plot

#rarefaction curves in Pinus stands
Sub_1<- as.matrix(t(df[c(36,42,48), ] ))
Sub_3<- as.matrix(t(df[c(37,43,49), ]))
Sub_6 <- as.matrix(t(df[c(38,44,50),]))
Sub_10 <-  as.matrix(t(df[c(39,45,51), ]))
Sub_15 <-  as.matrix(t(df[c(40,46,52),]))
Sub_20 <-  as.matrix(t(df[c(41,47,53),]))
Pinus_hill <-  list(s.1 = Sub_1, s.3 = Sub_3, s.6 = Sub_6 , s_10 = Sub_10, s_15 = Sub_15 , s_20 =Sub_20  )
typePinus <- lapply(Pinus_hill, as.abucount)
curve_Pinus <- iNEXT((typePinus), q=c(0), datatype = "abundance", endpoint = 140000)
curve_Pinus$DataInfo
Pinus_plot <- ggiNEXT(curve_Pinus , type=1,  color="site");Pinus_plot 
  

#rarefaction curves in Quercus stands
Sub_1<- as.matrix(t(df[c(19,25,31), ]) )
Sub_3<- as.matrix(t(df[c(20,26,32), ]))
Sub_6 <- as.matrix(t(df[c(21,27), ]))
Sub_10 <-  as.matrix(t(df[c(22,28), ]))
Sub_15 <-  as.matrix(t(df[c(23,29,34), ]))
Sub_20 <-  as.matrix(t(df[c(24,30,35), ]))
Quercus_hill<-  list(ss.1 = Sub_1, ss.3 = Sub_3, ss.6 = Sub_6 , ss_10 = Sub_10, ss_15 = Sub_15, ss_20 =Sub_20  )
typeQuercus = lapply(Quercus_hill, as.abucount)

curve_Quercus<- iNEXT((typeQuercus), q=0, datatype = "abundance", endpoint = 150000)
curve_Quercus$DataInfo
Quercus_plot <- ggiNEXT(curve_Quercus , type=1, color="site");Quercus_plot

ggarrange(Pinus_plot, Quercus_plot,mixed1_plot, labels = c("a)", "b)", "c)"), nrow = 1, ncol = 3, common.legend= TRUE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Betadisper to test the variance of Bray-Curtis matrix between the number of sample pools (Aim 2)####

par(mfrow=c(1,1))

#Mixed forest
Mixed <-  df[1:18,]
Mixed.hell <- decostand(Mixed, 'hell')
Mixed.m_beta <- vegdist(Mixed.hell, method = "bray")
mod_mix <- betadisper(Mixed.m_beta, Mixed_att , type = "centroid")
plot(mod_mix)
anova(mod_mix)

dfm <- data.frame(Distance_to_centroid=mod_mix$distances,Group=mod_mix$group)
groups <- mod_mix$group
m<- ggplot(data=dfm,aes(x=groups,y=Distance_to_centroid, fill= groups))+
  geom_boxplot(alpha=0.5)+
  scale_fill_brewer(palette = "Set1", name = "N. of sample pools", labels = c("1", "3", "6", "10", "15", "20")) +
  ggtitle("Mixed" ) +
  xlab("N. of sample pools") +
  scale_x_discrete(labels=c("1", "3", "6", "10", "15", "20"))+
  ylab("Distance to centroid") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position="top")+
  theme(panel.grid = element_blank());m


#now we do the same for Pinus.
Pinus_bd <- df[36:53,]
Pinus_bd.hell <- decostand(Pinus_bd , 'hell')
dist.Pinus_bd <- vegdist(Pinus_bd.hell, method = "bray")
modt_P <- betadisper(dist.Pinus_bd , P_env, type = "centroid")

spp_Pinus_bd <- data.frame(Distance_to_centroid=modt_P$distances,Group=modt_P$group)
hsd = TukeyHSD(modt_P)
groups <- modt_P$group


betadisp_Pinus<- ggplot(data=spp_Pinus_bd,aes(x=Group,y=Distance_to_centroid, fill= Group))+
  geom_boxplot(alpha=0.5)+
  scale_fill_brewer(palette = "Set1", name = "N. of sample pools", labels = c("1", "3", "6", "10", "15", "20")) +
  xlab("N. of sample pools") +
  scale_x_discrete(labels=c("1", "3", "6", "10", "15", "20"))+
  ggtitle("Pinus_s")+
  ylab("Distance to centroid") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face= "italic"))+
  theme(legend.position="top")+
  theme(panel.grid = element_blank());betadisp_Pinus
anova(modt_P)
TukeyHSD(modt_P)

# and now we perform the betadisper for Quercus
Quercus_bd <- df[19:35,]
Quercus_bd_hell <- decostand(Quercus_bd , 'hell')

dist.Quercus_bd   <- vegdist(Quercus_bd_hell, method = "bray")
modt_Q <- betadisper(dist.Quercus_bd, Q_env, type = "centroid")
anova(modt_Q)
TukeyHSD(modt_Q)
spp_Quercus_bd <- data.frame(Distance_to_centroid=modt_Q$distances,Group=modt_Q$group)
groups <- modt_Q$group
Quercus_betadisp<- ggplot(data=spp_Quercus_bd,aes(x=Group,y=Distance_to_centroid, fill= Group))+
  geom_boxplot(alpha=0.5)+
  scale_fill_brewer(palette = "Set1", name = "N. of sample pools", labels = c("1", "3", "6", "10", "15", "20")) +
  xlab("N. of sample pools") +
  ggtitle("Quercus_r") +
  scale_x_discrete(labels=c("1", "3", "6", "10", "15", "20"))+
  ylab("Distance to centroid") +

  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face= "italic"))+
  theme(legend.position="top")+
  theme(panel.grid = element_blank());Quercus_betadisp# we plotted the bestadisper results using ggplot2

# NMDS to diplay the lack of compositional differences between  number of soil sample pools  (aim 2) ####

# We use nmds to assess that there no difference in  species composition between sample pools (between the composite sample)
par(mfrow=c(1,1))

spnumfrec<- specnumber(df ,MARGIN=2)
spnumfrec

m_spp<-df[,spnumfrec>6]#we look at the species present in more than 10% of the sites
m_spp2_h <- decostand(m_spp, "hell")
nmds <- metaMDS(m_spp2_h, distance = "bray")

plot(nmds, type = "n")
points(nmds, pch=20, col=as.numeric(m_env$Subsample))


ordiellipse(nmds,m_env$Subsample ,show.groups="1",kind="sd",conf=0.95, col=6,lwd=2,lty=1,font=2,label = T)
ordiellipse(nmds, m_env$Subsample,show.groups="3",kind="sd",conf=0.95,col=2,lwd=2,lty=2,font=2,label = T)
ordiellipse(nmds, m_env$Subsample,show.groups="6",kind="sd",conf=0.95,col=3,lwd=2,font=2,label = T)
ordiellipse(nmds, m_env$Subsample,show.groups="10",kind="sd",conf=0.95, col=5,lwd=2,lty=2,font=2,label = T)
ordiellipse(nmds, m_env$Subsample,show.groups="15",kind="sd",conf=0.95, col=4,lwd=2,lty=3,font=2,label = T)
ordiellipse(nmds, m_env$Subsample,show.groups="20",kind="sd",conf=0.95, col=9,lwd=2,lty=3,font=2,label = T)
adonis(m_spp2_h~m_env$Subsample)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Assessing species losses (B) and species gains (C) using the beta-indices (tbi). Composite sample 1 used as reference and comparing the values with the other composite samples (aim 3)#### 

#Mixed forest
T1 <-  as.data.frame((df [c(1,7,13),]))
M3_T2 <-as.data.frame((df [c(2,8,14),]))
M6_T2 <- as.data.frame(df [c(3,9,15),])
M10_T2 <-as.data.frame(df [c(4,10,16), ])
M15_T2 <- as.data.frame(df[c(5,11,17), ])
M20_T2 <- as.data.frame(df [c(6,12,18), ])

#to get the permutation to work the function TBI must be changed to randomise it with n >2
#comparing composite sample 1 with composite sample 3
M3.TBI <- TBI(T1, M3_T2, method = "%diff", nperm = 999,  test.t.perm = TRUE)
M3.TBI$t.test_B.C #non significant p-value because  of the low number of samples used in the permutations
#comparing composite sample 1 with composite sample 6
M6.TBI <- TBI(T1,M6_T2, method = "%diff", nperm = 999, test.t.perm = TRUE)
M6prova.TBI$t.test_B.C
#comparing composite sample 1 with composite sample 10
M10.TBI <- TBI(T1,M10_T2, method = "%diff", nperm = 999, test.t.perm = TRUE)
M10.TBI$t.test_B.C
#comparing composite sample 1 with composite sample 15
M15.TBI <- TBI(T1,M15_T2, method = "%diff", nperm = 999, test.t.perm = TRUE)
M15.TBI$t.test_B.C
#comparing composite sample 1 with composite sample 20
M20.TBI <- TBI(T1,M20_T2, method = "%diff", nperm = 999, test.t.perm = TRUE)
M20.TBI$t.test_B.C

PT1 <- as.data.frame(df[c(36,42,48),])
P3_T2 <- as.data.frame(df[c(37,43,49),])
P6_T2 <-  as.data.frame(df[c(39,44,50),])
P10_T2 <- as.data.frame(df[c(40,45,51),])
P15_T2 <- as.data.frame(df[c(41,46,52),])
P20_T2 <- as.data.frame(df[c(42,47,53),])

#Pinus comparing composite sample 1 with composite sample 3
P_3.TBI <- TBI(ST1, S3_T2, method = "%diff", nperm = 999, test.t.perm = TRUE)
P_3.tbi <- as.data.frame(sp1_3.TBI$BCD.mat)

#Pinus comparing composite sample 1 with composite sample 6
P_6.TBI <- TBI(ST1,S6_T2, method = "%diff", nperm = 999, test.t.perm = TRUE)
P_6.tbi <- as.data.frame(sp1_6.TBI$BCD.mat)
#Pinus comparing composite sample 1 with composite sample 10
P_10.TBI <- TBI(ST1,S10_T2, method = "%diff", nperm = 999, test.t.perm = TRUE)
P_10.TBI$t.test_B.C
#Pinus comparing composite sample 1 with composite sample 15
P_15.TBI <- TBI(ST1,S15_T2, method = "%diff", nperm = 999, test.t.perm = TRUE)
P_15.TBI$t.test_B.C

P_20.TBI <- TBI(ST1, S20_T2, method = "%diff", nperm = 999, test.t.perm = TRUE)
P_20.TBI$t.test_B.C

#Quercus
QT1 <- as.data.frame(df[c(19,25,31),])
Q3_T2 <- as.data.frame(df[c(20,26,32), ])
Q6_T2 <-  as.data.frame(df[c(21,33), ])
QT1.1 <-df[c(1,7), ]#otherwise is not possible to compare it with Q10
Q10_T2 <- as.data.frame(df[c(22,28), ])
Q15_T2 <-  as.data.frame(df[c(23,29,34), ])
Q20_T2 <-  as.data.frame(df[c(24,30,35),])
#tree sp2 comparing composite sample 1 with composite sample 3
Q_3.TBI <- TBI(Q=RT1,R3_T2, method = "%diff", nperm = 999, test.t.perm = TRUE)
Q_3.tbi <- as.data.frame(sp2_3.TBI$BCD.mat)

#comparing composite sample 1 with composite sample 6
Q_6.TBI <- TBI(RT1,R6_T2, method = "%diff", nperm = 999, test.t.perm = TRUE)

#comparing composite sample 1 with composite sample 10
Q_10.TBI <- TBI(RT1,R10_T2, method = "%diff", nperm = 999, test.t.perm = TRUE)

#comparing composite sample 1 with composite sample 15
Q_15.TBI <- TBI(RT1,R15_T2, method = "%diff", nperm = 999, test.t.perm = TRUE)
Q_15.tbi <- sp2_15.TBI$BCD.mat

#comparing composite sample 1 with composite sample 20
Q_20.TBI <- TBI(TT1,R20_T2, method = "%diff", nperm = 999, test.t.perm = TRUE)



