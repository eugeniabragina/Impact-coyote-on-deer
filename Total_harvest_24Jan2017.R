# This code analyses the impact of coyote arrival on white-tailed deer 
# harvest in the Eastern US in 1981-2010 
# You will need 'Total_harvest.csv to run this code
library(MuMIn)
library(lmerTest)
library(nlme)
setwd('')

periods = c("1981", '1985', '1985', '1990', '1990', '1995', '1996', '2000', 
            '2000', '2005', '2005', '2010', '2011', '2014')
order = c(7,7,6,6,5,5,4,4,3,3,2,2,1)
i=1
for (i in c(1,3,5,7,9,11,13)) {
  print(c('!!!!!! NEW TIME PERIOD !!!!!!!!!'))
  print(c(periods[i], periods[i+1]))
  cd = read.table('Total_harvest.csv', sep = ',', header = T, check.names = F)
  cd=cd[2-dim(cd)[1],]
  cd$lam = (cd[,periods[i+1]]/cd[,periods[i]])^(1/(as.numeric(periods[i+1])-as.numeric(periods[i])))
  cd = cd[!is.na(cd$lam) & log(cd$lam) != -Inf & log(cd$lam) != Inf,] #get rid of NAs, Inf and -Inf
  cd$pop.growth =cd[,i+order[i] + 9]/cd[,i+order[i] + 8]
  
  #for (j in 1:length(cd$coyote_year_recl)) {                                                    
  #  if (is.na(as.numeric(periods[i+1]) - cd$coyote_year_recl[j])) {
  #    cd$year.since.coyotes[j] = NA}
  #  else if ((as.numeric(periods[i+1]) - cd$coyote_year_recl[j]) < 0) {
  #    cd$year.since.coyotes[j] = 0 }
  #  else
  #    cd$year.since.coyotes[j] = as.numeric(periods[i+1]) - cd$coyote_year_recl[j] }
  #cd = cd[!is.na(cd$year.since.coyote),]
  if (i<7) {
    m = lme(log(lam) ~ coyote_map + human.pop.dens.2010 + BIO12_median + BIO15_median + 
            BIO17_median + snow_mean + MIXED_FORE.1992 +	tree.canopy_mean.2011 +
            devel.21.22.2011 + NPPmedian	+ NPPstdev, random = ~1|state, 
            correlation = corExp(form =~xcoord_centroid+ycoord_centroid), data = cd, na.action = "na.fail")
    # print(summary(m))
  } else 
    m = lme(log(lam) ~ coyote_map + human.pop.dens.2010 + BIO12_median + BIO15_median + 
            BIO17_median + snow_mean + mixed.f.2011 +	tree.canopy_mean.2011 +
            devel.21.22.2011 + NPPmedian	+ NPPstdev, random = ~1|state, 
            correlation = corExp(form =~xcoord_centroid+ycoord_centroid), data = cd, na.action = "na.fail")
  print(summary(m))
  #dredge_m <- dredge(m)
  #print(dredge_m[1:10,])
  print('THE PERIOD IS OVER')
  print('/////////////////////////////////////////////////////////////////////////////////////////////')
  print('')
}
############# The best models#################
##############################################

#################################### 1981-1985
#The best model
m1 = lme(log(lam) ~ coyote_map + human.pop.dens.2010 + BIO12_median + BIO15_median + 
           BIO17_median + snow_mean + MIXED_FORE.1992 +	tree.canopy_mean.2011 +
           devel.21.22.2011 + NPPmedian	+ NPPstdev, random = ~1|state, data = cd, na.action = "na.fail")
r.squaredGLMM(m1) # conditional Rsq, see Nakagawa and Schielzeth 2013
m2 = lme(log(lam) ~ coyote_map + human.pop.dens.2010 + BIO12_median + BIO15_median + 
           BIO17_median + snow_mean + MIXED_FORE.1992 +	tree.canopy_mean.2011 +
           devel.21.22.2011 + NPPmedian	+ NPPstdev, 
         correlation = corExp(form =~xcoord_centroid+ycoord_centroid),
         random = ~1|state, data = cd, na.action = "na.fail")
r.squaredGLMM(m2)
m3 = lme(log(lam) ~ 1, random = ~1|state, correlation = 
  corExp(form =~xcoord_centroid+ycoord_centroid), data = cd, na.action = "na.fail")
r.squaredGLMM(m3)
#R2m        R2c
#[1,]   0 0.08874621
m1_1 = lme(log(lam) ~ snow_mean + NPPstdev, random = ~1|state, data = cd, na.action = "na.fail")
r.squaredGLMM(m1_1)

############################ end of 1981-1985

################################### 1985-1990
m1 = lme(log(lam) ~ 1, random = ~1|state, correlation = corExp(form =~xcoord_centroid+ycoord_centroid), 
         data = cd, na.action = "na.fail")
m2 = lme(log(lam) ~ pop.growth + 1, random = ~1|state, correlation = 
           corExp(form =~xcoord_centroid+ycoord_centroid), data = cd, na.action = "na.fail")
m3 = lme(log(lam) ~ WOODY_WETL.1992 + 1, random = ~1|state, correlation = 
           corExp(form =~xcoord_centroid+ycoord_centroid), data = cd, na.action = "na.fail")
AIC(m1, m2, m3)
# df       AIC
# 
# go with m
summary(m1)
# Fixed effects: log(lam) ~ 
# Value  Std.Error  DF   t-value p-value
############################ end of 1987-1990


################################### 1990-1995
m = lme(log(lam) ~ year.since.coyotes + pop.den + pop.growth + BIO5_median + BIO12_median + BIO15_median +  
          OPEN_WATER.1992 + DECIDUOUS_.1992 + EVERGREEN.1992 + grass.pastures.1992 + WOODY_WETL.1992 + HERBACEOUS_WETL.1992 + 
          LOW_INTENS.1992  + tree.canopy_mean.2011, random = ~1|state, correlation = corExp(form =~xcoord+ycoord), data = cd, na.action = "na.fail")

m1 = lme(log(lam) ~ pop.growth + LOW_INTENS.1992, random = ~1|state, correlation = corExp(form =~xcoord+ycoord), data = cd, na.action = "na.fail")
m2 = lme(log(lam) ~ pop.growth + LOW_INTENS.1992 + BIO15_median, random = ~1|state, correlation = corExp(form =~xcoord+ycoord), data = cd, na.action = "na.fail")
m3 = lme(log(lam) ~ pop.den + pop.growth, random = ~1|state, correlation = corExp(form =~xcoord_centroid+ycoord_centroid), data = cd, na.action = "na.fail")

AIC(m1, m2, m3)
# df       AIC
# m1  6 -442.6181
# m2  7 -441.1669
# m3  
# go with m1
summary(m1)
# Fixed effects: log(lam) ~  
# Value  Std.Error  DF   t-value p-value
############################ end of 1990-1995


################################### 1996-2000
m1 = lme(log(lam) ~ 1, random = ~1|state, correlation = corExp(form =~xcoord+ycoord), data = cd, na.action = "na.fail")
m2 = lme(log(lam) ~ pop.growth, random = ~1|state, correlation = corExp(form =~xcoord+ycoord), data = cd, na.action = "na.fail")
m3 = lme(log(lam) ~ devel.21.22.2001, random = ~1|state, correlation = corExp(form =~xcoord+ycoord), data = cd, na.action = "na.fail")
AIC(m1, m2, m3)
# df       AIC

# go with m1
summary(m1)
# Fixed effects: log(lam) ~ 
# Value  Std.Error  DF   t-value p-value
############################ end of 1996-2000


################################### 2000-2005
m1 = lme(log(lam) ~ pop.growth + 1, random = ~1|state, correlation = corExp(form =~xcoord+ycoord), data = cd, na.action = "na.fail")
m2 = lme(log(lam) ~ 1, random = ~1|state, correlation = corExp(form =~xcoord+ycoord), data = cd, na.action = "na.fail")
m3 = lme(log(lam) ~ pop.den + 1, random = ~1|state, correlation = corExp(form =~xcoord+ycoord), data = cd, na.action = "na.fail")
AIC(m1, m2, m3)
# df       AIC

# go with m1
summary(m1)
# Fixed effects: log(lam) ~  
# Value  Std.Error  DF   t-value p-value
############################ end of 1996-2000


################################### 2005-2010
m1 = lme(log(lam) ~ evergreen.2006 + 1, random = ~1|state, correlation = corExp(form =~xcoord+ycoord), data = cd, na.action = "na.fail")
m2 = lme(log(lam) ~ evergreen.2006 + pop.growth + 1, random = ~1|state, correlation = corExp(form =~xcoord+ycoord), data = cd, na.action = "na.fail")
m3 = lme(log(lam) ~ 1, random = ~1|state, correlation = corExp(form =~xcoord+ycoord), data = cd, na.action = "na.fail")
AIC(m1, m2, m3)
# df       AIC

# go with m1
summary(m1)
# Fixed effects: log(lam) ~ 
# Value  Std.Error  DF   t-value p-value
############################ end of 2005-2010

################################### 2011-2014 
m1 = lme(log(lam) ~ 1, random = ~1|state, correlation = corExp(form =~xcoord+ycoord), data = cd, na.action = "na.fail")
m2 = lme(log(lam) ~ pop.growth + 1, random = ~1|state, correlation = corExp(form =~xcoord+ycoord), data = cd, na.action = "na.fail")
m3 = lme(log(lam) ~ herbaceous.wetlands.2011 + 1, random = ~1|state, correlation = corExp(form =~xcoord+ycoord), data = cd, na.action = "na.fail")
AIC(m1, m2, m3)

# df       AIC
# m1  4 -599.6673
# m2  5 -596.0454
# m3  5 -592.2950
# go with m1
summary(m1)
############################ end of 2011-2014

############################################################################
########## Comparison of models with and without spatial correlation #######
# i=1 i.e. first period 1981-1985
i=3
m.spat = lme(log(lam) ~ coyote_James  + pop.growth + Pop.dens.2010 + BIO12_median + BIO15_median + 
               BIO17_median + snow_mean + MIXED_FORE.1992 +	shrub.2011 +	grass.pastures.1992	+ tree.canopy_mean.2011 +
               decidious.f.2011	+ evergreen.2011	+ crops.2011	+ woody.wetl.2011 +	herbaceous.wetlands.2011	+ 
               devel.21.22.2011 + NPPmedian	+ NPPstdev, random = ~1|state, 
             correlation = corExp(form =~xcoord+ycoord), data = cd, na.action = "na.fail")

m.nospat = lme(log(lam) ~ coyote_James  + pop.growth + Pop.dens.2010 + BIO12_median + BIO15_median + 
                 BIO17_median + snow_mean + MIXED_FORE.1992 +	shrub.2011 +	grass.pastures.1992	+ tree.canopy_mean.2011 +
                 decidious.f.2011	+ evergreen.2011	+ crops.2011	+ woody.wetl.2011 +	herbaceous.wetlands.2011	+ 
                 devel.21.22.2011 + NPPmedian	+ NPPstdev, random = ~1|state, data = cd, na.action = "na.fail")

m.spat.selected = lme(log(lam) ~ coyote_year_recl + BIO17_median, random = ~1|state, 
                      correlation = corExp(form =~xcoord+ycoord), data = cd, na.action = "na.fail")

AIC(m.spat.selected)
summary(m.spat)
dredge_m <- dredge(m.spat.selected)
dredge_m[1:20,]

