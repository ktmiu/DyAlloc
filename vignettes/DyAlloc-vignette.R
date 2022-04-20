## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(DyAlloc)

## ---- example without site----------------------------------------------------
#generate a data frame of 50 subjects in a 4 treatment group study
df1 <- data.frame("gender" = sample(c("female", "male"),50, TRUE),
                 "age" = sample(c("0-35", ">35"),50, TRUE),
                 "race" = sample(c("black", "white"),50,TRUE),
                 "trauma severity"= sample(c("good", "bad"),50,TRUE),
                 "treatment"=sample(c("1","2","3","4"),50, TRUE,c(1,1,5,10)),
                  stringsAsFactors = TRUE)
#covariates of the study
covars <-c("gender","age","race","trauma severity")
#new subject's profile
newsub <-c("male","0-35","black","bad")
#weight of the overall study, the stratum, and the 4 factors 
weight <-c(1,1,rep(1,4))
#print out the final weighted imbalance scores
imb <-WeiImb(df1,4,covars,weight,newsub,FALSE)
imb
#probabilities for the new subject to be assigned to the each of the treatment group 
Trt_Prob(imb,0.2)

#Using CAR function to get the final treatment assignment
DyAlloc(df1,4,covars,weight,newsub,FALSE)

## ---- example with site-------------------------------------------------------
df <- data.frame("gender" = sample(c("female", "male"),50, TRUE),
                 "age" = sample(c("0-35", ">35"), 50, TRUE),
                 "race" = sample(c("black", "white"),50,TRUE),
                 "trauma severity"= sample(c("good", "bad"),50,TRUE),
                 "treatment"=sample(c("1","2","3","4"),50, TRUE),
                 "site"=sample(c("1","2"),50, TRUE),
                 stringsAsFactors = TRUE)
covars <-c("gender","age","race","trauma severity")
#new subject's profile with the first element being site
newsub<-c("2","male","0-35","black","bad")
#weight of the overall study, the stratum, the site, and the 4 factors 
weight <-c(1,1,1,rep(1,4))
imb<-WeiImb(df,4,covars,weight,newsub,TRUE)
imb
Trt_Prob(imb,0.2)
DyAlloc(df,4,covars,weight,newsub,TRUE)

## ---- include=FALSE-----------------------------------------------------------
set.seed(12345)
mydf <- data.frame("gender" = sample(c("female", "male"),200, TRUE,c(.2,.8)),
                   "age" = sample(c("0-35", ">35"),200, TRUE),
                   "race" = sample(c("black", "white"),200,TRUE,c(0.7,0.3)),
                   "trauma severity"= sample(c("good", "bad"),200,TRUE), 
                   stringsAsFactors = TRUE)
covars <-c("gender","age","race","trauma severity")
weight <-c(1,1,rep(1,4))
treatment <- vector()
newdf <- data.frame("gender"=NA,"age"=NA,"race"=NA,"trauma severity"=NA,"treatment"=NA)
for (i in 1:200){
  newsub<-c(as.character(mydf[i,1]),as.character(mydf[i,2]),as.character(mydf[i,3]),
            as.character(mydf[i,4]))
  newtreatment <- DyAlloc(newdf,4,covars,weight,newsub,FALSE)
  treatment <- append(treatment,newtreatment)
  newdf <- cbind(mydf[1:i,],treatment)
}
newdf
newdf$srs = sample(c("1", "2", "3", "4"), 200, replace = TRUE)

## ---- echo=FALSE--------------------------------------------------------------
#number of female and male in the newdf
table(newdf$gender)
#proportion of each gender in each treatment 
prop.table(table(newdf$gender,newdf$treatment),1)
#proportion of each gender in each treatment 
prop.table(table(newdf$gender,newdf$srs),1)

## ---- echo=FALSE--------------------------------------------------------------
#disproportional race example
table(newdf$race)
prop.table(table(newdf$race,newdf$treatment),1)
prop.table(table(newdf$race,newdf$srs),1)

## ---- echo=FALSE--------------------------------------------------------------
#disproportional race example
table(newdf$race, newdf$gender)

ctab = rbind(prop.table(table(newdf$race[newdf$gender=="female"],newdf$treatment[newdf$gender=="female"]),1), prop.table(table(newdf$race[newdf$gender=="male"],newdf$treatment[newdf$gender=="male"]),1))
dimnames(ctab)[[1]] = c("black female", "white female", "black male", "white male")
ctab


ctab.srs = rbind(prop.table(table(newdf$race[newdf$gender=="female"],newdf$srs[newdf$gender=="female"]),1), prop.table(table(newdf$race[newdf$gender=="male"],newdf$srs[newdf$gender=="male"]),1))
dimnames(ctab.srs)[[1]] = c("black female", "white female", "black male", "white male")
ctab.srs

