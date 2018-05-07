#' prepData function estimate phenotypic, agronomic, and genetic trends using the control pop method
#'
#' @import emmeans
#' @import nlme
#' @param dat
#' @return tab
#' @export

prepData<- function(dat){
#get per year adjusted means
uyr<- unique(dat$year)
if(length(is.na(dat$phenoValue))>0){
  dat<- dat[which(!is.na(dat$phenoValue)),]
}
for(i in 1:length(uyr)){
  sub<- droplevels.data.frame(dat[which(dat$year==uyr[i]),])
  sub$site<- as.character(sub$site)
  sub$rep<- as.character(sub$rep)
  sub$lineID<- as.character(sub$lineID)
  sub$lineID_site<- paste(sub$lineID, sub$site)
  mod<- lme(phenoValue~ 1 + lineID + site + rep:site, random= ~ 1|lineID_site, data=sub)
  mod2<- update(mod, weights=varIdent(form~1|site))
  lsm<- emmeans(mod2, specs='lineID')
  mns<- as.data.frame(summary(lsm))
  df<- data.frame(year=uyr[i], mns)
  if(i==1){
    singlyr<- df
  }else{
    singlyr<- rbind(singlyr, df)
  }
}
colnames(singlyr)[1:4]<- c('year', 'gid', 'blue', 'se')
singlyr<- singlyr[,c('gid', 'blue', 'se', 'year')]
return(singlyr)
}

