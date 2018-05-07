#' controlGain function estimate phenotypic, agronomic, and genetic trends using the control pop method
#'
#' @import ggplot2
#' @import grid
#' @import lme4
#' @param dat a dataframe with the multi-year data
#' @param label a character with the plot label appendage
#' @param tunit a character naming the units
#' @return list of objects containing the results
#' @export

controlGain2<- function(dat, label='', tunit='units',minNumb=7,minNumbCk=1){
  if('phenoValue' %in% colnames(dat)){
    dat<- prepData(dat)
  }
  #Check the column names
  if(any(colnames(dat)!= c("gid", "blue","se","year")))
    stop("data should contain gid, blue, se, and year columns")
  dat$gid<- as.character(dat$gid)
  dat$year<- as.numeric(dat$year)
  dat$blue<- as.numeric(dat$blue)
  dat$se<- as.numeric(dat$se)

  #check the dataset for duplicate gids within a year
  id_sea<- paste(dat$gid, dat$year)
  if(!length(unique(id_sea))==length(id_sea))
    stop("Error: There should be only one value per genotype per year")

  #identify which are the potential checks
  tb<- table(dat$gid)
  ckcand<- tb[which(tb>=minNumb)] #checks should be in at least 10 years
  if(length(ckcand)==0)
    stop("Error: There are no checks in the dataset that have been used for the minimum number of years")
  ckcand<- names(ckcand)

  #discard earlier years where the checks are not sufficient
  dck<-dat[which(dat$gid %in% ckcand),]
  octab<- reshape::cast(dck, gid~year, value='blue')
  octab2<- octab[,-1]
  tally<- apply(octab2, 2, function(x)length(na.omit(x)))
  colrm<- which(tally<minNumbCk)
  if(length(colrm)!=0){
    start_num<- max(which(tally<minNumbCk))+1
  }else{
    start_num<- 1
  }

  if(start_num>ncol(octab2))
    stop("Error: There needs to be at least common 3 checks used for the past 10 years by default")
  sea_sel<- colnames(octab2)[start_num:ncol(octab2)]
  nsea<- length(sea_sel)

  if(nsea<minNumb)
    stop("Error: There needs to be at least common 3 checks used for the past 10 years by default")

  #select the years with sufficent checks for the analysis
  sea_sel<- as.numeric(sea_sel)
  dat<- dat[which(dat$year %in% sea_sel),]

  #identify the checks
  tb<- table(dat$gid)
  tbsel<- tb[which(tb==length(unique(dat$year)))]
  cks<- names(tbsel)

  #add the check indicator column
  dat<- data.frame(dat, Population='Selected', stringsAsFactors = FALSE)
  dat[which(dat$gid %in% cks),'Population']<- 'Control'
  dat$Population<- factor(dat$Population, levels=c('Control', 'Selected'))

  #keep earliest occurance of varieties only
  dat<- dat[order(dat$year),] #order by year
  selgids<- unique(dat$gid[which(dat$Population=='Selected')])
  ixkp<- c(which(dat$Population=='Control'), match(selgids, dat$gid))
  dat<- dat[ixkp, ]

  #start time at 0
  mnsea<- min(dat$year)
  dat$year<- dat$year-mnsea

  #fit model to get points for genetic value (stage 1)
  dat$year<- as.character(dat$year)
  wtMod0<- 1/dat$se^2
  year_num<- as.numeric(dat$year)
  dat<- data.frame(dat, year_num)

  mod0<- lmer(blue~ Population+year_num+Population:year_num+(1|year)+(1|Population:year),
              weights=wtMod0, data=dat)
  modnull<- lmer(blue~ Population+year_num+(1|year)+
                   (1|Population:year),
                 weights=wtMod0, data=dat)
  av<- anova(modnull, mod0)
  pval<- as.data.frame(av)[2,8]
  est<- summary(mod0)$coefficients[4,1]
  se<- summary(mod0)$coefficients[4,2]


  out<- list(Estimate=est, SE=se, Pvalue=pval)
  return(out)
}


