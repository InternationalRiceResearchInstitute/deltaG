#' controlGain function estimate phenotypic, agronomic, and genetic trends using the control pop method
#'
#' @import ggplot2
#' @import grid
#' @import emmeans
#' @param dat a dataframe with the multi-year data
#' @param label a character with the plot label appendage
#' @param tunit a character naming the units
#' @return list of objects containing the results
#' @export

controlGain<- function(dat, label='', tunit='units',minNumb=7,minNumbCk=1){
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
  mod0<- lm(blue~ Population+year+Population:year, weights=wtMod0, data=dat)
  lsm<- emmeans::emmeans(mod0, specs='Population', by='year', cov.reduce=FALSE)
  lsm<- emmeans::contrast(lsm, method='trt.vs.ctrl')
  pts<- as.data.frame(summary(lsm))

  #genetic trend model and fitted values (stage 2)
  pts$year<- as.numeric(as.character(pts$year))
  wt<- 1/pts$SE^2
  mdG<-lm(estimate~year, weights=wt, data=pts)
  lsm<- emmeans::emmeans(mdG, specs='year', cov.reduce=FALSE)
  genEst<- as.data.frame(summary(lsm))
  pts$year<-pts$year+mnsea #Gen-val tab

  #get the points for the other plots (stage 1)
  lsm<- emmeans::emmeans(mod0, specs='year', by='Population', cov.reduce=FALSE)
  ptsPop<- as.data.frame(summary(lsm))

  #Phenotypic trends (stage 2)
  ptsPop$year<- as.numeric(as.character(ptsPop$year))

  wt<- 1/ptsPop$SE^2
  mdPS<-lm(emmean~year+Population+Population:year, data=ptsPop, weights=wt)
  lsm<- emmeans::emmeans(mdPS, specs='Population', by='year' ,cov.reduce=FALSE)
  fitted<- as.data.frame(summary(lsm))
  ptsPop$year<- ptsPop$year+mnsea
  fitted$year<- fitted$year+mnsea

  #create the phenotypic trend plot
  dat$year<- as.numeric(dat$year)+mnsea
  p1<- ggplot(data=ptsPop, aes(x=year, y=emmean, group=Population)) +
    theme_minimal()+
    geom_point(size=2,stroke=1,aes(shape=Population, color=Population))+
    scale_color_manual(values=c('slategray4', 'darkorange'))+
    scale_shape_manual(values=c(0, 1))+
    geom_errorbar(size=1, aes(ymin=emmean-SE, ymax=emmean+SE, width=0.2, color=Population))+
    geom_line(data=fitted, size=0.5, aes(linetype=Population,  color=Population))+
    scale_linetype_manual(values=c("longdash", "twodash"))+
    geom_ribbon(data=fitted, aes(ymin=lower.CL, ymax=upper.CL, fill=Population), alpha=0.2)+
    scale_fill_manual(values=c('slategray4', 'darkorange'), name="fill")+
    labs(y = paste("Phenotypic value in", tunit, sep=" "), x='Year')+
    ggtitle(paste("Phenotypic trends", label, sep=""))+
    theme(legend.position="top",plot.title = element_text(hjust = 0.5))+
    guides(fill = "none")

  #create the genetic trend plot
  #geom_point(color='slategray4', size=2, shape=0)+geom_line(color='black')+
  genEst$year<- genEst$year+mnsea
  colnames(pts)[3]<- 'emmean'
  p2<- ggplot(pts,aes(x=year, y=emmean)) +
    geom_point(shape=1, size=2,stroke=1, color='slategray4', data=pts, aes(x=year,y=emmean))+
    geom_line(data=genEst, color='slategray4', linetype='longdash')+
    geom_errorbar(color='slategray4', data=pts, size=1, aes(ymin=emmean-SE, ymax=emmean+SE, width=0.2))+
    theme_minimal()+
    geom_ribbon(data=genEst, aes(ymin=emmean-SE, ymax=emmean+SE), fill='slategray4', alpha=0.2)+
    labs(y = paste("Predicted average genetic value in", tunit,sep=" "), x='Year')+
    ggtitle(paste("Predicted genetic trend", label, sep=""))+
    theme(plot.title = element_text(hjust = 0.5))

  #get the model information
  genTab<- as.data.frame(summary(mdG)$coefficients)
  psTab<- as.data.frame(summary(mdPS)$coefficients)

  cfnms<-c(paste('Control baseline in', tunit, sep=" "),
                paste('Control trend in', tunit, "per year", sep=" "),
                paste('Selected baseline in', tunit, sep=" "),
                paste('Selected-control trend in', tunit, "per year", sep=" "),
                paste('Predicted genetic value baseline in', tunit, sep=" "),
                paste('Predicted genetic trend in', tunit, "per year", sep=" "))
  rslts<- data.frame(Model=c(rep('P',4),
             rep('G',2)),cfnms, rbind(psTab, genTab))
  row.names(rslts)<- c(1:nrow(rslts))
  colnames(rslts)<- c( "Model","Parameter", "Estimate", "Standard error",
    't-value', 'p-value')
  rslts<-format(rslts, digits=2)


  #make anova table into a dataframe
  avP<- as.data.frame(anova(mdPS))
  avG<- as.data.frame(anova(mdG))
  varnms<- c('year', 'Population', 'year x population',
             'Residuals', 'year', 'Residuals')

  av<- data.frame(Model=c(rep('P',4),rep('G',2)), varnms, rbind(avP, avG))
  row.names(av)<- c(1:nrow(av))
  av<- format(av, digits=2)
  av[which(trimws(av[,6])=='NA'),6]<-""
  av[which(trimws(av[,7])=='NA'),7]<-""

  which(is.na(trimws(av[,6])))
  colnames(av)<- c('Model', 'Variable', 'df',
                   'Sum of squares', 'Mean square',
                   'F-value', 'p-value')

  out<- list(p1=p1, p2=p2, rslts=rslts, av=av, Gmod=mdG, Gvals=pts)
  return(out)
}


