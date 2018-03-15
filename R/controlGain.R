#' controlGain function estimate phenotypic, agronomic, and genetic trends using the control pop method
#'
#' @import ggplot2
#' @import grid
#' @import lme4
#' @import lsmeans
#' @param filenm a csv file name string
#' @param label a character with the plot label appendage
#' @param tunit a character naming the units
#' @return list of objects containing the results
#' @export
#'

controlGain<- function(dat, label='', tunit='units'){
  #Check the column names
  if(any(colnames(dat)!= c("gid", "blue","se","season_number")))
    stop("data should contain gid, blue, se, and season_number columns")
  dat$gid<- as.character(dat$gid)
  dat$season_number<- as.numeric(dat$season_number)
  dat$blue<- as.numeric(dat$blue)
  dat$se<- as.numeric(dat$se)

  #check the dataset for duplicate gids within a season
  id_sea<- paste(dat$gid, dat$season_number)
  if(!length(unique(id_sea))==length(id_sea))
    stop("Error: There should be only one value per genotype per season number")

  #identify which are the potential checks
  tb<- table(dat$gid)
  ckcand<- tb[which(tb>=10)] #checks should be in at least 10 years
  if(length(ckcand)==0)
    stop("Error: There are no checks in the dataset that have been used for 10+ years")
  ckcand<- names(ckcand)

  #discard earlier seasons where the checks are not sufficient
  dck<-dat[which(dat$gid %in% ckcand),]
  octab<- reshape::cast(dck, gid~season_number, value='blue')
  octab2<- octab[,-1]
  tally<- apply(octab2, 2, function(x)length(na.omit(x)))
  colrm<- which(tally<3)
  if(length(colrm)!=0){
    start_num<- max(which(tally<3))+1
  }else{
    start_num<- 1
  }

  if(start_num>ncol(octab2))
    stop("Error: There needs to be at least common 3 checks used for the past 10 seasons (minimum)")
  sea_sel<- colnames(octab2)[start_num:ncol(octab2)]
  nsea<- length(sea_sel)

  if(nsea<10)
    stop("Error: There needs to be at least common 3 checks used for the past 10 seasons (minimum)")

  #select the seasons with sufficent checks for the analysis
  sea_sel<- as.numeric(sea_sel)
  dat<- dat[which(dat$season_number %in% sea_sel),]

  #identify the checks
  tb<- table(dat$gid)
  tbsel<- tb[which(tb==length(unique(dat$season_number)))]
  cks<- names(tbsel)

  #add the check indicator column
  dat<- data.frame(dat, Population='Selected', stringsAsFactors = FALSE)
  dat[which(dat$gid %in% cks),'Population']<- 'Control'
  dat$Population<- factor(dat$Population, levels=c('Control', 'Selected'))

  #start time at 0
  mnsea<- min(dat$season_number)
  dat$season_number<- dat$season_number-mnsea

  #fit model to get points for genetic value (stage 1)
  dat$season_number<- as.character(dat$season_number)
  wtMod0<- 1/(dat$se/max(dat$se))
  mod0<- lm(blue~ Population+season_number+Population:season_number, weights=wtMod0, data=dat)
  lsm<- lsmeans::lsmeans(mod0, specs='Population', by='season_number', cov.reduce=FALSE)
  lsm<- lsmeans::contrast(lsm, method='trt.vs.ctrl')
  pts<- as.data.frame(summary(lsm))

  #genetic trend model and fitted values (stage 2)
  pts$season_number<- as.numeric(as.character(pts$season_number))
  wt<- 1/(pts$SE/max(pts$SE))
  mdG<-lm(estimate~season_number, weights=wt, data=pts)
  lsm<- lsmeans::lsmeans(mdG, specs='season_number', cov.reduce=FALSE)
  genEst<- as.data.frame(summary(lsm))
  pts$season_number<-pts$season_number+mnsea #Gen-val tab

  #get the points for the other plots (stage 1)
  lsm<- lsmeans::lsmeans(mod0, specs='season_number', by='Population', cov.reduce=FALSE)
  ptsPop<- as.data.frame(summary(lsm))

  #Phenotypic trends (stage 2)
  ptsPop$season_number<- as.numeric(as.character(ptsPop$season_number))
  wt<- 1/(ptsPop$SE/max(ptsPop$SE))
  mdPS<-lm(lsmean~season_number+Population+Population:season_number, data=ptsPop, weights=wt)
  lsm<- lsmeans::lsmeans(mdPS, specs='Population', by='season_number' ,cov.reduce=FALSE)
  fitted<- as.data.frame(summary(lsm))
  ptsPop$season_number<- ptsPop$season_number+mnsea
  fitted$season_number<- fitted$season_number+mnsea

  #create the phenotypic trend plot
  dat$season_number<- as.numeric(dat$season_number)+mnsea
  p1<- ggplot(data=ptsPop, aes(x=season_number, y=lsmean, group=Population)) +
    theme_minimal()+
    geom_point(size=2,stroke=1,aes(shape=Population, color=Population))+
    scale_color_manual(values=c('slategray4', 'darkorange'))+
    scale_shape_manual(values=c(0, 1))+
    geom_errorbar(size=1, aes(ymin=lsmean-SE, ymax=lsmean+SE, width=0.2, color=Population))+
    geom_line(data=fitted, size=0.5, aes(linetype=Population,  color=Population))+
    scale_linetype_manual(values=c("longdash", "twodash"))+
    geom_ribbon(data=fitted, aes(ymin=lower.CL, ymax=upper.CL, fill=Population), alpha=0.2)+
    scale_fill_manual(values=c('slategray4', 'darkorange'), name="fill")+
    labs(y = paste("Phenotypic value in", tunit, sep=" "), x='Season number')+
    ggtitle(paste("Phenotypic trends", label, sep=""))+
    theme(legend.position="top",plot.title = element_text(hjust = 0.5))+
    guides(fill = "none")

  #create the genetic trend plot
  #geom_point(color='slategray4', size=2, shape=0)+geom_line(color='black')+
  genEst$season_number<- genEst$season_number+mnsea
  colnames(pts)[3]<- 'lsmean'
  p2<- ggplot(pts,aes(x=season_number, y=lsmean)) +
    geom_point(shape=1, size=2,stroke=1, color='slategray4', data=pts, aes(x=season_number,y=lsmean))+
    geom_line(data=genEst, color='slategray4', linetype='longdash')+
    geom_errorbar(color='slategray4', data=pts, size=1, aes(ymin=lsmean-SE, ymax=lsmean+SE, width=0.2))+
    theme_minimal()+
    geom_ribbon(data=genEst, aes(ymin=lsmean-SE, ymax=lsmean+SE), fill='slategray4', alpha=0.2)+
    labs(y = paste("Predicted average genetic value in", tunit,sep=" "), x='Season number')+
    ggtitle(paste("Predicted genetic trend", label, sep=""))+
    theme(plot.title = element_text(hjust = 0.5))

  #get the model information
  genTab<- as.data.frame(summary(mdG)$coefficients)
  psTab<- as.data.frame(summary(mdPS)$coefficients)

  cfnms<-c(paste('Control baseline in', tunit, sep=" "),
                paste('Control trend in', tunit, "per season number", sep=" "),
                paste('Selected baseline in', tunit, sep=" "),
                paste('Selected-control trend in', tunit, "per season number", sep=" "),
                paste('Predicted genetic value baseline in', tunit, sep=" "),
                paste('Predicted genetic trend in', tunit, "per season number", sep=" "))
  rslts<- data.frame(Model=c(rep('P',4),
             rep('G',2)),cfnms, rbind(psTab, genTab))
  row.names(rslts)<- c(1:nrow(rslts))
  colnames(rslts)<- c("Parameter", "Model", "Estimate", "Standard error",
    't-value', 'p-value')
  rslts<-format(rslts, digits=2)


  #make anova table into a dataframe
  avP<- as.data.frame(anova(mdPS))
  avG<- as.data.frame(anova(mdG))
  varnms<- c('Season number', 'Population', 'Season number x population',
             'Residuals', 'Season number', 'Residuals')

  av<- data.frame(Model=c(rep('P',4),rep('G',2)), varnms, rbind(avP, avG))
  row.names(av)<- c(1:nrow(av))
  av<- format(av, digits=2)
  av[which(trimws(av[,6])=='NA'),6]<-""
  av[which(trimws(av[,7])=='NA'),7]<-""

  which(is.na(trimws(av[,6])))
  colnames(av)<- c('Model', 'Variable', 'df',
                   'Sum of squares', 'Mean square',
                   'F-value', 'p-value')

  out<- list(p1=p1, p2=p2, rslts=rslts, av=av)
  return(out)
}


