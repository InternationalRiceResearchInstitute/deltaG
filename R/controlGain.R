#' controlGain function estimate phenotypic, agronomic, and genetic trends using the control pop method
#'
#' @import ggplot2
#' @import grid
#' @import lme4
#' @import lsmeans
#' @param filenm a csv file name string
#' @param label a character with the plot label appendage
#' @param tunit a character naming the units
#' @param x1 lower range of the ylimit of the first plot
#' @param y1 upper range of the ylimit of the first plot
#' @param x2 lower range of the ylimit of the second plot
#' @param y2 upper range of the ylimit of the second plot
#' @return list of objects containing the results
#' @export
#'

controlGain<- function(dat, label='', tunit='units', x1=NULL, y1=NULL, x2=NULL, y2=NULL){
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
  dat<- data.frame(dat, ISCK=FALSE)
  dat[which(dat$gid %in% cks),'ISCK']<- TRUE
  dat$ISCK<- factor(dat$ISCK, levels=c('TRUE', 'FALSE'))

  #run the analysis
  mnsea<- min(dat$season_number)
  dat$season_number<- dat$season_number-mnsea
  wt<- 1/(dat$se/max(dat$se))
  dat<- data.frame(dat, seafac=as.character(dat$season_number))
  if(var(wt)>0){
    mod0<- lmer(blue~ (1|seafac)+ISCK+season_number, weight=wt, data=dat)
    mod1<- lmer(blue~ (1|seafac)+ISCK+season_number+ISCK:season_number, weight=wt, data=dat)
  }else{
    mod0<- lmer(blue~ (1|seafac)+ISCK+season_number, data=dat)
    mod1<- lmer(blue~ (1|seafac)+ISCK+season_number+ISCK:season_number, data=dat)
  }
  modcomp<- anova(mod0, mod1)
  pest<- modcomp$`Pr(>Chisq)`[2]
  cfs<- fixef(mod1)
  secoef<- summary(mod1)$coefficients[,2]

  #Get the genetic trend estimate
  RateEst<- cfs['ISCKFALSE:season_number']
  seEst<- secoef[4]

  #Get the fitted values
  lsm<- lsmeans::lsmeans(mod1, specs='season_number', by='ISCK', cov.reduce=FALSE)
  fitted<- as.data.frame(summary(lsm))
  colnames(fitted)[c(3:4)]<- c('blue', 'se')

  #Get the contrast means (predicted genetic trend)
  lsm<- lsmeans::lsmeans(mod1, specs='ISCK', by='season_number', cov.reduce=FALSE)
  lsm<- lsmeans::contrast(lsm, method='trt.vs.ctrl')
  genEst<- data.frame(summary(lsm))

  #phenotypic and agronomic trends
  dat$season_number<- dat$season_number+mnsea
  fitted$season_number<- fitted$season_number+mnsea
  dat<- data.frame(dat, group2=paste(dat$ISCK, dat$season_number))

  #Change labeling of populations
  colnames(dat)[5]<- 'Population'
  dat$Population<- as.character(dat$Population)
  dat[which(dat$Population==TRUE),'Population']<- "Control"
  dat[which(dat$Population==FALSE),'Population']<- "Selected"
  colnames(fitted)[2]<- 'Population'
  fitted$Population<- as.character(fitted$Population)
  fitted[which(fitted$Population==TRUE),'Population']<- "Control"
  fitted[which(fitted$Population==FALSE),'Population']<- "Selected"

  #create the phenotypic trend plot
  p1<- ggplot(data=dat, aes(x=season_number, y=blue, group=Population)) +
    geom_jitter(alpha=0.5, width=0.25, aes(color=Population, shape=Population))+
    stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1),
                 geom="crossbar", width=0.4, size=0.7, aes(color=Population))+
    scale_shape_manual(values=c(1, 0))+
    theme_minimal()+
    scale_color_manual(values=c('slategray4', 'darkorange'))+
    geom_line(data=fitted, size=0.5, aes(linetype=Population,  color=Population))+
    scale_linetype_manual(values=c("longdash", "twodash"))+
    geom_ribbon(data=fitted, aes(ymin=lower.CL, ymax=upper.CL, fill=Population), alpha=0.2)+
    scale_fill_manual(values=c('slategray4', 'darkorange'), name="fill")+
    labs(y = paste("Phenotypic value in", tunit, sep=" "), x='Season number')+
    ggtitle(paste("Phenotypic trends", label, sep=""))+
    theme(legend.position="top",plot.title = element_text(hjust = 0.5))+
    guides(fill = "none")
  if(!is.null(x1) & !is.null(y1))
    p1<- p1+ylim(x=x1, y=y1)


  #fit model to get points for genetic value plot
  dat$season_number<- as.character(dat$season_number)
  mod2<- lm(blue~ Population+season_number+Population:season_number, data=dat)
  lsm<- lsmeans::lsmeans(mod2, specs='Population', by='season_number', cov.reduce=FALSE)
  lsm<- lsmeans::contrast(lsm, method='trt.vs.ctrl')
  pts<- as.data.frame(summary(lsm))
  pts$season_number<- as.numeric(pts$season_number)+mnsea

  #create the genetic trend plot
  #geom_point(color='slategray4', size=2, shape=0)+geom_line(color='black')+
  genEst$season_number<- genEst$season_number+mnsea
  pts$SE<- pts$SE*2
  p2<- ggplot(genEst,aes(x=season_number,y=estimate)) +
    geom_point(shape=1, size=4.75,stroke=1, color='slategray4', data=pts, aes(x=season_number,y=estimate))+
    geom_line(color='slategray4', linetype='longdash')+
    geom_errorbar(color='slategray4', data=pts, size=1.15, aes(ymin=estimate-SE, ymax=estimate+SE, width=0.0))+
    theme_minimal()+
    geom_ribbon(aes(ymin=estimate-SE, ymax=estimate+SE), fill='slategray4', alpha=0.2)+
    labs(y = paste("Predicted average genetic value in", tunit,sep=" "), x='Season number')+
    ggtitle(paste("Predicted genetic trend", label, sep=""))+
    theme(plot.title = element_text(hjust = 0.5))

  if(!is.null(x2) & !is.null(y2))
    p2<- p2+ylim(x=x2, y=y2)

  #make results summary table
  rslts<- summary(mod1)$coefficients
  rslts<- data.frame(Parameter=c(paste('Control phenotypic baseline in', tunit, sep=" "),
    paste('Genetic value baseline in', tunit, sep=" "),
    paste('Agronomic trend,', tunit, "per season number", sep=" "),
    paste('Genetic trend,', tunit, "per season number", sep=" ")), rslts)
  row.names(rslts)<- c(1:4)
  colnames(rslts)[3:4]<- c("Standard_Error", 't-value')
  rslts<- format(rslts, digits=2)

  #make anova table into a dataframe
  av<- as.data.frame(anova(mod1))
  av<- data.frame(Parameter=c('Population',
  'Season number','Season number x population'), av)
  row.names(av)<- c(1:3)
  av<- format(av, digits=2)
  colnames(av)[c(2:5)]<- c('df', 'Sum of squares', 'Mean square', 'F-value')

  out<- list(p1=p1, p2=p2, rslts=rslts, av=av, r2)
  return(out)
}


