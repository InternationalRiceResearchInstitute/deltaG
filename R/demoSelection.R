#' demoSelection function to demonstrate phenotypic selection outcome
#' @import ggplot2
#' @param pop0min population min
#' @param pop0max population max
#' @param herit heritability in the TPE
#' @param popsize number of potential parents
#' @param numparents number of parents selected
#' @param ncycles number of cycles of selection
#' @param cycledur time required to complete 1 cycle
#' @export

demoSelection<- function(pop0min= 0, pop0max= 7,
                         herit= 0.1, popsize= 100,
                        numparents= 20, ncycles= 3,
                        cycledur=7, rnseed=99){

if(numparents>popsize){
  numparents=popsize
}

#variables
pop0mean<- mean(pop0min:pop0max)
popStd<- c(pop0max-pop0min)/8

#phenotypic variance
varP<- popStd^2

#proportion selected
p<- numparents/popsize

#selection intensity
i<- dnorm(qnorm(p))/p

#selection accuracy
selacc<- sqrt(herit)

#additive genetic variance
varA<- herit*varP

#Response per generation
Rpergen_avg<- i*sqrt(varA)*selacc

#cycle zero vec
set.seed(seed= rnseed)
cycle0 <- data.frame(Phenotypic_Value =
                       rnorm(popsize, pop0mean, sqrt(varP)))


#seed for after cycle 0
if(ncycles==1){
  newseed<- rnseed*popsize*numparents*123
}else{
  newseed<- rnseed*ncycles*herit*popsize*numparents*123
}

set.seed(seed=newseed)

#Make cycle 1 based on cycle 0
  #get response based on R= h2S
  cycle1<- cycle0
  k<-0
  totR<- 0
  Rvec<- c()
  while(k<ncycles){
    ordvec<- cycle1[order(-cycle1[,1]),]
    mnSel<- mean(ordvec[1:c(length(ordvec)*p)])
    mnTot<- mean(ordvec)
    S<- mnSel-mnTot
    Gain<- S*herit
    totR<- totR+Gain
    Rvec<- append(Rvec, totR)
    k<- k+1
    newMean<- mean(cycle1[,1])+Gain
    cycle1 <- data.frame(Phenotypic_Value =
                         rnorm(popsize, newMean, sqrt(varP)))
  }

meansVec<- c(pop0mean, Rvec+pop0mean)
dfMns<- rbind(data.frame(means=meansVec, year=c(0:ncycles)*cycledur, gain='In current program'),
        data.frame(means=(Rpergen_avg*c(0:ncycles))+pop0mean, year=c(0:ncycles)*cycledur, gain='On average'))


#population 1 mean
pop1mean<- pop0mean+totR

#Now, combine your two dataframes into one.  First make a new column in each that will be a variable to identify where they came from later.
cycle0$Population <- 'Before selection'
cycle1$Population <- 'After selection'

#and combine into your new data frame vegLengths
cycVecs <- rbind(cycle0, cycle1)
cycVecs <- cycVecs[order(-cycVecs[,1]),]
cycVecs<- data.frame(cycVecs, ranking=c(1:nrow(cycVecs)))

#get xlim
popMin<- min(cycVecs[,1])
popMax<- max(cycVecs[,1])
rg<- popMax-popMin
xmin<-popMin- rg *0.5
xmx<- popMax+ rg *0.5

#expected gain per year
Rperyear<- totR/(ncycles*cycledur)

#caption
Caption<- paste(paste(round(round(totR,3), 2), 'trait units are gained after', ncycles*cycledur, 'years. '),
    paste("This is equal to", round(Rperyear/sqrt(varA),2), 'genetic standard deviations per year'),
    paste(" or ", round(Rperyear/pop0mean *100, 2), ' percent per year in this scenario.', sep=""), sep="")


#plot title
plotTit<- paste("Gain from selection after ", ncycles*cycledur, ' years', sep="")


#make plot
plt<- ggplot2::ggplot(cycVecs, aes(Phenotypic_Value, fill = Population)) +
  geom_density(alpha = 0.3)+
  scale_x_continuous(limits = c(xmin, xmx))+
  xlab("Phenotypic value") +
  ylab("Density") +
  ggtitle(plotTit)+
  scale_fill_manual( values = c("orange","grey0"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))

plt2<- ggplot2::ggplot(cycVecs, aes(x = reorder(ranking, Phenotypic_Value),
                    y = Phenotypic_Value, fill = Population)) +
  geom_bar(stat='identity', position='identity',
           color='black', size= 4/popsize)+
  scale_fill_manual(values = c("peachpuff","grey65"))+
  xlab("Individual") +
  ylab("Phenotypic value") +
  ggtitle("Phenotypic values by population")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

rg<- max(dfMns$means)-min(dfMns$means)
#ymin<-min(dfMns$means)- rg*0.5
#ymx<- max(dfMns$means)+ rg *0.5
ymin<-pop0min
ymx<- pop0max+totR
plt3<- ggplot(data=dfMns, aes(x=year, y=means, colour=gain)) +
  geom_line(lwd=1)+
  geom_point()+
  scale_y_continuous(limits = c(ymin, ymx))+
  xlab("Year") +
  ylab("Average phenotypic value") +
  ggtitle("Average phenotypic values over time")+
  scale_colour_manual(name="Genetic gain",
                      values=c("orange", "blue"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), legend.key=element_blank(),
        axis.line = element_line(colour = "black"))


#variable table
tab<- data.frame(Variable=c('Percent selected', 'Selection intensty',
                           'Heritability', 'Selection accuracy',
                           'Additive genetic variance', 'Phenotypic variance', 'Starting population mean',
                           'Number of cycles', 'Expected genetic gain per cycle on average',
                           'Genetic gain per cycle', 'Total response',
                           'Number of years elapsed', 'Genetic gain per year'),
                 Value=round(c(p*100, i, herit, selacc, varA, varP, pop0mean, ncycles, Rpergen_avg,
                         totR/ncycles, totR, ncycles*cycledur, totR/(ncycles*cycledur)), 5))
return(list(plt=plt, plt2=plt2, tab=tab, phenos=cycVecs, plt3=plt3, Caption=Caption))
}


