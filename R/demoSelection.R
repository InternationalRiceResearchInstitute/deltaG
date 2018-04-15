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
                        numparents= 20, ncycles= 2, 
                        cycledur=7){

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
Rpergen<- i*sqrt(varA)*selacc

#cycle zero vec
cycle0 <- data.frame(Phenotypic_Value = 
                       rnorm(popsize, pop0mean, sqrt(varP)))

#Make cycle 1 based on cycle 0
if(ncycles==1){
  #get response based on R= h2S
  ordvec<- cycle0[order(-cycle0[,1]),]
  mnSel<- mean(ordvec[1:c(length(ordvec)*p)])
  mnTot<- mean(ordvec)
  S<- mnSel-mnTot
  totR<- S*herit
}else{
  #expected response per generation
  #based on R= i h stdA
  totR<- Rpergen*ncycles
}

#population 1 mean
pop1mean<- pop0mean+totR

#cycle 1
cycle1 <- data.frame(Phenotypic_Value = 
                       rnorm(popsize, pop1mean, sqrt(varP)))

#Now, combine your two dataframes into one.  First make a new column in each that will be a variable to identify where they came from later.
cycle0$Population <- 'Before selection'
cycle1$Population <- 'After selection'

#and combine into your new data frame vegLengths
cycVecs <- rbind(cycle0, cycle1)
cycVecs <- cycVecs[order(-cycVecs[,1]),]
cycVecs<- data.frame(cycVecs, id=paste('id', row.names(cycVecs)))

#get xlim
popMin<- min(cycVecs[,1])
popMax<- max(cycVecs[,1])
rg<- popMax-popMin
xmin<-popMin- rg *0.5
xmx<- popMax+ rg *0.5

#total gain
totGain<- pop1mean-pop0mean

#plot title
plotTit<- paste("Expected gain from selection = ", round(totGain,3), 
"after", ncycles*cycledur, 'years')

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

plt2<- ggplot2::ggplot(cycVecs, aes(x = reorder(id, Phenotypic_Value), 
                    y = Phenotypic_Value, fill = Population)) + 
  geom_bar(stat='identity', position='identity',
           color='black', size= 4/popsize)+  
  scale_fill_manual(values = c("peachpuff","grey65"))+
  xlab("Individual") +
  ylab("Phenotypic value") +
  ggtitle("Phenotypic values by population")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

#variable table
tab<- data.frame(Varible=c('Percent selected', 'Selection intensty',
                           'Heritability', 'Selection accuracy', 
                           'Additive genetic variance', 'Phenotypic variance', 'Starting population mean',
                           'Number of cycles', 'Expected genetic gain per cycle', 'Total response',
                           'Number of years elapsed', 'Expected genetic gain per year'), 
                 Value=c(p*100, i, herit, selacc, varA, varP, pop0mean, ncycles, 
                         Rpergen, totR, ncycles*cycledur, totR/(ncycles*cycledur)))
tab<- format(tab, digits=3)
return(list(plt=plt, plt2=plt2, tab=tab, phenos=cycVecs))
}






