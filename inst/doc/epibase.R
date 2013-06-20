
## @knitr eval=TRUE,echo=FALSE,results='hide'
library(knitr)
opts_chunk$set(fig.path='figs/epibase-', fig.keep='last', dev='pdf', fig.width=7, fig.height=7,
               tidy=FALSE, warning=FALSE, fig.show='asis', fig.align='center', cache=FALSE,
               out.width=".6\\textwidth")
options(width=80)


## @knitr 
library(epibase)
getClassDef("obkData")


## @knitr 
new("obkData")


## @knitr 
data(ToyOutbreak)
class(ToyOutbreak)
slotNames(ToyOutbreak)
head(ToyOutbreak)
summary(ToyOutbreak)


## @knitr 
head(ToyOutbreak@individuals)
head(ToyOutbreak@samples)
ToyOutbreak@trees


## @knitr 
class(ToyOutbreak@dna)
ToyOutbreak@dna
slotNames(ToyOutbreak@dna)
is.list(ToyOutbreak@dna@dna)
names(ToyOutbreak@dna@dna)
ToyOutbreak@dna@dna$gene1
class(ToyOutbreak@dna@dna$gene1)


## @knitr 
cf <- c("a", "b", "a", "c", "d")
ct <- c("b", "c", "c", "d", "b")
oc.static <- new("obkContacts", cf, ct, directed=FALSE)
slotNames(oc.static)
oc.static


## @knitr graphstat,out.width=".7\\textwidth"
plot(oc.static, main="Static contact network")


## @knitr 
onset <- c(1, 2, 3, 4, 5)
terminus <- c(1.2, 4, 3.5, 4.1, 6)
oc.dynamic <- new("obkContacts",cf,ct, directed=FALSE,
                  contactStart=onset, contactEnd=terminus)
slotNames(oc.dynamic)
oc.dynamic


## @knitr 
data.frame(onset,terminus,ct,cf)


## @knitr dynNet,out.width=".9\\textwidth"
par(mfrow=c(2,2))
plot(oc.dynamic@contacts,main="oc.dynamic - collapsed graph",
     displaylabels=TRUE)
plot(get.contacts(oc.dynamic, from=0, to=2),
     main="oc.dynamic - time 0--2", displaylabels=TRUE)
plot(get.contacts(oc.dynamic, from=2, to=4),
     main="oc.dynamic - time 2--4", displaylabels=TRUE)
plot(get.contacts(oc.dynamic, from=4, to=6),
     main="oc.dynamic - time 4--6", displaylabels=TRUE)


## @knitr 
new("obkData")


## @knitr 
data(ToyOutbreakRaw)
class(ToyOutbreakRaw)
names(ToyOutbreakRaw)


## @knitr 
head(ToyOutbreakRaw$samples)
head(ToyOutbreakRaw$individuals)
x <- new("obkData", samples=ToyOutbreakRaw$samples,
         individuals=ToyOutbreakRaw$individuals)
head(x)


## @knitr 
head(ToyOutbreakRaw$contacts)
head(ToyOutbreakRaw$contacts.start)
head(ToyOutbreakRaw$contacts.end)


## @knitr 
head(ToyOutbreakRaw$clinical$Fever)


## @knitr 
ToyOutbreakRaw$dna


## @knitr 
ToyOutbreakRaw$trees


## @knitr 
x <- new ("obkData", individuals=ToyOutbreakRaw$individuals,
          samples=ToyOutbreakRaw$samples,
          clinical=ToyOutbreakRaw$clinical, contacts=ToyOutbreakRaw$contacts,
          contacts.start=ToyOutbreakRaw$contacts.start,
          contacts.end=ToyOutbreakRaw$contacts.end,
          dna=ToyOutbreakRaw$dna, trees=ToyOutbreakRaw$trees)
head(x)
summary(x)


## @knitr eval=TRUE,echo=FALSE,results='hide'
opts_chunk$set(eval=FALSE, echo=FALSE)


## @knitr eval=TRUE,echo=FALSE,results='hide'
opts_chunk$set(eval=TRUE, echo=TRUE)


## @knitr eval=FALSE
## myFunction(x, y="foo")


## @knitr eval=FALSE
## myFunction(x, y="bar")


## @knitr 
data(ToyOutbreak)
set.seed(1)
get.nsamples(ToyOutbreak)
toKeep <- sample(1:nrow(ToyOutbreak@samples), 5)
toKeep
x <- subset(ToyOutbreak, row.samples=toKeep)
summary(x)


## @knitr 
get.nindividuals(x)
get.nindividuals(x, "contacts")
get.individuals(x)



## @knitr 
get.individuals(ToyOutbreak, "contacts")


## @knitr 
get.nsamples(x)
get.samples(x)


## @knitr 
get.nlocus(x)
get.locus(x)


## @knitr 
get.nsequences(x)
get.nsequences(x, "bylocus")
get.sequences(x)


## @knitr 
get.trees(x)


## @knitr 
get.dna(x)


## @knitr 
get.dna(x, locus=2)


## @knitr 
get.dna(x, id=c("311","222"))


## @knitr 
get.sequences(x)
identical(get.dna(x, id=c("311","222")), get.dna(x, id=c(2,1)))


## @knitr 
get.ncontacts(ToyOutbreak)
get.individuals(ToyOutbreak@contacts)
get.individuals(x)
get.ncontacts(x)


## @knitr getData1
get.data(x,"temperature")


## @knitr 
get.data(x, "Sex")


## @knitr 
get.data(x, c("Sex","Age","infector"))


## @knitr 
get.data(x, c("Sex","Age","infector"), showSource=TRUE)


## @knitr 
get.data(x, "date")


## @knitr 
get.data(x, "date", showSource=TRUE)


## @knitr 
get.data(x, "date", where="clinical", showSource=TRUE)


## @knitr sugarman, warning=TRUE
get.data(x, "sugarman")


## @knitr warning=TRUE
x@clinical <- NULL
get.data(x, "date", where="clinical")


## @knitr out.width==".8\\textwidth"
x <- subset(ToyOutbreak, individuals=1:10)
get.ncontacts(x)
plot(x@contacts, main="Contacts in x", label.cex=1.25, vertex.cex=2)


## @knitr 
as.matrix(x@contacts)
as.matrix(x@contacts, "edgelist")


## @knitr eval=FALSE, tidy=FALSE
## subset(x, individuals=NULL, samples=NULL, locus=NULL, sequences=NULL,
##        date.from=NULL, date.to=NULL, date.format=NULL,
##        row.individuals=NULL, row.samples=NULL,...)


## @knitr subset1
data(ToyOutbreak)
set.seed(1)
get.nsamples(ToyOutbreak)
toKeep <- sample(1:nrow(ToyOutbreak@samples), 10)
toKeep
x <- subset(ToyOutbreak, row.samples=toKeep)
summary(x)
get.individuals(x)


## @knitr 
x1 <- subset(x, indiv=c("60","168"))


## @knitr 
x2 <- subset(x, indiv=c(3,5))
identical(x1,x2)


## @knitr 
data(FluH1N1pdm2009)
x <- new("obkData", individuals = FluH1N1pdm2009$individuals, samples =
         FluH1N1pdm2009$samples, dna = FluH1N1pdm2009$dna, trees =
         FluH1N1pdm2009$trees)
range(get.data(x, "date", where="samples"))


## @knitr subsetdate
min.date <- min(get.data(x, "date", where="samples"))
min.date
min.date+31
x1 <- subset(x, date.to=min.date+31)
summary(x)
summary(x1)


## @knitr 
toKeep <- get.data(x, "location")=="Mexico"
sum(toKeep)
x.mex <- subset(x, row.individuals=toKeep)
summary(x.mex)
head(x.mex)


## @knitr lastsubset
x.summerEur <- subset(x, date.from="01/06/2009", date.to="31/08/2009",
                      row.indiv=get.data(x, "location")=="Europe")
summary(x.summerEur)
head(x.summerEur)


## @knitr 
x.summerEur@trees <- NULL
get.nsequences(x.summerEur)


## @knitr makephylo, fig.keep="all"
x2 <- make.phylo(x.summerEur)
summary(x2)


## @knitr out.width=".75\\textwidth"
library(ape)
plot(get.trees(x2)[[1]])
axisPhylo()


## @knitr tree1,fig.keep="last",out.width=".75\\textwidth"
tree1 <- make.phylo(x.summerEur, locus=1, ask=FALSE, model="K80",
                    plot=TRUE, color.by="dat")
axisPhylo()


## @knitr fig.width=10, out.width="\\textwidth"
set.seed(1)
x <- simuEpi(N=50,beta=0.01,showPlots=TRUE,makePhylo=TRUE)
summary(x)


## @knitr 
data(HorseFlu)
summary(HorseFlu)


## @knitr plottime1
plot(HorseFlu,'timeline')


## @knitr plotfirst20
plot(HorseFlu,selection=1:20)


## @knitr 
plot(HorseFlu,selection=1:20,colorBy='yardID')


## @knitr 
plot(HorseFlu,selection=1:20,colorBy='yardID',orderBy='yardID')


## @knitr 
data(ToyOutbreak)
head(ToyOutbreak@individuals)


## @knitr plotgeo1,dev='png',results='hide',message=FALSE
plot(ToyOutbreak,'geo', location=c('lon','lat'), isLonLat=TRUE, zoom=14)


## @knitr dev='png',results='hide',message=FALSE
plot(ToyOutbreak,'geo', location=c('lon','lat'), isLonLat=TRUE, zoom=15,
     colorBy='Sex', center='11')


## @knitr plotmst1, cache=TRUE, fig.keep="last"
data(HorseFlu)
plot(HorseFlu,'mst')


## @knitr 
plot(HorseFlu,'mst',individualID=42)


## @knitr 
data(FluH1N1pdm2009)
x <- new("obkData", individuals = FluH1N1pdm2009$individuals,
         samples = FluH1N1pdm2009$samples, dna = FluH1N1pdm2009$dna,
         trees = FluH1N1pdm2009$trees)
head(x)


## @knitr 
get.trees(x)
tre <- get.trees(x)[[1]]
tre


## @knitr pdh1n1tree1, out.width="0.8\\textwidth"
plot(get.trees(x)[[1]], show.tip=FALSE)


## @knitr pdh1n1,dev='png'
plot(x, colorBy="location", orderBy="location")


## @knitr fig.keep="last", out.width="0.8\\textwidth"
plotggphy(x)


## @knitr pdh1n1tree2, out.width="0.8\\textwidth"
p <- plotggphy(x, ladderize = TRUE, build.tip.attribute = TRUE,
               branch.unit = "year", tip.dates = "date")


## @knitr 
head(x@individuals)


## @knitr pdh1n1tree3, out.width="\\textwidth"
p <- plotggphy(x, ladderize = TRUE, build.tip.attribute = TRUE,
               branch.unit = "year", tip.dates = "date", tip.colour = "location",
               tip.size = 3, tip.alpha = 0.75)


