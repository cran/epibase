
## @knitr eval=TRUE,echo=FALSE,results='hide'
library(knitr)
opts_chunk$set(fig.path='figs/epibase-', fig.keep='last', dev='pdf', fig.width=7, fig.height=7,
               tidy=FALSE, warning=FALSE, fig.show='asis', fig.align='center', cache=FALSE,
               out.width=".6\\textwidth")
options(width=80)


## @knitr results='hide'
library(epibase)


## @knitr 
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
head(ToyOutbreak@records$Fever)
ToyOutbreak@trees


## @knitr 
class(ToyOutbreak@dna)
ToyOutbreak@dna
slotNames(ToyOutbreak@dna)
is.list(ToyOutbreak@dna@dna)
names(ToyOutbreak@dna@dna)
ToyOutbreak@dna@dna$gene1
class(ToyOutbreak@dna@dna$gene1)
class(ToyOutbreak@dna@meta)
head(ToyOutbreak@dna@meta)


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
                  start=onset, end=terminus)
slotNames(oc.dynamic)
oc.dynamic


## @knitr 
as.data.frame(oc.dynamic)


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
head(ToyOutbreakRaw$individuals)


## @knitr 
lapply(ToyOutbreakRaw$records, head)


## @knitr 
head(ToyOutbreakRaw$contacts)
head(ToyOutbreakRaw$contacts.start)
head(ToyOutbreakRaw$contacts.end)


## @knitr 
ToyOutbreakRaw$dna


## @knitr 
ToyOutbreakRaw$trees


## @knitr 
attach(ToyOutbreakRaw)

x <- new ("obkData", individuals=individuals, records=records,
          contacts=contacts, contacts.start=contacts.start,
          contacts.end=contacts.end, dna=dna,
          dna.individualID=dna.info$individualID,
          dna.date=dna.info$date, sample=dna.info$sample, trees=trees)

detach(ToyOutbreakRaw)

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
toKeep <- sample(get.nindividuals(ToyOutbreak),5)
toKeep
x <- subset(ToyOutbreak, individuals=toKeep)
summary(x)


## @knitr 
get.nindividuals(x)
get.nindividuals(x, "records")
get.nindividuals(x, "dna")
get.nindividuals(x, "contacts")


## @knitr 
get.individuals(ToyOutbreak, "contacts")


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
get.data(x,"temperature", showSource=TRUE)


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
get.data(x, "date", where="records", showSource=TRUE)


## @knitr sugarman, warning=TRUE
get.data(x, "sugarman")


## @knitr out.width==".8\\textwidth"
x <- subset(ToyOutbreak, individuals=1:10)
get.ncontacts(x)
plot(x@contacts, main="Contacts in x", label.cex=1.25, vertex.cex=2)


## @knitr 
as.matrix(x@contacts)
as.matrix(x@contacts, "edgelist")


## @knitr 
as.data.frame(x@contacts)


## @knitr eval=FALSE, tidy=FALSE
## subset(x, individuals=NULL, locus=NULL, sequences=NULL,
##        date.from=NULL, date.to=NULL, date.format=NULL, ...)


## @knitr subset1
data(ToyOutbreak)
x1 <- subset(ToyOutbreak, individuals=1:10)
x2 <- subset(ToyOutbreak, get.individuals(ToyOutbreak)[1:10])
identical(x1,x2)


## @knitr 
data(FluH1N1pdm2009)
attach(FluH1N1pdm2009)

x <- new("obkData", individuals = individuals, dna = FluH1N1pdm2009$dna,
      dna.individualID = samples$individualID, dna.date = samples$date,
      trees = FluH1N1pdm2009$trees)

detach(FluH1N1pdm2009)

range(get.data(x, "date"))


## @knitr subsetdate
min.date <- min(get.dates(x))
min.date
min.date+31
x1 <- subset(x, date.to=min.date+31)
summary(x)
summary(x1)


## @knitr lastsubset
temp <- get.data(x, "location", showSource=TRUE)
head(temp)
toKeep <- temp$individualID[temp$location=="Europe"]
x.summerEur <- subset(x, date.from="01/06/2009", date.to="31/08/2009",
                      indiv=toKeep)
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


## @knitr 
plot(x2, "phylo")


## @knitr tree1,fig.keep="last",out.width=".75\\textwidth"
x3 <- make.phylo(x.summerEur, locus=1, ask=FALSE, model="K80")
plot(get.trees(x3)[[1]])
axisPhylo()


## @knitr fig.width=10, out.width="\\textwidth"
set.seed(1)
x <- simuEpi(N=50, D=20, beta=0.01,plot=TRUE,makePhylo=TRUE)
summary(x)
x$dynamics
summary(x$x)


## @knitr 
plot(x$x, "contacts", main="Transmission tree")


## @knitr 
plot(x$x, "phylo")


## @knitr 
data(HorseFlu)
summary(HorseFlu)


## @knitr plottime1
plot(HorseFlu,'timeline')


## @knitr 
args(plotIndividualTimeline)


## @knitr 
plot(HorseFlu,'timeline', what="Vac")


## @knitr 
plotIndividualTimeline(HorseFlu, what="dna", colorBy="yardID", orderBy="yardID",plotNames=TRUE)


## @knitr plotfirst20
plot(HorseFlu,selection=1:20, colorBy="yardID", orderBy="yardID", size=5)


## @knitr 
data(ToyOutbreak)
head(ToyOutbreak@individuals)


## @knitr plotgeo1,dev='png',results='hide',message=FALSE
plot(ToyOutbreak,'geo', location=c('lon','lat'), zoom=14)


## @knitr dev='png',results='hide',message=FALSE
plot(ToyOutbreak,'geo', location=c('lon','lat'), zoom=15,
     colorBy='Sex', center='11')


## @knitr plotmst1, cache=TRUE, fig.keep="last"
data(HorseFlu)
plot(HorseFlu,'mst')


## @knitr 
plot(HorseFlu,'mst',individualID=42)


## @knitr 
data(FluH1N1pdm2009)
attach(FluH1N1pdm2009)

x <- new("obkData", individuals = individuals, dna = FluH1N1pdm2009$dna,
      dna.individualID = samples$individualID, dna.date = samples$date,
      trees = FluH1N1pdm2009$trees)

detach(FluH1N1pdm2009)

summary(x)


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
p <- plotggphy(x, ladderize = TRUE,  branch.unit = "year")


## @knitr 
head(x@individuals)


## @knitr pdh1n1tree3, out.width="\\textwidth"
p <- plotggphy(x, ladderize = TRUE, branch.unit = "year",
               tip.color = "location", tip.size = 3, tip.alpha = 0.75)


