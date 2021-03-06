
############################
####  CLASS DEFINITION  ####
############################

## CLASS DESCRIPTION:
## - instance of obkContacts store contacts between individuals
## - employs networkDynamic

setClass("obkContacts",representation(contacts="networkDynamicOrNetworkOrNULL", origin="DateOrNULL"),prototype(contacts=NULL))

setClassUnion("obkContactsOrNULL", c("obkContacts", "NULL"))




######################
####  CONSTRUCTOR ####
######################

## INPUT DESCRIPTION:
## contactFrom: a vector of characters indicating IDs of "sender" individuals
## contactTo: a vector of characters indicating IDs of "receiver" individuals
## directed: should we consider the network as directed
## contactStart: if present, a vector of dates of beginning of contact
## contactEnd: if present, a vector of dates of end of contact
## duration: another way to specify contactEnd, as duration of contact; if dates (not numbers)
## are provided in contactFrom, 'duration' is expected to be days
setMethod("initialize", "obkContacts", function(.Object, from=NULL, to=NULL, directed=FALSE,
                                                start=NULL, end=NULL,duration=NULL) {

    ## RETRIEVE PROTOTYPED OBJECT ##
    x <- .Object

    ## escape if the minimum information is not provided ##
    if(is.null(from) || is.null(to)) return(x)

    ## escape of obkContacts is provided ##
    if(inherits(from, "obkContacts")) return(from)
    if(inherits(to, "obkContacts")) return(to)


    ## PROCESS ARGUMENTS
    if(is.list(from)) from <- unlist(from)
    if(is.list(to)) from <- unlist(to)

    ## fill slots
    from <- as.character(from)
    to <- as.character(to)
    uniqueIDs <- unique(c(from,to))
    numIDs <- length(uniqueIDs)
    numedges <- length(from)
    y <- network.initialize(numIDs,directed=directed,multiple=TRUE)
    network.vertex.names(y) <- uniqueIDs
    ## static network
    if(is.null(start)){
      for(i in 1:numedges){
        v1 <- match(from[i],uniqueIDs)
        v2 <- match(to[i],uniqueIDs)
        add.edge(y,v1,v2)
        }
    }
    ## dynamic network
    if(!is.null(start)){
      if(is.null(end)){
        if(is.null(duration)) stop("Need to specify duration if end is missing")
        # single timestamps
        end <- start+duration
      }

      ## handle 'POSIXct' dates
      if(inherits(start, "POSIXct")){
          start <- as.Date(start)
      }
      if(inherits(end, "POSIXct")){
          end <- as.Date(end)
      }

      ## handle 'Date' dates
      if(inherits(start, "Date")){
          x@origin <- min(start)
          start <- as.numeric(start - x@origin)
          end <- as.numeric(end - x@origin)
      } else {
          x@origin <- NULL
      }
      for(i in 1:numedges){
        v1 <- match(from[i],uniqueIDs)
        v2 <- match(to[i],uniqueIDs)
        add.edge(y,v1,v2)
        activate.edges(y,onset=start[i],terminus=end[i],e=get.edgeIDs(y,v=v1,alter=v2,neighborhood="out"))
      }
    }
    x@contacts <- y

    return(x)
}) # end obkContacts constructor






####################
####  ACCESSORS ####
####################

######################
## get.nindividuals ##
######################
setMethod("get.nindividuals","obkContacts", function(x, ...){
	if(is.null(x@contacts)) return(0)
    return(x@contacts%n%"n")
})


######################
## get.individuals ##
######################
setMethod("get.individuals","obkContacts", function(x, ...){
    if(is.null(x@contacts)) return(NULL)
    return(network.vertex.names(x@contacts))
})


###############
## get.dates ##
###############
setMethod("get.dates","obkContacts", function(x, ...){
    if(is.null(x@contacts)) return(NULL)
    out <- x@origin + get.change.times(x@contacts)
    return(out)
})


################
## get.ndates ##
################
setMethod("get.ndates","obkContacts", function(x, ...){
    if(is.null(x@contacts)) return(0)
    return(length(get.dates(x)))
})



##################
## get.contacts ##
##################
setMethod("get.contacts","obkContacts", function(x, from=NULL, to=NULL, ...){
    if(is.null(x@contacts)) return(0)

    ## handle 'POSIXct' dates
    if(inherits(from, "POSIXct")){
        from <- as.Date(from)
        to <- as.Date(to)
    }

    ## handle character dates
    if(is.character(from)) from <- .process.Date(from)
    if(is.character(to)) to <- .process.Date(to)

    ## handle 'Date' dates
    if(inherits(from, "Date")) from <- as.numeric(from - x@origin)
    if(inherits(to, "Date")) to <- as.numeric(to - x@origin)

    ## extract network
    if(is.null(from)) from <- -1
    if(is.null(to)) to <- Inf
    res <- network.extract(x@contacts, onset=from, terminus=to)

    return(res)
})


######################
#### get.ncontacts ###
######################
setMethod("get.ncontacts","obkContacts", function(x, from=NULL, to=NULL, ...){
    if(is.null(x@contacts)) return(0)
    return(network.edgecount(get.contacts(x, from=from, to=to)))
})


######################
####  SHOW METHOD ####
######################

setMethod ("show", "obkContacts", function(object){
    nindividuals <- get.nindividuals(object)
    ncontacts <- get.ncontacts(object)
    #individualword <- ifelse(nindividuals>1, "individuals", "individual")
    #contactword <- ifelse(nindividuals>1, "individuals", "individual")
    contacts <- get.contacts(object)
    if(class(contacts)[1]=="network"){
      contacttype <- " Contacts = fixed"
    }
    else{
      contacttype <- " Contacts = dynamic"
    }
    cat(paste(" Number of individuals = ", nindividuals, "\n"," Number of contacts = ",ncontacts,"\n",contacttype,"\n",sep=""))
    if(ncontacts>0) print(object@contacts)
    if(contacttype==" Contacts = dynamic"){
        cat("\nDate of origin: ")
        print(object@origin)
    }
})



######################
####  PLOT METHOD ####
######################
## hack to remove the NOTE in R CMD check about:
## plot,obkContacts: no visible binding for global variable ‘y’
if(getRversion() >= "2.15.1")  utils::globalVariables("y")

setMethod ("plot", "obkContacts", function(x, y=NULL, labels=get.individuals(x), ...){
    plot(x@contacts, label=labels, ...)
    return(invisible())
})




##########################
#### AS.MATRIX METHOD ####
##########################
setMethod ("as.matrix", "obkContacts", function(x, matrix.type=c("adjacency","incidence","edgelist"),
                                                use.labels=TRUE, ...){
    g <- x@contacts
    if(is.null(g)) return(NULL)
    matrix.type <- match.arg(matrix.type)
    if(is.networkDynamic(g)) g <- network.collapse(g)
    set.network.attribute(g, "multiple", FALSE)
    out <- as.matrix(g, matrix.type=matrix.type)
    if(use.labels){
        if(matrix.type=="edgelist"){
            lab <- attr(out, "vnames")
            out <- matrix(lab[out], ncol=2)
        }
        if(matrix.type=="incidence"){
            temp <- as.matrix(g, matrix.type="edgelist")
            v.lab <- attr(temp, "vnames")
            e.lab <- apply(matrix(v.lab[temp],ncol=2), 1, paste, collapse="-")
            rownames(out) <- v.lab
            colnames(out) <- e.lab
        }
    }
    return(out)
})



##############################
#### AS.DATA.FRAME METHOD ####
##############################
setMethod ("as.data.frame", "obkContacts", function(x, row.names = NULL, optional = FALSE,
                                                    use.labels=TRUE, ...){
    g <- x@contacts
    if(is.null(g)) return(NULL)
    if(is.networkDynamic(g)) {
        out <- as.data.frame(g, row.names=row.names, optional=optional, ...)
        if(!is.null(x@origin)){
            out$onset <- x@origin + out$onset
            out$terminus <- x@origin + out$terminus
        }
        if(use.labels){
            lab <- get.individuals(x)
            out$tail <- lab[out$tail]
            out$head <- lab[out$head]
        }
    } else {
        out <- as.matrix(x, matrix.type="edgelist", use.labels=use.labels)
    }

    return(out)
})




##################
####  TESTING ####
##################
## NOTE: THIS MUST BE COMMENTED WHEN COMPILING/INSTALLING THE PACKAGE

## ## test constructor / show
## cf <- c("a", "b", "a", "c", "d")
## ct <- c("b", "c", "c", "d", "b")
## onset <- c(1, 2, 3, 4, 5)
## terminus <- c(1.2, 4, 3.5, 4.1, 6)
## oc.static <- new("obkContacts",cf,ct,FALSE) # static network
## oc.dynamic <- new("obkContacts",cf,ct,FALSE,onset,terminus)
## oc.static
## oc.dynamic
