datacache <- new.env(hash=TRUE, parent=emptyenv())

org.Cquinquefasciatus.eg <- function() showQCData("org.Cquinquefasciatus.eg", datacache)
org.Cquinquefasciatus.eg_dbconn <- function() dbconn(datacache)
org.Cquinquefasciatus.eg_dbfile <- function() dbfile(datacache)
org.Cquinquefasciatus.eg_dbschema <- function(file="", show.indices=FALSE) dbschema(datacache, file=file, show.indices=show.indices)
org.Cquinquefasciatus.eg_dbInfo <- function() dbInfo(datacache)

org.Cquinquefasciatus.egORGANISM <- "Culex quinquefasciatus"

.onLoad <- function(libname, pkgname)
{
    ## Connect to the SQLite DB
    dbfile <- system.file("extdata", "org.Cquinquefasciatus.eg.sqlite", package=pkgname, lib.loc=libname)
    assign("dbfile", dbfile, envir=datacache)
    dbconn <- dbFileConnect(dbfile)
    assign("dbconn", dbconn, envir=datacache)

    ## Create the OrgDb object
    sPkgname <- sub(".db$","",pkgname)
    db <- loadDb(system.file("extdata", paste(sPkgname,
      ".sqlite",sep=""), package=pkgname, lib.loc=libname),
                   packageName=pkgname)    
    dbNewname <- AnnotationDbi:::dbObjectName(pkgname,"OrgDb")
    ns <- asNamespace(pkgname)
    assign(dbNewname, db, envir=ns)
    namespaceExport(ns, dbNewname)
        
    packageStartupMessage(AnnotationDbi:::annoStartupMessages("org.Cquinquefasciatus.eg.db"))
}

.onUnload <- function(libpath)
{
    dbFileDisconnect(org.Cquinquefasciatus.eg_dbconn())
}

