library(RMySQL)

# Close outstanding results & connections.
conns <- dbListConnections(MySQL())
lapply(conns, function(x) {
  results <- dbListResults(x)
  lapply(results, function(y) dbClearResult(y))
  dbDisconnect(x)
})

connect <- function(dbname=NULL) {
  if(is.null(dbname)) {
    stop("Must give dbname to connect to mysql!")
  }

  if (Sys.getenv('USER') == 'gj1') {
    host = 'ens-research'
    port=3306
    user='ensadmin'
    password='ensembl'
    userpass='ensadmin:ensembl'
  } else {
    host = 'mysql-greg.ebi.ac.uk'
    port=4134
    user='slrsim'
    password='slrsim'
    userpass='slrsim:slrsim'
  }

  con <- dbConnect(MySQL(), host=host, port=port, user=user, password=password, dbname=dbname)
  dbURL = paste("mysql://",userpass,"@",host,":",port,"/",dbname,sep="")
  print(paste("Connected to:",user,"@",host,":",port,"/",dbname))
  print(paste("[",dbURL,"]"))
  return(con)
}
