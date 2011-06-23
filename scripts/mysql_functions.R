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

dbUpdateVars <- function(conn, dbtable, dataframe=NULL, primary, vars=names(dataframe)) { 
  if (!dbExistsTable(conn, dbtable)) { 
    stop("The target table \"", dbtable, "\" doesn't exist in the database \"", dbGetInfo(conn)$dbname, "\"\n\n", call. = FALSE)
  }
  if (is.null(dataframe)) { 
    stop("The source dataframe is missing, with no default\n\n", call. = FALSE) 
  } 
  if (!(toupper(primary) %in% toupper(names(dataframe)))) {
    stop("The primary key variable doesn't exist in the source dataframe\n\n", call. = FALSE) 
  } 
  if (!all(toupper(vars) %in% toupper(names(dataframe)))) { 
    stop("One or more variables don't exist in the source dataframe\n\n", call. = FALSE) 
  } 
  if (!(toupper(primary) %in% toupper(dbListFields(con, dbtable)))) { 
    stop("The primary key variable doesn't exist in the target table\n\n", call. = FALSE) 
  } 
  if (!all(toupper(vars) %in% toupper(dbListFields(con, dbtable)))) { 
    stop("One or more variables don't exist in the target table\n\n", call. = FALSE) 
  } 

  all.vars <- unique(c(primary,vars))

  if(length(all.vars) > 1) { 
    pastedvars <- paste("'", apply(dataframe[, all.vars], 1, paste, collapse="', '"), "'", sep="") 
  } else { 
    pastedvars <- paste("'", dataframe[, all.vars], "'", sep="") 
  }

  varlist <- paste(dbtable, "(", paste(all.vars, collapse=", "), ")", sep="") 
  datastring <- paste("(", paste(pastedvars, collapse="), ("), ")", sep="")
  toupdate <- paste(paste(vars, "=VALUES(", vars, ")", sep=""), collapse=", ")
  sqlstring <- paste("INSERT INTO", varlist, "VALUES", datastring, "ON DUPLICATE KEY UPDATE", toupdate)
  dbSendQuery(conn, sqlstring) 
} 