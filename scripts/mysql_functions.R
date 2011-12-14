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
#  print(paste("Connected to:",user,"@",host,":",port,"/",dbname))
#  print(paste("[",dbURL,"]"))
  return(con)
}

disconnect <- function(con) {
  info <- mysqlConnectionInfo(con)

  userpass <- info$user
  host <- info$host
  port <- 'port'
  dbname <- info$dbname  

  dbDisconnect(con)
#  print(paste("Disconnected from:",userpass,"@",host,":",port,"/",dbname))
}

connect.livemirror <- function(dbname=NULL) {
  if(is.null(dbname)) {
    stop("Must give dbname to connect to mysql!")
  }

  if (Sys.getenv('USER') == 'gj1') {
    host = 'ens-livemirror'
    port=3306
    user='ensadmin'
    password='ensembl'
    userpass='ensadmin:ensembl'
  } else {
    host = 'ensembldb.ensembl.org'
    port = 5306
    user = 'anonymous'
    password = NULL
    userpass='anonymous'
  }

  con <- dbConnect(MySQL(), host=host, port=port, user=user, password=password, dbname=dbname)
  dbURL = paste("mysql://",userpass,"@",host,":",port,"/",dbname,sep="")
  #print(paste("Connected to:",user,"@",host,":",port,"/",dbname))
  #print(paste("[",dbURL,"]"))
  return(con)
}


dbUpdateVars <- function(conn, dbtable, dataframe=NULL, primary, vars=colnames(dataframe)) { 
  #print(paste("dbupdateVars for table", dbtable))
  if (!dbExistsTable(conn, dbtable)) { 
    stop("The target table \"", dbtable, "\" does not exist in the database \"", dbGetInfo(conn)$dbname, "\"\n\n", call. = FALSE)
  }
  if (is.null(dataframe)) { 
    stop("The source dataframe is missing, with no default\n\n", call. = FALSE) 
  } 
  if (!(toupper(primary) %in% toupper(names(dataframe)))) {
    #stop("The primary key variable does not exist in the source dataframe\n\n", call. = FALSE) 
  } 
  if (!all(toupper(vars) %in% toupper(names(dataframe)))) { 
    stop("One or more variables do not exist in the source dataframe\n\n", call. = FALSE) 
  } 
  if (!(toupper(primary) %in% toupper(dbListFields(conn, dbtable)))) { 
    #stop("The primary key variable does not exist in the target table\n\n", call. = FALSE) 
  } 

  # Make the variable names OK.
  subst.vars <- gsub('[\\. ]', '_', vars)
  colnames(dataframe) <- subst.vars
  vars <- subst.vars

  if (!all(toupper(vars) %in% toupper(dbListFields(conn, dbtable)))) { 
    print(setdiff(tolower(subst.vars), tolower(dbListFields(conn, dbtable))))
    stop("One or more variables don't exist in the target table\n\n", call. = FALSE) 
  } 

  all.vars <- unique(c(primary,vars))

  for (i in 1:length(all.vars)) {
    cur.var <- all.vars[i]
    if (is.character(dataframe[, cur.var])) {
      dataframe[, cur.var] <- gsub("'", '', dataframe[, cur.var])
    }
  }

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

write.or.update <- function(df, tbl, con, primary) {
  if (dbExistsTable(con, tbl)) {
    #print(head(df))
    dbUpdateVars(con, tbl, df, primary)
  } else {
    print("Table doesn't exist! Writing first row...")
    dbWriteTable(con, tbl, df, row.names=F)
    print("Creating primary key...")
    primary.str <- primary
    if (is.character(df[1, primary]) || is.factor(df[1, primary])) {
      primary.str <- paste(primary, '(64)', sep='')
    }
    db.str <- paste('ALTER TABLE ', tbl, ' ADD UNIQUE (', primary.str, ')', sep='')
    #print(db.str)
    dbSendQuery(con, db.str)

    db.str <- paste('ALTER TABLE ', tbl, ' engine=InnoDB', sep='')
    #print(db.str)
    dbSendQuery(con, db.str)

    print("Done!")
  }
}