
duplicated2 <- function(x){
  if (sum(dup <- duplicated(x))==0)
    return(dup)
  if (class(x) %in% c("data.frame","matrix"))
    duplicated(rbind(x[dup,],x))[-(1:sum(dup))]
  else duplicated(c(x[dup],x))[-(1:sum(dup))]
}

protcov <- function(x,y){
  
  tablebyrun <- select(x,Protein.Ids,Run,Stripped.Sequence)
  table1 <- select(x,Protein.Ids,Stripped.Sequence)
  
  table1a <- filter(table1,!grepl(";",Protein.Ids)) #tables des ID unique (pas de ;)
  table1b <- filter(table1,grepl(";",Protein.Ids)) #Table des groupes d'ID (avec un ou plusieurs ;)
  table1b$Protein.Idsbis <- table1b$Protein.Ids
  table1bb <-table1b %>%
    mutate(Protein.Ids= strsplit(Protein.Ids, ";")) %>%
    unnest(Protein.Ids)
  #Pour séparer les différents ID (tout en gardant l'info du groupe) afin de pouvoir les utiliser contre la protein database
  
  
  
  table2b <- table1bb[!duplicated(table1bb),]
  
  table2 <- table1a[!duplicated(table1a),]
  
  MyListeAccession <- table2$Protein.Ids
  
  
  MyListeAccessionb <- table2b$Protein.Ids
  
  #table1b <- select(x,Protein.Ids,Run,Stripped.Sequence)
  #table1b$Protein.Ids<-trimws(table1$Protein.Ids,whitespace = ";.*")
  runs <- as.factor(x$Run)
  runs2 <- levels(runs)
  
  
  #accessiononly <- str_split_fixed(names(y),"[|]",n=3)
  #names(y) <- accessiononly[,2]
  #startTime <- Sys.time()
  i <- for (i in 1:length(MyListeAccession))
  {
    #access <- grep(MyListeAccession[i],names(y)) 1.82 min
    #access <- which(str_detect(names(y),MyListeAccession[i])) 1.57 min
    #access <- str_which(names(y),MyListeAccession[i]) 1.29 min
    access <- which(grepl(pattern=MyListeAccession[i],names(y),fixed = TRUE,useBytes = TRUE)) #38sec
    table2$sequencet[i]<-y[[access]]
    
  }
  
  i <- for (i in 1:length(MyListeAccessionb))
  {
    #access <- grep(MyListeAccession[i],names(y)) 1.82 min
    #access <- which(str_detect(names(y),MyListeAccession[i])) 1.57 min
    #access <- str_which(names(y),MyListeAccession[i]) 1.29 min
    access <- which(grepl(pattern=MyListeAccessionb[i],names(y),fixed = TRUE,useBytes = TRUE)) #38sec
    table2b$sequencet[i]<-y[[access]]
    
  }
  
  #sleep_func()
  #endTime <- Sys.time()
  #print(endTime - startTime)
  
  #Calcul du pourcentage de couverture totale     
  total_cova <-calculate_sequence_coverage(table2,protein_sequence=sequencet,peptides = Stripped.Sequence)
  total_covb <-calculate_sequence_coverage(table2b,protein_sequence=sequencet,peptides = Stripped.Sequence) 
  total_cova$coverage <- round(total_cova$coverage,2)
  total_covb$coverage <- round(total_covb$coverage,2)
  
  #Mis en forme de l'info pour chaque Id dans le cas du groupe
  total_covb2 <- select(total_covb,Protein.Idsbis,coverage)
  total_covb2b <- total_covb2[!duplicated(total_covb2),]
  total_covb3 <- total_covb2b %>% 
    group_by(Protein.Idsbis) %>% 
    dplyr::summarise(across(everything(), str_c, collapse="/")) 
  
  colnames(total_covb3) <- c("Protein.Ids","coverage")
  
  total_cova2 <- select(total_cova,Protein.Ids,coverage)
  total_cova3 <- total_cova2[!duplicated(total_cova2),]

  total_coverage <- rbind(total_cova3,total_covb3)
  
  
  table2byrun <- select(table2,-Stripped.Sequence) 
  table2byrun$Protein.Idsbis <- table2byrun$Protein.Ids
  table2byrunb <- select(table2b,-Stripped.Sequence)
  table2byrun <- rbind(table2byrun,table2byrunb)
  
  tablebyrun2 <- merge (table2byrun,tablebyrun, by.x="Protein.Idsbis",by.y="Protein.Ids")
  
  table2byrun <- table2byrun[!duplicated(table2byrun),]
  coveragebyrun = c()
  tablecover = vector("list",length(runs2))
  
  j <- for(j in 1:length(runs2))
  {
    coverrun <- filter(tablebyrun2,Run==runs2[j])
    
    #coverrun2 <- merge(coverrun,table2byrun,by="Protein.Ids")
    coveragebyrun <-  calculate_sequence_coverage(coverrun,protein_sequence =sequencet,peptides = Stripped.Sequence )
    coveragebyrun <- select(coveragebyrun,Protein.Ids,Protein.Idsbis,coverage)
    coveragebyrun$coverage <- round(coveragebyrun$coverage,2)
    colnames(coveragebyrun) <- c("Protein.Ids","Protein.Idsbis",paste0(runs2[j],"_coverage"))
    coveragebyrun <- coveragebyrun[!duplicated(coveragebyrun),]
    tablecover [[j]] <- coveragebyrun
    
  }
  
  
  tablecover2 <- merge(tablecover[[1]],tablecover[[2]],by=c("Protein.Ids","Protein.Idsbis"),all=TRUE)
  
  j <- for(j in 3:length(runs2))
  {
    
    tablecover2 <- merge(tablecover2,tablecover[j],by=c("Protein.Ids","Protein.Idsbis"), all=TRUE)
    
  }
  
  colnames(tablecover2) <- gsub("X","",colnames(tablecover2))
  
  #separer à nouveau les groupes et les uniques
  tablecover2a <- filter(tablecover2,!grepl(";",Protein.Idsbis)) #tables des ID unique (pas de ;)
  tablecover2a$Protein.Idsbis <- NULL
  tablecover2b <- filter(tablecover2,grepl(";",Protein.Idsbis))
  tablecover3b <- tablecover2b %>% 
    group_by(Protein.Idsbis) %>% 
    dplyr::summarise(across(everything(), str_c, collapse="/")) 
  tablecover3b$Protein.Ids <- NULL
  names(tablecover3b)[names(tablecover3b)=='Protein.Idsbis'] <- "Protein.Ids"
    
  tablecover3 <- rbind(tablecover2a,tablecover3b)
  
  total_coverage <- merge(tablecover3,total_coverage,by="Protein.Ids")
  
  
  return(total_coverage)
  
}

formatDIANN <- function(x,y){ #x=data.table and y=total_coverage
  table1 <- select(x,Protein.Ids,Protein.Names,Run,Precursor.Id,Precursor.Charge,contains("PG"),contains("Protein"),Precursor.Normalised)
  #table1$Protein.Ids <- trimws(table1$Protein.Ids,whitespace = ";.*")
  #table1$Protein.Names <- trimws(table1$Protein.Names,whitespace = ";.*")
  
  table1$difference <- table1$PG.Normalised-table1$Precursor.Normalised
  table2 <- filter(table1,difference==0)
  table3 <- select(table2,Protein.Ids,Protein.Names,Run,Precursor.Id) 
  table4 <- table3[duplicated2(table3[1:3]),]
  if (nrow(table4)!=0)
  {
  table5 <- unite(table4,col="merge",Protein.Ids,Protein.Names,Run,sep="_")
  table5$Precursor.Id <- paste0(table5$Precursor.Id,"/")
  table6 <- reshape2::dcast(data=table5,formula=merge~Precursor.Id,value.var = 'Precursor.Id')
  table6 [is.na(table6)] <-""
  namecolonne <- colnames(table6)
  namecolonne <- namecolonne[-1]
  table7 <- unite(table6,col="Sequence(s)",c(namecolonne),sep="")
  table7b <- str_split_fixed(table7$merge,"[_]",n=4)
  table7b <- as.data.frame(table7b)
  table7b<-unite(table7b,col="Protein.Names",V2,V3,sep="_")
  table7 <- cbind(table7b,table7)
  table7$merge <- NULL
  colnames(table7)<- c("Protein.Ids","Protein.Names","Run","Precursor.Id")
  
  table3b <- anti_join(table3,table4)
  
  table3 <- rbind(table3b,table7)
  }
  
  #Je compte le nb de peptide ou de run
  aggregat1<- reshape2::dcast(data = table1,formula = Protein.Ids+Protein.Names ~ Run,value.var = 'Precursor.Id',fun.aggregate = function(x) length(unique(x)))
  
  #Changer titre colonne run en nb de peptide
  aggregat1r <-rename_with(aggregat1,~paste0("nb.peptide.",.), -c(Protein.Ids,Protein.Names))
  
  
  options(digits = 3, scipen = -2)
  
  #Aggregation en fonction de PG.MaxLFQ, j'ai fait mean, puisque la valeur est partout la même
  
  #aggregat2 <- reshape2::dcast(data = table1,formula = Protein.Ids+Protein.Names ~ Run,value.var = 'PG.MaxLFQ',fun.aggregate = mean)
  aggregat2 <- reshape2::dcast(data = table1,formula = Protein.Ids+Protein.Names ~ Run,value.var = 'PG.MaxLFQ',fun.aggregate = mean)
  aggregat2r <-rename_with(aggregat2,~paste0("MaxLFQ.",.), -c(Protein.Ids,Protein.Names))
  
  aggregat3 <- reshape2::dcast(data = table1,formula = Protein.Ids+Protein.Names ~ Run,value.var = 'Protein.Q.Value',fun.aggregate = mean)
  aggregat3r <-rename_with(aggregat3,~paste0("Protein.Q.Value.",.), -Protein.Ids)
  
  #check if you have an unique sequence by protein
  aggregat4r <- reshape2::dcast(data=table3,formula = Protein.Ids+Protein.Names~Run ,value.var = 'Precursor.Id',fun.aggregate = function(x) length(unique(x)))
  
  aggregat4 <- reshape2::dcast(data=table3,formula = Protein.Ids+Protein.Names~Run ,value.var = 'Precursor.Id')
  
  
  
  
  final <- Merge (aggregat1r,aggregat2r,aggregat3r,aggregat4,y,id=~Protein.Ids,all=TRUE)
  
  
  
  final
}








count_condition <- function(x){
  
  cond <- as.data.frame(table(x$Run))
  cond$Var1 <- gsub("_[0-9]","",cond$Var1)
  recap <- ddply(cond,.(Var1),c("nrow"))
  colnames(recap) <- c("Conditions","Nb of replicat")
  recap
}
ggpairsformat <- function(x){
  df <- select(x,starts_with("MaxLFQ"))
  df[is.na(df)] <- 0
  dflog <- log2(df) 
  dflog
}

boxplotformat <- function(x){
  
  #coldata <- filter(x,Protein.Ids == y)
  coldata <- select(x,starts_with("MaxLFQ"))
  coldata[is.na(coldata)] <- 0
  coldata <- log2(coldata)
  a <- ncol(coldata)
  titrecol <- colnames(coldata)
  groupe <- gsub("MaxLFQ.","",titrecol)
  groupe <- gsub("[_][0-9]","",groupe)
  colnames(coldata)<-groupe
  IDProteins <- x$Protein.Ids
  boxplotdata <- cbind (IDProteins,coldata)
  boxplotdata2 <- pivot_longer(boxplotdata,cols = -"IDProteins")
  boxplotdata2
  
  
  
}

Barformat <- function(x){
  
 
  coldata2 <- select(x,starts_with("nb.peptide."))
  a <- ncol(coldata2)
  titrecol <- colnames(coldata2)
  groupe <- gsub("nb.peptide.","",titrecol)
  groupe <- gsub("[_][0-9]","",groupe)
  colnames(coldata2)<-groupe
  nbProt <- colSums(coldata2!=0)
  formatbarplot <- cbind(groupe,nbProt)
  #formatbarplot2 <- t(formatbarplot)
  formatbarplot2 <- as.data.frame(formatbarplot)
  formatbarplot2$nbProt <- as.numeric(formatbarplot2$nbProt)
  formatbarplot2
}


data_summary <- function(data, varname, groupnames){
  
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, summarize,mean=round(mean(nbProt),2),sd=round(sd(nbProt),2))
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}



### Upset format and extract

formatUpset <- function(x) {
  
 
  coldata <- select(x,starts_with("MaxLFQ"))
  titrecol <- colnames(coldata)
  groupe <- gsub("MaxLFQ.","",titrecol)
  groupe <- gsub(".raw","",groupe)
  groupe <- gsub("[_][0-9]","",groupe)
  #groupe <- gsub("[0-9]","",groupe)
  groupe <- as.factor(groupe)
  nomgroupe <- levels(groupe)
  
  
  
  condition_list = vector(mode = "list",nlevels(groupe))
  
  i <- for (i in 1:nlevels(groupe))
  {
    name = nomgroupe[i]
    data <-  select(x,starts_with(paste("MaxLFQ.",nomgroupe[i],sep="")))
    replicat<-ncol(data)
    data[is.na(data)] <- 0
    data [data>1] <- 1
    colnames(data) <- gsub ("MaxLFQ.","",colnames(data))
    zeropartout <- apply(data,1,sum)
    titre<-colnames(data)
    data <- cbind(data,select(x,Protein.Ids),zeropartout)
    data <- subset(data,data$zeropartout>0)
    data$zeropartout <- 1
    colnames(data)[colnames(data) == 'zeropartout'] <- name
    data <- data[,-(1:replicat)]
    
    groupe1 <- paste("groupe_",name,sep="")
    assign (groupe1,data)
    condition_list[[i]]=data
    
    
    
    
  }
  
  upsetdata <- merge(condition_list[[1]],condition_list[[2]],by="Protein.Ids",all=TRUE)
 
  
  i <- for (i in 3:nlevels(groupe))
  {
    upsetdata <- merge(upsetdata,condition_list[[i]],by="Protein.Ids",all=TRUE)
   
    
  }
  upsetdata[is.na(upsetdata)]<- 0
  
  
 
  
  
  
  
  return(upsetdata)
  
}

#### To create the commun list data

formatListUpsetData <- function(x) {
  
  group <- colnames(x)
  
  group <- group[-1]
  
  group <- as.factor(group)
  
  
  nameintersection <- paste(levels(group),collapse = "_")
  
  
  
  
  data_with_intersection <- x %>%
    unite(col = "intersection", -c("Protein.Ids"), sep = "")
  
  
  
  allintersection <- as.factor(data_with_intersection$intersection)
  
  levelintersection <- levels(allintersection)
  
  nbintersection <- nlevels(allintersection)
  
  versus <- str_split(nameintersection,"_")
  
  versus <- unlist(versus)
  
  com=vector(mode = "list",nbintersection)
  
  i <- for (i in 1:nbintersection)
  {
  
    com[[i]] <- filter(data_with_intersection,intersection==levelintersection[i])
    toto <- as.numeric(strsplit(as.character(levelintersection[i]),"")[[1]])
    listintersec <- as.data.frame(cbind(versus,toto))
    listintersec <- filter(listintersec,toto!=0)
    listintersec2 <- paste(listintersec$versus,collapse = "_&_")
    colnames(com[[i]]) <- c(listintersec2,"intersection")
    com[[i]] <- com[[i]][-2]
    
    
    
  }
 
  dataintersection <- rbind.fill(com) 
  return(dataintersection)
  
  
}


#### format data pirateplot

formatpirateplot <- function(x) {
  
  coldata <- select(x,starts_with("MaxLFQ"))
  titrecol <- colnames(coldata)
  groupe <- gsub("MaxLFQ.","",titrecol)
  groupe <- gsub(".raw","",groupe)
  groupe <- gsub("[_][0-9]","",groupe)
  #groupe <- gsub("[0-9]","",groupe)
  groupe <- as.factor(groupe)
  nomgroupe <- levels(groupe)
  
  
  final=c()
  
  i <- for (i in 1:nlevels(groupe))
  {
    name = nomgroupe[i]
    data <-  select(x,starts_with(paste("MaxLFQ.",nomgroupe[i],sep="")))
    sd_prot <- apply(data,1,sd)
    moy_prot <- apply(data,1,mean)
    cv <- sd_prot/moy_prot
    colnames(data) <- gsub ("MaxLFQ.","",colnames(data))
    data <- cbind(select(x,Protein.Ids),moy_prot,cv)
    data$condition <- name
    
    final <- rbind(final,data)
    final <- na.omit(final)
    
    
    
  }
  return(final)
}