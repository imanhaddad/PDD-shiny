formatDIANN <- function(x){
  table1 <- select(x,Protein.Ids,Protein.Names,Run,Modified.Sequence,Precursor.Charge,contains("PG"),contains("Protein"),Precursor.Normalised)
  table1$difference <- table1$PG.Normalised-table1$Precursor.Normalised
  table2 <- filter(table1,difference==0)
  table3 <- select(table2,Protein.Ids,Run,Modified.Sequence)
  
  
  #Je compte le nb de peptide ou de run
  aggregat1<- reshape2::dcast(data = table1,formula = Protein.Ids+Protein.Names ~ Run,value.var = 'Modified.Sequence',fun.aggregate = function(x) length(unique(x)))
  
  #Changer titre colonne run en nb de peptide
  aggregat1r <-rename_with(aggregat1,~paste0("nb.peptide.",.), -c(Protein.Ids,Protein.Names))
  
  
  options(digits = 3, scipen = -2)
  
  #Aggregation en fonction de PG.MaxLFQ, j'ai fait mean, puisque la valeur est partout la même
  
  aggregat2 <- reshape2::dcast(data = table1,formula = Protein.Ids ~ Run,value.var = 'PG.MaxLFQ',fun.aggregate = mean)
  aggregat2r <-rename_with(aggregat2,~paste0("MaxLFQ.",.), -Protein.Ids)
  
  aggregat3 <- reshape2::dcast(data = table1,formula = Protein.Ids ~ Run,value.var = 'Protein.Q.Value',fun.aggregate = mean)
  aggregat3r <-rename_with(aggregat3,~paste0("Protein.Q.Value.",.), -Protein.Ids)
  
  aggregat4 <- reshape2::dcast(data=table3,formula = Protein.Ids~Run ,value.var = 'Modified.Sequence')
  
  
  final <- Merge (aggregat1r,aggregat2r,aggregat3r,aggregat4,id=~Protein.Ids,all=TRUE)
  final
}



ggpairsformat <- function(x){
  df <- select(x,starts_with("MaxLFQ"))
  #df[is.na(df)] <- 0
  dflog <- log2(df) 
  dflog
}

boxplotformat <- function(x){
  
  #coldata <- filter(x,Protein.Ids == y)
  coldata <- select(x,starts_with("MaxLFQ"))
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
  
  
  #data_with_intersection <- upsetdata %>%
   # unite(col = "intersection", -c("Protein.Ids"), sep = "")
  
 # test <- data_with_intersection %>%
  #  group_by(intersection)%>%
   # nest()%>%
    #select(data)%>%
    #unlist(recursive = F)
  
  
  
  
  return(upsetdata)
  
}

#### To create the commun list




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
    
    
    
    
  }
  return(final)
}