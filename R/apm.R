apm<-function(mix, mono, ry=NULL, method="loreau"){
  
  #Error interceptions
  if (is.data.frame(mix)==TRUE|is.matrix(mix)==TRUE) {} else {stop("mix must be a matrix or a data frame")}
  if (is.data.frame(mono)==TRUE|is.matrix(mono)==TRUE) {} else {stop("mono must be a matrix or a data frame")}
  if (is.numeric(mix)==TRUE|is.numeric(mono)==TRUE) {} else {"mix and mono must contain numeric values"}
  if (method=="loreau"|method=="fox") {} else {stop("method must be eihter loreau or fox")}

  if (is.null(colnames(mix))==TRUE|is.null(colnames(mono))==TRUE){stop("Provide species names/codes as column names in mix and mono")}
  for (i in 1:ncol(mix)){if (colnames(mix)[i] %in% colnames(mono)) {} else {stop(paste(colnames(mix)[i], " was not found in mono", sep=""))}}
  
  if (is.data.frame(mix)==TRUE){mix<-as.matrix(mix)}
  if (is.data.frame(mono)==TRUE){mono<-as.matrix(mono)}
  
  #Calculate expected relative yield for each species
  if (is.null(ry)==TRUE){
    ry<-matrix(1/ncol(mix), ncol=ncol(mix), nrow=nrow(mix))
    colnames(ry)<-colnames(mix)
    rownames(ry)<-rownames(mix)}
  else {
    if (is.data.frame(ry)==TRUE|is.matrix(ry)==TRUE) {} else {stop("ry must be a matrix or a data frame")}
    if (is.numeric(ry)==FALSE){stop("ry must contain numeric values")}
    if (nrow(ry)!=nrow(mix)|ncol(ry)!=ncol(mix)){stop("mix and ry must have the same dimension")}
    if (is.null(colnames(ry))==TRUE){stop("Provide species names/codes as column names in ry")}
    for (i in 1:ncol(ry)){if (colnames(ry)[i] %in% colnames(mix)) {} else {stop(paste(colnames(ry)[i], " was not found in mix", sep=""))}}
    sumry<-apply(ry, 1, sum)
    if (length(which(sumry>1))>0){stop(paste("Sum of expected relative yields greater than 1"))}
    if (is.data.frame(ry)==TRUE){ry<-as.matrix(ry)}}
  
  #Detect zero or NA values in monocultures
  if (TRUE %in% is.na(apply(mono, 2, mean, na.rm=TRUE))){
    index<-which(is.na(apply(mono, 2, mean, na.rm=TRUE))==TRUE)
    for (i in 1:length(index)){if (colnames(mono)[index[i]] %in% colnames(mix)) {stop(paste("No monoculture biomass for ", colnames(mono)[index[i]], sep=""))}}}
  
  if (0 %in% apply(mono, 2, mean, na.rm=TRUE)){
    index<-which(apply(mono, 2, mean, na.rm=TRUE)==0)
    for (i in 1:length(index)){if (colnames(mono)[index[i]] %in% colnames(mix)) {stop(paste("Monoculture biomass for ", colnames(mono)[index[i]], " is equal to zero", sep=""))}}}
  
  mono<-t(as.matrix(apply(mono, 2, mean, na.rm=TRUE)))
  
  #Additive partitioning
  N<-apply(ry, 1, function(x){sum(x>0)}) #Sown species richness in each plot
  obs.tot.y<-apply(mix, 1, sum) #Total observed yield of the mixtures
  
  obs.ry<-mix
  exp.y<-ry
  for (i in 1:ncol(mix)){
    obs.ry[,i]<-obs.ry[,i]/mono[,colnames(obs.ry)[i]] #Observed relative yield of each species in mixtures
    exp.y[,i]<-exp.y[,i]*mono[,colnames(exp.y)[i]]} #Expected yield of each species in mixtures
  
  exp.tot.y<-apply(exp.y, 1, sum) #Total expected yield of the mixtures
  delta.ry<-obs.ry-ry[,match(colnames(ry), colnames(obs.ry))] #Deviation from expected relative yield of each species in the mixtures
  
  mean.delta.ry<-c()
  mean.mono<-c()
  covar<-c()
  
  for (i in 1:nrow(ry)){
    mean.delta.ry[i]<-mean(delta.ry[i, colnames(ry)[which(ry[i,]>0)]])
    mean.mono[i]<-mean(mono[,colnames(ry)[which(ry[i,]>0)]])
    covar[i]<-covariance(delta.ry[i, colnames(ry)[which(ry[i,]>0)]], mono[, colnames(ry)[which(ry[i,]>0)]])}
  
  if (method=="loreau"){
    results<-matrix(ncol=3, nrow=nrow(mix))
    results[,1]<-obs.tot.y-exp.tot.y #Net biodiversity effect
    results[,2]<-N*mean.delta.ry*mean.mono #Complementarity effect
    results[,3]<-N*covar #Selection effect
    colnames(results)<-c("NBE", "CE", "SE")
    rownames(results)<-rownames(mix)}
  
  if (method=="fox"){
    obs.ry.tot<-apply(obs.ry, 1, sum) #Observed relative yield total
    obs.freq<-obs.ry
    for (i in 1:ncol(obs.freq)){obs.freq[,i]<-obs.freq[,i]/obs.ry.tot} #Observed frequency of each species in mixtures
    
    covar1<-c()
    covar2<-c()
    
    for (i in 1:nrow(ry)){
      
      col<-colnames(ry)[which(ry[i,]>0)]
      covar1[i]<-covariance(obs.ry[i, col]-obs.freq[i, col], mono[, col])
      covar2[i]<-covariance(obs.freq[i, col]-ry[i, col], mono[, col])}
    
    results<-matrix(ncol=4, nrow=nrow(mix))
    results[,1]<-obs.tot.y-exp.tot.y #Net biodiversity effect
    results[,2]<-N*mean.delta.ry*mean.mono #Trait-independent complementarity effect
    results[,3]<-N*covar1 #Trait-dependent complementarity effect
    results[,4]<-N*covar2 #Dominance effect
    colnames(results)<-c("NBE", "TICE", "TDCE", "DE")
    rownames(results)<-rownames(mix)}
  
  return(results)}