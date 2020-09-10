library(SummarizedExperiment)
computeConfInt <- function(txiImp, perc = 95)
{
    if(is(txiImp, "SummarizedExperiment"))
    {
        infReps <- assays(txiImp)[grep("infRep",assayNames(txiImp))]
        infReps <- abind::abind(as.list(infReps), along=3)
        if(!("means" %in% assayNames(txiImp)))
            assays(txiImp, withDimnames = F)[["means"]] <- apply(infReps, 1:2, mean)    
        if(!("variance" %in% assayNames(txiImp)))
            assays(txiImp, withDimnames = F)[["variance"]] <- apply(infReps, 1:2, var)
        txiImp <- computeInfRV(txiImp)
    
        N <- dim(infReps)[3]
        k = as.integer(N*(1-perc/100)/2)
        
        assays(txiImp, withDimnames = F)[["LowC"]] <- apply(infReps, 1:2, function(x) quantile(x, probs = c(0.025)))
        assays(txiImp, withDimnames = F)[["HighC"]] <- apply(infReps, 1:2, function(x) quantile(x, probs = c(0.975)))

        return(txiImp)
    }
    infLi <- txiImp[["infReps"]]
    N <- ncol(infLi[[1]])
    k = as.integer(N*(1-perc/100)/2)
    conf <- vector(mode = "list", length(infLi))
    for(i in seq_along(conf))
    {
        conf[[i]] <- matrix(0, nrow = nrow(infLi[[1]]), ncol = 4)
        conf[[i]] <- as.matrix(t(apply(infLi[[i]], 1, function(x)
        {
            sortX <- sort(x)
            c(sortX[k], sortX[length(x)-k+1], mean(x), var(x))
            #c(Rfast::nth(x, k, descending = F), Rfast::nth(x, k, descending = T), mean(x), var(x))
        })))
        dimnames(conf[[i]]) = list(rownames(txiImp[["counts"]]), c("LowC", "HighC", "Mean", "Var"))
        conf[[i]] <- as.data.frame(conf[[i]])
        conf[[i]] <- cbind(conf[[i]], "Width" = conf[[i]][,2] - conf[[i]][,1])
    }
    txiImp$conf <- conf
    return(txiImp)
}

countReads <- function(df)
{
    countDf <- list()
    for(i in seq(1, nrow(df),2))
    {
        ens <- strsplit(strsplit(df[i,1,drop=T], split = ";", fixed= T)[[1]][1], "/", fixed = T)[[1]][2]
        if(ens %in% names(countDf))
            countDf[[ens]] = countDf[[ens]] + 1
        else
            countDf[[ens]] = 1
    }
    return(data.frame(countDf))
}

createPlotDf <- function(matList, indsList = NULL, namesMat = NULL)
{
    if(is.null(indsList))
    {
        indsList <- vector(mode = "list", length(matList))
        names(indsList) <- names(matList)
        for(i in seq_along(matList))
            indsList[[i]] = c(1:nrow(matList[[i]]))
    }
    if(length(matList) != length(indsList))
        stop("Lengths are not same")
    if(is.null(names(matList)) | is.null(names(indsList)))
        stop("Names are null")
    if(sum(names(matList) != names(indsList)) > 0)
        stop("Names are not same")
    
    i = 1
    matList[[i]] <- as.matrix(matList[[i]][indsList[[i]],])
    if(!is.null(namesMat))
    {
        if(ncol(matList[[i]]) != length(namesMat))
            stop("More columns provided")
        colnames(matList[[i]]) <- namesMat
    }
        
    df <- cbind(matList[[i]], type = names(matList)[i])
    
    if(length(indsList) > 1)
    {
        for(i in seq(2, length(indsList)))
        {
            matList[[i]] <- as.matrix(matList[[i]][indsList[[i]],])
            colnames(matList[[i]]) <- namesMat
            df <- rbind(df, cbind(matList[[i]], type = names(matList)[i]))
        }
    }
    
    return(as.data.frame(df))
}

computeInfRep <- function(matList, pc = 5, shift = .01)
{
    infVars <- lapply(matList, function(infMat) matrixStats::rowVars(as.matrix(infMat)))
    infMeans <- lapply(matList, function(infMat) rowMeans(as.matrix(infMat)))
    InfRVs <- lapply(seq_along(matList), function(x) pmax(infVars[[x]] - infMeans[[x]], 0)/(infMeans[[x]] + pc) + shift)
    return(InfRVs)
}

computeInfRV <- function(y, pc=5, shift=.01, meanVariance) {
    if (missing(meanVariance)) {
        meanVariance <- all(c("mean","variance") %in% assayNames(y))
    }
    if (meanVariance) {
        stopifnot(all(c("mean","variance") %in% assayNames(y)))
        infVar <- assays(y)[["variance"]]
        mu <- assays(y)[["mean"]]
    } else {
        infReps <- assays(y)[grep("infRep",assayNames(y))]
        infReps <- abind::abind(as.list(infReps), along=3)
        infVar <- apply(infReps, 1:2, var)
        mu <- apply(infReps, 1:2, mean)
    }
    # the InfRV computation:
    InfRV <- pmax(infVar - mu, 0)/(mu + pc) + shift
    mcols(y)$meanInfRV <- rowMeans(InfRV)
    assays(y, withDimnames = F)[["infRV"]] <- InfRV
    y
}

plotHist <- function(plotDf, var)
{
    library(ggpubr)
    library(dplyr)
    if(!var %in% colnames(plotDf))
        stop("Var not in colnames")
    pList <- list()
    for(t in unique(plotDf$type))
    {
        pList[[t]] <- ggplot(plotDf %>% filter(type == t), aes_string(x = var)) + geom_histogram(bins=200) + labs(title = t)
    }
    figure <- ggarrange(plotlist = pList)
    return(annotate_figure(figure, top = paste("Plotting histogram across", var)))
}

plotViol <- function(plotDf, var, type = "Viol")
{
    library(ggpubr)
    library(dplyr)
    if(!var %in% colnames(plotDf))
        stop("Var not in colnames")
    
    p <- ggplot(plotDf , aes_string(x = "type", y = var))
    if(type == "Viol")
        p <- p + geom_violin()
    if(type == "Box")
        p <- p + geom_boxplot()
    if(type == "Bar")
        p <- p + geom_bar()
    p <- p + labs(title = paste("Plotting Violin Plot across", var))
    return(p)
}

extractBinInds <- function(counts, breaks = 100)
{
    library(binr)
    uniqueCounts <- unique(counts)
    if(is.null(breaks))
        breaks <- length(uniqueCounts)
    binCounts <- bins(uniqueCounts, target.bins = breaks, max.breaks = breaks)
    binNumbers = cut(counts, bins.getvals(binCounts), labels = names(binCounts$binct))
    inds <- lapply(levels(binNumbers), function(b) which(binNumbers == b))
    upLev <- binCounts$binhi
    #gsub("\\]$","",gsub(pattern="\\[\\d+[,]\\s", "", levels(binNumbers))
    names(inds) <- upLev
    return(inds)
}

computeCoverage <- function(counts, infRep, indsList, prop = T)
{
    if(class(indsList) != "list")
        stop("Indexes should be list")
    if(nrow(infRep) != length(counts))
    {
        print(paste(nrow(infRep), length(counts)))
        stop("rows are not same as length")
    }
    
    if(is(infRep, "SummarizedExperiment"))
    {
        highC <- assays(infRep)[["HighC"]]
        lowC <- assays(infRep)[["LowC"]]
        if(prop)
        {
            df <- matrix(0, ncol = ncol(infRep), nrow = length(indsList), dimnames = list(seq(length(indsList)), colnames(infRep)))
            for(i in seq_along(indsList))
            {
                inds = indsList[[i]]
                for(j in seq(ncol(infRep)))
                    df[i,j] <- sum(lowC[inds, j] <= counts[inds] & counts[inds] <= highC[inds, j])/length(inds)
            }
        }
        return(df)
    }
    
    if(prop)
    {
        rat <- sapply(indsList, function(inds)
        {
            sum(infRep[inds, "LowC"] <= counts[inds] & counts[inds] <= infRep[inds, "HighC"])/length(inds)
        })    
    }
    else
    {
        rat <- sapply(indsList, function(inds)
        {
            ifelse(infRep[inds, "LowC"] <= counts[inds] & counts[inds] <= infRep[inds, "HighC"], 1, 0)
        })
    }
    
    names(rat) <- names(indsList)
    return(rat)
}

plotCovDf <- function(covDf, line = F)
{
    library(ggplot2)
    p <- ggplot(covDf, aes(x = log2Counts, y = CI_Coverage, color = Type, shape = Type)) + 
        geom_point() +  scale_x_continuous()
    if(line)
        p <- p + geom_line()
    return(p)
}

createCovDf <- function(confList, counts, cInds)
{
    if(is.null(names(confList)))
        stop("Names cannot be NULL")
    l <- length(confList)
    covDf <- data.frame(Type = rep(names(confList), each = length(cInds)), CI_Coverage = rep(0, l*length(cInds)), 
                        log2Counts = rep(round(log2(as.numeric(names(cInds))), 3), l))
    for(i in seq_along(confList))
        covDf[1:length(cInds) + (i-1)*length(cInds), "CI_Coverage"] <- computeCoverage(counts, confList[[i]], cInds)
    return(covDf)
}

plotBarplot <- function(dfList, indsList, title = F)
{
    library(data.table)
    library(ggplot2)
    if(sum(names(dfList) != names(dfList)) != 0)
        stop("Names dont match")
    i=1
    df <- lapply(seq_along(dfList), function(x) data.frame(table(dfList[[x]][indsList[[x]]]), Type = names(dfList)[x]))
    lengths <- paste(sapply(indsList, length), collapse = " ")
    
    df <- rbindlist(df)
    df <- df[df$Var1==1,]
    p <- ggplot(df, aes(x=Type, y=Freq)) + geom_bar(stat = "identity") 
    if(title)
        p <- p + labs(title = paste("Total transcripts", lengths)) + theme(plot.title = element_text(hjust = 1))
    return(p)
    
}