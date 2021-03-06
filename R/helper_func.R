library(SummarizedExperiment)
library(ggplot2)
library(ggpubr)

getInf <- function(se)
{
    library(SummarizedExperiment)
    infReps <- assays(se)[grep("infRep",assayNames(se))]
    infReps <- abind::abind(as.list(infReps), along=3)   
    return(infReps)
}

appTrueCounts <- function(se, trueCounts)
{
    if(length(trueCounts) != nrow(se))
        stop("Number of transcripts in single cell experiment not same as the length of true counts")
    if(is.null(names(trueCounts)))
        stop("True Counts vector is unnamed")
    
    metadata(se)[["trueCounts"]] <- trueCounts[rownames(se)]
    return(se)
}

##Scaled true counts
computeScaledTrueCounts <- function(se, tCounts = NULL, sf = T, type = "DESeq2", log = F)
{
    if(is.null(tCounts))
    {
        if(!"trueCounts" %in% names(metadata(se)))
            stop("True counts not in se")
        scaledTrueCounts <- matrix(rep(metadata(se)[["trueCounts"]], ncol(se)), ncol = ncol(se), 
            dimnames= list(list(rownames(se)), list(colnames(se))))
    }
    else
    {
        if(is.null(names(tCounts)))
            stop("tCounts cannot have null names")
        tCounts <- tCounts[rownames(se)]
        scaledTrueCounts <- matrix(rep(tCounts, ncol(se)), ncol = ncol(se), 
            dimnames= list(rownames(se), colnames(se)))

    }
        
    if(sf)
    {
        if(!"sf" %in% names(metadata(se)))
            stop("Size Factor needs to be computed")
        metadata(se)[["type"]] <- type
        sizeFac <- metadata(se)[["sf"]][[type]]
        for(i in seq(ncol(scaledTrueCounts)))
            scaledTrueCounts[,i] <- scaledTrueCounts[,i]/sizeFac[1,i]
    }
    if(log)
        scaledTrueCounts <- log(scaledTrueCounts)
    return(scaledTrueCounts)
}

##Scaled Inferential Replicated
computeScaledInfReplicates <- function(se, sf = T, type = "DESeq2", log = F)
{
    infReps <- assays(se)[grep("infRep",assayNames(se))]
    infReps <- abind::abind(as.list(infReps), along=3)
    if(sf)
    {
        if(! "sf" %in% names(metadata(se)))
            stop("Size Factor needs to be computed")
        
        metadata(se)[["type"]] <- type
        sizeFac <- metadata(se)[["sf"]][[type]]
        for(i in seq(dim(sizeFac)[2]))
            infReps[,i,] <- t(t(infReps[,i,])/sizeFac[3:dim(sizeFac)[1],i])
    }
    if(log)
        infReps <- log2(infReps)
    return(infReps)
}

computeSizeFactors <- function(se)
{
    if(!is(se, "SummarizedExperiment"))
        stop("Not a summarized experiment")
    if(!"trueCounts" %in% names(metadata(se)))
        stop("True counts not in metadata")
    
    trueCounts <- metadata(se)[["trueCounts"]]
    infReps <- assays(se)[grep("infRep",assayNames(se))]
    infReps <- abind::abind(as.list(infReps), along=3)
    
    sizeFac <- vector(mode = "list", length = 3)
    names(sizeFac) <- c("DESeq2", "medScale", "depthScale")
    
    for(i in seq_along(sizeFac))
        sizeFac[[i]] <- matrix(1, nrow = dim(infReps)[3] + 2, ncol = ncol(se), 
                                      dimnames = list(c("True", "Sim", dimnames(infReps)[[3]]), colnames(se)))    
    
    
    for(i in seq(ncol(se)))
    {
        print(i)
        sizeFac[["DESeq2"]][,i] <- DESeq2::estimateSizeFactorsForMatrix(cbind(trueCounts, assays(se)[["counts"]][,i], infReps[,i,]))
        
        sfs <- exp(median(log(trueCounts) - log(assays(se)[["counts"]][,i])))
        sfs <- c(sfs, apply(infReps[,i,], 2, function(x) exp(median(log(trueCounts) - log(x)))))
        sizeFac[["medScale"]][2:nrow(sizeFac[[i]]),i] <- sfs
        sizeFac[["medScale"]]["True",i] <- mean(sfs)
        
        sizeFac[["depthScale"]][1,i] <- sum(trueCounts)/sum(assays(se)[["counts"]][,i])
    }

    metadata(se)[["sf"]] <- sizeFac
    return(se)
}

computeConfInt <- function(se, perc = 95, sf = T, type = "DESeq2", log = F)
{
    if(is(se, "SummarizedExperiment"))
    {
        infReps <- assays(se)[grep("infRep",assayNames(se))]
        infReps <- abind::abind(as.list(infReps), along=3)
        if("type" %in% names(metadata(se)))
            metadata(se)[["type"]] <- NULL
        metadata(se)[["log"]] <- log
        sizeFac <- NULL
        if(sf)
        {
            if(! "sf" %in% names(metadata(se)))
                stop("Size Factor needs to be computed")
            
            metadata(se)[["type"]] <- type
            sizeFac <- metadata(se)[["sf"]][[type]]            
            for(i in seq(dim(sizeFac)[2]))
                infReps[,i,] <- t(t(infReps[,i,])/sizeFac[3:dim(sizeFac)[1],i])
            
        }
        if(log)
            infReps <- log(infReps)
        
        assays(se, withDimnames = F)[["mean"]] <- apply(infReps, 1:2, mean)
        assays(se, withDimnames = F)[["median"]] <- apply(infReps, 1:2, median)
        assays(se, withDimnames = F)[["variance"]] <- apply(infReps, 1:2, var)
        
        se <- computeInfRV(se, meanVariance = T)
    
        N <- dim(infReps)[3]
        k = (100-perc)/200
        
        assays(se, withDimnames = F)[["LowC"]] <- apply(infReps, 1:2, function(x) quantile(x, probs = c(k)))
        assays(se, withDimnames = F)[["HighC"]] <- apply(infReps, 1:2, function(x) quantile(x, probs = c(1-k)))
        assays(se, withDimnames = F)[["Width"]] <- assays(se, withDimnames = F)[["HighC"]] - assays(se, withDimnames = F)[["LowC"]]
        if("trueCounts" %in% names(metadata(se)))
        {
            sf <- metadata(se)[["sf"]][[type]][1,]
            cM  <- metadata(se)[["trueCounts"]]
            tC <- matrix(0, length(cM), length(sf))
            for(i in seq_along(sf))
                tC[,i] <- cM/sf[i]
            assays(se, withDimnames = F)[["bias"]] <- tC - assays(se)[["median"]]
        }
        return(se)
    }
    infLi <- se[["infReps"]]
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
        dimnames(conf[[i]]) = list(rownames(se[["counts"]]), c("LowC", "HighC", "Mean", "Var"))
        conf[[i]] <- as.data.frame(conf[[i]])
        conf[[i]] <- cbind(conf[[i]], "Width" = conf[[i]][,2] - conf[[i]][,1])
    }
    se$conf <- conf
    return(se)
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

computeRatio <- function(counts, type = "P", pseudoCount = 1) ##No checking of the rows, assuming there are checked before
    #type stands for paternal or maternal
{
    if(!type %in% c("P", "M"))
        stop("Only P or M needed")
    
    if(is.null(rownames(counts)) & is.null(names(counts)))
        stop("names not provided for txps")
    
    if(!is.null(rownames(counts)))
    {
        tx <- rownames(counts)[1]
        midInds <- nrow(counts)/2
        ratio <- matrix(0, nrow = nrow(counts)/2, ncol = ncol(counts))
    }
    else
    {
        tx <- names(counts)[1]
        midInds <- length(counts)/2
        ratio <- rep(0, length(counts)/2)
    }
    
    txType <- unlist(strsplit(unlist(strsplit(tx, split = "|", fixed = T))[1], "_", fixed=T))[2]
    
    
    reqInds <- c(1:midInds)
    remInds <- c(midInds+1:midInds)
    
    if(txType != type)
    {
        reqInds <- remInds
        remInds <- c(1:midInds)
    }
        
    
    if(length(dim(counts)) > 1)
    {
        for(i in seq(ncol(counts)))
            ratio[,i] <- (counts[reqInds,i] + pseudoCount)/(counts[reqInds,i] + counts[remInds,i] + 2*pseudoCount)
    }
    else
        ratio <- (counts[reqInds] + pseudoCount)/(counts[reqInds] + counts[remInds] + 2*pseudoCount)
    return(ratio)
        
}

computeAlleleConfInt <- function(se, pseudoCount = 1, sf = T, perc = 95, type = "DESeq2", ratioType = "P", log = F)
{
    if(!is(se, "SummarizedExperiment"))
        stop("Summarized Experiment not passed as input")
    ### Checking if the rownames are in the order assumed
    rnames <- rownames(se)
    midInd <- nrow(se)/2
    rnamesM <- sapply(strsplit(rownames(se)[1:midInd], "_", fixed = T), function(x) x[1]) ##1:nrow/2
    rnamesP <- sapply(strsplit(rownames(se)[midInd+1:midInd], "_", fixed = T), function(x) x[1]) ##nrow/2+1:nrow
    if(sum(rnamesM != rnamesP) != 0)
        stop("Maternal and Paternal not following the same order")
    
    infReps <- computeScaledInfReplicates(se, sf = sf, type = type, log = log)
    
    tnames <- rnamesM
    counts <- assays(se)[["counts"]][1:midInd,]
    rownames(counts) <- tnames
    for(i in seq(ncol(counts)))
        counts[,i] <- assays(se)[["counts"]][1:midInd,i] + assays(se)[["counts"]][midInd+1:midInd]
        
    trSe <- SummarizedExperiment(assays = list("counts" = counts), colData = colData(se))
    metadata(trSe)[["log"]] = log
    metadata(trSe)[["type"]] = type
    for(i in seq(dim(infReps)[3]))
        assays(trSe, withDimnames = F)[[paste("infRep", i, sep = "")]] <- computeRatio(infReps[,,i])
    
    infReps <- assays(trSe)[grep("infRep",assayNames(trSe))]
    infReps <- abind::abind(as.list(infReps), along=3)
    
    if(log)
        infReps <- log(infReps)
    
    k = (100-perc)/200
    assays(trSe, withDimnames = F)[["LowC"]] <- apply(infReps, 1:2, function(x) quantile(x, probs = c(k)))
    assays(trSe, withDimnames = F)[["HighC"]] <- apply(infReps, 1:2, function(x) quantile(x, probs = c(1-k)))
    assays(trSe, withDimnames = F)[["Width"]] <- assays(trSe, withDimnames = F)[["HighC"]] - assays(trSe, withDimnames = F)[["LowC"]]
    return(trSe)
}

computeCoverage <- function(counts, se, indsList, trSe = NULL, allele = F, prop = T, pseudoCount = 1, ratioType = "P")
{
    if(class(indsList) != "list")
        stop("Indexes should be list")
    if(nrow(se) != length(counts))
    {
        print(paste(nrow(se), length(counts)))
        stop("rows are not same as length")
    }
    
    if(is(se, "SummarizedExperiment"))
    {
        sf <- rep(1, ncol(se)) ##Sizefactor for true counts
        type <- NULL
        log <- F
        if("type" %in% names(metadata(se)))
        {
            type <- metadata(se)[["type"]]
            sf <- metadata(se)[["sf"]][[type]][1,]
            log <- metadata(se)[["log"]]
        }
        sfFlag = ifelse(is.null(type), F, T)
        scTVals <- computeScaledTrueCounts(se, tCounts = counts, sf = sfFlag, type = type, log = log)
        #print(head(scTVals))
        if(allele)
        {
            scTVals <- computeRatio(scTVals, type = ratioType, pseudoCount = pseudoCount)
            
            if(is.null(trSe))
                trSe <- computeAlleleConfInt(se, pseudoCount = pseudoCount, sf = sfFlag, type = type, log = log, ratioType = ratioType)
            
            highC <- assays(trSe)[["HighC"]]
            lowC <- assays(trSe)[["LowC"]]
        }
        else
        {
            if(!"LowC" %in% assayNames(se))
                stop("LowC not in assayNames")
            highC <- assays(se)[["HighC"]]
            lowC <- assays(se)[["LowC"]]
        }
        df <- matrix(0,1,1)
        covInds <- list()
        if(prop)
            df <- matrix(0, ncol = ncol(se), nrow = length(indsList), dimnames = list(seq(length(indsList)), colnames(se)))
        else
        {
            covInds <- vector(mode = "list", length(indsList))
            names(covInds) <- names(indsList)
            for(i in seq_along(covInds))
            {
                covInds[[i]] <- vector(mode = "list", length = ncol(se))
                names(covInds[[i]]) <- colnames(se)
            }    
        }
        
        for(i in seq_along(indsList))
        {
            inds = indsList[[i]]
            for(j in seq(ncol(se)))
            {
    
                if(prop)
                    df[i,j] <- sum(lowC[inds, j] <= scTVals[inds,j] & scTVals[inds,j] <= highC[inds, j])/length(inds)
                else
                    covInds[[i]][[j]] <- which(lowC[inds, j] <= scTVals[inds,j] & scTVals[inds,j] <= highC[inds, j])
            }
        }
        if(prop)
            return(df)
        else
            return(covInds)
    }
    
    # if(prop)
    # {
    #     rat <- sapply(indsList, function(inds)
    #     {
    #         sum(infRep[inds, "LowC"] <= counts[inds] & counts[inds] <= infRep[inds, "HighC"])/length(inds)
    #     })    
    # }
    # else
    # {
    #     rat <- sapply(indsList, function(inds)
    #     {
    #         ifelse(infRep[inds, "LowC"] <= counts[inds] & counts[inds] <= infRep[inds, "HighC"], 1, 0)
    #     })
    # }
    # 
    # names(rat) <- names(indsList)
    # return(rat)
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
    #gsub("\\]$","",gsub(pattern="\\[\\d+[,]\\s", "", levels(binNumbers))
    names(inds) <- sapply(inds, function(x) as.character(max(counts[x])))
    return(inds)
}

plotCovDf <- function(covDf, line = F, ratio = F)
{
    library(ggplot2)
    if(ratio)
        p <- ggplot(covDf, aes(x = Ratio, y = Coverage, color = Type, shape = Type)) + 
            geom_point() +  scale_x_continuous()
    else
        p <- ggplot(covDf, aes(x = log2Counts, y = Coverage, color = Type, shape = Type)) + 
            geom_point() +  scale_x_continuous()
    if(line)
        p <- p + geom_line()
    p <- p + theme_grey(base_size = 12) + theme(legend.position = "bottom") 
    return(p)
}

plotDensity <- function(se, transcript, cols = NULL)
{
    library(ggplot2)
    inf <- getInf(se)
    if(!transcript %in% dimnames(se)[[1]])
        stop("Transcript not present in dimnames")
    df <- data.frame(t(inf[transcript,,]))
    if(!is.null(cols))
    {
        if(ncol(df) != length(cols))
            stop("Columns are not equal")
        colnames(df) <- cols
    }
    
    df <- reshape::melt(df)
    
    colnames(df) <- c("Type", "infCounts")
    countT <- data.frame(counts = assays(se)[["counts"]][transcript,], type = unique(df[,"Type"]))
    
    p <- ggplot(df, aes(x = infCounts, color=Type)) + geom_density() +
        geom_vline(data = countT, aes(xintercept=counts, color = type), linetype = "dashed", size = 1)
    return(p)
}

plotCounts <- function(se, trueCounts, log = T, inds = 1)
{
    library(ggplot2)
    df <- data.frame(cbind(estCounts = assays(se)[["counts"]][,inds], trueCounts = trueCounts))
    tit <- ""
    if(log)
    {
        tit <- "log2 "
        df <- log2(df+1)
    }
    
    p <- ggplot(df, aes(x = estCounts, y = trueCounts)) +
        geom_point() + geom_abline(slope = 1, intercept = 0) +
         xlab(paste(tit, "Estimated Counts", sep = "")) + 
         ylab(paste(tit, "True Counts", sep = ""))
    return(p)    
}

createCovDf <- function(confList, counts, cInds, cols = NULL, logC = T, allele = F, trSe = NULL, ratio = F)
{
    library(reshape)
    if(is(confList, "SummarizedExperiment"))
    {
        
        df <- computeCoverage(counts, confList, cInds, allele = allele, trSe = trSe)
        if(!is.null(cols))
        {
            if(length(cols) != ncol(df))
                stop("Not enough columns")
            colnames(df) <- cols
        }
        counts = as.numeric(names(cInds))
        if(logC)
            counts=log2(counts)
        dfMelt <- melt(df)
        dfMelt <- cbind(dfMelt, rep(counts, ncol(df)))
        colnames(dfMelt) <- c("inds", "Type", "Coverage", "log2Counts")
        if(ratio)
            colnames(dfMelt)[4] <- "Ratio"
        return(dfMelt)
        
    }
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

createWidthDf <- function(se, counts = NULL)
{
    if(is.null(counts))
    {
        if(! "trueCounts" %in% names(metadata(se)))
            stop("True counts not in the summ Experiment and also not provided by the user")
        counts <- metadata(se)[["trueCounts"]]
    }
    sfT <- metadata(se)[["sf"]][["DESeq2"]][1,] #size factor for true counts with Boots being first followed by GS and/or Polee
    df <- data.frame(widthB = assays(se)[["Width"]][,1], widthGS = assays(se)[["Width"]][,2])
    df[["wiBProp"]] = df[["widthB"]]/(counts/sfT[1])
    df[["wiGSProp"]] = df[["widthGS"]]/(counts/sfT[2])
    df[["BCov"]] <- ifelse(assays(se)[["LowC"]][,1] <= counts/sfT[1] & counts/sfT[1] <= assays(se)[["HighC"]][,1], T, F)
    df[["GSCov"]] <- ifelse(assays(se)[["LowC"]][,2] <= counts/sfT[2] & counts/sfT[2] <= assays(se)[["HighC"]][,2], T, F)
    
    if(dim(se)[2] == 3)
    {
        df <- cbind(df, widthPolee = assays(se)[["Width"]][,3])
        df[["wiPoleeCount"]] = df[["widthPolee"]]/(counts/sfT[3])
        df[["PoleeCov"]] <- ifelse(assays(se)[["LowC"]][,3] <= (counts/sfT[3]) & (counts/sfT[3]) <= assays(se)[["HighC"]][,3], T, F)
        df <- df[,c(1:2,7,3:4,8,5:6,9)]
    }
    df <- cbind(df, counts = counts)
    df <- df[order(df[["counts"]]),]
    return(df)
}

plotWidthDf <- function(df, widthCols = c(1:2), widthPropCols = c(3:4), log = T, hex = T)
{
    if(log)
        df[,c(widthCols, widthPropCols, ncol(df))] <- log2(df[,c(widthCols, widthPropCols, ncol(df))] + 1)
    
    pWidth <- vector(mode="list", length(widthCols))
    j=1
    
    for(i in widthCols)
    {
        pWidth[[j]] <- ggplot(df, aes_string(x="counts", y =colnames(df)[i]))
        if(hex)
            pWidth[[j]] <- pWidth[[j]] + geom_hex()
        else
            pWidth[[j]] <- pWidth[[j]] + geom_point()
        pWidth[[j]] <- pWidth[[j]] + geom_smooth(method = "lm") + 
            xlab(paste("log2", colnames(df)[i])) + ylab(paste("log2", "counts"))
        j = j + 1
    }
    ncol=2
    nrow=1
    if(length(widthCols) > 2)
        nrow=2
    pWidth <- ggarrange(plotlist = pWidth, ncol = ncol, nrow = nrow)
    
    pWidthProp <- list()
    j=1
    for(i in widthPropCols)
    {
        pWidthProp[[j]] <- ggplot(df, aes_string(x = "counts", y = colnames(df)[i]))
        if(hex)
            pWidthProp[[j]] <- pWidthProp[[j]]  + geom_hex()
        else
            pWidthProp[[j]] <- pWidthProp[[j]]  + geom_point()
            
        pWidthProp[[j]] <- pWidthProp[[j]] + geom_smooth(method = "lm") + 
            ylab(paste("log2", colnames(df)[i])) + xlab(paste("log2", "counts"))
        j=j+1
    }
    pWidthProp <- ggarrange(plotlist = pWidthProp, nrow =nrow, ncol=ncol)
    
    return(list(pWidth, pWidthProp))
}

plotSummary <- function(se, tCounts = NULL, type = "Scatter", col_inds = c(1,2), summQuant = "mean", nbreaks = 50, log=T)
{
    if(!is.null(tCounts))
        cInds <- extractBinInds(tCounts, breaks = nbreaks)
    else
        cInds <- extractBinInds(assays(se)[["counts"]][,1], breaks = nbreaks)
    pList <- vector(mode = "list", length(cInds))
    summDf <- data.frame(assays(se)[[summQuant]][,col_inds])
    if(summQuant == "bias")
        summDf <- abs(summDf)
    for(i in seq_along(cInds))
    {
        inds = cInds[[i]]
        if(!summQuant %in% assayNames(se))
            stop(paste(summQuant, "does not exist"))
        if(log)
            df <- log2(summDf[inds,]+1)
        if(type == "Scatter")
        {
            pList[[i]] <- ggplot(log2(df+1), aes_string(x=colnames(df)[1], y=colnames(df)[2])) + 
            geom_point() + geom_smooth(method = "lm", color = "red") +
            geom_abline(intercept = 0, slope = 1, color = "blue") + 
            ggtitle(round(log2(as.numeric(names(cInds)[i])),2))
        }
        else
        {
            df <- reshape::melt(df)
            pList[[i]] <- ggplot(df, aes(x=variable, y=value)) + geom_boxplot() +
                ggtitle(round(log2(as.numeric(names(cInds)[i])),2))
        }
    }
    return(pList)
}