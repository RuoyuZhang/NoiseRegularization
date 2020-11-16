# helper function
library(dplyr)
library(tibble)
#library(data.table)

# permute matrix by each column
shuffle = function(m,seed = 1){
    ncell = nrow(m)
    ngene = ncol(m)
    set.seed(seed)
    for (i in 1:ngene){
        m[,i] = m[,i][sample(1:ncell,ncell)]
    }
    return(m)
}


# calculate quantile matrix
calculate_quantile_matrix = function(m,expr.gene.list,cell.names=NULL,quantile.list){
    if (is.null(cell.names)){cell.names = names(expr.gene.list)}
    
    quantile.mat = NULL
    for (cell in cell.names){
        expr.genes = m[cell,expr.gene.list[[cell]]]
        if (min(expr.genes) < 0){
          expr.genes = expr.genes + abs(min(expr.genes))
        }
        cell.quantile = quantile(expr.genes,probs = quantile.list)
        quantile.mat = cbind(quantile.mat, cell.quantile)
    }
    colnames(quantile.mat) = cell.names
    return(quantile.mat)
}


# calculate quantile matrix gene by gene
calculate_quantile_matrix.bygene = function(m,expr.cell.list=NULL,gene.names=NULL,quantile.list){
    if (is.null(gene.names)){gene.names = colnames(m)}
    
    quantile.mat = NULL
    for (gene in gene.names){
        if (is.null(expr.cell.list)){
          expr.cell = m[,gene]
        }else{
          expr.cell = m[expr.cell.list[[gene]],gene]
        }
        if (min(expr.cell) < 0){
          expr.cell = expr.cell + abs(min(expr.cell))
        }
        gene.quantile = quantile(expr.cell,probs = quantile.list)
        quantile.mat = cbind(quantile.mat, gene.quantile)
    }
    colnames(quantile.mat) = gene.names
    return(quantile.mat)
}



# add random number to the matrix
add.random = function(m,min,max,seed = 1){
    max = as.numeric(max)
    if (max < 0){
      min.temp = min
      min = max
      max = min.temp
    }
    #print(max)
    ncell = nrow(m)
    ngene = ncol(m)
    set.seed(seed)
    for (i in 1:ngene){
        m[,i] = m[,i]+runif(ncell,min,max)
    }
    return(m)
}


# add random number to the matrix gene by gene
add.random.bygene = function(m,min,max.vec,max.low=1,seed = 1){
    ncell = nrow(m)
    ngene = ncol(m)
    gene.names = colnames(m)
    set.seed(seed)
    for (gene in gene.names){
        max = max(max.vec[gene],max.low)
        m[,gene] = m[,gene]+runif(ncell,min,max)
    }
    return(m)
}


# calculate correlation
cal.cor = function(m,method){
    m.cor = cor(m,method = method)
    m.cor.pairwise = lower.tri(m.cor,diag = F)
    m.cor.pairwise = m.cor[m.cor.pairwise]
    return(m.cor.pairwise)
}


# calculate pairwise correlation
pair_cor = function(m,gene.list,ncore=4,nblocks=20){
    print('start')
    print(Sys.time())
    gene.list = sort(gene.list)
    m.take = m[,gene.list]
    n = length(gene.list)
    message(Sys.time())
    message(n)
    # calculate correlation
    # m.cor.spearman = cor(m.take,method = 'spearman',nThreads=4)
    m.cor = bigcorPar(m.take, nblocks = nblocks, verbose = F, ncore=ncore)
    m.cor.spearman = m.cor[['corMAT']]
    m.cor.p = m.cor[['pvalMAT']]

    #m.cor.pearson = cor(m.take,method = 'pearson')
    
    #print(m.cor.spearman[1:5,1:5])
    #print(m.cor.pearson[1:5,1:5])
    message(Sys.time())
    message('finish correlation')
    
    # create naming matrix
    df.name = expand.grid(gene.list,gene.list)
    #df.name$Var3 = paste0(df.name$Var1,'_',df.name$Var2)
    df.name.vec = paste0(df.name$Var1,'_',df.name$Var2)

    #m.name = matrix(df.name$Var3,nrow = n)
    m.name = matrix(df.name.vec,nrow = n)
    
    message(Sys.time())
    message('finish renaming')
    
    m.cor.spearman.vec = m.cor.spearman[upper.tri(m.cor.spearman,diag = F)]
    m.cor.p.vec = m.cor.p[upper.tri(m.cor.p,diag = F)]

    #m.cor.pearson.vec = m.cor.pearson[upper.tri(m.cor.pearson,diag = F)]
    m.cor.names = m.name[upper.tri(m.name,diag = F)]
    
    print(Sys.time())
    print('finish taking vector')
    
    cor.df = tibble(pair_name = m.cor.names, spearman = m.cor.spearman.vec, p = m.cor.p.vec)
    
    #message(Sys.time())
    #message('start sorting')
    
    #cor.df = cor.df[order(cor.df$spearman, decreasing = T),]
    #message(Sys.time())
    #message('finish sorting')
    
    return(cor.df)
}


# calculate the fraction of correlation in ppi 
ppi_frac_plot = function(cor.df,cuts,string.data,title){
    
    in.ppi.frac.vec = c()
    in.ppi.num = c()
    for (cut in cuts){
        cor.df.sub = cor.df[abs(as.numeric(cor.df[,'spearman']))>=cut,]
        in.ppi = length(intersect(cor.df.sub[,'pair_name'], string.data$pair_name))
        in.ppi.frac = in.ppi/nrow(cor.df.sub)
        in.ppi.frac.vec = c(in.ppi.frac.vec,in.ppi.frac)
        in.ppi.num = c(in.ppi.num,paste0(in.ppi,'_',nrow(cor.df.sub)))
    }
    fraction.df = data.frame(cor.cut = cuts, in.ppi.frac = in.ppi.frac.vec,in.ppi.num=in.ppi.num)
    print(ggplot(fraction.df,aes(cor.cut, in.ppi.frac))+geom_point()+ ggtitle(title)+ ylim(0,1)+
          geom_text_repel(aes(label = in.ppi.num)) )
    #return(fraction.df)
}


# calculate the fraction of correlation in ppi by rank of gene pair
ppi_frac_plot_rank = function(cor.df,cuts,string.data,title,type,do.plot=T){
    cor.df = cor.df[order(cor.df$spearman, decreasing = T),]
    
    in.ppi.frac.vec = c()
    in.ppi.num = c()
    for (cut in cuts){
        cor.df.sub = cor.df[1:cut,]
        in.ppi = length(intersect(cor.df.sub$pair_name, string.data$pair_name))
        in.ppi.frac = in.ppi/nrow(cor.df.sub)
        in.ppi.frac.vec = c(in.ppi.frac.vec,in.ppi.frac)
        in.ppi.num = c(in.ppi.num,paste0(in.ppi,'_',nrow(cor.df.sub)))
    }
    fraction.df = tibble(cor.cut = cuts, in.ppi.frac = in.ppi.frac.vec,in.ppi.num=in.ppi.num,type = type)
    if (do.plot){
      print(ggplot(fraction.df,aes(cor.cut, in.ppi.frac))+geom_point()+ ggtitle(title)+ ylim(0,1)+
          geom_text_repel(aes(label = in.ppi.num)) )
    }
    return(fraction.df)
}


# plot ppi plot with different random number
ppi_frac_plot_dif_quantile = function(m,gene.list, cuts = seq(100,10000,500), quantile.mat, quantile.colname, quantile.list,string.data,title,out.dir,seed=1){
  #add.cor.df = pair_cor(m=m, gene.list = gene.list)
  #frac.df.all = ppi_frac_plot_rank(cor.df = add.cor.df, cuts = cuts,string.data = string.data,title = '',type =0,do.plot = F)
  frac.df.all = NULL
  for (q in quantile.list){
      m.use=add.random(m[,gene.list],min=0,max=quantile.mat[q,quantile.colname],seed=seed)
      add.cor.df = pair_cor(m=m.use,gene.list = gene.list)
      frac.df = ppi_frac_plot_rank(cor.df = add.cor.df, cuts = cuts,string.data = string.data,title = '',type = q,do.plot = F)
      frac.df.all = rbind(frac.df.all,frac.df)
      q = sub('%','p',q)
      saveRDS(frac.df,file=paste0(out.dir,'/',title,'.',q,'.frac.df.rds'))
      saveRDS(add.cor.df,file=paste0(out.dir,'/',title,'.',q,'.cor.df.rds'))
  }
  
  p=ggplot(frac.df.all, aes(cor.cut,in.ppi.frac,col = factor(type,levels = c(0,quantile.list)))) + geom_point() + geom_path() + scale_color_discrete(name = '') +
    ggtitle(title)
  #print(p)
  saveRDS(frac.df.all,file=paste0(out.dir,'/',title,'.all.frac.df.rds'))
  return(list(df=frac.df.all,p=p))
}


# plot ppi plot with different random number (gene by gene)
noise.regularization = function(m,gene.list, cuts = seq(100,10000,500), max.low=1, quantile.mat, quantile.list,string.data,title,out.dir,ncore=4,nblocks=20,seed=1){
  #add.cor.df = pair_cor(m=m, gene.list = gene.list)
  #frac.df.all = ppi_frac_plot_rank(cor.df = add.cor.df, cuts = cuts,string.data = string.data,title = '',type =0,do.plot = F)
  dir.create(out.dir,showWarnings = F,recursive = T)
  frac.df.all = NULL
  for (q in quantile.list){
    print(seed)
      if (q != 'add0'){
        m.use=add.random.bygene(m[,gene.list],min=0,max=quantile.mat[q,gene.list],max.low=max.low,seed=seed)
      } else{
        m.use = m
      }
      add.cor.df = pair_cor(m=m.use,gene.list = gene.list,ncore=ncore,nblocks=nblocks)
      #frac.df = ppi_frac_plot_rank(cor.df = add.cor.df, cuts = cuts,string.data = string.data,title = '',type = q,do.plot = F)
      #frac.df.all = rbind(frac.df.all,frac.df)
      q = sub('%','p',q)
      #saveRDS(frac.df,file=paste0(out.dir,'/',title,'.',q,'.frac.df.rds'))
      saveRDS(add.cor.df,file=paste0(out.dir,'/',title,'.',q,'.cor.df.rds'))
      remove(m.use)
      remove(add.cor.df)
  }
  
  #p=ggplot(frac.df.all, aes(cor.cut,in.ppi.frac,col = factor(type,levels = c(0,quantile.list)))) + geom_point() + geom_path() + scale_color_discrete(name = '') +
  #  ggtitle(title)
  #print(p)
  #saveRDS(frac.df.all,file=paste0(out.dir,'/',title,'.all.frac.df.rds'))
  #return(list(df=frac.df.all))
}

# calculate and plot pairwise correlation
pairwise.cor = function(gene.list,do.shuffle=F,do.random=F,seed=1,title='',min=0,max=1,
                       umi = umi, NormUMI=NormUMI, logNormUMI=logNormUMI, scale.data = scale.data,
                       nbr=nbr,nb=nb,poisson = poisson, magic = magic, dca = dca, do.return = F, do.plot=T,method='spearman'){
  
  umi=umi[,gene.list]
  NormUMI=NormUMI[,gene.list]
  logNormUMI=logNormUMI[,gene.list]
  scale.data=scale.data[,gene.list]
  nbr=nbr[,gene.list]
  nb=nb[,gene.list]
  poisson=poisson[,gene.list]
  magic = magic[,gene.list]
  dca = dca[,gene.list]
  
  #print(sum(umi[,1]))
  #print(sum(nbr[,1]))
  
  # adding a random number
  if (do.random){
    if (length(max) == 1){
        max = rep(max, 9)
        names(max) = c('umi','NormUMI','logNormUMI','scale.data','nbr','nb','poisson','magic','dca')
    }
    
    umi = add.random(umi,min,max['umi'],seed)
    NormUMI = add.random(NormUMI,min,max['NormUMI'],seed)
    logNormUMI = add.random(logNormUMI,min,max['logNormUMI'],seed)
    scale.data = add.random(scale.data,min,max['scale.data'],seed)
    nbr = add.random(nbr,min,max['nbr'],seed)
    nb = add.random(nb,min,max['nb'],seed)
    poisson = add.random(poisson,min,max['poisson'],seed)
    magic = add.random(magic,min,max['magic'],seed)
    dca = add.random(dca,min,max['dca'],seed)
  }  
  
    
  if (do.shuffle){
    umi = shuffle(umi,seed)
    NormUMI = shuffle(NormUMI,seed)
    logNormUMI = shuffle(logNormUMI,seed)
    scale.data = shuffle(scale.data,seed)
    nbr = shuffle(nbr,seed)
    nb = shuffle(nb,seed)
    poisson = shuffle(poisson,seed)
    magic = shuffle(magic,seed)
    dca = shuffle(dca,seed)
  }
  
  #print(sum(umi[,1]))
  #print(sum(nbr[,1]))
  
  umi.cor.pairwise = cal.cor(umi,method = method)
  NormUMI.cor.pairwise = cal.cor(NormUMI,method = method)
  logNormUMI.cor.pairwise = cal.cor(logNormUMI,method = method)
  scale.cor.pairwise = cal.cor(scale.data,method = method)
  nbr.cor.pairwise = cal.cor(nbr,method = method)
  nb.cor.pairwise = cal.cor(nb,method = method)
  poisson.cor.pairwise = cal.cor(poisson,method = method)
  magic.cor.pairwise = cal.cor(magic,method = method)
  dca.cor.pairwise = cal.cor(dca,method = method)
  
  
  # data.frame
  umi.cor.df = data.frame(cor = umi.cor.pairwise, type = 'UMI')
  NormUMI.cor.df = data.frame(cor = NormUMI.cor.pairwise, type = 'Norm UMI')
  logNormUMI.cor.df = data.frame(cor = logNormUMI.cor.pairwise, type = 'logNorm UMI')
  scale.cor.df = data.frame(cor = scale.cor.pairwise, type = 'Scale Data')
  nbr.cor.df = data.frame(cor = nbr.cor.pairwise, type = 'NBR')
  nb.cor.df = data.frame(cor = nb.cor.pairwise, type = 'NB')
  poisson.cor.df = data.frame(cor = poisson.cor.pairwise, type = 'Poisson')
  magic.cor.df = data.frame(cor = magic.cor.pairwise, type = 'MAGIC')
  dca.cor.df = data.frame(cor = dca.cor.pairwise, type = 'DCA')

  cor.df = rbind(umi.cor.df, NormUMI.cor.df, logNormUMI.cor.df,scale.cor.df,
                    nbr.cor.df,nb.cor.df,poisson.cor.df,magic.cor.df,dca.cor.df)
  
  p = ggplot(cor.df,aes(type,cor,col=type)) + geom_violin()+ geom_boxplot(width=0.1) + ggtitle(title) + geom_hline(yintercept=c(0.2,0.4))
  if (do.plot){
    return(p)
  }
  
  if (do.return){
    return(cor.df)
  }
  
}




# plot histogram
plot_multi_histogram <- function(df, feature, label_column) {
    
    plt <- ggplot(df, aes(x=eval(parse(text=feature)), fill=eval(parse(text=label_column)))) +
            geom_histogram(breaks = seq(-1,1,0.05),aes(y = stat(width*density)), alpha = 0.3, color = 'black') +
              facet_wrap(~ type,nrow=2) + scale_y_continuous(labels = percent_format()) + 
                geom_hline(yintercept = 0.1)+ geom_vline(xintercept = 0,linetype = "dashed") +
              theme_bw() + guides(fill=guide_legend(title=label_column))
    print(plt)
}


# plot single gene pairs
plot_gene_pair = function(type,m,g1,g2){
    v1 = m[,g1]
    v2 = m[,g2]
    df = data.frame(x=v1, y=v2)
    cor.value1 = signif(cor(v1,v2,method='spearman'), digits = 3)
    cor.value2 = signif(cor(v1,v2,method='pearson'), digits = 3)
    #title = paste(type, g1, g2, '\nSpearman correlation:', cor.value1, '\nPearson corrleation:',cor.value2)
    title = paste(type, g1, g2, '\nSpearman correlation:', cor.value1)

    p = ggplot(df,aes(x,y)) + geom_point(alpha=0.1) + ggtitle(title) +
        xlab(g1)+ylab(g2)
    return(p)
}


# cor to p value
# borrow some code from psych::corr.test, but much faster than corr.test
cor2pvalue = function(r, n) {
  t <- (r*sqrt(n-2))/sqrt(1-r^2)
  p <- 2*(1 - pt(abs(t),(n-2)))
  se <- sqrt((1-r*r)/(n-2))
  out <- list(r, n, t, p, se)
  names(out) <- c("r", "n", "t", "p", "se")
  return(out)
}



# 
# multi core from https://gist.github.com/bobthecat/5024079
# did some modifications
bigcorPar <- function(x, nblocks = 10, verbose = TRUE, ncore="all", ...){
  require(doMC)
  require(psych)
	if(ncore=="all"){
		ncore = parallel::detectCores()
    #cl <- makeCluster(ncore)
    #registerDoParallel(cl)
		registerDoMC(cores = ncore)
	} else{
		#cl <- makeCluster(ncore)
    #registerDoParallel(cl)
    registerDoMC(cores = ncore)

	}

	NCOL <- ncol(x)
  print(NCOL)

	## test if ncol(x) %% nblocks gives remainder 0
	#if (NCOL %% nblocks != 0){stop("Choose different 'nblocks' so that ncol(x) %% nblocks = 0!")}

	## preallocate square matrix of dimension
	## ncol(x) in 'ff' single format
	#corMAT <- ff(vmode = "single", dim = c(NCOL, NCOL))
  corMAT = matrix(rep(0,NCOL*NCOL),ncol = NCOL)
  gene.names = colnames(x)
  colnames(corMAT) = gene.names
  rownames(corMAT) = gene.names
  
  pvalMAT = corMAT
  
	## split column numbers into 'nblocks' groups
  each_block = ceiling(NCOL/nblocks)
  left = ifelse(NCOL%%nblocks == 0, each_block, NCOL - each_block*(nblocks-1) )
	SPLIT <- split(1:NCOL, c(rep(1:(nblocks-1), each = each_block), rep(nblocks, left)))
  #print(SPLIT)

	## create all unique combinations of blocks
	COMBS <- expand.grid(1:length(SPLIT), 1:length(SPLIT))
	COMBS <- t(apply(COMBS, 1, sort))
	COMBS <- unique(COMBS)

	## iterate through each block combination, calculate correlation matrix
	## between blocks and store them in the preallocated matrix on both
	## symmetric sides of the diagonal
  gc()
	results <- foreach(i = 1:nrow(COMBS)) %dopar% {
		COMB <- COMBS[i, ]
		G1 <- SPLIT[[COMB[1]]]
		G2 <- SPLIT[[COMB[2]]]
		if (verbose) cat("Block", COMB[1], "with Block", COMB[2], "\n")
		flush.console()
		COR.cor <- cor(x[, G1], x[, G2], method = 'spearman', ...)
    #COR <- corr.test(x[, G1], x[, G2], method = 'spearman', adjust = 'none', ci = FALSE, ...)
    
    #r <- cor(test.data, use = "pairwise.complete.obs", method= 'spearman') # MUCH MUCH faster than corr.test()
    n <- t(!is.na(x[, G1])) %*% (!is.na(x[, G2])) # same as count.pairwise(x,y) from psych package  
    
    COR = cor2pvalue(COR.cor,n)
    
    COR.cor = COR$r
    COR.p = COR$p
    
    #print(COR[1:5,1:5])
    
    rownames(COR.cor) = gene.names[G1]
    colnames(COR.cor) = gene.names[G2]
    
    rownames(COR.p) = gene.names[G1]
    colnames(COR.p) = gene.names[G2]
    
    gc()
    return(list(cor = COR.cor, p = COR.p))
	}
  

  while(sum(sapply(results, is.null)) >0){
    print('rerun')

    gc()
    results <- foreach(i = 1:nrow(COMBS)) %dopar% {
      COMB <- COMBS[i, ]
      G1 <- SPLIT[[COMB[1]]]
      G2 <- SPLIT[[COMB[2]]]
      if (verbose) cat("Block", COMB[1], "with Block", COMB[2], "\n")
      flush.console()
      COR.cor <- cor(x[, G1], x[, G2], method = 'spearman', ...)
    
    n <- t(!is.na(x[, G1])) %*% (!is.na(x[, G2])) # same as count.pairwise(x,y) from psych package  
    
    COR = cor2pvalue(COR.cor,n)
    
    COR.cor = COR$r
    COR.p = COR$p
    
    #print(COR[1:5,1:5])
    
    rownames(COR.cor) = gene.names[G1]
    colnames(COR.cor) = gene.names[G2]
    
    rownames(COR.p) = gene.names[G1]
    colnames(COR.p) = gene.names[G2]
    
    gc()
    return(list(cor = COR.cor, p = COR.p))
    }
}
  
  print('block done')
  print(length(results))
  for (small in 1:length(results)){
    COR = results[[small]]
    G1 = rownames(COR[['cor']])
    G2 = colnames(COR[['cor']])
    corMAT[G1,G2] = COR[['cor']]
    corMAT[G2, G1] <- t(as.matrix(COR[['cor']]))
    
    pvalMAT[G1,G2] = COR[['p']]
    pvalMAT[G2, G1] <- t(as.matrix(COR[['p']]))
  }
  
  #print(results)
	gc()
	return(list(corMAT = corMAT, pvalMAT = pvalMAT))
}


# adding expression level to each pair
calculate.mean.pair = function(df, gene.mean.expr){
    #df = logNormUMI.cor.add0[1:10,]
    gene.names = do.call(rbind,(strsplit(as.character(df[,'pair_name']),split = '_')))
    df = cbind(df,gene.names)
    colnames(df)[3:4] = c('gene1','gene2')
    gene1.expr = gene.mean.expr[as.character(df$gene1)]
    gene2.expr = gene.mean.expr[as.character(df$gene2)]
    gene.mean.expr = sqrt(gene1.expr * gene2.expr)
    df = cbind(df,gene1.expr,gene2.expr,gene.mean.expr)
    return(df)
}


# compile function
library('compiler')
noise.regularization.cmp = cmpfun(noise.regularization)


### pick max correlation among clusters
pick.max.cluster = function(cor.dir,prefix,add='add0',clusters = 0:9, cell.num, cell.per.cut = 0.01, min.cell.cut = 50,gene.count,prefix2=NULL){
  require(matrixStats)
  require(tibble)
  
  if(is.null(prefix2)){prefix2=prefix}
  #all.cor.df = readRDS(paste0(cor.dir,'/',prefix,'/all/', prefix2,'.bygene.',add,'.cor.df.rds'))
  #this.cor.df = readRDS(paste0(cor.dir,'/',prefix,'/cluster',clusters[1], '/', prefix2,'.bygene.',add,'.cor.df.rds'))
  #cor.df = tibble(pair_name = this.cor.df$pair_name,all = this.cor.df$spearman)
  
  cor.df = NULL
  for (i in clusters){
      print(i)
      this.cor.df = readRDS(paste0(cor.dir,'/','/cluster',i, '/', prefix2,'.bygene.',add,'.cor.df.rds'))
      
      if (i==0){cor.df = tibble(pair_name = this.cor.df$pair_name)}

      if(!identical(cor.df$pair_name,this.cor.df$pair_name)){message('pair name not match!')}
      #gene.cluster = colnames(gene.expr.per.by.cluster)[gene.expr.per.by.cluster[as.character(i),]>0.002]
      gene.index = gene.count[as.character(i),] > max(min.cell.cut,cell.num[as.character(i)]*cell.per.cut)
      print(sum(gene.index))
      gene.cluster = colnames(gene.count)[gene.index]
  
      gene.cluster = sort(gene.cluster)
      take.pair.name = expand.grid(gene.cluster,gene.cluster)
      #df.name$Var3 = paste0(df.name$Var1,'_',df.name$Var2)
      take.pair.name.vec = paste0(take.pair.name$Var1,'_',take.pair.name$Var2)
      new.value = this.cor.df$spearman
      new.value[!(this.cor.df$pair_name %in% take.pair.name.vec)] = -2
      
      this.col = paste0('c',i)
      if (this.col %in% colnames(cor.df)){
          cor.df= cor.df %>% select(-!!(this.col))
      }
      
      cor.df = add_column(cor.df, !!(this.col) := new.value)
      remove(this.cor.df)
  }
  
  cor.df = as.data.frame(cor.df,row.names = cor.df$pair_name) 
  rownames(cor.df) = cor.df$pair_name
  
  cluster.col = paste0('c',clusters)
  max.col = rowMaxs(as.matrix(cor.df[,cluster.col]))
  cor.df = add_column(cor.df, max = max.col)
  
  return(cor.df)
}


