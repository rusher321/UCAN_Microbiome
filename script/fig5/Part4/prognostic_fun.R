########## ---------------- pro ---------------- ###########

HR <- function(data,x,time,censor,adj=NULL,outputAll=FALSE){
  data <- data[!is.na(data[,x]),,drop=F]
  if(is.factor(data[,x])){
    data[,x] <- droplevels(data[,x])
  }else if(is.numeric(data[,x])){
    data[,x] <- scale(data[,x])
  }
  lab.x <- paste(x, collapse = "+")
  formula <- paste0("Surv(", time,",", censor,")~",lab.x)
  if(!is.null(adj)){
    lab.adj <- paste(adj, collapse = "+")
    formula <- paste0(formula,"+",lab.adj)
  }
  formula <- as.formula(formula)
  #print(formula)
  #x <- unlist(strsplit(as.character(formula)[3]," "))[1]
  fit <- coxph(formula, data=data)
  a <- summary(fit)
  # output
  out <- cbind(
    a$coefficients[,c(1,4,5),drop=F],
    a$conf.int[,c(1,3,4),drop=F]
  )
  colnames(out) <- c("coef","z","P","HR","lower.95","upper.95")
  out <- out %>%
    as.data.frame() %>%
    rownames_to_column(var = "x") %>%
    mutate(
      HR2 = paste0(round(HR,2)," (",round(lower.95,2)," to ",round(upper.95,2),")")
    )
  if(!outputAll){
    tmp <- lapply(x,function(i,j) grep(i,j),j=out$x)
    tmp <- unlist(tmp)
    out <- out[tmp,]
  }
  out
}

HR_Forest <- function(dat){
  dat$` ` <- paste(rep(" ", 20), collapse = " ")
  colnames(dat)[colnames(dat)=="HR2"] <- "HR (95% CI)"
  forest(
    data = dat[,c(1,9,8,4),drop=F],
    est=dat$HR,
    lower = dat$lower.95,
    upper = dat$upper.95,
    ci_column = 2,
    ref_line = 1
  )
}

multiHR <- function(dat,x,time,censor,adj=NULL){
  res <- lapply(x,HR,dat=dat,time=time,censor=censor,adj=adj)
  res <- do.call(rbind,res)
  out <- as.data.frame(res)
  out$fdr <- p.adjust(out$P, "BH")
  out
}

LogT <- function(x){
  x[x==0] <- min(x[x>0])/2
  scale(log10(x))
}

combindfile2 <- function (filelist) 
{
  tmp <- filelist[[1]]
  tmp$name <- tmp[,1]
  for (i in 2:length(filelist)) {
    file2 <- filelist[[i]]
    file2$name <- filelist[[i]][,1]
    tmp <- merge(tmp, file2, by = "name", all = T)
  }
  rownames(tmp) <- tmp$name
  out <- tmp[, -which(colnames(tmp) == "name")]
  #out[is.na(out)] <- 0
  return(out)
}

Index <- function(dat, pos , neg){
  a1 <- dat[, pos]
  a2 <- dat[, neg]
  
  datx1 <- a1
  datx2 <- a2
  datx1 <- t(datx1)
  datx2 <- t(datx2)
  
  MH_species <- rownames(datx1)
  MN_species <- rownames(datx2)
  # Extracting Health-prevalent species present in metagenome
  # Extracting Health-scarce species present in metagenome
  MH_species_metagenome <- datx1
  MN_species_metagenome <- datx2
  
  # Diversity among Health-prevalent species
  # Diversity among Health-scarce species
  alpha <- function(x){sum((log(x[x>0]))*(x[x>0]))*(-1)}
  MH_shannon <- apply((MH_species_metagenome), 2, alpha) 
  MN_shannon <- apply((MN_species_metagenome), 2, alpha) 
  
  fun <- function(x){-log(x) * x}
  a <- apply(MH_species_metagenome, 2, fun) 
  b <- apply(MN_species_metagenome, 2, fun) 
  a <- apply(a,1,function(x) mean(x,na.rm=T))
  b <- apply(b,1,function(x) mean(x,na.rm=T))
  a <- data.frame(Mean=a)
  b <- data.frame(Mean=b)
  
  # Richness of Health-prevalent species
  # Richness of Health-scarce species
  R_MH <- apply(MH_species_metagenome, 2, function(i) (sum(i > 0))) 
  R_MN <- apply(MN_species_metagenome, 2, function(i) (sum(i > 0)))
  
  a <- data.frame(RMH=R_MH,RMN=R_MN)
  a <- a[order(a$RMN),]
  a <- a[order(a$RMH,decreasing = T),]
  
  # Median RMH from 1% of the top-ranked samples (see Methods)
  # Median RMN from 1% of the bottom-ranked samples (see Methods)
  MH_prime <- median(a[1:round(nrow(a)*0.01),"RMH"])
  MN_prime <- median(a[(round(nrow(a)*0.01):nrow(a)),"RMN"])
  
  
  # Collective abundance of Health-prevalent species
  # Collective abundance of Health-scarce species
  psi_MH <- ((R_MH/MH_prime)*MH_shannon) 
  psi_MN <- ((R_MN/MN_prime)*MN_shannon)
  
  # final OS-index
  dat.index <- data.frame(
    Index = log10((psi_MH+0.00001)/(psi_MN+0.00001))
  )
  rownames(dat.index) <- rownames(dat)
  dat.index
  
}

library(pec)
cindexFun_rfs <- function(dat,x,add){
  x <- union(x,add)
  x <- paste(x,collapse = "+")
  f <- paste0("Surv(RFS_5Years_with_censored, RFS_Status_5Years_with_censored)~",x)
  print(f)
  f <- as.formula(f)
  dynpred::cindex(f, data=dat)$cindex
}

cindexFun_os <- function(dat,x,add){
  x <- union(x,add)
  x <- paste(x,collapse = "+")
  f <- paste0("Surv(OS_5Years_with_censored, OS_Status_5Years_with_censored)~",x)
  print(f)
  f <- as.formula(f)
  dynpred::cindex(f, data=dat)$cindex
}

cindexFun.multi <- function(dat,x,add=NULL, class = "OS"){
  if(class == "OS"){
    a <- sapply(x,cindexFun_os,dat=dat,add=add)
  }else{
    a <- sapply(x,cindexFun_rfs,dat=dat,add=add)
  }
  a <- data.frame(Cindex=a)
}

plot_cindx <- function(a, b, col_name = c("raw","raw_index"), label = c("Raw","+Index")){
  
  a <- cbind(a,b)
  colnames(a) <- col_name
  a$x <- rownames(a)
  #a <- a[order(a$raw,decreasing = T),]
  a$x <- factor(a$x,levels = rev(a$x))
  a <- a %>% gather(type,Cindex,-x)
  a$type <- factor(a$type,levels = col_name, labels = label)
  ggplot(a,aes(x,Cindex,color=type))+
    geom_line(aes(group=x),color="black",arrow=arrow(length = unit(0.08, "inches")))+
    geom_point()+
    coord_flip()+
    theme_pubr()+
    theme(legend.position = "right")+
    xlab("")+ylab("C-index")+scale_color_manual(values = c("#235789", "#C1292E"))
  
}




