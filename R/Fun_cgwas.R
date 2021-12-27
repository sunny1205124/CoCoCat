#### to make sure that var(x)>0
f0 <- function(X0){
  sum(X0)>=1
}

#### The conditional regression model
regression <- function(m_time, data_m, fmly, p_value, name_y){
  del_causal <- NULL
  name_causal <- colnames(data_m)[1:(ncol(data_m)-1)]
  idx_m <- combn(length(name_causal),m_time)
  if(length(name_causal)-m_time>0){
    for (cc in 1:(length(name_causal)-m_time)) {
      idx_m[,idx_m[1,]==cc] <- idx_m[,idx_m[1,]==cc][,order(idx_m[m_time,idx_m[1,]==cc],decreasing = TRUE)]
    }
  }else{
    idx_m <- idx_m
  }
  for (aa in 1:ncol(idx_m)) {
    mul_snp <- name_causal[idx_m[,aa]]
    if(sum(mul_snp %in% del_causal)>0) next
    forml <- paste0(name_y,'~',paste(mul_snp,collapse = '+'))
    fit0 <- glm(as.formula(forml),data_m,family = fmly)
    res <- summary(fit0)$coef[-1,]
    if(sum(is.na(fit0$coefficients))>=(m_time-1)){
      del_causal <- del_causal
    }else if(sum(res[,4]>p_value,na.rm = T)==nrow(res)){
      del_causal <- del_causal
    }else{
      del_causal <- c(del_causal,names(which(res[,4]>p_value)))
    }
  }
  return(del_causal)
}

#### multi-stage GWAS fine mapping process
regression_cgwas <- function(
    data_snp, # the imported data of genomic and outcome data
    name_snps, # the column name of the imported data
    m_times, # the number of elements in multiple regression models
    Pvalue, # the threshold of significance of SNPs
    family0, # 'gaussian' for quantitative traits, 'binary' for qualitative diseases
    Ny # the name of the outcome y
    ){
  if (length(name_snps)>2){
    del_snp <- name_snps
    while (length(name_snps)>m_times & length(del_snp)!=0){
      cat(m_times,'\n')
      del_snp <- regression(m_time=m_times, data_m=data_snp, fmly=family0, p_value=Pvalue, name_y = Ny)
      name_snps <- name_snps[!name_snps %in% del_snp]
      data_snp <- data_snp[,name_snps]
      m_times=m_times+1
    }
    name_snp <- name_snps[1:(length(name_snps)-1)]
    m_time <- m_times-1
  }else{
    name_snp <- name_snps[1:(length(name_snps)-1)]
    m_time <- m_times-1
  }
  return(name_snp)
}


#### run the CoCoCat function
CGWAS <- function(data0, Pvalue, family0){
  name_snp0 <- colnames(data0)
  NameY <- name_snp0[length(name_snp0)]
  n_snp <- ncol(data0)-1
  ##### univariable regression model
  coef1 <- rep(NA,n_snp)
  p1 <- rep(NA,n_snp)

  for (m in 1:n_snp) {
    uni_snp <- data0[,m]
    outcomeY <- data0[,ncol(data0)]
    fit1 <- glm(outcomeY~uni_snp, family=family0)
    if(!is.na((fit1)$coef[2])){
      coef1[m] <- summary(fit1)$coef[2,1]
      p1[m] <- summary(fit1)$coef[2,4]
    }else{
      coef1[m] <- NA
      p1[m] <- NA
    }
  }

  data_s0 <- data0[,1:n_snp][,!is.na(p1)]
  data_s1 <- data_s0[,p1[!is.na(p1)]<Pvalue]
  data_s2 <- data_s1[,order(p1[!is.na(p1)][p1[!is.na(p1)]<Pvalue])]
  data_st <- data.frame(data_s2,outcomeY)
  name_snp1 <- colnames(data_s0)[p1[!is.na(p1)]<Pvalue]

  # continue to apply binary regression models
  m_time0=2
  name_csnp <- regression_cgwas(data_snp=data_st, name_snps=name_snp1,
                                m_times=m_time0, Pvalue=Pvalue, family0=family0, Ny='outcomeY')
  return(name_csnp)
}
