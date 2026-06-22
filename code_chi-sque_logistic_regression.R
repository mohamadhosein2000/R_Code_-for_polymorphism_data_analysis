### code for chi-squred test for all 3 SNP between Healthy vs Patients

genes <- list(SCN1A=SCN1A, ABCB1=ABCB1, GABRG2=GABRG2)

lapply(genes, function(x) {
  chisq.test(table(ifelse(x$group == 1, 1, 2), x$allele))
})

##### code for the output table

lapply(genes, function(x) {
  tbl <- table(ifelse(x$group == 1, 1, 2), x$allele)
  list(table = tbl, test = chisq.test(tbl))
})


### save full reslt (p-value + expected + table)
genes <- list(SCN1A, ABCB1, GABRG2)
names(genes) <- c("SCN1A", "ABCB1", "GABRG2")

chi_results <- lapply(genes, function(x) {
  tbl <- table(ifelse(x$group == 1, 1, 2), x$allele)
  test <- chisq.test(tbl)
  
  list(
    observed = tbl,
    expected = test$expected,
    statistic = unname(test$statistic),
    df = test$parameter,
    p.value = test$p.value
  )
})

##### summary of results
chi_summary <- data.frame(
  Gene = names(chi_results),
  Chi_square = sapply(chi_results, function(x) x$statistic),
  df = sapply(chi_results, function(x) x$df),
  p_value = sapply(chi_results, function(x) x$p.value)
)

chi_summary



##################################################################
##############################################################
########## allele analysis (healthy vs patients)
genes <- list(SCN1A, ABCB1, GABRG2)
names(genes) <- c("SCN1A", "ABCB1", "GABRG2")

allele_results <- lapply(genes, function(x) {
  
  grp  <- x$group[1:129]
  alle <- x$allele[1:129]
  
  # group 1
  T1 <- sum(alle[grp == 1] == 1) * 2 +
    sum(alle[grp == 1] == 2)
  
  C1 <- sum(alle[grp == 1] == 3) * 2 +
    sum(alle[grp == 1] == 2)
  
  # pooled group 2 + 3
  T23 <- sum(alle[grp != 1] == 1) * 2 +
    sum(alle[grp != 1] == 2)
  
  C23 <- sum(alle[grp != 1] == 3) * 2 +
    sum(alle[grp != 1] == 2)
  
  tbl <- matrix(c(T1, C1, T23, C23), nrow = 2, byrow = TRUE)
  colnames(tbl) <- c("T", "C")
  rownames(tbl) <- c("Group1", "Group2_3")
  
  test <- chisq.test(tbl)
  
  list(
    table = tbl,
    chi_square = unname(test$statistic),
    df = test$parameter,
    p.value = test$p.value
  )
})


#### summary 
allele_summary <- data.frame(
  Gene = names(allele_results),
  Chi_square = sapply(allele_results, function(x) x$chi_square),
  df = sapply(allele_results, function(x) x$df),
  p_value = sapply(allele_results, function(x) x$p.value)
)

allele_summary


x


#####################################
###############################################
###### genotype analysis (responsive vs resistant)
genes <- list(SCN1A, ABCB1, GABRG2)
names(genes) <- c("SCN1A", "ABCB1", "GABRG2")

genotype_results <- lapply(genes, function(x) {
  
  grp  <- x$group[1:129]
  alle <- x$allele[1:129]
  
  # select only group 2 and 3
  sel <- grp %in% c(2,3)
  grp23 <- grp[sel]
  alle23 <- alle[sel]
  
  # create 2 x 3 contingency table (group 2 vs 3, allele 1/2/3)
  tbl <- table(grp23, alle23)
  rownames(tbl) <- c("Group2", "Group3")
  colnames(tbl) <- c("Allele1", "Allele2", "Allele3")
  
  # chi-square test
  test <- chisq.test(tbl)
  
  list(
    table = tbl,
    chi_square = unname(test$statistic),
    df = test$parameter,
    p.value = test$p.value
  )
})

#### summary of genotype analysis (responsive vs resistant)
genotype_summary <- data.frame(
  Gene = names(genotype_results),
  Chi_square = sapply(genotype_results, function(x) x$chi_square),
  df = sapply(genotype_results, function(x) x$df),
  p_value = sapply(genotype_results, function(x) x$p.value)
)

genotype_summary




###################################
###################################
###################################
###################################
#### allele analysis (responsive vs resistant)
genes <- list(SCN1A, ABCB1, GABRG2)
names(genes) <- c("SCN1A", "ABCB1", "GABRG2")

allele_group23_results <- lapply(genes, function(x) {
  
  grp  <- x$group[1:129]
  alle <- x$allele[1:129]
  
  # select only group 2 and 3
  sel <- grp %in% c(2,3)
  grp23 <- grp[sel]
  alle23 <- alle[sel]
  
  # Compute allele counts per group
  # T = 2*TT + TC
  # C = 2*CC + TC
  T2 <- sum(alle23[grp23 == 2] == 1) * 2 + sum(alle23[grp23 == 2] == 2)
  C2 <- sum(alle23[grp23 == 2] == 3) * 2 + sum(alle23[grp23 == 2] == 2)
  
  T3 <- sum(alle23[grp23 == 3] == 1) * 2 + sum(alle23[grp23 == 3] == 2)
  C3 <- sum(alle23[grp23 == 3] == 3) * 2 + sum(alle23[grp23 == 3] == 2)
  
  # create 2x2 table
  tbl <- matrix(c(T2, C2, T3, C3), nrow = 2, byrow = TRUE)
  rownames(tbl) <- c("Group2", "Group3")
  colnames(tbl) <- c("T", "C")
  
  test <- chisq.test(tbl)
  
  list(
    table = tbl,
    chi_square = unname(test$statistic),
    df = test$parameter,
    p.value = test$p.value
  )
})

###################### summary of allele analysis (responsive vs resistant)
allele_group23_summary <- data.frame(
  Gene = names(allele_group23_results),
  Chi_square = sapply(allele_group23_results, function(x) x$chi_square),
  df = sapply(allele_group23_results, function(x) x$df),
  p_value = sapply(allele_group23_results, function(x) x$p.value)
)

allele_group23_summary


##################################33
###############################
################3 CODE FOR logistic regression (healthy vs patients) adjusted by sex and gender .....for SCN1A
## ===============================

# ===============================
# 1. SUBSET DATA
# ===============================
SCN1A_sub <- SCN1A[1:129, ]

# Outcome: 0 = group1, 1 = group2+3
SCN1A_sub$group_bin <- ifelse(SCN1A_sub$group == 1, 0, 1)

geno <- SCN1A_sub$allele
group <- SCN1A_sub$group_bin

# Age and sex (handle missing)
age <- if("age" %in% colnames(SCN1A_sub)) SCN1A_sub$age else rep(NA, nrow(SCN1A_sub))
sex <- if("gender" %in% colnames(SCN1A_sub)) SCN1A_sub$gender else rep(NA, nrow(SCN1A_sub))

# ===============================
# 2. ALLELIC MODEL (pooled alleles)
# ===============================

# Compute number of C and T alleles per person
C_count <- ifelse(geno == 1, 2, ifelse(geno == 2, 1, 0))
T_count <- ifelse(geno == 3, 2, ifelse(geno == 2, 1, 0))

# Make a data frame with group and allele counts
allelic_data <- data.frame(
  group = group,
  C = C_count,
  T = T_count,
  age = age,
  sex = sex
)

# For allelic logistic regression, expand each person into **alleles**
allele_long <- data.frame(
  allele = c(rep("C", sum(allelic_data$C)), rep("T", sum(allelic_data$T))),
  group = c(rep(0, sum(allelic_data$C[group==0])) , rep(1, sum(allelic_data$C[group==1])),
            rep(0, sum(allelic_data$T[group==0])) , rep(1, sum(allelic_data$T[group==1])))
)

# Logistic regression: allele ~ group
allelic_model <- glm(group ~ allele, data = allele_long, family = binomial())

# Extract OR, CI95%, p-value
confint_allelic <- confint(allelic_model)
colnames(confint_allelic) <- c("CI95_lower", "CI95_upper")

allelic_res <- cbind(
  OR = exp(coef(allelic_model)),
  CI95_lower = exp(confint_allelic[,1]),
  CI95_upper = exp(confint_allelic[,2]),
  p.value = summary(allelic_model)$coefficients[,4]
)

# Frequency table
allelic_freq <- rbind(
  Group1 = c(C = sum(allelic_data$C[group==0]), T = sum(allelic_data$T[group==0])),
  Group2_3 = c(C = sum(allelic_data$C[group==1]), T = sum(allelic_data$T[group==1]))
)

# ===============================
# 3. OTHER 6 GENETIC MODELS (genotype-based)
# ===============================
geno_het <- ifelse(geno == 2, 1, ifelse(geno == 1, 0, NA))
geno_hom <- ifelse(geno == 3, 1, ifelse(geno == 1, 0, NA))
geno_dom <- ifelse(geno == 1, 0, 1)
geno_rec <- ifelse(geno == 3, 1, 0)
geno_over <- ifelse(geno == 2, 1, 0)
geno_add <- ifelse(geno == 1, 0, ifelse(geno == 2, 1, 2))

pred_list <- list(
  Heterozygous = geno_het,
  Homozygous = geno_hom,
  Dominant = geno_dom,
  Recessive = geno_rec,
  Overdominant = geno_over,
  LogAdditive = geno_add
)

fit_model <- function(pred, group, age, sex){
  df <- data.frame(group_bin = group, pred = pred, age = age, sex = sex)
  m <- glm(group_bin ~ pred + age + sex, data = df, family = binomial())
  confint_m <- confint(m)
  colnames(confint_m) <- c("CI95_lower", "CI95_upper")
  res <- cbind(
    OR = exp(coef(m)),
    CI95_lower = exp(confint_m[,1]),
    CI95_upper = exp(confint_m[,2]),
    p.value = summary(m)$coefficients[,4]
  )
  return(res)
}

results_6 <- lapply(pred_list, function(p) fit_model(p, group, age, sex))

all_results <- c(list(Allelic = allelic_res), results_6)

# Frequency tables for verification
freq_list <- list(
  Allelic = allelic_freq,
  Heterozygous = table(group[!is.na(geno_het)], geno_het[!is.na(geno_het)]),
  Homozygous = table(group[!is.na(geno_hom)], geno_hom[!is.na(geno_hom)]),
  Dominant = table(group, geno_dom),
  Recessive = table(group, geno_rec),
  Overdominant = table(group, geno_over),
  LogAdditive = table(group, geno_add)
)

rownames(freq_list$Heterozygous) <- rownames(freq_list$Homozygous) <- 
  rownames(freq_list$Dominant) <- rownames(freq_list$Recessive) <- 
  rownames(freq_list$Overdominant) <- rownames(freq_list$LogAdditive) <- 
  c("Group1", "Group2+3")

colnames(freq_list$Heterozygous) <- c("CC", "CT")
colnames(freq_list$Homozygous) <- c("CC", "TT")
colnames(freq_list$Dominant) <- c("CC", "CT+TT")
colnames(freq_list$Recessive) <- c("CC+CT", "TT")
colnames(freq_list$Overdominant) <- c("CC+TT", "CT")
colnames(freq_list$LogAdditive) <- c("CC", "CT", "TT")

# ===============================
# 4. OUTPUT
# ===============================
all_results   # OR, CI95%, p-value for all 7 models
freq_list     # Frequency tables for verification


###########################################
#########################################
### one code for all SNP (SCN1A, ABCB1, GABRG2) logistic regression for 7 genetic models
# healthy(1) vs patients(2+3)
# ===============================
# FUNCTION TO RUN 7 GENETIC MODELS PER GENE
# ===============================
run_genetic_models <- function(gene_data){
  
  # Keep first 129 rows
  gene_sub <- gene_data[1:129, ]
  
  # Outcome: 0 = group1, 1 = group2+3
  gene_sub$group_bin <- ifelse(gene_sub$group == 1, 0, 1)
  group <- gene_sub$group_bin
  geno <- gene_sub$allele
  
  # Age and sex
  age <- if("age" %in% colnames(gene_sub)) gene_sub$age else rep(NA, nrow(gene_sub))
  sex <- if("gender" %in% colnames(gene_sub)) gene_sub$gender else rep(NA, nrow(gene_sub))
  
  # ===============================
  # 1. ALLELIC MODEL
  # ===============================
  allele_df <- data.frame(
    allele = unlist(lapply(geno, function(g) {
      if(g==1) return(c("C","C"))
      if(g==2) return(c("C","T"))
      if(g==3) return(c("T","T"))
    })),
    group = rep(group, each=2),
    age = rep(age, each=2),
    sex = rep(sex, each=2)
  )
  
  allele_df$allele_bin <- ifelse(allele_df$allele=="C", 0, 1)
  allelic_model <- glm(group ~ allele_bin + age + sex,
                       data=allele_df, family=binomial())
  
  confint_allelic <- confint(allelic_model)
  colnames(confint_allelic) <- c("CI95_lower", "CI95_upper")
  
  allelic_res <- cbind(
    OR = exp(coef(allelic_model)),
    CI95_lower = exp(confint_allelic[,1]),
    CI95_upper = exp(confint_allelic[,2]),
    p.value = summary(allelic_model)$coefficients[,4],
    AIC = AIC(allelic_model),
    BIC = BIC(allelic_model)
  )
  
  allelic_freq <- table(allele_df$group, allele_df$allele)
  rownames(allelic_freq) <- c("Group1","Group2+3")
  
  # ===============================
  # 2. OTHER 6 GENETIC MODELS
  # ===============================
  geno_het <- ifelse(geno == 2, 1, ifelse(geno == 1, 0, NA))
  geno_hom <- ifelse(geno == 3, 1, ifelse(geno == 1, 0, NA))
  geno_dom <- ifelse(geno == 1, 0, 1)
  geno_rec <- ifelse(geno == 3, 1, 0)
  geno_over <- ifelse(geno == 2, 1, 0)
  geno_add <- geno
  
  pred_list <- list(
    Heterozygous = geno_het,
    Homozygous = geno_hom,
    Dominant = geno_dom,
    Recessive = geno_rec,
    Overdominant = geno_over,
    LogAdditive = geno_add
  )
  
  fit_model <- function(pred, group, age, sex){
    df <- data.frame(group_bin = group, pred = pred, age = age, sex = sex)
    m <- glm(group_bin ~ pred + age + sex, data = df, family = binomial())
    
    confint_m <- confint(m)
    colnames(confint_m) <- c("CI95_lower", "CI95_upper")
    
    res <- cbind(
      OR = exp(coef(m)),
      CI95_lower = exp(confint_m[,1]),
      CI95_upper = exp(confint_m[,2]),
      p.value = summary(m)$coefficients[,4],
      AIC = AIC(m),
      BIC = BIC(m)
    )
    return(res)
  }
  
  results_6 <- lapply(pred_list, function(p) fit_model(p, group, age, sex))
  
  all_results <- c(list(Allelic = allelic_res), results_6)
  
  # ===============================
  # Frequency tables
  # ===============================
  freq_list <- list(
    Allelic = allelic_freq,
    Heterozygous = table(group[!is.na(geno_het)], geno_het[!is.na(geno_het)]),
    Homozygous = table(group[!is.na(geno_hom)], geno_hom[!is.na(geno_hom)]),
    Dominant = table(group, geno_dom),
    Recessive = table(group, geno_rec),
    Overdominant = table(group, geno_over),
    LogAdditive = table(group, geno_add)
  )
  
  rownames(freq_list$Heterozygous) <- rownames(freq_list$Homozygous) <- 
    rownames(freq_list$Dominant) <- rownames(freq_list$Recessive) <- 
    rownames(freq_list$Overdominant) <- rownames(freq_list$LogAdditive) <- 
    c("Group1", "Group2+3")
  
  colnames(freq_list$Heterozygous) <- c("CC", "CT")
  colnames(freq_list$Homozygous) <- c("CC", "TT")
  colnames(freq_list$Dominant) <- c("CC", "CT+TT")
  colnames(freq_list$Recessive) <- c("CC+CT", "TT")
  colnames(freq_list$Overdominant) <- c("CC+TT", "CT")
  colnames(freq_list$LogAdditive) <- c("CC", "CT", "TT")
  
  return(list(Results = all_results, Frequencies = freq_list))
}

# ===============================
# APPLY TO ALL 3 GENES
# ===============================
genes <- list(SCN1A = SCN1A, ABCB1 = ABCB1, GABRG2 = GABRG2)

gene_analysis <- lapply(genes, run_genetic_models)
names(gene_analysis) <- names(genes)

# ===============================
# AUTOMATIC DISPLAY IN CONSOLE
# ===============================
for(g in names(gene_analysis)){
  cat("\n==============================\n")
  cat("GENE:", g, "\n")
  cat("==============================\n\n")
  
  cat("LOGISTIC REGRESSION RESULTS (OR, CI, p, AIC, BIC):\n")
  print(gene_analysis[[g]]$Results)
  
  cat("\nFREQUENCY TABLES:\n")
  print(gene_analysis[[g]]$Frequencies)
}









##################################################################################################
##################################################################################################
##################################################################################################
##############################################################################################
############code for logistic regression for all SNP (adjusted by age and gender) ... responsive vs resistant 
# ===============================
# FUNCTION TO RUN 7 GENETIC MODELS BETWEEN GROUP2 AND GROUP3
# ===============================

run_genetic_models_group2vs3 <- function(gene_data){
  
  # Keep first 129 rows
  gene_sub <- gene_data[1:129, ]
  
  # Subset only group 2 and 3
  gene_sub <- subset(gene_sub, group %in% c(2,3))
  
  # Outcome: 0 = group2, 1 = group3
  gene_sub$group_bin <- ifelse(gene_sub$group == 2, 0, 1)
  group <- gene_sub$group_bin
  geno <- gene_sub$allele
  
  # Age and sex
  age <- if("age" %in% colnames(gene_sub)) gene_sub$age else rep(NA, nrow(gene_sub))
  sex <- if("gender" %in% colnames(gene_sub)) gene_sub$gender else rep(NA, nrow(gene_sub))
  
  # ===============================
  # 1. ALLELIC MODEL
  # ===============================
  allele_df <- data.frame(
    allele = unlist(lapply(geno, function(g) {
      if(g==1) return(c("C","C"))
      if(g==2) return(c("C","T"))
      if(g==3) return(c("T","T"))
    })),
    group = rep(group, each=2),
    age = rep(age, each=2),
    sex = rep(sex, each=2)
  )
  
  allele_df$allele_bin <- ifelse(allele_df$allele=="C", 0, 1)
  allelic_model <- glm(group ~ allele_bin + age + sex, data=allele_df, family=binomial())
  
  confint_allelic <- confint(allelic_model)
  colnames(confint_allelic) <- c("CI95_lower", "CI95_upper")
  
  allelic_res <- cbind(
    OR = exp(coef(allelic_model)),
    CI95_lower = exp(confint_allelic[,1]),
    CI95_upper = exp(confint_allelic[,2]),
    p.value = summary(allelic_model)$coefficients[,4],
    AIC = AIC(allelic_model),
    BIC = BIC(allelic_model)   # >>> AIC/BIC ADDED <<<
  )
  
  allelic_freq <- table(allele_df$group, allele_df$allele)
  rownames(allelic_freq) <- c("Group2","Group3")
  
  # ===============================
  # 2. OTHER 6 GENETIC MODELS
  # ===============================
  geno_het <- ifelse(geno == 2, 1, ifelse(geno == 1, 0, NA))
  geno_hom <- ifelse(geno == 3, 1, ifelse(geno == 1, 0, NA))
  geno_dom <- ifelse(geno == 1, 0, 1)
  geno_rec <- ifelse(geno == 3, 1, 0)
  geno_over <- ifelse(geno == 2, 1, 0)
  geno_add <- geno
  
  pred_list <- list(
    Heterozygous = geno_het,
    Homozygous = geno_hom,
    Dominant = geno_dom,
    Recessive = geno_rec,
    Overdominant = geno_over,
    LogAdditive = geno_add
  )
  
  fit_model <- function(pred, group, age, sex){
    df <- data.frame(group_bin = group, pred = pred, age = age, sex = sex)
    m <- glm(group_bin ~ pred + age + sex, data = df, family = binomial())
    confint_m <- confint(m)
    colnames(confint_m) <- c("CI95_lower", "CI95_upper")
    
    res <- cbind(
      OR = exp(coef(m)),
      CI95_lower = exp(confint_m[,1]),
      CI95_upper = exp(confint_m[,2]),
      p.value = summary(m)$coefficients[,4],
      AIC = AIC(m),
      BIC = BIC(m)   # >>> AIC/BIC ADDED <<<
    )
    return(res)
  }
  
  results_6 <- lapply(pred_list, function(p) fit_model(p, group, age, sex))
  all_results <- c(list(Allelic = allelic_res), results_6)
  
  # ===============================
  # Frequency tables
  # ===============================
  freq_list <- list(
    Allelic = allelic_freq,
    Heterozygous = table(group[!is.na(geno_het)], geno_het[!is.na(geno_het)]),
    Homozygous = table(group[!is.na(geno_hom)], geno_hom[!is.na(geno_hom)]),
    Dominant = table(group, geno_dom),
    Recessive = table(group, geno_rec),
    Overdominant = table(group, geno_over),
    LogAdditive = table(group, geno_add)
  )
  
  rownames(freq_list$Heterozygous) <- rownames(freq_list$Homozygous) <- 
    rownames(freq_list$Dominant) <- rownames(freq_list$Recessive) <- 
    rownames(freq_list$Overdominant) <- rownames(freq_list$LogAdditive) <- 
    c("Group2", "Group3")
  
  colnames(freq_list$Heterozygous) <- c("CC", "CT")
  colnames(freq_list$Homozygous) <- c("CC", "TT")
  colnames(freq_list$Dominant) <- c("CC", "CT+TT")
  colnames(freq_list$Recessive) <- c("CC+CT", "TT")
  colnames(freq_list$Overdominant) <- c("CC+TT", "CT")
  colnames(freq_list$LogAdditive) <- c("CC", "CT", "TT")
  
  return(list(Results = all_results, Frequencies = freq_list))
}

# ===============================
# APPLY TO ALL 3 GENES
# ===============================
genes <- list(SCN1A = SCN1A, ABCB1 = ABCB1, GABRG2 = GABRG2)

gene_analysis_2vs3 <- lapply(genes, run_genetic_models_group2vs3)
names(gene_analysis_2vs3) <- names(genes)

# ===============================
# AUTOMATIC DISPLAY IN CONSOLE
# ===============================
for(g in names(gene_analysis_2vs3)){
  cat("\n==============================\n")
  cat("GENE:", g, "\n")
  cat("==============================\n\n")
  
  cat("LOGISTIC REGRESSION RESULTS (OR, CI, p, AIC, BIC):\n")
  print(gene_analysis_2vs3[[g]]$Results)
  
  cat("\nFREQUENCY TABLES:\n")
  print(gene_analysis_2vs3[[g]]$Frequencies)
}



###############################################
###############################################
###########################################
############# logistic regression for gabrg2 (allelic model) for asm response adjusted by multiple factor

# ===============================
# ALLELIC LOGISTIC REGRESSION
# GABRG2 | Group 2 vs Group 3
# ===============================

# Use GABRG2 only
df <- GABRG2[1:129, ]

# Keep only group 2 and 3
df <- subset(df, group %in% c(2,3))

# Outcome: 0 = Group2, 1 = Group3
df$group_bin <- ifelse(df$group == 2, 0, 1)

# ===============================
# Expand genotypes into alleles
# ===============================
allele_df <- data.frame(
  allele = unlist(lapply(df$allele, function(g) {
    if (g == 1) return(c("C","C"))
    if (g == 2) return(c("C","T"))
    if (g == 3) return(c("T","T"))
  })),
  
  group = rep(df$group_bin, each = 2),
  age = rep(df$age, each = 2),
  gender = rep(df$gender, each = 2),
  intelectual_disability = rep(df$`intelectual disability`, each = 2),
  FC = rep(df$FC, each = 2),
  infantile_spasm = rep(df$`infantile spasm`, each = 2),
  Familial_Hx = rep(df$`Familial Hx`, each = 2),
  abnormal_EEG = rep(df$`abnormal EEG`, each = 2)
)

# Allele coding: C = 0, T = 1
allele_df$allele_bin <- ifelse(allele_df$allele == "C", 0, 1)

# ===============================
# Logistic regression
# ===============================
allelic_model <- glm(
  group ~ allele_bin +
    age +
    gender +
    intelectual_disability +
    FC +
    infantile_spasm +
    Familial_Hx +
    abnormal_EEG,
  data = allele_df,
  family = binomial()
)

# ===============================
# Results (OR, CI, p-value)
# ===============================
confint_m <- confint(allelic_model)

results <- data.frame(
  OR = exp(coef(allelic_model)),
  CI95_lower = exp(confint_m[,1]),
  CI95_upper = exp(confint_m[,2]),
  p.value = summary(allelic_model)$coefficients[,4]
)

print(results)

# ===============================
# Model fit: AIC & BIC
# ===============================
cat("\nAIC:", AIC(allelic_model), "\n")
cat("BIC:", BIC(allelic_model), "\n")

# ===============================
# Allele frequency check
# ===============================
allele_freq <- table(allele_df$group, allele_df$allele)
rownames(allele_freq) <- c("Group2","Group3")
print(allele_freq)



###############################################
###############################################
###########################################
############# logistic regression for gabrg2 (homozygouos model) for asm response adjusted by multiple factor


# ===============================
# HOMOZYGOUS LOGISTIC REGRESSION
# GABRG2 | Group 2 vs Group 3
# ===============================
# ===============================
# HOMOZYGOUS LOGISTIC REGRESSION
# GABRG2 | Group 2 vs Group 3
# ===============================

# Use GABRG2 only
df <- GABRG2[1:129, ]

# Keep only group 2 and 3
df <- subset(df, group %in% c(2,3))

# Keep only homozygous genotypes: CC (1) or TT (3)
df_homo <- subset(df, allele %in% c(1,3))

# Outcome: 0 = Group2, 1 = Group3
df_homo$group_bin <- ifelse(df_homo$group == 2, 0, 1)

# Predictor: homozygous genotype
# CC = 0, TT = 1
df_homo$geno_bin <- ifelse(df_homo$allele == 3, 1, 0)

# ===============================
# Logistic regression
# ===============================
# If needed, clean column names for safety
names(df_homo) <- make.names(names(df_homo))

homo_model <- glm(
  group_bin ~ geno_bin +
    age +
    gender +
    intelectual.disability +
    FC +
    infantile.spasm +
    Familial.Hx +
    abnormal.EEG,
  data = df_homo,
  family = binomial()
)

# ===============================
# Results (OR, CI, p-value)
# ===============================
confint_m <- confint(homo_model)

results <- data.frame(
  OR = exp(coef(homo_model)),
  CI95_lower = exp(confint_m[,1]),
  CI95_upper = exp(confint_m[,2]),
  p.value = summary(homo_model)$coefficients[,4]
)

print(results)

# ===============================
# Model fit: AIC & BIC
# ===============================
cat("\nAIC:", AIC(homo_model), "\n")
cat("BIC:", BIC(homo_model), "\n")

# ===============================
# Genotype frequency table
# ===============================
geno_freq <- table(df_homo$group_bin, df_homo$geno_bin)
rownames(geno_freq) <- c("Group2","Group3")
colnames(geno_freq) <- c("CC","TT")
print(geno_freq)


##############################33
#################################################3
############################################################

############# code for logistic regression (log-additive model) .... group2vs group3
####adjusted by factors

# ===============================
# LOG-ADDITIVE LOGISTIC REGRESSION
# GABRG2 | Group 2 vs Group 3
# ===============================

# Use GABRG2 only
df <- GABRG2[1:129, ]

# Keep only group 2 and 3
df <- subset(df, group %in% c(2,3))

# Outcome: 0 = Group2, 1 = Group3
df$group_bin <- ifelse(df$group == 2, 0, 1)

# Predictor: log-additive genotype coding
# CC=1, CT=2, TT=3
df$geno_add <- df$allele

# ===============================
# Clean column names for safety (spaces → dots)
# ===============================
names(df) <- make.names(names(df))

# ===============================
# Logistic regression
# ===============================
logadd_model <- glm(
  group_bin ~ geno_add +
    age +
    gender +
    intelectual.disability +
    FC +
    infantile.spasm +
    Familial.Hx +
    abnormal.EEG,
  data = df,
  family = binomial()
)

# ===============================
# Results (OR per allele, CI, p-value)
# ===============================
confint_m <- confint(logadd_model)

results <- data.frame(
  OR = exp(coef(logadd_model)),
  CI95_lower = exp(confint_m[,1]),
  CI95_upper = exp(confint_m[,2]),
  p.value = summary(logadd_model)$coefficients[,4]
)

print(results)

# ===============================
# Model fit: AIC & BIC
# ===============================
cat("\nAIC:", AIC(logadd_model), "\n")
cat("BIC:", BIC(logadd_model), "\n")

# ===============================
# Genotype frequency table
# ===============================
geno_freq <- table(df$group_bin, df$geno_add)
rownames(geno_freq) <- c("Group2","Group3")
colnames(geno_freq) <- c("CC=1","CT=2","TT=3")
print(geno_freq)


#######################################################
##########################################################
############################# code for logistic regression adjusted by factors
######## group2 vs 3..... recessive model
# ===============================
# RECESSIVE LOGISTIC REGRESSION
# GABRG2 | Group 2 vs Group 3
# ===============================

# Use GABRG2 only
df <- GABRG2[1:129, ]

# Keep only group 2 and 3
df <- subset(df, group %in% c(2,3))

# Outcome: 0 = Group2, 1 = Group3
df$group_bin <- ifelse(df$group == 2, 0, 1)

# Predictor: recessive genotype coding
# TT = 1, CC+CT = 0
df$geno_rec <- ifelse(df$allele == 3, 1, 0)

# ===============================
# Clean column names for safety
# ===============================
names(df) <- make.names(names(df))

# ===============================
# Logistic regression
# ===============================
recessive_model <- glm(
  group_bin ~ geno_rec +
    age +
    gender +
    intelectual.disability +
    FC +
    infantile.spasm +
    Familial.Hx +
    abnormal.EEG,
  data = df,
  family = binomial()
)

# ===============================
# Results (OR, CI, p-value)
# ===============================
confint_m <- confint(recessive_model)

results <- data.frame(
  OR = exp(coef(recessive_model)),
  CI95_lower = exp(confint_m[,1]),
  CI95_upper = exp(confint_m[,2]),
  p.value = summary(recessive_model)$coefficients[,4]
)

print(results)

# ===============================
# Model fit: AIC & BIC
# ===============================
cat("\nAIC:", AIC(recessive_model), "\n")
cat("BIC:", BIC(recessive_model), "\n")

# ===============================
# Genotype frequency table
# ===============================
geno_freq <- table(df$group_bin, df$geno_rec)
rownames(geno_freq) <- c("Group2","Group3")
colnames(geno_freq) <- c("CC+CT=0","TT=1")
print(geno_freq)










#############################
##############################
################################## HWE test for control participants (SCN1A)


################ 
## =========================================
## Hardy–Weinberg Equilibrium (HWE) Report
## Gene: SCN1A
## Controls only (group == 1)
## Genotype coding: 1 = homo1, 2 = hetero, 3 = homo2
## =========================================

# STEP 1: Subset controls and remove missing genotypes
controls <- subset(SCN1A, group == 1)
controls <- controls[!is.na(controls$allele), ]

# STEP 2: Count observed genotypes
geno_counts <- table(controls$allele)

# STEP 3: Create plain numeric vector (length 3) for HWE
obs <- c(
  as.numeric(geno_counts["1"]),
  as.numeric(geno_counts["2"]),
  as.numeric(geno_counts["3"])
)

# Replace NA with 0 if some genotype is missing
obs[is.na(obs)] <- 0

# Convert to numeric vector and remove names
obs <- as.numeric(obs)  # CRITICAL: HWExact needs a plain numeric vector

# STEP 4: Allele frequencies
N <- sum(obs)
p <- (2*obs[1] + obs[2]) / (2*N)
q <- 1 - p

# STEP 5: Expected counts under HWE
expected <- c(p^2, 2*p*q, q^2) * N

# STEP 6: Chi-square test (if expected counts ≥ 5)
chi_valid <- all(expected >= 5)
if(chi_valid){
  chi_test <- chisq.test(obs, p = expected / N, correct = FALSE)
}

# STEP 7: Exact test (always run)
if(!requireNamespace("HardyWeinberg", quietly = TRUE)){
  install.packages("HardyWeinberg")
}
library(HardyWeinberg)

# HWExact requires plain numeric vector
exact_test <- HWExact(obs)

# STEP 8: Reporting table
report_table <- data.frame(
  Genotype = c("Homozygote_1", "Heterozygote", "Homozygote_2"),
  Observed = obs,
  Expected = round(expected, 3)
)

cat("\n=========================================\n")
cat("Hardy–Weinberg Equilibrium Analysis\n")
cat("Gene: SCN1A (Controls only)\n")
cat("=========================================\n")

cat("\nObserved vs Expected genotype counts:\n")
print(report_table)

cat("\nChi-square HWE test:\n")
if(chi_valid){
  cat("Chi-square =", round(chi_test$statistic, 3), "\n")
  cat("df =", chi_test$parameter, "\n")
  cat("p-value =", chi_test$p.value, "\n")
} else {
  cat("Chi-square test NOT performed (expected counts < 5)\n")
}

cat("\nExact HWE test:\n")
cat("p-HWE (exact) =", exact_test$pval, "\n")

cat("\nFinal interpretation:\n")
if(exact_test$pval >= 0.05){
  cat("✔ Genotype distribution is consistent with Hardy–Weinberg equilibrium.\n")
} else {
  cat("✘ Genotype distribution deviates from Hardy–Weinberg equilibrium.\n")
}
cat("=========================================\n")



###################################################
##############################3
############ HWE for patients
## =========================================
## Hardy–Weinberg Equilibrium (HWE) Report
## Gene: SCN1A
## Pooled cases (group == 2 or 3)
## Genotype coding: 1 = homo1, 2 = hetero, 3 = homo2
## =========================================

# STEP 1: Subset pooled cases and remove missing genotypes
cases <- subset(SCN1A, group %in% c(2,3))
cases <- cases[!is.na(cases$allele), ]

# STEP 2: Count observed genotypes
geno_counts <- table(cases$allele)

# STEP 3: Create plain numeric vector (length 3) for HWE
obs <- c(
  as.numeric(geno_counts["1"]),
  as.numeric(geno_counts["2"]),
  as.numeric(geno_counts["3"])
)

# Replace NA with 0 if some genotype is missing
obs[is.na(obs)] <- 0

# Convert to numeric vector and remove names
obs <- as.numeric(obs)  # HWExact needs a plain numeric vector

# STEP 4: Allele frequencies
N <- sum(obs)
p <- (2*obs[1] + obs[2]) / (2*N)
q <- 1 - p

# STEP 5: Expected counts under HWE
expected <- c(p^2, 2*p*q, q^2) * N

# STEP 6: Chi-square test (if expected counts ≥ 5 and valid)
chi_valid <- all(!is.na(expected)) && all(expected >= 5)
if(chi_valid){
  chi_test <- chisq.test(obs, p = expected / N, correct = FALSE)
}

# STEP 7: Exact test (always run)
if(!requireNamespace("HardyWeinberg", quietly = TRUE)){
  install.packages("HardyWeinberg")
}
library(HardyWeinberg)

# HWExact requires plain numeric vector
exact_test <- HWExact(obs)

# STEP 8: Reporting table
report_table <- data.frame(
  Genotype = c("Homozygote_1", "Heterozygote", "Homozygote_2"),
  Observed = obs,
  Expected = round(expected, 3)
)

cat("\n=========================================\n")
cat("Hardy–Weinberg Equilibrium Analysis\n")
cat("Gene: SCN1A (Pooled cases: group 2+3)\n")
cat("=========================================\n")

cat("\nObserved vs Expected genotype counts:\n")
print(report_table)

cat("\nChi-square HWE test:\n")
if(chi_valid){
  cat("Chi-square =", round(chi_test$statistic, 3), "\n")
  cat("df =", chi_test$parameter, "\n")
  cat("p-value =", chi_test$p.value, "\n")
} else {
  cat("Chi-square test NOT performed (expected counts < 5 or missing data)\n")
}

cat("\nExact HWE test:\n")
cat("p-HWE (exact) =", exact_test$pval, "\n")

cat("\nFinal interpretation:\n")
if(exact_test$pval >= 0.05){
  cat("✔ Genotype distribution is consistent with Hardy–Weinberg equilibrium.\n")
} else {
  cat("✘ Genotype distribution deviates from Hardy–Weinberg equilibrium.\n")
}
cat("=========================================\n")


