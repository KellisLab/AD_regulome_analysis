library(speckle)
library(limma)
library(edgeR)
library(pheatmap)
library(gt)


#### scATAC, TSS>6
# major cell type (Figure 5C); non vs early vs late
info = readRDS("scATAC.meta.20230312.rds")
info.uniq = unique(info[,c("projid","msex","age_death","pmi","Pathology")])
info.full = info.uniq
counts <- table(info$Celltype4, info$projid)
trueprops <- rowSums(counts)/sum(rowSums(counts))

projid.level = info.uniq$projid
sample <- factor(info$projid,levels=projid.level)
group <- factor(info$Pathology, levels=c("nonAD","earlyAD","lateAD"))

grp <- factor(info.uniq$Pathology, levels=c("nonAD","earlyAD","lateAD"))
sex <- factor(info.uniq$msex,levels=c(0,1))
pmi = info.uniq$pmi
age_death = info.uniq$age_death

prop.logit <- getTransformedProps(clusters = info$Celltype4, sample=sample,
                                  transform = "logit")
design.anova <- model.matrix(~0+grp+sex+pmi+age_death)

propeller.anova(prop.logit,design = design.anova, coef=c(1,2,3), robust=TRUE, 
                trend = FALSE, sort=F)

# major cell type (Figure 5C); non = 1 vs early = 2 vs late =3
info = readRDS("scATAC.meta.20230312.rds")
info$Pathology.cvt = 1
info$Pathology.cvt[info$Pathology=="earlyAD"] = 2
info$Pathology.cvt[info$Pathology=="lateAD"] = 3
info.uniq = unique(info[,c("projid","msex","age_death","pmi","Pathology","Pathology.cvt")])
info.full = info.uniq
counts <- table(info$Celltype4, info$projid)
trueprops <- rowSums(counts)/sum(rowSums(counts))

projid.level = info.uniq$projid
sample <- factor(info$projid,levels=projid.level)
group <- info$Pathology.cvt

grp <- info.uniq$Pathology.cvt
sex <- factor(info.uniq$msex,levels=c(0,1))
pmi = info.uniq$pmi
age_death = info.uniq$age_death

prop.logit <- getTransformedProps(clusters = info$Celltype4, sample=sample,
                                  transform = "logit")
des.dose <- model.matrix(~grp+sex+pmi+age_death)

fit <- lmFit(prop.logit$TransformedProps,des.dose)
fit <- eBayes(fit, robust=TRUE)
topTable(fit,coef=2)

# sub cell type
info.major = readRDS("scATAC.meta.20230312.rds")

# number of cells to be considered for fractional analysis
ncut = 100
ncells = data.frame(table(info.major$pred_cell_type_high_resolution))
ncells = ncells[ncells$Freq>ncut,]
info.major = droplevels(info.major[info.major$pred_cell_type_high_resolution %in% ncells$Var1,])
major.cts = unique(info.major$Celltype4)

for(ct in major.cts){
info = droplevels(info.major[info.major$Celltype4==ct,])
info.uniq = unique(info[,c("projid","msex","age_death","pmi","Pathology")])
counts <- table(info$pred_cell_type_high_resolution, info$projid)
trueprops <- rowSums(counts)/sum(rowSums(counts))

projid.level = info.uniq$projid
sample <- factor(info$projid,levels=projid.level)
group <- factor(info$Pathology, levels=c("nonAD","earlyAD","lateAD"))

grp <- factor(info.uniq$Pathology, levels=c("nonAD","earlyAD","lateAD"))
sex <- factor(info.uniq$msex,levels=c(0,1))
pmi = info.uniq$pmi
age_death = info.uniq$age_death

prop.logit <- getTransformedProps(clusters = info$pred_cell_type_high_resolution, sample=sample,
                                  transform = "logit")
design.anova <- model.matrix(~0+grp+sex+pmi+age_death)

out = propeller.anova(prop.logit,design = design.anova, coef=c(1,2,3), robust=TRUE, 
                trend = FALSE, sort=TRUE)
print(out)
}

# sub-cell type (Figure 5D); non = 1 vs early = 2 vs late =3
# number of cells to be considered for fractional analysis
info.major = readRDS("scATAC.meta.20230312.rds")
info.major$Pathology.cvt = 1
info.major$Pathology.cvt[info.major$Pathology=="earlyAD"] = 2
info.major$Pathology.cvt[info.major$Pathology=="lateAD"] = 3

ncut = 100
ncells = data.frame(table(info.major$pred_cell_type_high_resolution))
ncells = ncells[ncells$Freq>ncut,]
info.major = droplevels(info.major[info.major$pred_cell_type_high_resolution %in% ncells$Var1,])
major.cts = unique(info.major$Celltype4)
major.cts = c("In","Ast","Ex","OPC","Microglia","PerEndo","Oligo")

for(ct in major.cts){
  info = droplevels(info.major[info.major$Celltype4==ct,])
  info.uniq = unique(info[,c("projid","msex","age_death","pmi","Pathology","Pathology.cvt")])
counts <- table(info$pred_cell_type_high_resolution, info$projid)
trueprops <- rowSums(counts)/sum(rowSums(counts))

projid.level = info.uniq$projid
sample <- factor(info$projid,levels=projid.level)
group <- info$Pathology.cvt

grp <- info.uniq$Pathology.cvt
sex <- factor(info.uniq$msex,levels=c(0,1))
pmi = info.uniq$pmi
age_death = info.uniq$age_death

prop.logit <- getTransformedProps(clusters = info$pred_cell_type_high_resolution, sample=sample,
                                  transform = "logit")
des.dose <- model.matrix(~grp+sex+pmi+age_death)

fit <- lmFit(prop.logit$TransformedProps,des.dose)
fit <- eBayes(fit, robust=TRUE)
out = topTable(fit,coef=2)
print(ct)
print(out)
}

###### snRNA cell composition
info = readRDS("scRNA.meta.0816.rds")
info.uniq = info.full
info = merge(info,info.uniq,by="projid")

counts <- table(info$major.celltype, info$projid)
trueprops <- rowSums(counts)/sum(rowSums(counts))

projid.level = info.uniq$projid
sample <- factor(info$projid,levels=projid.level)
group <- factor(info$Pathology, levels=c("nonAD","earlyAD","lateAD"))

grp <- factor(info.uniq$Pathology, levels=c("nonAD","earlyAD","lateAD"))
sex <- factor(info.uniq$msex,levels=c(0,1))
pmi = info.uniq$pmi
age_death = info.uniq$age_death

prop.logit <- getTransformedProps(clusters = info$major.celltype, sample=sample,
                                  transform = "logit")
design.anova <- model.matrix(~0+grp+sex+pmi+age_death)

propeller.anova(prop.logit,design = design.anova, coef=c(1,2,3), robust=TRUE, 
                trend = FALSE, sort=F)


# scRNA major cell type (Figure 5B); non = 1 vs early = 2 vs late =3
info = readRDS("scRNA.meta.0816.rds")
info.uniq = info.full 
info = merge(info,info.uniq,by="projid")
counts <- table(info$major.celltype , info$projid)
trueprops <- rowSums(counts)/sum(rowSums(counts))

projid.level = info.uniq$projid
sample <- factor(info$projid,levels=projid.level)
group <- info$Pathology.cvt

grp <- info.uniq$Pathology.cvt
sex <- factor(info.uniq$msex,levels=c(0,1))
pmi = info.uniq$pmi
age_death = info.uniq$age_death

prop.logit <- getTransformedProps(clusters = info$major.celltype, sample=sample,
                                  transform = "logit")
des.dose <- model.matrix(~grp+sex+pmi+age_death)

fit <- lmFit(prop.logit$TransformedProps,des.dose)
fit <- eBayes(fit, robust=TRUE)
topTable(fit,coef=2)

# scRNA sub-cell type (Figure S5D); non = 1 vs early = 2 vs late =3
# number of cells to be considered for fractional analysis
info.major = readRDS("scRNA.meta.0816.rds")
info.uniq = info.full
info.major = merge(info.major,info.uniq,by="projid")

ncut = 100
ncells = data.frame(table(info.major$cell_type_high_resolution))
ncells = ncells[ncells$Freq>ncut,]
info.major = droplevels(info.major[info.major$cell_type_high_resolution %in% ncells$Var1,])
major.cts = unique(info.major$major.celltype)
major.cts = c("Inh","Ast","Exc","Opc","Mic","Vas","Oli")

for(ct in major.cts){
  info = droplevels(info.major[info.major$major.celltype==ct,])
  info.uniq = unique(info[,c("projid","msex","age_death","pmi","Pathology","Pathology.cvt")])
  counts <- table(info$cell_type_high_resolution, info$projid)
  trueprops <- rowSums(counts)/sum(rowSums(counts))
  
  projid.level = info.uniq$projid
  sample <- factor(info$projid,levels=projid.level)
  group <- info$Pathology.cvt
  
  grp <- info.uniq$Pathology.cvt
  sex <- factor(info.uniq$msex,levels=c(0,1))
  pmi = info.uniq$pmi
  age_death = info.uniq$age_death
  
  prop.logit <- getTransformedProps(clusters = info$cell_type_high_resolution, sample=sample,
                                    transform = "logit")
  des.dose <- model.matrix(~grp+sex+pmi+age_death)
  
  fit <- lmFit(prop.logit$TransformedProps,des.dose)
  fit <- eBayes(fit, robust=TRUE)
  out = topTable(fit,coef=2)
  print(ct)
  print(out)
}

### TSS>1 eroded cells included
library("statmod")
info = readRDS("PFC.scATAC.erosion.20220925.rds")
info$age_death = as.numeric(as.character(info$age_death))
info$projid = info$ProjId
#info.uniq = unique(info[,c("projid","msex","age_death","Pathology")])
info.uniq = info.full
counts <- table(info$Celltype5, info$projid)
trueprops <- rowSums(counts)/sum(rowSums(counts))

projid.level = info.uniq$projid
sample <- factor(info$projid,levels=projid.level)
group <- factor(info$Pathology, levels=c("nonAD","earlyAD","lateAD"))

grp <- factor(info.uniq$Pathology, levels=c("nonAD","earlyAD","lateAD"))
sex <- factor(info.uniq$msex,levels=c(0,1))
pmi = info.uniq$pmi
age_death = info.uniq$age_death

prop.logit <- getTransformedProps(clusters = info$Celltype5, sample=sample,
                                  transform = "logit")
design.anova <- model.matrix(~0+grp+sex+age_death+pmi)

propeller.anova(prop.logit,design = design.anova, coef=c(1,2,3), robust=TRUE, 
                trend = FALSE, sort=TRUE)

## late vs. non+early
install.packages("statmod")
info = readRDS("PFC.scATAC.erosion.20220925.rds")
info$age_death = as.numeric(as.character(info$age_death))
info$projid = info$ProjId
info$Pathology2 = "non_early_AD"
info$Pathology2[info$Pathology=="lateAD"] = "lateAD"
info$Pathology = info$Pathology2
info.uniq = unique(info[,c("projid","msex","age_death","Pathology")])
counts <- table(info$Celltype5, info$projid)
trueprops <- rowSums(counts)/sum(rowSums(counts))

projid.level = info.uniq$projid
sample <- factor(info$projid,levels=projid.level)
group <- factor(info$Pathology, levels=c("non_early_AD","lateAD"))

grp <- factor(info.uniq$Pathology, levels=c("non_early_AD","lateAD"))
sex <- factor(info.uniq$msex,levels=c(0,1))
age_death = info.uniq$age_death

prop.logit <- getTransformedProps(clusters = info$Celltype5, sample=sample,
                                  transform = "logit")
design.anova <- model.matrix(~0+grp+sex+age_death)

propeller.anova(prop.logit,design = design.anova, coef=c(1,2), robust=TRUE, 
                trend = FALSE, sort=TRUE)








