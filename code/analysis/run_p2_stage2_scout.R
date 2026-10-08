options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(pheatmap); library(scales)})

script_arg <- grep("^--file=",commandArgs(trailingOnly=FALSE),value=TRUE); script_path <- if(length(script_arg)) sub("^--file=","",script_arg[1]) else "."
root <- normalizePath(Sys.getenv("ORGANOID_PROJECT_ROOT",unset=file.path(dirname(script_path),"../..")),mustWork=FALSE)
lfc_file <- file.path(root, "00_source/nature_supplementary_tables/supplementary_table_6_revision.csv")
bin_file <- file.path(root, "00_source/nature_supplementary_tables/supplementary_table_5_revision.csv")
universe_file <- file.path(root, "03_gene_sets/METABOLIC_GENE_UNIVERSE.tsv")
meta_file <- file.path(root, "01_metadata/ORGANOID_METADATA.tsv")
control_file <- file.path(root, "00_source/figshare_28339340/control-genes_organoids.csv")
common_file <- file.path(root, "00_source/figshare_28339340/CRISPRInferredCommonEssentials.csv")
dir.create(file.path(root, "06_tables"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(root, "05_figures"), recursive=TRUE, showWarnings=FALSE)

write_tsv <- function(x, name) fwrite(x, file.path(root, name), sep="\t", na="NA", quote=FALSE)
safe_mad <- function(x) mad(x, constant=1, na.rm=TRUE)
iqr <- function(x) IQR(x, na.rm=TRUE, type=7)
collapse_unique <- function(x) paste(sort(unique(x[!is.na(x) & x != ""])), collapse=";")

# P2-1.5: membership tiering. GO is retained as broad ontology Tier 3. Curated
# Reactome/KEGG/Hallmark memberships are Tier 1 unless their set label explicitly
# describes regulation/response/signalling/disease-level context, which is Tier 2.
u <- fread(universe_file)
broad_pattern <- "REGULAT|RESPONSE|SIGNAL|CANCER|DISEASE|TRANSCRIPTION|EXPRESSION|HOMEOSTASIS|DEVELOP|DIFFERENTIAT|ORGANISM|CELLULAR_RESPONSE|HISTONE|METHYLAT.*PROTEIN|PROTEIN.*METHYLAT|CELL_CYCLE|APOPTOSIS|IMMUNE"
u[, tier := fifelse(source_database == "GO", "Tier 3",
                    fifelse(source_database == "MSigDB" | grepl(broad_pattern, source_gene_set, ignore.case=TRUE), "Tier 2", "Tier 1"))]
u[, tier_rule := fifelse(tier == "Tier 1", "curated Reactome/KEGG/Hallmark canonical metabolic or transport membership",
  fifelse(tier == "Tier 2", "Hallmark signature or curated extended/regulatory/response/disease-context membership",
                         "GO Biological Process ontology membership; retained but excluded from primary scores"))]
setcolorder(u, c("gene_symbol","ensembl_id","pathway_major","pathway_minor","source_database","source_gene_set","evidence_level","tier","tier_rule"))
write_tsv(u, "03_gene_sets/METABOLIC_GENE_UNIVERSE_TIERED.tsv")

qc <- u[, .(membership_rows=.N, unique_genes=uniqueN(gene_symbol)), by=.(pathway_major,tier)]
allpaths <- CJ(pathway_major=sort(unique(u$pathway_major)), tier=c("Tier 1","Tier 2","Tier 3"))
qc <- merge(allpaths, qc, by=c("pathway_major","tier"), all.x=TRUE)
qc[is.na(membership_rows), `:=`(membership_rows=0L, unique_genes=0L)]
qc[, tier1_size := unique_genes[tier == "Tier 1"], by=pathway_major]
qc[, tier1_size_flag := fifelse(tier1_size < 5, "LT5", fifelse(tier1_size > 150, "GT150", "OK"))]
write_tsv(qc, "06_tables/PATHWAY_TIER_QC.tsv")

met_members <- unique(u[tier %in% c("Tier 1","Tier 2"), .(gene=gene_symbol, pathway_major, tier)])
met_genes <- sort(unique(met_members$gene))
meta <- fread(meta_file)
meta <- meta[CRISPR_available == "Yes", .(sample_ID, cancer_type=primary_tumour_type)]

# Read the full matrices once; six columns and 2.95M rows each.
lfc <- fread(lfc_file, sep=" ", header=TRUE, showProgress=FALSE)
bin <- fread(bin_file, sep=" ", header=TRUE, showProgress=FALSE)

# Control-informed strong-dependency threshold: midpoint between the pooled median
# LFC of supplied curated BAGEL essential and non-essential controls.
controls <- fread(control_file)
control_labels <- controls[essentiality %in% c("curated_BAGEL_essential","curated_BAGEL_non_essential"), .(gene, essentiality)]
control_lfc <- merge(lfc[, .(gene,LFC)], control_labels, by="gene")
control_medians <- control_lfc[, .(median_LFC=median(LFC)), by=essentiality]
ess_med <- control_medians[essentiality == "curated_BAGEL_essential", median_LFC]
non_med <- control_medians[essentiality == "curated_BAGEL_non_essential", median_LFC]
strong_cutoff <- mean(c(ess_med, non_med))

# Restrict to metabolic Tier 1/2 genes. Library adjustment is gene-wise:
# subtract the gene×library mean and add the gene global mean. Then average the
# adjusted LFC across libraries for each organoid–gene pair. This preserves model
# as the analysis unit and prevents 16 dual-library models from being double-counted.
lfcm <- lfc[gene %in% met_genes]
lfcm[, gene_global_mean := mean(LFC), by=gene]
lfcm[, library_gene_mean := mean(LFC), by=.(gene,library)]
lfcm[, LFC_library_adjusted := LFC - library_gene_mean + gene_global_mean]
cons <- lfcm[, .(
  consensus_LFC=mean(LFC), adjusted_consensus_LFC=mean(LFC_library_adjusted),
  n_libraries=uniqueN(library), libraries=collapse_unique(library),
  minLib_LFC=if(any(library=="minLib")) mean(LFC[library=="minLib"]) else NA_real_,
  v1_1_LFC=if(any(library=="v1_1")) mean(LFC[library=="v1_1"]) else NA_real_
), by=.(sample_ID,gene)]
binm <- bin[gene %in% met_genes, .(binary_depleted_probability=mean(is_depleted), binary_calls=collapse_unique(as.character(is_depleted))), by=.(sample_ID,gene)]
cons <- merge(cons, binm, by=c("sample_ID","gene"), all.x=TRUE)
cons <- merge(cons, meta, by="sample_ID", all.x=TRUE)
stopifnot(uniqueN(cons$sample_ID)==162, max(cons[, uniqueN(sample_ID), by=gene]$V1)<=162)
write_tsv(cons, "02_crispr/METABOLIC_ORGANOID_GENE_CONSENSUS.tsv")

# Gene-level landscape.
gene_info <- u[tier %in% c("Tier 1","Tier 2"), .(
  pathways=collapse_unique(pathway_major), tiers=collapse_unique(tier),
  ensembl_id=collapse_unique(ensembl_id)
), by=.(gene=gene_symbol)]
gene_sum <- cons[, .(
  N_organoids=.N, median_LFC=median(adjusted_consensus_LFC), mean_LFC=mean(adjusted_consensus_LFC),
  IQR=iqr(adjusted_consensus_LFC), MAD=safe_mad(adjusted_consensus_LFC),
  variance=var(adjusted_consensus_LFC), min_LFC=min(adjusted_consensus_LFC), max_LFC=max(adjusted_consensus_LFC),
  dependency_range=max(adjusted_consensus_LFC)-min(adjusted_consensus_LFC),
  binary_depleted_fraction=mean(binary_depleted_probability),
  strongly_dependent_fraction=mean(adjusted_consensus_LFC <= strong_cutoff)
), by=gene]
for(ct in c("Colorectal","Oesophageal","Ovarian","Pancreatic","Gastric")) {
  z <- cons[cancer_type==ct, .(
    n=.N, median=median(adjusted_consensus_LFC), mean=mean(adjusted_consensus_LFC),
    depleted_fraction=mean(binary_depleted_probability), strong_fraction=mean(adjusted_consensus_LFC <= strong_cutoff)
  ), by=gene]
  setnames(z, c("n","median","mean","depleted_fraction","strong_fraction"), paste0(c("N_","median_LFC_","mean_LFC_","depleted_fraction_","strong_fraction_"),ct))
  gene_sum <- merge(gene_sum,z,by="gene",all.x=TRUE)
}
gene_sum <- merge(gene_info,gene_sum,by="gene",all.y=TRUE)
gene_sum[, strong_dependency_threshold := strong_cutoff]
write_tsv(gene_sum, "06_tables/METABOLIC_GENE_DEPENDENCY_SCOUT.tsv")

# Gene-wise standardization across organoids; positive scores mean stronger relative
# dependency. Rank sensitivity uses percentile rank of -LFC within each gene.
cons[, gene_sd := sd(adjusted_consensus_LFC), by=gene]
cons[, dependency_z := fifelse(gene_sd>0, -(adjusted_consensus_LFC-mean(adjusted_consensus_LFC))/gene_sd, 0), by=gene]
cons[, dependency_rank := frank(-adjusted_consensus_LFC, ties.method="average")/.N, by=gene]

score_membership <- function(tiers, zcol, outcol) {
  mem <- unique(met_members[tier %in% tiers, .(gene,pathway_major)])
  x <- merge(cons[, .(sample_ID,gene,cancer_type,value=get(zcol))], mem, by="gene", allow.cartesian=TRUE)
  x[, .(score=median(value), n_genes=uniqueN(gene)), by=.(sample_ID,cancer_type,pathway_major)][]
}
primary <- score_membership("Tier 1", "dependency_z", "primary")
rank_sens <- score_membership("Tier 1", "dependency_rank", "rank")
t12_sens <- score_membership(c("Tier 1","Tier 2"), "dependency_z", "t12")
setnames(primary,c("score","n_genes"),c("tier1_z_median_score","n_tier1_genes"))
setnames(rank_sens,c("score","n_genes"),c("tier1_rank_median_score","n_tier1_rank_genes"))
setnames(t12_sens,c("score","n_genes"),c("tier1plus2_z_median_score","n_tier1plus2_genes"))
pathmat <- Reduce(function(x,y) merge(x,y,by=c("sample_ID","cancer_type","pathway_major"),all=TRUE), list(primary,rank_sens,t12_sens))
pathmat <- pathmat[!is.na(tier1_z_median_score)]
write_tsv(pathmat, "06_tables/PATHWAY_DEPENDENCY_MATRIX.tsv")
pathsum <- pathmat[, .(
  N_organoids=.N, n_tier1_genes=as.integer(max(n_tier1_genes,na.rm=TRUE)),
  median_primary_score=median(tier1_z_median_score,na.rm=TRUE), mean_primary_score=mean(tier1_z_median_score,na.rm=TRUE),
  IQR_primary=iqr(tier1_z_median_score), MAD_primary=safe_mad(tier1_z_median_score),
  median_rank_sensitivity=median(tier1_rank_median_score,na.rm=TRUE),
  median_tier1plus2_sensitivity=median(tier1plus2_z_median_score,na.rm=TRUE),
  spearman_primary_vs_rank=cor(tier1_z_median_score,tier1_rank_median_score,use="pairwise.complete.obs",method="spearman"),
  spearman_primary_vs_tier1plus2=cor(tier1_z_median_score,tier1plus2_z_median_score,use="pairwise.complete.obs",method="spearman")
), by=pathway_major]
write_tsv(pathsum, "06_tables/PATHWAY_DEPENDENCY_SUMMARY.tsv")

# CRC versus ESCA: non-parametric independent-model test, library-adjusted primary
# effect and unadjusted/library-specific sensitivity estimates.
mw_stats <- function(x,y) {
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  if(length(x)<2 || length(y)<2) return(list(p=NA_real_,u=NA_real_,rbc=NA_real_,hl=NA_real_,lo=NA_real_,hi=NA_real_))
  wt <- suppressWarnings(wilcox.test(x,y,alternative="two.sided",exact=FALSE,conf.int=TRUE,conf.level=.95))
  U <- unname(wt$statistic)
  list(p=wt$p.value,u=U,rbc=2*U/(length(x)*length(y))-1,
       hl=if(!is.null(wt$estimate)) unname(wt$estimate) else NA_real_,
       lo=if(!is.null(wt$conf.int)) wt$conf.int[1] else NA_real_, hi=if(!is.null(wt$conf.int)) wt$conf.int[2] else NA_real_)
}
diff_one <- function(d, value_col, id_col) {
  ids <- unique(d[[id_col]])
  rbindlist(lapply(ids,function(id){
    a <- d[get(id_col)==id & cancer_type=="Colorectal", get(value_col)]
    b <- d[get(id_col)==id & cancer_type=="Oesophageal", get(value_col)]
    s <- mw_stats(a,b)
    data.table(id=id,n_CRC=length(a),n_ESCA=length(b),median_CRC=median(a),median_ESCA=median(b),
      median_difference_CRC_minus_ESCA=median(a)-median(b),hodges_lehmann_shift=s$hl,
      CI95_low=s$lo,CI95_high=s$hi,U=s$u,rank_biserial_CRC_vs_ESCA=s$rbc,raw_P=s$p)
  }))
}
gd <- diff_one(cons, "adjusted_consensus_LFC", "gene")
setnames(gd,"id","gene")
gu <- diff_one(cons, "consensus_LFC", "gene")[, .(gene=id,unadjusted_median_difference=median_difference_CRC_minus_ESCA,unadjusted_raw_P=raw_P)]
gd <- merge(gd,gu,by="gene",all.x=TRUE)
gmin <- diff_one(cons[!is.na(minLib_LFC)], "minLib_LFC", "gene")[, .(gene=id,minLib_n_CRC=n_CRC,minLib_n_ESCA=n_ESCA,minLib_median_difference=median_difference_CRC_minus_ESCA,minLib_raw_P=raw_P)]
gv1 <- diff_one(cons[!is.na(v1_1_LFC)], "v1_1_LFC", "gene")[, .(gene=id,v1_1_n_CRC=n_CRC,v1_1_n_ESCA=n_ESCA,v1_1_median_difference=median_difference_CRC_minus_ESCA,v1_1_raw_P=raw_P)]
gd <- Reduce(function(x,y) merge(x,y,by="gene",all.x=TRUE),list(gd,gmin,gv1))
gd[, BH_FDR := p.adjust(raw_P,method="BH")]
gd <- merge(gene_info,gd,by="gene",all.y=TRUE)
write_tsv(gd[order(BH_FDR,-abs(rank_biserial_CRC_vs_ESCA))], "06_tables/CRC_ESCA_GENE_DIFFERENTIAL_DEPENDENCY.tsv")

pdiff <- diff_one(pathmat, "tier1_z_median_score", "pathway_major")
setnames(pdiff,"id","pathway_major")
pdiff[, BH_FDR := p.adjust(raw_P,method="BH")]
write_tsv(pdiff[order(BH_FDR,-abs(rank_biserial_CRC_vs_ESCA))], "06_tables/CRC_ESCA_PATHWAY_DIFFERENTIAL_DEPENDENCY.tsv")

# Selectivity landscape and empirical classification.
shape <- cons[, {
  x <- adjusted_consensus_LFC; n <- length(x); m <- mean(x); s <- sd(x)
  skew <- if(s>0) mean(((x-m)/s)^3) else 0
  kurt <- if(s>0) mean(((x-m)/s)^4)-3 else 0
  bc <- (skew^2+1)/(kurt+3 + 3*(n-1)^2/((n-2)*(n-3)))
  .(skewness=skew,excess_kurtosis=kurt,bimodality_coefficient=bc,
    extreme_dependency_fraction=mean(x <= quantile(x,.10,type=7)))
},by=gene]
sel <- merge(gene_sum,shape,by="gene")
common <- sub(" \\(.*$", "", fread(common_file, sep=",")[[1]])
sel[, common_essential_flag := gene %in% common]
mad75 <- quantile(sel$MAD,.75,na.rm=TRUE); iqr75 <- quantile(sel$IQR,.75,na.rm=TRUE)
sel[, dependency_class := fifelse(binary_depleted_fraction>=.80 & median_LFC<=strong_cutoff,"A_common_metabolic_essential",
  fifelse(binary_depleted_fraction>=.10 & binary_depleted_fraction<.80 & (MAD>=mad75 | IQR>=iqr75),"B_selective_metabolic_dependency","C_weak_or_non_dependent"))]
sel[, bimodality_flag_descriptive := bimodality_coefficient > 5/9]
sel[, `:=`(MAD_75pct_threshold=mad75,IQR_75pct_threshold=iqr75,strong_dependency_threshold=strong_cutoff)]
write_tsv(sel[order(dependency_class,-MAD,median_LFC)], "06_tables/SELECTIVE_METABOLIC_DEPENDENCIES.tsv")

# Integrated data-driven shortlists. Union of top global strength, top selectivity,
# and FDR-significant cancer differences with at least upper-quartile absolute effect.
cand <- merge(sel,gd[,.(gene,CRC_ESCA_effect=rank_biserial_CRC_vs_ESCA,CRC_ESCA_median_difference=median_difference_CRC_minus_ESCA,FDR=BH_FDR)],by="gene",all.x=TRUE)
strength_cut <- quantile(-cand$median_LFC,.95,na.rm=TRUE)
select_cut <- quantile(cand$MAD,.90,na.rm=TRUE)
effect_cut <- quantile(abs(cand$CRC_ESCA_effect),.75,na.rm=TRUE)
cand[, reason_shortlisted := paste(
  fifelse(-median_LFC>=strength_cut,"top-5% global dependency strength",""),
  fifelse(MAD>=select_cut,"top-10% selectivity (MAD)",""),
  fifelse(!is.na(FDR)&FDR<.05&abs(CRC_ESCA_effect)>=effect_cut,"FDR-significant CRC-ESCA difference with upper-quartile effect",""), sep="; ")]
cand[, reason_shortlisted := gsub("(^; )|(; ; )|(; $)","",reason_shortlisted)]
cand <- cand[nchar(reason_shortlisted)>0]
cand[, cancer_specificity := fifelse(CRC_ESCA_median_difference<0,"CRC-leaning",fifelse(CRC_ESCA_median_difference>0,"ESCA-leaning","none/undetermined"))]
short <- cand[,.(gene,pathway=pathways,global_dependency_strength=-median_LFC,selectivity_MAD=MAD,
  depleted_fraction=binary_depleted_fraction,CRC_vs_ESCA_effect=CRC_ESCA_effect,FDR,cancer_specificity,
  Tier=tiers,common_essential_flag,reason_shortlisted)]
write_tsv(short[order(FDR,-global_dependency_strength,-selectivity_MAD)], "06_tables/P2_SCOUT_SHORTLIST.tsv")

pshort <- merge(pathsum,pdiff,by="pathway_major")
pshort[, abs_effect:=abs(rank_biserial_CRC_vs_ESCA)]
pshort <- pshort[BH_FDR<.10 | abs_effect>=quantile(abs_effect,.75,na.rm=TRUE) | MAD_primary>=quantile(MAD_primary,.75,na.rm=TRUE)]
pshort[, reason_shortlisted := paste0("data-driven pathway scout: FDR=",signif(BH_FDR,3),", |rank-biserial|=",signif(abs_effect,3),", MAD=",signif(MAD_primary,3))]
write_tsv(pshort[order(BH_FDR,-abs_effect)], "06_tables/P2_SCOUT_PATHWAY_SHORTLIST.tsv")

# Figures: colorblind-safe, PDF vector + 300-dpi PNG.
theme_set(theme_classic(base_size=8))
pal <- c(Colorectal="#0072B2",Oesophageal="#D55E00",Ovarian="#CC79A7",Pancreatic="#009E73",Gastric="#E69F00")
save_plot <- function(p,name,w=7.2,h=4.5){
  ggsave(file.path(root,"05_figures",paste0(name,".pdf")),p,width=w,height=h,units="in",device=cairo_pdf)
  ggsave(file.path(root,"05_figures",paste0(name,".png")),p,width=w,height=h,units="in",dpi=300)
}

# Heatmap sorted by cancer type; pathways clustered, samples not clustered.
hm <- dcast(pathmat, pathway_major~sample_ID, value.var="tier1_z_median_score")
rn <- hm$pathway_major; hm$pathway_major <- NULL; mat <- as.matrix(hm); rownames(mat)<-rn
ordmeta <- meta[match(colnames(mat),sample_ID)]; ord <- order(factor(ordmeta$cancer_type,levels=names(pal)),ordmeta$sample_ID)
mat <- mat[,ord,drop=FALSE]; ann <- data.frame(Cancer=ordmeta$cancer_type[ord]); rownames(ann)<-colnames(mat)
anncols <- list(Cancer=pal)
pdf(file.path(root,"05_figures/P2_F01_pathway_dependency_heatmap.pdf"),width=7.2,height=5.2,useDingbats=FALSE)
pheatmap(mat,cluster_cols=FALSE,cluster_rows=TRUE,show_colnames=FALSE,fontsize_row=6,
         annotation_col=ann,annotation_colors=anncols,color=colorRampPalette(c("#2166AC","white","#B2182B"))(101),
         main="Tier-1 pathway dependency (gene-wise z-score median)",border_color=NA)
dev.off()
png(file.path(root,"05_figures/P2_F01_pathway_dependency_heatmap.png"),width=2160,height=1560,res=300)
pheatmap(mat,cluster_cols=FALSE,cluster_rows=TRUE,show_colnames=FALSE,fontsize_row=6,
         annotation_col=ann,annotation_colors=anncols,color=colorRampPalette(c("#2166AC","white","#B2182B"))(101),
         main="Tier-1 pathway dependency (gene-wise z-score median)",border_color=NA)
dev.off()

topvar <- pathsum[order(-MAD_primary)][1:min(10,.N),pathway_major]
p2 <- ggplot(pathmat[pathway_major %in% topvar],aes(cancer_type,tier1_z_median_score,fill=cancer_type))+
  geom_violin(scale="width",trim=TRUE,color="grey35",linewidth=.25)+geom_boxplot(width=.16,outlier.shape=NA,fill="white",linewidth=.25)+
  facet_wrap(~pathway_major,scales="free_y",ncol=5)+scale_fill_manual(values=pal)+
  labs(x=NULL,y="Relative dependency score",title="Pathway dependency distributions (10 most variable)")+
  theme(legend.position="none",axis.text.x=element_text(angle=45,hjust=1),strip.text=element_text(size=6))
save_plot(p2,"P2_F02_pathway_dependency_distribution",7.2,5.2)

gdplot <- copy(gd); gdplot[, neglog10FDR := -log10(pmax(BH_FDR,1e-300))]
gdplot[, significant := BH_FDR<.05]
labg <- gdplot[order(BH_FDR,-abs(rank_biserial_CRC_vs_ESCA))][1:min(15,.N)]
p3 <- ggplot(gdplot,aes(rank_biserial_CRC_vs_ESCA,neglog10FDR,color=significant))+
  geom_point(size=.8,alpha=.65)+geom_hline(yintercept=-log10(.05),linetype=2,color="grey40")+
  geom_text(data=labg,aes(label=gene),size=2,vjust=-.5,check_overlap=TRUE)+
  scale_color_manual(values=c(`FALSE`="grey65",`TRUE`="#D55E00"))+
  labs(x="Rank-biserial effect (CRC vs ESCA; negative = stronger CRC dependency)",y="−log10(BH FDR)",title="Metabolic gene differential dependency")+
  theme(legend.position="none")
save_plot(p3,"P2_F03_CRC_ESCA_gene_volcano",6.2,4.2)

p4d <- pdiff[order(hodges_lehmann_shift)]; p4d[, pathway_major:=factor(pathway_major,levels=pathway_major)]
p4 <- ggplot(p4d,aes(hodges_lehmann_shift,pathway_major,color=BH_FDR<.05))+
  geom_vline(xintercept=0,color="grey60")+geom_errorbarh(aes(xmin=CI95_low,xmax=CI95_high),height=.18,linewidth=.35)+geom_point(size=2)+
  scale_color_manual(values=c(`FALSE`="grey55",`TRUE`="#0072B2"))+
  labs(x="Hodges-Lehmann shift: CRC - ESCA (95% CI)",y=NULL,title="Tier-1 pathway effects")+theme(legend.position="none")
save_plot(p4,"P2_F04_CRC_ESCA_pathway_effects",6.2,5.2)

p5 <- ggplot(sel,aes(median_LFC,MAD,color=dependency_class,size=binary_depleted_fraction))+
  geom_point(alpha=.7)+geom_vline(xintercept=strong_cutoff,linetype=2,color="grey45")+
  scale_color_manual(values=c(A_common_metabolic_essential="#D55E00",B_selective_metabolic_dependency="#0072B2",C_weak_or_non_dependent="grey70"))+
  scale_size(range=c(.6,3.2))+
  labs(x="Median library-adjusted LFC",y="MAD",color="Dependency class",size="Depleted fraction",title="Metabolic dependency strength and selectivity")+
  theme(legend.position="right")
save_plot(p5,"P2_F05_selective_dependency_landscape",6.8,4.5)

# Machine-readable run metadata.
dual <- cons[n_libraries==2 & is.finite(minLib_LFC) & is.finite(v1_1_LFC)]
dual_spearman <- cor(dual$minLib_LFC,dual$v1_1_LFC,method="spearman")
dual_median_abs_difference <- median(abs(dual$minLib_LFC-dual$v1_1_LFC))
run_qc <- data.table(metric=c("strong_dependency_threshold","essential_control_median","nonessential_control_median","consensus_rows","consensus_models","dual_library_models","dual_library_gene_pairs","dual_library_spearman","dual_library_median_absolute_LFC_difference","tier1_pathways_LT5","tier1_pathways_GT150","gene_tests","pathway_tests","random_seed"),
  value=c(strong_cutoff,ess_med,non_med,nrow(cons),uniqueN(cons$sample_ID),uniqueN(cons[n_libraries==2,sample_ID]),
          nrow(dual),dual_spearman,dual_median_abs_difference,
          uniqueN(qc[tier1_size_flag=="LT5",pathway_major]),uniqueN(qc[tier1_size_flag=="GT150",pathway_major]),nrow(gd),nrow(pdiff),20260901))
write_tsv(run_qc,"07_logs/P2_STAGE2_SCOUT_RUN_QC.tsv")

session <- capture.output(sessionInfo())
writeLines(session,file.path(root,"07_logs/P2_STAGE2_SESSION_INFO.txt"))
