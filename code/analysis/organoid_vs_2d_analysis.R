options(stringsAsFactors=FALSE, warn=1)
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(scales)})
script_arg <- grep("^--file=", commandArgs(trailingOnly=FALSE), value=TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else "."
repo_default <- normalizePath(file.path(dirname(script_path), "../.."), mustWork=FALSE)
root <- normalizePath(Sys.getenv("ORGANOID_PROJECT_ROOT", unset=repo_default), mustWork=FALSE)
tab <- normalizePath(Sys.getenv("ORGANOID_DERIVED_ROOT", unset=file.path(root,"data/derived/analysis_outputs")), mustWork=FALSE)
fig <- normalizePath(Sys.getenv("ORGANOID_FIGURE_ROOT", unset=file.path(root,"figures/generated")), mustWork=FALSE)
source_root <- normalizePath(Sys.getenv("ORGANOID_SOURCE_ROOT", unset=file.path(root,"data/source")), mustWork=FALSE)
metadata_root <- normalizePath(Sys.getenv("ORGANOID_METADATA_ROOT", unset=file.path(root,"data/external_metadata")), mustWork=FALSE)
dir.create(tab, recursive=TRUE, showWarnings=FALSE); dir.create(fig, recursive=TRUE, showWarnings=FALSE)

diff <- fread(file.path(tab,"CRC_ESCA_GENE_DIFFERENTIAL_DEPENDENCY.tsv"))
sel <- fread(file.path(tab,"SELECTIVE_METABOLIC_DEPENDENCIES.tsv"))
short <- fread(file.path(tab,"P2_SCOUT_SHORTLIST.tsv"))
cons <- fread(file.path(source_root,"02_crispr/METABOLIC_ORGANOID_GENE_CONSENSUS.tsv"))
tier <- fread(file.path(metadata_root,"METABOLIC_GENE_UNIVERSE_TIERED.tsv"))

mw <- function(x,y) {
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  if(length(x)<2 || length(y)<2) return(list(U=NA_real_, p=NA_real_, rbc=NA_real_, meddiff=NA_real_))
  z <- suppressWarnings(wilcox.test(x,y,exact=FALSE,correct=FALSE))
  U <- unname(z$statistic)
  list(U=U,p=z$p.value,rbc=2*U/(length(x)*length(y))-1,meddiff=median(x)-median(y))
}

# OXPHOS gene-level decomposition is required to identify drivers of the 299-gene aggregate.
oxgenes <- unique(tier[pathway_major=="OXPHOS" & tier=="Tier 1", gene_symbol])
oxd <- diff[gene %in% oxgenes]
setorder(oxd, BH_FDR, rank_biserial_CRC_vs_ESCA)
oxd[, `:=`(
  oxphos_alignment=ifelse(rank_biserial_CRC_vs_ESCA < 0,"aligned_CRC_stronger",ifelse(rank_biserial_CRC_vs_ESCA>0,"opposing_ESCA_stronger","neutral")),
  driver_score=pmax(0,-rank_biserial_CRC_vs_ESCA) * abs(median_difference_CRC_minus_ESCA),
  submodule="other_respiratory_or_linked"
)]
oxd[grepl("^NDUF|^NUBPL$|^TMEM126|^ECSIT$|^ACAD9$|^FOXRED1$|^TIMMDC1$",gene),submodule:="complex_I"]
oxd[grepl("^SDH",gene),submodule:="complex_II"]
oxd[grepl("^UQCR|^CYC1$|^BCS1L$|^TTC19$",gene),submodule:="complex_III"]
oxd[grepl("^COX|^SCO[12]$|^SURF1$|^COA[3-9]|^PET100$|^PET117$|^TACO1$",gene),submodule:="complex_IV"]
oxd[grepl("^ATP5|^ATPIF1$|^TMEM70$",gene),submodule:="complex_V_ATP_synthase"]
oxd[grepl("^MRPL|^MRPS|^MTIF|^MTERF|^TSFM$|^TUFM$|^GFM|^LRPPRC$|^FASTKD",gene),submodule:="mitochondrial_translation_expression"]
oxd[grepl("^SLC25|^TIMM|^TOMM|^PMPCA$|^PMPCB$|^AFG3L2$|^SPG7$",gene),submodule:="mitochondrial_transport_proteostasis"]
oxd[, driver_rank:=frank(-driver_score,ties.method="min")]
oxd[, driver_category:=fcase(BH_FDR<0.05 & rank_biserial_CRC_vs_ESCA<0,"primary_FDR_driver",
                             rank_biserial_CRC_vs_ESCA<=-0.20,"supporting_aligned_driver",
                             rank_biserial_CRC_vs_ESCA>=0.20,"opposing_gene",
                             default="weak_or_neutral")]

sig33 <- diff[BH_FDR<0.05,gene]
bshort <- intersect(sel[grepl("^B_",dependency_class),gene],short$gene)
ox_short <- intersect(short[grepl("OXPHOS",pathway),gene],oxd$gene)
topox <- oxd[gene %in% ox_short][order(-driver_score)][1:min(15,.N),gene]
candidates <- unique(c(sig33,bshort,topox))
basis <- data.table(gene=candidates)
basis[, candidate_basis:=vapply(gene,function(g) paste(c(if(g %in% sig33) "CRC_ESCA_FDR33",
                                                          if(g %in% bshort) "class_B_in_scout_shortlist",
                                                          if(g %in% topox) "top_OXPHOS_interpreter"),collapse=";"),character(1))]

load(file.path(source_root,"00_source/figshare_28339340/07_EssMatrix_qnorm_corrected_logFCs.RData"))
load(file.path(source_root,"00_source/figshare_28339340/07_EssMatrix_bDepletionsB2.Rdata"))
load(file.path(source_root,"00_source/figshare_28339340/10_PANCANCER_coreFitness_genes.RData"))
model <- fread(file.path(metadata_root,"DepMap_24Q2_Model.csv"))
model <- model[ModelID %in% colnames(qnorm_corrected_logFCs)]
crc_ids <- model[OncotreePrimaryDisease=="Colorectal Adenocarcinoma",ModelID]
esca_ids <- model[grepl("^Esophageal",OncotreePrimaryDisease) | OncotreeCode %in% c("ESCC","ESCA"),ModelID]
related_es_ids <- model[OncotreeLineage=="Esophagus/Stomach",ModelID]
annot_ids <- model$ModelID

ice <- fread(file.path(source_root,"00_source/figshare_28339340/CRISPRInferredCommonEssentials.csv"),header=TRUE)
icegenes <- sub(" \\(.*$","",ice[[1]])
all_genes <- intersect(candidates,rownames(qnorm_corrected_logFCs))

cmp <- vector("list",length(candidates))
for(i in seq_along(candidates)) {
  g <- candidates[i]
  drow <- diff[gene==g][1]
  srow <- sel[gene==g][1]
  direction <- if(is.finite(drow$median_difference_CRC_minus_ESCA) && drow$median_difference_CRC_minus_ESCA<0) "CRC-leaning" else "ESCA-leaning"
  org_all <- cons[gene==g,adjusted_consensus_LFC]
  org_crc <- cons[gene==g & cancer_type=="Colorectal",adjusted_consensus_LFC]
  org_es <- cons[gene==g & cancer_type=="Oesophageal",adjusted_consensus_LFC]
  if(g %in% rownames(qnorm_corrected_logFCs)) {
    ext_all <- as.numeric(qnorm_corrected_logFCs[g,])
    ext_crc <- as.numeric(qnorm_corrected_logFCs[g,intersect(crc_ids,colnames(qnorm_corrected_logFCs)),drop=TRUE])
    ext_es <- as.numeric(qnorm_corrected_logFCs[g,intersect(esca_ids,colnames(qnorm_corrected_logFCs)),drop=TRUE])
    ext_rel <- as.numeric(qnorm_corrected_logFCs[g,intersect(related_es_ids,colnames(qnorm_corrected_logFCs)),drop=TRUE])
    ext_ann <- as.numeric(qnorm_corrected_logFCs[g,intersect(annot_ids,colnames(qnorm_corrected_logFCs)),drop=TRUE])
    dep_all <- as.numeric(bDepletionsB2[g,])
    dep_match <- if(direction=="CRC-leaning") as.numeric(bDepletionsB2[g,intersect(crc_ids,colnames(bDepletionsB2)),drop=TRUE]) else as.numeric(bDepletionsB2[g,intersect(esca_ids,colnames(bDepletionsB2)),drop=TRUE])
  } else { ext_all<-ext_crc<-ext_es<-ext_rel<-ext_ann<-dep_all<-dep_match<-numeric() }
  org_match <- if(direction=="CRC-leaning") org_crc else org_es
  ext_match <- if(direction=="CRC-leaning") ext_crc else ext_es
  tmatch <- mw(org_match,ext_match); tglobal <- mw(org_all,ext_all)
  line_crc <- mw(ext_crc,ext_ann[!intersect(annot_ids,colnames(qnorm_corrected_logFCs)) %in% crc_ids])
  line_es <- mw(ext_es,ext_ann[!intersect(annot_ids,colnames(qnorm_corrected_logFCs)) %in% esca_ids])
  cmp[[i]] <- data.table(gene=g,pathway=if(nrow(short[gene==g])) short[gene==g,pathway][1] else drow$pathways,
    candidate_basis=basis[gene==g,candidate_basis],organoid_direction="negative_LFC_is_stronger_dependency",CRC_ESCA_direction=direction,
    organoid_effect_global_median=median(org_all,na.rm=TRUE),organoid_effect_lineage_median=median(org_match,na.rm=TRUE),
    cellline_effect_global_median=if(length(ext_all)) median(ext_all,na.rm=TRUE) else NA_real_,
    cellline_effect_matched_median=if(length(ext_match)) median(ext_match,na.rm=TRUE) else NA_real_,
    cellline_effect_related_ES_stomach_median=if(length(ext_rel)) median(ext_rel,na.rm=TRUE) else NA_real_,
    organoid_vs_cellline_difference=tmatch$meddiff,matched_U=tmatch$U,matched_rank_biserial=tmatch$rbc,matched_raw_P=tmatch$p,
    global_organoid_vs_cellline_difference=tglobal$meddiff,global_rank_biserial=tglobal$rbc,global_raw_P=tglobal$p,
    n_organoid_lineage=length(org_match),n_cellline_matched=length(ext_match),n_cellline_related_ES_stomach=length(ext_rel),n_cellline_global=length(ext_all),
    cellline_binary_depleted_fraction_global=if(length(dep_all)) mean(dep_all,na.rm=TRUE) else NA_real_,
    cellline_binary_depleted_fraction_matched=if(length(dep_match)) mean(dep_match,na.rm=TRUE) else NA_real_,
    common_essential_status=if(g %in% PanCancerCoreFitnessGenes && g %in% icegenes) "both_core_fitness_lists" else if(g %in% PanCancerCoreFitnessGenes) "Sanger_pan_cancer_core_fitness" else if(g %in% icegenes) "DepMap_inferred_common_essential" else "not_flagged",
    selectivity=if(nrow(srow)) srow$dependency_class else NA_character_,organoid_selectivity_MAD=if(nrow(srow)) srow$MAD else NA_real_,
    CRC_cellline_lineage_rank_biserial=line_crc$rbc,CRC_cellline_lineage_raw_P=line_crc$p,
    ESCA_cellline_lineage_rank_biserial=line_es$rbc,ESCA_cellline_lineage_raw_P=line_es$p,
    organoid_CRC_ESCA_rank_biserial=drow$rank_biserial_CRC_vs_ESCA,organoid_CRC_ESCA_FDR=drow$BH_FDR)
}
cmp <- rbindlist(cmp,fill=TRUE)
cmp[, matched_BH_FDR:=p.adjust(matched_raw_P,"BH")]
cmp[, global_BH_FDR:=p.adjust(global_raw_P,"BH")]
cmp[, CRC_cellline_lineage_BH_FDR:=p.adjust(CRC_cellline_lineage_raw_P,"BH")]
cmp[, ESCA_cellline_lineage_BH_FDR:=p.adjust(ESCA_cellline_lineage_raw_P,"BH")]
cmp[, comparison_category:=fcase(is.na(matched_rank_biserial),"not_testable",
  matched_BH_FDR<0.05 & matched_rank_biserial<=-0.20,"stronger_in_organoids",
  matched_BH_FDR<0.05 & matched_rank_biserial>=0.20,"weaker_in_organoids",
  default="similar_or_not_resolved")]
cmp[, lineage_external_aligned:=fcase(CRC_ESCA_direction=="CRC-leaning" & CRC_cellline_lineage_BH_FDR<0.05 & CRC_cellline_lineage_rank_biserial<0,TRUE,
                                     CRC_ESCA_direction=="ESCA-leaning" & ESCA_cellline_lineage_BH_FDR<0.05 & ESCA_cellline_lineage_rank_biserial<0,TRUE,default=FALSE)]
cmp[, external_support_category:=fcase(comparison_category=="stronger_in_organoids","candidate_organoid_enhanced_dependency",
  lineage_external_aligned,"lineage_associated_dependency",
  common_essential_status!="not_flagged" | cellline_binary_depleted_fraction_global>=0.8,"general_common_fitness_dependency",
  comparison_category=="weaker_in_organoids","cellline_stronger_dependency",
  default="limited_or_mixed_external_support")]
cmp[, QC_warning:=paste0("cross-platform LFC scales; matched cell-line n=",n_cellline_matched,
                         "; 576/930 ACH models annotated by DepMap 24Q2; non-significant classified as similar_or_not_resolved")]
setcolorder(cmp,c("gene","pathway","candidate_basis","organoid_direction","CRC_ESCA_direction","organoid_effect_global_median","organoid_effect_lineage_median","cellline_effect_global_median","cellline_effect_matched_median","organoid_vs_cellline_difference","matched_rank_biserial","matched_raw_P","matched_BH_FDR","comparison_category","common_essential_status","selectivity","external_support_category","QC_warning"))
cmp[, .abs_matched_rbc:=abs(matched_rank_biserial)]
setorder(cmp,matched_BH_FDR,-.abs_matched_rbc)
cmp[, .abs_matched_rbc:=NULL]
fwrite(cmp,file.path(tab,"P2_3A_DEPMAP_CANDIDATE_COMPARISON.tsv"),sep="\t",na="NA")

oxout <- merge(oxd[,.(gene,ensembl_id,submodule,n_CRC,n_ESCA,median_CRC,median_ESCA,median_difference_CRC_minus_ESCA,
                       rank_biserial_CRC_vs_ESCA,raw_P,BH_FDR,oxphos_alignment,driver_score,driver_rank,driver_category)],
               cmp[,.(gene,candidate_external_tested=TRUE,common_essential_status,cellline_effect_global_median,cellline_binary_depleted_fraction_global,
                       CRC_cellline_lineage_rank_biserial,CRC_cellline_lineage_BH_FDR,ESCA_cellline_lineage_rank_biserial,ESCA_cellline_lineage_BH_FDR)],
               by="gene",all.x=TRUE)
oxout[is.na(candidate_external_tested),candidate_external_tested:=FALSE]
setorder(oxout,driver_rank)
fwrite(oxout,file.path(tab,"P2_3A_OXPHOS_DRIVER_GENES.tsv"),sep="\t",na="NA")

# Conservative <=10 candidate shortlist: reward organoid enhancement, reproducible lineage direction,
# prior organoid FDR/selectivity and OXPHOS driver status; penalize generic common fitness.
prio <- copy(cmp)
prio[, ox_driver:=gene %in% oxout[driver_category %in% c("primary_FDR_driver","supporting_aligned_driver") & driver_rank<=20,gene]]
prio[, priority_score:=3*(comparison_category=="stronger_in_organoids") + 2*lineage_external_aligned +
       2*(organoid_CRC_ESCA_FDR<0.05) + 2*grepl("^B_",selectivity) + 1*ox_driver -
       2*(common_essential_status!="not_flagged")]
prio[, evidence_count:=(comparison_category=="stronger_in_organoids")+lineage_external_aligned+(organoid_CRC_ESCA_FDR<0.05)+grepl("^B_",selectivity)+ox_driver]
prio <- prio[evidence_count>=2]
prio[, .abs_org_rbc:=abs(organoid_CRC_ESCA_rank_biserial)]
setorder(prio,-priority_score,matched_BH_FDR,-.abs_org_rbc,gene)
prio[, .abs_org_rbc:=NULL]
prio <- prio[1:min(10,.N)]
prio[, priority:=ifelse(seq_len(.N)<=5,"high","medium")]
outshort <- prio[,.(gene,pathway,organoid_direction,CRC_ESCA_direction,
  organoid_effect=organoid_effect_lineage_median,cell_line_effect=cellline_effect_matched_median,
  organoid_vs_cell_line_difference=organoid_vs_cellline_difference,common_essential_status,selectivity,
  external_support_category,QC_warning,priority,priority_score)]
fwrite(outshort,file.path(tab,"P2_3A_EXTERNAL_VALIDATION_SHORTLIST.tsv"),sep="\t",na="NA")

# Figure: candidate matched-platform effects, lineage concordance, and OXPHOS gene drivers.
plotp <- prio
plotp[, gene:=factor(gene,levels=rev(gene))]
cols <- c("stronger_in_organoids"="#0072B2","similar_or_not_resolved"="#999999","weaker_in_organoids"="#D55E00","not_testable"="#000000")
p1 <- ggplot(plotp,aes(x=matched_rank_biserial,y=gene,color=comparison_category))+
  geom_vline(xintercept=c(-.2,0,.2),linetype=c(3,2,3),linewidth=.3,color="grey55")+geom_point(size=2.2)+
  scale_color_manual(values=cols,drop=FALSE)+labs(x="Organoid vs matched 2D rank-biserial\n(negative = organoid stronger)",y=NULL,color=NULL,title="a  Prioritized candidates")+
  theme_bw(base_size=8)+theme(legend.position="bottom",panel.grid.minor=element_blank())
plotall <- cmp[is.finite(organoid_CRC_ESCA_rank_biserial)]
plotall[, ext_lineage:=ifelse(CRC_ESCA_direction=="CRC-leaning",CRC_cellline_lineage_rank_biserial,-ESCA_cellline_lineage_rank_biserial)]
plotall[, label:=ifelse(gene %in% as.character(prio$gene),gene,"")]
plotall[, highlighted:=ifelse(label=="","other candidates","prioritized")]
p2 <- ggplot(plotall,aes(x=organoid_CRC_ESCA_rank_biserial,y=ext_lineage,color=highlighted))+
  geom_hline(yintercept=0,linewidth=.3,color="grey65")+geom_vline(xintercept=0,linewidth=.3,color="grey65")+
  geom_point(alpha=.65,size=1.5)+geom_text(aes(label=label),check_overlap=TRUE,size=2.2,vjust=-.5)+
  scale_color_manual(values=c("other candidates"="#999999","prioritized"="#D55E00"))+guides(color="none")+
  labs(x="Organoid CRC vs ESCA rank-biserial\n(negative = CRC stronger)",y="Aligned 2D lineage effect\n(negative = matched lineage stronger)",title="b  Lineage concordance")+
  theme_bw(base_size=8)+theme(panel.grid.minor=element_blank())
oxplot <- oxout[oxphos_alignment=="aligned_CRC_stronger"][order(driver_rank)][1:min(15,.N)]
oxplot[, gene:=factor(gene,levels=rev(gene))]
p3 <- ggplot(oxplot,aes(x=-rank_biserial_CRC_vs_ESCA,y=gene,color=submodule))+
  geom_segment(aes(x=0,xend=-rank_biserial_CRC_vs_ESCA,yend=gene),linewidth=.5)+geom_point(size=2)+
  labs(x="CRC-strength aligned effect (-rank-biserial)",y=NULL,color="Submodule",title="c  OXPHOS gene drivers")+
  theme_bw(base_size=8)+theme(legend.position="bottom",panel.grid.minor=element_blank())
suppressPackageStartupMessages(library(grid))
pngfile <- file.path(fig,"P2_F06_organoid_vs_cellline_candidates.png")
pdffile <- file.path(fig,"P2_F06_organoid_vs_cellline_candidates.pdf")
drawfig <- function(device=NULL){
  if(!is.null(device)) device()
  grid.newpage(); lay<-grid.layout(2,2,widths=unit(c(1,1),"null"),heights=unit(c(1,1.15),"null")); pushViewport(viewport(layout=lay))
  print(p1,vp=viewport(layout.pos.row=1,layout.pos.col=1)); print(p2,vp=viewport(layout.pos.row=1,layout.pos.col=2)); print(p3,vp=viewport(layout.pos.row=2,layout.pos.col=1:2)); popViewport()
  if(!is.null(device)) dev.off()
}
drawfig(function() png(pngfile,width=3600,height=3000,res=450,type="cairo"))
drawfig(function() cairo_pdf(pdffile,width=8,height=6.7))

cat(sprintf("candidates=%d sig33=%d Bshort=%d topOX=%d shortlist=%d ox=%d annotated=%d crcCL=%d escaCL=%d related=%d\n",
            nrow(cmp),length(sig33),length(bshort),length(topox),nrow(outshort),nrow(oxout),length(annot_ids),length(crc_ids),length(esca_ids),length(related_es_ids)))
cat("SHORTLIST\n"); print(outshort)
