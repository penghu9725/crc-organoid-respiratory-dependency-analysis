options(stringsAsFactors=FALSE,warn=1)
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(scales);library(pheatmap);library(grid)})
set.seed(20260901)
script_arg<-grep("^--file=",commandArgs(trailingOnly=FALSE),value=TRUE);script_path<-if(length(script_arg))sub("^--file=","",script_arg[1])else"."
repo_default<-normalizePath(file.path(dirname(script_path),"../.."),mustWork=FALSE)
root<-normalizePath(Sys.getenv("ORGANOID_PROJECT_ROOT",unset=repo_default),mustWork=FALSE)
tab<-normalizePath(Sys.getenv("ORGANOID_DERIVED_ROOT",unset=file.path(root,"data/derived/analysis_outputs")),mustWork=FALSE)
fig<-normalizePath(Sys.getenv("ORGANOID_FIGURE_ROOT",unset=file.path(root,"figures/generated")),mustWork=FALSE)
source_root<-normalizePath(Sys.getenv("ORGANOID_SOURCE_ROOT",unset=file.path(root,"data/source")),mustWork=FALSE)
metadata_root<-normalizePath(Sys.getenv("ORGANOID_METADATA_ROOT",unset=file.path(root,"data/external_metadata")),mustWork=FALSE)
dir.create(tab,recursive=TRUE,showWarnings=FALSE);dir.create(fig,recursive=TRUE,showWarnings=FALSE)
mods<-c("CI_structural","CI_assembly","CII","CIII_structural","CIII_assembly","CIV_structural","CIV_assembly","CV_structural","CV_assembly","FeS_biogenesis","CoQ_biosynthesis","mtFAS","mtDNA_maintenance_expression")
cands<-c("NDUFAF8","NDUFA6","COA5","COX6B1","TIMMDC1","MCAT","ISCA2")
map<-fread(file.path(metadata_root,"OXPHOS_SUBMODULE_MAP.tsv"))
map<-unique(map[submodule_primary%in%mods & !grepl("^MT-",gene_symbol),.(gene=gene_symbol,module=submodule_primary,source,evidence_level,notes)])
cons_all<-fread(file.path(source_root,"02_crispr/METABOLIC_ORGANOID_GENE_CONSENSUS.tsv"))
cons<-cons_all[cancer_type%in%c("Colorectal","Oesophageal")]
diff<-fread(file.path(tab,"CRC_ESCA_GENE_DIFFERENTIAL_DEPENDENCY.tsv"))
ext<-fread(file.path(tab,"P2_3A_DEPMAP_CANDIDATE_COMPARISON.tsv"))

mw<-function(x,y){x<-x[is.finite(x)];y<-y[is.finite(y)];if(length(x)<2||length(y)<2)return(list(U=NA,p=NA,rbc=NA,md=NA));z<-wilcox.test(x,y,exact=FALSE,correct=FALSE);U<-unname(z$statistic);list(U=U,p=z$p.value,rbc=2*U/(length(x)*length(y))-1,md=median(x)-median(y))}
boot_ci<-function(x,y,B=2000){x<-x[is.finite(x)];y<-y[is.finite(y)];if(length(x)<2||length(y)<2)return(c(NA,NA));z<-replicate(B,median(sample(x,length(x),TRUE))-median(sample(y,length(y),TRUE)));quantile(z,c(.025,.975),na.rm=TRUE,names=FALSE)}

# Higher dependency_z means stronger dependency. Standardization occurs gene-wise across the 144 primary-comparison organoids.
cg<-cons[gene%in%map$gene,.(sample_ID,gene,cancer_type,lfc=adjusted_consensus_LFC)]
cg[,dependency_z:=as.numeric(scale(-lfc)),by=gene]
cg<-cg[is.finite(dependency_z)]
cg<-merge(cg,map[,.(gene,module)],by="gene",allow.cartesian=TRUE)
cg[,dependency_rank:=frank(dependency_z,ties.method="average")/.N,by=gene]
scores<-cg[,.(module_score=median(dependency_z),rank_module_score=median(dependency_rank)),by=.(module,sample_ID,cancer_type)]

gene_stats<-cg[,{
 x<-dependency_z[cancer_type=="Colorectal"];y<-dependency_z[cancer_type=="Oesophageal"];z<-mw(x,y)
 .(n_CRC=length(x),n_ESCA=length(y),median_dependency_z_CRC=median(x),median_dependency_z_ESCA=median(y),median_difference_CRC_minus_ESCA=z$md,rank_biserial_CRC_vs_ESCA=z$rbc,raw_P=z$p)
},by=.(module,gene)]
gene_stats[,BH_FDR_within_all_mapped_genes:=p.adjust(raw_P,"BH")]
gene_stats<-merge(gene_stats,diff[,.(gene,gene_level_prior_BH_FDR=BH_FDR,gene_level_prior_rank_biserial=rank_biserial_CRC_vs_ESCA)],by="gene",all.x=TRUE)
gene_stats[,direction:=fcase(median_difference_CRC_minus_ESCA>0,"CRC-leaning",median_difference_CRC_minus_ESCA<0,"ESCA-leaning",default="neutral")]
gene_stats[,gene_level_FDR_driver:=!is.na(gene_level_prior_BH_FDR)&gene_level_prior_BH_FDR<.05&gene_level_prior_rank_biserial<0]

res<-scores[,{
 x<-module_score[cancer_type=="Colorectal"];y<-module_score[cancer_type=="Oesophageal"];z<-mw(x,y);ci<-boot_ci(x,y)
 xr<-rank_module_score[cancer_type=="Colorectal"];yr<-rank_module_score[cancer_type=="Oesophageal"];zr<-mw(xr,yr)
 gs<-gene_stats[module==.BY$module]
 .(n_evaluable_genes=uniqueN(cg[module==.BY$module,gene]),n_CRC=length(x),n_ESCA=length(y),CRC_median=median(x),ESCA_median=median(y),median_difference_CRC_minus_ESCA=z$md,effect_size_rank_biserial=z$rbc,CI95_low=ci[1],CI95_high=ci[2],U=z$U,raw_P=z$p,rank_sensitivity_effect=zr$rbc,rank_sensitivity_raw_P=zr$p,n_CRC_leaning_genes=sum(gs$direction=="CRC-leaning"),proportion_CRC_leaning_genes=mean(gs$direction=="CRC-leaning"),n_gene_level_FDR_drivers=sum(gs$gene_level_FDR_driver),proportion_gene_level_FDR_drivers=mean(gs$gene_level_FDR_driver),directional_consistency=mean(gs$direction=="CRC-leaning"))
},by=module]
res[,BH_FDR:=p.adjust(raw_P,"BH")];res[,rank_sensitivity_BH_FDR:=p.adjust(rank_sensitivity_raw_P,"BH")]
res[,interpretation:=fcase(BH_FDR<.05&median_difference_CRC_minus_ESCA>0&n_evaluable_genes>=5,"CRC_stronger_significant",
                           BH_FDR<.05&median_difference_CRC_minus_ESCA>0&n_evaluable_genes<5,"CRC_stronger_significant_low_gene_coverage",
                           BH_FDR<.05&median_difference_CRC_minus_ESCA<0&n_evaluable_genes>=5,"ESCA_stronger_significant",
                           BH_FDR<.05&median_difference_CRC_minus_ESCA<0&n_evaluable_genes<5,"ESCA_stronger_significant_low_gene_coverage",default="not_significant")]
setorder(res,BH_FDR)
fwrite(res,file.path(tab,"P2_3B_OXPHOS_SUBMODULE_RESULTS.tsv"),sep="\t",na="NA")
gene_stats[,abs_median_difference:=abs(median_difference_CRC_minus_ESCA)]
setorder(gene_stats,module,gene_level_prior_BH_FDR,-abs_median_difference)
gene_stats[,abs_median_difference:=NULL]
fwrite(gene_stats,file.path(tab,"P2_3B_OXPHOS_SUBMODULE_GENE_DRIVERS.tsv"),sep="\t",na="NA")

# Focused biomarker analysis
rna<-fread(file.path(source_root,"00_source/figshare_28339340/data_RNAseq_ALL-organoids.csv"),select=c("sample_ID",paste0(cands,"_exp")))
dep7<-cons_all[gene%in%cands,.(sample_ID,gene,cancer_type,dependency_score=-adjusted_consensus_LFC)]
expr_long<-melt(rna,id.vars="sample_ID",variable.name="feature",value.name="expression")
expr_long[,gene:=sub("_exp$","",feature)];expr_long[,feature:=NULL]
expr_long<-merge(dep7,expr_long,by=c("sample_ID","gene"),all.x=TRUE)
expr_long[,expression_rank_within_lineage:=frank(expression,ties.method="average")/.N,by=.(gene,cancer_type)]
expr_long[,dependency_rank_within_lineage:=frank(dependency_score,ties.method="average")/.N,by=.(gene,cancer_type)]
expr_res<-expr_long[,{
 ok<-is.finite(expression)&is.finite(dependency_score);zp<-if(sum(ok)>=3)cor.test(expression[ok],dependency_score[ok],method="spearman",exact=FALSE) else NULL
 za<-if(sum(ok)>=3)cor.test(expression_rank_within_lineage[ok],dependency_rank_within_lineage[ok],method="spearman",exact=FALSE) else NULL
 .(N=sum(ok),spearman_rho=if(is.null(zp))NA_real_ else unname(zp$estimate),raw_P=if(is.null(zp))NA_real_ else zp$p.value,
   lineage_adjusted_spearman_rho=if(is.null(za))NA_real_ else unname(za$estimate),lineage_adjusted_raw_P=if(is.null(za))NA_real_ else za$p.value,
   median_expression=median(expression[ok]),dependency_definition="higher score = stronger dependency")
},by=gene]
expr_res[,BH_FDR:=p.adjust(raw_P,"BH")];expr_res[,lineage_adjusted_BH_FDR:=p.adjust(lineage_adjusted_raw_P,"BH")]
expr_res[,evidence:=ifelse(lineage_adjusted_BH_FDR<.05,ifelse(lineage_adjusted_spearman_rho>0,"higher_expression_stronger_dependency_lineage_adjusted","higher_expression_weaker_dependency_lineage_adjusted"),"no_lineage_adjusted_FDR_support")]
fwrite(expr_res,file.path(tab,"P2_3B_CANDIDATE_EXPRESSION_DEPENDENCY.tsv"),sep="\t",na="NA")

mut<-fread(file.path(source_root,"00_source/figshare_28339340/data_driver-mutations_ALL-organoids.csv"),select=c("sample_ID",paste0(c("KRAS","APC","TP53","PIK3CA","BRAF"),"_mut")))
cnv<-fread(file.path(source_root,"00_source/figshare_28339340/data_CNV-genes_ALL-organoids.csv"),select=c("sample_ID","wholeGenomeDuplication","focal_amplifical"))
cov<-fread(file.path(source_root,"00_source/figshare_28339340/data_covariates_ALL-organoids.csv"),select=c("sample_ID","ploidy","msStatus"))
bio<-Reduce(function(x,y)merge(x,y,by="sample_ID",all=TRUE),list(mut,cnv,cov))
bio[,BRAF_KRAS_mut:=as.integer(BRAF_mut==1|KRAS_mut==1)]
catfeat<-c("KRAS_mut","APC_mut","TP53_mut","PIK3CA_mut","BRAF_mut","BRAF_KRAS_mut","msStatus","wholeGenomeDuplication","focal_amplifical")
gres<-list();k<-1
for(g in cands){d<-merge(dep7[gene==g],bio,by="sample_ID",all.x=TRUE);d[,dependency_lineage_residual:=dependency_score-median(dependency_score,na.rm=TRUE),by=cancer_type]
 for(f in catfeat){v<-d[[f]];ok<-is.finite(v)&is.finite(d$dependency_score);n1<-sum(v[ok]==1);n0<-sum(v[ok]==0);adequate<-n1>=5&n0>=5
  zp<-if(adequate)mw(d$dependency_score[ok&v==1],d$dependency_score[ok&v==0]) else list(U=NA,p=NA,rbc=NA,md=NA);za<-if(adequate)mw(d$dependency_lineage_residual[ok&v==1],d$dependency_lineage_residual[ok&v==0]) else list(U=NA,p=NA,rbc=NA,md=NA)
  gres[[k]]<-data.table(gene=g,biomarker=f,biomarker_type=if(f%in%c("wholeGenomeDuplication","focal_amplifical"))"CNV_global_categorical" else if(f=="msStatus")"MSI_categorical" else "driver_mutation_categorical",N=sum(ok),group1_N=n1,group0_N=n0,pooled_effect=zp$rbc,pooled_raw_P=zp$p,effect=za$rbc,effect_definition="lineage-adjusted rank-biserial: positive = altered group stronger dependency",median_difference=za$md,raw_P=za$p,underpowered=!adequate,availability="available");k<-k+1}
 d[,ploidy_rank_within_lineage:=frank(ploidy,ties.method="average")/.N,by=cancer_type];d[,dependency_rank_within_lineage:=frank(dependency_score,ties.method="average")/.N,by=cancer_type]
 ok<-is.finite(d$ploidy)&is.finite(d$dependency_score);zp<-if(sum(ok)>=10)cor.test(d$ploidy[ok],d$dependency_score[ok],method="spearman",exact=FALSE) else NULL;za<-if(sum(ok)>=10)cor.test(d$ploidy_rank_within_lineage[ok],d$dependency_rank_within_lineage[ok],method="spearman",exact=FALSE) else NULL
 gres[[k]]<-data.table(gene=g,biomarker="ploidy",biomarker_type="CNV_continuous",N=sum(ok),group1_N=NA_integer_,group0_N=NA_integer_,pooled_effect=if(is.null(zp))NA_real_ else unname(zp$estimate),pooled_raw_P=if(is.null(zp))NA_real_ else zp$p.value,effect=if(is.null(za))NA_real_ else unname(za$estimate),effect_definition="lineage-adjusted Spearman rho: positive = higher ploidy, stronger dependency",median_difference=NA_real_,raw_P=if(is.null(za))NA_real_ else za$p.value,underpowered=sum(ok)<10,availability="available");k<-k+1
 gres[[k]]<-data.table(gene=g,biomarker="candidate_locus_copy_number",biomarker_type="CNV_gene_level_continuous",N=0L,group1_N=NA_integer_,group0_N=NA_integer_,effect=NA_real_,effect_definition="not estimable",median_difference=NA_real_,raw_P=NA_real_,underpowered=TRUE,availability="unavailable: frozen gene-level CNV table has no locus field for this candidate");k<-k+1
}
gres<-rbindlist(gres,fill=TRUE);gres[,BH_FDR:=p.adjust(raw_P,"BH")];gres[,evidence:=fcase(is.na(raw_P),"not_tested",BH_FDR<.05,"FDR_supported",underpowered,"underpowered",default="no_FDR_support")]
fwrite(gres,file.path(tab,"P2_3B_CANDIDATE_GENOMIC_BIOMARKERS.tsv"),sep="\t",na="NA")

# Integrated ranking
ranktab<-data.table(gene=cands)
ranktab<-merge(ranktab,map[gene%in%cands,.(gene,functional_module=module)][,.(functional_module=paste(unique(functional_module),collapse=";")),by=gene],by="gene",all.x=TRUE)
ranktab<-merge(ranktab,diff[,.(gene,CRC_vs_ESCA_effect=rank_biserial_CRC_vs_ESCA,CRC_vs_ESCA_FDR=BH_FDR)],by="gene",all.x=TRUE)
ranktab<-merge(ranktab,ext[,.(gene,organoid_vs_2D_category=comparison_category,common_essential_status)],by="gene",all.x=TRUE)
ranktab<-merge(ranktab,expr_res[,.(gene,RNA_dependency_rho=lineage_adjusted_spearman_rho,RNA_dependency_FDR=lineage_adjusted_BH_FDR,RNA_dependency_evidence=evidence)],by="gene",all.x=TRUE)
ranktab[,genomic_biomarker_evidence:=vapply(gene,function(g){x<-gres[gene==g&availability=="available"];if(any(x$BH_FDR<.05,na.rm=TRUE))paste(x[BH_FDR<.05,paste0(biomarker,"(FDR=",signif(BH_FDR,3),")")],collapse=";") else "no_FDR_supported_internal_genomic_biomarker"},character(1))]
ranktab<-merge(ranktab,res[,.(functional_module=module,submodule_effect=effect_size_rank_biserial,submodule_FDR=BH_FDR,submodule_interpretation=interpretation)],by="functional_module",all.x=TRUE)
ranktab[,novelty_risk:=ifelse(common_essential_status!="not_flagged","not_literature_assessed; elevated_common_fitness_confounding","not_literature_assessed; lower_common_fitness_confounding")]
ranktab[,submodule_support:=paste0(functional_module,": ",submodule_interpretation,"; FDR=",signif(submodule_FDR,3))]
ranktab[,QC_limitation:=vapply(seq_len(.N),function(i)paste(c("gene-level locus CNV unavailable",if(RNA_dependency_FDR[i]>=.05)"RNA association exploratory/non-significant",if(common_essential_status[i]!="not_flagged")"common-fitness confounding"),collapse="; "),character(1))]
ranktab[,evidence_score:=2*(CRC_vs_ESCA_FDR<.05)+2*(submodule_FDR<.05)+2*(organoid_vs_2D_category=="stronger_in_organoids")+1*(RNA_dependency_FDR<.05)+1*(genomic_biomarker_evidence!="no_FDR_supported_internal_genomic_biomarker")-1*(common_essential_status!="not_flagged")]
setorder(ranktab,-evidence_score,CRC_vs_ESCA_FDR,gene);ranktab[,final_priority:=fcase(evidence_score>=5,"high",evidence_score>=3,"medium",default="exploratory")]
fwrite(ranktab,file.path(tab,"P2_3B_FINAL_CANDIDATE_RANKING.tsv"),sep="\t",na="NA")

# F07 forest plot
rp<-copy(res);rp[,module:=factor(module,levels=rev(mods))]
p7<-ggplot(rp,aes(x=median_difference_CRC_minus_ESCA,y=module,color=interpretation))+geom_vline(xintercept=0,color="grey55",linewidth=.4)+geom_errorbarh(aes(xmin=CI95_low,xmax=CI95_high),height=.18,linewidth=.55)+geom_point(size=2.4)+scale_color_manual(values=c(CRC_stronger_significant="#0072B2",CRC_stronger_significant_low_gene_coverage="#E69F00",ESCA_stronger_significant="#D55E00",ESCA_stronger_significant_low_gene_coverage="#CC79A7",not_significant="#888888"))+labs(x="CRC - ESCA median module dependency z-score (95% bootstrap CI)",y=NULL,color=NULL,title="OXPHOS submodule effects")+theme_bw(base_size=8)+theme(legend.position="bottom",panel.grid.minor=element_blank())
ggsave(file.path(fig,"P2_F07_OXPHOS_submodule_effects.png"),p7,width=7.2,height=5.2,dpi=450);ggsave(file.path(fig,"P2_F07_OXPHOS_submodule_effects.pdf"),p7,width=7.2,height=5.2,device=cairo_pdf)

# F08 all evaluable CI/CIV structural/assembly genes, columns ordered by cancer type then module-score hierarchy.
hmods<-c("CI_structural","CI_assembly","CIV_structural","CIV_assembly");hd<-cg[module%in%hmods]
hwide<-dcast(hd,gene+module~sample_ID,value.var="dependency_z",fun.aggregate=median)
rn<-paste(hwide$module,hwide$gene,sep=" | ");hm<-as.matrix(hwide[,-c("gene","module")]);rownames(hm)<-rn
ann<-unique(cons[,.(sample_ID,cancer_type)]);ann<-ann[match(colnames(hm),sample_ID)];ord<-order(ann$cancer_type,colMeans(hm,na.rm=TRUE),decreasing=FALSE);hm<-hm[,ord,drop=FALSE];annd<-data.frame(Cancer=ann$cancer_type[ord],row.names=colnames(hm));ac<-list(Cancer=c(Colorectal="#0072B2",Oesophageal="#D55E00"))
png(file.path(fig,"P2_F08_CI_CIV_dependency_heatmap.png"),width=3600,height=4200,res=450,type="cairo");pheatmap(hm,cluster_rows=TRUE,cluster_cols=FALSE,show_colnames=FALSE,fontsize_row=4.5,annotation_col=annd,annotation_colors=ac,color=colorRampPalette(c("#D55E00","white","#0072B2"))(101),breaks=seq(-3,3,length.out=102),main="CI/CIV dependency (gene-wise z-score; blue = stronger)");dev.off()
cairo_pdf(file.path(fig,"P2_F08_CI_CIV_dependency_heatmap.pdf"),width=8,height=9.3);pheatmap(hm,cluster_rows=TRUE,cluster_cols=FALSE,show_colnames=FALSE,fontsize_row=4.5,annotation_col=annd,annotation_colors=ac,color=colorRampPalette(c("#D55E00","white","#0072B2"))(101),breaks=seq(-3,3,length.out=102),main="CI/CIV dependency (gene-wise z-score; blue = stronger)");dev.off()

# F09 expression scatter facets + focused genomic effect tiles.
expr_long[,cancer_group:=ifelse(cancer_type%in%c("Colorectal","Oesophageal"),cancer_type,"Other")]
p9a<-ggplot(expr_long,aes(x=expression,y=dependency_score,color=cancer_group))+geom_point(alpha=.55,size=.75)+geom_smooth(method="lm",se=FALSE,linewidth=.45)+facet_wrap(~gene,scales="free_x",ncol=4)+scale_color_manual(values=c(Colorectal="#0072B2",Oesophageal="#D55E00",Other="#777777"))+labs(x="Baseline RNA expression (log2 scale)",y="CRISPR dependency score (-LFC)",color=NULL,title="a  RNA-dependency associations")+theme_bw(base_size=7)+theme(legend.position="bottom",panel.grid.minor=element_blank())
tile<-gres[availability=="available" & biomarker!="ploidy",.(gene,biomarker,effect,BH_FDR)];tile[,sig:=ifelse(BH_FDR<.05,"*","")]
p9b<-ggplot(tile,aes(x=biomarker,y=gene,fill=effect))+geom_tile(color="white",linewidth=.2)+geom_text(aes(label=sig),size=3)+scale_fill_gradient2(low="#D55E00",mid="white",high="#0072B2",midpoint=0,limits=c(-1,1))+labs(x=NULL,y=NULL,fill="Rank-biserial",title="b  Genomic biomarker effects (blue = altered group stronger)")+theme_bw(base_size=7)+theme(axis.text.x=element_text(angle=45,hjust=1),panel.grid=element_blank())
draw9<-function(dev){dev();grid.newpage();lay<-grid.layout(2,1,heights=unit(c(1.7,1),"null"));pushViewport(viewport(layout=lay));print(p9a,vp=viewport(layout.pos.row=1));print(p9b,vp=viewport(layout.pos.row=2));popViewport();dev.off()}
draw9(function()png(file.path(fig,"P2_F09_candidate_biomarkers.png"),width=3600,height=3600,res=450,type="cairo"));draw9(function()cairo_pdf(file.path(fig,"P2_F09_candidate_biomarkers.pdf"),width=8,height=8))

cat("modules\n");print(res);cat("expression\n");print(expr_res);cat("ranking\n");print(ranktab)
