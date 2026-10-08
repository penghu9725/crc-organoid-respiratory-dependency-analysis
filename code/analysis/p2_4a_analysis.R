options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(survival)})
script_arg <- grep('^--file=',commandArgs(trailingOnly=FALSE),value=TRUE); script_path <- if(length(script_arg)) sub('^--file=','',script_arg[1]) else '.'
root <- normalizePath(Sys.getenv('ORGANOID_PROJECT_ROOT',unset=file.path(dirname(script_path),'../..')),mustWork=FALSE)
td <- normalizePath(Sys.getenv('ORGANOID_P2_4A_TEMP',unset=file.path(root,'data/source/.tmp_p2_4a')),mustWork=FALSE)
out <- normalizePath(Sys.getenv('ORGANOID_DERIVED_ROOT',unset=file.path(root,'data/derived/analysis_outputs')),mustWork=FALSE)
fig <- normalizePath(Sys.getenv('ORGANOID_FIGURE_ROOT',unset=file.path(root,'figures/generated')),mustWork=FALSE)
mods <- c('CI_structural','CI_assembly','CIV_structural','CIV_assembly','FeS_biogenesis')
tc <- fread(file.path(td,'tcga_scores.tsv')); tc[, module:=factor(module,levels=mods)]
rb_unpaired <- function(x,y){ w<-wilcox.test(x,y,exact=FALSE)$statistic; as.numeric(2*w/(length(x)*length(y))-1) }
tumor_normal <- rbindlist(lapply(c('COADREAD','COAD','READ'),function(pr) rbindlist(lapply(mods,function(m){
 d<-tc[module==m & (if(pr=='COADREAD') TRUE else project==pr)]; x<-d[status=='Tumor',module_score]; y<-d[status=='Normal',module_score]
 wt<-wilcox.test(x,y,exact=FALSE); data.table(comparison=pr,module=m,tumor_N=length(x),normal_N=length(y),tumor_median=median(x),normal_median=median(y),tumor_mean=mean(x),normal_mean=mean(y),median_difference_tumor_minus_normal=median(x)-median(y),effect_size_rank_biserial=rb_unpaired(x,y),raw_P=wt$p.value)
}))))
tumor_normal[,BH_FDR:=p.adjust(raw_P,'BH')]
fwrite(tumor_normal,file.path(out,'P2_4A_TCGA_MODULE_SUMMARY.tsv'),sep='\t',na='NA')

tt<-tc[status=='Tumor']; tt[,stage_major:=sub('Stage ([IV]+).*','\\1',stage)]; tt[!stage_major %chin% c('I','II','III','IV'),stage_major:=NA_character_]; tt[,stage_num:=match(stage_major,c('I','II','III','IV'))]
clin <- list(); k<-1
for(m in mods){
 d<-tt[module==m & !is.na(MSI) & MSI!='' & MSI %chin% c('MSS','MSI-L','MSI-H')]; ns<-d[,.(N=.N,median=median(module_score)),by=MSI]
 if(nrow(ns)>=2 && all(ns$N>=5)){ z<-kruskal.test(module_score~factor(MSI),d); eps<-max(0,(as.numeric(z$statistic)-length(unique(d$MSI))+1)/(nrow(d)-length(unique(d$MSI)))); clin[[k]]<-data.table(analysis='MSI_Kruskal',module=m,N=nrow(d),groups=paste(ns$MSI,ns$N,sep=':',collapse=';'),group_medians=paste(ns$MSI,sprintf('%.4f',ns$median),sep=':',collapse=';'),effect=eps,effect_label='epsilon_squared',raw_P=z$p.value,underpowered=FALSE); k<-k+1 } else { clin[[k]]<-data.table(analysis='MSI_Kruskal',module=m,N=nrow(d),groups=paste(ns$MSI,ns$N,sep=':',collapse=';'),group_medians=paste(ns$MSI,sprintf('%.4f',ns$median),sep=':',collapse=';'),effect=NA_real_,effect_label='epsilon_squared',raw_P=NA_real_,underpowered=TRUE); k<-k+1 }
 d<-tt[module==m & !is.na(stage_num)]; ns<-d[,.(N=.N,median=median(module_score)),by=stage_major][order(stage_major)]; z<-kruskal.test(module_score~factor(stage_major),d); eps<-max(0,(as.numeric(z$statistic)-length(unique(d$stage_major))+1)/(nrow(d)-length(unique(d$stage_major)))); clin[[k]]<-data.table(analysis='Stage_Kruskal',module=m,N=nrow(d),groups=paste(ns$stage_major,ns$N,sep=':',collapse=';'),group_medians=paste(ns$stage_major,sprintf('%.4f',ns$median),sep=':',collapse=';'),effect=eps,effect_label='epsilon_squared',raw_P=z$p.value,underpowered=any(ns$N<5)); k<-k+1
 sp<-cor.test(d$module_score,d$stage_num,method='spearman',exact=FALSE); clin[[k]]<-data.table(analysis='Stage_Spearman',module=m,N=nrow(d),groups=paste(ns$stage_major,ns$N,sep=':',collapse=';'),group_medians='',effect=unname(sp$estimate),effect_label='spearman_rho',raw_P=sp$p.value,underpowered=any(ns$N<5)); k<-k+1
}
clin<-rbindlist(clin,fill=TRUE); clin[,BH_FDR:=p.adjust(raw_P,'BH'),by=analysis]; fwrite(clin,file.path(out,'P2_4A_TCGA_MSI_STAGE.tsv'),sep='\t',na='NA')

coxrows<-list(); k<-1
for(m in mods){
 d<-copy(tt[module==m & !is.na(OS.time) & !is.na(OS)]); d[,score_z:=as.numeric(scale(module_score))]
 for(model in c('univariable','adjusted_age_stage_MSI')){
  if(model=='univariable') dd<-d[complete.cases(OS.time,OS,score_z)] else dd<-d[complete.cases(OS.time,OS,score_z,age,stage_major,MSI) & MSI %chin% c('MSS','MSI-L','MSI-H')]
  form<-if(model=='univariable') Surv(OS.time,OS)~score_z else Surv(OS.time,OS)~score_z+age+factor(stage_major)+factor(MSI)
  fit<-coxph(form,dd); sm<-summary(fit); cf<-sm$coefficients['score_z',]; ci<-sm$conf.int['score_z',]; php<-tryCatch(cox.zph(fit)$table['score_z','p'],error=function(e) NA_real_)
  coxrows[[k]]<-data.table(model=model,module=m,N=nrow(dd),events=sum(dd$OS==1),HR_per_1SD=ci['exp(coef)'],CI95_low=ci['lower .95'],CI95_high=ci['upper .95'],coefficient=cf['coef'],raw_P=cf['Pr(>|z|)'],PH_assumption_P=php); k<-k+1
 }
}
coxres<-rbindlist(coxrows); coxres[,BH_FDR:=p.adjust(raw_P,'BH'),by=model]; fwrite(coxres,file.path(out,'P2_4A_TCGA_SURVIVAL.tsv'),sep='\t',na='NA')
cov<-fread(file.path(td,'gene_coverage.tsv')); cov[,coverage_fraction:=detected_genes/requested_genes]; fwrite(cov[grepl('TCGA',dataset)],file.path(out,'P2_4A_TCGA_GENE_COVERAGE.tsv'),sep='\t',na='NA'); fwrite(cov[dataset=='GSE132465'],file.path(out,'P2_4A_SC_GENE_COVERAGE.tsv'),sep='\t',na='NA')

p10<-ggplot(tc,aes(status,module_score,fill=status))+geom_violin(scale='width',trim=TRUE,alpha=.7)+geom_boxplot(width=.15,outlier.shape=NA)+facet_grid(project~module,scales='free_y')+scale_fill_manual(values=c(Normal='#56B4E9',Tumor='#D55E00'))+labs(x=NULL,y='Within-sample rank module score',title='TCGA COAD/READ module expression programs',subtitle='Expression activity; not CRISPR dependency')+theme_bw(base_size=10)+theme(legend.position='none',axis.text.x=element_text(angle=25,hjust=1),strip.text.x=element_text(size=8))
ggsave(file.path(fig,'P2_F10_TCGA_module_scores.png'),p10,width=12,height=6,dpi=300); ggsave(file.path(fig,'P2_F10_TCGA_module_scores.pdf'),p10,width=12,height=6)
pe<-rbind(clin[analysis %chin% c('MSI_Kruskal','Stage_Kruskal'),.(family=analysis,module,effect,lo=NA_real_,hi=NA_real_,FDR=BH_FDR)],coxres[,.(family=paste0('Cox_',model),module,effect=log(HR_per_1SD),lo=log(CI95_low),hi=log(CI95_high),FDR=BH_FDR)])
p11<-ggplot(pe,aes(module,effect,color=FDR<.05))+geom_hline(yintercept=0,lty=2,color='grey50')+geom_point(size=2)+geom_errorbar(aes(ymin=lo,ymax=hi),width=.15,na.rm=TRUE)+facet_wrap(~family,scales='free_y',ncol=2)+coord_flip()+scale_color_manual(values=c('FALSE'='grey40','TRUE'='#D55E00'))+labs(x=NULL,y='Effect (epsilon² or log HR)',title='TCGA clinical associations',color='FDR < 0.05')+theme_bw(base_size=10)
ggsave(file.path(fig,'P2_F11_TCGA_clinical_associations.png'),p11,width=10,height=8,dpi=300); ggsave(file.path(fig,'P2_F11_TCGA_clinical_associations.pdf'),p11,width=10,height=8)

sc<-fread(file.path(td,'sc_scores.tsv')); sc[,module:=factor(module,levels=mods)]
ct<-sc[level=='patient_class_celltype']; cts<-ct[,.(patient_N=uniqueN(Patient),total_cells=sum(n_cells),median_program_score=median(module_score),mean_program_score=mean(module_score),IQR_program_score=IQR(module_score)),by=.(Class,Cell_type,module)]; fwrite(cts,file.path(out,'P2_4A_SC_CELLTYPE_PROGRAM.tsv'),sep='\t',na='NA')
epi<-ct[Cell_type=='Epithelial cells']; wide<-dcast(epi,Patient+module~Class,value.var='module_score'); nwide<-dcast(epi,Patient+module~Class,value.var='n_cells'); setnames(nwide,c('Normal','Tumor'),c('Normal_cells','Tumor_cells')); wide<-merge(wide,nwide,by=c('Patient','module'),all=TRUE)
tests<-rbindlist(lapply(mods,function(m){d<-wide[module==m & !is.na(Normal)&!is.na(Tumor)]; dif<-d$Tumor-d$Normal; wt<-wilcox.test(d$Tumor,d$Normal,paired=TRUE,exact=FALSE); rr<-rank(abs(dif)); rb<-if(sum(rr)>0) sum(rr*sign(dif))/sum(rr) else 0; data.table(row_type='paired_test',Patient=NA_character_,Class='Tumor_minus_Normal',module=m,n_cells=NA_integer_,module_score=median(dif),paired_N=nrow(d),effect_size_matched_rank_biserial=rb,raw_P=wt$p.value)})); tests[,BH_FDR:=p.adjust(raw_P,'BH')]
ps<-epi[,.(row_type='patient_score',Patient,Class,module,n_cells,module_score,paired_N=NA_integer_,effect_size_matched_rank_biserial=NA_real_,raw_P=NA_real_,BH_FDR=NA_real_)]; fwrite(rbindlist(list(ps,tests),fill=TRUE),file.path(out,'P2_4A_SC_PSEUDOBULK.tsv'),sep='\t',na='NA')

p12<-ggplot(cts,aes(module,paste(Class,Cell_type,sep=': '),fill=median_program_score))+geom_tile(color='white')+scale_fill_gradient2(low='#2166AC',mid='white',high='#B2182B')+labs(x=NULL,y=NULL,fill='Median score',title='GSE132465 patient-level cell-type programs',subtitle='Descriptive; groups require ≥100 cells per patient/class/cell type')+theme_minimal(base_size=10)+theme(axis.text.x=element_text(angle=35,hjust=1))
ggsave(file.path(fig,'P2_F12_SC_celltype_program.png'),p12,width=10,height=6,dpi=300); ggsave(file.path(fig,'P2_F12_SC_celltype_program.pdf'),p12,width=10,height=6)
pairedp<-epi[Patient %chin% wide[!is.na(Normal)&!is.na(Tumor),unique(Patient)]]
p13<-ggplot(pairedp,aes(Class,module_score,group=Patient,color=Patient))+geom_line(alpha=.55)+geom_point(size=2)+facet_wrap(~module,scales='free_y')+labs(x=NULL,y='Patient pseudobulk module score',title='Paired malignant vs normal epithelial programs',subtitle='Patient is the inferential unit')+theme_bw(base_size=10)+theme(legend.position='none')
ggsave(file.path(fig,'P2_F13_SC_malignant_vs_normal.png'),p13,width=11,height=6,dpi=300); ggsave(file.path(fig,'P2_F13_SC_malignant_vs_normal.pdf'),p13,width=11,height=6)

org<-fread(file.path(out,'P2_3B_OXPHOS_SUBMODULE_RESULTS.tsv')); org<-org[module %chin% mods]
integ<-merge(data.table(module=mods),org[,.(module,organoid_CRISPR_evidence=paste0(interpretation,'; CRC-ESCA effect=',signif(median_difference_CRC_minus_ESCA,3),'; FDR=',signif(BH_FDR,3)))],by='module',all.x=TRUE)
integ<-merge(integ,tumor_normal[comparison=='COADREAD',.(module,TCGA_tumor_normal_direction=ifelse(median_difference_tumor_minus_normal>0,'tumor_higher','tumor_lower'),TCGA_tumor_normal_effect=effect_size_rank_biserial,TCGA_tumor_normal_FDR=BH_FDR)],by='module',all.x=TRUE)
clbest<-clin[order(BH_FDR),.SD[1],by=module]; svbest<-coxres[order(BH_FDR),.SD[1],by=module]; integ<-merge(integ,clbest[,.(module,TCGA_clinical_association=paste0(analysis,' effect=',signif(effect,3),'; FDR=',signif(BH_FDR,3)))],by='module',all.x=TRUE); integ<-merge(integ,svbest[,.(module,TCGA_survival_association=paste0(model,' HR=',signif(HR_per_1SD,3),'; FDR=',signif(BH_FDR,3)))],by='module',all.x=TRUE)
mal<-cts[Class=='Tumor'&Cell_type=='Epithelial cells',.(module,malignant_epithelial_median=median_program_score)]; oth<-cts[!(Class=='Tumor'&Cell_type=='Epithelial cells'),.(other_median=median(median_program_score)),by=module]; integ<-merge(integ,merge(mal,oth,by='module'),by='module',all.x=TRUE); integ[,single_cell_malignant_localization:=ifelse(malignant_epithelial_median>other_median,'higher_than_cross_celltype_median','not_higher_than_cross_celltype_median')]
integ<-merge(integ,tests[,.(module,SC_paired_N=paired_N,SC_malignant_vs_normal_effect=module_score,SC_malignant_vs_normal_FDR=BH_FDR)],by='module',all.x=TRUE)
integ[,agreement_disagreement:=ifelse((TCGA_tumor_normal_effect>0)==(SC_malignant_vs_normal_effect>0),'TCGA_and_SC_direction_agree','TCGA_and_SC_direction_disagree')]; integ[,QC_limitation:='Expression programs do not measure dependency; TCGA legacy expression has incomplete assembly-factor coverage; SC paired epithelial analysis has six patients and uses published tumor epithelial annotation.']
fwrite(integ,file.path(out,'P2_4A_EXTERNAL_VALIDATION_SUMMARY.tsv'),sep='\t',na='NA')
cat('DONE\n')
