#!/usr/bin/env Rscript

# Summarize the 14 final GWASs and reproduce the Manhattan and QQ figures.
# All variants are used for inference and QQ plots. Manhattan plots retain
# every P < 1e-4 variant plus a deterministic sample of background variants.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

gwas_dir <- Sys.getenv("GWAS_DIR")
output_dir <- Sys.getenv("GWAS_REPORT_DIR", file.path(gwas_dir,"report"))
if (!nzchar(gwas_dir)) stop("Set GWAS_DIR to the merged GWAS directory.")
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)

networks <- c("VN","SMN","DAN","VAN","LN","FPN","DMN")
modalities <- c("SFC","DSFC")
phenotypes <- unlist(lapply(modalities,function(x)paste0(x,"_",networks)))
network_colours <- c(VN="#781286",SMN="#4682B4",DAN="#00760E",VAN="#C43AFA",
                     LN="#B7D67A",FPN="#E69422",DMN="#CD3E4E")
p_genomewide <- 5e-8
p_studywide <- p_genomewide/length(phenotypes)
retain_below <- 1e-4
max_background <- 250000L
lead_window_bp <- 1e6
set.seed(20261108)

safe_p <- function(x) pmax(x,.Machine$double.xmin)

find_file <- function(phenotype) {
  candidates <- file.path(gwas_dir,paste0("gwas.",phenotype,".glm.linear",c(".gz","")))
  found <- candidates[file.exists(candidates)]
  if (!length(found)) stop("Missing merged GWAS file for ",phenotype)
  found[[1L]]
}

read_gwas <- function(file) {
  x <- fread(file,showProgress=FALSE)
  chr_name <- intersect(c("#CHROM","CHROM","CHR"),names(x))[1L]
  required <- c(chr_name,"POS","ID","P","OBS_CT")
  if (anyNA(required) || !all(required %in% names(x)))
    stop("Required PLINK columns are missing from ",file)
  if (chr_name != "CHR") setnames(x,chr_name,"CHR")
  if ("TEST" %in% names(x)) x <- x[TEST=="ADD"]
  x[,`:=`(CHR=as.integer(CHR),POS=as.numeric(POS),P=as.numeric(P),
          OBS_CT=as.numeric(OBS_CT))]
  x <- x[CHR %between% c(1L,22L) & is.finite(P) & P>=0 & P<=1]
  x[P==0,P:=.Machine$double.xmin]
  for (name in c("REF","ALT","A1","OMITTED","A1_FREQ","BETA","SE","T_STAT"))
    if (!name %in% names(x)) x[,(name):=NA]
  x
}

lambda_gc <- function(p) {
  p <- p[is.finite(p) & p>0 & p<=1]
  median(qchisq(1-p,df=1),na.rm=TRUE)/qchisq(0.5,df=1)
}

select_leads <- function(hits,window_bp=lead_window_bp) {
  if (!nrow(hits)) return(hits)
  remaining <- copy(hits)[order(P)]
  leads <- list()
  while (nrow(remaining)) {
    lead <- remaining[1L]
    leads[[length(leads)+1L]] <- lead
    remaining <- remaining[!(CHR==lead$CHR & abs(POS-lead$POS)<=window_bp)]
  }
  rbindlist(leads,fill=TRUE)
}

summary_list <- hit_list <- lead_list <- manhattan_list <- qq_list <- chr_list <-
  vector("list",length(phenotypes))

for (i in seq_along(phenotypes)) {
  phenotype <- phenotypes[[i]]
  modality <- sub("_.*$","",phenotype)
  network <- sub("^[^_]+_","",phenotype)
  message("Processing ",phenotype)
  x <- read_gwas(find_file(phenotype))
  if (!nrow(x)) stop("No valid additive results for ",phenotype)
  top <- x[which.min(P)]

  summary_list[[i]] <- data.table(
    Phenotype=phenotype,Modality=modality,Network=network,
    N_variants=nrow(x),Median_N=median(x$OBS_CT,na.rm=TRUE),
    Lambda_GC=lambda_gc(x$P),Top_SNP=top$ID,Top_CHR=top$CHR,
    Top_POS=top$POS,Minimum_P=top$P,
    N_P_lt_5e_8=sum(x$P<p_genomewide),
    N_P_lt_3_57e_9=sum(x$P<p_studywide)
  )

  hits <- x[P<p_genomewide,.(Phenotype=phenotype,Modality=modality,Network=network,
    CHR,POS,ID,REF,ALT,A1,OMITTED,A1_FREQ,BETA,SE,T_STAT,OBS_CT,P)]
  hits[,Significance:=ifelse(P<p_studywide,"Study-wide","Genome-wide")]
  hit_list[[i]] <- hits
  lead_list[[i]] <- select_leads(hits)

  signal <- which(x$P<retain_below)
  background <- setdiff(seq_len(nrow(x)),signal)
  if (length(background)>max_background)
    background <- sample(background,max_background)
  manhattan_list[[i]] <- x[sort(unique(c(signal,background))),
    .(Phenotype=phenotype,Modality=modality,Network=network,CHR,POS,ID,P)]

  ordered_p <- sort(x$P)
  qq_index <- unique(round(seq(1,length(ordered_p),length.out=min(length(ordered_p),15000L))))
  qq_list[[i]] <- data.table(
    Phenotype=phenotype,Modality=modality,Network=network,
    Expected=-log10((qq_index-0.5)/length(ordered_p)),
    Observed=-log10(safe_p(ordered_p[qq_index])))
  chr_list[[i]] <- x[,.(CHR_length=max(POS,na.rm=TRUE)),by=CHR]
  rm(x); invisible(gc(FALSE))
}

summary_table <- rbindlist(summary_list)
hit_table <- rbindlist(hit_list,fill=TRUE)
lead_table <- rbindlist(lead_list,fill=TRUE)
manhattan <- rbindlist(manhattan_list)
qq <- rbindlist(qq_list)

fwrite(summary_table,file.path(output_dir,"GWAS_summary_by_phenotype.csv"))
fwrite(hit_table,file.path(output_dir,"GWAS_variants_P_lt_5e-8.csv"))
if (nrow(hit_table)) {
  fwrite(hit_table[P<p_studywide],file.path(output_dir,"GWAS_variants_P_lt_3.57e-9.csv"))
} else {
  fwrite(hit_table,file.path(output_dir,"GWAS_variants_P_lt_3.57e-9.csv"))
}
fwrite(lead_table,file.path(output_dir,"GWAS_display_lead_variants.csv"))

chr_lengths <- rbindlist(chr_list)[,.(CHR_length=max(CHR_length)),by=CHR][order(CHR)]
chr_lengths[,offset:=shift(cumsum(CHR_length),fill=0)]
chr_lengths[,midpoint:=offset+CHR_length/2]
manhattan <- merge(manhattan,chr_lengths[,.(CHR,offset)],by="CHR")
manhattan[,BP_cum:=POS+offset]
if (nrow(lead_table)) {
  lead_table <- merge(lead_table,chr_lengths[,.(CHR,offset)],by="CHR")
  lead_table[,BP_cum:=POS+offset]
}

theme_set(theme_classic(base_size=7,base_family="sans")+
  theme(axis.line=element_line(linewidth=0.3),axis.ticks=element_line(linewidth=0.3),
        strip.background=element_blank(),strip.text=element_text(face="bold"),
        panel.grid=element_blank()))

make_manhattan <- function(modality) {
  d <- manhattan[Modality==modality]
  labels <- lead_table[Modality==modality]
  ggplot(d,aes(BP_cum,-log10(safe_p(P)),colour=Network))+
    geom_point(size=0.42,alpha=0.55,stroke=0)+
    geom_hline(yintercept=-log10(p_genomewide),linetype="dashed",linewidth=0.3)+
    geom_hline(yintercept=-log10(p_studywide),linetype="solid",linewidth=0.3)+
    geom_text_repel(data=labels,aes(BP_cum,-log10(safe_p(P)),label=ID,colour=Network),
                    inherit.aes=FALSE,size=1.8,min.segment.length=0,
                    box.padding=0.15,point.padding=0.10,max.overlaps=Inf,show.legend=FALSE)+
    facet_grid(Network~.,scales="free_y")+
    scale_colour_manual(values=network_colours,drop=FALSE)+
    scale_x_continuous(breaks=chr_lengths$midpoint,labels=chr_lengths$CHR,
                       expand=expansion(mult=c(0.005,0.005)))+
    labs(title=modality,x="Chromosome",y=expression(-log[10](italic(P))))+
    theme(legend.position="none",panel.spacing.y=grid::unit(1.5,"mm"),
          strip.text.y=element_text(angle=0),axis.text.x=element_text(size=5.5))
}

p_manhattan <- make_manhattan("SFC") | make_manhattan("DSFC")

lambda_labels <- summary_table[,.(Modality,Network,
  label=sprintf("lambda[GC] == %.3f",Lambda_GC),x=0.10,y=Inf)]
p_qq <- ggplot(qq,aes(Expected,Observed,colour=Network))+
  geom_abline(intercept=0,slope=1,colour="grey35",linewidth=0.3)+
  geom_point(size=0.28,alpha=0.6,stroke=0)+
  geom_text(data=lambda_labels,aes(x=x,y=y,label=label),inherit.aes=FALSE,
            parse=TRUE,hjust=0,vjust=1.2,size=2.1)+
  facet_grid(Modality~Network)+
  scale_colour_manual(values=network_colours,drop=FALSE)+
  coord_equal()+
  labs(x=expression(Expected~~-log[10](italic(P))),
       y=expression(Observed~~-log[10](italic(P))))+
  theme(legend.position="none")

ggsave(file.path(output_dir,"GWAS_Manhattan_SFC_DSFC.pdf"),p_manhattan,
       width=183,height=215,units="mm",device=cairo_pdf)
ggsave(file.path(output_dir,"GWAS_Manhattan_SFC_DSFC.tiff"),p_manhattan,
       width=183,height=215,units="mm",dpi=600,compression="lzw")
ggsave(file.path(output_dir,"GWAS_QQ_all_phenotypes.pdf"),p_qq,
       width=183,height=62,units="mm",device=cairo_pdf)
ggsave(file.path(output_dir,"GWAS_QQ_all_phenotypes.tiff"),p_qq,
       width=183,height=62,units="mm",dpi=600,compression="lzw")

message("GWAS summaries and figures written to ",output_dir)
