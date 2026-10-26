# 首次运行如缺少依赖：install.packages(c("data.table", "readxl", "ggplot2"))
# 输入为全细菌基因丰度；输出表示按分类群汇总的基因相对丰度。
if (.Platform$OS.type == "windows" && !isTRUE(l10n_info()[["UTF-8"]])) {
  suppressWarnings(Sys.setlocale("LC_CTYPE", "Chinese_China.utf8"))
}
library(data.table)
setDTthreads(2L)
library(readxl)
library(ggplot2)

cli <- commandArgs(trailingOnly=TRUE)
if(length(cli)==0 || "--help" %in% cli) {
  cat("Rscript correlate_taxa.R --bacteria genes.tsv --pollutants concentrations.xlsx --output result_dir --pollutant-columns A,B [--permutations 9999 --min-n 8]\n")
  quit(status=0)
}
option <- function(name,default=NULL) {
  i<-match(name,cli)
  if(is.na(i)) {if(is.null(default)) stop("Missing option: ",name);return(default)}
  if(i==length(cli)) stop("Missing value: ",name)
  cli[i+1]
}
bacteria_file <- option("--bacteria")
pua_file <- option("--pollutants")
output_dir <- option("--output")
puas <- trimws(strsplit(option("--pollutant-columns"),",",fixed=TRUE)[[1]])
if(!length(puas)||any(puas=="")||anyDuplicated(puas)) stop("Invalid pollutant columns")
min_n <- as.integer(option("--min-n","8"))
n_perm <- as.integer(option("--permutations","9999"))
if(is.na(min_n)||min_n<3||is.na(n_perm)||n_perm<99) stop("Invalid test parameters")
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
samples <- NULL
levels <- c(Family="f__", Genus="g__", Species="s__")
set.seed(20261017)

numeric_checked <- function(x, label, allow_na=FALSE) {
  text <- trimws(as.character(x))
  missing <- is.na(x) | text == ""
  value <- suppressWarnings(as.numeric(text))
  if (any(!missing & is.na(value)) || any(!is.na(value) & !is.finite(value)) ||
      any(value < 0, na.rm=TRUE) || (!allow_na && anyNA(value)))
    stop("数值缺失、负数或格式异常：", label)
  value
}

message("读取 PUA 并核对样本……")
pua <- if(tolower(tools::file_ext(pua_file)) %in% c("xlsx","xls")) {
  as.data.table(read_excel(pua_file))
} else {
  fread(pua_file)
}
setnames(pua, trimws(gsub("[\r\n\t]", " ", names(pua))))
if (anyDuplicated(names(pua))) stop("PUA 存在重复表头")
find_column <- function(pattern) {
  hit <- grep(pattern, names(pua), ignore.case=TRUE, perl=TRUE, value=TRUE)
  if(length(hit)!=1L) stop("无法唯一识别表头：",pattern)
  hit
}
old <- c(find_column("^(Station|SampleID)(?=$|\\s|\\[|\\(|（)"), vapply(puas, function(x) {
  hit<-names(pua)[startsWith(tolower(names(pua)),tolower(x))]
  hit<-hit[vapply(hit,function(v) {tail<-substring(v,nchar(x)+1);tail==""||grepl("^[[:space:]\\[（(]",tail)},logical(1))]
  if(length(hit)!=1L) stop("Cannot identify pollutant column: ",x)
  hit
}, character(1)))
fwrite(data.table(Original=old, Analysis=c("SampleID",puas)), file.path(output_dir,"PUA_column_mapping.csv"))
setnames(pua,old,c("SampleID",puas))
pua[,SampleID:=trimws(as.character(SampleID))]
available <- names(fread(bacteria_file,nrows=0))
samples <- pua$SampleID
if(!length(samples)||!all(samples %in% available)) stop("Pollutant samples must match gene abundance columns")
if(anyNA(pua$SampleID)||any(pua$SampleID=="")||anyDuplicated(pua$SampleID)) stop("PUA 样本名缺失或重复")
if(!all(samples %in% pua$SampleID)) stop("PUA 缺少样本：",paste(setdiff(samples,pua$SampleID),collapse=","))
pua <- pua[match(samples,SampleID)]
if(!"Layer" %in% names(pua)) pua[,Layer:="Unspecified"]
for(x in puas) set(pua,j=x,value=numeric_checked(pua[[x]],x,TRUE))
fwrite(pua[,c("SampleID","Layer",puas),with=FALSE],file.path(output_dir,"PUA_analysis_input.csv"))

message("读取并检查全细菌表……")
dat <- fread(bacteria_file, select=c("Reference_ID",samples,"Taxonomy"),
             colClasses=list(character=c("Reference_ID","Taxonomy")), na.strings=c("","NA"))
if(anyNA(dat$Reference_ID)||any(trimws(dat$Reference_ID)=="")||anyDuplicated(dat$Reference_ID)) stop("Reference_ID 缺失或重复")
for(x in samples) set(dat,j=x,value=numeric_checked(dat[[x]],x))
if(anyNA(dat$Taxonomy)||any(!grepl("(^|;)\\s*[kd]__Bacteria\\s*(;|$)",dat$Taxonomy))) stop("输入包含非细菌或缺失分类记录")
totals <- vapply(dat[, samples, with=FALSE],sum,numeric(1))
if(any(totals<=0)) stop("样本总丰度为零")
fwrite(data.table(SampleID=samples,Total=totals),file.path(output_dir,"sample_totals.csv"))

# 先按完整 Taxonomy 加总，减少重复解析；分母仍是全部细菌记录。
tax <- dat[,lapply(.SD,sum),by=Taxonomy,.SDcols=samples]
rm(dat); invisible(gc())
extract_rank <- function(strings,prefix) {
  pattern <- paste0("(?:^|;)\\s*",prefix,"([^;]*)")
  hit <- regexec(pattern,strings,perl=TRUE)
  pieces <- regmatches(strings,hit)
  vapply(pieces,function(z) if(length(z)>=2) trimws(z[2]) else NA_character_,character(1))
}
prefixes <- c("k__","p__","c__","o__","f__","g__","s__")
for(prefix in prefixes) set(tax,j=prefix,value=extract_rank(tax$Taxonomy,prefix))
tax[, k__ := fifelse(is.na(k__),extract_rank(Taxonomy,"d__"),k__)]
unknown <- function(x) is.na(x)|trimws(x)==""|grepl("unclassified|unassigned|unknown|uncultured|unidentified",x,ignore.case=TRUE)

# 所有候选复用相同的随机置换索引；每次均完整保存有效样本数。
perm_cache <- new.env(parent=emptyenv())
test_pair <- function(x,y) {
  ok <- is.finite(x)&is.finite(y); x<-x[ok]; y<-y[ok]; n<-length(x)
  if(n<min_n) return(list(N=n,Rho=NA_real_,P_value=NA_real_,Status="Insufficient_N"))
  if(length(unique(x))<2L||length(unique(y))<2L)
    return(list(N=n,Rho=NA_real_,P_value=NA_real_,Status="Constant"))
  a<-rank(x); a<-a-mean(a); b<-rank(y); b<-b-mean(b)
  denominator<-sqrt(sum(a*a)*sum(b*b)); rho<-sum(a*b)/denominator
  key<-as.character(n)
  if(!exists(key,envir=perm_cache,inherits=FALSE))
    assign(key,replicate(n_perm,sample.int(n)),envir=perm_cache)
  indexes<-get(key,envir=perm_cache)
  randomized<-matrix(b[indexes],nrow=n)
  null<-as.vector(crossprod(a,randomized))/denominator
  list(N=n,Rho=rho,P_value=(1+sum(abs(null)>=abs(rho)-1e-12))/(n_perm+1),Status="OK")
}

result_list<-list(); abundance_list<-list(); coverage_list<-list()
for(level in names(levels)) {
  message("汇总并分析：",level)
  prefix<-levels[[level]]
  lineage_cols<-prefixes[seq_len(match(prefix,prefixes))]
  usable<-!unknown(tax[[prefix]])
  classified_totals<-vapply(tax[usable, samples, with=FALSE],sum,numeric(1))
  coverage_list[[level]]<-data.table(Level=level,SampleID=samples,Classified_fraction=classified_totals/totals)
  aggregated<-tax[usable,lapply(.SD,sum),by=lineage_cols,.SDcols=samples]
  if(!nrow(aggregated)) next
  aggregated[,Taxon:=get(prefix)]
  aggregated[,Taxon_key:=do.call(paste,c(.SD,sep=" | ")),.SDcols=lineage_cols]
  for(s in samples) set(aggregated,j=s,value=aggregated[[s]]/totals[[s]])
  mat<-as.matrix(aggregated[, samples, with=FALSE])
  aggregated[,`:=`(Level=level,Mean_relative_abundance=rowMeans(mat),
    Nonzero_samples=rowSums(mat>0),Prevalence=rowMeans(mat>0))]
  aggregated[,Mean_abundance_percent:=100*Mean_relative_abundance]
  # 标记 MAG / bacterium 等标签；保留在探索分析中供人工判断。
  aggregated[,Provisional_label:=grepl("\\bbacterium\\b|\\bMAG\\b|metagenome|\\bsp\\.",Taxon,ignore.case=TRUE)]
  abundance_list[[level]]<-aggregated
  selected<-aggregated[Mean_relative_abundance>=0.001]
  if(!nrow(selected)) next
  rows<-vector("list",nrow(selected)*length(puas)); k<-0L
  for(i in seq_len(nrow(selected))) {
    x<-as.numeric(unlist(selected[i, samples, with=FALSE],use.names=FALSE))
    for(p in puas) {
      k<-k+1L
      metadata<-selected[i,.(Level,Taxon_key,Taxon,Mean_relative_abundance,
        Mean_abundance_percent,Nonzero_samples,Prevalence,Provisional_label)]
      rows[[k]]<-cbind(metadata,data.table(PUA=p),as.data.table(test_pair(x,pua[[p]])))
    }
  }
  result_list[[level]]<-rbindlist(rows)
}
fwrite(rbindlist(coverage_list),file.path(output_dir,"classified_abundance_coverage.csv"))
abundances<-rbindlist(abundance_list,fill=TRUE)
fwrite(abundances,file.path(output_dir,"taxon_relative_abundance.csv"))
if(!length(result_list)) stop("没有分类群达到平均相对丰度 0.1%；丰度表已保存")
results<-rbindlist(result_list)
results[,Q_BH_level:=p.adjust(P_value,"BH"),by=Level]
results[,Q_BH_global:=p.adjust(P_value,"BH")]
results[,`:=`(Direction=fifelse(Rho>=0,"Positive","Negative"),
  Candidate_raw_p=!is.na(P_value)&abs(Rho)>=0.4&P_value<0.05,
  Candidate_FDR_level=!is.na(Q_BH_level)&abs(Rho)>=0.4&Q_BH_level<0.05)]
setorder(results,Level,P_value)
fwrite(results,file.path(output_dir,"correlation_all.csv"))
fwrite(results[Candidate_raw_p==TRUE],file.path(output_dir,"candidates_raw_p.csv"))
fwrite(results[Candidate_FDR_level==TRUE],file.path(output_dir,"candidates_FDR_level.csv"))
threshold_results<-rbindlist(lapply(c(.001,.005,.01),function(t) {
  z<-copy(results[Mean_relative_abundance>=t])
  z[,Threshold_percent:=100*t]
  z[,Q_BH_threshold_level:=p.adjust(P_value,"BH"),by=Level]
  z
}))
fwrite(threshold_results,file.path(output_dir,"correlation_by_threshold.csv"))
summary<-threshold_results[,.(Taxa_tested=sum(Status=="OK"),
  Significant_pairs_raw_p=sum(Candidate_raw_p),
  Positive_pairs=sum(Candidate_raw_p & Direction=="Positive",na.rm=TRUE),
  Negative_pairs=sum(Candidate_raw_p & Direction=="Negative",na.rm=TRUE),
  Significant_pairs_FDR=sum(abs(Rho)>=.4 & Q_BH_threshold_level<.05,na.rm=TRUE)),
  by=.(Level,Threshold_percent,PUA)]
fwrite(summary,file.path(output_dir,"screening_summary.csv"))

candidates<-results[Candidate_raw_p==TRUE]
plot_dir<-file.path(output_dir,"candidate_plots")
dir.create(plot_dir,showWarnings=FALSE)
for(i in seq_len(nrow(candidates))) {
  item<-candidates[i]
  row<-abundances[Level==item$Level & Taxon_key==item$Taxon_key]
  pd<-copy(pua)
  pd[,Abundance_percent:=100*as.numeric(unlist(row[, samples, with=FALSE],use.names=FALSE))]
  pd[,Concentration:=get(item$PUA)]
  g<-ggplot(pd,aes(Concentration,Abundance_percent,color=Layer))+
    geom_point(size=3,na.rm=TRUE)+
    geom_text(aes(label=SampleID),vjust=-.7,size=3,na.rm=TRUE)+
    scale_color_brewer(palette="Dark2",na.value="#888888")+
    labs(title=paste(item$Level,item$Taxon,sep=": "),
      subtitle=sprintf("%s | rho=%.3f | p=%.4g | q(level)=%.4g | n=%d",item$PUA,item$Rho,item$P_value,item$Q_BH_level,item$N),
      x=paste(item$PUA,"concentration (input units)"),y="Relative abundance (%)")+
    theme_classic(base_size=12)
  ggsave(file.path(plot_dir,sprintf("candidate_%04d.png",i)),g,width=8,height=5.5,dpi=200)
}
candidates[,Plot_file:=sprintf("candidate_%04d.png",seq_len(.N))]
fwrite(candidates,file.path(output_dir,"candidate_plot_index.csv"))
writeLines(c("Mean abundance includes zeros across all input samples.",
  "Denominator includes all input bacterial records, including unclassified taxa.",
  "Taxon abundance is summed gene abundance, not a direct cell count.",
  "Unrestricted permutation assumes independent, exchangeable samples.",
  "Tests assume independent exchangeable samples; not adjusted for layer or other covariates.",
  "Unknown labels excluded at target rank; provisional species labels flagged.",
  paste("Permutation count:",n_perm),"Seed: 20261017"),file.path(output_dir,"analysis_settings.txt"))
capture.output(sessionInfo(),file=file.path(output_dir,"sessionInfo.txt"))
message("完成：",output_dir)
