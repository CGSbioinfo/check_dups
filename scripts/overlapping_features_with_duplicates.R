#!/usr/local/bin/Rscript

suppressMessages(library(reshape))
suppressMessages(library(ggplot2))
suppressMessages(library(GenomicRanges))
suppressMessages(library(rtracklayer))
suppressMessages(library(Rsamtools))
suppressMessages(library(GenomicAlignments))

indir=commandArgs(TRUE)[1]
patt=commandArgs(TRUE)[2]
gtf.file=commandArgs(TRUE)[3]
outdir=commandArgs(TRUE)[4]
bam_files_dir=commandArgs(TRUE)[5]
bam_files_patt=commandArgs(TRUE)[6]

# Load GTF file
gtf.file=import.gff(gtf.file, format='gtf', feature.type='gene')
gtf.file=gtf.file[,c(5,7,9)]
interg=gaps(gtf.file) # Get intergenic region
interg$gene_id=paste('intergenic',seq(1:length(interg)),sep='_')
interg$gene_name=paste('intergenic',seq(1:length(interg)),sep='_')
interg$gene_biotype='intergenic'
gtf.file=sort(c(gtf.file,interg))


# Load files
files=list.files(indir,pattern=patt)
for (i in 1:length(files)){
	name=gsub(patt,'',files[i])
	name=gsub('-','',name)
	assign(name, overlap_with_feature(sample=files[i],gtf.file=gtf.file, patt=patt))
}

mprop=data.frame(matrix(ncol=2))
colnames(mprop)=c('vals','sample')
for (i in 1:length(files)){
	name=gsub(patt,'',files[i])
	name=gsub('-','',name)
	vals=data.frame(cbind(vals=eval(as.name(name))$nreads_sum), sample=name)
	mprop=rbind(mprop,vals)
}
mprop=mprop[-1,]
pdf(paste0(outdir,'genome_proportion_dups.pdf'), width=12)
ggplot(data=mprop, aes(x=sample,y=vals)) + geom_boxplot()
dev.off()

# Intersect gene
name=gsub(patt,'',files[1])
name1=gsub('-','',name)
name=gsub(patt,'',files[2])
name2=gsub('-','',name)
intersect_genes=c(as.character(eval(as.name(name1))$gene),as.character(eval(as.name(name2))$gene))
for (i in 2:length(files)){
	name=gsub(patt,'',files[i])
	name=gsub('-','',name)
	intersect_genes=intersect(intersect_genes,as.character(eval(as.name(name))$gene))
}

mprop=data.frame(matrix(nrow=length(intersect_genes), ncol=length(files)))
rownames(mprop)=intersect_genes
colnames(mprop)=gsub('-','',gsub(patt,'',files))
for (i in 1:length(files)){
	name=gsub(patt,'',files[i])
	name=gsub('-','',name)
	subset_temp=eval(as.name(name))
	subset_temp=subset_temp[which(subset_temp$gene %in% intersect_genes),]
	rownames(subset_temp)=subset_temp$gene
	mprop[,name]=subset_temp[rownames(mprop),]$nreads_sum
}
mprop=cbind(gene=rownames(mprop),mprop)
mprop.melt=melt(mprop)
pdf(paste0(outdir,'intersect_dups_features.pdf'), width=12)
ggplot(data=mprop.melt, aes(x=variable,y=value, group=gene,colour=gene)) + geom_point(size=3)
dev.off()


# Plot by chromosome #
#--------------------#
# Get gene list of all included in all samples
all_genes=c()
for (i in 1:length(files)){
	name=gsub(patt,'',files[i])
	name=gsub('-','',name)
	temp_genes=as.character(eval(as.name(name))$gene)
	all_genes=c(all_genes,temp_genes)
}
all_genes=unique(all_genes)
#matrix_genes=data.frame(matrix(nrow=(length(all_genes)*length(files)), ncol=6))
matrix_genes=data.frame(matrix( ncol=6))
colnames(matrix_genes)=c('gene','chr','nreads_sum','biotype','sample','shared')
for (i in 1:length(files)){
	name=gsub(patt,'',files[i])
	name=gsub('-','',name)
	temp_sample=eval(as.name(name))
	#emp_sample=temp_sample[which(temp_sample$gene %in% all_genes),]
	temp_sample=subset(temp_sample,select=c(gene, chr, nreads_sum, biotype))
	temp_sample$sample=name
	temp_sample$shared=NA
	matrix_genes=rbind(matrix_genes,temp_sample)
}
matrix_genes=matrix_genes[-1,]

pdf(paste0(outdir,'dups_by_chr_by_sample.pdf'), width=14)
ggplot(matrix_genes, aes(x=chr, y=nreads_sum, colour=biotype)) + geom_point(shape=19, alpha=0.6) + facet_wrap(~sample, ncol=3) + theme()
dev.off()

top10_genes=head(unique(matrix_genes[order(matrix_genes$nreads_sum, decreasing=TRUE),]$gene),10)
for (i in 1:length(top10_genes)){
	temp_gene=top10_genes[i]
	matrix_genes[which(matrix_genes$gene == temp_gene),]$shared=temp_gene
}
matrix_genes[is.na(matrix_genes$shared),]$shared='unique'
submatrix_genes=matrix_genes[!matrix_genes$shared=='unique',]
ggplot(submatrix_genes, aes(x=chr, y=nreads_sum, colour=shared)) + geom_jitter(position=position_jitter(width=.2), size=3, shape=19) + facet_wrap(~sample, ncol=3)

# Refine quantification of top 10 duplications #
#----------------------------------------------#

gene_matrix=as.data.frame(matrix(nrow=length(top10_genes),ncol=4))
rownames(gene_matrix)=top10_genes
colnames(gene_matrix)=c('gene','chr','min_start','max_end')
for (i in 1:length(top10_genes)){
	start_min=c()
	for (j in 1:length(files)){
		name=gsub(patt,'',files[j])
		name=gsub('-','',name)
		temp_sample=eval(as.name(name))
		temp_start_min=temp_sample[temp_sample$gene==top10_genes[i],]$start_min
		if (!length(temp_start_min)==0){
			names(temp_start_min)=name
			start_min=c(start_min,temp_start_min)
		}
	}
	start_min_index=which(start_min==min(start_min))[1]
	start_min=start_min[start_min_index]

	end_max=c()
	for (j in 1:length(files)){
		name=gsub(patt,'',files[j])
		name=gsub('-','',name)
		temp_sample=eval(as.name(name))
		temp_end_max=temp_sample[temp_sample$gene==top10_genes[i],]$end_max
		if (!length(temp_end_max)==0){
			names(temp_end_max)=name
			end_max=c(end_max,temp_end_max)
		}
	}
	end_max_index=which(end_max==max(end_max))[1]
	end_max=end_max[end_max_index]

	gene_matrix[top10_genes[i],]$gene=top10_genes[i]
	gene_matrix[top10_genes[i],]$min_start=start_min
	gene_matrix[top10_genes[i],]$max_end=end_max
	gene_matrix[top10_genes[i],]$chr=as.character(seqnames(gtf.file[gtf.file$gene_name==top10_genes[i]]))
}

bam_files=list.files(bam_files_dir, pattern=bam_files_patt)
quantified_regions<-as.data.frame(matrix(ncol=10))
colnames(quantified_regions)<-c('gene','chr','start','end', 'strand_prop','nreads','cigar_prop','mapq_prop','qwidth_prop', 'sample')
for (i in 1:length(bam_files)){
	for(j in 1:nrow(gene_matrix)){
		info=quantify_regions(bam_files[i], gene_matrix[j,]$chr, gene_matrix[j,]$min_start, gene_matrix[j,]$max_end)
		quantified_regions=rbind(quantified_regions,c(gene=gene_matrix[j,]$gene,chr=gene_matrix[j,]$chr,info, sample=bam_files[i]))
	}
}
quantified_regions=quantified_regions[-1,]
quantified_regions$sample=gsub(bam_files_patt,'',quantified_regions$sample)
quantified_regions$nreads=as.numeric(quantified_regions$nreads)
quantified_regions$start=as.numeric(quantified_regions$start)
quantified_regions$end=as.numeric(quantified_regions$end)

gene_matrix_gene_size=as.data.frame(matrix(nrow=length(top10_genes),ncol=4))
rownames(gene_matrix_gene_size)=top10_genes
colnames(gene_matrix_gene_size)=c('gene','chr','min_start','max_end')
for (i in 1:length(top10_genes)){
	gene_matrix_gene_size[top10_genes[i],]$gene=top10_genes[i]
	gene_matrix_gene_size[top10_genes[i],]$min_start=start(gtf.file[gtf.file$gene_name==top10_genes[i]])-75
	gene_matrix_gene_size[top10_genes[i],]$max_end=end(gtf.file[gtf.file$gene_name==top10_genes[i]])+75
	gene_matrix_gene_size[top10_genes[i],]$chr=as.character(seqnames(gtf.file[gtf.file$gene_name==top10_genes[i]]))
}
gene_size=abs(gene_matrix_gene_size$max_end-gene_matrix_gene_size$min_start)
names(gene_size)=rownames(gene_matrix_gene_size)
quantified_regions$size_dup_region=abs(quantified_regions$end-quantified_regions$start)
quantified_regions$dup_region_size_prop=0
for (i in 1:length(gene_size)){
	temp_index=which(quantified_regions$gene == names(gene_size)[i])
	quantified_regions$dup_region_size_prop[temp_index]=quantified_regions$size_dup_region[temp_index]/gene_size[i]
}
pdf(paste0(outdir,'quantified_regions.pdf'), width=14)
ggplot(quantified_regions, aes(x=chr, y=as.numeric(nreads), colour=gene)) + geom_jitter(position=position_jitter(width=.2), size=3, shape=19, alpha=0.5) + facet_wrap(~sample, ncol=3)
dev.off()
pdf(paste0(outdir,'quantified_regions_metrics.pdf'), width=14)
ggplot(quantified_regions, aes(x=as.numeric(mapq_prop), y=as.numeric(cigar_prop), colour=sample, size=as.numeric(qwidth_prop), alpha=as.numeric(dup_region_size_prop))) + geom_point(shape=19) + facet_wrap(~gene, ncol=3)
dev.off()







############ functions to be moved to independent rscript
overlap_with_feature=function(sample,gtf.file,patt){
	outname=gsub(patt,'',sample)
	sample=read.csv(sample)
	#sample$chr=gsub('chr_','',sample$chr)
	#sample$F_loci=gsub('^1:','',sample$F_loci)
	#sample$F_strand=gsub('^1:','',sample$F_strand)

	sample_gr=GRanges(seqnames=sample$seqnames, 
		ranges=IRanges(start=as.numeric(sample$start), end=(as.numeric(sample$end_summary))), value=as.numeric(sample$dup_score), strand=sample$strand, median=as.numeric(sample$median), mapq_summary=sample$mapq_summary, cigar_summary=sample$cigar_summary, read=sample$read, pair=sample$pair, prop=as.numeric(sample$prop))

	overlap=findOverlaps(query=sample_gr, subject=gtf.file)
	temp=cbind(as.data.frame(sample_gr[queryHits(overlap)]), as.data.frame(gtf.file[subjectHits(overlap)]))
	colnames(temp)[1:3]=c('seqnames_dup','start_dup','end_dup')
	temp=data.frame(cbind(chr=as.character(temp$seqnames_dup), start=temp$start_dup, end=temp$end_dup, strand=as.character(temp$strand), nreads=temp$value, median=temp$median, mapq_summary=as.character(temp$mapq_summary), cigar_summary=as.character(temp$cigar_summary), read=as.character(temp$read), pair=as.character(temp$pair), gene=temp$gene_name, biotype=temp$gene_biotype, prop=temp$prop))
	write.csv(temp, paste0(outdir,'/',outname,'_dup_overlaps.csv'), row.names=FALSE)
	temp_summarized=summarize_overlap(temp)
	write.csv(temp_summarized, paste0(outdir,'/',outname,'_dup_overlaps_summary.csv'), row.names=FALSE)
	return(temp_summarized)
}

summarize_overlap=function(df){
	df$start=as.numeric(as.character(df$start))
	df$end=as.numeric(as.character(df$end))
	df$nreads=as.numeric(as.character(df$nreads))
	df$prop=as.numeric(as.character(df$prop))
	df$median_prop=as.numeric(as.character(df$median))/df$nreads
	df_temp=do.call('rbind',as.list(by(df, df['gene'], transform, start_sd=sd_norm(start))))
	df_temp=do.call('rbind',as.list(by(df_temp, df_temp['gene'], transform, end_max=max(end))))
	df_temp=do.call('rbind',as.list(by(df_temp, df_temp['gene'], transform, start_min=min(start))))
	df_temp=do.call('rbind',as.list(by(df_temp, df_temp['gene'], transform, nreads_sum=sum(nreads))))
	df_temp=do.call('rbind',as.list(by(df_temp, df_temp['gene'], transform, prop_sum=sum(prop))))
	df_temp=do.call('rbind',as.list(by(df_temp, df_temp['gene'], transform, strand_unique=paste_unique_vals(strand))))
	df_temp=do.call('rbind',as.list(by(df_temp, df_temp['gene'], transform, read_unique=paste_unique_vals(read))))
	df_temp=do.call('rbind',as.list(by(df_temp, df_temp['gene'], transform, pair_unique=paste_unique_vals(pair))))
	df_temp=do.call('rbind',as.list(by(df_temp, df_temp['gene'], transform, cigar_sum=mapqs_summary(cigar_summary))))
	df_temp=do.call('rbind',as.list(by(df_temp, df_temp['gene'], transform, mapq_sum=mapqs_summary(mapq_summary))))
	df_temp=do.call('rbind',as.list(by(df_temp, df_temp['gene'], transform, median_median_prop=median(median_prop))))
	rownames(df_temp)=1:nrow(df_temp)
	df_temp=unique(subset(df_temp,select=c(chr,start_min,end_max,start_sd,strand_unique,read_unique,pair_unique,median_median_prop,nreads_sum,prop_sum,gene,biotype, cigar_sum, mapq_sum)))
	return(df_temp)
}

sd_norm<-function(x){
	sd_normalised=sd(x)/mean(x)
	return(sd_normalised)
}

paste_unique_vals<-function(x){
	temp_unique=paste0(as.character(unique(x)),collapse=';')
	return(temp_unique)
}

mapqs_summary<-function(x){
	x=as.character(x)
	x=gsub('^ ','',unlist(strsplit(x,split=';')))
	existing_vals=strsplit(x,split='=')
	existing_vals=unique(unlist(lapply(existing_vals,'[[',1)))
	m_temp=c(rep(0,length(existing_vals)))
	names(m_temp)=paste0('m',existing_vals)
	for (i in 1:length(existing_vals)){
		x_temp=grep(paste0('^',existing_vals[i],'='),x, value=TRUE)
		x_temp=sum(as.numeric(gsub(paste0('^',existing_vals[i],'='),'',x_temp)))
		temp_index=paste0('m',existing_vals[i])
		m_temp[temp_index]=x_temp
	}
	mapq_sum=sort(gsub('^m','',paste(names(m_temp),m_temp,sep='=')))
	mapq_sum=paste(mapq_sum,collapse='; ')
	return(mapq_sum)
}


quantify_regions<-function(file, chr, start, end){
	bamFile=BamFile(paste0(bam_files_dir,'/',file))
	what=c('rname','strand', 'pos', 'qwidth', 'qname', 'mapq', 'cigar')
	which=GRanges(seqnames=chr, IRanges(start=start,end=end))
	flag=scanBamFlag(isUnmappedQuery=FALSE, isPaired=TRUE, isProperPair=TRUE, isSecondaryAlignment=FALSE)
	param<-ScanBamParam(what=what, which=which, flag=flag)
	x=scanBam(bamFile, param=param)
	lst=lapply(names(x[[1]]), function(elt){do.call(c,unname(lapply(x,'[[',elt)))})
	names(lst)=names(x[[1]])
	x=do.call('DataFrame',lst)
	x=x[x$pos>=start,]
	x=x[x$pos<=end,]
	# nreads
	nreads=nrow(x)
	# cigar summary
	subset_cigar=paste0(70:max(x$qwidth), 'M')
	tcigar=table(x$cigar)
	subset_cigar=sum(tcigar[grep(paste0(subset_cigar,collapse='|'),names(tcigar))])
	cigar_prop=subset_cigar/nreads
	# mapq
	mapq_prop=table(x$mapq)['255']/nreads
	# qwidth
	subset_qwidth=70:max(x$qwidth)
	qwidth_prop=sum(table(x$qwidth)[which(names(table(x$qwidth))%in%subset_qwidth)])/nreads
	#strand
	strand_prop=table(x$strand)[1]/sum(table(x$strand))
	info=c(start=min(x$pos), end=max(x$pos), strand_prop=strand_prop, nreads=nrow(x), cigar_prop=cigar_prop, mapq_prop=mapq_prop, qwidth_prop=qwidth_prop)
	return(info)
}





