summaryFunction<-function(seqname, file, what, seqinfo, flag){
	bamFile=BamFile(file)
	end_coord=seqlengths(seqinfo[seqname])
	which=GRanges(seqnames=seqname, IRanges(1,end_coord))
	param<-ScanBamParam(what=what, which=which, flag=flag)
	print(paste0('Reading chromosome ', seqname, ' sample ', file))
	alignmentPairs=readGAlignmentPairs(file, param=param)
	print(paste0('Calculating dup matrix chromosome ', seqname, ' sample ', file))
	dup_mat_temp<-get_dup_matrix(alignmentPairs=alignmentPairs)
	return(dup_mat_temp)
}

get_dup_matrix<-function(alignmentPairs){
	# Possible pairs: 
	# F1R2 - forward in positive strand, reverse in negative strand
	# F2R1 - forward in negative strand, reverse in positive strand

	# Forward alignments. As explained above, some will be in + and some in -
	f_alignment=first(alignmentPairs)

	# Matrix with dups, from pairs f1r2 anf f2r1
	str<-c('+','-')
	mat<-as.data.frame(matrix(0,ncol=10))
	colnames(mat)<-c('seqnames','strand','start','end_summary','dup_score','median','mapq_summary','cigar_summary','read','pair')
	for (i in 1:length(str)){
		print(i)
		temp_mat=get_dup_pair(f_alignment,alignmentPairs=alignmentPairs,str[i])
		mat=rbind(mat,temp_mat)
	}
	mat=mat[-1,]
	return(mat)
}

get_dup_pair<-function(f_alignment,alignmentPairs,str){
	if(str=='+'){
		pair='f1r2'
	} else if(str=='-'){
		pair='f2r1'
	}

	tmp_alignmentPairs=alignmentPairs[which(as.vector(strand(f_alignment))==str),]
	
	# Forward alignments from forward in str
	tmp_f_alignment=first(tmp_alignmentPairs)
	df_first=as.data.frame(tmp_f_alignment)
	df_first=df_first[order(df_first$start),]
	dup_counts=c(table(df_first$start))
	dup_vals=rep(dup_counts,dup_counts)
	df_first=cbind(df_first,dup_first=dup_vals)
	rownames(df_first)=df_first$qname
	
	# Reverse alignments from forward in str
	tmp_r_alignment=last(tmp_alignmentPairs)
	df_last=as.data.frame(tmp_r_alignment)
	df_last=df_last[order(df_last$start),]
	dup_counts=c(table(df_last$start))
	dup_vals=rep(dup_counts,dup_counts)
	df_last=cbind(df_last,dup_last=dup_vals)
	rownames(df_last)=df_last$qname

	# Obtain duplication level in complementaries
	df=as.data.frame(cbind(qn1=df_first[rownames(df_first),]$qname, dup_first=df_first[rownames(df_first),]$dup_first, qn2=df_last[rownames(df_first),]$qname, dup_last=df_last[rownames(df_first),]$dup_last))
	df$dup_first=as.numeric(as.character(df$dup_first))
	df$dup_last=as.numeric(as.character(df$dup_last))
	rownames(df)=df$qn1
	df_first=cbind(df_first[rownames(df),], dup_last=df$dup_last, diff=abs(df$dup_first-df$dup_last))
	df_last=cbind(df_last[rownames(df),], dup_first=df$dup_first, diff=abs(df$dup_first-df$dup_last))

	# Keep dups that meet the threshold level 
	tmpPair_f=df_first[df_first$dup_first>=dup_threshold,]
	if (nrow(tmpPair_f)!=0){
		tmpPair_f=do.call('rbind',as.list(by(tmpPair_f, tmpPair_f['start'], transform, median=median(diff))))
		tmpPair_f=do.call('rbind',as.list(by(tmpPair_f, tmpPair_f['start'], transform, mapq_summary=paste(paste(names(table(mapq)),table(mapq), sep='='), collapse='; '))))
		tmpPair_f=do.call('rbind',as.list(by(tmpPair_f, tmpPair_f['start'], transform, cigar_summary=get_cigar(cigar))))
		tmpPair_f=do.call('rbind',as.list(by(tmpPair_f, tmpPair_f['start'], transform, end_summary=max(end))))
		tmpPair_f=unique(subset(tmpPair_f, select=c(seqnames, strand, start, end_summary, dup_first, median, mapq_summary, cigar_summary)))
		colnames(tmpPair_f)[which(colnames(tmpPair_f)=='dup_first')]='dup_score'
		if(pair=='f1r2'){
			tmpPair_f$read='f1'
		} else if(pair=='f2r1'){
			tmpPair_f$read='f2'
		}
		tmpPair_f$pair=pair
	}
	tmpPair_r=df_last[df_last$dup_last>=dup_threshold,]
	if (nrow(tmpPair_r)!=0){
		tmpPair_r=do.call('rbind',as.list(by(tmpPair_r, tmpPair_r['start'], transform, median=median(diff))))
		tmpPair_r=do.call('rbind',as.list(by(tmpPair_r, tmpPair_r['start'], transform, mapq_summary=paste(paste(names(table(mapq)),table(mapq), sep='='),collapse='; '))))
		tmpPair_r=do.call('rbind',as.list(by(tmpPair_r, tmpPair_r['start'], transform, cigar_summary=get_cigar(cigar))))
		tmpPair_r=do.call('rbind',as.list(by(tmpPair_r, tmpPair_r['start'], transform, end_summary=max(end))))
		tmpPair_r=unique(subset(tmpPair_r, select=c(seqnames, strand, start, end_summary, dup_last, median, mapq_summary, cigar_summary)))
		colnames(tmpPair_r)[which(colnames(tmpPair_r)=='dup_last')]='dup_score'
		if(pair=='f1r2'){
			tmpPair_r$read='r2'
		} else if(pair=='f2r1'){
			tmpPair_r$read='r1'
		}
		tmpPair_r$pair=pair
	}
	return(rbind(tmpPair_f,tmpPair_r))
}

get_cigar=function(col){
	temp_table=sort(table(col))
	if (length(temp_table)>4){
		temp_table=temp_table[(length(temp_table)-4):(length(temp_table))]
	}
	temp_table=paste(paste(names(temp_table),temp_table,sep='='),collapse='; ')
}


n_pairs_aligned<-function(seqname, file, what, seqinfo, flag){
	bamFile=BamFile(file)
	end_coord=seqlengths(seqinfo[seqname])
	which=GRanges(seqnames=seqname, IRanges(1,end_coord))
	param<-ScanBamParam(what=what, which=which, flag=flag)
	print(paste0('Reading chromosome ', seqname, ' sample ', file))
	alignmentPairs=readGAlignmentPairs(file, param=param)
	n_pairs=length(alignmentPairs)
	print(n_pairs)
	#print(paste0('Calculating dup matrix chromosome ', seqname, ' sample ', file))
	#dup_mat_temp<-get_dup_matrix(alignmentPairs=alignmentPairs)
	return(n_pairs)
}







