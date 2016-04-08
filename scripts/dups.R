#!/usr/local/bin/Rscript
suppressMessages(library(Rsamtools))
suppressMessages(library(GenomicAlignments))
source('dups_functions.R')

file=commandArgs(TRUE)[1]
out=commandArgs(TRUE)[2]
alignments=commandArgs(TRUE)[3]
readType=commandArgs(TRUE)[4]
file_all_alignments=commandArgs(TRUE)[5]
dup_threshold=as.numeric(commandArgs(TRUE)[6])
print(dup_threshold)

seqinfo<-seqinfo(BamFile(file))
what=c('rname','strand', 'pos', 'qwidth', 'qname', 'mapq')
if (readType=='pairedEnd'){
	if (alignments=='secondary'){
		flag=scanBamFlag(isUnmappedQuery=FALSE, isPaired=TRUE, isProperPair=TRUE)
	} else if (alignments=='primary'){
		flag=scanBamFlag(isUnmappedQuery=FALSE, isPaired=TRUE, isProperPair=TRUE, isSecondaryAlignment=FALSE)
	}
} else if (readType=='singleEnd'){
	if (alignments=='secondary'){
                flag=scanBamFlag(isUnmappedQuery=FALSE)
        } else if (alignments=='primary'){
                flag=scanBamFlag(isUnmappedQuery=FALSE, isSecondaryAlignment=FALSE)
        }
}

duplications=as.data.frame(matrix(0,ncol=10))
		colnames(duplications)=c('seqnames','strand','start','end_summary','dup_score','median','mapq_summary','cigar_summary','read','pair')
for (i in 1:length(seqinfo)){
	print(seqnames(seqinfo)[i])
	new=summaryFunction(seqname=seqnames(seqinfo)[i], file=file, what=what, seqinfo=seqinfo, flag=flag)
	duplications=rbind(duplications,new)
}
duplications=duplications[-1,]

n_pairs=c()
for (i in 1:length(seqinfo)){
	print(seqnames(seqinfo)[i])
	new=n_pairs_aligned(seqname=seqnames(seqinfo)[i], file=file_all_alignments, what=what, seqinfo=seqinfo, flag=flag)
	n_pairs=c(n_pairs,new)
}

n_reads_total=sum(n_pairs)*2

duplications=cbind(duplications,total_nreads=n_reads_total,prop=duplications$dup_score/n_reads_total)

write.csv(duplications, paste0(out,'_duplications.csv'), row.names=FALSE)




########## functions to be moved to independent file
