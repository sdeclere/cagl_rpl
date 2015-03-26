library(GenomicRanges)
library(GenomicAlignments)

# listing BAM files 
bams = list.files(pattern="*sorted*.bam$")

## Definition of ARS intervals 
for (bam in bams) { 

# results in file 	
sink(paste ("ARS_", bam, ".txt" ,sep=""))

aln = readGAlignments(bam)
chrs=names(genome(aln))
cvg = coverage(aln)
lib.mean = mean(mean(coverage(aln)))

for (k in chrs){
	cov.k = cvg[k][[1]]
	cov.vector = as.vector(cov.k) 
	pos.ars = which(cov.vector >= 5)
	min.ars = pos.ars[1]
	last.pos = pos.ars[1]
	count =0
	# iter on all position 
	for (pos in pos.ars){
		if (pos-last.pos>1){
			cat(bam, k, count, min.ars, last.pos, round(mean(cov.vector[min.ars:last.pos])),last.pos-min.ars,"\n", sep="\t")
			min.ars = pos
			count = count+1
		}
		last.pos=pos
	}
}
sink()

}

## ARS post-traitement
txt = list.files(pattern="*.txt$")

for (t in txt){
    ars <- read.table(t) 
    names(ars) <- c("lib", "k", "ars_id", "start", "end", "cov", "len")
    sink(paste(t, "_fused", sep=""))
    last <- ars[1,] 
    for (i in 2:nrow(ars)) { 
        cur=ars[i,]
        if (cur$k != last$k) { last=cur; next};  
        d = (cur$start - last$end)
        if (d <= 20){
            new_id = paste(last$ars_id, "+", cur$ars_id, sep="");
            #cat(cur$lib, cur$k, new_id , last$start, cur$end, ((cur$cov+last$cov)/2), (cur$len+last$len),"\n", sep="\t");
            cur = data.frame(lib=cur$lib, k=cur$k, ars_id=new_id , start=last$start, end=cur$end, cov=((cur$cov+last$cov)/2), len=(cur$len+last$len));
        } else {
            cat( as.character(last$lib), as.character(last$k), as.character(last$ars_id), last$start, last$end, last$cov,last$len,"\n", sep="\t")
        } 
        last = cur;
    }
    sink() 
}
