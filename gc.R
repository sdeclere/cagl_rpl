##
## Compute GC curves. 
##

library(seqinr)
chroms = read.fasta(file = "CAGL.faa")
chr = c("Cagl0A", "Cagl0B", "Cagl0C", "Cagl0D", "Cagl0E", "Cagl0F", "Cagl0G", "Cagl0H", "Cagl0I", "Cagl0J", "Cagl0K", "Cagl0L", "Cagl0M")

for (k in 1:13) {
	ka=chroms[[k]]
	chr.name = chr[k]
	cov = subset(T1, regexpr(chr.name, T1$chr)>0)
	starts <- seq(1, length(ka)-50000, by = 1000)
	n <- length(starts)  
	chunkGCs <- numeric(n)

	for (i in 1:n) {
        chunk <- ka[starts[i]:(starts[i]+49999)]
        chunkGC <- GC(chunk)
        chunkGCs[i] <- chunkGC
	}

	png( paste(chr.name,"_covgc.png", sep=""), width=8, height=6, units="in", res=300, bg="white")
	par(mfrow=c(2,1))
	plot(cov$pos, cov$dp, type='l', col='green', xlab="coverage")
	plot(starts, chunkGCs, type='l', col='red', xlab="GC")
	dev.off()
}

