##
#Copyright (C) DECLERE 2015
#This program is free software; you can redistribute it and/or
#modify it under the terms of the GNU General Public License
#as published by the Free Software Foundation; either version 2
#of the License, or (at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program; if not, write to the Free Software
#Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#Compute GC curves.
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

