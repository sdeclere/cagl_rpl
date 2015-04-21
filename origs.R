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

###
### R code related to origins detection
###

# -----------------------------------------------------------------------------
require(zoo)
require(ggplot2)
require(parallel)


# -----------------------------------------------------------------------------

# extract pos & depth from vcf files
read.vcf <- function (filename) {
	unix.cmd = paste ("cat", filename, "| awk '{split($8, aa, \";\"); split(aa[1], bb,\"=\"); print $1, $2, bb[2] + 0 }' ");
	read.table( pipe(unix.cmd),  comment.char = "#", col.names=c("chr", "pos", "dp") )
}

# median filter :
# filter out too big or too tiny values of coverage
median.filter <- function (t) {
	m = median(t$dp)
	t=subset(t, t$dp > m/2)
	subset(t, t$dp < m*2)
}

# 4P logistic regression
# F(x) = ((A-D)/(1+((x/C)^B))) + D
exp.fit <- function (T) {
	time  = c(0, 55, 60, 65, 70, 75, 80)
	dp = as.numeric(T[2:8])
	a = min(dp)
	d = max(dp)
	dat = data.frame(y=dp, x=time)
	f = function(x,b,c) { ((a-d) / (1+ ((x/c)^b)))+d } 
	fm <- nls(y ~ f(x,b,c), data = dat, start=c(b=10, c=70), control = list(maxiter = 500)) 
	co <- coef(fm) 
	co
}

# used by mclapply to // computation
# return NA in case of failure 
loopfun <- function(i) {
	 return(tryCatch(exp.fit(ka[i,]), error=function(e) NA))
}

# find peaks by comparing a local maximum filter to the smooth
# http://stats.stackexchange.com/questions/36309/how-do-i-find-peaks-in-a-dataset
my.peaks <- function(y,w=50) {
	n <- length(y)
	y.max <- rollapply(zoo(y), 2*w+1, max, align="center")
	delta <- y.max - y[-c(1:w, n+1-1:w)]
	i.max <- which(delta <= 0) + w
	i.max
}

# find valleys by comparing a local maximum filter to the smooth
my.valleys <- function(y,w=50) {
	n <- length(y)
	y.min <- rollapply(zoo(y), 2*w+1, min, align="center")
	delta <- y.min - y[-c(1:w, n+1-1:w)]
	i.min <- which(delta >= 0) + w
	i.min
}

## Same function but using "pastecs" package instead of the local maxima method 
## -> this version was used to produced the results presented in the paper 

# new version of my.peaks/valleys
require(pastecs)
	my.peaks <- function(y) {
	tp=turnpoints(y)
	tp$peaks
}

my.valleys <- function(y) {
	tp=turnpoints(y)
	tp$pit
}


# -----------------------------------------------------------------------------
# Analysis Workflow 
#

# load CVF file (pileups) 
T0 <- read.vcf("s_1_GF-t0_nodup.vcf")
T1 <- read.vcf("s_2_GF-t1_nodup.vcf")
T2 <- read.vcf("s_3_GF-t2_nodup.vcf")
T3 <- read.vcf("s_4_GF-t3_nodup.vcf")
T4 <- read.vcf("s_5_GF-t4_nodup.vcf")
T5 <- read.vcf("s_6_GF-t5_nodup.vcf")
T6 <- read.vcf("s_7_GF-t6_nodup.vcf")

# filter artefacs 
T0 = median.filter(T0)
T1 = median.filter(T1)
T2 = median.filter(T2)
T3 = median.filter(T3)
T4 = median.filter(T4)
T5 = median.filter(T5)
T6 = median.filter(T6)

# DNA Q correction (given by FACS) 
T0$dp.Q <- T0$dp * 1.15
T1$dp.Q <- T1$dp * 1.22
T2$dp.Q <- T2$dp * 1.37
T3$dp.Q <- T3$dp * 1.51
T4$dp.Q <- T4$dp * 1.65
T5$dp.Q <- T5$dp * 1.74
T6$dp.Q <- T6$dp * 1.75

# Q of reads correction 
ref.cov.med = median(T0$dp.Q)
T0$dp.C <- T0$dp.Q * ( (median(T0$dp.Q) / ref.cov.med ) )
T1$dp.C <- T1$dp.Q * ( (median(T1$dp.Q) / ref.cov.med ) )
T2$dp.C <- T2$dp.Q * ( (median(T2$dp.Q) / ref.cov.med ) )
T3$dp.C <- T3$dp.Q * ( (median(T3$dp.Q) / ref.cov.med ) )
T4$dp.C <- T4$dp.Q * ( (median(T4$dp.Q) / ref.cov.med ) )
T5$dp.C <- T5$dp.Q * ( (median(T5$dp.Q) / ref.cov.med ) )
T6$dp.C <- T6$dp.Q * ( (median(T6$dp.Q) / ref.cov.med ) )

# set dataset annotations 
comment(T0) <- "T0"
comment(T1) <- "T1"
comment(T2) <- "T2"
comment(T3) <- "T3"
comment(T4) <- "T4"
comment(T5) <- "T5"
comment(T6) <- "T6"

# list of chr ids 
chr = c("Cagl0A", "Cagl0B", "Cagl0C", "Cagl0D", "Cagl0E", "Cagl0F", "Cagl0G", "Cagl0H", "Cagl0I", "Cagl0J", "Cagl0K", "Cagl0L", "Cagl0M")

T50 <- new.env()
for  (k in chr) {
	# for each chr extract coverage by time 
	kt0 = subset(T0, regexpr(k, T0$chr)>0, select=c(pos,dp.C))
	names(kt0) <- c("pos", "0")
	kt1 = subset(T1, regexpr(k, T1$chr)>0, select=c(pos,dp.C))
	names(kt1) <- c("pos", "1")
	kt2 = subset(T2, regexpr(k, T2$chr)>0, select=c(pos,dp.C))
	names(kt2) <- c("pos", "2")
	kt3 = subset(T3, regexpr(k, T3$chr)>0, select=c(pos,dp.C))
	names(kt3) <- c("pos", "3")
	kt4 = subset(T4, regexpr(k, T4$chr)>0, select=c(pos,dp.C))
	names(kt4) <- c( "pos", "4")
	kt5 = subset(T5, regexpr(k, T5$chr)>0, select=c(pos,dp.C))
	names(kt5) <- c("pos", "5")
	kt6 = subset(T6, regexpr(k, T6$chr)>0, select=c(pos,dp.C))
	names(kt6) <- c("pos", "6")
	# merged data 
	ka=Reduce(function(x,y) {merge(x,y, by="pos")}, list(kt0, kt1,kt2,kt3,kt4,kt5,kt6) )
	# do compute the T50
	co <- mclapply(1:nrow(ka),loopfun, mc.cores=6L)
	result = data.frame(do.call("rbind", co))
	assign (k, data.frame(pos=ka$pos, t50=result$c, s=result$b), envir=T50)
}


# -----------------------------------------------------------------------------
## -- Search for Peaks & Valleys 

LOESS <- new.env() 
MODEL <- new.env()

# smoothing T50 curves 
for  (k in chr) {
	
	# loess on T50
	df = get(k, env=T50)
	x = df$"pos"
	y = df$"t50"
	y.loess <- loess(y ~ x, span=0.04, data.frame(x=x, y=y))
	y.predict <- predict(y.loess, data.frame(x=x))
	assign (k, data.frame(pos=x, predict=y.predict), envir=LOESS)
	assign (k, y.loess, envir=MODEL)
}

TERM <- new.env() 

# peaks calling  
for (k in chr) {
	df = get(k, env=LOESS)
	p=my.peaks(as.numeric(df$predict))
	pos = df$pos[p]
	peaks = df$predict[p]
	assign (k, data.frame(x.pos=pos, y.peaks=peaks), envir=TERM)
}

ORIG <- new.env()

# valleys calling  
for (k in chr) {
	df = get(k, env=LOESS)
	p=my.valleys(as.numeric(df$predict))
	pos = df$pos[p]
	peaks = df$predict[p]
	assign (k, data.frame(x.pos=pos, y.valleys=peaks), envir=ORIG)
}

# save results in file 
for (k in chr) {
	term = get(k, env=TERM)
	orig = get(k, env=ORIG)
	write.csv(orig, file= paste(k,"_orig.csv", sep=""), row.names = F)
	write.csv(term, file= paste(k,"_term.csv", sep=""), row.names = F)
}

# save results in file 
df=data.frame()
for (k in chr) {
	orig = get(k, env=ORIG)
	chr.name = rep(k,length(orig$x.pos))
	head(chr.name)
	df = rbind(df, data.frame(chr=chr.name, pos=orig$x.pos, time=orig$y.valleys))
}
write.csv(df, file="orig.csv", sep="", row.names = F, append=TRUE)


# save results in file 
df=data.frame()
for (k in chr) {
	orig = get(k, env=TERM)
	chr.name = rep(k,length(orig$x.pos))
	head(chr.name)
	df = rbind(df, data.frame(chr=chr.name, pos=orig$x.pos, time=orig$y.peaks))
}
write.csv(df, file="term.csv", row.names = F)

## -- END 
