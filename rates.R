##
#Copyright (C) DECLERE 2015
#
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
#

# This script contains instruction used to compute migration forks rates 
# see article for more explanation. 

require (zoo)
require ( ggplot2 )

# new operator defined for conveniance 
`%between%`<-function(x,rng) x>rng[1] & x<rng[2]

ret0 = matrix(ncol=50)
ret4 = matrix(ncol=50)

for (k in chr) {
	# this contains bonafides origins 
    orig.all = get(k, env=BORIG)
	
    # early / lates filter
    orig = orig.all[(orig.all$y.valleys  > 65.5) && (orig.all$y.valleys  < 68.5),]

	for (i in 1:nrow(orig)) {
        o = orig[i, 1]
		# get raw values bw o -> o +2kb 
        ds0 = T0[T0$pos %between% c(o, o+20000) & T0$chr==k,]$dp.C
		ds4 = T4[T4$pos %between% c(o, o+20000) & T4$chr==k,]$dp.C
        
		# apply a "floating" mean on those values 
        slope0 = rollapply(ds0, FUN=mean,width=1000)
		slope4 = rollapply(ds4, FUN=mean,width=1000)
		
		# extract a sample a points 
        sample0 = seq(1, length(slope0), length(slope0)/50)
		sample4 = seq(1, length(slope4), length(slope4)/50)
        x0 = sapply(sample0, function(x){slope0[x]})
		x4 = sapply(sample4, function(x){slope4[x]})
		
		# accumulate res 
        ret0 = rbind(ret0, x0)
		ret4 = rbind(ret4, x4)
	}
}

#compute the regresssion 
rc = matrix(ncol=7)
for (c in seq(1,50)) {
    ry = c(mean(ret4[,c]/ret0[,c], na.rm=T), mean(ret5[,c]/ret0[,c], na.rm=T))
    rx = c(1,2)
    rxy = data.frame(x=rx, y=ry)  
    fit <- lm(y ~ x, rxy)

    ry.min = c(min(ret4[,c]/ret0[,c], na.rm=T), min(ret5[,c]/ret0[,c], na.rm=T))  
    ry.max = c(max(ret4[,c]/ret0[,c], na.rm=T), max(ret5[,c]/ret0[,c], na.rm=T))
 
    rxy.min = data.frame(x=rx, y=ry.min)
    rxy.max = data.frame(x=rx, y=ry.max)  
    fit.min <- lm(y ~ x, rxy.min)
    fit.max <- lm(y ~ x, rxy.max)

    toto = cbind(c, fit$coefficients[1], fit$coefficients[2], fit.min$coefficients[1], fit.min$coefficients[2], fit.max$coefficients[1], fit.max$coefficients[2])
    rc = rbind(rc, toto)
}


## plots 
pdf()

line1=data.frame(x=1:15, y=rc[,3][2:16])
model.line1 <- lm(y ~ x, line1)

line3=data.frame(x=1:5, y=rc[,3][2:6])
model.line3 <- lm(y ~ x, line3)

line2=data.frame(x=c(which.min(rc[,3]), which.min(rc[,3])-1), y=c(min(rc[,3], na.rm=T), min(rc[,3], na.rm=T)))
model.line2 <- lm(y ~ x, line2)

p=qplot(1:51, rc[,3])+geom_point(shape=1)
p= p + geom_abline(intercept=coef(model.line1)[1],slope=coef(model.line1)[2],colour="red")
p= p + geom_abline(intercept=coef(model.line2)[1],slope=coef(model.line2)[2],colour="blue")
p= p + geom_abline(intercept=coef(model.line3)[1],slope=coef(model.line3)[2],colour="orange")

print(p)
dev.off()

# used to search intersect between the two lines 
c(-solve(cbind(cm[,2],-1)) %*% cm[,1])



