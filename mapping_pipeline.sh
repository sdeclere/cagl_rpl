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


## Mapping pipeline.  
## bwa used is verson =  0.7.8-r455
##


# fastx doesn't handle properly compressed files
for f in *.bz2; do bunzip2 $f; done

# path to ngs executable 
export PATH=/Users/stef/pavane/bioports/ports/Darwin\@i386/bin/:$PATH

# Remove 9 first bases of each reads 
for f in *.txt; do 
    fr=`basename $f _sequence.txt`.fq.gz ; 
    echo "fastx_trimmer -f 9 -z -i ${f} -o ${fr}" ;
    fastx_trimmer -f 9 -z -i $f -o $fr
done 

# build index 
bwa index -a is CAGL_masked.faa

# gzip all fq (because i need to make room for analysis)
for f in *.fq; do gzip $f; done

# bwa aln 
# q trimming N at end end of the reads 
# I quality encoding is illumina 1.3 (look like 1.5+)
for f in *.fq; do
    fr=`basename $f .fq.gz`.sai
    echo "bwa aln -q 20 -t 2 CAGL_masked.faa $f > $fr; "
    bwa aln -q 20 -t 2 CAGL_masked.fasta $f > $fr; 
done 
  
# produce the sam file 
for f in *.sai; do 
    fq=`basename $f .sai`
    sam=`basename $f .sai`.fq.sam
	
    echo "bwa samse CAGL_masked.fasta $f $fq -f $sam"
    bwa samse -f $sam  CAGL_masked.fasta $f $fq 
done 

# produce bam file 
for f in *.sam; do
    ff=`basename ${f} .sam`
    samtools view -bq 2 -S -F 4 -T  CAGL_masked.fasta -o ${ff}_usorted.bam $f
    samtools sort ${ff}_usorted.bam ${ff}
    samtools index ${ff}.bam
    samtools rmdup ${ff}.bam `basename $ff .bam`_nodup.bam
done 

# RMDUP
for f in *.bam; do samtools rmdup -s $f `basename $f .bam`_nodup.bam ; done

# vcf 
for f in *nodup.bam; do 
    samtools mpileup -6 -uf CAGL_masked.fasta $f | bcftools view  - > `basename $f .fq.fq_nodup.bam`.vcf
done 
