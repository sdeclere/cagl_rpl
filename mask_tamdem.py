#!/usr/bin/env python

import sys
import textwrap

""" Mask a region in a fasta formated sequence"""

__author__ = "<stef>"
__date__ = "Thu Oct  4, 2012  5:58 PM"
__version__ = "1"
__email__ = "sdeclere@gmail.com"


# char to be replaced 
CHAR='g'

# location disct 
LOCS = {
    "Cagl0A":[(484269,490829),],
    "Cagl0B":[(500190,500919),],
    "Cagl0C":[(20100,22024),],
    "Cagl0D":[(588091,588386),],
    "Cagl0E":[(2670,4186),(17388,19192),(162644,163491),(671322,673549),],
    "Cagl0F":[(6193,7200),],
    "Cagl0G":[(989909,990543),],
    "Cagl0H":[(5206,5625),(254146,254547),(861098,861448),(1045107,1045432),],
    "Cagl0I":[(3588,5022),(12815,13416),(705198,707332),(962102,962826),(972648,982076),(993571,997131),(1003538,1005113),(1014891,1015638),(1021936,1024307),(1098588,1099005),(1099895,1100349),],
    "Cagl0J":[(166314,170552),(175004,175349),(490474,492446),(1166634,1168634),],
    "Cagl0K":[(15424,17414),(1299827,1300509),],
    "Cagl0L":[(1958,4642),(21502,21851),(1066365,1067898),(1085927,1087868),(1427915,1431339),(1437129,1438370),],
    "Cagl0M":[(11952,12277),(764830,765090),],
    }

def fastadb(fileobj):
    '''
    return a dict object that could iterate over a sequence of fasta formated entries.
    This will load the whole file into memory which can be a "very bad idea"(c)
    depending on the length of your sequence(s). 
    '''

    ret = {}
    lines = fileobj.readlines() # load all im memory 
    seqs = ( ''.join(lines) ).split('>')
    seqs.pop(0)# remove first null element
    
    for s in seqs:
        elines=s.split('\n')
        title = elines[0].strip().split()[0]
        seq = map(lambda x: str.strip(x), elines[1:]) 
        ret[title] = ''.join(seq)
    return ret

def replace_subseq(start, end, lst_seq):
    for i in range (start, end):
        lst_seq[i]=CHAR
    return lst_seq

def seq2fasta (seq, length=60):
    splited_seq = [seq[i:i+length] for i in range(0, len(seq), length)]
    return '\n'.join(splited_seq)

if __name__ == "__main__":
    fp = open (sys.argv[1], 'r')
    db = fastadb(fp)

    
    for key, value  in LOCS.items():
        #print 'Replace %s with %s in %s ' %  (value, CHAR, key)
        seq = db[key]
        lst_seq=list(seq)

        for slice in value:
            start = slice[0]
            end = slice[1]
            replace_subseq(start, end, lst_seq)
    
        rseq= ''.join(lst_seq)
        print '> %s' % key 
        print seq2fasta(rseq)
