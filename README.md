# bedCons

Given a BED6 file with genomic regions, and a bigWig file with some genomic signal, this program calculates mean, median, and percent of bins with greater than `minscore` signal over the given regions. 

## Install

Download the binary and make it executable: 
```
wget https://github.com/jakevc/bedCons/blob/master/bedCons
chmod +x bedCons
./bedCons -h
Usage of bedCons:
  -bedfile string
    	BED6 format file with regions to summarize wig signal over
  -bigwig string
    	big wig file
  -binsize int
    	bin size to summarize bigwig values over (default 50)
  -minscore float
    	the output column percentgt calcualtes the percentage of bigWig bins greater than or equal to this score (default 0.5)
  -outfile string
    	output filename
```


## Example

Some example bedfile with regions from the human genome: 
```
head data/hg38_n_boundaries.bed
chr1	955000	967500	hg38_n_B00012	0.624106357531	.
chr1	1307500	1340000	hg38_n_B00040	0.448953308449	.
chr1	1515000	1527500	hg38_n_B00056	0.629629903428	.
chr1	1695000	1740000	hg38_n_B00067	0.109148377678	.
chr1	1915000	1925000	hg38_n_B00084	-0.582895935256	.
chr1	2175000	2185000	hg38_n_B00109	0.132596293231	.
```

Summarize phyloP score over these regions:

```
bedCons -bigwig hg38.phyloP100way.bw -bedfile hg38_n_boundaries.bed -minscore 0.6 -outfile test.bed
```

```
head test.bed
chrom	start	end	name	score	strand	mean	median	percentgt
chr1	955000	967500	hg38_n_B00012	0	.	0.314180	-0.368830	0.232000
chr1	1307500	1340000	hg38_n_B00040	0	.	0.193881	-0.321280	0.215385
chr1	1515000	1527500	hg38_n_B00056	0	.	0.028473	-0.358220	0.128000
chr1	1695000	1740000	hg38_n_B00067	0	.	NaN	-0.130650	0.126667
chr1	1915000	1925000	hg38_n_B00084	0	.	-0.359669	-0.531190	0.130000
chr1	2175000	2185000	hg38_n_B00109	0	.	-0.040980	-0.164530	0.055000
chr1	2415000	2515000	hg38_n_B00132	0	.	-0.531748	-0.732310	0.068000
chr1	3015000	3025000	hg38_n_B00154	0	.	-0.070420	-0.334660	0.120000
chr1	3435000	3445000	hg38_n_B00190	0	.	0.113821	-0.250460	0.185000
```

