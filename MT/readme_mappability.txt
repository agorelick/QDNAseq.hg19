Steps taken to generate 50-mer mappability data for MT. Generally this follows instructions here:
- https://github.com/asntech/QDNAseq.hg38
- https://github.com/cpockrandt/genmap 

1. I downloaded chrM.fa.gz from https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chrM.fa.gz, then I extracted and renamed chrM to MT to keep convention with hg38.
> wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chrM.fa.gz
> mv chrM.fa.gz MT.fa.gz
> gunzip MT.fa.gz
> (replace chrM with MT in the header)

2. Next I used genmap to create an index for MT.fa
> genmap index -F MT.fa -I MT_index

3. Next use genmap to generate the 50mer mappability file with 2-mismatches.
> genmap map -K 50 -E 2 -I MT_index -O MT_mappability -w

4. Convert the wig file to bigwig
> wigToBigWig MT_mappability.wig MT_mappability.chrom.sizes MT_mappability.genmap.50mer.bigwig

10/19/2023
Alex Gorelick


