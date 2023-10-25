# This script loads pre-generated QDNAseq rda objects (with hg38 bin annotations for various binsizes)
# and then generates and adds the appropriate bins for MT/chrM.
#
# The pre-generated QDNAseq objects are located in MT/rda_orig
# The updated versions with MT/chrM data added are saved to data/
#
# Be sure to download binary version of bigWigAverageOverBed at https://hgdownload.soe.ucsc.edu/admin/exe/ and add it to your $PATH

rm(list=ls())
library(Biobase)
library(BSgenome.Hsapiens.UCSC.hg19)
library(QDNAseq)
library(future)
library(data.table)

# set working directory
setwd('~/lab_repos/QDNAseq.hg19')

## load chrM sequence
mt <- read.table('MT/MT.fa',sep='\n',header=T)[[1]]
mt <- paste(mt, collapse='')
mt_chr <- strsplit(mt,'')[[1]]
mt_len <- length(mt_chr)
mt_dat <- data.table(nt=mt_chr)
mt_dat$pos <- 1:nrow(mt_dat)
mt_dat[,pos2:=pos+1]
mt_dat$chr <- 'chrM'
setkey(mt_dat,'chr','pos','pos2')


for (binsize in c(1000, 500, 100, 50, 30, 15, 10, 5, 1)) {
    message(binsize)

    bin_chr <- function(binsize, dat, len, label) {
        starts <- seq(1, len, by=binsize*1000)
        ends <- c(tail(starts, -1) - 1, len)
        bins <- data.table(chr=label, start=starts, end=ends)
        bins[,region:=paste0(chr,':',start,'-',end)]
        setkey(bins,'chr','start','end')
        dat_binned <- foverlaps(dat, bins, type='any')
        dat_binned <- dat_binned[!duplicated(pos),]
        dat_binned
    }
    mt_dat_binned <- bin_chr(binsize, mt_dat, mt_len, 'chrM')

    summarize_bin <- function(mt_dat_binned) { 
        len <- nrow(mt_dat_binned)
        nonN_bases <- sum(toupper(mt_dat_binned$nt) %in% c('A','C','G','T'))
        bases <- round((nonN_bases / len)*100, 4)
        n_gc <- sum(toupper(mt_dat_binned$nt) %in% c('C','G'))
        gc <- round((n_gc / nonN_bases)*100,5)
        data.table(bases=bases, gc=gc, blacklist=0, residual=NA, use=F, len=len)
    }
    mt_res <- mt_dat_binned[,summarize_bin(.SD), by=c('chr','start','end','region')]
    setnames(mt_res,'chr','chromosome')

    add_mappability <- function(res, map_file) {
        res[,region:=paste0(chromosome,':',start,'-',end)] 
        res <- as.data.frame(res)
        rownames(res) <- res$region
        res$region <- NULL

        ## calculate average mappability for MT
        mappability <- calculateMappability(res, bigWigFile=map_file, bigWigAverageOverBed="bigWigAverageOverBed")
        res$mappability <- mappability

        res
    }
    ## duplicate chrM for compatibility with GRCh38 and hg38
    res <- add_mappability(mt_res, 'MT/chrM_map.bw')
    res$chromosome <- 'MT'
    rownames(res) <- gsub('chrM','MT',rownames(res))

    ## load the rdata object for the given binsize
    rda_file <- paste0('MT/rda_orig/hg19.',binsize,'kbp.SR50.rda')
    load(rda_file)
    obj <- eval(parse(text=paste0('hg19.',binsize,'kbp.SR50')))
    bins <- obj@data
    res <- res[,names(bins)]
    bins <- bins[!bins$chromosome %in% c('Y','MT'),]
    bins <- rbind(bins, res)
    bins$residual[is.nan(bins$residual)] <- NA
    bins$gc[is.nan(bins$gc)] <- NA


    bins <- AnnotatedDataFrame(bins,
        varMetadata=data.frame(labelDescription=c(
        "Chromosome name",
        "Base pair start position",
        "Base pair end position",
        "Percentage of non-N nucleotides (of full bin size)",
        "Percentage of C and G nucleotides (of non-N nucleotides)",
        "Average mappability of 50mers with a maximum of 2 mismatches",
        "Percent overlap with ENCODE blacklisted regions",
        "Median loess residual from 1000 Genomes (50mers)",
        "Whether the bin should be used in subsequent analysis steps"),
        row.names=colnames(bins)))

    QDNAseqInfo <- list(
                        author="Modified by Alex Gorelick",
                        date=Sys.time(),
                        organism='Hsapiens',
                        build='hg19',
                        version=packageVersion("QDNAseq"),
                        url=paste0("https://github.com/agorelick/QDNAseq.hg19/tree/main/data/hg19.",binsize,"kbp.SR50.rda"),
                        md5=digest::digest(bins@data),
                        sessionInfo=sessionInfo())

    attr(bins, "QDNAseq") <- QDNAseqInfo
   
    ## hacky solution to save the rda object with the expected name format e.g. hg38.1000kbp.SR50
    eval(parse(text=paste0('hg19.',binsize,'kbp.SR50 <- bins'))) 
    save(list=paste0('hg19.',binsize,'kbp.SR50'), file=paste0("data/hg19.", binsize, "kbp.SR50.rda"), compress='xz')
}


#rm(list=ls())
#load('data/hg38.1000kbp.SR50.rda')
#qc <- as.data.table(hg38.1000kbp.SR50@data)
#qc[chromosome=='Y']
#qc[start <= 2000001 & chromosome=='Y',]
#qc[chromosome=='MT',]

## after installing the package, try these commands. rowNames should include MT bins
#library(QDNAseq)
#library(QDNAseq.hg38)
#bins <- getBinAnnotations(binSize=50, genome="hg19")
#qc <- as.data.table(bins@data)


