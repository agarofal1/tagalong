suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(Ckmeans.1d.dp))
suppressMessages(library(gtools))
suppressMessages(library(pbmcapply))

args <- commandArgs(T)

samp.id <- args[1]
left.clip.path <- args[2]
right.clip.path <- args[3]
clip.length.cutoff <- as.numeric(args[4])
map.q.cutoff <- as.numeric(args[5])
cluster.supp.cutoff <- as.numeric(args[6])
outdir <- args[7]

get.clusters <- function(reads){
  chrs <- unique(reads$chr)
  cluster.chr <- function(chrom){
    chr.sub <- subset(reads, chr == chrom)
    if(length(unique(chr.sub$break.pos)) != 1){
      clusters <- paste(chrom, Cksegs.1d.dp(chr.sub$break.pos, c(2, nrow(chr.sub)))$cluster, sep="_")
    } else{
      clusters <- paste(chrom, rep(1, nrow(chr.sub)), sep="_")
    }
    return(clusters)
  }
  
  clusters <- unlist(pbmclapply(chrs, function(x) cluster.chr(x), mc.cores=length(chrs)))
}

mixedrank = function(x) order(gtools::mixedorder(x))

us.clips <- fread(right.clip.path, header=F, sep="\t", stringsAsFactors=F) %>% dplyr::rename(name=1, flag=2, chr=3, start.pos=4, map.q=5, cigar=6, seq=7) %>% group_by(name, flag) %>% mutate(read.length=nchar(seq)) %>% ungroup()
us.clips$clip.length <- as.numeric(sapply(strsplit(us.clips$cigar, "[A-Z]"), tail, 1))
us.clips <- us.clips %>% mutate(map.length = read.length - clip.length, break.pos = start.pos + map.length - 1, chr.numeric = case_when(grepl("X", chr)~23, grepl("Y", chr)~24, TRUE~as.numeric(gsub("chr", "", chr)))) %>% filter(clip.length >= clip.length.cutoff & map.q >= map.q.cutoff) %>% arrange(chr.numeric, break.pos) %>% dplyr::select(-chr.numeric)
us.clips$break.cluster <- get.clusters(us.clips)
us.clips <- us.clips %>% mutate(host.seq=substr(seq, 1, map.length), vir.seq=substr(seq, map.length + 1, read.length), stream="up") %>% group_by(break.cluster) %>% dplyr::mutate(cluster.support = n()) %>% ungroup() %>% filter(cluster.support >= cluster.supp.cutoff) %>% group_by(name) %>% filter(clip.length == max(clip.length)) %>% filter(flag == min(flag)) %>% ungroup()

ds.clips <- fread(left.clip.path, header=F, sep="\t", stringsAsFactors=F) %>% dplyr::rename(name=1, flag=2, chr=3, start.pos=4, map.q=5, cigar=6, seq=7) %>% group_by(name, flag) %>% mutate(read.length=nchar(seq)) %>% ungroup()
ds.clips$clip.length <- as.numeric(sapply(strsplit(ds.clips$cigar, "[A-Z]"), head, 1))
ds.clips <- ds.clips %>% mutate(map.length = read.length - clip.length, break.pos = start.pos, chr.numeric = case_when(grepl("X", chr)~23, grepl("Y", chr)~24, TRUE~as.numeric(gsub("chr", "", chr)))) %>% filter(clip.length >= clip.length.cutoff & map.q >= map.q.cutoff) %>% arrange(chr.numeric, break.pos) %>% dplyr::select(-chr.numeric)
ds.clips$break.cluster <- get.clusters(ds.clips)
ds.clips <- ds.clips %>% mutate(host.seq=substr(seq, clip.length + 1, read.length), vir.seq=substr(seq, 1, clip.length), stream="down") %>% group_by(break.cluster) %>% dplyr::mutate(cluster.support = n()) %>% ungroup() %>% filter(cluster.support >= cluster.supp.cutoff) %>% group_by(name) %>% filter(clip.length == max(clip.length)) %>% filter(flag == min(flag)) %>% ungroup() 

all.clips <- bind_rows(us.clips, ds.clips) %>% group_by(name) %>% filter(n() == 1 | !((stream == "down" & start.pos != max(start.pos)) | (stream == "up" & start.pos != min(start.pos))))
  
print.fastas <- function(clust, reads){
  cluster.sub <- subset(reads, paste(break.cluster, stream, sep="_") == clust)
  stream <- cluster.sub$stream[1]
  fasta.headers <- paste(">", cluster.sub$name, ":", cluster.sub$chr, ":", cluster.sub$start.pos, ":cluster-", cluster.sub$break.cluster, "-", stream, sep="")
  fasta.seqs <- cluster.sub$vir.seq
  cluster.fasta <- data.frame(as.vector(rbind(fasta.headers, fasta.seqs)), row.names=NULL, stringsAsFactors=F)
  write.table(cluster.fasta, file=paste(outdir, "/", samp.id, "_cluster_", clust, "_nonhuman_clips.fa", sep=""), col.names=F, row.names=F, sep="\t", quote=F)
}

all.clusters <- unique(paste(all.clips$break.cluster, all.clips$stream, sep="_"))
sapply(all.clusters, function(x) print.fastas(x, all.clips))

write.table(all.clips, file=paste(outdir, "/", samp.id, ".all.passed.softclip.reads.txt", sep=""), col.names=T, row.names=F, sep="\t", quote=F)
