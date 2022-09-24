suppressMessages(library(data.table))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(GenomicRanges))
suppressMessages(library(TxDb.Hsapiens.UCSC.hg19.knownGene))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(biomartr))
suppressMessages(library(ape))
suppressMessages(library(taxonomizr))
suppressMessages(library(VariantAnnotation))

args <- commandArgs(T)

sample <- args[1]
blast.path <- args[2]
softclips.path <- args[3]
qcov.cutoff <- as.numeric(args[4])
e.cutoff <- as.numeric(args[5])
ident.cutoff <- 100*(as.numeric(args[6]))
clust.perc <- as.numeric(args[7])
outdir <- args[8]

virus.dict <- fread("~/CNV_integration_method/reference/viral_id_map_table.txt", header=F, sep="\t", stringsAsFactors = F) %>% dplyr::rename(id=1, name=2)

refdir <- "/drive3/staging/agarofal/CNV_integration_method/reference/"
taxa.table <- file.path(refdir, "accessionTaxa.sql")

blast.frame <- tryCatch(fread(blast.path, header = F, sep = '\t', stringsAsFactors = F), error=function(e) NULL)
blast.frame <- blast.frame %>% dplyr::rename(qseqid=1, sseqid=2, pident=3, qlen=4, alen=5, mismatch=6, gapopen=7, qstart=8, qend=9, sstart=10, send=11, evalue=12, bitscore=13) %>% separate(qseqid, into=c("name", "cluster", "stream"), sep="-", remove=F) %>% mutate(name = gsub(":chr..*", "", name)) %>% unite("cluster", c("cluster", "stream"), sep="_") %>% mutate(qcov = alen / qlen) %>% filter(qcov >= qcov.cutoff & evalue < e.cutoff)

sc.reads <- fread(softclips.path, header=T, stringsAsFactors=F, sep="\t")

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

txdb.hg19 <- TxDb.Hsapiens.UCSC.hg19.knownGene

extract.breaks <- function(clust){
    cluster.sub <- subset(blast.frame, cluster == clust)
    virus.id <- names(table(cluster.sub$sseqid))[which.max(table(cluster.sub$sseqid))]
    virus.name <- virus.dict$name[virus.dict$id == virus.id]
    virus.taxid <- accessionToTaxa(virus.id, taxa.table)
    consensus.sub <- subset(cluster.sub, sseqid == virus.id)
    virus.breakpos <- ifelse(grepl("up", clust), getmode(consensus.sub$sstart), getmode(consensus.sub$send))
    
    gff.path <- Sys.glob(paste("/drive3/staging/agarofal/CNV_integration_method/reference/*", virus.taxid, "*gff.gz", sep=""))
    
    if(length(gff.path) == 0){
      getGFF(db = "refseq", organism = virus.taxid, path=refdir)
      gff.path <- Sys.glob(paste(refdir, "/*", virus.taxid, "*gff.gz", sep=""))
      rm.command <- paste("rm ", refdir, "/doc* ", refdir, "/*md5*", sep="")
      system(rm.command)
    }
    
    if(length(gff.path) != 0){
      gff <- read_gff(gff.path)
      gff.sub <- gff %>% filter(type == "gene", start <= virus.breakpos & end >= virus.breakpos)
      if(nrow(gff.sub) != 0){
        viral.gene <- paste(sapply(strsplit(sapply(strsplit(gff.sub$attribute, ";"), "[[", 3), "="), "[[", 2), collapse=";")
      } else{
        gene.sub <- gff %>% filter(type == "gene")
        closest.start.sub <- gene.sub[which.min(abs(gene.sub$start - virus.breakpos)),]
        closest.start <- closest.start.sub$start[1]
        
        closest.end.sub <- gene.sub[which.min(abs(gene.sub$end - virus.breakpos)),]
        closest.end <- closest.start.sub$end[1]
        
        closest <- c(closest.start, closest.end)[which.min(abs(c(closest.start, closest.end) - virus.breakpos))]
        
        if(closest == closest.start){
          gff.sub <- closest.start.sub
        } else{
          gff.sub <- closest.end.sub
        }
        if(nrow(gff.sub) != 0){
          viral.gene <- paste(sapply(strsplit(sapply(strsplit(gff.sub$attribute, ";"), "[[", 3), "="), "[[", 2), collapse=";")
        } else{
          viral.gene <- "NA"
        }
      }
    } else{
      viral.gene <- "NA"
    }
    
    if(grepl("up", clust)){
      final.sub <- subset(consensus.sub, sstart == virus.breakpos)
    } else{
      final.sub <- subset(consensus.sub, send == virus.breakpos)
    }
    cluster.name <- paste(unlist(strsplit(clust, "_"))[1:2], collapse="_")
    strm <- unlist(strsplit(clust, "_"))[3]
    
    sc.sub <- sc.reads %>% filter(break.cluster == cluster.name & stream == strm & name %in% final.sub$name)
    
    human.breakchr <- sc.sub$chr[1]
    human.breakpos <- getmode(sc.sub$break.pos)
    human.coord <- paste(paste(gsub("chr", "", human.breakchr), human.breakpos, sep=":"), human.breakpos, sep="-")
    read.support <- nrow(sc.sub)
    
    human.gr <- GRanges(seqnames = human.breakchr, range = IRanges(start=human.breakpos))
    human.annot <- locateVariants(human.gr, txdb.hg19, AllVariants())
    human.type <- as.character(human.annot$LOCATION[human.annot$LOCATION == names(which.max(table(human.annot$LOCATION)))][1])
    cols <- c("SYMBOL")
    keys <- na.omit(unique(human.annot$GENEID))
    human.gene <- ifelse(human.type == "intergenic", "NA", select(org.Hs.eg.db, keys, cols, keytype="ENTREZID")$SYMBOL[1])
    
    readnames <- sc.sub$name
    human.seqs <- sc.sub$host.seq
    viral.seqs <- sc.sub$vir.seq
    
    fasta.headers <- paste(">", readnames, sep="")
    cluster.human.fa <- data.frame(as.vector(rbind(fasta.headers, human.seqs)), row.names=NULL, stringsAsFactors=F)
    cluster.viral.fa <- data.frame(as.vector(rbind(fasta.headers, viral.seqs)), row.names=NULL, stringsAsFactors=F)
    
    human.fa.path <- paste(outdir, "/", sample, "_cluster_", clust, "_human_seqs.fa", sep="")
    viral.fa.path <- paste(outdir, "/", sample, "_cluster_", clust, "_viral_seqs.fa", sep="")
    
    human.cdhit.path <- paste(outdir, "/", sample, "_", clust, "_human.cdhit.out.txt", sep="")
    viral.cdhit.path <- paste(outdir, "/", sample, "_", clust, "_viral.cdhit.out.txt", sep="")
    
    write.table(cluster.human.fa, file=human.fa.path, col.names=F, row.names=F, sep="\t", quote=F)
    write.table(cluster.viral.fa, file=viral.fa.path, col.names=F, row.names=F, sep="\t", quote=F)
    
    cdhit.command <- paste("sh /drive3/staging/agarofal/CNV_integration_method/scripts/run-cdhit.sh", sample, human.fa.path, outdir, clust.perc, sep=" ")
    system(cdhit.command)
    cdhit.command2 <- paste("sh /drive3/staging/agarofal/CNV_integration_method/scripts/run-cdhit.sh", sample, viral.fa.path, outdir, clust.perc, sep=" ")
    system(cdhit.command2)
    
    if(file.exists(human.cdhit.path)){
      human.cdhit <- fread(human.cdhit.path, header=F, stringsAsFactors = F, sep="\t") %>% filter(!grepl(">", V1))
      human.seq <- ifelse(nrow(human.cdhit) == 1, human.cdhit$V1[1], "NA")
    } else{
      human.seq <- "NA"
    }
    
    if(file.exists(viral.cdhit.path)){
      viral.cdhit <- fread(viral.cdhit.path, header=F, stringsAsFactors = F, sep="\t") %>% filter(!grepl(">", V1))
      viral.seq <- viral.cdhit$V1[which.max(nchar(viral.cdhit$V1))]
    } else{
      viral.seq <- "NA"
    }
    
    rm.command <- paste("rm ", outdir, "/*cdhit*", sep="")
    system(rm.command)
    
    if(human.seq != "NA" & viral.seq != "NA"){
      junc.seq <- ifelse(strm == "up", paste(human.seq, viral.seq, sep="|"), paste(viral.seq, human.seq, sep="|"))
    } else{
      junc.seq <- "NA"
    }
    
    output.vec <- c(clust, human.breakchr, human.breakpos, human.type, human.gene, human.seq, virus.id, virus.name, virus.breakpos, viral.gene, viral.seq, strm, junc.seq, read.support)
    return(output.vec)
}

output.table <- do.call(rbind.data.frame, unname(lapply(unique(blast.frame$cluster), function(x) extract.breaks(x)))) %>% dplyr::rename(cluster_name=1, human_chr=2, human_breakpos=3, human_region=4, human_gene=5, human_seq=6, virus_id=7, virus_name=8, virus_breakpos=9, virus_gene=10, virus_seq=11, break_ori=12, breakpt_seq=13, frag_support=14)

write.table(output.table, file=paste(outdir, "/", sample, ".integration.calls.txt", sep=""), col.names=T, row.names=F, sep="\t", quote=F)


