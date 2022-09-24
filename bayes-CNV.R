suppressMessages(library(dplyr))
suppressMessages(library(tools))
suppressMessages(library(Rcpp))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(gtools))
suppressMessages(library(scales))
suppressMessages(library(viridis))
suppressMessages(library(tidyr))
suppressMessages(library(fitdistrplus))
suppressMessages(library(plyr))
suppressMessages(library(gamlss))
suppressMessages(library(VGAM))
suppressMessages(library(foreach))
suppressMessages(library(doSNOW))
suppressMessages(library(mcreplicate))
suppressMessages(library(brms))
suppressMessages(library(tidybayes))
suppressMessages(library(akmedoids))
suppressMessages(library(progress))
suppressMessages(library(pbmcapply))

options(scipen=999)
options(warn=-1)

args <- commandArgs(TRUE)

info.bed.path <- args[1] 
ctr.cov.path <- args[2] 
sample.cov.path <- args[3]
bin.size <- as.numeric(args[4])
bins.pool <- as.numeric(args[5])
num.cores <- as.numeric(args[6])
p.thresh <- as.numeric(args[7])
outdir <- args[8] 
seed <- as.numeric(args[9])
ds.factor <- 3e4*bin.size/3e9

if(is.na(seed)){
  seed <- 123
}

print("Starting off-target CNV analysis...")

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

info.bed <- fread(info.bed.path, header=T, sep="\t", stringsAsFactors=F) %>% dplyr::rename(chr=1, start=2, end=3, AT=4, GC=5, numA=6, numC=7, numG=8, numT=9, numN=10, numOther=11, length=12)
info.bed$bin.num <- seq(1, nrow(info.bed), 1)
info.bed <- info.bed %>% mutate(bad.bin=ifelse(numN >= 100, 1, 0))
info.bed.good <- info.bed %>% filter(bad.bin == 0)

mixedrank <- function(x) order(gtools::mixedorder(x))

if(ctr.cov.path != "NA"){
  print("Reading in CTR off-target coverage data...")

  if(file_ext(ctr.cov.path) != "rda"){
    ctr.read.counts <- fread(ctr.cov.path, header=T, stringsAsFactors=F, sep="\t") %>% filter(chr != "chrM")
    ctr.read.counts <- merge(ctr.read.counts, info.bed.good %>% dplyr::select(chr, start, end, bin.num), by=c("chr", "start", "end")) #%>% arrange(sample, bin.num)
    
    ctr.stats <- ctr.read.counts %>% group_by(sample) %>% dplyr::summarise(RC_median = median(RC), RC_mean = mean(RC), RC_sd = sd(RC))
    
    print("Averaging CTR samples to generate CTR metasample...")
    
    ctr.count.average <- ctr.read.counts %>% group_by(bin.num, chr, start, end) %>% dplyr::summarise(GC = GC[1], RC=mean(RC)) %>% ungroup() %>% dplyr::select(bin.num, chr, start, end, GC, RC) %>% mutate(RC_norm=RC/median(RC)) %>% filter(RC_norm != 0) %>% arrange(bin.num)
    
    c1 <- makeCluster(num.cores)
    registerDoSNOW(c1)
    splits <- sort(rank(1:nrow(ctr.count.average)) %% num.cores)
    ctr.count.average$density <- foreach(i=unique(splits), .combine=c) %dopar% {as.numeric(get_density(ctr.count.average$GC[splits == i], log10(ctr.count.average$RC[splits == i])))}
    stopCluster(c1)
    
    print("GC-correcting CTR metasample...")
    
    set.seed(seed)
    ctr.gc.loess <- loess(RC~GC,data=ctr.count.average %>% sample_n(ds.factor*nrow(ctr.count.average)), control = loess.control(surface = "direct"), degree=2, span=0.1)
    
    c1 <- makeCluster(num.cores)
    registerDoSNOW(c1)
    splits <- sort(rank(1:nrow(ctr.count.average)) %% num.cores)
    preds <- foreach(i=unique(splits), .combine=c) %dopar% {as.numeric(predict(ctr.gc.loess, ctr.count.average$GC[splits == i]))}
    stopCluster(c1)
    
    ctr.count.average <- ctr.count.average %>% mutate(RC_GCcorr = ctr.count.average$RC / preds, RC_GCcorr = ifelse(RC_GCcorr < 0, 0, RC_GCcorr), RC_percentile = ntile(RC_GCcorr, 1000))
  } else{
    load(ctr.cov.path)
  }
  
  print("Generating read count vs GC scatter plot for CTR metasample...")
  jpeg(file=paste(outdir, "/CTR_metasample_RC_GC_scatterplot.jpeg", sep=""), height=3, width=4, res=600, units="in")
  p <- ctr.count.average %>% sample_n(0.05*nrow(ctr.count.average)) %>% ggplot(aes(x=GC, y=RC, color=density)) + geom_point(alpha=0.5, size=0.5) + xlab("% GC") + ylab("Read count") + theme_classic() + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) + annotation_logticks(sides="l") + geom_smooth(method="loess", span=0.05, colour="red", se=F) + scale_color_viridis()
  print(p)
  graphics.off()
  
  write.table(ctr.count.average %>% dplyr::select(-density), file=paste(outdir, "/CTR_metasample_bin_stats.txt", sep=""), col.names=T, row.names=F, sep="\t", quote=F)
  
  blacklist.perc <- 1
  
  bad.bins <- ctr.count.average %>% filter(RC_percentile > (1000 - blacklist.perc*1000/100))
  ctr.no.bad <- ctr.count.average %>% filter(!(bin.num %in% bad.bins$bin.num))
}

print("Reading in sample coverage data and censoring blacklisted bins...")

samp.read.counts <- fread(sample.cov.path, header=T, stringsAsFactors=F, sep="\t") %>% filter(chr != "chrM") %>% mutate(sample=gsub("-T1", "", sample))

if(ctr.cov.path != "NA"){
  samp.read.counts <- merge(samp.read.counts, info.bed.good %>% dplyr::select(chr, start, end, bin.num), by=c("chr", "start", "end")) %>% arrange(bin.num) %>% filter(bin.num %in% ctr.no.bad$bin.num)
} else{
  samp.read.counts <- merge(samp.read.counts, info.bed.good %>% dplyr::select(chr, start, end, bin.num), by=c("chr", "start", "end")) %>% arrange(bin.num)
}

samp.id <- samp.read.counts$sample[1]

samp.stats <- samp.read.counts %>% dplyr::summarise(RC_median = median(RC), RC_mean = mean(RC), RC_sd = sd(RC)) %>% ungroup()

samp.read.counts <- samp.read.counts %>% dplyr::select(sample, bin.num, chr, start, end, GC, RC) %>% mutate(RC_norm = RC/mean(RC))

print("Generating read count vs GC scatter plot for sample...")

samp.read.counts.density <- samp.read.counts %>% filter(RC != 0)

c1 <- makeCluster(num.cores)
registerDoSNOW(c1)
splits <- sort(rank(1:nrow(samp.read.counts.density)) %% num.cores)
samp.read.counts.density$density <- foreach(i=unique(splits), .combine=c) %dopar% {as.numeric(get_density(samp.read.counts.density$GC[splits == i], log10(samp.read.counts.density$RC[splits == i]), n=100))}
stopCluster(c1)

jpeg(file=paste(outdir, "/", samp.id, "_RC_GC_scatterplot.jpeg", sep=""), height=3, width=4, res=600, units="in")
p <- samp.read.counts.density %>% sample_n(0.05*length(sample)) %>% ggplot(aes(x=GC, y=RC, colour=density)) + geom_point(alpha=0.5, size=0.5) + xlab("% GC") + ylab("Read count") + geom_smooth(method="loess", span=0.05, colour="red", se=F) + theme_classic() + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) + annotation_logticks(sides="l") + scale_color_viridis()
print(p)
graphics.off()

print("GC-correcting case sample...")
set.seed(seed)
samp.loess <- loess(RC ~ GC, data = samp.read.counts %>% sample_n(ds.factor*nrow(samp.read.counts)), control = loess.control(surface = "direct"), degree = 2, span = 0.05) 
c1 <- makeCluster(num.cores)
registerDoSNOW(c1)
splits <- sort(rank(1:nrow(samp.read.counts)) %% num.cores)
samp.preds <- foreach(i=unique(splits), .combine=c) %dopar% {as.numeric(predict(samp.loess, samp.read.counts$GC[splits == i]))}
stopCluster(c1)

samp.rc.adj <- samp.read.counts
samp.rc.adj$gc_pred <- samp.preds
samp.rc.adj <- samp.rc.adj %>% mutate(RC_GCcorr = RC / gc_pred) %>% mutate(RC_GCcorr = ifelse(RC_GCcorr < 0 | is.na(RC_GCcorr), 0, RC_GCcorr), RC_norm = ifelse(is.na(RC_norm), 0, RC_norm)) %>% mutate(RC_GCcorr_count=round(RC_GCcorr, 0)) %>% dplyr::select(-gc_pred)

print("Plotting read count histograms for case sample pre- and post-GC correction...")
hist.frame <- samp.rc.adj %>% dplyr::select(sample, RC, RC_GCcorr_count) %>% gather(state, RC, 2:3) %>% dplyr::mutate(state = ifelse(state == "RC", "Raw", "GC-corrected"), state = factor(state, levels=c("Raw", "GC-corrected")), RC=ifelse(RC == 0, 0.5, RC))

zero.counts <- hist.frame %>% group_by(state) %>% dplyr::summarise(zero_perc=round(100*length(sample[RC == 0.5])/length(sample), 2))

samp.rc.nozero <- samp.rc.adj %>% filter(RC_GCcorr != 0 & RC != 0)

jpeg(file=paste(outdir, "/", samp.id, "_RC_histogram_before_after_GC_correction.jpeg", sep=""), height=3, width=4, res=600, units="in")
p <- ggplot(hist.frame, aes(x=RC, fill=state)) + xlab("Per-bin read count") + ylab("Density") + geom_histogram(alpha=0.8, colour="black", position="dodge") + theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.text=element_text(size=12), axis.title=element_text(size=14), legend.text=element_text(size=10), legend.position="top", legend.title=element_blank()) + scale_fill_brewer(palette="Set1") + annotation_logticks(sides="b") + scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) 
print(p)
graphics.off()
#
## Bayesian analysis for updating CNV in local neighborhoods #
#
generate.prior <- function(read.counts, reps, n, dist, sd){
  seeds <- seq(sd, sd + reps - 1, 1)
  sample.prior <- function(sd2){
    set.seed(sd2)
    return(sample(read.counts$RC_GCcorr_count, n, replace=T))
  }
  sim.mat <- pbmcmapply(function(x) sample.prior(x), seeds, mc.cores=num.cores, ignore.interactive = T)
  
  compute.params <- function(sim.data, i){
    column <- sim.data[,i]
    fit <- fitdist(column, dist)
    return(fit$estimate)
  }
  
  index.vec <- seq(1, reps, 1)
  sim.params <- pbmclapply(index.vec, function(x) compute.params(sim.mat, x), mc.cores = num.cores, ignore.interactive=T)
  sim.sizes <- sapply(sim.params, "[[", 1)
  sim.mus <- sapply(sim.params, "[[", 2)
  
  size.weibull.fit <- fitdist(sim.sizes, "weibull")
  size.weibull.shape <- size.weibull.fit$estimate[1]
  size.weibull.scale <- size.weibull.fit$estimate[2]
  
  mu.lnorm.fit <- fitdist(sim.mus, "lnorm")
  mu.lnorm.meanlog <- mu.lnorm.fit$estimate[1]
  mu.lnorm.meansd <- mu.lnorm.fit$estimate[2]
  
  prior.params <- list(size=c(size.weibull.shape, size.weibull.scale), mu=c(mu.lnorm.meanlog, mu.lnorm.meansd))
  
  return(prior.params)
}

print("Generating prior...")
prior <- generate.prior(samp.rc.adj, 1e4, 1e4, "nbinom", seed)

group.bins <- function(read.counts, window.size){
  rc.grouped <- read.counts %>% group_by(chr) %>% dplyr::mutate(grp = rep(row_number(), length.out = n(), each = window.size), window.id = paste(chr, grp, sep=":"))
  window.id <- unique(rc.grouped$window.id)
  window <- seq(1, length(window.id), 1)
  map.frame <- data.frame(cbind(window.id, window), stringsAsFactors=F, row.names=NULL) %>% mutate(window = as.numeric(window))
  rc.grouped <- merge(rc.grouped, map.frame, by="window.id") %>% dplyr::select(-grp, -window.id) %>% arrange(window, start)
  return(rc.grouped)
}
samp.rc.grouped <- group.bins(samp.rc.adj, bins.pool)

window.summary <- samp.rc.grouped %>% group_by(window) %>% dplyr::summarise(window.chr=chr[1], window.start=min(start), window.end=max(end), window.length=window.end-window.start, window.min=min(RC_GCcorr_count), window.max=max(RC_GCcorr_count), window.median=median(RC_GCcorr_count), window.75=unname(quantile(RC_GCcorr_count))[4], window.mean=mean(RC_GCcorr_count), window.sd=sd(RC_GCcorr_count), window.zeros=length(RC_GCcorr_count[RC_GCcorr_count == 0])/length(RC_GCcorr_count))

zero.windows <- subset(window.summary, window.zeros == 1)$window

write.table(window.summary, file=paste(outdir, "/", samp.id, "_segment_summary_stats.txt", sep=""), col.names=T, row.names=F, sep="\t", quote=F)

print("Calculating posteriors...")

generate.posterior <- function(grouped.bin.counts, pr, reps, nchain, sd){  
  single.post <- function(window.num, sd2){
    
    dat <- subset(grouped.bin.counts, window == window.num)
    if(nrow(dat) < bins.pool){
      set.seed(sd2)
      dat <- dat %>% sample_n(bins.pool, replace=T)
    }
    
    stanvars <- stanvar(unname(pr$mu[1]), name='m1') + stanvar(unname(pr$mu[2]), name='m2') + stanvar(unname(pr$size[1]), name='s1') + stanvar(unname(pr$size[2]), name='s2')
    
    stan.prior <- c(prior(lognormal(m1, m2), class="Intercept", lb=0), prior(weibull(s1, s2), class="shape"))
    
    stan.model <- brm(data = dat, family=negbinomial, RC_GCcorr_count ~ 1, prior = stan.prior, cores=1, iter=1.2*reps, warmup=0.2*reps, chains=nchain, stanvars=stanvars, seed=sd2, silent=2, refresh=0, open_progress=F)
    
    set.seed(sd2)
    stan.pp <- as.numeric(posterior_predict(stan.model, ndraws = 1, cores = 1))
    return(stan.pp)
  }
  
  winds <- unique(grouped.bin.counts$window)
  c1 <- makeCluster(num.cores)
  registerDoSNOW(c1)
  pb <- txtProgressBar(min=1, max=length(winds), style=3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
  window.posteriors <- foreach(i=sort(winds), .packages=c("brms", "rstan", "StanHeaders", "Rcpp", "dplyr"), .options.snow=opts, .combine='rbind', .export='bins.pool', .inorder=T) %dopar% {as.numeric(single.post(i, sd))}
  close(pb)
  stopCluster(c1)
  
  rownames(window.posteriors) <- sort(winds)
 
  return(window.posteriors)
}

posteriors <- generate.posterior(samp.rc.grouped, prior, 1000, 3, seed)

plot.all.post <- function(pr, pst, n, sd){
  set.seed(sd)
  prior.sizes <- rweibull(n, shape=pr$size[1], scale=pr$size[2])
  set.seed(sd)
  prior.mus <- rlnorm(n, meanlog=pr$mu[1], sdlog=pr$mu[2])
  
  prior.sim <- matrix(cbind(prior.mus, prior.sizes), ncol=2)
  prior.pred <- rep(rnbinom(nrow(prior.sim), mu=prior.sim[, 1], size=prior.sim[, 2]), length(unique(window.summary$window.chr)))
  prior.dist <- rep("prior", length(prior.pred))
  prior.chr <- rep(unique(window.summary$window.chr), each=n)
  prior.frame <- data.frame(cbind(prior.dist, prior.chr, prior.pred, prior.dist), stringsAsFactors=F, row.names = NULL) %>% dplyr::rename(window=1, window.chr=2, value=3, dist=4) %>% mutate(value=as.numeric(value))
  
  window.chr.sub <- window.summary %>% dplyr::select(window, window.chr)
  
  post.frame <- data.frame(t(pst)) 
  names(post.frame) <- gsub("^X", "", names(post.frame))
  post.frame <- merge(post.frame %>% gather(window, value) %>% mutate(window = as.numeric(window), dist="posterior") %>% arrange(window), window.chr.sub, by="window") %>% dplyr::select(window, window.chr, value, dist)
  
  bayes.frame <- rbind(prior.frame, post.frame) %>% mutate(value = as.numeric(value)) %>% mutate(dist=factor(dist, levels=c("prior", "posterior")), col=factor(ifelse(dist == "prior", "WG", window.chr)))
  
  cols <- rainbow(length(unique(bayes.frame$col)) - 1)
  
  bayes.frame.plot <- bayes.frame %>% mutate(value = ifelse(value == 0, 0.5, value), window.chr = factor(window.chr, levels=unique(window.summary$window.chr)), col = factor(col, levels=c(unique(window.summary$window.chr), "WG")))

  p <- ggplot(bayes.frame.plot, aes(x=value, group=window, colour=col, size=dist)) + xlab("GC-corrected read count") + ylab("Density") + geom_density() + theme_classic() + scale_size_manual(values = c("prior" = .6, "posterior"=0.08), guide="none") + theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.text=element_text(size=12), axis.title=element_text(size=14), legend.text=element_text(size=10), legend.position="none", legend.title=element_blank()) + scale_colour_manual(values=c(cols, "black")) + annotation_logticks(sides="b") + scale_x_log10(limits=c(0.4, .01*max(bayes.frame.plot$value)), breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) + facet_wrap(~window.chr, ncol=3, scales="free_y")
  
  return(p)
}

print("Plotting posterior distributions...")
jpeg(file=paste(outdir, "/", samp.id, "_all_posterior_distributions.jpeg", sep=""), height=10, width=8, res=600, units="in")
all.post <- plot.all.post(prior, posteriors, 1000, seed)
plot(all.post)
graphics.off()

print("Generating probabilities from posterior distributions...")

convert.post.to.params <- function(pst){
  windows <- rownames(pst)
  get.params <- function(wind){
    posterior.vec <- pst[rownames(pst) == wind, ]
    post.fit <- fitdist(posterior.vec, "nbinom", lower=c(0,0), start=list(size=1, mu=1))
    post.fit.params <- post.fit$estimate
    return(post.fit.params)
  }
  window.params <- do.call(rbind.data.frame, pbmclapply(windows, function(x) get.params(x), ignore.interactive=T, mc.cores=num.cores)) %>% dplyr::rename(size=1, mu=2)
  window.params <- data.frame(cbind(windows, window.params), stringsAsFactors=F, row.names=NULL) %>% dplyr::rename(window=1)
  return(window.params)
}

window.params <- convert.post.to.params(posteriors)

samp.rc.grouped <- merge(samp.rc.grouped, window.params, by="window") %>% group_by(bin.num) %>% mutate(bin.p=(1-pnbinom(RC_GCcorr_count, size=size, mu=mu))) %>% ungroup() %>% dplyr::select(-c("size", "mu"))

print("Plotting histogram of predictive posterior p-values...")
jpeg(file=paste(outdir, "/", samp.id, "_posterior_p_histogram.jpeg", sep=""), width=5, height=4, res=600, units="in")
plot <- ggplot(samp.rc.grouped, aes(x=bin.p)) + geom_density(colour="black", fill="navyblue", alpha=0.5, size=1) + xlab("Posterior predictive p-value") + ylab("Density") + theme_classic() + theme(axis.title = element_text(size=14), axis.text=element_text(size=12), panel.border=element_rect(colour="black", fill=NA, size=1.5), legend.background=element_blank())
print(plot)
dev.off()

print("Trying different p-value thresholds...")
p.increments <- seq(0,0.1, 0.00001)

cut.p <- function(p){
  n <- ifelse(p == 0, nrow(subset(samp.rc.grouped, bin.p <= p)), nrow(subset(samp.rc.grouped, bin.p < p)))
  return(n)
}
pos.counts <- unlist(pbmclapply(p.increments, function(x) cut.p(x), mc.cores = num.cores, ignore.interactive=T))
thresh.hits <- data.frame(cbind(p.increments, pos.counts), row.names=NULL, stringsAsFactors = F) %>% dplyr::rename(p.cutoff=1, hits=2) %>% mutate(gen.perc=hits*bin.size/3e9)

elbow.x <- elbow_point(thresh.hits$p.cutoff, thresh.hits$gen.perc)$x
elbow.y <- elbow_point(thresh.hits$p.cutoff, thresh.hits$gen.perc)$y

elbow.x.log <- elbow_point(thresh.hits$p.cutoff, log10(thresh.hits$gen.perc))$x
elbow.y.log <- elbow_point(thresh.hits$p.cutoff, log10(thresh.hits$gen.perc))$y

jpeg(file=paste(outdir, "/", samp.id, "_posterior_p_cutoffs_vs_perc_genome_hits.jpeg", sep=""), height=3, width=4, res=600, units="in")
p <- thresh.hits %>% ggplot(aes(x=p.cutoff, y=gen.perc)) + geom_point(size=0.3, colour="navyblue", alpha=0.4) + geom_vline(xintercept=elbow.x, colour="black", linetype="dashed", size=0.5) + geom_hline(yintercept=elbow.y, colour="black", linetype="dashed", size=0.5) + stat_smooth(method="lm", formula = y ~ poly(x, 25), se=F, inherit.aes=T, size=0.5, colour="red") + xlab("Posterior p-value cutoff") + ylab("Fraction of genome\namp-positive") + theme_bw() + scale_x_continuous(breaks = seq(0, 0.1, 0.01))
print(p)
graphics.off()

jpeg(file=paste(outdir, "/", samp.id, "_posterior_p_cutoffs_vs_log_perc_genome_hits.jpeg", sep=""), height=3, width=4, res=600, units="in")
p <- thresh.hits %>% ggplot(aes(x=p.cutoff, y=gen.perc)) + geom_point(size=0.3, colour="navyblue", alpha=0.4) + geom_vline(xintercept=elbow.x.log, colour="black", linetype="dashed", size=0.5) + geom_hline(yintercept=10^elbow.y.log, colour="black", linetype="dashed", size=0.5) + stat_smooth(method="lm", formula = y ~ poly(x, 25), se=F, inherit.aes=T, size=0.5, colour="red") + xlab("Posterior p-value cutoff") + ylab("Fraction of genome\namp-positive") + theme_bw() + scale_x_continuous(breaks = seq(0, 0.1, 0.01)) + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) + annotation_logticks(sides="l")
print(p)
graphics.off()

print("Writing out bed file...")
sig.bed <- samp.rc.grouped %>% filter(bin.p < p.thresh) %>% dplyr::select(chr, start, end, bin.p)
write.table(sig.bed, file=paste(outdir, "/", samp.id, ".sig_bins.bed", sep=""), col.names=F, row.names=F, sep="\t", quote=F)

write.table(samp.rc.grouped, file=paste(outdir, "/", samp.id, "_bin_stats.txt", sep=""), col.names=T, row.names=F, sep="\t", quote=F)

print("Done.")
