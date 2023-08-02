options(width=180)

load("raw_tox.rda")

# sample names are x_Sx_L00y where x is the sample number and y is the lane number in [1,4]
# there are 120 ASPA samples, 32 on toxicity and 88 VAMP samples

isamples  = seq(32)
dates     = rep(seq(4), each=8)
time      = rep(c(0,0,5,5,7,7,9,9), 4)
tech_rep  = rep(seq(2), 16)

samples = data.frame(date=dates, time=time, tech_rep=tech_rep, sample_id=sprintf("%d_S%d",isamples,isamples))

# ===========================================================================================================================================
#
# Lane merged (lm)
#
raw_lm = data.frame(barcode=raw$barcode, var=raw$var_aa)
lane_rp = c()
lane_nreads = c()
# loop over all samples
for (i in seq(nrow(samples))) {
    cn_lanes = paste0(samples[i,"sample_id"],"_L00",seq(4))
    raw_lm[,samples[i,"sample_id"]] = apply(raw[,cn_lanes], MARGIN=1, sum)

    # Correlations
    rp = c(cor(raw[,cn_lanes[1]], raw[,cn_lanes[2]], method="pearson"))
    rp = append(rp, cor(raw[,cn_lanes[1]], raw[,cn_lanes[3]], method="pearson"))
    rp = append(rp, cor(raw[,cn_lanes[1]], raw[,cn_lanes[4]], method="pearson"))
    rp = append(rp, cor(raw[,cn_lanes[2]], raw[,cn_lanes[3]], method="pearson"))
    rp = append(rp, cor(raw[,cn_lanes[2]], raw[,cn_lanes[4]], method="pearson"))
    rp = append(rp, cor(raw[,cn_lanes[3]], raw[,cn_lanes[4]], method="pearson"))
    print(sprintf("Lane correlations of %6s: %s", samples[i,"sample_id"], paste(rp,collapse=", ")))
    lane_rp = c(lane_rp, rp)

    # Read counts per lane
    nreads = apply(raw[,cn_lanes], MARGIN=2, sum)
    lane_nreads = c(lane_nreads, nreads)
    print(sprintf("Number of reads for lanes: %s", paste(nreads, collapse=", ")))
}
print(sprintf("Lane correlations %.2f sd %.2f min %.2f max %.2f",mean(lane_rp),sd(lane_rp),min(lane_rp),max(lane_rp)))
print(sprintf("Average reads per lane %.0f sd %.0f min %d max %d",mean(lane_nreads),sd(lane_nreads),min(lane_nreads),max(lane_nreads)))


# ===========================================================================================================================================
#
# Technical replica merged (rm)
#
rep_nreads = c()
rep_rp_raw = c()
raw_rm = data.frame(barcode=raw$barcode, var_aa=raw$var_aa, var_dna=raw$var_dna)
for (i in which(samples$tech_rep==1)) {
    # sample row of other tech_rep 
    i2 = which(samples$date     == samples[i,"date"]      &  samples$time     == samples[i,"time"] &
               samples$tech_rep == 2)
    stopifnot(length(i2) == 1)
    cn_rep = c(samples[i,"sample_id"], samples[i2,"sample_id"])
    cn = sprintf("run%d_day%d", samples[i,"date"], samples[i,"time"])
    raw_rm[,cn] = apply(raw_lm[,cn_rep], MARGIN=1, sum)

    # Correlations
    rp_raw = cor(raw_lm[,cn_rep[1]], raw_lm[,cn_rep[2]], method="pearson")
    rep_rp_raw = c(rep_rp_raw, rp_raw)
    
    # Read counts per replica
    nreads = apply(raw_lm[,cn_rep], MARGIN=2, sum)
    rep_nreads = c(rep_nreads, nreads)
    print(sprintf("Pearson between technical replica %6s and %6s %.2f. Reads %7d and %7d",cn_rep[1],cn_rep[2],rp_raw,nreads[1],nreads[2]))
}
print(sprintf("Replica raw correlations %.2f sd %.2f min %.2f max %.2f",mean(rep_rp_raw),sd(rep_rp_raw),min(rep_rp_raw),max(rep_rp_raw)))
print(sprintf("Average reads per replica %.0f sd %.0f min %d max %d",mean(rep_nreads),sd(rep_nreads),min(rep_nreads),max(rep_nreads)))

print(sprintf("Average coverage of %d barcodes is x%.2f in each technical replicate", nrow(raw), mean(rep_nreads)/nrow(raw)))

# # rename and remove barcodes without reads
# raw_rm$sum_all = apply(raw_rm[,sprintf("run%d_day%d",rep(seq(4),each=4),rep(c(0,5,7,9),4))], MARGIN=1, sum)
# print(sprintf("Remove %d barcodes without reads observed in any replica", sum(raw_rm$sum_all==0)))
# raw_bc = raw_rm[which(raw_rm$sum_all > 0),]

# keep variants without reads, it gives the information at a barcode should be present in the lib
raw_bc = raw_rm

# ===========================================================================================================================================
#
# Technical replicates have poor correlations and much worse than the vamp data?
#

# # Look for patterns per date or time
# print("Pearson correlation of raw counts per growth timepoint (averaged over 4 dates/experiments)")
# rep_rp_timeavg = aggregate(rep_rp_raw, by=list(samples[samples$tech_rep==1,"time"]), mean)
# colnames(rep_rp_timeavg) = c("time","mean")
# print(rep_rp_timeavg)

# print("Pearson correlation of raw counts per date/experiment (averaged over growth timepoints)")
# rep_rp_dateavg = aggregate(rep_rp_raw, by=list(samples[samples$tech_rep==1,"date"]), mean)
# colnames(rep_rp_dateavg) = c("date","mean")
# print(rep_rp_dateavg)

# Perhaps a small dependency on growth timepoint with growth time zero having the worst coorlation.
# Growth point "Day 0" should have the highest complexity. Calculate number of genotypes above a threshold (that relates to the distribution)
#   and perhaps an effective complexity like log-differential entropy of the distribution.

# Does reads-per-barcode correlate with replica correlation?
cns = paste0("run",rep(unique(samples$date),each=4),"_day",rep(unique(samples$time),4))
names(rep_rp_raw) = cns
read_stats = data.frame(name=cns, reads=apply(raw_bc[,cns], MARGIN=2, sum), rp_raw=rep_rp_raw)
# number of barcodes with non-zero number of reads
read_stats$nbc_nz = apply(raw_bc[,cns], MARGIN=2, function(v){ sum(v>0) })
# effective number of barcodes with reads, based on p*log(p)
read_stats$nbc_eff = apply(raw_bc[,cns], MARGIN=2, function(v){ p=v/sum(v)+1e-12; exp(-sum(p*log(p))) })
# average reads per barcode (non-zero)
read_stats$avg_rpb = read_stats$reads / read_stats$nbc_eff
# correlate everything
quartz(width=8, height=8)
plot(read_stats[,c("rp_raw","nbc_nz","nbc_eff","avg_rpb")])
quartz.save("aspa_tox_complexity_and_reads.png", type="png")


# Complexity as a function of growth time
read_stats$run = rep(seq(4), each=4)
read_stats$day = rep(c(0,5,7,9), 4)
quartz(width=8, height=6)
plot(0,0,col=0, xlim=c(0,9), ylim=c(min(read_stats$nbc_eff),max(read_stats$nbc_eff)), xlab="Growth time [days]", ylab="Complexity [effective number of barcodes]")
for (run in seq(4)) { i=which(read_stats$run==run); lines(read_stats[i,"day"], read_stats[i,"nbc_eff"], col=run, lwd=2) }
legend("topright", sprintf("Experiment %d",seq(4)), lty=1, lwd=2, col=seq(4))
quartz.save("aspa_tox_complexity_and_time.png", type="png")


# # Histograms of reads per barcode
# col = rep(seq(4), each=4); names(col) = cns
# lty = rep(seq(4), 4); names(lty) = cns
# max_counts = max(raw_bc[,cns])
# quartz(width=12, height=8)
# plot(0,0,col=0, xlim=c(0,100), ylim=c(0,6000), xlab="Read counts", ylab="Barcodes")
# for (cn in cns) {
#     h = hist(raw_bc[,cn], breaks=seq(0,max_counts), plot=F)
#     lines(h$mids, h$counts, lty=lty[cn], col=col[cn])
#     nbc = sum(raw_bc[,cn]>0)
#     p = raw_bc[,cn]/sum(raw_bc[,cn])
#     nbc_eff = exp(-sum(p*log(p+1e-12)))
#     print(sprintf("%s: Barcodes with counts: %6d. Mean reads per barcode %4.1f. Pearson %4.2f. Effective number of barcodes %7.1f",
#                    cn,nbc,sum(raw_bc[,cn])/nbc,rep_rp_raw[cn],nbc_eff))
# }
# legend("topright",cns,lty=lty,col=col)
# quartz.save("tox_read_distributions.png", type="png")


# ===========================================================================================================================================
#
# Fit growth profiles per barcode
#
# The read frequency should be scaled by the complexity of the library at a given time point, e.g. if the complexity is 100 at
#   time 0 but 10 at time 1, frequencies will on average increase a factor of 10. 

# Function that normalize all columns to RPM
reads_per_million = function(df) {
    col_sums = apply(df, MARGIN=2, sum)
    df * matrix(10^6/col_sums, nrow=dim(df)[1], ncol=dim(df)[2], byrow=TRUE)
}

line_fit_fast = function(y,x,sum_x,sum_xx,delta) {
    # # normalize to first time point
    # y = y/y[1]
    # # y = log( y/y[1] )
    
    # # normalize as a distribution - this seems to give a robust slope on a relevant scale to compare between variants
    # y = y/sum(y)
    # sum_y = 1.0
    
    sum_y = sum(y)
    sum_x_y=sum(x*y)
    slope = (length(x)*sum_x_y - sum_x*sum_y)/delta
    intersept = (sum_xx*sum_y - sum_x*sum_x_y)/delta
    # return( c(intersept,slope) )
    residuals = y - (intersept + slope*x)
    r2 = 1-sum(residuals^2)/sum((y-mean(y))^2)
    return( c(intersept,slope,r2) )    
}

line_fit_lm = function(y, x) {
    # if (y[1] > 0) {
    #     # normalize to first point
    #     y = y/y[1]
        fit = lm(y ~ x)
        r2 = 1-sum(fit$residuals^2)/sum((y-mean(y))^2)
        return( c(coef(fit), r2) )
    # } else {
    #     return( c(NA,NA,NA) )
    # }
}

line_fit_lm_w = function(yw, x) {
    n = length(x)
    w = unlist( yw[(n+1):(2*n)] )
    y = unlist( yw[1:n] )
    # if (y[1] > 0) {
        # # normalize to first point
        # y = y/y[1]
        fit = lm(y ~x, weights=w)
        r2 = 1-sum(fit$residuals^2)/sum((y-mean(y))^2)
        return( c(coef(fit), r2) )
    # } else {
    #     return( c(NA,NA,NA) )
    # }
}

growth_profile = function(df, x, name, wdf=NULL) {
    stopifnot(length(x) == ncol(df))

    if (is.null(wdf)) {
        # a simple line - fast
        sum_x = sum(x)
        sum_x_sq = sum(x*x) #namespace check! this is sum_xx
        delta = length(x)*sum(x*x) - sum_x*sum_x    
        linefit = t(apply(df, MARGIN=1, line_fit_fast, x=x, sum_x=sum_x, sum_xx=sum_x_sq, delta=delta))
        # # linear regression using lm
        # linefit = t(apply(df, MARGIN=1, line_fit_lm, x=x))
    } else {
        # weighted least squares
	# There's something odd about always having low weights on all low read counts, e.g. y=c(1000,1000,0,0) or c(1,1,0,0) would have a
	#   slope of 0 and excellent fit but I would consider the first dying with some certainty and the other highly uncertain	
	# So don't use weighted least squares!
        linefit = t(apply(cbind(df,wdf), MARGIN=1, line_fit_lm_w, x=x))
    }
    
    # put in a data frame and give column name
    ret_df = data.frame(linefit)
    
    # ret_df$rs = apply(df, MARGIN=1, function(v){ stddev=sd(v); if (stddev < 1e-3) NA else cor(seq(length(v)),v,method="spearman") })

    # colnames(ret_df) = paste0(name,c(".slope",".intersept",".rs"))
    colnames(ret_df) = paste0(name,c("_intersept","_slope","_r2"))
    return(ret_df)
}

# Calculate count frequencies with psudo_counts
psudo_counts = 0
print(sprintf("Normalizing barcode counts to frequencies using psudo count %.3f",psudo_counts))
rpm_bc = data.frame(barcode=raw_bc$barcode, var_aa=raw_bc$var_aa)
cns = sprintf("run%d_day%d", rep(seq(4),each=4), rep(c(0,5,7,9),4))
rpm_bc = cbind(rpm_bc, reads_per_million(raw_bc[,cns] + psudo_counts))

# # NA out barcode counts below a threshold - note that these still count in the frequencies
# read_cut = 10
# for (cn in cns) { i=which(raw_bc[,cn] < read_cut); rpm_bc[i,cn] = NA}

# Growth per profile calculated per barcode
growth_bc = data.frame(barcode=raw_bc$barcode, var_aa=raw_bc$var_aa, var_dna=raw_bc$var_dna)
# add number of aa substitutions
splitlist_aa = strsplit(growth_bc$var_aa, ":")
growth_bc$n_aa = sapply(splitlist_aa, length)
# add number of DNA substitutions
splitlist_dna = strsplit(growth_bc$var_dna, ":")
growth_bc$n_dna = sapply(splitlist_dna, length)

for (si in which(samples$time==0 & samples$tech_rep==1)) {
    cn_times = sprintf("run%d_day%d",samples[si,"date"],c(0,5,7,9))
    cn = sprintf("run%d", samples[si,"date"])
    
    print(sprintf("Calculate growth profiles per barcode for %s",cn))
    row_sums = data.frame(rs = apply(rpm_bc[,cn_times], MARGIN=1, sum))
    df = rpm_bc[,cn_times] / row_sums[,rep(1,length(cn_times))]
    growth_bc = cbind(growth_bc, growth_profile(df, x=c(0,5,7,9), name=cn))
    
    # # Weighted least squares
    # wdf = raw_bc[,cn_times]
    # growth_bc = cbind(growth_bc, growth_profile(df, x=c(0,5,7,9), name=cn, wdf=wdf))
}

# Sum raw reads per growth curve and calc an unceratinty from this
for (run in unique(samples$date)) {
    cn_times = sprintf("run%d_day%d", run, c(0,5,7,9))
    cn_raw = sprintf("run%d_sum",run)
    raw_bc[,cn_raw] = apply(raw_bc[,cn_times], MARGIN=1, sum)
    growth_bc[,sprintf("run%d_reads",run)] = raw_bc[,cn_raw]
    cn_grow = sprintf("run%d_nu",run)
    growth_bc[,cn_grow] = 1 - 1/sqrt(raw_bc[,cn_raw]/4)
    growth_bc[which(is.infinite(growth_bc[,cn_grow])), cn_grow] = NA
}
cns = sprintf("run%d_day%d",rep(seq(4),each=4),rep(c(0,5,7,9),4))
growth_bc$all_reads = apply(raw_bc[,cns], MARGIN=1, sum)
stopifnot( all( growth_bc$all_sum == apply(growth_bc[,sprintf("run%d_reads",seq(4))], MARGIN=1, sum) ) )

# Calculate a mean of replica
cn = sprintf("run%d_slope",seq(4))
growth_bc$mean.slope = apply(growth_bc[,cn], MARGIN=1, mean, na.rm=T)
growth_bc$sd.slope = apply(growth_bc[,cn], MARGIN=1, sd, na.rm=T)

# Min-max normalization
var_tox_ctrl = c("C152W")
bc_slope_tox = mean( growth_bc[which(growth_bc$var_aa %in% var_tox_ctrl),"mean.slope"], na.rm=T )
bc_slope_wt = mean( growth_bc[which(growth_bc$var_dna==""),"mean.slope"], na.rm=T )
bc_scale = 1/(bc_slope_tox - bc_slope_wt)
growth_bc$mean.slope.norm = (growth_bc$mean.slope - bc_slope_wt) * bc_scale
growth_bc$sd.slope.norm = growth_bc$sd.slope * bc_scale

# Tox-score correlations between biological replicates
i = which(growth_bc$n_aa == 1 )
print("Barcode slope correlation of single AA variants between biological replicates (not normalized)")
rps = c(); maes = c()
for ( run1 in seq(1,3) ) {
    x1 = (growth_bc[i,sprintf("run%d_slope",run1)] - bc_slope_wt) * bc_scale
    for ( run2 in seq(run1+1,4) ) {
        x2 = (growth_bc[i,sprintf("run%d_slope",run2)] - bc_slope_wt) * bc_scale
        rp = cor(x1, x2, method="pearson", use="complete.obs")
	mae = mean(abs(x1-x2), na.rm=T)
	rps = c(rps, rp)
	maes = c(maes, mae)
	print(sprintf("Run %d versus run %d: rp = %.3f, MAE = %.3f",run1,run2,rp,mae))
    }
}
print(sprintf("Average of %d pairs: rp %.3f (range %.3f-%.3f), MAE %.3f (range %.3f-%.3f)",
              length(rps),mean(rps),min(rps),max(rps),mean(maes),min(maes),max(maes)))

# # Distribution of slopes for WT barcodes
# i_wt = which(growth_bc$var_aa == "")
# quartz(height=6, width=6)
# plot(0, 0, col=0, xlim=c(-1.5,0.5), ylim=c(0,4), xlab="Slope", ylab="Density", main=sprintf("%d wild-type barcodes",length(i_wt)))
# breaks = seq(-3,13,0.1)
# for (run in seq(4)) {
#     cn = sprintf("run%d_slope",run)
#     h = hist(growth_bc[i_wt,cn], breaks=breaks, plot=F)
#     lines(h$mids, h$density, col=run)
# }
# legend("topright", sprintf("run %d",seq(4)), lty=1, col=seq(4))


# # check calculations by plotting a few random curves
# i_ex = sample(seq(nrow(growth_bc)), size=4, replace=F)
# i_ex = which(raw_bc$var_aa == "C152W")
# i_ex = sample(i_wt, size=5, replace=F)
# quartz(height=6, width=8)
# plot(0, 0, col=0, xlim=c(0,9), ylim=c(0,5), xlab="Growth time [days]", ylab="", main=sprintf("Plotting growth curves for %d barcodes", length(i_ex)))
# for (ii in seq_along(i_ex)) {
#     i = i_ex[ii]
#     for (run in seq(4)) {
#         cns = sprintf("run%d_day%d",run,c(0,5,7,9))
#         # if (rpm_bc[i,cns[1]] > 10) {
#         if (all(! is.na(rpm_bc[i,cns])) & rpm_bc[i,cns[1]] > 0) {
# 	    print(log10(raw_bc[i,cns]+2))
#             points(c(0,5,7,9), rpm_bc[i,cns] / rpm_bc[i,cns[1]], type="b", col=ii, lty=run, cex=unlist(log10(raw_bc[i,cns]+2)))
#             abline(growth_bc[i,sprintf("run%d_intersept",run)], growth_bc[i,sprintf("run%d_slope",run)], lwd=2, col=ii, lty=run) # col=run, lty=ii)
# 	}
#     }
# }
# legend("topleft", c(sprintf("run %d",seq(4)), sprintf("barcode %d",i_ex)), col=c(rep(1,4),seq_along(i_ex)), lty=c(seq(4),rep(1,length(i_ex))), lwd=2, ncol=1)


# ===========================================================================================================================================
#
# Aggregate per amino acid variant
#
print("Aggregate barcodes per amino acid variant")
agg = aggregate(seq(nrow(rpm_bc)), list(rpm_bc$var_aa), paste, collapse=":")
var_indices = lapply( strsplit(agg[,2], ":"), as.numeric)
names(var_indices) = agg[,1]

raw_aa = data.frame(var=names(var_indices), n_bc=sapply(var_indices,length))
for (i in which(samples$tech_rep==1)) {
    cn = sprintf("run%d_day%d",samples[i,"date"],samples[i,"time"])
    raw_aa[,cn] = sapply(var_indices, function(v) { sum(raw_bc[v,cn], na.rm=T) })
}

# Total number of reads per replica
for (run in unique(samples$date)) {
    cns = sprintf("run%d_day%d",run,c(0,5,7,9))
    cn = sprintf("run%1d_sum",run)
    raw_aa[,cn] = apply(raw_aa[,cns], MARGIN=1, sum)
}


# Compare, replica correlation in raw counts
rps_bc = c()
for (day in unique(samples$time)) {
    print(sprintf("Replica correlations of raw counts per barcode at day %d",day))
    # log counts to avoid inflated correlations due to jackpotting
    # df = log( raw_bc[,sprintf("run%d_day%d",seq(4),day)] + 1 )
    df = raw_bc[,sprintf("run%d_day%d",seq(4),day)]

    # # plot to look for jackpotting
    # quartz(width=8, height=8)
    # plot(df, pch=".")
    # # quartz.save(sprintf("replica_bc_cor_day%d.png",day), type="png")
    
    # correlations
    cm = cor(df, method="pearson")
    # print(cm)
    rps_bc = c(rps_bc, cm[upper.tri(cm)])
    print(sprintf("Average Pearson correlation %.4f range %.4f - %.4f",mean(cm[upper.tri(cm)]),min(cm[upper.tri(cm)]),max(cm[upper.tri(cm)])))
}
print(sprintf("Average of %d Pearson barcode correlations %.4f range %.4f - %.4f",length(rps_bc),mean(rps_bc),min(rps_bc),max(rps_bc)))

rps_aa = c()
for (day in unique(samples$time)) {
    print(sprintf("Replica correlations of raw counts per amino acid variant at day %d",day))
    # exclude WT with very high read counts and/or log counts to avoid inflated correlations
    # df = log( raw_aa[,sprintf("run%d_day%d",seq(4),day)] + 1 )
    df = raw_aa[2:nrow(raw_aa),sprintf("run%d_day%d",seq(4),day)]
    
    # # plot to look for jackpotting
    # quartz(width=8, height=8)
    # plot(df)
    # # quartz.save(sprintf("replica_aa_cor_day%d.png",day), type="png")
    
    # correlations
    cm = cor(df, method="pearson")
    rps_aa = c(rps_aa, cm[upper.tri(cm)])
    # print(cm)
    # print(sprintf("Average Pearson correlation %.4f range %.4f - %.4f",mean(cm[upper.tri(cm)]),min(cm[upper.tri(cm)]),max(cm[upper.tri(cm)])))
}
print(sprintf("Average of %d Pearson AA correlations %.4f range %.4f - %.4f",length(rps_aa),mean(rps_aa),min(rps_aa),max(rps_aa)))


# Now, sequencing coverage is much better
cn = paste0("run",rep(unique(samples$date),each=4),"_day",rep(unique(samples$time),4))
sample_coverage = apply(raw_bc[,cn], MARGIN=2, sum, na.rm=T) / nrow(raw_bc)
print(sprintf("Per-barcode sequencing coverage average x%.2f (range %.2f-%.2f) [reads per barcode]",
              mean(sample_coverage), min(sample_coverage), max(sample_coverage)))
print(sample_coverage)

cn = paste0("run",rep(unique(samples$date),each=4),"_day",rep(unique(samples$time),4))
sample_coverage = apply(raw_aa[,cn], MARGIN=2, sum, na.rm=T) / nrow(raw_aa)
print(sprintf("Per AA variant sequencing coverage average x%.2f (range %.2f-%.2f) [reads per variant]",
              mean(sample_coverage), min(sample_coverage), max(sample_coverage)))
print(sample_coverage)


# Number of reads
raw_aa_sum = apply(raw_aa[,sprintf("run%d_day%d",rep(seq(4),each=4),rep(c(0,5,7,9)))], MARGIN=2, sum)
print(sprintf("Average number of matched reads per AA variant %.0f range %d - %d", mean(raw_aa_sum), min(raw_aa_sum), max(raw_aa_sum)))


# ===========================================================================================================================================
#
# Fit growth profiles per amino acid variant
#

# Calculate RPM with psudo_counts defined above
# psudo_counts = 0
print(sprintf("Normalizing amino acid counts to frequencies using psudo count %.3f",psudo_counts))
rpm_aa = data.frame(var=raw_aa$var)
cn = paste0("run",rep(unique(samples$date),each=4),"_day",rep(unique(samples$time),4))
rpm_aa = cbind(rpm_aa, reads_per_million(raw_aa[,cn] + psudo_counts))

# Growth profile calculated per AA variant
# Its a linear fit to (log) frequencies but not normalized to OD or similar. So if a barcode has a dieing phenotype but many others are dieing faster, if may appear enriched because 
growth_aa = data.frame(var=raw_aa$var, n_bc=raw_aa$n_bc)

# add number of aa substitutions
splitlist = strsplit(growth_aa$var, ":")
growth_aa$n_aa = sapply(splitlist, length)

# # add raw read counts for error estimation
# for (run in seq(4)) { cn=sprintf("run%d_sum",run); growth_aa[,cn] = raw_aa[,cn] }
# growth_aa$sumsum = apply(growth_aa[,sprintf("run%d_sum",seq(4))], MARGIN=1, sum)

# calculate line fit 
for (i in which(samples$time==0 & samples$tech_rep==1)) {
    cn_times = sprintf("run%d_day%d",samples[i,"date"],c(0,5,7,9))
    cn = sprintf("run%d", samples[i,"date"])
    
    # Normalize to sum of frequencies (distribution like)
    print(sprintf("Calculate distribution-style growth profiles per AA variant for %s",cn))
    row_sums = data.frame(rs = apply(rpm_aa[,cn_times], MARGIN=1, sum))
    df = rpm_aa[,cn_times] / row_sums[,rep(1,length(cn_times))]
    growth_aa = cbind(growth_aa, growth_profile(df, x=c(0,5,7,9), name=cn))
    
    # # Normalize to first timepoint, need to add psudo counts here
    # print(sprintf("Calculate t0 normalized growth profiles per AA variant for %s",cn))
    # df = rpm_aa[,cn_times] / rpm_aa[,rep(cn_times[1],length(cn_times))]
    # growth_aa = cbind(growth_aa, growth_profile(df, x=c(0,5,7,9), name=cn))    

    # # Enrich2 method with weighted least squares
    # print(sprintf("Calculate Enrich2-style growth profiles per AA variant for %s",cn))
    # df = log( (raw_aa[,cn_times]+0.5) / (raw_aa[rep(1,nrow(raw_aa)),cn_times]+0.5) )
    # wdf = raw_aa[,cn_times]+0.5
    # growth_aa = cbind(growth_aa, growth_profile(df, x=c(0,5,7,9), name=cn, wdf=wdf))
}

reads_threshold = 20

# Sum raw reads per growth curve and calc an unceratinty from this
for (run in unique(samples$date)) {
    cn_times = sprintf("run%d_day%d", run, c(0,5,7,9))
    cn_raw = sprintf("run%d_sum",run)
    raw_aa[,cn_raw] = apply(raw_aa[,cn_times], MARGIN=1, sum)
    growth_aa[,sprintf("run%d_reads",run)] = raw_aa[,cn_raw]

    # remove scores based on few counts
    i_low = which(growth_aa[,sprintf("run%d_reads",run)] < reads_threshold)
    growth_aa[i_low,sprintf("run%d_intersept",run)] = NA
    growth_aa[i_low,sprintf("run%d_slope",run)] = NA
    # growth_aa[i_low,sprintf("run%d_r2",run)] = NA
    growth_aa[i_low,sprintf("run%d_reads",run)] = 0
    print(sprintf("Removed %d slopes with fewer than %d reads from run %d", length(i_low), reads_threshold, run))

    cn_grow = sprintf("run%d_nu",run)
    growth_aa[,cn_grow] = 1 - 1/sqrt(growth_aa[,sprintf("run%d_reads",run)]/4)
    growth_aa[which(is.infinite(growth_aa[,cn_grow])), cn_grow] = NA
}
growth_aa$all_reads = apply(growth_aa[,sprintf("run%d_reads",seq(4))], MARGIN=1, sum)

i_grow0 = which(growth_aa$all_reads == 0)
raw_sums = apply(raw_aa[,sprintf("run%d_sum",seq(4))], MARGIN=1, sum)
i_raw0 = which(raw_sums == 0)
print(sprintf("AA variants without reads in any replica before filtering %d, after %d", length(i_raw0), length(i_grow0)))
i1 = which(growth_aa$n_aa == 1)
print(sprintf("Single AA variants without reads in any replica before filtering %d, after %d",
              length(intersect(i1,i_raw0)), length(intersect(i1,i_grow0))))
print(raw_aa[intersect(i1,i_grow0),])

# order variants
growth_aa$resi = sapply(splitlist, function(v){ if (length(v)>0) as.numeric(substr(v[1],2,nchar(v[1])-1)) else 0 })
growth_aa$mut = sapply(splitlist, function(v){ if (length(v)>0) substr(v[1],nchar(v[1]),nchar(v[1])) else "" })
growth_aa = growth_aa[order(growth_aa$n_aa, growth_aa$resi, growth_aa$mut, decreasing=F),]

# Calculate a mean of replica
cn = sprintf("run%d_slope",seq(4))
growth_aa$mean.slope = apply(growth_aa[,cn], MARGIN=1, mean, na.rm=T)
growth_aa$sd.slope = apply(growth_aa[,cn], MARGIN=1, sd, na.rm=T)
growth_aa$n_rep = apply(growth_aa[,cn], MARGIN=1, function(v){ sum(! is.na(v)) })

# Min-max normalization
# slope_tox = mean( growth_aa[c(which(growth_aa$var == "C152W"),which(growth_aa$var == "A305E")),"mean.slope"] )
# scale = 1/(growth_aa[1,"mean.slope"] - slope_tox)
# growth_aa$mean.slope.norm = (growth_aa$mean.slope - slope_tox) *scale
# growth_aa$sd.slope.norm = growth_aa$sd.slope *scale

# var_tox_ctrl = c("C152W")
slope_tox = mean( growth_aa[which(growth_aa$var %in% var_tox_ctrl),"mean.slope"] )
slope_wt = growth_aa[1,"mean.slope"] 
scale = 1/(slope_tox - slope_wt)
growth_aa$mean.slope.norm = (growth_aa$mean.slope - slope_wt) * scale
growth_aa$sd.slope.norm = growth_aa$sd.slope * abs(scale)
print(sprintf("Normalized growth slopes to WT = %.3f and toxic = %.3f", slope_wt, slope_tox))

# Tox-score correlations between biological replicates
i = which(growth_aa$n_aa == 1 & growth_aa$all_reads > 500)
print("Amino acid slope correlation of single AA variants between biological replicates (normalized scores)")
# print( cor( growth_aa[i,sprintf("run%d_slope",seq(4))], method="pearson", use="complete.obs"))
rps = c(); maes = c()
for ( run1 in seq(1,3) ) {
    x1 = (growth_aa[i,sprintf("run%d_slope",run1)]-slope_wt)*scale
    for ( run2 in seq(run1+1,4) ) {
        x2 = (growth_aa[i,sprintf("run%d_slope",run2)]-slope_wt)*scale
        rp = cor(x1, x2, method="pearson", use="complete.obs")
	mae = mean(abs(x1-x2), na.rm=T)
	rps = c(rps, rp)
	maes = c(maes, mae)
	print(sprintf("Run %d versus run %d: rp = %.3f, MAE = %.3f",run1,run2,rp,mae))
    }
}
tox_err = mean(maes)
print(sprintf("Average of %d pairs: rp %.3f (range %.3f-%.3f), MAE %.3f (range %.3f-%.3f)",
              length(rps),mean(rps),min(rps),max(rps),tox_err,min(maes),max(maes)))


# scatter plots of normalized scores
quartz(width=8, height=8)
panel = function(x,y,...) {
    points(x, y, ...)
    abline(0, 1)
    points(x[1], y[1], pch=20, col=2)
}
plot( (growth_aa[i,sprintf("run%d_slope",seq(4))] - slope_wt) *scale, labels=sprintf("Replica %d",seq(4)), panel=panel)
title(xlab="Toxicity score", line=4)
title(ylab="Toxicity score", line=3)

# Growth curves for controls - looks like read levels are indeed affected between day 0 and day 5
iplot = c(1, which(raw_aa$var=="C152W"), which(raw_aa$var=="A305E"))
# iplot = c(1, which(raw_aa$var=="C152W"), which(raw_aa$var=="E9I"), which(raw_aa$var=="Q12V"), which(raw_aa$var=="N37V"))
# iplot = c(1, which(raw_aa$var=="I191A"), which(raw_aa$var=="W252A"))
quartz(width=8, height=6)
plot(0, 0, col=0, xlim=c(0,9), ylim=c(0,0.8), xlab="Growth time [days]", ylab="Reads frequency [%]")
x = c(0,5,7,9)
for (ii in seq_along(iplot)) {
    i = iplot[ii]
    for (run in seq(4)) {
        y = raw_aa[i,sprintf("run%d_day%d",run,c(0,5,7,9))]/raw_aa[i,sprintf("run%d_sum",run)]
        # y = rpm_aa[i,sprintf("run%d_day%d",run,c(0,5,7,9))]
	# y = y/y[1,1] *100
	points(x, y, type="b", col=ii)
	abline(growth_aa[i,sprintf("run%d_intersept",run)], growth_aa[i,sprintf("run%d_slope",run)], lty=2, col=ii)
    }
}
legend("topright", c("WT", "C152W", "A305E", "Fit"), lty=c(1,1,1,2), lwd=2, col=c(seq_along(iplot),1))
# legend("topright", c("WT","C152W","E9I","Q12V","N37V"), lty=1, lwd=2, col=seq_along(iplot))
# legend("topright", c("WT","I191A","W252A"), lty=1, lwd=2, col=seq_along(iplot))
quartz.save("aspa_tox_ctrl_growth.png", type="png")


# ===========================================================================================================================================
#
#   Poorly estimated variants
#
# we dont use a threshold but give all error measures in the final file and the it is up the user to sort
# # quality cutoff on number of reads per AA variant across all replica
# reads_cut = 100


# poorly reproduced slopes are those with few reads - but few reads may also reproduce well
quartz(width=8, height=6)
plot(growth_aa$all_reads, growth_aa$sd.slope.norm, cex=0.6, log="x", xlab="Reads summed across all replica", ylab="Slope SD between replica")
i = which(growth_aa$n_aa==1)
points(growth_aa[i,"all_reads"], growth_aa[i,"sd.slope.norm"], pch=20, col=2, cex=0.6)
legend("topright", c("All AA variants","Single AA variants"), pch=c(1,20), col=c(1,2))
# abline(v=reads_cut, h=tox_err)
quartz.save("aspa_tox_reads_sd.png", type="png")


# # There is no correlation between the number of reads and the r2 of the line fit
# quartz(width=6, height=6)
# plot(1,1,col=0, xlim=c(1,10000), ylim=c(0,1), xlab="reads", ylab="r2", log="x", main="Single amino acid variants per replica")
# for (i in seq(4)) {
#     # print(summary(growth_aa[,sprintf("run%d_%s",i,c("r2","reads"))]))
#     points(growth_aa[,sprintf("run%d_%s",i,c("reads","r2"))], col=i, pch=".")
# }
# quartz.save("aspa_tox_reads_r2_cor.png", type="png")


# ===========================================================================================================================================
#
# A curated and nice data frame 
#
# make a nice data frame of single mutants with good read counts
i = which(growth_aa$n_aa < 2 )  # & growth_aa$sumsum > 100)
tox = data.frame(var  = growth_aa[i,"var"])
nc_var = nchar(tox$var)
tox$wt    = substr(tox$var, 1, 1)
tox$resi  = as.numeric(substr(tox$var, 2, nc_var-1))
tox$mut   = substr(tox$var, nc_var, nc_var)
tox$barcodes = growth_aa[i,"n_bc"]
tox$reads = growth_aa[i,"all_reads"]
tox$score = growth_aa[i,"mean.slope.norm"]
tox$sd    = growth_aa[i,"sd.slope.norm"]
tox[1,"resi"] = 0
tox[1,"var"] = "WT"
tox = tox[order(tox$resi, tox$mut),]

# coverage, which single substitutions are present
i1 = which(growth_aa$n_aa==1)   # & growth_aa$mut != "*")
agg = aggregate(growth_aa[i1,"mut"], list(growth_aa[i1,"resi"]), paste0, collapse="")

residue = data.frame(resi=seq(nchar(wt[["aa"]])), wt=strsplit(wt[["aa"]],"")[[1]])
residue$mut = ""
residue[agg[,1],"mut"] = agg[,2]
residue$coverage = sapply(residue$mut, nchar)
i_no_nons = which(growth_aa$mut != "*")
residue$median_tox = sapply(residue$resi, function(resi) { i = i1[which(growth_aa[i1,"resi"]==resi)];
                                                           median(growth_aa[intersect(i,i_no_nons),"mean.slope.norm"], na.rm=T) })
residue$tox_sd_prop = sapply(residue$resi, function(resi) { i = i1[which(growth_aa[i1,"resi"]==resi)];
                                                            sds = growth_aa[intersect(i,i_no_nons),"sd.slope.norm"];
		    					    sqrt(sum(sds^2)) })


# ===========================================================================================================================================
#
#   Aggregate per DNA variant
#
# Idea here is to make a growth curve per DNA variant and treat these as independent replica of an amino acid variant (in addition to the
#   biological replicates)

print("Aggregate barcodes per DNA variant")
agg = aggregate(seq(nrow(raw_bc)), list(raw_bc$var_dna), paste, collapse=":")
dna_indices = lapply( strsplit(agg[,2], ":"), as.numeric)
names(dna_indices) = agg[,1]

raw_dna = data.frame(var_dna=names(dna_indices), var_aa=raw_bc[sapply(dna_indices,"[[",1),"var_aa"], n_bc=sapply(dna_indices,length))
dna_split_list = strsplit(raw_dna$var_dna, ":")
raw_dna$n_dna = sapply(dna_split_list, length)
dna_split_list = strsplit(raw_dna$var_dna, ":")
raw_dna$n_dna = sapply(dna_split_list, length)
aa_split_list = strsplit(raw_dna$var_aa, ":")
raw_dna$n_aa = sapply(aa_split_list, length)
for (i in which(samples$tech_rep==1)) {
    cn = sprintf("run%d_day%d",samples[i,"date"],samples[i,"time"])
    raw_dna[,cn] = sapply(dna_indices, function(v) { sum(raw_bc[v,cn], na.rm=T) })
}
# Total number of reads per replica
for (run in unique(samples$date)) {
    cns = sprintf("run%d_day%d",run,c(0,5,7,9))
    cn = sprintf("run%1d_sum",run)
    raw_dna[,cn] = apply(raw_dna[,cns], MARGIN=1, sum)
}
# check for NA counts
stopifnot(sum(is.na(raw_dna[,sprintf("run%d_sum",unique(samples$date))])) == 0)

# # histogram of counts per dna variant and run, summed over days in growth tracejtory
# quartz(width=8, height=6)
# hist(unlist(raw_dna[,sprintf("run%d_sum",unique(samples$date))]), breaks=seq(0,500000,10), xlim=c(0,3000), xlab="Total read counts per DNA variant and biological replica")

for (day in unique(samples$time)) {
    print(sprintf("Replica correlations of raw counts per DNA variant at day %d",day))
    print(cor(raw_dna[,paste("run",seq(4),"_day",5, sep="")]))
}

cn = paste0("run",rep(unique(samples$date),each=4),"_day",rep(unique(samples$time),4))
sample_coverage = apply(raw_bc[,cn], MARGIN=2, sum, na.rm=T) / nrow(raw_bc)
sample_coverage = apply(raw_dna[,cn], MARGIN=2, sum, na.rm=T) / nrow(raw_dna)
print(sprintf("Per DNA variant sequencing coverage average x%.2f (range %.2f-%.2f) [reads per variant]",
              mean(sample_coverage), min(sample_coverage), max(sample_coverage)))
print(sample_coverage)

# Calculate RPM with psudo_counts defined above
# psudo_counts = 0
print(sprintf("Normalizing DNA variant counts to frequencies using psudo count %.3f",psudo_counts))
rpm_dna = data.frame(var_dna=raw_dna$var_dna)
cn = paste0("run",rep(unique(samples$date),each=4),"_day",rep(unique(samples$time),4))
rpm_dna = cbind(rpm_dna, reads_per_million(raw_dna[,cn] + psudo_counts))

# Growth profile calculated per DNA variant
growth_dna = data.frame(var_dna=raw_dna$var_dna, var_aa=raw_dna$var_aa, n_bc=raw_dna$n_bc, n_dna=raw_dna$n_dna, n_aa=raw_dna$n_aa)
for (i in which(samples$time==0 & samples$tech_rep==1)) {
    cn_times = sprintf("run%d_day%d",samples[i,"date"],c(0,5,7,9))
    cn = sprintf("run%d", samples[i,"date"])
    
    print(sprintf("Calculate distribution-style growth profiles per DNA variant for %s",cn))
    row_sums = data.frame(rs = apply(rpm_dna[,cn_times], MARGIN=1, sum))
    df = rpm_dna[,cn_times] / row_sums[,rep(1,length(cn_times))]
    growth_dna = cbind(growth_dna, growth_profile(df, x=c(0,5,7,9), name=cn))    

    # print(sprintf("Calculate t0 normalized growth profiles per DNA variant for %s",cn))
    # df = rpm_dna[,cn_times] / rpm_dna[,rep(cn_times[1],length(cn_times))]
    # growth_dna = cbind(growth_dna, growth_profile(df, x=c(0,5,7,9), name=cn))
    
    # # Enrich2 method
    # print(sprintf("Calculate Enrich2-style growth profiles per DNA variant for %s",cn))
    # df = log( (raw_dna[,cn_times]+0.5) / (raw_dna[rep(1,nrow(raw_dna)),cn_times]+0.5) )
    # wdf = raw_dna[,cn_times]+0.5
    # growth_dna = cbind(growth_dna, growth_profile(df, x=c(0,5,7,9), name=cn, wdf=wdf))

}

# Sum raw reads per growth curve and calc an unceratinty from this
for (run in unique(samples$date)) {
    cn_times = sprintf("run%d_day%d", run, c(0,5,7,9))
    cn_raw = sprintf("run%d_sum",run)
    raw_dna[,cn_raw] = apply(raw_dna[,cn_times], MARGIN=1, sum)
    growth_dna[,sprintf("run%d_reads",run)] = raw_dna[,cn_raw]
    cn_grow = sprintf("run%d_nu",run)
    growth_dna[,cn_grow] = 1 - 1/sqrt(raw_dna[,cn_raw]/4)
    growth_dna[which(is.infinite(growth_dna[,cn_grow])), cn_grow] = NA
}
cns = sprintf("run%d_day%d",rep(seq(4),each=4),rep(c(0,5,7,9),4))
growth_dna$all_reads = apply(raw_dna[,cns], MARGIN=1, sum)
stopifnot( all( growth_dna$all_sum == apply(growth_dna[,sprintf("run%d_reads",seq(4))], MARGIN=1, sum) ) ) 

# order variants
growth_dna$resi = sapply(aa_split_list, function(v){ if (length(v)>0) as.numeric(substr(v[1],2,nchar(v[1])-1)) else 0 })
growth_dna$mut = sapply(aa_split_list, function(v){ if (length(v)>0) substr(v[1],nchar(v[1]),nchar(v[1])) else "" })
growth_dna = growth_dna[order(growth_dna$n_aa, growth_dna$resi, growth_dna$mut, growth_dna$n_dna, decreasing=F),]

# Calculate a mean of replica
cn = sprintf("run%d_slope",seq(4))
growth_dna$mean.slope = apply(growth_dna[,cn], MARGIN=1, mean, na.rm=T)
growth_dna$sd.slope = apply(growth_dna[,cn], MARGIN=1, sd, na.rm=T)

# Min-max normalization based on DNA variants
# dna_slope_tox = mean( growth_dna[c(which(growth_dna$var_aa == "C152W"),which(growth_dna$var_aa == "A305E")),"mean.slope"] )
# dna_scale = 1/(growth_dna[1,"mean.slope"] - dna_slope_tox)
# growth_dna$mean.slope.norm = (growth_dna$mean.slope - dna_slope_tox) *dna_scale
# growth_dna$sd.slope.norm = growth_dna$sd.slope *dna_scale
dna_slope_tox = mean( growth_dna[which(growth_dna$var_aa %in% var_tox_ctrl),"mean.slope"] )
dna_slope_wt = growth_dna[1,"mean.slope"] 
dna_scale = 1/(dna_slope_tox - dna_slope_wt)
growth_dna$mean.slope.norm = (growth_dna$mean.slope - dna_slope_wt) * dna_scale
growth_dna$sd.slope.norm = growth_dna$sd.slope * abs(dna_scale)
print(sprintf("Normalized growth slope of DNA variants to WT = %.3f and toxic = %.3f", dna_slope_wt, dna_slope_tox))
print(sprintf("Difference in normalization DNA-AA:  wt %.5f  and  tox %.5f", dna_slope_wt-slope_wt, dna_slope_tox-slope_tox))

# Tox-score correlations between biological replicates
i = which(growth_dna$n_aa == 1 )  # & growth_dna$sumsum > 100)
print("DNA slope correlation between biological replicates (DNA normalized scores)")
rps = c(); maes = c()
for ( run1 in seq(1,3) ) {
    x1 = (growth_dna[i,sprintf("run%d_slope",run1)]-dna_slope_wt)*dna_scale
    for ( run2 in seq(run1+1,4) ) {
        x2 = (growth_dna[i,sprintf("run%d_slope",run2)]-dna_slope_wt)*dna_scale
        rp = cor(x1, x2, method="pearson", use="complete.obs")
	mae = mean(abs(x1-x2), na.rm=T)
	rps = c(rps, rp)
	maes = c(maes, mae)
	print(sprintf("Run %d versus run %d: rp = %.3f, MAE = %.3f",run1,run2,rp,mae))
    }
}
print(sprintf("Average of %d pairs: rp %.3f (range %.3f-%.3f), MAE %.3f (range %.3f-%.3f)",
              length(rps),mean(rps),min(rps),max(rps),mean(maes),min(maes),max(maes)))

# Is there a difference in the distribution of normalized slopes between AA and DNA variants
quartz(width=8, height=6)
breaks = seq(-10,90,0.1)
i1_aa = which(growth_aa$n_aa == 1)
h_aa = hist(growth_aa[i1_aa,"mean.slope.norm"], breaks=breaks, plot=F)
i1_dna = which(growth_dna$n_aa == 1)
h_dna = hist(growth_dna[i1_dna,"mean.slope.norm"], breaks=breaks, plot=F)
plot(0, 0, col=0, xlim=c(-1,2), ylim=c(0,3), xlab="Normalized slope", ylab="Density", main="Only variants with one AA substitution")
lines(h_aa$mids, h_aa$density, lwd=2, col=2)
lines(h_dna$mids, h_dna$density, lwd=2, col=3)
legend("topleft", c(sprintf("AA variants %d",length(i1_aa)),sprintf("DNA variants %d",length(i1_dna))), lty=1, lwd=2, col=c(2,3))

# make a data frame of synonymous WT (and WT) to be compared to AA variants - growth_aa and growth_dna are normalized in the same way but possible with different numbers
i_wt_syn = which(growth_dna$var_aa=="" & growth_dna$all_reads > 0)
tox_wt_syn = data.frame(var=growth_dna[i_wt_syn,"var_dna"], n_bc=growth_dna[i_wt_syn,"n_bc"], all_reads=growth_dna[i_wt_syn,"all_reads"])
# Use normalization of AA variants here
tox_wt_syn$score = (growth_dna[i_wt_syn,"mean.slope"] - slope_wt) *scale
tox_wt_syn$sd = growth_dna[i_wt_syn,"sd.slope"] * abs(scale)
stopifnot(tox_wt_syn[1,"var"] == "")
tox_wt_syn[1,"var"] = "WT"

# Signal, noise and coverage for different threshold of raw read counts
i_wt = which(rpm_dna$var_dna=="")
slopes_wt = unlist(growth_dna[i_wt,sprintf("run%d_slope",seq(4))])
i1 = which(raw_dna$n_aa == 1)
print(sprintf("WT     slope %8.5f +- %.5f coverage %d AA variant", mean(slopes_wt), sd(slopes_wt), length(unique(growth_dna[i1,"var_aa"]))))

# To kontroller til ASPA tox assay: C152W og A305E
for (var_tox in c("C152W","A305E")) {
    i_tox = which(raw_dna$var_aa == var_tox)
    slopes_tox = unlist(growth_dna[i_tox,sprintf("run%d_slope",seq(4))])
    run_reads = raw_dna[i_tox,sprintf("run%d_sum",seq(4))]
    print(sprintf("%5s  slope %8.5f +- %.5f with %d read counts (per run counts %s)", var_tox, mean(slopes_tox), sd(slopes_tox), sum(run_reads), paste(run_reads, collapse=", ")))
}

for (raw_sumsum_threshold in c(0,10,20,50,100,200,300,400,500,1000,2000,5000)) {
    raw_sumsum = apply(raw_dna[,sprintf("run%d_sum",seq(4))], MARGIN=1, sum)
    i_syn = which(raw_dna$n_dna>0 & raw_dna$var_aa=="" & raw_sumsum>raw_sumsum_threshold )
    slopes_syn = unlist(growth_dna[i_syn,sprintf("run%d_slope",seq(4))])

    i = which(raw_dna$n_aa == 1 & raw_sumsum > raw_sumsum_threshold)
    print(sprintf("WT syn slope %8.5f +- %.5f coverage %d AA variants using threshold %d",
                  mean(slopes_syn, na.rm=T), sd(slopes_syn, na.rm=T), length(unique(growth_dna[i,"var_aa"])), raw_sumsum_threshold))
}
# a threshold on sumsum of 100 seems good to give an unceratinty of 0.01 and good coverage

# # Plot growth profile for WT and syn WT
# quartz(width=8, height=6)
# plot(0, 0, col=0, xlim=c(0,9), ylim=c(0,1), xlab="Growth time [days]", ylab="Relative populations", main="Normalized RPM growth curves for DNA variants")
# i = which(rpm_dna$var_dna=="")
# for (run in seq(4)) {
#     x = c(0,5,7,9)
#     y = unlist( rpm_dna[i,sprintf("run%d_day%d",run,c(0,5,7,9))] )
#     lines(x, y/sum(y), col=run, lwd=3)
#     # lines(x, y/y[1], col=run, lwd=3)
# }
# raw_sum_threshold = 100
# i_syn = which(raw_dna$n_dna>0 & raw_dna$var_aa=="")
# for (i in i_syn) {
#     for (run in seq(4)) {
#         if (raw_dna[i,sprintf("run%d_sum",run)] > raw_sum_threshold) {
#             x = c(0,5,7,9)
#             y = unlist(rpm_dna[i,sprintf("run%d_day%d",run,c(0,5,7,9))])
#             lines(x, y/sum(y), col=run, lwd=1)
# 	    # if (y[1] > 0) {
#             #     lines(x, y/y[1], col=run, lwd=1)
# 	    # } else {
# 	    #     print(sprintf("Skip %s run %d with %d reads but %d at time 0", rpm_dna[i,"var_dna"], run, raw_dna[i,sprintf("run%d_sum",run)], y[1]))
# 	    # }
# 	} else {
# 	    print(sprintf("Skip %s run %d with %3d reads < %d", rpm_dna[i,"var_dna"], run, raw_dna[i,sprintf("run%d_sum",run)], raw_sum_threshold))
# 	}
#     }
# }


# ===========================================================================================================================================
#
# Dump 
#
# save everything in R data file
save(tox, tox_err, tox_wt_syn, growth_aa, growth_dna, growth_bc, raw_aa, raw_dna, raw_bc, samples, residue, wt,
     slope_tox, slope_wt, scale, dna_slope_tox, dna_slope_wt, dna_scale, file="aspa_tox.rda")

# dump tox scores as csv
write.table(tox[,c("var","barcodes","reads","score","sd")], col.names=c("variant","barcodes","reads","tox_score","tox_std"),
            file="aspa_toxicity.csv", sep=";", row.names=F, quote=F)

# dump tox residue frame
write.table(residue, file="aspa_tox_residues.csv", sep=";", row.names=F, quote=F)

# Dump synonymous WT
write.table(data.frame(var=tox_wt_syn$var, barcodes=tox_wt_syn$n_bc, tox_score=tox_wt_syn$score, tox_sd=tox_wt_syn$sd),
            file="aspa_tox_wt_syn.csv", sep=";", row.names=F, quote=F)

# dump prism file
f = file("prism_mave_161_ASPA_tox.txt", "wt")
write("# --------------------", f)
write("# version: 1", f)
write("# protein:", f)
write("#     name: ASPA", f)
write(sprintf("#     sequence: %s",wt[["aa"]]), f)
write("#     organism: Homo sapiens (Human)", f)
write("#     uniprot: P45381", f)
write("# mave:", f)
write("#     organism: Homo sapiens (Human)", f)
write("#     cloning: Landing pad", f)
write("#     expression: Overexpression", f)
write("#     technology: growth", f)
write("#     doi: unpublished", f)
write("#     year: 2022", f)
write("# variants:", f)
write(sprintf("#     number: %d",nrow(tox)), f)
write(sprintf("#     coverage: %.2f",mean(residue$coverage > 0)), f)
write(sprintf("#     depth: %.2f",mean(nchar(gsub("*", "", residue[which(residue$coverage > 0),"mut"], fixed=T)))), f)
write("#     width: single mutants", f)
write("# columns:", f)
write("#     barcodes: Number of different barcodes observed for variant", f)
write("#     reads: Number of matched illumina reads observed across all replica and barcodes", f)
write("#     tox_score: Toxicity score is the slope of a growth curve normalized to wild-type and the toxic variant C152W", f)
write("#     tox_std: Toxicity score standard deviation between 4 biological replicates", f)
write("# --------------------", f)
write("#", f)
write("# Unpublished version of Sep 2022 - kristoffer.johansson@bio.ku.dk", f)
write("# ", f)
write(sprintf("%5s  %8s  %8s  %10s  %8s", "var", "barcodes", "reads", "tox_score", "tox_std"), f)
write(sprintf("%5s  %8d  %8d  %10.4f  %8.4f",tox$var, tox$barcodes, tox$reads, tox$score, tox$sd), f)
close(f)


# Analysis
breaks = seq(-1,2,0.05)
ia_a1 = which(growth_aa$n_aa == 1)
h1 = hist(growth_aa[ia_a1,"mean.slope.norm"], breaks=breaks, plot=F)

# id_wt_syn = which(growth_dna$n_dna > 0 & growth_dna$n_aa == 0)
# hsyn = hist(growth_dna[id_wt_syn,"mean.slope.norm"], breaks=breaks, plot=F)
hsyn = hist(tox_wt_syn[2:nrow(tox_wt_syn),"score"], breaks=breaks, plot=F)

plot(0, 0, col=0, xlim=c(-0.5, 1.5), xlab="Toxicity score", ylim=c(0,9), ylab="Density")
lines(h1$mids, h1$density, col=8, lwd=2)
lines(hsyn$mids, hsyn$density, col=2, lwd=2)

# # plot WT from both AA and DNA data frames
# points(c(growth_aa[1,"mean.slope.norm"], growth_dna[1,"mean.slope.norm"]), c(0,0.5), pch="|", cex=2, col=c(8,2))

# plot toxic controls from AA data frame
ia = c(1, which(growth_aa$var %in% c("C152W","A305E")))
points(growth_aa[ia,"mean.slope.norm"], rep(0,length(ia)), pch="|", cex=2, col=8)

# plot error bars on all controls
arrows(growth_aa[ia,"mean.slope.norm"]-growth_aa[ia,"sd.slope.norm"], 0.03*seq(0,length(ia)-1), x1=growth_aa[ia,"mean.slope.norm"]+growth_aa[ia,"sd.slope.norm"],
       code=3, angle=90, length=.1, col=8)

# plot toxic controls from DNA data frame
id = c(1, which(growth_dna$var_aa %in% c("C152W","A305E")))
points(growth_dna[id,"mean.slope.norm"], rep(0.5,length(id)), pch="|", cex=2, col=2)

arrows(growth_dna[id,"mean.slope.norm"]-growth_dna[id,"sd.slope.norm"], 0.5+0.03*seq(0,length(id)-1), x1=growth_dna[id,"mean.slope.norm"]+growth_dna[id,"sd.slope.norm"],
       code=3, angle=90, length=.1, col=2)

# plot an errorbar of replica MAE
arrows(growth_aa[1,"mean.slope.norm"]-tox_err, 0.25, x1=growth_aa[1,"mean.slope.norm"]+tox_err,
       code=3, angle=90, length=.1, col=1)

legend("topright", c("Single amino acid variants", "Synonymous WT",
                    sprintf("Replica MAE of all AA variants: %.3f",tox_err),
                    sprintf("WT (AA) replica SD: %.3f",growth_aa[1,"sd.slope.norm"]),
		    sprintf("WT (DNA) replica SD: %.3f", growth_dna[1,"sd.slope.norm"])),
       lty=c(1,1,NA,NA,NA), pch=c(NA,NA,"|","|","|"), col=c(8,2,1,8,2))

quartz.save("aspa_tox_distribution.png", type="png")



