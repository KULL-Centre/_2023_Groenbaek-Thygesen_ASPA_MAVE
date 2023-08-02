options(width=180)

load("raw_vamp.rda")

isamples  = seq(33,120)
dates     = rep(seq(4), c(16,24,24,24))               # Date 1, FACS replica 3 failed
facs_rep  = c(rep(seq(2), each=8), rep(rep(seq(3), each=8), 3))
gates     = rep(rep(c("A","B","C","D"), each=2), 11)
tech_rep  = rep(seq(2), 44)

samples = data.frame(date=dates, gate=gates, facs_rep=facs_rep, tech_rep=tech_rep, sample_id=sprintf("%d_S%d",isamples,isamples))

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
    print(sprintf("Lane correlations of %s: %s", samples[i,"sample_id"], paste(rp,collapse=", ")))
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
raw_rm = data.frame(barcode=raw$barcode, var_dna=raw$var_dna, var_aa=raw$var_aa)
subst_list_aa = strsplit(raw_rm$var_aa, ":")
raw_rm$n_aa = sapply(subst_list_aa, length)
subst_list_dna = strsplit(raw_rm$var_dna, ":")
raw_rm$n_dna = sapply(subst_list_dna, length)
for (i in which(samples$tech_rep==1)) {
    # sample row of other tech_rep 
    i2 = which(samples$date     == samples[i,"date"]      &  samples$gate     == samples[i,"gate"] &
               samples$facs_rep == samples[i,"facs_rep"]  &  samples$tech_rep == 2)
    stopifnot(length(i2) == 1)
    cn_rep = c(samples[i,"sample_id"], samples[i2,"sample_id"])
    cn = sprintf("%s-%s%d", samples[i,"date"], samples[i,"gate"], samples[i,"facs_rep"])
    raw_rm[,cn] = apply(raw_lm[,cn_rep], MARGIN=1, sum)

    # Correlations
    rp_raw = cor(raw_lm[,cn_rep[1]], raw_lm[,cn_rep[2]], method="pearson")
    rep_rp_raw = c(rep_rp_raw, rp_raw)
    print(sprintf("Pearson between technical replica %s and %s is %.2f",cn_rep[1],cn_rep[2],rp_raw))
    
    # Read counts per replica
    nreads = apply(raw_lm[,cn_rep], MARGIN=2, sum)
    rep_nreads = c(rep_nreads, nreads)
    print(sprintf("Number of reads for technical replica: %s", paste(nreads, collapse=", ")))
}
print(sprintf("Replica raw correlations %.2f sd %.2f min %.2f max %.2f",mean(rep_rp_raw),sd(rep_rp_raw),min(rep_rp_raw),max(rep_rp_raw)))
print(sprintf("Average reads per replica %.0f sd %.0f min %d max %d",mean(rep_nreads),sd(rep_nreads),min(rep_nreads),max(rep_nreads)))

raw_bc = raw_rm

# Does reads-per-barcode correlate with replica correlation?
# cns = colnames(raw_bc[,3:ncol(raw_bc)])
is = which(samples$tech_rep==1)
cns = paste0(samples[is,"date"],"-",samples[is,"gate"],samples[is,"facs_rep"])
names(rep_rp_raw) = cns
read_stats = data.frame(name=cns, reads=apply(raw_bc[,cns], MARGIN=2, sum), rp_raw=rep_rp_raw)
read_stats$nbc_nz = apply(raw_bc[,cns], MARGIN=2, function(v){ sum(v>0) })
read_stats$nbc_eff = apply(raw_bc[,cns], MARGIN=2, function(v){ p=v/sum(v); exp(-sum(p*log(p+1e-12))) })
read_stats$avg_rpb = read_stats$reads / read_stats$nbc_eff
quartz(width=8, height=8)
plot(read_stats[,c("rp_raw","nbc_nz","nbc_eff","avg_rpb")])
quartz.save("vamp_complexity_and_reads.png", type="png")


# Function that normalize all columns to RPM
reads_per_million = function(df) {
    col_sums = apply(df, MARGIN=2, sum)
    df * matrix(10^6/col_sums, nrow=dim(df)[1], ncol=dim(df)[2], byrow=TRUE)
}

# ===========================================================================================================================================
#
# PSI per barcode
#
# Calculate count frequencies with psudo_counts
psudo_counts = 0
rpm_bc = data.frame(barcode=raw_bc$barcode, var=raw_bc$var_aa)
is = which(samples$tech_rep==1)
cns = paste0(samples[is,"date"],"-",samples[is,"gate"],samples[is,"facs_rep"])
rpm_bc = cbind(rpm_bc, reads_per_million(raw_bc[,cns] + psudo_counts))

# Function to calculate PSI based on all columns of a data frame
protein_stability_index = function(df, name) {
    # First gate (gate 1) is stable with index 4, and the late (gate 4) unstable with index 1
    psi = t(apply(df, MARGIN=1, function(v){vsum=sum(v); c(vsum, sum( v/vsum * c(1,2,3,4) ))}))
    # put in a data frame and give column name
    ret_df = data.frame(psi)
    colnames(ret_df) = paste0(name,c("_rpm_sum","_psi"))
    return(ret_df)
}

# PSI calculated per barcode
psi_bc = data.frame(barcode=raw_bc$barcode, var_aa=raw_bc$var_aa, var_dna=raw_bc$var_dna)
for (i in which(samples$gate=="A" & samples$tech_rep==1)) {
    cn_gates = sprintf("%s-%s%d", samples[i,"date"], c("A","B","C","D"), samples[i,"facs_rep"])
    cn = sprintf("run%s_facs%d",samples[i,"date"],samples[i,"facs_rep"])
    print(sprintf("Calculate PSI for %s",cn))
    psi_bc = cbind(psi_bc, protein_stability_index(rpm_bc[,cn_gates], cn))
}


# ===========================================================================================================================================
#
# Aggregate raw counts per AA variant
#
print("Aggregate barcodes per AA variant")
agg = aggregate(seq(nrow(raw_bc)), list(raw_bc$var_aa), paste, collapse=":")
var_indices = lapply( strsplit(agg[,2], ":"), as.numeric)
names(var_indices) = agg[,1]

raw_aa = data.frame(var=names(var_indices), n_bc=sapply(var_indices,length))
subst_list = strsplit(raw_aa$var, ":")
raw_aa$n_subst = sapply(subst_list, length)
for (is in which(samples$tech_rep==1)) {
    cns = sprintf("%d-%s%d",samples[is,"date"],samples[is,"gate"],samples[is,"facs_rep"])
    cn = sprintf("run%s_facs%d_gate%s",samples[is,"date"],samples[is,"facs_rep"],samples[is,"gate"])
    print(sprintf("Aggregate %s",cn))
    raw_aa[,cn] = sapply(var_indices, function(v) { sum(raw_bc[v,cns], na.rm=T) })
}

# Total number of reads per replica
for (i in which(samples$tech_rep==1 & samples$gate=="A")) {
    cns = sprintf("run%d_facs%d_gate%s",samples[i,"date"],samples[i,"facs_rep"],c("A","B","C","D"))
    cn = sprintf("run%d_facs%d_sum",samples[i,"date"],samples[i,"facs_rep"])
    raw_aa[,cn] = apply(raw_aa[,cns], MARGIN=1, sum)
}

# order variants
raw_aa$resi = sapply(subst_list, function(v){ if (length(v)>0) as.numeric(substr(v[1],2,nchar(v[1])-1)) else 0 })
raw_aa$mut = sapply(subst_list, function(v){ if (length(v)>0) substr(v[1],nchar(v[1]),nchar(v[1])) else "" })
raw_aa = raw_aa[order(raw_aa$n_subst, raw_aa$resi, raw_aa$mut, decreasing=F),]

# calculate correlations of replica raw reads per gate 
rps_bc = c()
rps_aa = c()
for (gate in c("A","B","C","D")) {
    print(sprintf("Replica correlations of raw counts per amino acid variant at gate %s",gate))
    
    # exclude WT with very high read counts and/or log counts to avoid inflated correlations
    is = which(samples$tech_rep==1 & samples$gate==gate)
    
    cns_bc = sprintf("%d-%s%d",samples[is,"date"],gate,samples[is,"facs_rep"])
    i_bc = 2:nrow(raw_bc)
    # i_bc = which(raw_bc$n_aa == 1)
    df_bc = raw_bc[i_bc,cns_bc]
    
    cns_aa = sprintf("run%d_facs%d_gate%s",samples[is,"date"],samples[is,"facs_rep"],gate)
    i_aa = 2:nrow(raw_aa)
    # i_aa = which(raw_aa$n_subst == 1)
    df_aa = raw_aa[i_aa,cns_aa]
    
    # # plot to look for jackpotting
    # quartz(width=8, height=8)
    # plot(df)
    # # quartz.save(sprintf("aspa_vamp_raw_cor_gate%s.png",gate), type="png")
    
    # correlations
    cm_bc = cor(df_bc, method="pearson")
    rps_bc = c(rps_bc, cm_bc[upper.tri(cm_bc)])
    
    cm_aa = cor(df_aa, method="pearson")
    rps_aa = c(rps_aa, cm_aa[upper.tri(cm_aa)])
    
    print(sprintf("Average Pearson correlation for barcodes gate %s %.4f range %.4f - %.4f", gate,
                  mean(cm_bc[upper.tri(cm_bc)]), min(cm_bc[upper.tri(cm_bc)]), max(cm_bc[upper.tri(cm_bc)])))
    print(sprintf("Average Pearson correlation for AA var   gate %s %.4f range %.4f - %.4f", gate,
                  mean(cm_aa[upper.tri(cm_aa)]), min(cm_aa[upper.tri(cm_aa)]), max(cm_aa[upper.tri(cm_aa)])))
}
print(sprintf("Average of %d Pearson correlations for barcodes    %.4f range %.4f - %.4f",length(rps_bc),mean(rps_bc),min(rps_bc),max(rps_bc)))
print(sprintf("Average of %d Pearson correlations for AA variants %.4f range %.4f - %.4f",length(rps_aa),mean(rps_aa),min(rps_aa),max(rps_aa)))


# # Sequencing coverage
# cn = paste0("run",rep(unique(samples$date),each=4),"_day",rep(unique(samples$time),4))
# sample_coverage = apply(raw_bc[,cn], MARGIN=2, sum, na.rm=T) / nrow(raw_bc)
# print(sprintf("Per-barcode sequencing coverage average x%.2f (range %.2f-%.2f) [reads per barcode]",
#               mean(sample_coverage), min(sample_coverage), max(sample_coverage)))
# print(sample_coverage)

# cn = paste0("run",rep(unique(samples$date),each=4),"_day",rep(unique(samples$time),4))
# sample_coverage = apply(raw_aa[,cn], MARGIN=2, sum, na.rm=T) / nrow(raw_aa)
# print(sprintf("Per AA variant sequencing coverage average x%.2f (range %.2f-%.2f) [reads per variant]",
#               mean(sample_coverage), min(sample_coverage), max(sample_coverage)))
# print(sample_coverage)


# # Number of reads
# raw_aa_sum = apply(raw_aa[,sprintf("run%d_day%d",rep(seq(4),each=4),rep(c(0,5,7,9)))], MARGIN=2, sum)
# print(sprintf("Average number of matched reads per AA variant %.0f range %d - %d", mean(raw_aa_sum), min(raw_aa_sum), max(raw_aa_sum)))


# ===========================================================================================================================================
#
# PSI per amino acid variant
#
# Calculate count frequencies with psudo_counts
psudo_counts = 0
rpm_aa = data.frame(var=raw_aa$var)
is = which(samples$tech_rep==1)
cns = sprintf("run%d_facs%d_gate%s",samples[is,"date"],samples[is,"facs_rep"],samples[is,"gate"])
rpm_aa = cbind(rpm_aa, reads_per_million(raw_aa[,cns] + psudo_counts))

# PSI calculated per variant
psi_aa = data.frame(var=raw_aa$var, n_subst=raw_aa$n_subst, resi=raw_aa$resi, mut=raw_aa$mut, n_bc=raw_aa$n_bc)
for (is in which(samples$gate=="A" & samples$tech_rep==1)) {
    cn_gates = sprintf("run%d_facs%d_gate%s",samples[is,"date"],samples[is,"facs_rep"],c("A","B","C","D"))
    cn = sprintf("run%s_facs%d",samples[is,"date"],samples[is,"facs_rep"])
    print(sprintf("Calculate PSI for %s",cn))
    psi_aa = cbind(psi_aa, protein_stability_index(rpm_aa[,cn_gates], cn))
}

reads_threshold = 20

# Total number of reads per replica
for (is in which(samples$tech_rep==1 & samples$gate=="A")) {
    run = samples[is,"date"]
    facs = samples[is,"facs_rep"]
    raw_reads_cn = sprintf("run%d_facs%d_sum",run,facs)
    psi_reads_cn = sprintf("run%d_facs%d_reads",run,facs)
    psi_aa[,psi_reads_cn] = raw_aa[,raw_reads_cn]

    # remove scores based on few counts
    i_low = which(psi_aa[,psi_reads_cn] < reads_threshold)
    psi_aa[i_low,sprintf("run%d_facs%d_psi",run,facs)] = NA
    psi_aa[i_low,sprintf("run%d_facs%d_rpm_sum",run,facs)] = NA
    psi_aa[i_low,psi_reads_cn] = 0
    print(sprintf("Removed %d PSI's with fewer than %d reads from run %d facs %d",
                  length(i_low), reads_threshold, run, facs))

    psi_nu_cn = sprintf("run%d_facs%d_nu",run,facs)
    psi_aa[,psi_nu_cn] = 1 - 1/sqrt(raw_aa[,raw_reads_cn]/4)
    psi_aa[which(is.infinite(psi_aa[,psi_nu_cn])), psi_nu_cn] = NA
}
is = which(samples$gate=="A" & samples$tech_rep==1)
cns_psi = sprintf("run%s_facs%d_reads",samples[is,"date"],samples[is,"facs_rep"])
psi_aa$all_reads = apply(psi_aa[,cns_psi], MARGIN=1, sum)

i_psi0 = which(psi_aa$all_reads == 0)
cns_raw = sprintf("run%s_facs%d_sum",samples[is,"date"],samples[is,"facs_rep"])
raw_sums = apply(raw_aa[,cns_raw], MARGIN=1, sum)
i_raw0 = which(raw_sums == 0)
print(sprintf("AA variants without reads in any replica before filtering %d, after %d", length(i_raw0), length(i_psi0)))
i1 = which(psi_aa$n_subst == 1)
print(sprintf("Single AA variants without reads in any replica before filtering %d, after %d",
              length(intersect(i1,i_raw0)), length(intersect(i1,i_psi0))))
# print(raw_aa[intersect(i1,i_psi0),])

# Calculate a mean of replica
is = which(samples$gate=="A" & samples$tech_rep==1)
cns = sprintf("run%s_facs%d_psi",samples[is,"date"],samples[is,"facs_rep"])
psi_aa$mean_psi = apply(psi_aa[,cns], MARGIN=1, mean, na.rm=T)
psi_aa$sd_psi = apply(psi_aa[,cns], MARGIN=1, sd, na.rm=T)

# Min-max normalization
i_nons = which(psi_aa$mut == "*" & psi_aa$n_subst == 1)
psi_nonsense = median(psi_aa[i_nons,"mean_psi"], na.rm=T)
scale = 1/(psi_aa[1,"mean_psi"] - psi_nonsense)
psi_aa$mean_psi_norm = (psi_aa$mean_psi - psi_nonsense) * scale
psi_aa$sd_psi_norm = psi_aa$sd_psi * scale
print(sprintf("Normalized PSI to WT = %.3f and nonsense = %.3f",psi_aa[1,"mean_psi"],psi_nonsense))


# ===========================================================================================================================================
#
# Aggregate PSI per AA variant
#
print("Aggregate barcodes per variant")
agg = aggregate(seq(nrow(psi_bc)), list(psi_bc$var_aa), paste, collapse=":")
var_indices = lapply( strsplit(agg[,2], ":"), as.numeric)
names(var_indices) = agg[,1]

psi_bc_aa = data.frame(var=names(var_indices), n_bc=sapply(var_indices,length))
subst_list = strsplit(psi_bc_aa$var, ":")
psi_bc_aa$n_subst = sapply(subst_list, length)
for (i in which(samples$gate=="A" & samples$tech_rep==1)) {
    cn = sprintf("run%s_facs%d",samples[i,"date"],samples[i,"facs_rep"])
    print(sprintf("Aggregate %s",cn))
    psi_bc_aa[,cn] = sapply(var_indices, function(v) { mean(psi_bc[v,paste0(cn,"_psi")], na.rm=T) })
    psi_bc_aa[,paste0("sd_",cn)] = sapply(var_indices, function(v) { sd(psi_bc[v,paste0(cn,"_psi")], na.rm=T) })
}

# order variants
psi_bc_aa$resi = sapply(subst_list, function(v){ if (length(v)>0) as.numeric(substr(v[1],2,nchar(v[1])-1)) else 0 })
psi_bc_aa$mut = sapply(subst_list, function(v){ if (length(v)>0) substr(v[1],nchar(v[1]),nchar(v[1])) else "" })
psi_bc_aa = psi_bc_aa[order(psi_bc_aa$n_subst, psi_bc_aa$resi, psi_bc_aa$mut, decreasing=F),]

# Calculate a mean of replica
cn = unique( sprintf("run%s_facs%d",dates,facs_rep) )
psi_bc_aa$mean_psi = apply(psi_bc_aa[,cn], MARGIN=1, mean, na.rm=T)      # hvorfor er F209L=NA med?
psi_bc_aa$sd_psi = apply(psi_bc_aa[,cn], MARGIN=1, sd, na.rm=T)

# Min-max normalization
psi_nonsense_bcaa = median(psi_bc_aa[which(psi_bc_aa$mut == "*"),"mean_psi"], na.rm=T)
scale_bcaa = 1/(psi_bc_aa[1,"mean_psi"] - psi_nonsense_bcaa)
psi_bc_aa$mean_psi_norm = (psi_bc_aa$mean_psi - psi_nonsense_bcaa) * scale_bcaa
psi_bc_aa$sd_psi_norm = psi_bc_aa$sd_psi * scale_bcaa
print(sprintf("Normalized PSI to WT = %.3f and nonsense = %.3f",psi_bc_aa[1,"mean_psi"],psi_nonsense_bcaa))


# ===========================================================================================================================================
#
# A curated and nice data frame 
#
# make a nice data frame of single mutants with good read counts
i = which(psi_aa$n_subst < 2 & ! is.na(psi_aa$mean_psi))
vamp = data.frame(var  = psi_aa[i,"var"])
nc_var = nchar(vamp$var)
vamp$wt    = substr(vamp$var, 1, 1)
vamp$resi  = as.numeric(substr(vamp$var, 2, nc_var-1))
vamp$mut   = substr(vamp$var, nc_var, nc_var)
vamp$barcodes   = psi_aa[i,"n_bc"]
vamp$reads = psi_aa[i,"all_reads"]
vamp$score = psi_aa[i,"mean_psi_norm"]
vamp$sd    = psi_aa[i,"sd_psi_norm"]
vamp[1,"resi"] = 0
vamp[1,"var"] = "WT"

# coverage, which single substitutions are present
i1 = which(psi_aa$n_subst==1)   # & psi_aa$mut != "*")
agg = aggregate(psi_aa[i1,"mut"], list(psi_aa[i1,"resi"]), paste0, collapse="")

residue = data.frame(resi=seq(nchar(wt[["aa"]])), wt=strsplit(wt[["aa"]],"")[[1]])
residue$mut = ""
residue[agg[,1],"mut"] = agg[,2]
residue$coverage = sapply(residue$mut, nchar)
no_nonsense = psi_aa$mut != "*"
i_no_nons = which(psi_aa$mut != "*")
residue$median_psi  = sapply(residue$resi, function(resi) { i = i1[which(psi_aa[i1,"resi"]==resi)];
                                                            median(psi_aa[intersect(i,i_no_nons),"mean_psi"], na.rm=T) })
residue$median_vamp = sapply(residue$resi, function(resi) { i = i1[which(psi_aa[i1,"resi"]==resi)]
                                                            median(psi_aa[intersect(i,i_no_nons),"mean_psi_norm"], na.rm=T) })
residue$vamp_sd_prop = sapply(residue$resi, function(resi) { i = i1[which(psi_aa[i1,"resi"]==resi)];
                                                             sds = psi_aa[intersect(i,i_no_nons),"sd_psi_norm"];
							     sqrt(sum(sds^2)) })


# ===========================================================================================================================================
#
# Synonymous WT
#
# aggragate DNA variants of WT
wt_i = which(psi_bc$var_aa=="")
agg = aggregate(wt_i, list(psi_bc[wt_i,"var_dna"]), paste, collapse=":")
var_indices = lapply( strsplit(agg[,2], ":"), as.numeric)
names(var_indices) = agg[,1]

# new data frame with WT DNA variants and PSI of all replicates
vamp_wt_syn = data.frame(var=names(var_indices), n_bc=sapply(var_indices,length))
subst_list = strsplit(vamp_wt_syn$var, ":")
vamp_wt_syn$n_subst = sapply(subst_list, length)
for (i in which(samples$gate=="A" & samples$tech_rep==1)) {
    cn = sprintf("run%s_facs%d",samples[i,"date"],samples[i,"facs_rep"])
    print(sprintf("Aggregate %s",cn))
    vamp_wt_syn[,cn] = sapply(var_indices, function(v) { mean(psi_bc[v,paste0(cn,"_psi")], na.rm=T) })
    vamp_wt_syn[,paste0("sd_",cn)] = sapply(var_indices, function(v) { sd(psi_bc[v,paste0(cn,"_psi")], na.rm=T) })
}

# Average over replicates
cn = unique( sprintf("run%s_facs%d",dates,facs_rep) )
vamp_wt_syn$mean_psi = apply(vamp_wt_syn[,cn], MARGIN=1, mean, na.rm=T)
vamp_wt_syn$sd_psi = apply(vamp_wt_syn[,cn], MARGIN=1, sd, na.rm=T)

# Normalize using same values as for full vamp
vamp_wt_syn$mean_psi_norm = (vamp_wt_syn$mean_psi - psi_nonsense) * scale
vamp_wt_syn$sd_psi_norm = vamp_wt_syn$sd_psi * scale

vamp_wt_syn[1,1] = "WT"

# plot
quartz(width=8, height=6)
breaks = seq(-0.3, 1.6, 0.02)
h_wt = hist(vamp_wt_syn[,"mean_psi_norm"], breaks=breaks, plot=F)
h_nons = hist(psi_aa[which(psi_aa$mut=="*"),"mean_psi_norm"], breaks=breaks, plot=F)
h_lib = hist(psi_aa[which(psi_aa$mut!="*" & psi_aa$n_subst==1),"mean_psi_norm"], breaks=breaks, plot=F)
plot(0,0,col=0, xlim=c(-.1,1.3), ylim=c(0,20.0), xlab="VAMP score", ylab="Density")
lines(h_lib$mids, h_lib$density, col=1, lwd=2)
lines(h_nons$mids, h_nons$density, col=2, lwd=2)
lines(h_wt$mids, h_wt$density, col=3, lwd=2)
legend("top", c("WT syn.","Nonsense","Library"), lty=1, lwd=2, col=c(3,2,1))
quartz.save("vamp_distributions.png", type="png")


# ===========================================================================================================================================
#
# Dump 
#
# save everything in R data file
save(vamp, vamp_wt_syn, psi_aa, psi_bc, raw_bc, samples, residue, wt, psi_nonsense, scale, file="aspa_vamp.rda")

# dump vamp scores as csv
write.table(vamp[,c("var","barcodes","reads","score","sd")], col.names=c("variant","barcodes","reads","vamp_score","vamp_std"),
            file="aspa_abundance.csv", sep=";", row.names=F, quote=F)

# dump vamp residue frame
write.table(residue, file="aspa_vamp_residues.csv", sep=";", row.names=F, quote=F)

# Dump synonymous WT
write.table(data.frame(var=vamp_wt_syn$var, barcodes=vamp_wt_syn$n_bc, vamp_score=vamp_wt_syn$mean_psi_norm, vamp_std=vamp_wt_syn$sd_psi_norm),
            file="aspa_vamp_wt_syn.csv", sep=";", row.names=F, quote=F)

# dump prism file
f = file("prism_mave_159_ASPA_vamp.txt", "wt")
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
write("#     technology: VAMPseq", f)
write("#     doi: unpublished", f)
write("#     year: 2022", f)
write("# variants:", f)
write(sprintf("#     number: %d",nrow(vamp)), f)
write(sprintf("#     coverage: %.2f",mean(residue$coverage > 0)), f)
write(sprintf("#     depth: %.2f",mean(nchar(gsub("*", "", residue[which(residue$coverage > 0),"mut"], fixed=T)))), f)
write("#     width: single mutants", f)
write("# columns:", f)
write("#     barcodes: Number of different barcodes observed for variant", f)
write("#     reads: Number of matched illumina reads observed across all replica and barcodes", f)
write("#     vamp_score: VAMPseq score is protein stability index (PSI) normalized to wild-type and nonsense variants", f)
write("#     vamp_std: VAMPseq score standard deviation between 11 replicates (4 biological)", f)
write("# --------------------", f)
write("#", f)
write("# Unpublished version of Sep 2022 - kristoffer.johansson@bio.ku.dk", f)
write("# ", f)
write(sprintf("%5s  %8s  %8s  %10s  %8s", "var", "barcodes", "reads", "vamp_score", "vamp_std"), f)
write(sprintf("%5s  %8d  %8d  %10.4f  %8.4f",vamp$var, vamp$barcodes, vamp$reads, vamp$score, vamp$sd), f)
close(f)


# ===========================================================================================================================================
#
# Analysis
#
date_rp = c()
date_mae = c()
for (date in unique(dates)) {
    cn = sprintf("run%s_facs%d_psi",date,unique(samples[which(samples$date==date),"facs_rep"]))
    
    rp = cor(psi_aa[,cn[1]], psi_aa[,cn[2]], method="pearson", use="complete.obs")
    mae = mean( abs(psi_aa[,cn[1]]-psi_aa[,cn[2]]), na.rm=T)
    if (length(cn) > 2) {
        rp = c(rp, cor(psi_aa[,cn[1]], psi_aa[,cn[3]], method="pearson", use="complete.obs"))
        rp = c(rp, cor(psi_aa[,cn[2]], psi_aa[,cn[3]], method="pearson", use="complete.obs"))
        mae = c(mae, mean( abs(psi_aa[,cn[1]]-psi_aa[,cn[3]]), na.rm=T))
        mae = c(mae, mean( abs(psi_aa[,cn[2]]-psi_aa[,cn[3]]), na.rm=T))
    }
    print(sprintf("Pearson correlations of variant PSI from date %d: %s",date,paste(rp,collapse=", ")))
    date_rp = c(date_rp, rp)
    
    # print(sprintf("Mean absolute error of variant PSI from %s: %s",samples[i,"date"],paste(mae,collapse=", "))) # what is i?
    print(sprintf("Mean absolute error of variant PSI from date %d: %s",date,paste(mae,collapse=", ")))
    date_mae = c(date_mae, mae)
}
print(sprintf("Replica correlations %.2f sd %.2f min %.2f max %.2f",mean(date_rp),sd(date_rp),min(date_rp),max(date_rp)))
print(sprintf("Replica MAE %.2f sd %.2f min %.2f max %.2f",mean(date_mae),sd(date_mae),min(date_mae),max(date_mae)))

facs_rp = c()
facs_mae = c()
for (facs in unique(facs_rep)) {
    cn = sprintf("run%s_facs%d_psi",unique(samples[which(samples$facs_rep==facs),"date"]),facs)
    
    rp = cor(psi_aa[,cn[1]], psi_aa[,cn[2]], method="pearson", use="complete.obs")
    mae = mean( abs(psi_aa[,cn[1]]-psi_aa[,cn[2]]), na.rm=T)    
    if (length(cn) > 2) {
        rp = c(rp, cor(psi_aa[,cn[1]], psi_aa[,cn[3]], method="pearson", use="complete.obs"))
        rp = c(rp, cor(psi_aa[,cn[2]], psi_aa[,cn[3]], method="pearson", use="complete.obs"))
        mae = c(mae, mean( abs(psi_aa[,cn[1]]-psi_aa[,cn[3]]), na.rm=T))
        mae = c(mae, mean( abs(psi_aa[,cn[2]]-psi_aa[,cn[3]]), na.rm=T))
    }
    # print(sprintf("Pearson correlations of variant PSI for FACS rep %d: %s",samples[i,"facs_rep"],paste(rp,collapse=", ")))
    print(sprintf("Pearson correlations of variant PSI for FACS rep %d: %s",facs,paste(rp,collapse=", ")))
    facs_rp = c(facs_rp, rp)

    # print(sprintf("Mean absolute error of variant PSI from %s: %s",samples[i,"date"],paste(mae,collapse=", ")))
    print(sprintf("Mean absolute error of variant PSI for FACS rep %d: %s",facs,paste(mae,collapse=", ")))
    facs_mae = c(date_mae, mae)
}
print(sprintf("Replica correlations %.2f sd %.2f min %.2f max %.2f",mean(facs_rp),sd(facs_rp),min(facs_rp),max(facs_rp)))
print(sprintf("Replica MAE %.2f sd %.2f min %.2f max %.2f",mean(facs_mae),sd(facs_mae),min(facs_mae),max(facs_mae)))


psi_bc_aa$psi_aa = psi_aa[match(psi_bc_aa$var,psi_aa$var),"mean_psi_norm"]
psi_bc_aa$psi_aa_sd = psi_aa[match(psi_bc_aa$var,psi_aa$var),"sd_psi_norm"]

# Compare AA-aggregation at raw reads level or PSI level
quartz(width=6, height=6)
plot(psi_bc_aa$psi_aa, psi_bc_aa$mean_psi_norm, xlab="New Vamp", ylab="Old Vamp", xlim=c(-0.05,0.2), ylim=c(-0.05,0.2))
abline(0,1)
i = which(psi_bc_aa$n_subst==1)
points(psi_bc_aa[i,"psi_aa"], psi_bc_aa[i,"mean_psi_norm"], col=2, pch=20)
inons = intersect(i, which(grepl('*', psi_bc_aa$var, fixed=T)))
points(psi_bc_aa[inons,"psi_aa"], psi_bc_aa[inons,"mean_psi_norm"], col=3, pch=20)
legend("topleft", c("All AA variants","Single AA var","Single nonsense"), pch=c(1,20,20), col=c(1,2,3))

breaks = seq(-0.1, 1.3, 0.05)
hnew = hist(psi_bc_aa[i,"psi_aa"], breaks=breaks, plot=F)
hold = hist(psi_bc_aa[i,"mean_psi_norm"], breaks=breaks, plot=F)
quartz(width=8, height=6)
plot(0, 0, col=0, xlim=c(-.1,1.3), ylim=c(0,5), xlab="VAMP score", ylab="Density")
lines(hnew$mids, hnew$density, col=2)
lines(hold$mids, hold$density, col=3)
legend("topright", c("New Vamp","Old Vamp"), lty=1, col=c(2,3))