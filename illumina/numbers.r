
##
##   Replica correlations
##
stat_rep = function(df) {
    # cor_mat = cor(df, method="pearson", use="complete.obs")
    # rps = cor_mat[upper.tri(cor_mat)]
    # for (ic in seq(ncol(df))) {
    #     print(sprintf("Replica %2d summary:",ic))
    # 	print(summary(df[,ic]))
    # }
    rss = c()
    rps = c()
    maes = c()
    for (ic1 in seq(ncol(df)-1)) {
        for ( ic2 in seq(ic1+1,ncol(df)) ) {
            rs = cor(df[,ic1], df[,ic2], method="spearman", use="complete.obs")
	    rss = c(rss, rs)
            rp = cor(df[,ic1], df[,ic2], method="pearson", use="complete.obs")
	    rps = c(rps, rp)
	    mae = mean( abs(df[,ic2] - df[,ic1]), na.rm=T)
	    maes = c(maes, mae)
	}
    }
    print(sprintf("Average Pearson %.4f, range %.4f - %.4f", mean(rps), min(rps), max(rps)))
    print(sprintf("Average Spearman %.4f, range %.4f - %.4f", mean(rss), min(rss), max(rss)))
    print(sprintf("Average MAE %.4f, range %.4f - %.4f", mean(maes), min(maes), max(maes)))
    return( c(mean(rps), mean(maes)) )
}


# abundance data
load("aspa_vamp.rda")
samples_vamp = samples
scale_vamp = scale
res_vamp = residue

# nonsense mutations not included in residue data frame
psi_aa$median = residue[match(psi_aa$resi,residue$resi), "median_vamp"]
i = which(psi_aa$n_subst==1 & psi_aa$mut != '*' & (! is.na(psi_aa$mean_psi_norm) ) & (! psi_aa$resi %in% c(1,188,189,242,243)) )

var_data = var(psi_aa[i,"mean_psi_norm"])
var_model = var(psi_aa[i,"mean_psi_norm"] - psi_aa[i,"median"])
ftest = (1 - var_model / var_data )
print(sprintf("Median position explains %.2f %% of data variance (%.2f -> %.2f) in abundance score units (data sd = %.2f)",
              ftest*100, var_data, var_model, sqrt(var_data)))

# crude histogram of variance
breaks = seq(-2,2,0.1)
hist(psi_aa[i,"mean_psi_norm"], breaks=breaks)
h_model = hist(psi_aa[i,"mean_psi_norm"]-psi_aa[i,"median"], breaks=breaks, plot=F)
lines(h_model$mids, h_model$counts, col=2)


n_nons_tot = nrow(res_vamp)-1  # stop at first res dosn't measure anything
n_nons_meas = sum(grepl('*', res_vamp$mut, fixed=T))
print(sprintf("VAMP nonsense scores determined for %d of %d positions (%.2f%%)", n_nons_meas, n_nons_tot, n_nons_meas/n_nons_tot*100))

aa_one = strsplit("ACDEFGHIKLMNPQRSTVWY","")[[1]]
n_mis_tot = nrow(res_vamp)*19
n_mis_meas = sum(vamp$mut %in% aa_one)
print(sprintf("VAMP missense scores determined for %d of %d positions (%.2f%%)", n_mis_meas, n_mis_tot, n_mis_meas/n_mis_tot*100))

n_tot = n_mis_tot + n_nons_tot
n_meas = n_mis_meas + n_nons_meas
print(sprintf("Total VAMP scores for %d out of %d (%.2f%%) variants measured (excl WT)", n_meas, n_tot, n_meas/n_tot*100))

print("================= Stat abundance scores =================")
n_nons_tot = nrow(res_vamp)-1  # stop at first res dosn't measure anything
n_nons_meas = sum(grepl('*', res_vamp$mut, fixed=T))
print(sprintf("VAMP nonsense scores determined for %d of %d positions (%.2f%%)", n_nons_meas, n_nons_tot, n_nons_meas/n_nons_tot*100))

n_mis_tot = nrow(res_vamp)*19
n_mis_meas = sum(vamp$mut %in% aa_one)
print(sprintf("VAMP missense scores determined for %d of %d positions (%.2f%%)", n_mis_meas, n_mis_tot, n_mis_meas/n_mis_tot*100))

n_tot = n_mis_tot + n_nons_tot
n_meas = n_mis_meas + n_nons_meas
print(sprintf("Total VAMP scores for %d out of %d (%.2f%%) variants measured (excl WT)", n_meas, n_tot, n_meas/n_tot*100))

# i1_psi = which(psi_aa$n_subst == 1)
is = which(samples_vamp$gate=="A" & samples_vamp$tech_rep==1)
vamp_replica_cn = sprintf("run%s_facs%d_psi",samples[is,"date"],samples[is,"facs_rep"])
# stat_rep(psi_aa[i1_psi,vamp_replica_cn])
dfv = (psi_aa[match(vamp$var,psi_aa$var),vamp_replica_cn] - psi_nonsense) * scale_vamp   # WT not matched
rownames(dfv) = vamp$var
stat_rep(dfv)
print("=========================================================")
print("")

# quartz(height=8, width=8)
jpeg("vamp_replica.jpg", width=19, heigh=19, units="cm", res=300, pointsize=9, quality=90)
plot(dfv, pch=".", labels=sprintf("Run %d\nFACS %d",samples[is,"date"],samples[is,"facs_rep"]), upper.panel=NULL)
title(xlab="Abundance score", line=4)
title(ylab="Abundance score", line=3)
# quartz.save("vamp_replica.png", type="png")
dev.off()


load("aspa_tox.rda")
samples_tox = samples
scale_tox = scale
res_tox = residue

n_nons_tot = nrow(res_tox)-1  # stop at first res dosn't measure anything
n_nons_meas = sum(grepl('*', res_tox$mut, fixed=T))
print(sprintf("Tox nonsense scores determined for %d of %d positions (%.2f%%)", n_nons_meas, n_nons_tot, n_nons_meas/n_nons_tot*100))

aa_one = strsplit("ACDEFGHIKLMNPQRSTVWY","")[[1]]
n_mis_tot = nrow(res_tox)*19
n_mis_meas = sum(tox$mut %in% aa_one)
print(sprintf("Tox missense scores determined for %d of %d positions (%.2f%%)", n_mis_meas, n_mis_tot, n_mis_meas/n_mis_tot*100))

n_tot = n_mis_tot + n_nons_tot
n_meas = n_mis_meas + n_nons_meas
print(sprintf("Total tox scores for %d out of %d (%.2f%%) variants measured (excl WT)", n_meas, n_tot, n_meas/n_tot*100))


print("================= Stat toxicity scores  =================")

n_nons_tot = nrow(res_tox)-1  # stop at first res dosn't measure anything
n_nons_meas = sum(grepl('*', res_tox$mut, fixed=T))
print(sprintf("Tox nonsense scores determined for %d of %d positions (%.2f%%)", n_nons_meas, n_nons_tot, n_nons_meas/n_nons_tot*100))

aa_one = strsplit("ACDEFGHIKLMNPQRSTVWY","")[[1]]
n_mis_tot = nrow(res_tox)*19
n_mis_meas = sum(tox$mut %in% aa_one)
print(sprintf("Tox missense scores determined for %d of %d positions (%.2f%%)", n_mis_meas, n_mis_tot, n_mis_meas/n_mis_tot*100))

n_tot = n_mis_tot + n_nons_tot
n_meas = n_mis_meas + n_nons_meas
print(sprintf("Total tox scores for %d out of %d (%.2f%%) variants measured (excl WT)", n_meas, n_tot, n_meas/n_tot*100))

# i1_tox = which(growth_aa$n_aa == 1)
tox_replica_cn = sprintf("run%d_slope",seq(4))
# stat_rep(growth_aa[i1_tox,tox_replica_cn])
dft = (growth_aa[match(tox$var,growth_aa$var),tox_replica_cn] - slope_wt) * scale_tox
rownames(dft) = tox$var
stat_rep(dft)
print("=========================================================")
print("")

# quartz(height=6, width=6)
jpeg("tox_replica.jpg", width=8, heigh=8, units="cm", res=300, pointsize=9, quality=90)
plot(dft, pch=".", labels=sprintf("Replica %d",seq(4)), upper.panel=NULL)
title(xlab="Toxicity score", line=4)
title(ylab="Toxicity score", line=3)
# quartz.save("tox_replica.png", type="png")
dev.off()




