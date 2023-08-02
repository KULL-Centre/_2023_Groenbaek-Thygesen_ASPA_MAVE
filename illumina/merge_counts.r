options(width=160, digits=4, stringsAsFactors=F)

bc_file = "barcode_map.rda"
load(bc_file)
print(sprintf("Loaded barcode map for %s from %s",wt[["name"]],bc_file))

raw = data.frame(barcode = bc_map$barcode, var_aa = bc_map$var_aa, var_dna = bc_map$var_dna)

for (file in files) {
    d = read.table(file)
    fl = strsplit(file, "/")[[1]]
    stopifnot(substr(fl[length(fl)], nchar(fl[length(fl)])-10, nchar(fl[length(fl)])) == "_counts.txt")
    file_id = substr(fl[length(fl)], 1, nchar(fl[length(fl)])-11)
    raw[,file_id] = d[match(raw$barcode,d$V1),"V2"]
    na_mask = is.na(raw[,file_id])
    raw[na_mask,file_id] = 0
    # Report
    mapped_bc = sum(! na_mask)
    mapped_reads = sum(raw[which(! na_mask),file_id])
    print(sprintf("Mapped %.2f%% of barcodes and %.2f%% of reads from %s (%d barcodes and %d reads)",
                  mapped_bc/nrow(d)*100, mapped_reads/sum(d$V2)*100, file_id, mapped_bc, mapped_reads))
    
}

save(raw, wt, file="raw.rda")
