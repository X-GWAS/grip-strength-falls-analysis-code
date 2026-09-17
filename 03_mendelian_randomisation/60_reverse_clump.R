suppressMessages({ library(data.table); library(ieugwasr) })
Sys.setenv(OPENGWAS_JWT = trimws(readLines("/tmp/opengwas_jwt", warn = FALSE)[1]))
setwd("/Users/mac/Documents/Codex/2026-08-31/xian/work/反向MR-跌倒到握力")

fin <- fread("finn_FALLS_gws.tsv", sep = "\t")
setnames(fin, gsub("^#", "", names(fin)))
cat("全基因组显著 SNP:", nrow(fin), "\n")
cat("p 值范围:", sprintf("%.2e ~ %.2e", min(fin$pval), max(fin$pval)), "\n")
fin[, rsid := sub(",.*$", "", rsids)]
cat("有 rsID 的:", sum(fin$rsid != "NA" & fin$rsid != ""), "\n")
fin <- fin[rsid != "NA" & rsid != ""]
cat("去重前 rsID 数:", nrow(fin), " 唯一 rsID:", uniqueN(fin$rsid), "\n")
fin <- fin[!duplicated(rsid)]
cat("\n最小的 20 个位点:\n")
print(head(fin[order(pval), .(rsid, nearest_genes, pval, beta, sebeta, af_alt)], 20), row.names = FALSE)

cat("\n===== LD 剪枝 (r2<0.001, 10Mb) =====\n")
cl <- ld_clump(dat = fin[, .(rsid, pval, id = "finn_R12_FALLS")],
               clump_kb = 10000, clump_r2 = 0.001,
               clump_p = 5e-8, pop = "EUR")
cat("剪枝后独立工具变量:", nrow(cl), "\n")
print(cl, row.names = FALSE)
fwrite(cl, "01_工具变量_剪枝后.csv")
