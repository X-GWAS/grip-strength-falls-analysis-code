## 30_figure.R —— 反向 MR 图（ragg，普通字重，勿用粗体）
suppressMessages({ library(data.table); library(ggplot2); library(ragg); library(gridExtra) })
source("/Users/mac/Documents/Codex/2026-08-31/xian/work/tools/fig_lang.R")   # 中英图文字开关
setwd("/Users/mac/Documents/Codex/2026-08-31/xian/work/反向MR-跌倒到握力")
dir.create("res", showWarnings = FALSE)
## ---- 图内文字的英文映射（中文版文章出中文图，英文版文章出英文图）----
TRE <- c(
  "名义显著 (p<0.05)" = "Nominally significant (p<0.05)", "不显著" = "Not significant",
  "方向相反（跌倒多→握力高）" = "Opposite direction (more falls -> higher grip)",
  "方向一致（跌倒多→握力低）" = "Concordant direction (more falls -> lower grip)",
  "主分析" = "Primary", "工具变量拆分" = "Instrument subsets",
  "骨密度直接检验" = "Direct bone-density test", "重叠对照臂" = "Overlapping control arm",
  "阳性对照" = "Positive control", "弱阳性" = "Weakly positive", "阴性对照" = "Negative control"
)
trw <- function(x) if (EN) { y <- TRE[as.character(x)]; ifelse(is.na(y), as.character(x), unname(y)) } else as.character(x)
TPL <- function(x) {
  if (!EN) return(as.character(x))
  x <- as.character(x)
  x <- gsub("^跌倒 → ", "Falls -> ", x)
  x <- gsub("^骨密度 eBMD → ", "eBMD -> ", x)
  x <- gsub("^UKB 跌倒 → UKB 握力", "UKB falls -> UKB grip", x)
  x <- gsub("^身高 → 握力", "Height -> grip", x)
  x <- gsub("^四肢肌肉量 → 握力", "Appendicular lean mass -> grip", x)
  x <- gsub("^BMI → 握力", "BMI -> grip", x)
  x <- gsub("^LDL → 握力", "LDL -> grip", x)
  x <- gsub("握力（右手）", "grip (right hand)", x)
  x <- gsub("握力（左手）", "grip (left hand)", x)
  x <- gsub("低握力 EWGSOP", "low grip strength (EWGSOP)", x)
  x <- gsub("跌倒 → 握力", "falls -> grip", x)
  x <- gsub("排除骨密度位点", "excluding bone-density loci", x)
  x <- gsub("仅骨密度位点", "bone-density loci only", x)
  x <- gsub("加权中位数", "weighted median", x)
  x <- gsub("（阳性对照）", " (positive control)", x, fixed = TRUE)
  x <- gsub("（阴性对照）", " (negative control)", x, fixed = TRUE)
  x <- gsub("（弱阳性）", " (weakly positive)", x, fixed = TRUE)
  gsub("骨密度", "bone density", x)
}
main_dt <- as.data.table(fread("res/01_反向MR_主结果.csv"))
set4    <- as.data.table(fread("res/14_工具变量集拆分.csv"))
ebmd    <- as.data.table(fread("res/21_eBMD到握力.csv"))
ovl     <- as.data.table(fread("res/16_重叠对照臂.csv"))
ctl     <- as.data.table(fread("res/15b_对照_补充.csv"))
ivw <- function(x) x[方法 == "Inverse variance weighted"]

## ---------- A 逐 SNP ----------
d <- as.data.table(fread("res/20_逐SNP_Wald比.csv"))
g <- unique(as.data.table(fread("res/12_工具变量表型画像.csv"))[, .(SNP, nearest_genes)])
d <- merge(d, g, by = "SNP", all.x = TRUE)
d[, lab := sprintf("%s (%s)", SNP, ifelse(is.na(nearest_genes), "", nearest_genes))]
d[, sig := ifelse(p_wald < 0.05, "名义显著 (p<0.05)", "不显著")]
d[, dir := ifelse(wald > 0, "方向相反（跌倒多→握力高）", "方向一致（跌倒多→握力低）")]
d[, sig_d := trw(sig)]; d[, dir_d := trw(dir)]
d[, `:=`(lo = wald - 1.96*wald_se, hi = wald + 1.96*wald_se)]
setorder(d, wald); d[, lab := factor(lab, levels = lab)]
ivw_b <- ivw(main_dt[结局 == "握力（右手，UKB）"])$b

A <- ggplot(d, aes(x = wald, y = lab, colour = sig_d, shape = dir_d)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey45") +
  geom_vline(xintercept = ivw_b, linetype = 3, colour = "#2E74B5") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.22, linewidth = 0.55) +
  geom_point(size = 2.4) +
  scale_colour_manual(values = setNames(c("#C0504D", "grey45"), trw(c("名义显著 (p<0.05)", "不显著")))) +
  scale_shape_manual(values = setNames(c(16, 17), trw(c("方向一致（跌倒多→握力低）", "方向相反（跌倒多→握力高）")))) +
  labs(title = TT("A  逐个工具变量的 Wald 比：遗传预测跌倒 → 握力（右手，SD）",
                  "A  Wald ratio for each instrument: genetically predicted falls and grip strength (right hand, SD)"),
       subtitle = TT("蓝虚线＝合并 IVW。14 个位点中 6 个名义显著，但方向相反，合并后互相抵消",
                     "Blue dashed line = pooled inverse-variance weighted estimate. Six of 14 loci were nominally significant, but in opposing directions, so they cancel on pooling"),
       x = TT("握力变化（SD）／每 1 单位 log-OR 跌倒遗传易感性",
              "Change in grip strength (SD) per 1-unit higher log-odds genetic liability to falls"),
       y = NULL, colour = NULL, shape = NULL) +
  theme_minimal() +
  theme(plot.title = element_text(size = 12, colour = "#1F3A5F"),
        plot.subtitle = element_text(size = 9.2, colour = "#5A5A5A"),
        axis.text.y = element_text(size = 8.4), axis.text.x = element_text(size = 9),
        legend.position = "bottom", legend.text = element_text(size = 8.4),
        legend.title = element_blank(), legend.box = "vertical",
        panel.grid.minor = element_blank(), plot.margin = margin(6, 12, 4, 6))

## ---------- B 反向 MR 汇总（不含对照） ----------
B <- rbind(
  data.table(项 = "跌倒 → 握力（右手）", b = ivw(main_dt[结局=="握力（右手，UKB）"])$b,
             se = ivw(main_dt[结局=="握力（右手，UKB）"])$se, grp = "主分析"),
  data.table(项 = "跌倒 → 握力（右手）·加权中位数", b = main_dt[结局=="握力（右手，UKB）" & 方法=="Weighted median", b],
             se = main_dt[结局=="握力（右手，UKB）" & 方法=="Weighted median", se], grp = "主分析"),
  data.table(项 = "跌倒 → 握力（左手）", b = ivw(main_dt[结局=="握力（左手，UKB）"])$b,
             se = ivw(main_dt[结局=="握力（左手，UKB）"])$se, grp = "主分析"),
  data.table(项 = "跌倒 → 低握力 EWGSOP", b = ivw(main_dt[结局=="低握力 EWGSOP（UKB）"])$b,
             se = ivw(main_dt[结局=="低握力 EWGSOP（UKB）"])$se, grp = "主分析"),
  data.table(项 = sprintf("跌倒 → 握力：排除骨密度位点（n=%d）", ivw(set4[工具变量集=="排除骨密度位点"])$n_snp),
             b = ivw(set4[工具变量集=="排除骨密度位点"])$b, se = ivw(set4[工具变量集=="排除骨密度位点"])$se, grp = "工具变量拆分"),
  data.table(项 = sprintf("跌倒 → 握力：仅骨密度位点（n=%d）", ivw(set4[工具变量集=="仅骨密度位点"])$n_snp),
             b = ivw(set4[工具变量集=="仅骨密度位点"])$b, se = ivw(set4[工具变量集=="仅骨密度位点"])$se, grp = "工具变量拆分"),
  data.table(项 = sprintf("骨密度 eBMD → 握力（n=%d）", ebmd[method=="Inverse variance weighted", nsnp]),
             b = ebmd[method=="Inverse variance weighted", b], se = ebmd[method=="Inverse variance weighted", se],
             grp = "骨密度直接检验"),
  data.table(项 = sprintf("UKB 跌倒 → UKB 握力（n=%d）", ivw(ovl[暴露=="ukb-b-2535"])$n_snp),
             b = ivw(ovl[暴露=="ukb-b-2535"])$b, se = ivw(ovl[暴露=="ukb-b-2535"])$se, grp = "重叠对照臂")
)
B[, `:=`(lo = b - 1.96*se, hi = b + 1.96*se)]
B[, lab_d := TPL(项)]
B[, lab := factor(lab_d, levels = rev(lab_d))]
pal <- c("主分析"="#2E74B5","工具变量拆分"="#C0504D","骨密度直接检验"="#8C6D1F","重叠对照臂"="#7030A0")
names(pal) <- trw(names(pal))

pB <- ggplot(B, aes(x = b, y = lab, colour = grp)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey45") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.22, linewidth = 0.6) +
  geom_point(size = 2.6) +
  scale_colour_manual(values = pal) +
  labs(title = TT("B  反向 MR 汇总（不含对照）", "B  Reverse-direction MR, summary estimates (controls excluded)"),
       subtitle = TT("8 个估计全部跨越 0；拆分工具变量后点估计方向翻转",
                     "All eight estimates cross zero; the point estimate flips direction when the instrumental variables are split"),
       x = TT("握力变化（SD）／每 1 单位 log-OR 跌倒遗传易感性",
              "Change in grip strength (SD) per 1-unit higher log-odds genetic liability to falls"),
       y = NULL, colour = NULL) +
  theme_minimal() +
  theme(plot.title = element_text(size = 12, colour = "#1F3A5F"),
        plot.subtitle = element_text(size = 9.2, colour = "#5A5A5A"),
        axis.text.y = element_text(size = 8.6), axis.text.x = element_text(size = 9),
        legend.position = "bottom", legend.text = element_text(size = 8.4),
        legend.title = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = margin(6, 12, 4, 6))

## ---------- C 对照 ----------
C <- data.table(
  项 = c("身高 → 握力（阳性对照）", "四肢肌肉量 → 握力（阳性对照）",
         "BMI → 握力（弱阳性）", "LDL → 握力（阴性对照）"),
  b  = c(ctl[grepl("身高", 对照), b_IVW], ctl[grepl("四肢肌肉量", 对照), b_IVW],
         ctl[grepl("^BMI", 对照), b_IVW], ctl[grepl("LDL", 对照), b_IVW]),
  se = c(ctl[grepl("身高", 对照), se_IVW], ctl[grepl("四肢肌肉量", 对照), se_IVW],
         ctl[grepl("^BMI", 对照), se_IVW], ctl[grepl("LDL", 对照), se_IVW])
)
C[, `:=`(lo = b - 1.96*se, hi = b + 1.96*se)]
C[, lab_d := TPL(项)]
C[, lab := factor(lab_d, levels = rev(lab_d))]
C[, 类 := ifelse(grepl("阳性", 项), "阳性对照", ifelse(grepl("阴性", 项), "阴性对照", "弱阳性"))]
C[, 类_d := trw(类)]

pC <- ggplot(C, aes(x = b, y = lab, colour = 类_d)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey45") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2, linewidth = 0.6) +
  geom_point(size = 2.8) +
  scale_colour_manual(values = setNames(c("#4E7B4E","#B07A2A","grey45"), trw(c("阳性对照","弱阳性","阴性对照")))) +
  labs(title = TT("C  对照检验（结局同为 UKB 握力）", "C  Control analyses with UK Biobank grip strength as the outcome"),
       subtitle = TT("身高、四肢肌肉量极显著为正，LDL 精确为零：流程与结局数据可靠",
                     "Height and appendicular lean mass were strongly positive and LDL was precisely null: the pipeline and the outcome data behave as expected"),
       x = TT("握力变化（SD）／每 1 单位暴露", "Change in grip strength (SD) per 1-unit higher exposure"),
       y = NULL, colour = NULL) +
  theme_minimal() +
  theme(plot.title = element_text(size = 12, colour = "#1F3A5F"),
        plot.subtitle = element_text(size = 9.2, colour = "#5A5A5A"),
        axis.text.y = element_text(size = 8.6), axis.text.x = element_text(size = 9),
        legend.position = "bottom", legend.text = element_text(size = 8.4),
        legend.title = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = margin(6, 12, 4, 6))

agg_png(paste0("res/30_反向MR图", FSFX, ".png"), width = 3000, height = 2600, res = 210, background = "white")
grid.arrange(A, pB, pC, layout_matrix = rbind(c(1,1), c(2,3)), heights = c(1, 1))
dev.off()
cat("图已生成: res/30_反向MR图", FSFX, ".png\n", sep = "")
