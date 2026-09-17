## 82_Figure3_稳健性.R —— Figure 3：阳性对照与窗口稳健性
## Panel A 逐窗口 OR；Panel B 阳性对照（年龄/ADL/BMI → 跌倒）；Panel C 关键敏感性
suppressMessages({ library(data.table); library(survey); library(ggplot2); library(ragg); library(gridExtra) })
options(survey.lonely.psu = "adjust")
source("/Users/mac/Documents/Codex/2026-08-31/xian/work/tools/fig_lang.R")   # 中英图文字开关
B <- "/Users/mac/Documents/Codex/2026-08-31/xian/work"
T1 <- file.path(B, "外部复现-ELSA与HRS/t1")
dir.create(file.path(B, "外部复现-ELSA与HRS/res"), showWarnings = FALSE, recursive = TRUE)

prep <- function(rds) {
  d <- as.data.table(readRDS(rds))
  if ("ID" %in% names(d) && !"id" %in% names(d)) setnames(d, "ID", "id")
  d
}
coh <- list(
  CHARLS = list(d = prep(file.path(B, "新暴露-体力活动与握力轨迹/t1/样本_CHARLS.rds")),
                ids = "communityID", strata = "stratum", wt = "wt3", lab = "CHARLS 中国", fv = "fall_1518"),
  ELSA   = list(d = prep(file.path(T1, "样本_ELSA.rds")), ids = "id", strata = NULL, wt = "wt", lab = "ELSA 英国", fv = "fall"),
  HRS    = list(d = prep(file.path(T1, "样本_HRS.rds")), ids = "psu", strata = "strat", wt = "wt", lab = "HRS 美国", fv = "fall"),
  SHARE  = list(d = prep(file.path(T1, "样本_SHARE.rds")), ids = "psu", strata = "strat", wt = "wt", lab = "SHARE 欧洲", fv = "fall")
)
des_of <- function(x, ids, strata, wt) {
  x <- as.data.frame(x)
  if (is.null(strata)) svydesign(ids = as.formula(paste0("~", ids)), weights = as.formula(paste0("~", wt)), data = x)
  else svydesign(ids = as.formula(paste0("~", ids)), strata = as.formula(paste0("~", strata)),
                 weights = as.formula(paste0("~", wt)), data = x, nest = TRUE)
}

## ================= Panel B：阳性对照（年龄 / ADL / BMI → 跌倒）=================
cat("\n########## Panel B 阳性对照 ##########\n")
pc_rows <- list()
for (cn in names(coh)) {
  z <- coh[[cn]]; d <- copy(z$d)
  setnames(d, z$fv, "fall")
  d[, age10 := age/10][, bmi5 := bmi/5]
  vv <- c("fall","age10","bmi5","adl","edu")
  if (cn == "SHARE") vv <- c(vv, "eurod")
  dd <- d[complete.cases(d[, ..vv])]
  des <- des_of(dd, z$ids, z$strata, z$wt)
  specs <- list(list("年龄（每 +10 岁）", "age10", c("bmi5","adl","edu")),
                list("ADL 受限（每 +1 项）", "adl", c("age10","bmi5","edu")),
                list("BMI（每 +5 kg/m²）", "bmi5", c("age10","adl","edu")))
  for (s in specs) {
    f <- as.formula(paste0("fall ~ ", s[[2]], " + ", paste(s[[3]], collapse = " + ")))
    m <- svyglm(f, design = des, family = quasibinomial())
    co <- summary(m)$coefficients
    b <- co[s[[2]], 1]; se <- co[s[[2]], 2]; p <- co[s[[2]], 4]
    cat(sprintf("%-14s %-22s OR=%.3f (%.3f-%.3f) p=%.3g  n=%d\n",
                z$lab, s[[1]], exp(b), exp(b-1.96*se), exp(b+1.96*se), p, length(m$y)))
    pc_rows[[paste(cn, s[[1]])]] <- data.table(队列 = z$lab, 对照 = s[[1]], n = length(m$y),
                                               OR = exp(b), 下限 = exp(b-1.96*se), 上限 = exp(b+1.96*se), p = p)
  }
}
PC <- rbindlist(pc_rows)
fwrite(PC, "res/50_阳性对照.csv")

## ================= Panel A：逐窗口 =================
cat("\n########## Panel A 逐窗口 ##########\n")
w1 <- fread("res/02_ELSA_分窗口.csv")[, .(队列="ELSA 英国", 窗口)]
w2 <- fread("res/04_HRS_分窗口.csv")[, .(队列="HRS 美国", 窗口)]
w3 <- fread("res/12_SHARE_分窗口.csv")[, .(队列="SHARE 欧洲", 窗口)]
ch <- fread(file.path(B, "新暴露-体力活动与握力轨迹/res/08_2011到2015窗口.csv"))
ch <- ch[grepl("握力水平", 项) & grepl("跌倒", 项)]
WIN <- rbind(
  cbind(队列="ELSA 英国", 窗口=w1$窗口, fread("res/02_ELSA_分窗口.csv")[, .(n, 事件, OR, 下限, 上限, p)]),
  cbind(队列="HRS 美国",  窗口=w2$窗口, fread("res/04_HRS_分窗口.csv")[, .(n, 事件, OR, 下限, 上限, p)]),
  cbind(队列="SHARE 欧洲",窗口=w3$窗口, fread("res/12_SHARE_分窗口.csv")[, .(n, 事件, OR, 下限, 上限, p)]),
  data.table(队列="CHARLS 中国", 窗口="2011→2015（替代窗口）", n = NA_integer_, 事件 = NA_integer_,
             OR = ch$OR, 下限 = ch$下限, 上限 = ch$上限, p = ch$p))
WIN[, 标签 := paste0(队列, "  ", 窗口)]
for (cn in unique(WIN$队列)) {
  s <- WIN[队列 == cn]
  if (nrow(s) >= 3) {
    y <- log(s$OR); se <- (log(s$上限) - log(s$下限))/(2*1.96)
    w <- 1/se^2; Q <- sum(w*(y - weighted.mean(y, w))^2)
    cat(sprintf("%-12s %d 个窗口；OR 范围 %.3f–%.3f；同一队列内 Q=%.2f, df=%d, P=%.3f\n",
                cn, nrow(s), min(s$OR), max(s$OR), Q, nrow(s)-1, pchisq(Q, nrow(s)-1, lower.tail=FALSE)))
  } else cat(sprintf("%-12s 仅 1 个可用窗口（主分析用 2015+2018 合并）\n", cn))
}
fwrite(WIN, "res/51_逐窗口.csv")

## ================= Panel C：关键敏感性 =================
cat("\n########## Panel C 关键敏感性 ##########\n")
## 【2026-09-17 v2】原来这张表的 OR/CI/p 全部硬编码在脚本里，ELSA 设计修正后不会传导，
##   违反"数字必须由脚本现算"（台账 C104）。现全部改为从各队列结果文件读取。
SENS1 <- function(p, pat, orc = "OR", pc = "p", nf = "n") {
  d <- as.data.table(fread(p, encoding = "UTF-8"))
  key <- if ("项" %in% names(d)) d$项 else d$分析      # 各结果文件的标识列名不统一
  i <- grep(pat, key)
  if (length(i) == 0) stop("未匹配到: ", pat, " in ", p)
  r <- d[i[1]]
  list(OR = as.numeric(r[[orc]]), L = as.numeric(r[["下限"]]), U = as.numeric(r[["上限"]]),
       p = as.numeric(r[[pc]]), n = if (nf %in% names(r)) as.numeric(r[[nf]]) else NA_real_)
}
CH   <- file.path(B, "新暴露-体力活动与握力轨迹/res")
c_main <- SENS1(file.path(CH, "12_CHARLS_设计效应.csv"), "^CHARLS 跌倒（主模型）", orc = "OR_完整设计", pc = "p_完整设计")
c_alt  <- list(OR = ch$OR, L = ch$下限, U = ch$上限, p = ch$p, n = NA_real_)
c_sl   <- SENS1(file.path(CH, "09_三波斜率稳健性.csv"), "^三波斜率")
e_main <- SENS1("res/01_ELSA_主结果.csv", "^ELSA：.*跌倒（任意）")
e_exb  <- SENS1("res/06_ELSA_敏感性.csv", "仅基线未跌倒者")
e_wt   <- SENS1("res/06_ELSA_敏感性.csv", "结局波权重")
h_main <- SENS1("res/03_HRS_主结果.csv", "^HRS：.*跌倒（任意）")
h_exb  <- SENS1("res/05_HRS_敏感性.csv", "仅基线未跌倒者")
s_main <- SENS1("res/11_SHARE_主结果.csv", "^SHARE：.*跌倒")
s_exb  <- SENS1("res/13_SHARE_敏感性.csv", "仅基线未跌倒者")
s_wt   <- SENS1("res/33_SHARE_权重标准化敏感性.csv", "权重国内标准化")
sens <- rbind(
  data.table(队列="CHARLS 中国", 检验="主分析（2011→2015+2018）", n=c_main$n, OR=c_main$OR, 下限=c_main$L, 上限=c_main$U, p=c_main$p),
  data.table(队列="CHARLS 中国", 检验="替代窗口（2011→2015）", n=c_alt$n, OR=c_alt$OR, 下限=c_alt$L, 上限=c_alt$U, p=c_alt$p),
  data.table(队列="CHARLS 中国", 检验="握力三波斜率（每 SD）", n=c_sl$n, OR=c_sl$OR, 下限=c_sl$L, 上限=c_sl$U, p=c_sl$p),
  data.table(队列="ELSA 英国", 检验="主分析", n=e_main$n, OR=e_main$OR, 下限=e_main$L, 上限=e_main$U, p=e_main$p),
  data.table(队列="ELSA 英国", 检验="排除基线跌倒者", n=e_exb$n, OR=e_exb$OR, 下限=e_exb$L, 上限=e_exb$U, p=e_exb$p),
  data.table(队列="ELSA 英国", 检验="改用结局波权重", n=e_wt$n, OR=e_wt$OR, 下限=e_wt$L, 上限=e_wt$U, p=e_wt$p),
  data.table(队列="HRS 美国", 检验="主分析", n=h_main$n, OR=h_main$OR, 下限=h_main$L, 上限=h_main$U, p=h_main$p),
  data.table(队列="HRS 美国", 检验="排除基线跌倒者", n=h_exb$n, OR=h_exb$OR, 下限=h_exb$L, 上限=h_exb$U, p=h_exb$p),
  data.table(队列="SHARE 欧洲", 检验="主分析", n=s_main$n, OR=s_main$OR, 下限=s_main$L, 上限=s_main$U, p=s_main$p),
  data.table(队列="SHARE 欧洲", 检验="排除基线跌倒者", n=s_exb$n, OR=s_exb$OR, 下限=s_exb$L, 上限=s_exb$U, p=s_exb$p),
  data.table(队列="SHARE 欧洲", 检验="权重国内标准化", n=s_wt$n, OR=s_wt$OR, 下限=s_wt$L, 上限=s_wt$U, p=s_wt$p)
)
fwrite(sens, "res/52_关键敏感性.csv")
cat(sprintf("敏感性检验共 %d 条，全部 OR<1 的比例：%d/%d\n", nrow(sens), sum(sens$OR<1), nrow(sens)))

## ================= 作图 =================
## ---- 图内文字的英文映射（中文版文章出中文图，英文版文章出英文图）----
TR <- c(
  "年龄（每 +10 岁）"      = "Age (per +10 years)",
  "ADL 受限（每 +1 项）"   = "ADL limitation (per +1 item)",
  "BMI（每 +5 kg/m²）"     = "BMI (per +5 kg/m²)",
  "主分析（2011→2015+2018）" = "Primary (2011→2015+2018)",
  "替代窗口（2011→2015）"  = "Alternative window (2011→2015)",
  "握力三波斜率（每 SD）"  = "Grip slope across waves (per SD)",
  "主分析"                 = "Primary analysis",
  "排除基线跌倒者"         = "Excluding baseline fallers",
  "改用结局波权重"         = "Outcome-wave weights",
  "权重国内标准化"         = "Weights standardised within country",
  "2011→2015（替代窗口）"  = "2011→2015 (alternative)"
)
tr <- function(x) if (EN) { y <- TR[as.character(x)]; ifelse(is.na(y), as.character(x), unname(y)) } else as.character(x)
## 合并估计由结果文件现算（台账 C104：图内数字不得手写）
POOL_FALL <- tryCatch({
  m <- fread(file.path(B, "外部复现-ELSA与HRS/res/08_多队列合并.csv"))
  as.numeric(m[结局 == "跌倒（任意）", 合并OR][1])
}, error = function(e) NA_real_)
if (is.na(POOL_FALL)) stop("缺少 res/08_多队列合并.csv，请先运行 30_多队列对照与合并.R")

WIN[, 签名 := ifelse(p<0.05, "P<0.05", "P≥0.05")]
WIN[, 标签 := paste0(CLAB(队列), "  ", tr(窗口))]
WIN[, 标签 := factor(标签, levels = rev(标签))]
A <- ggplot(WIN, aes(x = OR, y = 标签, colour = 签名)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey45") +
  geom_vline(xintercept = POOL_FALL, linetype = 3, colour = "#2E74B5") +
  geom_errorbarh(aes(xmin = 下限, xmax = 上限), height = 0.22, linewidth = 0.5) +
  geom_point(size = 2.1) +
  scale_colour_manual(values = c("P<0.05"="#C0504D","P≥0.05"="grey45")) +
  labs(title = TT(sprintf("A  逐随访窗口的握力—跌倒关联（%d 个窗口）", nrow(WIN)),
                  sprintf("A  Grip strength and falls in each follow-up window (%d windows)", nrow(WIN))),
       subtitle = TT(sprintf("蓝虚线＝四队列合并估计 %.3f。全部 %d 个窗口点估计均 <1，方向无一例外", POOL_FALL, nrow(WIN)),
                     sprintf("Blue dashed line = pooled estimate %.3f. All %d window-specific point estimates are below 1.",
                             POOL_FALL, nrow(WIN))),
       x = TT("OR（每 +1 SD 握力）", "OR per 1-SD higher grip strength"), y = NULL, colour = NULL) +
  theme_minimal() +
  theme(plot.title = element_text(size = 11.5, colour = "#1F3A5F"),
        plot.subtitle = element_text(size = 9, colour = "#5A5A5A"),
        axis.text.y = element_text(size = 7.4), axis.text.x = element_text(size = 9),
        legend.position = "bottom", legend.text = element_text(size = 8.2),
        legend.title = element_blank(), panel.grid.minor = element_blank())

PC[, 标签 := paste0(CLAB(队列), "  ", tr(对照))]
PC[, 标签 := factor(标签, levels = rev(unique(标签)))]
Bp <- ggplot(PC, aes(x = OR, y = 标签, colour = 对照)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey45") +
  geom_errorbarh(aes(xmin = 下限, xmax = 上限), height = 0.2, linewidth = 0.55) +
  geom_point(size = 2.4) +
  scale_colour_manual(values = c("年龄（每 +10 岁）"="#2E74B5","ADL 受限（每 +1 项）"="#C0504D","BMI（每 +5 kg/m²）"="#4E7B4E"),
                      labels = tr) +
  scale_x_log10(breaks = c(1, 1.2, 1.5, 2, 3, 5)) +
  labs(title = TT("B  阳性对照：已知危险因素 → 跌倒", "B  Positive controls: established risk factors for falls"),
       subtitle = TT("年龄与 ADL 受限在四队列全部显著升高风险；\nBMI 在三队列显著、CHARLS 为零效应（真实结果，未删）",
                     "Age and ADL limitation raised the odds of falls in all four cohorts;\nBMI was significant in three cohorts, a true null in CHARLS"),
       x = TT("OR（对数刻度）", "OR (log scale)"), y = NULL, colour = NULL) +
  theme_minimal() +
  theme(plot.title = element_text(size = 11.5, colour = "#1F3A5F"),
        plot.subtitle = element_text(size = 9, colour = "#5A5A5A"),
        axis.text.y = element_text(size = 8.2), axis.text.x = element_text(size = 9),
        legend.position = "bottom", legend.text = element_text(size = 8),
        legend.title = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = margin(t = 6, r = 30, b = 6, l = 6))

sens[, 签名 := ifelse(p<0.05, "P<0.05", "P≥0.05")]
sens[, 标签 := paste0(CLAB(队列), "  ", tr(检验))]
sens[, 标签 := factor(标签, levels = rev(标签))]
Cp <- ggplot(sens, aes(x = OR, y = 标签, colour = 签名)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey45") +
  geom_errorbarh(aes(xmin = 下限, xmax = 上限), height = 0.2, linewidth = 0.55) +
  geom_point(size = 2.2) +
  scale_colour_manual(values = c("P<0.05"="#C0504D","P≥0.05"="grey45")) +
  labs(title = TT(sprintf("C  关键敏感性检验（%d 条）", nrow(sens)),
                  sprintf("C  Key sensitivity analyses (%d)", nrow(sens))),
       subtitle = TT("替代窗口、排除基线跌倒者、换权重口径后方向全部不变",
                     "Direction unchanged across alternative windows,\nexclusion of baseline fallers and alternative weights"),
       x = TT("OR（每 +1 SD 握力）", "OR per 1-SD higher grip strength"), y = NULL, colour = NULL) +
  theme_minimal() +
  theme(plot.title = element_text(size = 11.5, colour = "#1F3A5F"),
        plot.subtitle = element_text(size = 9, colour = "#5A5A5A"),
        axis.text.y = element_text(size = 7.8), axis.text.x = element_text(size = 9),
        legend.position = "bottom", legend.text = element_text(size = 8.2),
        legend.title = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = margin(t = 6, r = 30, b = 6, l = 6))

agg_png(paste0("res/53_Figure3", FSFX, ".png"), width = 3000, height = 2500, res = 210, background = "white")
grid.arrange(A, Bp, Cp, layout_matrix = rbind(c(1,1), c(2,3)), heights = c(1.25, 1))
dev.off()
cat("\n图已生成: res/53_Figure3", FSFX, ".png\n", sep = "")
