## 85_Table2_森林图.R —— Table 2 亚组森林图（四队列逐层 + 合并）
suppressMessages({ library(data.table); library(ggplot2); library(ragg) })
source("/Users/mac/Documents/Codex/2026-08-31/xian/work/tools/fig_lang.R")   # 中英图文字开关
RES <- "/Users/mac/Documents/Codex/2026-08-31/xian/work/外部复现-ELSA与HRS/res"
R    <- fread(file.path(RES, "60_Table2_分层估计.csv"))
P    <- fread(file.path(RES, "62_Table2_层合并.csv"))
## 【2026-09-17 修订】层间交互 P 的范围原来硬编码在副标题里（写死 0.39–0.78），
##   主模型改动后会悄悄失配；改为每次从 63 号文件实时计算。
IP   <- fread(file.path(RES, "63_Table2_层间交互P.csv"))
ipc  <- grep("P$|P$|P$", names(IP), value = TRUE)
ipc  <- names(IP)[grepl("P", names(IP))][1]
ipr  <- range(as.numeric(IP[[ipc]]), na.rm = TRUE)
IP_TXT_ZH <- sprintf("%.2f–%.2f", ipr[1], ipr[2])
IP_TXT_EN <- sprintf("%.2f–%.2f", ipr[1], ipr[2])
cat("层间交互 P 范围：", IP_TXT_EN, "\n")

LV <- c("<65","65–74","≥75","<25","25–29.9","≥30","0","1","≥2",
        "无 ADL 受限","有 ADL 受限","从不吸烟","曾/现吸烟","中/高教育","低教育")
VN <- c(年龄="A 年龄（岁）", BMI="B BMI（kg/m²）", `心血管代谢共病（高血压/糖尿病）`="C 心血管代谢共病（高血压/糖尿病）",
        `ADL 受限`="D ADL 受限", 吸烟="E 吸烟", 教育="F 受教育程度")
D <- R[分层变量 != "总体", .(分层变量, 层, 来源 = 队列, n, 事件, OR, 下限, 上限, p)]
P2 <- P[, .(分层变量, 层, 来源 = "四队列合并", n = NA_integer_, 事件 = NA_integer_,
            OR = 合并OR, 下限, 上限, p)]
A <- rbind(D, P2, fill = TRUE)
A[, 变量标签 := VN[分层变量]]
ORD <- LV[LV %in% unique(A$层)]
A[, 层 := factor(层, levels = ORD)]
CLV   <- c("CHARLS 中国","ELSA 英国","HRS 美国","SHARE 欧洲","四队列合并")
CLV_D <- TT(CLV, c("CHARLS (China)","ELSA (England)","HRS (USA)","SHARE (Europe)","Pooled (4 cohorts)"))
A[, 来源 := factor(来源, levels = CLV)]
## ---- 图内文字语言切换（中文版文章出中文图，英文版文章出英文图）----
LV_EN <- c("<65"="<65", "65–74"="65–74", "≥75"="≥75", "<25"="<25", "25–29.9"="25–29.9", "≥30"="≥30",
           "0"="0", "1"="1", "≥2"="≥2",
           "无 ADL 受限"="No ADL limitation", "有 ADL 受限"="Any ADL limitation",
           "从不吸烟"="Never smoker", "曾/现吸烟"="Former or current smoker",
           "中/高教育"="Medium or high education", "低教育"="Low education")
VN_EN <- c("A 年龄（岁）"="A Age (years)", "B BMI（kg/m²）"="B BMI (kg/m²)",
           "C 心血管代谢共病（高血压/糖尿病）"="C Cardiometabolic multimorbidity (hypertension or diabetes)",
           "D ADL 受限"="D ADL limitation", "E 吸烟"="E Smoking", "F 受教育程度"="F Education")
if (EN) {
  A[, 层 := factor(unname(LV_EN[as.character(层)]), levels = unname(LV_EN[ORD]))]
  A[, 变量标签 := unname(VN_EN[变量标签])]
  A[, 来源 := factor(unname(CLV_D[match(as.character(来源), CLV)]), levels = CLV_D)]
}
A[, 显著 := p < 0.05]
setorder(A, 变量标签, 层, 来源)

COL <- setNames(c("#C0504D","#2E74B5","#4E7B4E","#8064A2","#1F3A5F"), CLV_D)

g <- ggplot(A, aes(x = OR, y = 层, colour = 来源, shape = 来源)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey45") +
  geom_errorbarh(aes(xmin = 下限, xmax = 上限), height = 0.28, linewidth = 0.42,
                 position = position_dodge(width = 0.78)) +
  geom_point(aes(size = 来源), position = position_dodge(width = 0.78)) +
  scale_colour_manual(values = COL) +
  scale_shape_manual(values = c(16, 16, 16, 16, 18)) +
  scale_size_manual(values = c(1.7, 1.7, 1.7, 1.7, 2.8)) +
  scale_x_log10(breaks = c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.5)) +
  facet_grid(变量标签 ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(title = TT("Table 2  握力（每 +1 SD）与跌倒：四队列亚组分析与合并估计",
                  "Table 2  Grip strength and falls: subgroup analyses and pooled estimates"),
       subtitle = TT(paste0("亚组模型统一使用核心协变量（年龄、教育、吸烟、高血压、糖尿病、BMI、ADL；分层变量本身不纳入校正）；\n四队列按 DerSimonian–Laird 随机效应合并。所有 6 个分层变量的层间交互检验均不显著（P = ", IP_TXT_ZH, "）"),
                    paste0("All subgroup models use the same core covariates (age, education, smoking, hypertension,\ndiabetes, BMI, ADL; the stratifying variable itself is not adjusted for); cohorts are pooled by\nDerSimonian-Laird random effects. No stratum-by-stratum interaction test was significant for any\nof the 6 stratifying variables (P = ", IP_TXT_EN, ")")),
       x = TT("OR（每 +1 SD 握力，对数刻度）", "OR per 1-SD higher grip strength (log scale)"),
       y = NULL, colour = NULL, shape = NULL, size = NULL) +
  theme_minimal() +
  theme(plot.title = element_text(size = 13, colour = "#1F3A5F"),
        plot.subtitle = element_text(size = 8.6, colour = "#5A5A5A", lineheight = 1.15),
        strip.text.y.left = element_text(size = 9.5, colour = "#1F3A5F", angle = 0, hjust = 0),
        strip.placement = "outside",
        axis.text.y = element_text(size = 8.4),
        axis.text.x = element_text(size = 9),
        panel.spacing.y = unit(0.35, "lines"),
        legend.position = "bottom", legend.text = element_text(size = 9),
        panel.grid.minor = element_blank(), panel.grid.major.y = element_line(colour = "grey93"))

agg_png(file.path(RES, paste0("65_Table2_亚组森林图", FSFX, ".png")), width = 2600, height = 3000, res = 200, background = "white")
print(g)
dev.off()
cat("图已生成: res/65_Table2_亚组森林图", FSFX, ".png\n", sep = "")
cat("层数 =", nrow(P), " 逐队列估计 =", nrow(D), "\n")
cat("层p<0.05 数 =", sum(P$p < 0.05), "/", nrow(P), "\n")
