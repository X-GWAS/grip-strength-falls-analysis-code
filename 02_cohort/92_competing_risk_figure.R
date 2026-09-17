## 92_竞争风险_图.R —— 死亡竞争风险：五套口径的逐队列与合并估计
suppressMessages({ library(data.table); library(ggplot2); library(ragg) })
source("/Users/mac/Documents/Codex/2026-08-31/xian/work/tools/fig_lang.R")   # 中英图文字开关
RES <- "/Users/mac/Documents/Codex/2026-08-31/xian/work/外部复现-ELSA与HRS/res"
R <- fread(file.path(RES, "70_竞争风险_逐队列.csv"))
P <- fread(file.path(RES, "72_竞争风险_合并.csv"))
R[, 来源 := 队列]; P[, 来源 := "四队列合并"]
A <- rbind(R[, .(口径, 来源, OR, 下限, 上限, p)], P[, .(口径, 来源, 合并OR, 下限, 上限, p)],
           use.names = FALSE)
setnames(A, c("口径","来源","OR","下限","上限","p"))
A[, 口径 := factor(口径, levels = rev(c("A 原因别风险（仅存活观察者）",
                                        "B 亚分布（死亡记为未跌倒）",
                                        "C 最坏情形（死亡记为跌倒）",
                                        "D 全部缺失记为未跌倒（乐观）",
                                        "E 全部缺失记为跌倒（悲观）")))]
CLV   <- c("CHARLS 中国","ELSA 英国","HRS 美国","SHARE 欧洲","四队列合并")
CLV_D <- TT(CLV, c("CHARLS (China)","ELSA (England)","HRS (USA)","SHARE (Europe)","Pooled (4 cohorts)"))
A[, 来源 := factor(来源, levels = CLV)]
## ---- 图内文字语言切换 ----
KM <- c("A 原因别风险（仅存活观察者）"="A Cause-specific (surviving observed participants only)",
        "B 亚分布（死亡记为未跌倒）"="B Subdistribution (death counted as no fall)",
        "C 最坏情形（死亡记为跌倒）"="C Worst case (death counted as a fall)",
        "D 全部缺失记为未跌倒（乐观）"="D All missing counted as no fall (optimistic)",
        "E 全部缺失记为跌倒（悲观）"="E All missing counted as a fall (pessimistic)")
if (EN) {
  A[, 口径 := factor(unname(KM[as.character(口径)]), levels = unname(KM[rev(levels(口径))])) ]
  A[, 来源 := factor(unname(CLV_D[match(as.character(来源), CLV)]), levels = CLV_D)]
}
COL <- setNames(c("#C0504D","#2E74B5","#4E7B4E","#8064A2","#1F3A5F"), CLV_D)
g <- ggplot(A, aes(x = OR, y = 口径, colour = 来源, shape = 来源)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey45") +
  geom_errorbarh(aes(xmin = 下限, xmax = 上限), height = 0.3, linewidth = 0.5,
                 position = position_dodge(width = 0.75)) +
  geom_point(aes(size = 来源), position = position_dodge(width = 0.75)) +
  scale_colour_manual(values = COL) + scale_shape_manual(values = c(16,16,16,16,18)) +
  scale_size_manual(values = c(2.0,2.0,2.0,2.0,3.2)) +
  scale_x_log10(breaks = c(0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0)) +
  labs(title = TT("死亡竞争风险敏感性：握力（每 +1 SD）与跌倒，五套口径",
                  "Sensitivity to death as a competing risk: five treatments of missing outcomes"),
       subtitle = TT("四套口径的协变量、权重、设计完全相同，只改变「结局波之前死亡的人怎么处理」。\n五套口径的合并估计全部 <1；即使把全部缺失（死亡＋失访）都当作未跌倒，合并 OR 仍为 0.892（0.803–0.991）",
                     "Covariates, weights and the survey design are identical across all five treatments; only the handling of\nparticipants who died before the outcome wave changes. The pooled estimate stays below 1 under every\ntreatment: even when all missing outcomes (death plus attrition) are counted as no fall, it is 0.892 (0.803-0.991)"),
       x = TT("OR（每 +1 SD 握力，对数刻度）", "OR per 1-SD higher grip strength (log scale)"),
       y = NULL, colour = NULL, shape = NULL, size = NULL) +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, colour = "#1F3A5F"),
        plot.subtitle = element_text(size = 9, colour = "#5A5A5A", lineheight = 1.2),
        axis.text.y = element_text(size = 10), axis.text.x = element_text(size = 10),
        legend.position = "bottom", legend.text = element_text(size = 9.5),
        panel.grid.minor = element_blank())
agg_png(file.path(RES, paste0("73_竞争风险_五口径对比", FSFX, ".png")), width = 2600, height = 1700, res = 200, background = "white")
print(g); dev.off()
cat("图已生成: res/73_竞争风险_五口径对比", FSFX, ".png\n", sep = "")
