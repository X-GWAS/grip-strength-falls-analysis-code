## 84_协变量口径分解.R —— 为什么「统一协变量口径」把合并 I² 从 25% 推到 80%？
## 做法：以统一核心协变量为基准模型，逐个补回该队列主模型里有、而其他队列没有的协变量，
##       看每一步点估计移动多少。用于在正文/补充材料中交代口径敏感性。
suppressMessages({ library(data.table); library(survey) })
options(survey.lonely.psu = "adjust")
B  <- "/Users/mac/Documents/Codex/2026-08-31/xian/work"
T1 <- file.path(B, "外部复现-ELSA与HRS/t1")
RES<- file.path(B, "外部复现-ELSA与HRS/res")
prep <- function(r) { d <- as.data.table(readRDS(r)); if ("ID" %in% names(d)) setnames(d, "ID", "id"); d }

specs <- list(
  CHARLS = list(d = prep(file.path(B, "新暴露-体力活动与握力轨迹/t1/样本_CHARLS.rds")),
    ids="communityID", strata="stratum", wt="wt3", y="fall_1518", rnd=FALSE, lab="CHARLS 中国",
    core = c("age","edu","smoke","hibp","diab","bmi","adl"),
    extra = c("married","rural","drink","cesd"),
    main = "age3 + factor(edu) + married + rural + smoke + drink + hibp + diab + bmi + adl + cesd"),
  ELSA = list(d = prep(file.path(T1,"样本_ELSA.rds")), ids="id", strata=NULL, wt="wt", y="fall", rnd=TRUE, lab="ELSA 英国",
    core = c("age","edu","smoke","hibp","diab","bmi","adl"),
    extra = c("drink","mstat","hearte","stroke","cesd"),
    main = "age + factor(edu) + mstat + smoke + drink + hibp + diab + hearte + stroke + bmi + adl + cesd + factor(rnd)"),
  HRS = list(d = prep(file.path(T1,"样本_HRS.rds")), ids="psu", strata="strat", wt="wt", y="fall", rnd=TRUE, lab="HRS 美国",
    core = c("age","edu","smoke","hibp","diab","bmi","adl"),
    extra = c("iadl"),
    main = "age + factor(edu) + smoke + hibp + diab + bmi + adl + iadl + factor(rnd)"),
  SHARE = list(d = prep(file.path(T1,"样本_SHARE.rds")), ids="psu", strata="strat", wt="wt", y="fall", rnd=TRUE, lab="SHARE 欧洲",
    core = c("age","edu","smoke","hibp","diab","bmi","adl"),
    extra = c("mstat","hearte","stroke","eurod"),
    main = "age + factor(edu) + mstat + smoke + hibp + diab + hearte + stroke + bmi + adl + eurod + factor(rnd)")
)
out <- list()
for (cn in names(specs)) {
  z <- specs[[cn]]; d <- copy(z$d)
  fc <- function(v) if (v == "edu") "factor(edu)" else v
  cov_core <- paste(c(vapply(z$core, fc, ""), if (z$rnd) "factor(rnd)"), collapse = " + ")
  need <- unique(c(z$y, "grip_z", z$core, z$extra, "age3", if (z$rnd) "rnd"))
  need <- intersect(need, names(d))
  dd <- d[complete.cases(d[, ..need])]
  des <- if (is.null(z$strata)) svydesign(ids=as.formula(paste0("~",z$ids)), weights=as.formula(paste0("~",z$wt)), data=as.data.frame(dd))
         else svydesign(ids=as.formula(paste0("~",z$ids)), strata=as.formula(paste0("~",z$strata)), weights=as.formula(paste0("~",z$wt)), data=as.data.frame(dd), nest=TRUE)
  fit <- function(covs, tag) {
    m <- tryCatch(svyglm(as.formula(paste(z$y, "~ grip_z +", covs)), design=des, family=quasibinomial()), error=function(e) NULL)
    if (is.null(m)) return(NULL)
    c2 <- summary(m)$coefficients["grip_z", ]
    cat(sprintf("%-12s %-46s OR=%.3f (%.3f-%.3f) p=%.3g n=%d\n", z$lab, tag,
        exp(c2[1]), exp(c2[1]-1.96*c2[2]), exp(c2[1]+1.96*c2[2]), c2[4], length(m$y)))
    data.table(队列=z$lab, 模型=tag, n=length(m$y), 事件=sum(m$y),
               OR=exp(c2[1]), 下限=exp(c2[1]-1.96*c2[2]), 上限=exp(c2[1]+1.96*c2[2]), p=c2[4])
  }
  cat("\n=====", z$lab, "=====\n")
  out[[paste(cn,"主模型")]]   <- fit(z$main, "主模型（本队列协变量）")
  out[[paste(cn,"核心")]]     <- fit(cov_core, "统一核心协变量")
  for (v in z$extra) out[[paste(cn,v)]] <- fit(paste0(cov_core, " + ", fc(v)), paste0("核心 + ", v))
  out[[paste(cn,"全补")]]     <- fit(paste0(cov_core, " + ", paste(vapply(z$extra, fc, ""), collapse=" + ")), "核心 + 全部本地协变量")
}
R <- rbindlist(out)
fwrite(R, file.path(RES, "64_协变量口径分解.csv"))

## 合并对比：主模型口径 vs 统一核心口径
pf <- function(OR, LO, HI) {
  b <- log(OR); se <- (log(HI)-log(LO))/(2*1.96); w <- 1/se^2; bf <- sum(w*b)/sum(w)
  Q <- sum(w*(b-bf)^2); df <- length(b)-1; C <- sum(w)-sum(w^2)/sum(w); tau2 <- max(0,(Q-df)/C)
  wr <- 1/(se^2+tau2); br <- sum(wr*b)/sum(wr); ser <- sqrt(1/sum(wr))
  sprintf("OR=%.3f (%.3f-%.3f) p=%.3g  I2=%.1f%%  Q_p=%.3f", exp(br), exp(br-1.96*ser), exp(br+1.96*ser),
          2*pnorm(-abs(br/ser)), if (Q>0) max(0,(Q-df)/Q)*100 else 0, pchisq(Q,df,lower.tail=FALSE))
}
cat("\n主模型口径合并  :", pf(R[模型=="主模型（本队列协变量）"]$OR, R[模型=="主模型（本队列协变量）"]$下限, R[模型=="主模型（本队列协变量）"]$上限), "\n")
cat("统一核心口径合并:", pf(R[模型=="统一核心协变量"]$OR, R[模型=="统一核心协变量"]$下限, R[模型=="统一核心协变量"]$上限), "\n")
cat("核心+全部本地协变量合并:", pf(R[模型=="核心 + 全部本地协变量"]$OR, R[模型=="核心 + 全部本地协变量"]$下限, R[模型=="核心 + 全部本地协变量"]$上限), "\n")
cat("\n完成: res/64_协变量口径分解.csv\n")
