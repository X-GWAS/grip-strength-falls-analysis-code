## 86_RCS_握力与跌倒.R —— 握力—跌倒的剂量反应（限制性立方样条 RCS）
## 设计要点：
##   1) 四队列样本均为【女性】子样本（CHARLS 4724 / ELSA 8925 / HRS 13098 / SHARE 87099），
##      故 RCS 直接按绝对握力（kg）建模，无需再分性别；
##   2) 分析 A（临床曲线）：各队列用自己的握力分布取 4 结点（5/35/65/95 百分位），
##      参照点 = 该队列握力中位数，横轴 kg —— 直接回答「低于多少 kg 风险开始上升」；
##   3) 分析 B（可合并曲线）：把握力换算为队列内 z（(grip-均值)/队列SD），
##      四队列使用【同一组固定结点 z = -2.5/-1.5/-0.5/+0.5/+1.5】（5 结点），
##      参照 z = 0。只有结点位置完全一致，四队列的样条基函数才相同、系数才可合并；
##      合并用逐系数 DerSimonian–Laird 随机效应（对角 τ²，忽略队列内交叉协方差，已如实标注）；
##   4) 协变量与主分析一致：age + edu + smoke + hibp + diab + bmi + adl（+ factor(rnd)）；
##      全部模型经 svydesign + svyglm 设计校正（CHARLS 分层+PSU+权重；ELSA 仅权重；HRS/SHARE 分层+PSU+权重）。
suppressMessages({ library(data.table); library(survey); library(ggplot2); library(ragg) })
options(survey.lonely.psu = "adjust")
B   <- "/Users/mac/Documents/Codex/2026-08-31/xian/work"
T1  <- file.path(B, "外部复现-ELSA与HRS/t1")
RES <- file.path(B, "外部复现-ELSA与HRS/res")

## ---------- Harrell 截断幂参数化 RCS 基 ----------
rcs_basis <- function(x, knots) {
  k <- length(knots); out <- matrix(NA_real_, length(x), k - 1); out[, 1] <- x
  if (k > 2) for (j in seq_len(k - 2)) {
    out[, j + 1] <- pmax(x - knots[j], 0)^3 -
      pmax(x - knots[k - 1], 0)^3 * (knots[k] - knots[j]) / (knots[k] - knots[k - 1])
  }
  colnames(out) <- c("x", paste0("x", seq_len(k - 2))); out
}

## ---------- DL 随机效应（单变量） ----------
dl <- function(b, se) {
  w <- 1 / se^2; bf <- sum(w * b) / sum(w); Q <- sum(w * (b - bf)^2); df <- length(b) - 1
  C <- sum(w) - sum(w^2) / sum(w); tau2 <- if (df > 0) max(0, (Q - df) / C) else 0
  wr <- 1 / (se^2 + tau2); br <- sum(wr * b) / sum(wr); ser <- sqrt(1 / sum(wr))
  list(b = br, se = ser, lo = br - 1.96 * ser, hi = br + 1.96 * ser,
       p = 2 * pnorm(-abs(br / ser)), I2 = if (Q > 0 && df > 0) max(0, (Q - df) / Q) * 100 else 0,
       Q = Q, df = df, tau2 = tau2)
}

prep <- function(rds) { d <- as.data.table(readRDS(rds)); if ("ID" %in% names(d) && !"id" %in% names(d)) setnames(d, "ID", "id"); d }
coh <- list(
  CHARLS = list(d = prep(file.path(B, "新暴露-体力活动与握力轨迹/t1/样本_CHARLS.rds")),
                g = "grip_2011", ids = "communityID", strata = "stratum", wt = "wt3", lab = "CHARLS 中国", y = "fall_1518", rnd = FALSE),
  ELSA   = list(d = prep(file.path(T1, "样本_ELSA.rds")),   g = "grip", ids = "id",   strata = NULL,    wt = "wt", lab = "ELSA 英国",   y = "fall", rnd = TRUE),
  HRS    = list(d = prep(file.path(T1, "样本_HRS.rds")),    g = "grip", ids = "psu", strata = "strat", wt = "wt", lab = "HRS 美国",    y = "fall", rnd = TRUE),
  SHARE  = list(d = prep(file.path(T1, "样本_SHARE.rds")),  g = "grip", ids = "psu", strata = "strat", wt = "wt", lab = "SHARE 欧洲",  y = "fall", rnd = TRUE)
)
CORE <- c("age", "edu", "smoke", "hibp", "diab", "bmi", "adl")
ZK   <- c(-2.5, -1.5, -0.5, 0.5, 1.5)          # 分析 B 的四队列共用结点（z 尺度）

## ---------- 逐队列建模 ----------
fit_cohort <- function(cn, scale = c("kg", "z")) {
  scale <- match.arg(scale); z <- coh[[cn]]; d <- copy(z$d)
  d[, gz := (get(z$g) - mean(get(z$g), na.rm = TRUE)) / sd(get(z$g), na.rm = TRUE)]
  xv <- if (scale == "kg") z$g else "gz"
  keep <- c(z$y, xv, CORE, if (z$rnd) "rnd")
  dd <- d[complete.cases(d[, ..keep])]
  knots <- if (scale == "kg") as.numeric(quantile(dd[[xv]], c(.05, .35, .65, .95))) else ZK
  ref   <- if (scale == "kg") as.numeric(median(dd[[xv]])) else 0
  X <- rcs_basis(dd[[xv]], knots); dx <- cbind(as.data.frame(dd), as.data.frame(X))
  des <- if (is.null(z$strata)) {
    svydesign(ids = as.formula(paste0("~", z$ids)), weights = as.formula(paste0("~", z$wt)), data = dx)
  } else {
    svydesign(ids = as.formula(paste0("~", z$ids)), strata = as.formula(paste0("~", z$strata)),
              weights = as.formula(paste0("~", z$wt)), data = dx, nest = TRUE)
  }
  covf <- paste(c(CORE, if (z$rnd) "factor(rnd)"), collapse = " + ")
  nm <- colnames(X)                      # 4 结点 -> x,x1,x2；5 结点 -> x,x1,x2,x3
  m <- svyglm(as.formula(paste(z$y, "~", paste(nm, collapse = " + "), "+", covf)),
              design = des, family = quasibinomial())
  b <- coef(m)[nm]; V <- vcov(m)[nm, nm]
  nl <- nm[-1]; dnl <- length(nl)        # 非线性项个数 = 结点数 - 2
  rng   <- if (scale == "kg") as.numeric(quantile(dd[[xv]], c(.01, .99))) else as.numeric(quantile(dd[[xv]], c(.005, .995)))
  gridx <- seq(rng[1], rng[2], length.out = 400)
  Bx <- rcs_basis(c(ref, gridx), knots)
  gx <- Bx[-1, , drop = FALSE] - matrix(Bx[1, ], nrow(Bx) - 1, ncol(Bx), byrow = TRUE)
  rc <- as.numeric(gx %*% b); se <- sqrt(rowSums((gx %*% V) * gx))
  ## 非线性 Wald（结点项的联合检验，设计校正 vcov）
  W <- as.numeric(t(b[nl]) %*% solve(V[nl, nl]) %*% b[nl])
  pnl <- pchisq(W, dnl, lower.tail = FALSE)
  list(cn = cn, lab = z$lab, scale = scale, knots = knots, ref = ref, n = nrow(dd), ev = sum(dd[[z$y]]),
       grid = gridx,
       or = exp(rc), lo = exp(rc - 1.96 * se), hi = exp(rc + 1.96 * se),
       b = b, V = V, pnl = pnl, W = W, dnl = dnl, sdg = sd(dd[[z$g]]), medg = median(dd[[z$g]]),
       p_range = rng)
}

resA <- lapply(names(coh), fit_cohort, scale = "kg"); names(resA) <- names(coh)
resB <- lapply(names(coh), fit_cohort, scale = "z");  names(resB) <- names(coh)

## ---------- 分析 A：临床阈值表（相对各队列握力中位数） ----------
CUT <- c(14, 16, 18, 20, 22, 24, 26, 28)
thr <- rbindlist(lapply(resA, function(r) {
  cu <- CUT[abs(CUT - r$ref) > 1e-9]          # 参照点本身不重复列出
  Bx <- rcs_basis(c(r$ref, cu), r$knots)
  gx <- Bx[-1, , drop = FALSE] - matrix(Bx[1, ], length(cu), ncol(Bx), byrow = TRUE)
  rc <- as.numeric(gx %*% r$b); se <- sqrt(rowSums((gx %*% r$V) * gx))
  data.table(队列 = r$lab, n = r$n, 事件 = r$ev, 参照握力中位数 = round(r$ref, 1),
             握力kg = cu, OR = exp(rc), 下限 = exp(rc - 1.96 * se), 上限 = exp(rc + 1.96 * se),
             p = 2 * pnorm(-abs(rc / se)), 该值在数据范围内 = cu >= r$p_range[1] & cu <= r$p_range[2])
}))

## ---------- 分析 A：风险开始上升的握力（点估计上穿 1 的最低值） ----------
cross <- rbindlist(lapply(resA, function(r) {
  ord <- order(r$grid); g <- r$grid[ord]; o <- r$or[ord]; lo <- r$lo[ord]
  below <- g < r$ref - 0.05                             # 只看中位数以下（风险上升侧），排除参照点浮点伪影
  pick  <- function(cond) if (any(below & cond)) max(which(below & cond)) else NA_integer_
  ## 大样本下「CI 下限 > 1」在紧邻中位数处就已成立，不是有临床意义的阈值；
  ## 主阈值改用「OR ≥ 1.25 且 95%CI 下限 > 1」的最高握力（相对该队列中位数）
  i125  <- pick(o >= 1.25 & lo > 1); i150 <- pick(o >= 1.5)
  data.table(队列 = r$lab, n = r$n, 事件 = r$ev, 握力中位数 = round(r$medg, 1),
             主阈值kg_OR达1.25 = if (is.na(i125)) NA_real_ else round(g[i125], 1),
             该点OR = if (is.na(i125)) NA_real_ else round(o[i125], 3),
             该点CI下限 = if (is.na(i125)) NA_real_ else round(lo[i125], 3),
             该点CI上限 = if (is.na(i125)) NA_real_ else round(r$hi[ord][i125], 3),
             OR达1.5的握力kg = if (is.na(i150)) NA_real_ else round(g[i150], 1))
}))

## ---------- 分析 B：共用结点的系数 → 可合并 ----------
BB <- rbindlist(lapply(resB, function(r) data.table(队列 = r$lab, n = r$n, 事件 = r$ev,
  项 = names(r$b), 系数 = as.numeric(r$b), SE = sqrt(diag(r$V)))))
BB[, `:=`(OR每单位 = exp(系数), 下限 = exp(系数 - 1.96*SE), 上限 = exp(系数 + 1.96*SE))]

## ---------- 分析 B：逐 z 点 DL 合并 → 合并剂量反应曲线 ----------
ZG <- seq(-3, 2, by = 0.25)
pool_at <- function(zv) {
  o <- rbindlist(lapply(resB, function(r) {
    Bx <- rcs_basis(c(0, zv), r$knots)
    gx <- matrix(Bx[2, ] - Bx[1, ], nrow = 1)
    rc <- as.numeric(gx %*% r$b); se <- sqrt(as.numeric(gx %*% r$V %*% t(gx)))
    data.table(队列 = r$lab, b = rc, se = se)
  }))
  ## 参照点（z = 0）各队列对比恒为 0、SE 恒为 0，DL 权重为 Inf → 直接给定值 1
  if (abs(zv) < 1e-9) return(data.table(z = 0, OR = 1, 下限 = 1, 上限 = 1, p = 1, I2 = 0, k = nrow(o)))
  p <- dl(o$b, o$se)
  data.table(z = zv, OR = exp(p$b), 下限 = exp(p$lo), 上限 = exp(p$hi), p = p$p, I2 = p$I2, k = nrow(o))
}
PC <- rbindlist(lapply(ZG, pool_at))

## ---------- 分析 B：合并非线性检验（对角 τ² 的多变量 DL） ----------
P  <- length(resB[[1]]$b)                                  # 系数个数
tau2j <- sapply(1:P, function(j) dl(sapply(resB, function(r) r$b[j]), sqrt(sapply(resB, function(r) r$V[j, j])))$tau2)
Tm <- diag(tau2j)
Ai <- lapply(resB, function(r) solve(r$V + Tm))
A  <- Reduce(`+`, Ai)
bm <- solve(A, Reduce(`+`, lapply(seq_along(resB), function(i) Ai[[i]] %*% resB[[i]]$b)))
Vm <- solve(A)
idx <- 2:P
Wnl <- as.numeric(t(bm[idx]) %*% solve(Vm[idx, idx]) %*% bm[idx])
Pnl <- pchisq(Wnl, length(idx), lower.tail = FALSE)

## ---------- 落盘 ----------
fwrite(thr,   file.path(RES, "80_RCS_临床阈值表.csv"))
fwrite(cross, file.path(RES, "81_RCS_上穿1的握力.csv"))
fwrite(BB,    file.path(RES, "82_RCS_共用结点系数.csv"))
fwrite(PC,    file.path(RES, "83_RCS_合并曲线.csv"))
NL <- rbindlist(lapply(resB, function(r) data.table(队列 = r$lab, 尺度 = "队列内 z（共用结点）",
       Wald = round(r$W, 3), 自由度 = r$dnl, P非线性 = r$pnl)))
NL <- rbind(NL, rbindlist(lapply(resA, function(r) data.table(队列 = r$lab, 尺度 = "绝对 kg（本队列结点）",
       Wald = round(r$W, 3), 自由度 = r$dnl, P非线性 = r$pnl))))
NL <- rbind(NL, data.table(队列 = "四队列合并", 尺度 = "队列内 z（共用结点）", Wald = round(Wnl, 3),
       自由度 = length(idx), P非线性 = Pnl))
fwrite(NL, file.path(RES, "84_RCS_非线性检验.csv"))
saveRDS(list(A = resA, B = resB, PC = PC, tau2j = tau2j, Wnl = Wnl, dfnl = length(idx),
             Pnl = Pnl, bm = bm, Vm = Vm),
        file.path(RES, "85_RCS_对象.rds"))

cat("\n=========== 分析 A：绝对握力（kg）非线性检验 ===========\n"); print(resA |> lapply(\(r) c(队列 = r$lab, P = r$pnl)))
cat("\n=========== 分析 B：z 尺度非线性检验 ===========\n"); print(resB |> lapply(\(r) c(队列 = r$lab, P = r$pnl)))
cat(sprintf("\n合并非线性 Wald = %.3f, df = %d, P = %.4g\n", Wnl, length(idx), Pnl))
cat("\n=========== 风险上穿 1 的握力 ===========\n"); print(cross)
cat("\n=========== 阈值表（部分）===========\n"); print(thr[握力kg %in% c(16, 18, 20, 24)][order(队列, -握力kg)])
cat("\n完成：res/80–85_RCS_*.csv/rds\n")
