## 00_deff_helper.R —— 设计效应（design effect）统一计算
## 用法：source 本文件后调用 deff_block(des, data, formula, term, label)
##   对比三种方差来源：
##     1) 完整抽样设计（分层 + PSU + 权重）——本文主结果
##     2) 简单随机抽样（无权重的普通 glm）
##     3) 仅加权（ids = ~1）
##   DEFF_总体 = (SE_完整设计 / SE_简单随机)^2
##   DEFF_加权 = (SE_仅加权 / SE_简单随机)^2
##   DEFF_聚类 = (SE_完整设计 / SE_仅加权)^2   （剥离加权后的聚类部分）
##   有效样本量 = n / DEFF_总体

## 传进来的是【该模型的分析数据框】（已做完整案例筛选），函数内部构造三种方差来源
survey_deff <- function(dat, formula, term, label, wvar, psu = NULL, strata = NULL,
                        family = quasibinomial(), n_events = NULL, 事件列 = NULL) {
  dat <- as.data.frame(dat)
  dat <- dat[!is.na(dat[[wvar]]) & dat[[wvar]] > 0, ]
  if (!is.null(psu))    dat <- dat[!is.na(dat[[psu]]), ]
  if (!is.null(strata)) dat <- dat[!is.na(dat[[strata]]), ]
  f <- as.formula(formula)
  ## 【2026-09-17 修正】先按模型变量做完整案例筛选再建抽样设计，
  ##   否则设计对象建立在一份"筛后但未去缺失"的数据上，会与真正进入模型的样本错位
  ##   （表现为 n 与 PSU 数不一致：ELSA 曾出现 n=8925 而 PSU 数=5109，实际人数为 4153）。
  mvars <- all.vars(f)
  dat <- dat[stats::complete.cases(dat[, mvars, drop = FALSE]), , drop = FALSE]
  ## 【2026-09-17 修正】原来要求 PSU 与分层必须同时存在才算"含聚类"，
  ##   导致 ELSA（有参与者 id 可作最后阶段聚类单位、但无分层变量）被误判为"仅权重"，
  ##   设计效应表的 SE 与主模型（ids = "id"）不一致。现允许只有 PSU、没有分层。
  has_cluster <- !is.null(psu)
  des <- if (has_cluster && !is.null(strata))
    svydesign(ids = as.formula(paste0("~", psu)), strata = as.formula(paste0("~", strata)),
              weights = as.formula(paste0("~", wvar)), data = dat, nest = TRUE)
  else if (has_cluster)
    svydesign(ids = as.formula(paste0("~", psu)),
              weights = as.formula(paste0("~", wvar)), data = dat)
  else svydesign(ids = ~1, weights = as.formula(paste0("~", wvar)), data = dat)
  des_w <- svydesign(ids = ~1, weights = as.formula(paste0("~", wvar)), data = dat)
  suppressWarnings({
    m_svy <- svyglm(f, design = des,   family = family)
    m_wt  <- svyglm(f, design = des_w, family = family)
    m_glm <- glm(f, data = dat, family = family)
  })
  cs <- summary(m_svy)$coefficients
  if (!term %in% rownames(cs)) return(NULL)
  b <- cs[term, 1]; se <- cs[term, 2]
  b_srs <- coef(m_glm)[[term]]; se_srs <- summary(m_glm)$coefficients[term, 2]
  cs_w <- summary(m_wt)$coefficients
  b_wt <- cs_w[term, 1]; se_wt <- cs_w[term, 2]

  deff_tot <- (se / se_srs)^2
  deff_wt  <- (se_wt / se_srs)^2
  deff_cl  <- (se / se_wt)^2
  ## 关键：报告【模型实际使用】的样本量，而不是"符合条件但协变量有缺失"的样本量
  n <- length(m_glm$y)
  y_used <- m_glm$y
  ## ---- 独立核验：Kish 权重设计效应（只与权重分布有关，与模型无关）----
  ## 2026-09-17 修正：w 必须取【模型实际使用】的权重（m_svy$prior.weights），
  ## 否则 n 与 w 长度不一致会把"筛后样本"的权重分布当成模型样本的分布。
  ## 注意必须取 m_svy（含抽样权重）而不是 m_glm（无权重，prior.weights 恒为 1）。
  w <- if (!is.null(m_svy$prior.weights) && length(m_svy$prior.weights) == length(m_svy$y))
         as.numeric(m_svy$prior.weights) else dat[[wvar]]
  kish  <- length(w) * sum(w^2) / sum(w)^2  # DEFF_w = 1 + CV(w)^2
  cv_w  <- sqrt(kish - 1)
  data.table(
    项 = label, n = n,
    事件 = if (!is.null(事件列)) sum(y_used == 1, na.rm = TRUE) else n_events,
    OR_完整设计 = exp(b), 下限 = exp(b - 1.96 * se), 上限 = exp(b + 1.96 * se),
    SE_完整设计 = se, p_完整设计 = cs[term, 4],
    OR_仅加权 = exp(b_wt), SE_仅加权 = se_wt, p_仅加权 = cs_w[term, 4],
    OR_简单随机 = exp(b_srs), SE_简单随机 = se_srs, p_简单随机 = summary(m_glm)$coefficients[term, 4],
    DEFF_总体 = deff_tot, DEFF_加权 = deff_wt, DEFF_聚类 = deff_cl,
    权重CV = cv_w, DEFF_权重Kish = kish,
    有效样本量 = round(n / deff_tot),
    含聚类 = has_cluster,
    有分层 = !is.null(strata),
     PSU数 = if (has_cluster) length(unique(des$cluster[, 1])) else NA_integer_,
    层数 = if (is.null(des$strata)) NA_integer_ else length(unique(des$strata[[1]]))
  )
}
