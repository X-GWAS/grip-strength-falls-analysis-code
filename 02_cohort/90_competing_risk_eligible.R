## 90_重建符合条件人群.R —— 把"因死亡/失访而未观察到结局"的人重新拉回分析
## 背景：四个队列现有的分析样本都只保留「结局可观察到」的人，死亡者在建样本时已被静默排除。
##       要做竞争风险，必须先重建「符合条件人群」＝女性、满足年龄、暴露波握力非缺失、有暴露波权重，
##       然后按结局波状态分三类：观察到结局（obs）／结局前死亡（died）／其他缺失（othermiss）。
## 死亡判定：
##   CHARLS/ELSA/SHARE：结局波 rXiwstat ∈ {5,6}（5=本波死亡、6=上一波死亡）；3=首波前死亡
##   HRS：无 iwstat，改用 Cross-Wave Tracker 的 EXDEATHYR/EXDEATHMO 与结局波年份比较
suppressMessages({ source("/Users/mac/.codex/knowledge/数据库资产/db_paths.R"); library(haven); library(data.table) })
OUT <- "/Users/mac/Documents/Codex/2026-08-31/xian/work/竞争风险死亡/t1"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
## 取变量；变量名不存在时返回 NA 向量而不是 NULL
## （教训：data.table 里只要有 1 个长度为 0 的列，整张表会静默变成 0 行，
##   HRS 的慢病变量实际叫 rXrxhibp / rXrxdiab，写成 hibpe/diabe 就会中招）
gv <- function(dt, v) if (v %in% names(dt)) dt[[v]] else rep(NA_real_, nrow(dt))

## ---------- 1. CHARLS（暴露 w3=2015，结局 w4=2018）----------
cat("\n===== CHARLS =====\n")
cl <- c("ID","communityID","ragender","raeduc_c","r3agey","r3gripsum","r4iwstat",
        "r3mbmi","r3smokev","r3hibpe","r3diabe","r3adlwa","r3wtresp")
C <- as.data.table(read_dta(DB$file$charls_h, col_select = any_of(cl)))
## 抽样设计信息在 2011 的 psu.dta 里（省 / 城乡），Harmonized 文件没有
PSU <- unique(as.data.table(read_dta(file.path(DB$charls, "2011charls/psu.dta"),
                                    col_select = any_of(c("communityID","province","urban_nbs")))))
C <- merge(C, PSU, by = "communityID", all.x = TRUE)
## 跌倒题在原始问卷里（Harmonized CHARLS 无此项），2018 波变量名 da023_w4
FF <- as.data.table(read_dta(file.path(DB$charls, "2018charls/Health_Status_and_Functioning.dta"),
                             col_select = any_of(c("ID", "da023_w4"))))
setnames(FF, "da023_w4", "fail1518")
C <- merge(C, FF, by = "ID", all.x = TRUE)
C[, `:=`(female = as.integer(ragender == 2),
         age = as.numeric(r3agey),
         grip = as.numeric(r3gripsum),
         edu = as.numeric(raeduc_c),
         bmi = as.numeric(r3mbmi),
         smoke = as.integer(r3smokev == 1),
         hibp = as.integer(r3hibpe == 1),
         diab = as.integer(r3diabe == 1),
         adl = as.numeric(r3adlwa),
         wt = as.numeric(r3wtresp),
         fall = fifelse(is.na(fail1518), NA_integer_, as.integer(fail1518 == 1)),
         iws = as.numeric(r4iwstat))]
C[grip <= 0 | grip > 100, grip := NA_real_]
E <- C[female == 1 & !is.na(age) & age >= 45 & !is.na(grip) & !is.na(wt) & wt > 0]
E[, `:=`(died = as.integer(iws %in% c(5, 6)),
         obs  = as.integer(!is.na(fall)),
         othermiss = as.integer(is.na(fall) & !(iws %in% c(5, 6))))]
E[, `:=`(rnd = 1L, grip_wave = 3L, fall_wave = 4L, psu = communityID,
         strat = paste(province, urban_nbs, sep = "_"))]
cat(sprintf("符合条件 %d 人；观察到结局 %d；结局前死亡 %d；其他缺失 %d\n",
            nrow(E), sum(E$obs), sum(E$died), sum(E$othermiss)))
print(E[obs == 0, .N, by = .(iws)])
saveRDS(E[, .(id = ID, rnd, grip_wave, fall_wave, grip, fall, age, edu, smoke, hibp, diab, bmi, adl,
              wt, psu, strat, died, obs, othermiss)], file.path(OUT, "elig_CHARLS.rds"))

## ---------- 2. ELSA（暴露 w2/4/6/8，结局 w3/5/7/9）----------
cat("\n===== ELSA =====\n")
gw <- c(2,4,6,8); fw <- c(3,5,7,9)
el <- c("idauniq","ragender", paste0("r", gw, "agey"), paste0("r", gw, "lgrip1"), paste0("r", gw, "rgrip1"),
        paste0("r", fw, "fall"), paste0("r", fw, "iwstat"), "raeduc_e",
        paste0("r", gw, "mbmi"), paste0("r", gw, "smokev"), paste0("r", gw, "hibpe"),
        paste0("r", gw, "diabe"), paste0("r", gw, "adlwa"), paste0("r", gw, "cwtresp"))
E0 <- as.data.table(read_dta(DB$file$elsa_h, col_select = any_of(el)))
EL <- rbindlist(lapply(seq_along(gw), function(k) {
  g <- gw[k]; f <- fw[k]
  a <- suppressWarnings(as.numeric(E0[[paste0("r", g, "lgrip1")]]))
  b <- suppressWarnings(as.numeric(E0[[paste0("r", g, "rgrip1")]]))
  a[a <= 0 | a > 100] <- NA; b[b <= 0 | b > 100] <- NA
  gr <- rowMeans(cbind(a, b), na.rm = TRUE); gr[is.nan(gr)] <- NA
  iws <- as.numeric(E0[[paste0("r", f, "iwstat")]])
  fv  <- fifelse(is.na(E0[[paste0("r", f, "fall")]]), NA_integer_, as.integer(E0[[paste0("r", f, "fall")]] == 1))
  data.table(id = E0$idauniq, rnd = k, grip_wave = g, fall_wave = f, grip = gr, fall = fv,
             age = as.numeric(E0[[paste0("r", g, "agey")]]), female = as.integer(E0$ragender == 2),
             edu = as.numeric(E0$raeduc_e), bmi = as.numeric(E0[[paste0("r", g, "mbmi")]]),
             smoke = as.integer(E0[[paste0("r", g, "smokev")]] == 1),
             hibp = as.integer(E0[[paste0("r", g, "hibpe")]] == 1),
             diab = as.integer(E0[[paste0("r", g, "diabe")]] == 1),
             adl = as.numeric(E0[[paste0("r", g, "adlwa")]]),
             wt = as.numeric(E0[[paste0("r", g, "cwtresp")]]), iws = iws)
}))
EL[grip <= 0 | grip > 100, grip := NA_real_]
EL <- EL[female == 1 & !is.na(age) & age >= 50 & !is.na(grip) & !is.na(wt) & wt > 0]
EL[, `:=`(died = as.integer(iws %in% c(5, 6)), obs = as.integer(!is.na(fall)),
          othermiss = as.integer(is.na(fall) & !(iws %in% c(5, 6))))]
cat(sprintf("ELSA 符合条件窗口 %d；观察到结局 %d；结局前死亡 %d；其他缺失 %d\n",
            nrow(EL), sum(EL$obs), sum(EL$died), sum(EL$othermiss)))
print(EL[, .(符合条件 = .N, 观测 = sum(obs), 死亡 = sum(died), 其他缺失 = sum(othermiss)), by = .(窗口 = paste0(grip_wave, "→", fall_wave))])
saveRDS(EL[, .(id, rnd, grip_wave, fall_wave, grip, fall, age, edu, smoke, hibp, diab, bmi, adl,
               wt, died, obs, othermiss)], file.path(OUT, "elig_ELSA.rds"))

## ---------- 3. SHARE（暴露 w1/2/4/5/6，结局 w2/4/5/6/7）----------
cat("\n===== SHARE =====\n")
gw2 <- c(1,2,4,5,6); fw2 <- c(2,4,5,6,7)
sh <- c("mergeid","ragender","raeducl","country", paste0("r", gw2, "agey"),
        paste0("r", gw2, "lgrip1"), paste0("r", gw2, "rgrip1"),
        paste0("r", fw2, "fall_s"), paste0("r", fw2, "iwstat"),
        paste0("r", gw2, "bmi"), paste0("r", gw2, "smokev"), paste0("r", gw2, "hibpe"),
        paste0("r", gw2, "diabe"), paste0("r", gw2, "adlwa"), paste0("r", gw2, "wtresp"))
S0 <- as.data.table(read_dta(DB$file$share_h, col_select = any_of(sh)))
cat("SHARE 读入变量:", paste(intersect(sh, names(S0)), collapse=", "), "\n")
SH <- rbindlist(lapply(seq_along(gw2), function(k) {
  g <- gw2[k]; f <- fw2[k]
  a <- suppressWarnings(as.numeric(S0[[paste0("r", g, "lgrip1")]]))
  b <- suppressWarnings(as.numeric(S0[[paste0("r", g, "rgrip1")]]))
  a[a <= 0 | a > 100] <- NA; b[b <= 0 | b > 100] <- NA
  gr <- rowMeans(cbind(a, b), na.rm = TRUE); gr[is.nan(gr)] <- NA
  fcol <- paste0("r", f, "fall_s"); icol <- paste0("r", f, "iwstat")
  fv <- if (fcol %in% names(S0)) fifelse(is.na(S0[[fcol]]), NA_integer_, as.integer(S0[[fcol]] == 1)) else NA_integer_
  iws <- if (icol %in% names(S0)) as.numeric(S0[[icol]]) else NA_real_
  data.table(id = S0$mergeid, rnd = k, grip_wave = g, fall_wave = f, grip = gr, fall = fv,
             age = as.numeric(S0[[paste0("r", g, "agey")]]), female = as.integer(S0$ragender == 2),
             edu = as.numeric(S0$raeducl), bmi = as.numeric(S0[[paste0("r", g, "bmi")]]),
             smoke = as.integer(S0[[paste0("r", g, "smokev")]] == 1),
             hibp = as.integer(S0[[paste0("r", g, "hibpe")]] == 1),
             diab = as.integer(S0[[paste0("r", g, "diabe")]] == 1),
             adl = as.numeric(S0[[paste0("r", g, "adlwa")]]),
             wt = as.numeric(S0[[paste0("r", g, "wtresp")]]), iws = iws)
}))
SH[grip <= 0 | grip > 100, grip := NA_real_]
SH <- SH[female == 1 & !is.na(age) & age >= 50 & !is.na(grip) & !is.na(wt) & wt > 0]
SH[, `:=`(died = as.integer(iws %in% c(5, 6)), obs = as.integer(!is.na(fall)),
          othermiss = as.integer(is.na(fall) & !(iws %in% c(5, 6))))]
cat(sprintf("SHARE 符合条件窗口 %d；观察到结局 %d；结局前死亡 %d；其他缺失 %d\n",
            nrow(SH), sum(SH$obs), sum(SH$died), sum(SH$othermiss)))
print(SH[, .(符合条件 = .N, 观测 = sum(obs), 死亡 = sum(died), 其他缺失 = sum(othermiss)), by = .(窗口 = paste0(grip_wave, "→", fall_wave))])
saveRDS(SH[, .(id, rnd, grip_wave, fall_wave, grip, fall, age, edu, smoke, hibp, diab, bmi, adl,
               wt, died, obs, othermiss)], file.path(OUT, "elig_SHARE.rds"))

## ---------- 4. HRS（暴露 w8-13，结局 w9-14；用 Tracker 死亡日期）----------
cat("\n===== HRS =====\n")
gw3 <- 8:13; fw3 <- 9:14
wave_year <- function(w) 1990 + 2 * w
hr <- c("hhidpn","raeducl",
        paste0("r", gw3, "lgrip1"), paste0("r", gw3, "rgrip1"),
        paste0("r", fw3, "fall"), paste0("r", gw3, "mbmi"), paste0("r", gw3, "smokef"),
        paste0("r", gw3, "rxhibp"), paste0("r", gw3, "rxdiab"), paste0("r", gw3, "adlfive"),
        paste0("r", gw3, "nwtresp"))
H0 <- as.data.table(read_dta(DB$file$hrs_h, col_select = any_of(hr)))
TRK <- as.data.table(read_dta(file.path(DB$hrs, "HRS_美国/Raw_data/Cross-Wave Tracker File/trk2020tr_r.dta"),
                              col_select = any_of(c("HHID","PN","GENDER","BIRTHYR","SECU","STRATUM",
                                                    "EXDEATHYR","EXDEATHMO"))))
TRK[, key := as.numeric(HHID) * 1000 + as.numeric(PN)]
TRK <- TRK[!is.na(key) & !duplicated(key)]
TRK <- TRK[, .(key, GENDER, BIRTHYR, SECU, STRATUM,
               EXDEATHYR = as.numeric(EXDEATHYR), EXDEATHMO = as.numeric(EXDEATHMO))]
cat("Tracker 行数:", nrow(TRK), " 有死亡年月:", sum(!is.na(TRK$EXDEATHYR)), "\n")
cat("主键交集:", sum(H0$hhidpn %in% TRK$key), "/", nrow(H0), "\n")
H0 <- merge(H0, TRK, by.x = "hhidpn", by.y = "key", all.x = TRUE)
cat(sprintf("  [诊断] 合并后 H0 行数 %d；列数 %d；GENDER 非缺失 %d；gw3 长度 %d\n",
            nrow(H0), ncol(H0), sum(!is.na(H0$GENDER)), length(gw3)))
HR <- rbindlist(lapply(seq_along(gw3), function(k) {
  g <- gw3[k]; f <- fw3[k]
  a <- suppressWarnings(as.numeric(H0[[paste0("r", g, "lgrip1")]]))
  b <- suppressWarnings(as.numeric(H0[[paste0("r", g, "rgrip1")]]))
  a[a <= 0 | a > 100] <- NA; b[b <= 0 | b > 100] <- NA
  gr <- rowMeans(cbind(a, b), na.rm = TRUE); gr[is.nan(gr)] <- NA
  if (k == 1) cat(sprintf("   [内层] g=%d f=%d a=%d b=%d gr=%d id=%d\n",
                          g, f, length(a), length(b), length(gr), length(H0$hhidpn)))
  data.table(id = H0$hhidpn, rnd = k, grip_wave = g, fall_wave = f, grip = gr,
             fall = fifelse(is.na(H0[[paste0("r", f, "fall")]]), NA_integer_,
                            as.integer(H0[[paste0("r", f, "fall")]] == 1)),
             birthyr = as.numeric(H0$BIRTHYR), female = as.integer(H0$GENDER == 2),
             edu = as.numeric(H0$raeducl), bmi = as.numeric(H0[[paste0("r", g, "mbmi")]]),
             smoke = as.integer(H0[[paste0("r", g, "smokef")]] >= 1),
             hibp = as.integer(gv(H0, paste0("r", g, "rxhibp")) == 1),
             diab = as.integer(gv(H0, paste0("r", g, "rxdiab")) == 1),
             adl = as.numeric(gv(H0, paste0("r", g, "adlfive"))),
             wt = as.numeric(H0[[paste0("r", g, "nwtresp")]]),
             psu = as.numeric(H0$SECU), strat = as.numeric(H0$STRATUM),
             EXDEATHYR = as.numeric(H0$EXDEATHYR), EXDEATHMO = as.numeric(H0$EXDEATHMO),
             结局年 = wave_year(f))
}))
HR[, age := (1990 + 2 * grip_wave) - birthyr]
HR[grip <= 0 | grip > 100, grip := NA_real_]
cat(sprintf("  [诊断] HR 行数 %d；female 非缺失 %d（女性 %d）；age 非缺失 %d；grip 非缺失 %d；wt 有效 %d\n",
            nrow(HR), sum(!is.na(HR$female)), sum(HR$female == 1, na.rm = TRUE),
            sum(!is.na(HR$age)), sum(!is.na(HR$grip)), sum(!is.na(HR$wt) & HR$wt > 0, na.rm = TRUE)))
HR <- HR[female == 1 & !is.na(age) & age >= 50 & !is.na(grip) & !is.na(wt) & wt > 0]
## 结局波年份之前死亡＝结局前死亡；死亡年 = 结局年时按月份判断（年中 6 月前记为死亡）
HR[, died := as.integer(!is.na(EXDEATHYR) &
                          (EXDEATHYR < 结局年 | (EXDEATHYR == 结局年 & !is.na(EXDEATHMO) & EXDEATHMO <= 6)))]
HR[, `:=`(obs = as.integer(!is.na(fall)),
          othermiss = as.integer(is.na(fall) & died == 0))]
cat(sprintf("HRS 符合条件窗口 %d；观察到结局 %d；结局前死亡 %d；其他缺失 %d\n",
            nrow(HR), sum(HR$obs), sum(HR$died), sum(HR$othermiss)))
print(HR[, .(符合条件 = .N, 观测 = sum(obs), 死亡 = sum(died), 其他缺失 = sum(othermiss)), by = .(窗口 = paste0(grip_wave, "→", fall_wave))])
saveRDS(HR[, .(id, rnd, grip_wave, fall_wave, grip, fall, age, edu, smoke, hibp, diab, bmi, adl,
               wt, psu, strat, died, obs, othermiss)], file.path(OUT, "elig_HRS.rds"))

cat("\n完成：", OUT, "\n")
