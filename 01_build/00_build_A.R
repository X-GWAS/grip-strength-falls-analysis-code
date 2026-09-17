## 00_build_A.R —— 路径A 建库：累积体力活动 / 握力轨迹 → 跌倒·骨折·体能
## 设计：女性 ≥45 岁；暴露 2011–2015（三波）；结局 2015→2018（前瞻，暴露在前）
suppressMessages({ library(haven); library(data.table) })

## 数据路径：2026-09-16 起全部迁至 TSD302 移动盘，统一由 db_paths.R 解析
suppressWarnings(source("/Users/mac/.codex/knowledge/数据库资产/db_paths.R"))
root <- DB$charls
out  <- "/Users/mac/Documents/Codex/2026-08-31/xian/work/新暴露-体力活动与握力轨迹"
dir.create(out, showWarnings = FALSE, recursive = TRUE)

gv <- function(dt, v) if (v %in% names(dt)) dt[[v]] else rep(NA_real_, nrow(dt))

## ---------- 1) Harmonized CHARLS D ----------
cat("读取 Harmonized CHARLS D ...\n")
H <- as.data.table(read_dta(file.path(root, "Harmonized_CHARLS_D/H_CHARLS_D_Data.dta")))

keep <- c("ID","communityID","hhid","ragender","raeduc_c","raeducl",
          "r1agey","r2agey","r3agey","r4agey",
          "r1mstat","h1rural","r1rural2",
          "r1smokev","r1drinkev","r1drinkl","r1drinkn_c",
          "r1hibpe","r1diabe","r1cancre","r1lunge","r1hearte","r1stroke",
          "r1arthre","r1memrye","r1adla_c","r1adlwa","r1iadla_c","r1cesd10",
          "r1mheight","r3mheight","r1mweight","r3mweight",
          "r1vgact_c","r2vgact_c","r3vgact_c","r4vgact_c",
          "r1mdact_c","r2mdact_c","r3mdact_c","r4mdact_c",
          "r1ltact_c","r2ltact_c","r3ltact_c","r4ltact_c",
          "r1vgactx_c","r2vgactx_c","r3vgactx_c","r4vgactx_c",
          "r1mdactx_c","r2mdactx_c","r3mdactx_c","r4mdactx_c",
          "r1ltactx_c","r2ltactx_c","r3ltactx_c","r4ltactx_c",
          "r1gripsum","r2gripsum","r3gripsum",
          "r1wspeed","r2wspeed","r3wspeed","r1chaira","r3chaira",
          "r1wtresp","r3wtresp","r1wtrespb","r3wtrespb",
          "r1iwstat","r2iwstat","r3iwstat","r4iwstat")
keep <- intersect(keep, names(H))
H <- H[, ..keep]
cat("Harmonized 行数:", nrow(H), " 列数:", ncol(H), "\n")

## ---------- 2) 跌倒 / 髋骨骨折（原始问卷，Harmonized 无此项） ----------
read_ff <- function(path, fvars, tag) {
  d <- as.data.table(read_dta(path))
  d <- d[, c("ID", intersect(fvars, names(d))), with = FALSE]
  setnames(d, setdiff(names(d), "ID"), paste0(tag, "_", setdiff(names(d), "ID")))
  d
}
ff13 <- read_ff(file.path(root,"2013charls/Health_Status_and_Functioning.dta"),
                c("da023","da024","da025"), "w2")
ff15 <- read_ff(file.path(root,"2015charls/Health_Status_and_Functioning.dta"),
                c("da023","da024","da025"), "w3")
ff18 <- read_ff(file.path(root,"2018charls/Health_Status_and_Functioning.dta"),
                c("da023_w4","da024","da025_w4"), "w4")

## ---------- 3) 抽样设计信息（省 / 城乡） ----------
psu <- as.data.table(read_dta(file.path(root,"2011charls/psu.dta")))
psu <- unique(psu[, .(communityID, province, urban_nbs)])
psu[, `:=`(stratum = paste(province, urban_nbs, sep = "_"))]

## ---------- 4) 合并 ----------
D <- merge(H, psu, by = "communityID", all.x = TRUE)
D <- merge(D, ff13, by = "ID", all.x = TRUE)
D <- merge(D, ff15, by = "ID", all.x = TRUE)
D <- merge(D, ff18, by = "ID", all.x = TRUE)

## ---------- 5) 暴露：各波 MVPA（中高强度体力活动，天/周） ----------
mv <- function(v, m) pmin(7, fifelse(is.na(v), 0, v) + fifelse(is.na(m), 0, m))
D[, `:=`(
  mvpa_w1 = mv(gv(D,"r1vgactx_c"), gv(D,"r1mdactx_c")),
  mvpa_w2 = mv(gv(D,"r2vgactx_c"), gv(D,"r2mdactx_c")),
  mvpa_w3 = mv(gv(D,"r3vgactx_c"), gv(D,"r3mdactx_c")),
  mvpa_w4 = mv(gv(D,"r4vgactx_c"), gv(D,"r4mdactx_c")),
  ltpa_w1 = gv(D,"r1ltactx_c"), ltpa_w2 = gv(D,"r2ltactx_c"),
  ltpa_w3 = gv(D,"r3ltactx_c"), ltpa_w4 = gv(D,"r4ltactx_c")
)]
## 三波中高强度活动是否"有测量"
for (w in 1:3) {
  D[[paste0("has_mvpa", w)]] <- !is.na(D[[paste0("r", w, "vgactx_c")]]) |
                               !is.na(D[[paste0("r", w, "mdactx_c")]])
}
D[, n_mvpa_wave := has_mvpa1 + has_mvpa2 + has_mvpa3]
D[, cum_mvpa := rowMeans(.SD, na.rm = TRUE), .SDcols = c("mvpa_w1","mvpa_w2","mvpa_w3")]
D[n_mvpa_wave == 0, cum_mvpa := NA_real_]

## ---------- 6) 结局变量 ----------
D[, `:=`(
  fall_1113   = fifelse(!is.na(w2_da023), as.integer(w2_da023 == 1), NA_integer_),
  fall_1315   = fifelse(!is.na(w3_da023), as.integer(w3_da023 == 1), NA_integer_),
  fall_1518   = fifelse(!is.na(w4_da023_w4), as.integer(w4_da023_w4 == 1), NA_integer_),
  fx_1113     = fifelse(!is.na(w2_da025), as.integer(w2_da025 == 1), NA_integer_),
  fx_1315     = fifelse(!is.na(w3_da025), as.integer(w3_da025 == 1), NA_integer_),
  fx_1518     = fifelse(!is.na(w4_da025_w4), as.integer(w4_da025_w4 == 1), NA_integer_)
)]
## 跌倒受伤需治疗（DA024 次数 ≥1，仅对跌倒者问）
D[, fall_inj_1518 := fifelse(is.na(w4_da024), NA_integer_, as.integer(w4_da024 >= 1))]
D[, fall_inj_1315 := fifelse(is.na(w3_da024), NA_integer_, as.integer(w3_da024 >= 1))]

## ---------- 7) 体能结局（2015） ----------
D[, `:=`(
  grip_2015  = as.numeric(gv(D,"r3gripsum")),
  grip_2011  = as.numeric(gv(D,"r1gripsum")),
  speed_2015 = fifelse(!is.na(gv(D,"r3wspeed")) & gv(D,"r3wspeed") >= 1, 2.5/as.numeric(gv(D,"r3wspeed")), NA_real_),
  speed_2011 = fifelse(!is.na(gv(D,"r1wspeed")) & gv(D,"r1wspeed") >= 1, 2.5/as.numeric(gv(D,"r1wspeed")), NA_real_)
)]
## 握力异常值清理
D[grip_2015 <= 0 | grip_2015 > 100, grip_2015 := NA_real_]
D[grip_2011 <= 0 | grip_2011 > 100, grip_2011 := NA_real_]

## ---------- 8) 协变量 ----------
D[, `:=`(
  female   = as.integer(ragender == 2),
  age      = as.numeric(gv(D,"r1agey")),
  age3     = as.numeric(gv(D,"r3agey")),
  edu      = as.numeric(gv(D,"raeduc_c")),
  married  = fifelse(!is.na(gv(D,"r1mstat")), as.integer(gv(D,"r1mstat") %in% c(1,2)), NA_integer_),
  rural    = as.numeric(gv(D,"h1rural")),
  smoke    = fifelse(!is.na(gv(D,"r1smokev")), as.integer(gv(D,"r1smokev") == 1), NA_integer_),
  drink    = fifelse(!is.na(gv(D,"r1drinkev")), as.integer(gv(D,"r1drinkev") == 1), NA_integer_),
  hibp     = fifelse(!is.na(gv(D,"r1hibpe")), as.integer(gv(D,"r1hibpe") == 1), NA_integer_),
  diab     = fifelse(!is.na(gv(D,"r1diabe")), as.integer(gv(D,"r1diabe") == 1), NA_integer_),
  adl      = as.numeric(gv(D,"r1adla_c")),
  cesd     = as.numeric(gv(D,"r1cesd10")),
  height1  = as.numeric(gv(D,"r1mheight")),
  weight1  = as.numeric(gv(D,"r1mweight")),
  wt1      = as.numeric(gv(D,"r1wtresp")),
  wt3      = as.numeric(gv(D,"r3wtresp")),
  wt1b     = as.numeric(gv(D,"r1wtrespb")),
  wt3b     = as.numeric(gv(D,"r3wtrespb"))
)]
## 身高若以 m 记录则换算为 cm
if (mean(D$height1, na.rm = TRUE) < 3) D[, height1 := height1 * 100]
D[, bmi := fifelse(height1 > 100 & weight1 > 20, weight1/(height1/100)^2, NA_real_)]
D[bmi < 10 | bmi > 60, bmi := NA_real_]

## ---------- 9) 分析队列：女性 ≥45 岁 ----------
setorder(D, ID)
saveRDS(D, file.path(out, "build_D.rds"))

cat("\n===== 队列盘点（女性 ≥45 岁）=====\n")
W <- D[ragender == 2 & !is.na(age) & age >= 45]
cat("2011 年女性 ≥45 岁:", nrow(W), "\n")
cat("其中三波 MVPA 中位可得波数:", median(W$n_mvpa_wave, na.rm = TRUE), "\n")
cat("有 ≥1 波 MVPA:", sum(W$n_mvpa_wave >= 1, na.rm = TRUE), "\n")
cat("有 ≥2 波 MVPA:", sum(W$n_mvpa_wave >= 2, na.rm = TRUE), "\n")
cat("累积 MVPA 非缺失:", sum(!is.na(W$cum_mvpa)), " 均值=", round(mean(W$cum_mvpa, na.rm = TRUE), 2),
    " SD=", round(sd(W$cum_mvpa, na.rm = TRUE), 2), "\n")
cat("\n结局事件:\n")
for (v in c("fall_1113","fall_1315","fall_1518","fall_inj_1518","fx_1113","fx_1315","fx_1518")) {
  cat(sprintf("  %-14s 非缺失=%6d  事件=%5d\n", v, sum(!is.na(W[[v]])), sum(W[[v]] == 1, na.rm = TRUE)))
}
cat("\n体能结局:\n")
for (v in c("grip_2011","grip_2015","speed_2011","speed_2015")) {
  cat(sprintf("  %-12s 非缺失=%6d  均值=%.3f\n", v, sum(!is.na(W[[v]])), mean(W[[v]], na.rm = TRUE)))
}
cat("\n权重: wt1 非缺失", sum(!is.na(W$wt1)), " wt1b 非缺失", sum(!is.na(W$wt1b)),
    " wt3 非缺失", sum(!is.na(W$wt3)), "\n")
cat("社区数:", length(unique(W$communityID)), " 层数:", length(unique(W$stratum)), "\n")
