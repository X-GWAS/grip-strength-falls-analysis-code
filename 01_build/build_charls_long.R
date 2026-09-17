suppressMessages({library(haven); library(data.table)})

## 数据路径：2026-09-16 起全部迁至 TSD302 移动盘，统一由 db_paths.R 解析
suppressWarnings(source("/Users/mac/.codex/knowledge/数据库资产/db_paths.R"))
root <- paste0(DB$charls, "/")

# ---------- 1) harmonized CHARLS D: longitudinal backbone ----------
h <- as.data.table(read_dta(paste0(root, "Harmonized_CHARLS_D/H_CHARLS_D_Data.dta")))
h <- h[, .(ID, communityID, ID_w1, ragender,
           r1agey, r2agey, r3agey, r4agey,
           inw1, inw2, inw3, inw4,
           r1iwy, r1iwm, r2iwy, r2iwm, r3iwy, r3iwm, r4iwy, r4iwm,
           r1wtresp, r1diabe, r2diabe, r3diabe, r4diabe,
           r1rxdiab, r2rxdiab, r3rxdiab, r4rxdiab,
           raeduc_c, r1mstat, h1rural, r1mbmi, r1mwaist,
           r1smokev, r1smoken, r1drinkev, r1drinkl,
           r1hibpe, r1dyslipe, r1vgact_c, r1mdact_c)]

# ---------- 2) raw wave 1 female reproductive exposure ----------
hs1 <- as.data.table(read_dta(paste0(root, "2011charls/health_status_and_functioning.dta")))
hs1 <- hs1[, .(ID = as.character(ID),
               menarche_raw = da026_2,
               meno_start = da027,
               meno_age_raw = da028_2,
               dm_self_w1 = da007_3_)]

# ---------- 3) raw wave 1 blood (2011) ----------
bld1 <- as.data.table(read_dta(paste0(root, "2011charls/Blood_20140429.dta")))
bld1 <- bld1[, .(ID = as.character(ID), glu11 = newglu, hba11 = newhba1c)]

# ---------- 4) raw wave 3 blood (2015) for lab-augmented DM at w3 ----------
bld3 <- as.data.table(read_dta(paste0(root, "2015charls/Blood.dta")))
bld3 <- bld3[, .(ID = as.character(ID), fast15 = bl_fasting, glu15 = bl_glu, hba15 = bl_hbalc)]

# ---------- 5) PSU/strata ----------
psu <- as.data.table(read_dta(paste0(root, "2011charls/psu.dta")))
psu[, communityID := as.character(communityID)]
psu[, strata_id := paste(province, urban_nbs, sep = "_")]

# ---------- merge ----------
d <- merge(h, hs1, by.x = "ID_w1", by.y = "ID", all.x = TRUE)
d <- merge(d, bld1, by.x = "ID_w1", by.y = "ID", all.x = TRUE)
d <- merge(d, bld3, by = "ID", all.x = TRUE)
d <- merge(d, psu[, .(communityID, strata_id)], by = "communityID", all.x = TRUE)
d[, strata_prov := sub("_[01]$", "", strata_id)]   # province-only strata (robust to sparse cells)

# ---------- cleaning ----------
d[, sex_f := ragender == 2]
d[is.na(ragender), sex_f := NA]
d[, age1 := r1agey]

# exposures (clean implausible year-like entries)
d[, menarche := fifelse(!is.na(menarche_raw) & menarche_raw %in% 8:25, menarche_raw, NA_real_)]
d[, meno_ever := fifelse(meno_start == 1, 1L, 0L)]
d[is.na(meno_start), meno_ever := NA_integer_]
d[, meno_age := fifelse(meno_ever == 1L & !is.na(meno_age_raw) & meno_age_raw %in% 20:60, meno_age_raw, NA_real_)]
d[, early_menarche12 := fifelse(!is.na(menarche) & menarche < 12, 1L, 0L)]
d[, early_menarche13 := fifelse(!is.na(menarche) & menarche < 13, 1L, 0L)]
d[, early_meno45 := fifelse(!is.na(meno_age) & meno_age < 45, 1L, 0L)]
d[, early_meno40 := fifelse(!is.na(meno_age) & meno_age < 40, 1L, 0L)]

# diabetes at waves (self-report + meds)
for (k in 1:4) {
  d[, paste0("dm", k) := fifelse(get(paste0("r", k, "diabe")) == 1, 1L, 0L)]
  d[get(paste0("r", k, "diabe")) == 0, paste0("dm", k) := 0L]
  d[, paste0("rx", k) := fifelse(get(paste0("r", k, "rxdiab")) == 1, 1L, 0L)]
  d[get(paste0("r", k, "rxdiab")) == 0, paste0("rx", k) := 0L]
}

# lab-augmented any DM at w1 (2011) and w3 (2015)
d[, lab_dm11 := fifelse((!is.na(hba11) & hba11 >= 6.5) | (!is.na(glu11) & glu11 >= 126), 1L, 0L)]
d[is.na(hba11) & is.na(glu11), lab_dm11 := NA_integer_]
d[, anydm11 := fifelse(dm1 == 1L | lab_dm11 == 1L, 1L, 0L)]
d[is.na(dm1) & is.na(lab_dm11), anydm11 := NA_integer_]

d[, lab_dm15 := fifelse((!is.na(hba15) & hba15 >= 6.5) |
                          (!is.na(fast15) & fast15 == 1 & !is.na(glu15) & glu15 >= 126), 1L, 0L)]
d[is.na(hba15) & (is.na(fast15) | fast15 != 1 | is.na(glu15)), lab_dm15 := NA_integer_]
d[, anydm15 := fifelse(dm3 == 1L | lab_dm15 == 1L, 1L, 0L)]
d[is.na(dm3) & is.na(lab_dm15), anydm15 := NA_integer_]

# covariates recode
d[, educ4 := fcase(raeduc_c %in% 1:2, "1_文盲及未完成小学",
                   raeduc_c %in% 3:4, "2_小学",
                   raeduc_c == 5, "3_初中",
                   raeduc_c %in% 6:10, "4_高中及以上",
                   default = NA_character_)]
d[, married := fifelse(r1mstat %in% c(1, 3), 1L, 0L)]
d[is.na(r1mstat), married := NA_integer_]
d[, rural := fifelse(h1rural == 1, 1L, 0L)]
d[is.na(h1rural), rural := NA_integer_]
d[, bmi1 := r1mbmi]
d[, smoke_cat := fcase(r1smoken == 1, "now",
                       r1smokev == 1 & r1smoken == 0, "former",
                       r1smokev == 0, "never",
                       default = NA_character_)]
d[, drink1y := fifelse(r1drinkl == 1, 1L, 0L)]
d[is.na(r1drinkl), drink1y := NA_integer_]
d[, htn := fifelse(r1hibpe == 1, 1L, 0L)]
d[is.na(r1hibpe), htn := NA_integer_]
d[, dyslip := fifelse(r1dyslipe == 1, 1L, 0L)]
d[is.na(r1dyslipe), dyslip := NA_integer_]
d[, act_mod := fifelse(r1mdact_c == 1, 1L, 0L)]
d[is.na(r1mdact_c), act_mod := NA_integer_]

# date (decimal year) at each interview
d[, t1 := r1iwy + (r1iwm - 0.5) / 12]
d[, t2 := r2iwy + (r2iwm - 0.5) / 12]
d[, t3 := r3iwy + (r3iwm - 0.5) / 12]
d[, t4 := r4iwy + (r4iwm - 0.5) / 12]

fwrite(d, "/Users/mac/Documents/Codex/2026-08-31/xian/work/charls_long_build.csv")
cat("saved rows:", nrow(d), "\n")

# ---------- exploration ----------
base <- d[sex_f == TRUE & age1 >= 45]
cat("\nWomen >=45 w1 interviewed:", base[inw1 == 1, .N], "\n")
cat("  with menarche:", base[!is.na(menarche), .N], "\n")
cat("  with menopause info (ever/never):", base[!is.na(meno_ever), .N],
    "  post:", base[meno_ever == 1, .N], "\n")
cat("  with meno_age valid:", base[!is.na(meno_age), .N], "\n")
cat("  dm1 self:", base[dm1 == 1, .N], " anydm11 (lab+self):", base[anydm11 == 1, .N], "\n")
cat("  w2/w3/w4 interviewed & dm known:",
    base[!is.na(dm2), .N], "/", base[!is.na(dm3), .N], "/", base[!is.na(dm4), .N], "\n")
cat("  incident self DM w2/3/4 among dm1==0:",
    base[dm1 == 0 & dm2 == 1, .N], "/", base[dm1 == 0 & dm3 == 1, .N], "/",
    base[dm1 == 0 & dm4 == 1, .N], "\n")
cat("  lab DM w3 among dm3==0: ", base[dm3 == 0 & lab_dm15 == 1, .N], "\n", sep = "")

cat("\nmenarche dist (women>=45):\n")
print(quantile(base$menarche, probs = c(0.01, 0.1, 0.5, 0.9, 0.99), na.rm = TRUE))
cat("menopause age dist:\n")
print(quantile(base$meno_age, probs = c(0.01, 0.1, 0.5, 0.9, 0.99), na.rm = TRUE))
cat("\ncovariate completeness (non-missing / cohort):\n")
for (vv in c("bmi1", "educ4", "married", "rural", "smoke_cat", "drink1y", "htn", "dyslip", "act_mod")) {
  cat(" ", vv, ":", sum(!is.na(base[[vv]])), "\n")
}
