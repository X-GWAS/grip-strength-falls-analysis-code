## db_paths.R —— 七大公共数据库路径解析（跨对话复用）
## 2026-09-16 起：数据全部存放在 TSD302 移动盘，本地不再保存。
## 用法（R 脚本开头）：
##   source("/Users/mac/.codex/knowledge/数据库资产/db_paths.R")
##   H <- read_dta(file.path(DB$charls, "Harmonized_CHARLS_D/H_CHARLS_D_Data.dta"))

.find_root <- function(cands, label) {
  for (p in cands) if (dir.exists(p)) return(p)
  stop(sprintf("【%s】数据目录不存在。请确认 TSD302 移动盘已挂载（/Volumes/TSD302）。", label),
       call. = FALSE)
}

.DBROOT_TSD302 <- "/Volumes/TSD302/七大数据库"

DB <- list(
  charls = .find_root(c(file.path(.DBROOT_TSD302, "charls数据库"),
                        "/Users/mac/Downloads/charls数据库"), "CHARLS"),
  elsa   = .find_root(c(file.path(.DBROOT_TSD302, "2. ELSA 英国"),   "/Users/mac/Downloads/ELSA"),   "ELSA"),
  hrs    = .find_root(c(file.path(.DBROOT_TSD302, "3. HRS  美国"),  "/Users/mac/Downloads/HRS"),    "HRS"),
  klosa  = .find_root(c(file.path(.DBROOT_TSD302, "4. KLoSA 韩国"),  "/Users/mac/Downloads/KLoSA"),  "KLoSA"),
  lasi   = .find_root(c(file.path(.DBROOT_TSD302, "5. LASI 印度"),   "/Users/mac/Downloads/LASI"),   "LASI"),
  mhas   = .find_root(c(file.path(.DBROOT_TSD302, "6. MHAS 墨西哥"), "/Users/mac/Downloads/MHAS"),   "MHAS"),
  share  = .find_root(c(file.path(.DBROOT_TSD302, "7.SHARE 欧洲"),   "/Users/mac/Downloads/SHARE"),  "SHARE")
)

## harmonized 主数据文件（实测路径，2026-09-16 核对）
DB$file <- list(
  charls_h = file.path(DB$charls, "Harmonized_CHARLS_D/H_CHARLS_D_Data.dta"),
  elsa_h   = file.path(DB$elsa,   "Raw_data/Harmonized ELSA/h_elsa_g3.dta"),
  hrs_h    = file.path(DB$hrs,    "HRS_美国/Raw_data/Gateway Harmonized HRS/H_HRS_d.dta"),
  lasi_h   = file.path(DB$lasi,   "LASI_印度/Raw_data/Harmonized LASI (A.3)/H_LASI_a3.dta"),
  mhas_h   = file.path(DB$mhas,   "MHAS_墨西哥/Raw_data/Harmonized MHAS File/H_MHAS_c2.dta"),
  share_h  = file.path(DB$share,  "SHARE_欧洲/Raw_data/Harmonized SHARE/H_SHARE_f2.dta")
)
## 注意：KLoSA 无 harmonized 文件，只有原始波次（Raw_data/wave*/w0X.dta）与 Temp_data，
##       需自行 harmonize；SHARE 的跌倒变量带 "_s" 后缀（r1fall_s 等）。
