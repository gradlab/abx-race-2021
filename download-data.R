library(tidyverse)
library(haven)
library(readxl)

download_stata_data <- function(url) {
  tf <- tempfile()
  download.file(url, tf)
  fn <- unzip(tf, list = TRUE)$Name[1]
  con <- unz(tf, fn, open = "rb")
  data <- read_stata(con)
  close(con)
  data
}

namcs_nhamcs <- tribble(
  ~dir, ~base_filename,
  "namcs", "namcs",
  "nhamcs", "ed"
) %>%
  crossing(year = c(2016, 2018)) %>%
  mutate(
    url = str_glue("https://ftp.cdc.gov/pub/Health_Statistics/NCHS/dataset_documentation/{dir}/stata/{base_filename}{year}-stata.zip"),
    data = map(url, download_stata_data)
  )

write_rds(namcs_nhamcs, "data/namcs-nhamcs.rds", compress = "gz")

tf <- tempfile()
download.file("https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6334180/bin/chuk046645.ww2.xlsx", tf, mode = "wb")
icd_codes <- read_excel(tf, sheet = "2016 ICD-10-CM")
unlink(tf)
write_rds(icd_codes, "data/icd-codes.rds", compress = "gz")
