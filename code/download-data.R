library(tidyverse)
library(haven)
library(readxl)

# download a Stata file from a url, unzip it, and read it into memory
download_stata_data <- function(url) {
  # download
  tf <- tempfile()
  download.file(url, tf)

  # unzip
  fn <- unzip(tf, list = TRUE)$Name[1]
  con <- unz(tf, fn, open = "rb")

  # read in data
  data <- read_stata(con)

  # clean up
  close(con)
  unlink(tf)
  unlink(fn)

  # return data
  data
}

# download NAMCS and NHAMCS ED files for 2016 and 2018
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

# save NAMCS/NHAMCS data
write_rds(namcs_nhamcs, "data/namcs-nhamcs.rds", compress = "gz")

# download and read in the ICD-10 codes, then save that data
tf <- tempfile()
download.file("https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6334180/bin/chuk046645.ww2.xlsx", tf, mode = "wb")
icd_codes <- read_excel(tf, sheet = "2016 ICD-10-CM")
unlink(tf)
write_rds(icd_codes, "data/icd-codes.rds", compress = "gz")
