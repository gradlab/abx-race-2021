library(tidyverse)
library(survey)

options(survey.lonely.psu = "adjust")

raw_icd_codes <- read_rds("data/icd-codes.rds")

raw_data <- read_rds("data/namcs-nhamcs.rds") %>%
  pull(data) %>%
  bind_rows()

icd_to_diag <- function(x) {
  # all codes should be length 3 (eg A00) or length 5 (eg A00.9)
  stopifnot(all(str_length(x) %in% c(3, 5)))
  # which are length 5?
  long <- str_length(x) == 5
  
  # length-3 codes get a trailing dash
  x[!long] <- str_c(x[!long], "-")
  # all length-5 codes should have a dot in position 4
  stopifnot(all(str_sub(x[long], 4, 4) == "."))
  # length-5 codes get the period removed
  str_sub(x[long], 4, 4) <- ""
  
  # now all codes should be length 4
  stopifnot(all(str_length(x) == 4))
  x
}

# test: "I10" -> "I10-"; "Z00.1" -> "Z001"
# icd_to_diag(c("I10", "Z00.1"))
icd_codes <- raw_icd_codes %>%
  filter(str_length(ICD10_CODE) %in% c(3, 5)) %>%
  mutate(value = icd_to_diag(ICD10_CODE))

# from pg 275: https://ftp.cdc.gov/pub/Health_Statistics/NCHS/dataset_documentation/NHAMCS/doc18-ed-508.pdf
abx_codes <- tribble(
  ~code, ~abx_group,
  "013", "penicillins",
  "009", "cephalosporins",
  "011", "macrolide derivatives",
  "014", "quinolones",
  "240", "lincomycin derivatives",
  "015", "sulfonamides",
  "016", "tetracyclines",
  "017", "urinary antiinfectives",
  "486", "oxazolidinone antibiotics", # eg linezolid
  "002", "amebicides" # eg metronidazole
)

# from pg 18: https://ftp.cdc.gov/pub/Health_Statistics/NCHS/Dataset_Documentation/NAMCS/doc2016.pdf
race_populations <- tribble(
  ~YEAR, ~name, ~value,
  2016, "Hispanic", "56,809,443",
  2016, "Non-Hispanic all races", "261,205,874",
  2016, "Non-Hispanic White", "195,051,980",
  2016, "Non-Hispanic Black", "38,914,755",
  2018, "Hispanic", "59,214,173",
  2018, "Non-Hispanic all races", "262,890,676",
  2018, "Non-Hispanic White", "194,658,743",
  2018, "Non-Hispanic Black", "39,598,970"
) %>%
  mutate(across(value, parse_number)) %>%
  pivot_wider() %>%
  mutate(`Non-Hispanic Other` = `Non-Hispanic all races` - `Non-Hispanic White` - `Non-Hispanic Black`) %>%
  pivot_longer(!YEAR, names_to = "RACERETH", values_to = "race_pop") %>%
  mutate(across(RACERETH, ~ factor(., levels = names(attr(raw_data$RACERETH, "labels"))))) %>%
  filter(!is.na(RACERETH))

raw_data2 <- raw_data %>%
  mutate(across(c(RACERETH), as_factor)) %>%
  left_join(race_populations, by = c("YEAR", "RACERETH")) %>%
  mutate(
    # an ID will be helpful for merging in abx use and appropriateness
    id = 1:n(),
    # visit=1 is useful for counting #s of visits
    visit = 1,
    # cf slide 3: https://www.peppercenter.org/docs/StartWithData/Understanding_and_using_NAMCSandNHAMCS_data.pdf
    racewt = PATWT / race_pop * 1e3
  ) %>%
  select(
    # variables to keep:
    # race/ethnicity
    RACERETH,
    # prescriptions
    matches("RX\\d+CAT\\d+"),
    # diagnoses
    matches("^DIAG\\d+"),
    # survey details
    YEAR, CSTRATM, CPSUM, YEAR, PATWT,
    # derived values
    racewt, visit, id
  )

# for which visits were antibiotic prescribed?
abx_visit_ids <- raw_data2 %>%
  select(id, starts_with("RX")) %>%
  pivot_longer(!id) %>%
  filter(value %in% abx_codes$code) %>%
  pull(id) %>%
  unique()

abx_categories <- raw_data2 %>%
  filter(id %in% abx_visit_ids) %>%
  select(id, starts_with("DIAG")) %>%
  pivot_longer(!id) %>%
  # check here for missing codes
  left_join(icd_codes, by = "value") %>%
  group_by(id) %>%
  summarize(
    always = "A" %in% CATEGORY,
    sometimes = "S" %in% CATEGORY,
    never = "N" %in% CATEGORY
  ) %>%
  mutate(name = case_when(
    always ~ "appropriate",
    sometimes ~ "potentially",
    never ~ "inappropriate",
    TRUE ~ "not_associated"
  )) %>%
  select(id, name) %>%
  mutate(value = 1, category = name) %>%
  pivot_wider(values_fill = 0)

data <- raw_data2 %>%
  # remove the prescription and diagnosis columns
  select(!matches("^(RX|DIAG)")) %>%
  # note which visits are abx visits
  mutate(got_abx = as.integer(id %in% abx_visit_ids)) %>%
  # merge in appropriateness
  left_join(abx_categories, by = "id") %>%
  replace_na(list(appropriate = 0, potentially = 0, inappropriate = 0))

design <- svydesign(
  ids = ~CPSUM,
  strata = ~CSTRATM,
  weights = ~racewt,
  nest = TRUE,
  data = data
)


# healthcare and antibiotic use by race ----------------------------------------

tidy_svyby <- function(x, target) {
  as_tibble(x) %>%
    rename(estimate := !!target)
}

rates <- tibble(target = c("visit", "got_abx", "appropriate", "potentially", "inappropriate")) %>%
  mutate(
    formula = map(target, ~ as.formula(paste0("~", .))),
    by_object = map(formula, function(x) svyby(x, by = ~RACERETH, design, svytotal, vartype = c("se", "ci"))),
    result = map2(by_object, target, tidy_svyby)
  ) %>%
  select(target, result) %>%
  unnest(cols = result)

# make a pretty table
format_result <- function(estimate, ci_l, ci_u, se) {
  case_when(
    se / estimate > 0.3 ~ "*",
    TRUE ~ as.character(str_glue("{round(estimate)} ({round(ci_l)} to {round(ci_u)})"))
  )
}

rates %>%
  mutate(label = format_result(estimate, ci_l, ci_u, se)) %>%
  select(target, name = RACERETH, value = label) %>%
  pivot_wider()


# check for statistical significance: compare rates to whites
crossing(
  target = c("visit", "got_abx", "appropriate", "potentially", "inappropriate"),
  RACERETH = levels(data$RACERETH)
) %>%
  filter(RACERETH != "Non-Hispanic White") %>%
  mutate(
    # compare this race against whites
    rows = map2(RACERETH, target, ~ data$RACERETH %in% c("Non-Hispanic White", ..1) & data[[..2]] == 1),
    this_design = map(rows, ~ subset(design, .)),
    # create a formula for a t-test
    formula = map(RACERETH, function(x) as.formula(as.character('I(RACERETH == "{x}") ~ 1'))),
    # run the t-test
    test = map2(formula, this_design, svyttest),
    result = map(test, broom::tidy)
  ) %>%
  select(RACERETH, target, result) %>%
  unnest(cols = result) %>%
  # check for Benjamini-Hochberg-adjusted significance
  mutate(sig = p.adjust(p.value, "BH") < 0.01)

# proportion of visits with abx did not vary by race
svyciprop(~got_abx, design)
svychisq(~RACERETH + got_abx, design)

# appropriateness of an abx visit did not vary by race
svyciprop(~appropriate, subset(design, got_abx == 1))
svyciprop(~potentially, subset(design, got_abx == 1))
svyciprop(~inappropriate, subset(design, got_abx == 1))
svychisq(~RACERETH + category, subset(design, got_abx == 1))


# rate matching
rates %>%
  filter(target %in% c("visit", "got_abx")) %>%
  select(RACERETH, name = target, value = estimate) %>%
  pivot_wider() %>%
  mutate(
    abx_per_visit = got_abx / visit,
    match_visits = abx_per_visit * visit[RACERETH == "Non-Hispanic White"],
    match_apv = visit * abx_per_visit[RACERETH == "Non-Hispanic White"],
    across(c(got_abx, match_visits, match_apv), ~ .[RACERETH == "Non-Hispanic White"] - .),
    across(c(match_visits, match_apv), ~ 1 - . / got_abx)
  ) %>%
  select(RACERETH, match_visits, match_apv) %>%
  mutate(across(!RACERETH, ~ scales::percent(., accuracy = 1)))
