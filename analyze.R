library(tidyverse)
library(survey)

# allow for singleton PSUs
options(survey.lonely.psu = "adjust")

# read in raw ICD code data
raw_icd_codes <- read_rds("data/icd-codes.rds")

# read in raw NAMCS/NHAMCS data
raw_data <- read_rds("data/namcs-nhamcs.rds") %>%
  pull(data) %>%
  bind_rows()

# ICD codes in NAMCS/NHAMCS are coded in a special way
# this function converts from normal ICD codes to the NAMCS style
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

# convert the ICD codes
icd_codes <- raw_icd_codes %>%
  filter(str_length(ICD10_CODE) %in% c(3, 5)) %>%
  mutate(value = icd_to_diag(ICD10_CODE))

# get a list of antibiotics based on Multum categories
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
    racewt = PATWT / race_pop * 1e3,
    # a binary race variable will simplify comparisons later
    is_white = as.integer(RACERETH == "Non-Hispanic White")
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
    id, visit, racewt, is_white
  )

# for which visits were antibiotic prescribed?
abx_visit_ids <- raw_data2 %>%
  select(id, starts_with("RX")) %>%
  pivot_longer(!id) %>%
  filter(value %in% abx_codes$code) %>%
  pull(id) %>%
  unique()

# categorize antibiotic visits as appropriate, potentially appropriate, or
# inappropriate
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


# get rates of visits ----------------------------------------------------------

visit_rates <- svyby(~visit, by = ~RACERETH, design, svytotal, vartype = c("se", "ci")) %>%
  as_tibble() %>%
  mutate(target = "visit") %>%
  rename(estimate = visit)

# check for statistical significance: compare rates to whites
tibble(RACERETH = c("Non-Hispanic Black", "Hispanic", "Non-Hispanic Other")) %>%
  mutate(
    # compare this race against whites
    rows = map(RACERETH, ~ data$RACERETH %in% c("Non-Hispanic White", ..1)),
    this_design = map(rows, ~ subset(design, .)),
    # run the t-test
    test = map(this_design, ~ svyttest(is_white ~ 1, .)),
    result = map(test, broom::tidy)
  ) %>%
  select(RACERETH, result) %>%
  unnest(cols = result) %>%
  # check for Benjamini-Hochberg-adjusted significance
  mutate(sig = p.adjust(p.value, "BH") < 0.01)


# proportion of visits with abx ------------------------------------------------

tidy_svy_by_ciprop <- function(x, target) {
  as_tibble(x) %>%
    rename(estimate := !!target) %>%
    rename_with(function(x) str_replace(x, "^se\\..*", "se"))
}

abx_visits <- svyby(~got_abx, by = ~RACERETH, design, svyciprop, vartype = c("se", "ci")) %>%
  tidy_svy_by_ciprop("got_abx") %>%
  mutate(
    target = "got_abx",
    across(where(is.numeric), ~ . * 100)
  )

# proportion of visits with abx did not vary by race
svyciprop(~got_abx, design)
svychisq(~RACERETH + got_abx, design)

# appropriateness --------------------------------------------------------------
  
# proportion of abx visits that are appropriate, potentially appropriate, etc.
appropriateness <- tibble(target = c("appropriate", "potentially", "inappropriate")) %>%
  mutate(
    formula = map(target, ~ as.formula(paste0("~", .))),
    by_object = map(formula, function(x) svyby(x, by = ~RACERETH, subset(design, got_abx == 1), svyciprop, vartype = c("se", "ci"))),
    result = map2(by_object, target, tidy_svy_by_ciprop)
  ) %>%
  select(target, result) %>%
  unnest(cols = result) %>%
  mutate(across(where(is.numeric), ~ . * 100))

# appropriateness of an abx visit did not vary by race
svyciprop(~appropriate, subset(design, got_abx == 1))
svyciprop(~potentially, subset(design, got_abx == 1))
svyciprop(~inappropriate, subset(design, got_abx == 1))
svychisq(~RACERETH + category, subset(design, got_abx == 1))


# make a pretty table ----------------------------------------------------------

format_result <- function(estimate, ci_l, ci_u, se) {
  case_when(
    se / estimate > 0.3 ~ "*",
    TRUE ~ as.character(str_glue("{round(estimate)} ({round(ci_l)} to {round(ci_u)})"))
  )
}

bind_rows(visit_rates, abx_visits, appropriateness) %>%
  mutate(label = format_result(estimate, ci_l, ci_u, se)) %>%
  select(target, name = RACERETH, value = label) %>%
  pivot_wider()


# rate matching ----------------------------------------------------------------

bind_rows(
  visit = visit_rates,
  got_abx = abx_visits,
  .id = "target"
) %>%
  select(RACERETH, name = target, value = estimate) %>%
  pivot_wider() %>%
  mutate(
    is_white = RACERETH == "Non-Hispanic White",
    across(got_abx, ~ . / 100),
    observed_abx = visit * got_abx,
    match_abx = visit[is_white] * got_abx,
    disparity = observed_abx[is_white] - observed_abx,
    match_disparity = match_abx[is_white] - match_abx,
    explained = 1 - match_disparity / disparity
  ) %>%
  select(RACERETH, explained) %>%
  mutate(across(explained, ~ scales::percent(., accuracy = 1)))
