library(tidyverse)
library(survey)
library(cowplot)
library(scales)

# allow for singleton PSUs
options(survey.lonely.psu = "adjust")

# load data -----------------------------------------------------------

# read in raw ICD code data
raw_icd_codes <- read_rds("data/icd-codes.rds")
additional_codes <- read_csv("data/additional-codes.csv")

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
  mutate(value = icd_to_diag(ICD10_CODE)) %>%
  bind_rows(additional_codes)

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
  # compute non-hispanic other; drop non-hispanic all
  mutate(`Non-Hispanic Other` = `Non-Hispanic all races` - `Non-Hispanic White` - `Non-Hispanic Black`) %>%
  select(-`Non-Hispanic all races`) %>%
  pivot_longer(!YEAR, names_to = "RACERETH", values_to = "race_pop") %>%
  mutate(across(RACERETH, ~ factor(., levels = names(attr(raw_data$RACERETH, "labels")))))

raw_data2 <- raw_data %>%
  mutate(across(c(RACERETH), as_factor)) %>%
  left_join(race_populations, by = c("YEAR", "RACERETH")) %>%
  mutate(
    # an ID will be helpful for merging in abx use and appropriateness
    id = 1:n(),
    # visit=1 is useful for counting #s of visits
    visit = 1,
    # cf slide 3: https://www.peppercenter.org/docs/StartWithData/Understanding_and_using_NAMCSandNHAMCS_data.pdf
    racewt = PATWT / race_pop,
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
    id, visit, racewt
  )

# for which visits were antibiotic prescribed?
abx_visit_ids <- raw_data2 %>%
  select(id, starts_with("RX")) %>%
  pivot_longer(!id) %>%
  filter(value %in% abx_codes$code) %>%
  pull(id) %>%
  unique()

# categorize visits as appropriate, potentially appropriate, or
# inappropriate
abx_categories <- raw_data2 %>%
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
    never ~ "inappropriate"
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

visits_by_race <- svyby(~visit, by = ~RACERETH, design, svytotal, vartype = "ci") %>%
  as_tibble() %>%
  mutate(category = "total")

visits_by_race_and_cat <- svyby(~visit, by = ~RACERETH + category, design, svytotal, vartype = "ci") %>%
  as_tibble()

visits <- bind_rows(visits_by_race, visits_by_race_and_cat)


# test differences in visits ---------------------------------------------------

# y is -1 for whites and 1 for others. The t-test asks if the mean of that
# variable is zero, ie, whether, rates of visits are equal across two races
test_design <- update(design, y = if_else(RACERETH == "Non-Hispanic White", -1, 1))

test_total <- function(d) {
  st <- svytotal(~y, d)
  ci <- confint(st)
  
  list(
    estimate = coef(st),
    ci_l = ci[1],
    ci_u = ci[2],
    p.value = (2 * pt(-abs(coef(st) / SE(st)), df = degf(d) - 1))[, 1]
  ) %>%
    map(unname)
}

test_visits_by_race <- function(race) {
  test_total(subset(test_design, RACERETH %in% c(race, "Non-Hispanic White")))
}

test_visits_by_race_and_cat <- function(race, cat) {
  test_total(subset(test_design, RACERETH %in% c(race, "Non-Hispanic White") & category == cat))
}

test_races <- c("Non-Hispanic Black", "Hispanic", "Non-Hispanic Other")

tests_by_race <- tibble(race = test_races) %>%
  mutate(
    category = "total",
    result = map(race, test_visits_by_race)
  )

tests_by_race_and_cat <- crossing(
  race = test_races,
  category = c("appropriate", "potentially", "inappropriate")
) %>%
  mutate(result = map2(race, category, test_visits_by_race_and_cat))

test_results <- bind_rows(tests_by_race, tests_by_race_and_cat) %>%
  unnest_wider(result) %>%
  mutate(sig = p.value < 0.01)


# table of visit rate results --------------------------------------------------

category_levels <- c("total", "appropriate", "potentially", "inappropriate")
category_labels <- c("All visits", "Antibiotic-\nappropropriate\nvisits", "Sometimes\nantibiotic-appropriate\nvisits", "Antibiotic-\ninappropriate\nvisits")
race_levels <- c("Non-Hispanic White", "Non-Hispanic Black", "Hispanic", "Non-Hispanic Other")

visits %>%
  left_join(test_results, by = c("RACERETH" = "race", "category")) %>%
  mutate(
    # round everything to 2 sigfigs
    across(where(is.numeric), signif, digits = 2),
    # create a label for visit rates
    visit_label = str_glue("{visit} ({ci_l.x} to {ci_u.x})"),
    # create a label for differences in rates
    across(c(estimate, ci_l.y, ci_u.y), ~ -.),
    sig_label = if_else(sig, "*", "", missing = ""),
    diff_label = str_glue("{estimate} ({ci_l.y} to {ci_u.y}){sig_label}"),
    # get the category factors in order
    across(RACERETH, factor, levels = race_levels),
    across(category, factor, levels = category_levels)
  ) %>%
  select(category, RACERETH, visit = visit_label, diff = diff_label) %>%
  arrange(category, RACERETH)


# plot of visit rate results ---------------------------------------------------

my_palette  <- c("#FFFFFF", "#009E73", "#F0E442", "#E69F00", "#56B4E9", "#0072B2", "#D55E00", "#CC79A7")

visits_plot <- visits %>%
  left_join(select(test_results, RACERETH = race, category, sig), by = c("RACERETH", "category")) %>%
  mutate(
    across(category, ~ factor(., levels = category_levels, labels = category_labels)),
    across(RACERETH, ~ factor(., levels = race_levels)),
    sig_label = if_else(sig, "*", "", missing = ""),
    label = str_glue("{visit} ({ci_l} to {ci_u}){sig_label}")
  ) %>%
  ggplot(aes(category, visit, fill = RACERETH)) +
  geom_col(position = position_dodge(width = 0.9), color = "black") +
  geom_errorbar(aes(ymin = ci_l, ymax = ci_u), position = position_dodge(width = 0.9), width = 0.5) +
  geom_text(aes(label = sig_label, y = ci_u + 0.1), position = position_dodge(width = 0.9)) +
  labs(
    y = "Annual visits per capita",
    fill = "Race/ethnicity"
  ) +
  scale_y_continuous(limits = c(0, 8), expand = c(0, 0), breaks = 0:8) +
  scale_fill_manual(values = my_palette) +
  theme_cowplot() +
  theme(
    axis.title.x = element_blank(),
    legend.position = c(0.35, 0.75),
    axis.ticks.x = element_blank()
  )

show(visits_plot)


# proportion of visit categories with abx --------------------------------------

abx_by_race <- svyby(~got_abx, by = ~RACERETH, design, svyciprop, vartype = "ci") %>%
  as_tibble() %>%
  mutate(category = "total")

abx_by_race_and_cat <- svyby(~got_abx, by = ~RACERETH + category, design, svyciprop, vartype = "ci") %>%
  as_tibble()

bind_rows(abx_by_race, abx_by_race_and_cat) %>%
  mutate(
    across(category, ~ factor(., levels = category_levels)),
    across(where(is.numeric), ~ signif(., 2) * 100),
    label = str_glue("{got_abx}% ({ci_l}% to {ci_u}%)")
  ) %>%
  select(RACERETH, name = category, value = label) %>%
  pivot_wider() %>%
  select(all_of(c("RACERETH", category_levels)))

abx_plot <- bind_rows(abx_by_race, abx_by_race_and_cat) %>%
  mutate(
    across(category, ~ factor(., levels = category_levels, labels = category_labels)),
    across(RACERETH, ~ factor(., levels = race_levels)),
  ) %>%
  ggplot(aes(category, got_abx, fill = RACERETH)) +
  geom_col(position = position_dodge(width = 0.9), color = "black") +
  geom_errorbar(aes(ymin = ci_l, ymax = ci_u), position = position_dodge(width = 0.9), width = 0.5) +
  scale_y_continuous(
    labels = partial(scales::percent, accuracy = 1),
    limits = c(0, 0.55),
    expand = c(0, 0)
  ) +
  scale_fill_manual(values = my_palette) +
  labs(
    y = "% of visits with antibiotics",
    fill = "Race/ethnicity"
  ) +
  theme_cowplot() +
  theme(
    axis.title.x = element_blank(),
    legend.position = "none",
    axis.ticks.x = element_blank()
  )

show(abx_plot)

plot_grid(visits_plot, abx_plot, labels = "auto", ncol = 1)

ggsave("fig.png", width = 7, height = 7)


# test for proportions ---------------------------------------------------------

# baseline % of visits with abx
svyciprop(~got_abx, design)

# and by race
svyby(~got_abx, by = ~RACERETH, design, svyciprop, vartype = "ci") %>%
  as_tibble() %>%
  mutate(across(!RACERETH, ~ scales::percent(., accuracy = 0.1)))

# that baseline does not vary by race
svychisq(~RACERETH + got_abx, design)

# nor by race for any subset of visits
svychisq(~RACERETH + got_abx, subset(design, category == "appropriate"))
svychisq(~RACERETH + got_abx, subset(design, category == "potentially"))
svychisq(~RACERETH + got_abx, subset(design, category == "inappropriate"))


# rates of visits with antibiotics ------------------------------------

abx_visits_by_race <- svyby(~got_abx, by = ~RACERETH, design, svytotal, vartype = "ci") %>%
  as_tibble() %>%
  mutate(category = "total")

abx_visits_by_race_and_cat <- svyby(~got_abx, by = ~RACERETH + category, design, svytotal, vartype = "ci") %>%
  as_tibble()

test_abx_visits_by_race <- function(race) {
  test_total(subset(test_design, got_abx & RACERETH %in% c(race, "Non-Hispanic White")))
}

test_abx_visits_by_race_and_cat <- function(race, cat) {
  test_total(subset(test_design, got_abx & RACERETH %in% c(race, "Non-Hispanic White") & category == cat))
}

abx_tests_by_race <- tibble(race = test_races) %>%
  mutate(result = map(race, test_abx_visits_by_race)) %>%
  unnest_wider(result) %>%
  mutate(category = "total")

abx_tests_by_race_and_cat <- crossing(
  race = test_races,
  category = c("appropriate", "potentially", "inappropriate")
) %>%
  mutate(result = map2(race, category, test_abx_visits_by_race_and_cat)) %>%
  unnest_wider(result)

# put together a table of results
left_join(
  bind_rows(abx_visits_by_race, abx_visits_by_race_and_cat),
  bind_rows(abx_tests_by_race, abx_tests_by_race_and_cat),
  by = c("RACERETH" = "race", "category")
) %>%
  mutate(
    across(c(estimate, ci_l.y, ci_u.y), ~ -.),
    across(where(is.numeric), signif, digits = 2),
    sig_label = if_else(p.value < 0.01, "*", ""),
    rate_label = str_glue("{got_abx} ({ci_l.x} to {ci_u.x})"),
    diff_label = str_glue("{estimate} ({ci_u.y} to {ci_l.y}){sig_label}"),
    across(category, factor, levels = category_levels),
    across(RACERETH, factor, levels = race_levels)
  ) %>%
  select(category, RACERETH, rate = rate_label, diff = diff_label) %>%
  arrange(category, RACERETH)
