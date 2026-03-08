# packages, seed, contrasts ####

library(bayestestR)
library(brms)
library(irr)
library(janitor)
library(patchwork)
library(readxl)
library(rstan)
library(StanHeaders)
library(tidyverse)

set.seed(33322)

options(contrasts = c("contr.sum", "contr.poly"))

# custom functions ####

brm_func <- function(categorical, data, prior, formula, control){
  if(missing(categorical)){
    brm(formula = formula,
      data    = data,
      family  = bernoulli,
      prior   = prior,
      iter    = 10000,
      warmup  = 2500,
      chains  = 7,
      cores   = 7,
      save_pars = save_pars(all = TRUE),
      backend = "cmdstanr",
      refresh = 5000)
  }
  else{
    brm(formula = formula,
        data    = data,
        family  = categorical,
        prior   = prior,
        iter    = 10000,
        warmup  = 2500,
        chains  = 7,
        cores   = 7,
        save_pars = save_pars(all = TRUE),
        backend = "cmdstanr",
        refresh = 5000)
  }
}

plot_func <- function(byparts, df, condition, facet_x, facet_y){
  if(missing(byparts)){
  ggplot(data = df %>% filter(!is.na(prompt)), 
         aes(x = {{condition}}, fill = prompt)) +
    facet_grid(vars({{facet_y}}),
               vars({{facet_x}})) +
    geom_bar(position = "fill") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    geom_text(stat = "count", aes(y = after_stat(count), label = after_stat(count), vjust = ifelse(after_stat(count) < 4, -0.25, 1.25)), position = "fill", size = 8) +
    guides(fill = guide_legend(title = "Abandoned on")) +
    scale_fill_manual(values = c("lightblue", "limegreen"), labels = c("1st prompt", "2nd+ prompt")) +
    theme(text = element_text(size = 20, family = "Times New Roman")) +
    labs(y = "Percentage") 
  }
  else{
    ggplot(data = df, 
           aes(x = {{condition}}, fill = ratio)) +
      facet_grid(vars({{facet_y}}),
                 vars({{facet_x}})) +
      geom_bar(position = "fill") +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
      geom_text(stat = "count", aes(y = after_stat(count), label = after_stat(count), vjust = ifelse(after_stat(count) < 2, -0.25, 1.25)), position = "fill", size = 8) +
      guides(fill = guide_legend(title = "Abandoned on 2+")) +
      scale_fill_manual(values = c("lightblue", "limegreen", "forest green")) +
      theme(text = element_text(size = 20, family = "Times New Roman")) +
      labs(y = "Percentage")
  }
}

table_func <- function(mod){
  hypoth <- hypothesis(mod, paste(rownames(fixef(mod)), "> 0"))[[1]]
  fixef(mod) %>%
    as.data.frame() %>%
    rownames_to_column("Parameter") %>%
    rename(Error = "Est.Error") %>% 
    mutate(Estimate = as.numeric(map_estimate(mod)[2]), # MAP estimate
           prob = hypoth$Post.Prob,
           "Ev. Ratio" = hypoth$Evid.Ratio,
           across(.cols = 2:5, ~round(.x, 2)),
           across(.cols = 2:5, ~format(., nsmall = 2)), # keep trailing zeros
           across(.cols = 6:7, ~round(.x, 3))) %>% 
    unite(HDI, Q2.5:Q97.5, sep = ", ") %>%
    mutate("95% HDI" = paste0("[", HDI, "]")) %>%
    select(-HDI) %>%
    relocate(prob, .after = everything()) %>% 
    relocate("Ev. Ratio", .after = everything()) %>% 
    mutate("Evidence" = ifelse(prob > 0.100 & prob < 0.899, "weak",
                                   ifelse(prob >= 0.950 | prob <= 0.050, "strong", "moderate")),
           Parameter = str_replace_all(Parameter, pattern = "_", replacement = " "),
           across(.cols = 6:7, ~format(., nsmall = 3))) %>% # keep trailing zeros 
    rename("Post. Prob" = prob) 
}

# data ####

## import ####

weco_reli_prompt <- read_xlsx("data/weco_reli.xlsx") %>%
  clean_names() %>% 
  mutate_at(vars(everything()), as.character)

weco_demo <- read_xlsx("data/weco_demo.xlsx") %>%
  clean_names() %>%
  mutate(birth_date = as.numeric(birth_date),
         birth_date = as.Date(birth_date, origin = "1900-01-01"), # fix a read_xl import issue, birth_date incorrectly converted to charater var
         test_date  = as.Date(test_date, origin = "1900-01-01"),
         age_days   = test_date - birth_date,
         age_years  = as.double(age_days/365.25)) %>% 
  select(part_id, sex, race, income, age_years)

weco_data_dirty <- read_xlsx("data/weco_data.xlsx") %>% 
  clean_names() %>%
  rename("ref_exp" = "we") %>% 
  group_by(part_id) %>% 
  mutate(trial_num = row_number()) %>% 
  ungroup()

## wrangling ####

weco_data <- left_join(weco_data_dirty, weco_demo, by = "part_id") %>% 
  mutate(age = ifelse(age_years > 3.00, "3_yo", "2_yo")) # automate age group creation

weco_data_nopilot <- weco_data %>%
  filter(!str_detect(part_id, 'weco_pilot')) %>% # remove rows of pilot data
  select(-c(date, location, notes)) #%>%
   # distinct(part_id) # 111 participants

weco_data_nopilot %>%
  group_by(final_sample) %>% 
  summarize(n = n())

weco_data_nopilot %>%
  filter(final_sample == "0") %>%  # retain only final sample excluded trials
  group_by(part_id, age) %>% 
  summarise(trials_excluded = n()) %>%
  ungroup() %>% 
  group_by(age, trials_excluded) %>% 
  summarize(n = n()) # N = 5 participants 2 trials excluded (N = 101 participants total; there are 111 participants in dataset, these extra 10 were sampled just for security, if needed)

weco_data_final_bytrial <- weco_data_nopilot %>%
  mutate(sex = fct_recode(sex,
                          "Male" = "m",
                          "Female" = "f"),
         age = ifelse(age == "2_yo", "2-year-olds", 
                      ifelse(age == "3_yo", "3-year-olds", NA)), # weco_final74 never completed demographics, thus age neither "2_yo" nor "3_yo" but "NA", so removing this participant
         order_partner = fct_recode(order_partner,
                                    "Dyad first" = "individual_first",
                                    "Group first" = "group_first"),
         condition_partner = fct_recode(condition_partner,
                                        "Dyad" = "individual",
                                        "Group" = "group"),
         condition_framing = fct_recode(condition_framing,
                                        "We" = "we",
                                        "You" = "you"),
         race = fct_recode(race,
                           "Asian" = "asian",
                           "Asian" = "asian (indian)",
                           "Black" = "black",
                           "Bi- or Multiracial (Unsp.)" = "biracial/multiracial",
                           "Bi- or Multiracial (Unsp.)" = "Native, Black, White",
                           "Biracial (AW)" = "asian/white",
                           "Biracial (AW)" = "biracial (asian/white)",
                           "Biracial (BW)" = "biracial (white/black)",
                           "Biracial (BW)" = "black and white",
                           "Biracial (HW)" = "white, hispanic",
                           "Hispanic or Latinx" = "hispanic or latino",
                           "White" = "white"))

weco_data_final_bytrial %>% group_by(age, condition_framing) %>% distinct(part_id) %>% summarize(n = n()) # need 24 per age-condition

weco_data_fullexclude0 <- weco_data_final_bytrial %>% 
  group_by(part_id) %>% 
  mutate(fully_excluded = ifelse(sum(final_sample) == 0, "yes", "no")) %>%  # both trials of all final sample participants included, so this is appropriate
  ungroup()

weco_data_fullexclude0 %>% group_by(fully_excluded) %>% distinct(part_id) %>% summarize(n = n()) # n = 5 participants had trials fully excluded

weco_data_fullexclude1 <- weco_data_fullexclude0 %>% 
  filter(fully_excluded != "yes") 
  
weco_data_final_id <- weco_data_fullexclude1 %>% # cull excess participants from participants with at least 1 trial included
  mutate(across(!age_years, as.factor)) %>% 
  distinct(part_id, age, condition_framing) %>%
  group_by(age, condition_framing) %>% 
  reframe(n = row_number(),
          part_id = part_id) %>%  # cf. summarize (deprecated)
  filter(n < 25) %>%  # cull excess participants
  ungroup() %>% 
  drop_na() # remove 1 participant (weco_final74) who did not complete consent

weco_data_final_id %>% group_by(age, condition_framing) %>% summarize(n = n()) # 24 per age-condition

weco_data_final0 <- weco_data_final_bytrial %>% # 192 test trials
  filter(part_id %in% weco_data_final_id$part_id) #%>% # retain only non-excess participants

weco_data_final0 %>% group_by(final_sample) %>% summarize(n = n()) # 14 excluded final sample trials
weco_data_final0 %>% group_by(final_sample, age) %>% summarize(n = n()) # 6 excluded final sample trials from 2-year-olds
weco_data_final0 %>% group_by(final_sample, condition_framing) %>% summarize(n = n()) # 11 excluded final sample trials following "we"-framing
weco_data_final0 %>% group_by(final_sample, condition_partner) %>% summarize(n = n()) # 4 excluded final sample trials in group contexts

weco_data_final <- weco_data_final0 %>% # 178 final sample included test trials from 96 participants
  filter(final_sample == "1") %>% # retain only final sample included trials
  select(-final_sample)

weco_data_final %>% summarize(n = n()) # 178 final sample trials
weco_data_final %>% distinct(part_id, age, sex) %>% group_by(age, sex) %>% summarize(n = n()) # sex breakdown: n = 24 female 2-yo's, n = 16 female 3-yo's
weco_data_final %>% distinct(part_id, e1, e2, condition_framing) %>% group_by(e1, condition_framing) %>% summarize(n = n())
weco_data_final %>% distinct(part_id, e1, e2, condition_framing) %>% group_by(e2) %>% summarize(n = n())
weco_data_final %>% distinct(part_id, e1, e2, condition_framing) %>% group_by(e1, e2) %>% summarize(n = n())

weco_data_final %>% # age breakdown
  distinct(part_id, age, age_years) %>%
  group_by(age) %>%
  summarize(n = n(),
            min  = min(age_years),
            mean = mean(age_years),
            max = max(age_years))

weco_data_final %>% distinct(part_id) # 96 participants
weco_data_final %>% group_by(age)                                                  %>% summarize(n = n ())
weco_data_final %>% group_by(age)                            %>% distinct(part_id) %>% summarize(n = n ())
weco_data_final %>% group_by(age, abandon)                   %>% distinct(part_id) %>% summarize(n = n ())
weco_data_final %>% group_by(age, condition_framing, prompt)                       %>% summarize(n = n()) %>% ungroup() %>% complete(age, condition_framing, prompt, fill = list(n = 0))
weco_data_final %>% group_by(age, condition_framing)                               %>% summarize(n = n ())
weco_data_final %>% group_by(age, condition_framing)         %>% distinct(part_id) %>% summarize(n = n ())
weco_data_final %>% group_by(age, condition_partner)         %>% distinct(part_id) %>% summarize(n = n ())
weco_data_final %>% group_by(age, order_partner)             %>% distinct(part_id) %>% summarize(n = n ())
weco_data_final %>% group_by(age, order_game)                %>% distinct(part_id) %>% summarize(n = n ())
weco_data_final %>% group_by(age, pup)                       %>% distinct(part_id) %>% summarize(n = n ())
weco_data_final %>% group_by(age, sex)                       %>% distinct(part_id) %>% summarize(n = n ())
weco_data_final %>% group_by(income)                         %>% distinct(part_id) %>% summarize(n = n ())
weco_data_final %>% group_by(race)                           %>% distinct(part_id) %>% summarize(n = n ())

## reliability ####

weco_data_reli_prompt <- weco_data %>% 
  select(part_id, prompt) %>% 
  filter(part_id %in% weco_reli_prompt$part_id_jv)

reli_dat_prompt <- cbind(weco_reli_prompt, weco_data_reli_prompt) %>% 
  mutate(prompt = as.numeric(prompt),
         disagree = ifelse(prompt == prompt_jv, 0, 1))

weco_kappa_prompt <- reli_dat_prompt %>%
  select(prompt, prompt_jv) %>% 
  mutate_at(vars(everything()), as.double)

kappa2(weco_kappa_prompt)
agree(weco_kappa_prompt)

## create df ####

weco_data_final %>% group_by(prompt) %>% summarize(n = n())
weco_data_final %>% group_by(abandon) %>% summarize(n = n())

data_prompt <- weco_data_final %>% # 173 obs (5 excluded, see below)
  mutate(prompt = fct_recode(prompt, # aggregating to 2+ prompt factor (latency)
                             "2+" = "2",
                             "2+" = "3",
                             "2+" = "4"),
         trial_num = as.factor(trial_num)) %>% 
  mutate(prompt_recoded1 = ifelse(prompt == 95, "2+", paste(prompt)),
         prompt_recoded2 = ifelse(prompt == 95, "remained", paste(prompt)),
         prompt_recoded2 = fct_recode(prompt_recoded2,
                                      "remained" = "remained",
                                      "first" = "1",
                                      "second_plus" = "2+"),
         prompt = ifelse(prompt == 95, NA, paste(prompt)),
         e2_original = e2,
         e2 = fct_recode(e2, "A" = "arod",
                         "B" = "eleanor",
                         "B" = "isabelle",
                         "C" = "jadelyn",
                         "D" = "kaelin",
                         "B" = "sabrina")) %>%
  mutate_if(is.character, as.factor)

data_prompt %>% filter(prompt != "NA") %>% distinct(part_id)
data_prompt %>% filter(prompt != "NA") %>% group_by(prompt) %>% summarize(n = n())
data_prompt %>% filter(prompt != "NA") %>% group_by(prompt, trial_num) %>% summarize(n = n())

data_prompt %>%
  group_by(part_id) %>%
  mutate(trials = n()) %>%
  ungroup %>%
  filter(is.na(prompt)) %>%
  group_by(part_id) %>%
  reframe(trials_excluded = n(),
            ratio = trials_excluded/trials) %>%
  distinct() # n = 4 participants excluded with n = 5 trials (96 and 99 fully excluded)

data_prompt %>% filter_all(any_vars(!is.na(prompt))) # 173 trials

# analyses ####

## priors ####

prior_bytrial <- c(prior(normal(0, 1.50), class = Intercept),
                   prior(normal(0, 1.00), class = b),
                   prior(student_t(10, 0, 0.1), class = sd)) # lots of prior info or else divergences

prior_byparts <- c(prior(normal(0, 1.50), class = Intercept),
                   prior(normal(0, 1.00), class = b))

prior_intercept_bytrial <- c(prior(normal(0, 1.50), class = Intercept),
                             prior(student_t(10, 0, 0.1), class = sd))

prior_categorical_byparts <- c(prior(normal(0, 1.50), class = Intercept, dpar = muSometimes), # categorical(refcat = "Never")
                               prior(normal(0, 1.50), class = Intercept, dpar = muAlways),
                               prior(normal(0, 1.00), class = b, dpar = muSometimes),
                               prior(normal(0, 1.00), class = b, dpar = muAlways))

prior_categorical_bytrial <- c(prior(normal(0, 1.50), class = Intercept, dpar = mufirst), # categorical(refcat = "remained")
                               prior(normal(0, 1.50), class = Intercept, dpar = musecondplus),
                               prior(normal(0, 1.00), class = b, dpar = mufirst),
                               prior(normal(0, 1.00), class = b, dpar = musecondplus),
                               prior(student_t(10, 0, 0.1), class = sd, dpar = mufirst),
                               prior(student_t(10, 0, 0.1), class = sd, dpar = musecondplus))

## main analysis ####

modc <- brm_func(data = data_prompt, prior = prior_bytrial, formula = prompt ~ sex + e1 + e2 + trial_num + order_game + pup + game     + (condition_partner | part_id))

table_modc <- table_func(modc)

mod1 <- brm_func(data = data_prompt, prior = prior_bytrial, formula = prompt ~ trial_num + age                                         + (condition_partner | part_id))
mod2 <- brm_func(data = data_prompt, prior = prior_bytrial, formula = prompt ~ trial_num + age + condition_partner                     + (condition_partner | part_id))
mod3 <- brm_func(data = data_prompt, prior = prior_bytrial, formula = prompt ~ trial_num + age * condition_partner                     + (condition_partner | part_id))
mod4 <- brm_func(data = data_prompt, prior = prior_bytrial, formula = prompt ~ trial_num + age + condition_framing                     + (condition_partner | part_id))
mod5 <- brm_func(data = data_prompt, prior = prior_bytrial, formula = prompt ~ trial_num + age * condition_framing                     + (condition_partner | part_id))
mod6 <- brm_func(data = data_prompt, prior = prior_bytrial, formula = prompt ~ trial_num + age + condition_partner + condition_framing + (condition_partner | part_id))
mod7 <- brm_func(data = data_prompt, prior = prior_bytrial, formula = prompt ~ trial_num + age * condition_framing * condition_partner + (condition_partner | part_id))

list_mod <- list(mod1,
                 mod2,
                 mod3,
                 mod4,
                 mod5,
                 mod6,
                 mod7)

names_mod <- c("trial + age",
               "trial + age + partner",
               "trial + age * partner",
               "trial + age + framing",
               "trial + age * framing",
               "trial + age + partner + framing",
               "trial + age * framing * partner")

names(list_mod) <- names_mod

for(i in list_mod){
  plot(pp_check(i, type = "dens_overlay", ndraws = 100) +
         ggtitle(paste(i[[1]])))
  }

loo_weight_mod <- lapply(list_mod, loo, reloo = TRUE)
round(loo_model_weights(setNames(loo_weight_mod, c(names(loo_weight_mod))), weights = "stacking"), 2)

table_mod7 <- table_func(mod7)

## additional analyses ####

### first set ####

# within condition

mod8  <- brm_func(data = data_prompt %>% filter(condition_partner == "Group"), prior = prior_byparts, formula = prompt ~ trial_num + age * condition_framing)
mod9  <- brm_func(data = data_prompt %>% filter(condition_partner == "Dyad"),  prior = prior_byparts, formula = prompt ~ trial_num + age * condition_framing)
mod10 <- brm_func(data = data_prompt %>% filter(condition_framing == "You"),   prior = prior_bytrial, formula = prompt ~ trial_num + age * condition_partner + (condition_partner | part_id))
mod11 <- brm_func(data = data_prompt %>% filter(condition_framing == "We"),    prior = prior_bytrial, formula = prompt ~ trial_num + age * condition_partner + (condition_partner | part_id))

table_mod8  <- table_func(mod8)
table_mod9  <- table_func(mod9)
table_mod10 <- table_func(mod10)
table_mod11 <- table_func(mod11)

plot_func(df = data_prompt,
          condition = condition_partner,
          facet_x = condition_framing) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "Partner Condition", title = "Partner Effects by Framing")

plot_func(df = data_prompt,
          condition = condition_framing,
          facet_x = condition_partner) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "Framing Condition", title = "Framing Effects by Partner")

### second set ####

# within age

mod12  <- brm_func(data = data_prompt %>% filter(age == "2-year-olds"), prior = prior_intercept_bytrial, formula = prompt ~ trial_num +                                         (condition_partner | part_id))
mod13  <- brm_func(data = data_prompt %>% filter(age == "2-year-olds"), prior = prior_bytrial,           formula = prompt ~ trial_num +                     condition_partner + (condition_partner | part_id))
mod14  <- brm_func(data = data_prompt %>% filter(age == "2-year-olds"), prior = prior_bytrial,           formula = prompt ~ trial_num + condition_framing +                     (condition_partner | part_id))
mod15  <- brm_func(data = data_prompt %>% filter(age == "2-year-olds"), prior = prior_bytrial,           formula = prompt ~ trial_num + condition_framing + condition_partner + (condition_partner | part_id))
mod16  <- brm_func(data = data_prompt %>% filter(age == "2-year-olds"), prior = prior_bytrial,           formula = prompt ~ trial_num + condition_framing * condition_partner + (condition_partner | part_id))

mod17  <- brm_func(data = data_prompt %>% filter(age == "3-year-olds"), prior = prior_intercept_bytrial, formula = prompt ~ trial_num                                         + (condition_partner | part_id))
mod18  <- brm_func(data = data_prompt %>% filter(age == "3-year-olds"), prior = prior_bytrial,           formula = prompt ~ trial_num +                     condition_partner + (condition_partner | part_id))
mod19  <- brm_func(data = data_prompt %>% filter(age == "3-year-olds"), prior = prior_bytrial,           formula = prompt ~ trial_num + condition_framing +                     (condition_partner | part_id))
mod20  <- brm_func(data = data_prompt %>% filter(age == "3-year-olds"), prior = prior_bytrial,           formula = prompt ~ trial_num + condition_framing + condition_partner + (condition_partner | part_id))
mod21  <- brm_func(data = data_prompt %>% filter(age == "3-year-olds"), prior = prior_bytrial,           formula = prompt ~ trial_num + condition_framing * condition_partner + (condition_partner | part_id))

list_mod2 <- list(mod12,
                  mod13,
                  mod14,
                  mod15,
                  mod16)
list_mod3 <- list(mod17,
                  mod18,
                  mod19,
                  mod20,
                  mod21)

names_mod2 <- c("trial",
                "trial + partner",
                "trial + framing",
                "trial + framing + partner",
                "trial + framing * partner")

names(list_mod2) <- names_mod2
names(list_mod3) <- names_mod2

for(i in list_mod2){
  plot(pp_check(i, type = "dens_overlay", ndraws = 100) +
         ggtitle(paste(i[[1]])))
}
for(i in list_mod3){
  plot(pp_check(i, type = "dens_overlay", ndraws = 100) +
         ggtitle(paste(i[[1]])))
}

loo_weight_mod2 <- lapply(list_mod2, loo, reloo = TRUE)
round(loo_model_weights(setNames(loo_weight_mod2, c(names(loo_weight_mod2))), weights = "stacking"), 2)

loo_weight_mod3 <- lapply(list_mod3, loo, reloo = TRUE)
round(loo_model_weights(setNames(loo_weight_mod3, c(names(loo_weight_mod3))), weights = "stacking"), 2)

table_mod16 <- table_func(mod16)
table_mod21 <- table_func(mod21)

# within 2-yo's, within partner condition

mod14a  <- brm_func(data = data_prompt %>% filter(condition_partner == "Group", age == "2-year-olds"), prior = prior_byparts, formula = prompt ~ trial_num + condition_framing)
mod14b  <- brm_func(data = data_prompt %>% filter(condition_partner == "Dyad", age == "2-year-olds"), prior = prior_byparts, formula = prompt ~ trial_num + condition_framing)

table_func(mod14a)
table_func(mod14b)

### third set ####

# by participant

data_prompt_byparts <- data_prompt %>%
  filter(!is.na(prompt)) %>% 
  group_by(part_id) %>%
  mutate(n = n(),
         prompt_bypart = ifelse(prompt == "2+", 1, 0),
         sum = sum(prompt_bypart),
         ratio = sum/n,
         ratio = as.factor(ratio),
         ratio = fct_relevel(ratio, c("1", "0.5", "0")),
         ratio = fct_recode(ratio, "Always" = "1", "Sometimes" = "0.5", "Never" = "0")) %>% # participants always left 
  distinct(part_id, .keep_all = TRUE) %>%
  ungroup()

data_prompt_byparts %>% group_by(ratio, sum) %>% summarise(n = n())

plot_func(byparts, df = data_prompt_byparts, condition = age) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "Age")

plot_func(byparts, df = data_prompt_byparts, condition = condition_framing, facet_x = age) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "Framing")

mod22  <- brm_func(categorical(refcat = "Never"), data = data_prompt_byparts, prior = prior_categorical_byparts, formula = ratio ~ age)
mod23  <- brm_func(categorical(refcat = "Never"), data = data_prompt_byparts, prior = prior_categorical_byparts, formula = ratio ~ age + condition_framing)
mod24  <- brm_func(categorical(refcat = "Never"), data = data_prompt_byparts, prior = prior_categorical_byparts, formula = ratio ~ age * condition_framing)

list_mod4 <- list(mod22,
                  mod23,
                  mod24)

names_mod4 <- c("age",
                "age + framing",
                "age * framing")

names(list_mod4) <- names_mod4

for(i in list_mod4){
  plot(pp_check(i, type = "dens_overlay", ndraws = 100) +
         ggtitle(paste(i[[1]])))
}

loo_weight_mod4 <- lapply(list_mod4, loo, reloo = TRUE)
round(loo_model_weights(setNames(loo_weight_mod4, c(names(loo_weight_mod4))), weights = "stacking"), 2)

table_mod24 <- table_func(mod24)

# age models

mod25 <- brm_func(categorical(refcat = "Never"), data = data_prompt_byparts %>% filter(age == "2-year-olds"), prior = prior_categorical_byparts, formula = ratio ~ condition_framing)
table_mod25 <- table_func(mod25)

mod26  <- brm_func(categorical(refcat = "Never"), data = data_prompt_byparts %>% filter(age == "3-year-olds"), prior = prior_categorical_byparts, formula = ratio ~ condition_framing)
table_mod26 <- table_func(mod26)

### fourth set ####

# full dataset

# same model

# assign *maximum* value to those who remained with partner(s)

mod28 <- brm_func(data = data_prompt, prior = prior_bytrial, formula = prompt_recoded1 ~ trial_num + age * condition_framing * condition_partner + (condition_partner | part_id))

table_mod28 <- table_func(mod28)

# different model

# assign *unique* value to those who remained with partner(s)

mod29 <- brm_func(categorical(refcat = "remained"), data = data_prompt, prior = prior_categorical_bytrial, formula = prompt_recoded2 ~ trial_num + age * condition_framing * condition_partner + (condition_partner | part_id))

table_mod29 <- table_func(mod29)

# figures ####

## figure 2 ####

fig2a <- plot_func(df = data_prompt,
                   condition = condition_framing) +
  labs(x = "Framing Condition") +
  theme(axis.title.x = element_text(size = 19))
fig2b <- plot_func(df = data_prompt,
                   condition = condition_partner) +
  labs(x = "Partner Number Condition", y = NULL) +
  theme(axis.title.x = element_text(size = 19))

fig2a + fig2b +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'A',
                  title = "Effects of Framing and Partner Number",
                  theme = theme(plot.title = element_text(size = 24,
                                                          family = "Times New Roman",
                                                          hjust = .50)))

# ggsave("figures/figure2.jpeg", height = 4, width = 8)

## figure 3 ####

# first set of additional analyses

fig3a <- plot_func(df = data_prompt %>%
                     filter(condition_partner == "Group") %>%
                     mutate(condition_partner = fct_recode(condition_partner,
                                                            "Group Partner" = "Group")),
                   condition = condition_framing) +
  labs(x = NULL) +
  facet_grid(~condition_partner)
fig3b <- plot_func(df = data_prompt %>%
                     filter(condition_partner == "Dyad") %>%
                     mutate(condition_partner = fct_recode(condition_partner,
                                                           "Dyadic Partner" = "Dyad")),
                   condition = condition_framing) +
  labs(x = "Framing Condition") +
  theme(axis.title.x = element_text(size = 19)) +
  facet_grid(~condition_partner)
fig3c <- plot_func(df = data_prompt %>%
                     filter(condition_framing == "You") %>%
                     mutate(condition_framing = fct_recode(condition_framing,
                                                           "You-framing" = "You")),
                   condition = condition_partner) +
  labs(x = NULL, y = NULL) +
  facet_grid(~condition_framing)
fig3d <- plot_func(df = data_prompt %>%
                     filter(condition_framing == "We") %>%
                     mutate(condition_framing = fct_recode(condition_framing,
                                                           "We-framing" = "We")),
                   condition = condition_partner) +
  labs(x = "Partner Number Condition", y = NULL) +
  facet_grid(~condition_framing)

fig3a + fig3c + fig3b + fig3d +
  plot_layout(guides = "collect") +
  plot_annotation(title = "First Set of Additional Analyses",
                  tag_levels = 'A',
                  theme = theme(plot.title = element_text(size = 28, family = "Times New Roman", hjust = .50)))

# ggsave("figures/figure3.jpeg", height = 8, width = 9)

## figure 4 ####

plot_func(df = data_prompt %>% filter(age == "2-year-olds"),
                   condition = condition_framing,
          facet_x = condition_partner) +
  labs(x = "Framing Condition", title = "Effects of Framing on 2-year-olds")

# ggsave("figures/figure4.jpeg", height = 4, width = 6)

## figure 5 ####

plot_func(by_parts,
          df = data_prompt_byparts %>%
            filter(age == "2-year-olds"),
          condition = condition_framing) +
  labs(x = "Framing Condition", title = "Effects of Framing on 2-year-olds' Ratios")

# ggsave("figures/figure5.jpeg", height = 4, width = 7)

# appendix ####

## table 1 (formulas) ####

app_table1 <- tibble(Formula = c(as.character(mod1$formula[1]), # main analysis
                                 as.character(mod2$formula[1]), # main analysis
                                 as.character(mod3$formula[1]), # main analysis
                                 as.character(mod4$formula[1]), # main analysis
                                 as.character(mod5$formula[1]), # main analysis
                                 as.character(mod6$formula[1]), # main analysis
                                 as.character(mod7$formula[1]), # main analysis
                                 as.character(mod12$formula[1]), # by-age analysis
                                 as.character(mod13$formula[1]), # by-age analysis
                                 as.character(mod14$formula[1]), # by-age analysis
                                 as.character(mod15$formula[1]), # by-age analysis
                                 as.character(mod16$formula[1]), # by-age analysis
                                 as.character(mod22$formula[1]), # by-participant analysis
                                 as.character(mod23$formula[1]), # by-participant analysis
                                 as.character(mod24$formula[1])), # by-participant analysis
                     Analysis = c(rep("Main", 7),
                                  rep("2nd Set", 5),
                                  rep("3rd Set", 3))) %>% 
  mutate(Formula = str_replace_all(Formula, pattern = "_", replacement = " "),
         Formula = str_replace_all(Formula, pattern = "condition ", replacement = ""),
         Formula = str_replace_all(Formula, pattern = "part id", replacement = "participant"),
         Formula = as_tibble(str_split_fixed(Formula, pattern = "~", 2))) %>% # use as_tibble is necessary because str_split_fixed creates a character matrix (indicated by colname[,1], etc)
  unnest_wider(Formula) %>% # use as_tibble, above, enables unnesting of character matrix created by str_split_fixed
  rename("Outcome" = V1, "Formula" = V2) %>% 
  relocate(Outcome, .before = everything()) %>% 
  relocate(Analysis, .before = everything())

# write_csv(app_table1, "tables/app_table1.csv")

## table 2 (pars) ####

app_table2 <- rbind(table_modc,
                    table_mod7,
                    table_mod8,
                    table_mod9,
                    table_mod10,
                    table_mod11,
                    table_mod16,
                    table_mod21,
                    table_mod24,
                    table_mod25,
                    table_mod26,
                    table_mod28,
                    table_mod29) %>% 
  mutate(Analysis = c(rep("Control", dim(table_modc)[1]), # main analysis, control
                      rep("Main", dim(table_mod7)[1]), # main analysis
                      rep("Set 1 (Group)", dim(table_mod8)[1]), # 1st addit'l analyses
                      rep("Set 1 (Dyad)", dim(table_mod9)[1]), # 1st addit'l analyses
                      rep("Set 1 (You)", dim(table_mod10)[1]), # 1st addit'l analyses
                      rep("Set 1 (We)", dim(table_mod11)[1]), # 1st addit'l analyses
                      rep("Set 2 (2-yo)", dim(table_mod16)[1]), # 2nd addit'l analyses
                      rep("Set 2 (3-yo)", dim(table_mod21)[1]), # 2nd addit'l analyses
                      rep("Set 3", dim(table_mod24)[1]), # 3rd addit'l analyses
                      rep("Ser 3 (2-yo)", dim(table_mod25)[1]), # 3rd addit'l analyses (within 2-yo's)
                      rep("Set 3 (3-yo)", dim(table_mod26)[1]), # 3rd addit'l analyses (within 3-yo's)
                      rep("Set 4 (R1)", dim(table_mod28)[1]),
                      rep("Set 4 (R2)", dim(table_mod29)[1])), # 4th addit'l analyses
         Parameter = str_replace_all(Parameter, ":", " * "),
         distr_par = ifelse(str_detect(Parameter, "mu"), paste(word(Parameter, 1)), NA),
         Parameter = ifelse(str_starts(Parameter, "mu"), str_replace(Parameter, "^[^ ]+\\s", ""), Parameter), # regex for identifying start of string (^), one or more non-space charaters ([^ ]+), and first space (\\s)
         Parameter = ifelse(str_starts(Parameter, "e"), Parameter, str_replace_all(Parameter, "1", "")),
         Parameter = ifelse(Parameter  %in% c("e11", "e21"), Parameter, str_replace_all(Parameter, "1$", "")),
         Parameter = str_replace_all(Parameter, "condition ", ""),
         Parameter = str_to_lower(Parameter),
         distr_par = ifelse(str_starts(distr_par, "mu"), str_replace_all(distr_par, "mu", ""), distr_par)) %>% 
  relocate(distr_par, .before = everything()) %>% 
  rename("Dis. Par." = distr_par) %>% 
  relocate(Analysis, .before = everything()) %>% 
  select(-c("Error", "Ev. Ratio"))
    
# write_csv(app_table2, "tables/app_table2.csv")

## table 3 ####

# mods 1-7 / main analysis

postprob_mods17 <- post_prob(mod1,
                             mod2,
                             mod3,
                             mod4,
                             mod5,
                             mod6,
                             mod7)

round(postprob_mods17, 3)

table_postprob_mods17 <- tibble("postprob_mods17" = postprob_mods17,
                                "bf_mods17" = postprob_mods17/postprob_mods17[1],
                                Model = names_mod) %>% 
  mutate(across(where(is.numeric), \(x) round(x, 3)),
         across(where(is.numeric), \(x) ifelse(x == 1, NA, x))) %>% 
  relocate(Model, .before = everything()) %>% 
  rename("Post. Prob." = "postprob_mods17",
         "BF" = "bf_mods17")

# write_csv(table_postprob_mods17, "tables/app_table3.csv")

## table 4 ####

# mods 12-6 & 17-21 / second set of additional analyses

postprob_mods1216 <- post_prob(mod12,
                               mod13,
                               mod14,
                               mod15,
                               mod16)

postprob_mods1721 <- post_prob(mod17,
                               mod18,
                               mod19,
                               mod20,
                               mod21)

round(postprob_mods1216, 3)
round(postprob_mods1721, 3)

table_postprob_mods1221 <- tibble("postprob_mods1216" = postprob_mods1216,
                                  "bf_mods1216" = postprob_mods1216/postprob_mods1216[1],
                                  "postprob_mods1721" = postprob_mods1721,
                                  "bf_mods1721" = postprob_mods1721/postprob_mods1721[1],
                                  Model = names_mod2) %>% 
  mutate(across(where(is.numeric), \(x) round(x, 3)),
         across(where(is.numeric), \(x) ifelse(x == 1, NA, x))) %>% 
  relocate(Model, .before = everything()) %>% 
  rename("Post. Prob. (2yr)" = "postprob_mods1216",
         "BF (2yr)" = "bf_mods1216",
         "Post. Prob. (3yr)" = "postprob_mods1721",
         "BF (3yr)" = "bf_mods1721")

# write_csv(table_postprob_mods1221, "tables/app_table4.csv")

## table 5 ####

# mods 22-4 / third set of additional analyses

postprob_mods2224 <- post_prob(mod22,
                               mod23,
                               mod24)

round(postprob_mods2224, 3)

table_postprob_mods2224 <- tibble("postprob_mods2224" = postprob_mods2224,
                                  "bf_mods2224" = postprob_mods2224/postprob_mods2224[1],
                                  Model = names_mod4) %>% 
  mutate(across(where(is.numeric), \(x) round(x, 3)),
         across(where(is.numeric), \(x) ifelse(x == 1, NA, x))) %>% 
  relocate(Model, .before = everything()) %>% 
  rename("Post. Prob." = "postprob_mods2224",
         "BF" = "bf_mods2224")

# write_csv(table_postprob_mods2224, "tables/app_table5.csv")

