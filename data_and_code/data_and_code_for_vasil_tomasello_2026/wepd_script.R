
start <- Sys.time()

# globals ####
library(brms) # models, inference
library(flextable) # tables
library(ggtext) # plots
library(irr) # reliability
library(officer) # tables
library(parallel) # detectcores
library(patchwork) # plots
library(scales) # pretty_breaks
library(tidybayes) # inference
library(tidyverse) # wrangling

set.seed(10003) # rng
options(contrasts = c("contr.treatment", "contr.poly")) # treatment coded model

# data ####
wepd_data1 <- read_csv("data/wepd_data.csv") # import data

wepd_data2 <- wepd_data1 %>%  # clean data
  filter(exclude_reason == "included") %>% 
  mutate(ref_ling_recoded = fct_recode(ref_ling, "no_ref" = "0", "nonling_only" = "no_ling_ref", "proform" = "null"), # collapse null (n = 1) with proform
         form = fct_relevel(factor(ref_ling_recoded, ordered = T), c("no_ref", "nonling_only", "proform", "np_unmod", "np_mod")),
         age_standardized = ((age - mean(age))/(2*sd(age))),
         across(!c(age, age_standardized), as.factor)) %>% 
  group_by(participant_id) %>% 
  mutate(trial = as_factor(row_number())) %>% 
  ungroup()

wepd_data <- wepd_data2 %>% # cleaned df
  select(participant_id, age_standardized, condition_discourse, condition_framing, form, trial, location, sex) %>% 
  filter(form != "no_ref")

## fig 2 ####
data_df_for_gg <- wepd_data %>% # visualize observations
  mutate(across(c(condition_discourse, condition_framing), ~str_to_title(as.character(.))),
         form = fct_rev(form),
         form = case_when(form == "nonling_only" ~ "Nonlinguistic",
                          form == "proform" ~ "Proform",
                          form == "np_unmod" ~ "NP Unmod.",
                          form == "np_mod" ~ "NP Mod."),
         form = fct_relevel(form, "NP Mod.", "NP Unmod.", "Proform", "Nonlinguistic"),
         condition_framing = paste0("''", condition_framing, "''", "-Framing"),
         trial = paste0("Trial ", trial)) %>%
  rename(Form = form)

ggplot(data = data_df_for_gg,
       mapping = aes(x = condition_discourse, fill = Form)) +
  geom_bar(position = "fill") +
  scale_fill_grey(start = 0.85, end = 0.45) +
  facet_grid(trial ~ condition_framing) +
  scale_y_continuous(limits = c(NA, NA),
                     labels = function(x) sub("^0\\.", ".", format(x, nsmall = 2))) +
  geom_text(stat = "count",
            aes(y = after_stat(count),
                label = after_stat(count),
                vjust = ifelse(after_stat(count) < 5, -100, 1.2)),
            position = "fill", size = 4.5, family = "Times New Roman") +
  theme_classic(base_size = 14, base_family = "Times New Roman") +
  labs(x = "Discourse Condition",
       y = "Observed Proportion", 
       title = "Observations") +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_text(vjust = 0, size = 13),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
ggsave("figs/fig2.png", width = 7, height = 4)

# mod ####
mod <- brm(family  = cumulative(),
           data    = wepd_data,
           formula = form ~ cs(condition_discourse) + cs(condition_framing) + condition_discourse * condition_framing + age_standardized + location + trial + sex + (trial | participant_id),
           prior   = c(set_prior("normal(0, 1.00)", "Intercept"),
                       set_prior("normal(0, 0.80)", "b"),
                       set_prior("normal(0, 0.20)", "sd")),
           control = list(adapt_delta = .99),
           chains  = detectCores() - 1,
           cores   = detectCores() - 1,
           init    = 0,
           backend = "cmdstanr")

# ppc's ####
# graphical
pp_check(mod, "dens_overlay_grouped", group = "condition_discourse", ndraws = 1000)
pp_check(mod, "dens_overlay_grouped", group = "condition_framing", ndraws = 1000)
pp_check(mod, "ecdf_overlay_grouped", group = "condition_discourse", ndraws = 1000)
pp_check(mod, "ecdf_overlay_grouped", group = "condition_framing", ndraws = 1000)
pp_check(mod, "bars_grouped", group = "condition_discourse", ndraws = 1000)
pp_check(mod, "bars_grouped", group = "condition_framing", ndraws = 1000)

# p values
# overall
yrep <- posterior_predict(mod) - 1
yobs <- mod$data$form # observations
str(yobs) # factor outcome
yobs <- as.numeric(yobs) - 1 # convert from factor to numeric
T_obs <- mean(yobs) # observed proportion of stereotypical responses
T_obs_var <- var(yobs) # observed variance of stereotypical responses
T_rep <- apply(yrep, 1, mean) # simulated proportions
T_rep_var <- apply(yrep, 1, var) # simulated proportions
mean(T_rep >= T_obs) # proportion of simulated datasets with means that exceed that observed
mean(T_rep_var >= T_obs_var) # proportion of simulated datasets with sd's that exceed that observed

# each outcome
K <- length(levels(mod$data$form))
obs_prop <- sapply(0:(K-1), function(k) mean(yobs == k))
obs_prop
rep_prop <- t(apply(yrep, 1, function(row) {
  sapply(0:(K-1), function(k) mean(row == k))
}))
bayes_p <- sapply(1:K, function(k) mean(rep_prop[, k] >= obs_prop[k]))
bayes_p

# condition * outcome
mod$data$cond <- paste(mod$data$condition_discourse, mod$data$condition_framing, sep = "_")
conditions <- unique(mod$data$cond)
bayes_p_cond <- matrix(NA, nrow = length(conditions), ncol = K,
                       dimnames = list(conditions, paste0("cat_", 0:(K-1))))
for(cond in conditions){
  idx <- which(mod$data$cond == cond)
  yobs_sub <- yobs[idx]
  obs_prop_sub <- sapply(0:(K-1), function(k) mean(yobs_sub == k))
  rep_prop_sub <- t(apply(yrep[, idx, drop = FALSE], 1, function(row) {
    sapply(0:(K-1), function(k) mean(row == k))
  }))
  bayes_p_cond[cond, ] <- sapply(1:K, function(k) mean(rep_prop_sub[, k] >= obs_prop_sub[k]))
}
bayes_p_cond

# cumulative probabilities
cum_prob <- function(x) sapply(0:(K-1), function(k) mean(x <= k))
obs_cum <- cum_prob(yobs)
rep_cum <- t(apply(yrep, 1, cum_prob))
bayes_p_cum <- sapply(1:K, function(k) mean(rep_cum[, k] >= obs_cum[k]))
obs_cum
bayes_p_cum

# inference ####
newdata <- crossing(condition_discourse = levels(mod$data$condition_discourse), # prediction grid
                    condition_framing = levels(mod$data$condition_framing),
                    location = levels(mod$data$location),
                    trial = levels(mod$data$trial),
                    sex = levels(mod$data$sex),
                    age_standardized = seq(from = min(mod$data$age_standardized),
                                           to = max(mod$data$age_standardized),
                                           length.out = 20))

epreds <- epred_draws(mod, newdata, re_formula = NA, seed = 27707) # rows = 960 (cells) * 7000 (post warmup draws) * 4 (outcome categories)

## expected value ####
# ie., expected informativeness
epreds_expected_value <- epreds %>%
  mutate(.category_numeric = as.numeric(.category)) %>% # ranges 1-4
  group_by(.draw, condition_discourse, condition_framing, age_standardized, location, trial, sex) %>%
  summarise(expected_value = sum(.category_numeric * .epred), .groups = "drop") # expected informativeness per predictor combination per draw

### fig 3 ####
# discourse contrasts
contrasts_discourse_marginal <- epreds_expected_value %>% # draws of marginal comparison of discourse
  group_by(.draw, condition_discourse) %>% 
  summarize(mean_expected_value = mean(expected_value),
            .groups = "drop") %>% 
  mutate(condition_discourse = factor(condition_discourse, levels = c("idiosyncratic", "conventional"))) %>% 
  group_by(condition_discourse) %>%
  compare_levels(mean_expected_value, by = condition_discourse) %>%
  ungroup() %>% 
  rename(comparison = condition_discourse,
         diff_mean_expected_value = mean_expected_value)
contrasts_discourse_marginal_summary <- contrasts_discourse_marginal %>% # summary of marginal comparison
  group_by(comparison) %>% 
  summarize(hdi_lo = hdi(diff_mean_expected_value)[ncol(hdi(diff_mean_expected_value)) - 1], # rightmost HDI LB to safeguard against multimodality
            mean_expected_value = mean(diff_mean_expected_value),
            hdi_hi = hdi(diff_mean_expected_value)[ncol(hdi(diff_mean_expected_value))], # rightmost HDI UB to safeguard against multimodality
            prop_gt = mean(diff_mean_expected_value > 0),
            prop_gt_char = sprintf("%.0f%%", 100 * prop_gt),
            .groups = "drop")

contrasts_discourse_conditional <- epreds_expected_value %>% # draws of marginal comparison of discourse conditional on framing
  group_by(.draw, condition_discourse, condition_framing) %>% 
  summarize(mean_expected_value = mean(expected_value),
            .groups = "drop") %>% 
  mutate(condition_discourse = factor(condition_discourse, levels = c("idiosyncratic", "conventional"))) %>% 
  group_by(condition_framing) %>%
  compare_levels(mean_expected_value, by = condition_discourse) %>%
  ungroup() %>% 
  rename(comparison = condition_discourse,
         diff_mean_expected_value = mean_expected_value)
contrasts_discourse_conditional_summary <- contrasts_discourse_conditional %>% # summary of conditional marginal comparison
  group_by(condition_framing, comparison) %>% 
  summarize(hdi_lo = hdi(diff_mean_expected_value)[ncol(hdi(diff_mean_expected_value)) - 1], # rightmost HDI LB to safeguard against multimodality
            mean_expected_value = mean(diff_mean_expected_value),
            hdi_hi = hdi(diff_mean_expected_value)[ncol(hdi(diff_mean_expected_value))], # rightmost HDI UB to safeguard against multimodality
            prop_gt = mean(diff_mean_expected_value > 0),
            prop_gt_char = sprintf("%.0f%%", 100 * prop_gt),
            .groups = "drop")

contrasts_discourse_draws <- bind_rows(contrasts_discourse_marginal %>% mutate(condition_framing = "overall"), contrasts_discourse_conditional) %>% 
  mutate(across(starts_with("condition"), ~str_to_title(.)),
         condition_framing = ifelse(condition_framing == "Overall", condition_framing, paste0('"', condition_framing, '"-Framing')),
         condition_framing = fct_relevel(condition_framing, '"You"-Framing', '"We"-Framing', "Overall")) # overall column first
contrasts_discourse_summary <- bind_rows(contrasts_discourse_marginal_summary %>% mutate(condition_framing = "overall"), contrasts_discourse_conditional_summary) %>% 
  mutate(hjust_label = ifelse(hdi_hi > -0.010 & hdi_hi < 0.010, 0.2, NA),
         across(starts_with("condition"), ~str_to_title(.)),
         condition_framing = ifelse(condition_framing == "Overall", condition_framing, paste0('"', condition_framing, '"-Framing')),
         condition_framing = fct_relevel(condition_framing, '"You"-Framing', '"We"-Framing', "Overall")) # overall column first

ggplot(contrasts_discourse_draws,
       aes(x = diff_mean_expected_value, y = condition_framing)) +
  stat_halfeye(aes(fill = after_stat(x < 0)), .width = 0.95, slab_alpha = 0.6, slab_linewidth = 0.5, slab_color = "black") +
  scale_fill_manual(values = c("TRUE" = "grey70", "FALSE" = "grey30"), guide = "none") + 
  geom_text(data = contrasts_discourse_summary,
            aes(x = hdi_hi + .02, label = prop_gt_char, hjust = hjust_label),
            size = 4,
            vjust = -2.5,
            family = "Times New Roman",
            fontface = "bold",
            color = "grey40") +
  geom_vline(xintercept = 0, alpha = 0.3, linetype = 2) +
  labs(x = "Difference in Predicted Informativeness",
       title = "Discourse Contrasts",
       subtitle = "Are participants more informative after conventional or idiosyncratic discourse?",
       y = NULL) +
  scale_x_continuous(labels = function(x) format(x, trim = TRUE),
                     breaks = pretty_breaks(3)) +
  scale_y_discrete(expand = expansion(mult = c(0.05, .8))) +
  theme_classic(base_size = 15, base_family = "Times New Roman") +
  theme(axis.text.y = element_text(vjust = -2),
        axis.ticks.y = element_blank(),
        plot.title = element_text(face = "bold", hjust = .3),
        plot.subtitle = element_text(face = "italic", hjust = .95, size = 13)) +
  annotate("text", x = -.30, 
           y =  4.2, 
           label = "Idiosyncratic\nDiscourse", hjust = .45, color = "grey70",
           family = "Times New Roman", fontface = "bold", size = 4) +
  annotate("text", x = .50, 
           y = 4.2, 
           label = "Conventional\nDiscourse", hjust = .5, color = "grey50",
           family = "Times New Roman", fontface = "bold", size = 4)
# ggsave("figs/fig3.png", width = 6, height = 3.5)

# framing * discourse
# epreds_expected_value_discourse_framing <- epreds_expected_value %>% # framing * discourse draws df
#   group_by(.draw, condition_discourse, condition_framing) %>% # marginalize over covariates
#   summarize(expected_value = mean(mean), .groups = "drop") # expected informativeness per condition combination per draw
# epreds_expected_value_discourse_framing_summary <- epreds_expected_value_discourse_framing %>% # discourse * framing summary df
#   group_by(condition_discourse, condition_framing) %>% 
#   summarize(ci_lo = hdi(expected_value)[ncol(hdi(expected_value)) - 1], # rightmost HDI LB to safeguard against multimodality
#             mean_expected_value = mean(expected_value),
#             ci_hi = hdi(expected_value)[ncol(hdi(expected_value))], # rightmost HDI UB to safeguard against multimodality
#             .groups = "drop")

## expected probability ####
### fig 4 ####
# step 1: discourse * framing
draws_mean_prob_framing_discourse <- epreds %>%
  mutate(.category_numeric = as.numeric(.category)) %>% # ranges 1-4
  group_by(.draw, condition_framing, condition_discourse, .category) %>% # marginalize across control vars
  summarise(mean_probability = mean(.epred), 
            .groups = "drop")
diff_draws_mean_prob_framing_discourse <- draws_mean_prob_framing_discourse %>% 
  mutate(condition_discourse = factor(condition_discourse, levels = c("idiosyncratic", "conventional"))) %>% 
  group_by(condition_framing, .category) %>% 
  compare_levels(mean_probability, by = condition_discourse) %>% 
  ungroup() %>% 
  rename(comparison = condition_discourse,
         diff_mean_prob = mean_probability)
summary_diff_draws_mean_prob_framing_discourse <- diff_draws_mean_prob_framing_discourse %>% 
  group_by(comparison, condition_framing, .category) %>% 
  summarize(hdi_lo = hdi(diff_mean_prob)[ncol(hdi(diff_mean_prob)) - 1], # rightmost HDI LB to safeguard against multimodality
            mean_diff_mean_prob = mean(diff_mean_prob),
            hdi_hi = hdi(diff_mean_prob)[ncol(hdi(diff_mean_prob))], # rightmost HDI UB to safeguard against multimodality
            prop_gt = mean(diff_mean_prob > 0),
            prop_gt_char = sprintf("%.0f%%", 100 * prop_gt),
            .groups = "drop")

# step 2: discourse
draws_mean_prob_discourse <- epreds %>%
  mutate(.category_numeric = as.numeric(.category)) %>% # ranges 1-4
  group_by(.draw, condition_discourse, .category) %>% # marginalize across control vars
  summarise(mean_probability = mean(.epred), 
            .groups = "drop") 
diff_draws_mean_prob_discourse <- draws_mean_prob_discourse %>% 
  mutate(condition_discourse = factor(condition_discourse, levels = c("idiosyncratic", "conventional"))) %>% 
  group_by(condition_discourse, .category) %>% 
  compare_levels(mean_probability, by = condition_discourse) %>% 
  ungroup() %>% 
  rename(comparison = condition_discourse,
         diff_mean_prob = mean_probability)
summary_diff_draws_mean_prob_discourse <- diff_draws_mean_prob_discourse %>% 
  group_by(comparison, .category) %>% 
  summarize(hdi_lo = hdi(diff_mean_prob)[ncol(hdi(diff_mean_prob)) - 1], # rightmost HDI LB to safeguard against multimodality
            mean_diff_mean_prob = mean(diff_mean_prob),
            hdi_hi = hdi(diff_mean_prob)[ncol(hdi(diff_mean_prob))], # rightmost HDI UB to safeguard against multimodality
            prop_gt = mean(diff_mean_prob > 0),
            prop_gt_char = sprintf("%.0f%%", 100 * prop_gt),
            .groups = "drop")

# step 3: combine df's
contrasts_cat <- bind_rows(diff_draws_mean_prob_discourse %>% mutate(condition_framing = "overall"), diff_draws_mean_prob_framing_discourse) %>% 
  mutate(across(starts_with("condition"), ~str_to_title(.)),
         .category = case_when(.category == "nonling_only" ~ "Nonlinguistic",
                               .category == "proform" ~ "Proform",
                               .category == "np_unmod" ~ "Unmod. NP",
                               .category == "np_mod" ~ "Mod. NP"),
         .category = fct_relevel(.category, "Mod. NP", "Unmod. NP", "Proform", "Nonlinguistic"),
         condition_framing = ifelse(condition_framing == "Overall", condition_framing, paste0('"', condition_framing, '"-Framing')),
         condition_framing = fct_relevel(condition_framing, "Overall")) # overall column first
contrasts_cat_summary <- bind_rows(summary_diff_draws_mean_prob_discourse %>% mutate(condition_framing = "overall"), summary_diff_draws_mean_prob_framing_discourse) %>% 
  mutate(hjust_label = ifelse(hdi_hi > -0.010 & hdi_hi < 0.010, 0.2, NA),
         across(starts_with("condition"), ~str_to_title(.)),
         .category = case_when(.category == "nonling_only" ~ "Nonlinguistic",
                               .category == "proform" ~ "Proform",
                               .category == "np_unmod" ~ "Unmod. NP",
                               .category == "np_mod" ~ "Mod. NP"),
         .category = fct_relevel(.category, "Mod. NP", "Unmod. NP", "Proform", "Nonlinguistic"),
         condition_framing = ifelse(condition_framing == "Overall", condition_framing, paste0('"', condition_framing, '"-Framing')),
         condition_framing = fct_relevel(condition_framing, "Overall")) # overall column first

# step 4: plot
ggplot() + # first plot: left y-axis label
geom_richtext(aes(x = 0, y = 0, label = "<span style='color:grey70'>Idiosyncratic</span><br>
                 <span style='color:grey70'>Discourse</span><br><br>
                 <span style='color:black'>vs.</span><br><br>
                 <span style='color:grey50'>Conventional</span><br>
                 <span style='color:grey50'>Discourse</span>"),
                fill = NA, label.color = NA,
              size = 4.5,
                hjust = 0.5, vjust = 0.5, family = "Times New Roman", fontface = "bold") + 
  theme_void() +
  ggplot(contrasts_cat, aes(x = diff_mean_prob, y = comparison)) + # second plot: densities
  stat_halfeye(aes(fill = after_stat(x < 0)), .width = 0.95,
               slab_alpha = 0.6, slab_linewidth = 0.5, slab_color = "black") +
  scale_fill_manual(values = c("TRUE" = "grey70", "FALSE" = "grey30"), guide = "none") + 
  facet_grid(.category ~ condition_framing, scales = "free") +
  geom_text(data = contrasts_cat_summary,
            aes(x = ifelse(hdi_hi < .2, hdi_hi + .1, hdi_hi + .06), label = prop_gt_char, hjust = hjust_label),
            vjust = -2.5,
            size = 3.75,
            family = "Times New Roman",
            fontface = "bold",
            color = "grey40") +
  geom_vline(xintercept = 0, alpha = 0.3, linetype = 2) +
  theme_classic(base_size = 14, base_family = "Times New Roman") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(face = "italic", size = 13, hjust = .5)) +
  labs(y = NULL, x = "Difference in Predicted Probability of Use",
       subtitle = "Which referring forms are used more often after conventional discourse?") +
  scale_x_continuous(breaks = seq(-0.3, 0.3, by = 0.3)) +
  scale_y_discrete(expand = expansion(mult = c(0.1, 0))) +
  plot_annotation(
    title = "Discourse Contrasts by Referring Form",
    theme = theme(plot.title = element_text(face = "bold", family = "Times New Roman", hjust = .6, size = 17))) +
  ggplot() + # third plot: right y-axis label
  geom_richtext(aes(x = 0, y = 0, label = "Referring Form"),
                fill = NA, label.color = NA, 
                family = "Times New Roman",
                fontface = "bold",
                size = 4.5, angle = 270) +
  theme_void() +
  plot_layout(widths = c(1.2, 6, 0.4))
# ggsave("figs/fig4.png", width = 7.2, height = 5.6)

# supplement ####
## tables ####
### 1 ####
exclusions <- wepd_data1 %>% # exclusions by condition
  filter(exclude_reason != "included") %>% 
  distinct(participant_id, .keep_all = T)

sup_tab1 <- exclusions %>% # exclusions by condition
  group_by(condition_discourse, condition_framing) %>% 
  summarize(Excluded = n(), 
            .groups = "drop") %>% 
  mutate(Proportion = sub("^0", "", format(round(Excluded / sum(Excluded), 2), nsmall = 2)),  # remove leading zero and keep 2 decimals
         across(starts_with("condition"), ~str_to_title(.))) %>% 
  rename("Discourse Condition" = condition_discourse,
         "Framing Condition" = condition_framing)

func_sup_tab <- function(df){
  flextable(df) %>% 
    autofit() %>%                                 # shrink to fit content
    font(fontname = "Times New Roman", part = "all") %>%  # font
    fontsize(size = 12, part = "all") %>%        # font size
    set_table_properties(width = 1, layout = "autofit") %>% # scale to page width
    border_remove() %>%                           # remove all borders first
    hline(i = 1, border = fp_border_default(), part = "header") %>%  # line below header
    align(align = "center", part = "all")        # center all cells
  }

ft_sup_tab1 <- func_sup_tab(sup_tab1)
# save_as_docx(ft_sup_tab1, path = "tables/sup_tab1.docx")

### 2 ####
sup_tab2 <- wepd_data2 %>% 
  distinct(participant_id, .keep_all = T) %>% 
  group_by(location, condition_discourse, condition_framing) %>%
  summarize(n = n(), 
            .groups = "drop") %>% 
  mutate("Testing Location" = case_when(location == 0 ~ "Uni. Lab",
                                        location == 1 ~ "Museum Demo.",
                                        location == 2 ~ "Museum Lab"),
         location = as.character(location),
         condition_discourse = ifelse(condition_discourse == "conventional", "conv.", "idio.")) %>% 
  group_by(location) %>% 
  mutate(included = sum(n),
         across(starts_with("condition"), ~str_to_title(.))) %>%
  pivot_wider(names_from = c(condition_discourse, condition_framing), values_from = n) %>% 
  rename_with(~str_replace_all(.x, "_", "-"), matches("_")) %>% 
  relocate(included, .after = everything()) %>%
  left_join(exclusions %>%
              group_by(location) %>%
              summarize(Excluded = n(),
                        .groups = "drop") %>% 
              mutate(location = as.character(location)),
            by = "location") %>% 
  ungroup() %>% 
  mutate(ProportionExcluded = sub("^0", "", format(round(Excluded/included, 2), nsmall = 2))) %>% 
  rename("Prop. Excluded" = ProportionExcluded) %>% 
  rename("Total Included" = included) %>% 
  select(-location)

ft_sup_tab2 <- func_sup_tab(sup_tab2)
# save_as_docx(ft_sup_tab2, path = "tables/sup_tab2.docx")

### 3 ####
sup_tab3 <- exclusions %>% 
  group_by(exclude_reason) %>% 
  summarize(Excluded = n(),
            .groups = "drop") %>% 
  mutate(Proportion = round(Excluded / sum(Excluded), 2),
         Proportion = sub("^0", "", format(Proportion, nsmall = 2)),
         "Exclusion Reason" = case_when(exclude_reason == "exp_error" ~ "Experimenter Error",
                                    exclude_reason == "idiosyncratic" ~ "Idiosyncratic",
                                    exclude_reason == "nonenglish" ~ "Non-English",
                                    exclude_reason == "shy" ~ "Shy",
                                    exclude_reason == "wrong_age" ~ "Wrong Age")) %>% 
  select(-exclude_reason) %>% 
  relocate("Exclusion Reason", .before = everything())

ft_sup_tab3 <- func_sup_tab(sup_tab3)
# save_as_docx(ft_sup_tab3, path = "tables/sup_tab3.docx")

### 4 ####
sup_tab4 <- wepd_data2 %>% 
  distinct(participant_id, .keep_all = T) %>% 
  group_by(income_recoded) %>% 
  count(name = "income_n") %>% 
  ungroup() %>% 
  mutate(income_prop = sub("^0", "", format(round(income_n / sum(income_n), 2), nsmall = 2))) %>% 
  bind_cols(wepd_data2 %>% 
              distinct(participant_id, .keep_all = T) %>% 
              group_by(race_recoded) %>% 
              count(name = "race_n") %>% 
              ungroup() %>% 
              mutate(race_prop = sub("^0", "", format(round(race_n / sum(race_n), 2), nsmall = 2)))) %>% 
  rename(Income = income_recoded, "Income n" = income_n, "Income Proportion" = income_prop,
         Race = race_recoded, "Race n" = race_n, "Race Proportion" = race_prop)

ft_sup_tab4 <- func_sup_tab(sup_tab4)
# save_as_docx(ft_sup_tab4, path = "tables/sup_tab4.docx")

### 5 ####
sup_tab5 <- wepd_data2 %>% 
  distinct(participant_id, .keep_all = T) %>% 
  group_by(e1_recoded, condition_discourse, condition_framing) %>%
  summarize(n = n(), 
            .groups = "drop") %>% 
  mutate(E1 = case_when(e1_recoded == "alissa"  ~ "1",
                        e1_recoded == "jadelyn" ~ "2",
                        e1_recoded == "kaelin"  ~ "3",
                        e1_recoded == "maya&naomi" ~ "4*"),
         e1_recoded = as.character(e1_recoded),
         condition_discourse = ifelse(condition_discourse == "conventional", "conv.", "idio.")) %>% 
  group_by(e1_recoded) %>% 
  mutate(included = sum(n),
         across(starts_with("condition"), ~str_to_title(.))) %>%
  pivot_wider(names_from = c(condition_discourse, condition_framing), values_from = n) %>% 
  relocate(included, .after = everything()) %>%
  ungroup() %>% 
  left_join(exclusions %>%
              group_by(e1_recoded) %>%
              summarize(Excluded = n(),
                        .groups = "drop") %>% 
              mutate(e1_recoded = as.character(e1_recoded)),
            by = "e1_recoded") %>% 
  ungroup() %>% 
  mutate(ProportionExcluded = sub("^0", "", format(round(Excluded/included, 2), nsmall = 2)),
         "Proportion of Inclusions" = sub("^0", "", format(round(included/sum(included), 2), nsmall = 2))) %>% 
  rename("Prop. Excluded" = ProportionExcluded) %>% 
  rename("Total Included" = included) %>% 
  select(-e1_recoded) %>% 
  rename_with(~str_replace_all(.x, "_", "-"), matches("_")) %>% 
  relocate("Proportion of Inclusions", .after = "Total Included")

ft_sup_tab5 <- func_sup_tab(sup_tab5)
# save_as_docx(ft_sup_tab5, path = "tables/sup_tab5.docx")

### 6 ####
sup_tab6 <- as_draws_df(mod) %>% # model posterior fixed effects
  select(starts_with("b_")) %>% 
  pivot_longer(everything(), names_to = "parameter") %>%
  group_by(parameter) %>% 
  summarize(mean = mean(value),
            median = median(value),
            sd = sd(value),
            mad = mad(value),
            hdi_lo = hdi(value)[ncol(hdi(value)) - 1], # rightmost HDI LB to safeguard against multimodality
            hdi_hi = hdi(value)[ncol(hdi(value))], # rightmost HDI UB to safeguard against multimodality
            q2.5   = quantile(value, 0.025),
            q97.5   = quantile(value, 0.975),
            prop_gt_0 = mean(value > 0)) %>% 
  mutate(across(where(is.numeric), ~ formatC(.x, format = "f", digits = 2)),
         prop_gt_0 = sub("^0", "", prop_gt_0),
         parameter = str_to_title(str_replace_all(str_replace_all(str_remove(parameter, "b_"), "_", " "), ":", " : ")),
         parameter = case_when(parameter == "Condition Discourseidiosyncratic" ~ "Idiosyncratic Discourse",
                               str_detect(parameter, ":") ~ "Discourse : Framing",
                               parameter == "Condition Framingyou" ~ "You-Framing",
                               parameter == "Sexm" ~ "Boys",
                               T ~ parameter)) %>% 
  rename_with(~ str_replace_all(.x, "_", " "), everything()) %>% 
  rename_with(str_to_title, everything()) %>% 
  rename_with(toupper, c(Sd, Mad)) %>% 
  rename_with(~ str_replace(.x, "Hdi", "HDI"), everything()) %>% 
  slice(c(1:5, 7, 6, 8:n()))

ft <- func_sup_tab(sup_tab6) %>%  width(j = "Q2.5", width = 1)  # width in inches
# save_as_docx(ft, path = "tables/sup_tab6.docx")

## figures ####
### 1 ####
# marginal predicted informativeness
epreds_expected_value_discourse <- epreds_expected_value %>% # discourse draws df
  group_by(.draw, condition_discourse) %>% # marginalize over covariates
  summarize(expected_value = mean(expected_value), .groups = "drop") # expected informativeness per condition combination per draw
epreds_expected_value_framing <- epreds_expected_value %>% # framing draws df
  group_by(.draw, condition_framing) %>% # marginalize over covariates
  summarize(expected_value = mean(expected_value), .groups = "drop") # expected informativeness per condition combination per draw

draws_df_marginal <- bind_rows(epreds_expected_value_discourse, epreds_expected_value_framing) %>%
  mutate(condition = case_when(condition_discourse == "conventional" ~ "Conventional\nDiscourse",
                               condition_discourse == "idiosyncratic" ~ "Idiosyncratic\nDiscourse",
                               condition_framing == "we" ~ "\n''We''-Framing",
                               condition_framing == "you" ~ "\n''You''-Framing")) %>% 
  select(-starts_with("condition_"))

draws_df_marginal_summary <- draws_df_marginal %>% 
  mutate(grand_mean = mean(expected_value)) %>% 
  group_by(condition) %>%
  summarize(hdi_lo = hdi(expected_value)[ncol(hdi(expected_value)) - 1], # rightmost HDI LB to safeguard against multimodality
            mean_expected_value = mean(expected_value),
            hdi_hi = hdi(expected_value)[ncol(hdi(expected_value))], # rightmost HDI UB to safeguard against multimodality
            prop_gt_gm = mean(expected_value > grand_mean),
            prop_gt_gm_char = sprintf("%.0f%%", 100 * prop_gt_gm),
            .groups = "drop") %>% 
  mutate(grand_mean = mean(draws_df_marginal$expected_value))
  
ggplot(draws_df_marginal,
       aes(y = condition, x = expected_value)) +
  stat_halfeye(aes(fill = after_stat(x < draws_df_marginal_summary$grand_mean)), .width = 0.95, slab_alpha = 0.6, slab_color = "black", slab_linewidth = 0.5) +
  scale_fill_manual(values = c("TRUE" = "grey70", "FALSE" = "grey30"), guide = "none") + 
  theme_classic(base_family = "Times New Roman", base_size = 15) +
  labs(x = "Predicted Informativeness", title = "Predicted Informativeness", y = "Factor Level",
       subtitle = "How informative are participants, on average, in each factor level?") +
  geom_vline(xintercept = draws_df_marginal_summary$grand_mean, linetype = 2) +
  scale_x_continuous(labels = function(x) format(x, trim = TRUE),
                     breaks = pretty_breaks(4)) +
  scale_y_discrete(expand = expansion(mult = c(0.1, 0.02))) +
  geom_text(data = draws_df_marginal_summary,
            aes(x = hdi_hi + .02, label = prop_gt_gm_char),
            size = 4.5,
            vjust = -2,
            family = "Times New Roman",
            fontface = "bold",
            color = "grey30") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(face = "italic", hjust = .8, size = 12),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(vjust =-.5, hjust = .5),
        legend.position = "none")
# ggsave("figs/sup_fig1.png", width = 6, height = 4)

### 2 ####
# conditional predicted informativeness
epreds_expected_value_discourse_framing <- epreds_expected_value %>% # framing * discourse draws df
  group_by(.draw, condition_discourse, condition_framing) %>% # marginalize over covariates
  summarize(expected_value = mean(expected_value), .groups = "drop") %>% # expected informativeness per condition combination per draw
  mutate(condition_discourse = fct_rev(case_when(condition_discourse == "conventional" ~ "Conventional\nDiscourse",
                                         condition_discourse == "idiosyncratic" ~ "Idiosyncratic\nDiscourse")),
         condition_framing = case_when(condition_framing == "we" ~ '"We"-Framing',
                                       condition_framing == "you" ~ '"You"-Framing'))
         
ggplot(epreds_expected_value_discourse_framing,
       aes(y = condition_discourse, x = expected_value)) +
  stat_halfeye(fill = "forestgreen", .width = 0.95, slab_alpha = 0.6, slab_color = "black", slab_linewidth = 0.8) +
  facet_wrap(~ condition_framing) +
  theme_classic(base_family = "Times New Roman", base_size = 15) +
  labs(x = "Predicted Informativeness", title = "Conditional Predicted Informativeness", y = "Discourse Level",
       subtitle = "How informative are participants in each condition?") +
  scale_x_continuous(labels = function(x) format(x, trim = TRUE),
                     breaks = pretty_breaks(3)) +
  scale_y_discrete(expand = expansion(mult = c(0.1, 0.02))) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(face = "italic", hjust = .5, size = 12),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(vjust = -0.5, hjust = .5),
        legend.position = "none")
# ggsave("figs/sup_fig2.png", width = 7, height = 3)

 ### 3 ####
# framing contrasts
contrasts_framing_marginal <- epreds_expected_value %>% # draws of marginal comparison of discourse
  group_by(.draw, condition_framing) %>% 
  summarize(mean_expected_value = mean(expected_value),
            .groups = "drop") %>% 
  mutate(condition_framing = factor(condition_framing, levels = c("you", "we"))) %>% 
  group_by(condition_framing) %>%
  compare_levels(mean_expected_value, by = condition_framing) %>%
  ungroup() %>% 
  rename(comparison = condition_framing,
         diff_mean_expected_value = mean_expected_value)
contrasts_framing_marginal_summary <- contrasts_framing_marginal %>% # summary of marginal comparison
  group_by(comparison) %>% 
  summarize(hdi_lo = hdi(diff_mean_expected_value)[ncol(hdi(diff_mean_expected_value)) - 1], # rightmost HDI LB to safeguard against multimodality
            mean_expected_value = mean(diff_mean_expected_value),
            hdi_hi = hdi(diff_mean_expected_value)[ncol(hdi(diff_mean_expected_value))], # rightmost HDI UB to safeguard against multimodality
            prop_gt = mean(diff_mean_expected_value > 0),
            prop_gt_char = sprintf("%.0f%%", 100 * prop_gt),
            .groups = "drop")

contrasts_framing_conditional <- epreds_expected_value %>% # draws of marginal comparison of discourse conditional on framing
  group_by(.draw, condition_framing, condition_discourse) %>% 
  summarize(mean_expected_value = mean(expected_value),
            .groups = "drop") %>% 
  mutate(condition_framing = factor(condition_framing, levels = c("you", "we"))) %>% 
  group_by(condition_discourse) %>%
  compare_levels(mean_expected_value, by = condition_framing) %>%
  ungroup() %>% 
  rename(comparison = condition_framing ,
         diff_mean_expected_value = mean_expected_value)
contrasts_framing_conditional_summary <- contrasts_framing_conditional %>% # summary of conditional marginal comparison
  group_by(condition_discourse, comparison) %>% 
  summarize(hdi_lo = hdi(diff_mean_expected_value)[ncol(hdi(diff_mean_expected_value)) - 1], # rightmost HDI LB to safeguard against multimodality
            mean_expected_value = mean(diff_mean_expected_value),
            hdi_hi = hdi(diff_mean_expected_value)[ncol(hdi(diff_mean_expected_value))], # rightmost HDI UB to safeguard against multimodality
            prop_gt = mean(diff_mean_expected_value > 0),
            prop_gt_char = sprintf("%.0f%%", 100 * prop_gt),
            .groups = "drop")

contrasts_framing_draws <- bind_rows(contrasts_framing_marginal %>% mutate(condition_discourse = "overall"), contrasts_framing_conditional) %>% 
  mutate(across(starts_with("condition"), ~str_to_title(.)),
         condition_discourse = ifelse(condition_discourse == "Overall",
                                      "Overall",
                                      paste0(condition_discourse, " Discourse")),
         condition_discourse = fct_recode(condition_discourse,
                                          "Idiosyncratic\nDiscourse" = "Idiosyncratic Discourse",
                                          "Conventional\nDiscourse" = "Conventional Discourse",
                                          "Overall\n" = "Overall"),
         condition_discourse = fct_relevel(condition_discourse,
                                           "Idiosyncratic\nDiscourse",
                                           "Conventional\nDiscourse",
                                           "Overall\n"))
contrasts_framing_summary <- bind_rows(contrasts_framing_marginal_summary %>% mutate(condition_discourse = "overall"), contrasts_framing_conditional_summary) %>% 
  mutate(hjust_label = ifelse(hdi_hi > -0.010 & hdi_hi < 0.010, 0.2, NA),
         across(starts_with("condition"), ~str_to_title(.)),
         condition_discourse = ifelse(condition_discourse == "Overall",
                                      "Overall",
                                      paste0(condition_discourse, " Discourse")),
         condition_discourse = fct_recode(condition_discourse,
                                          "Idiosyncratic\nDiscourse" = "Idiosyncratic Discourse",
                                          "Conventional\nDiscourse" = "Conventional Discourse",
                                          "Overall\n" = "Overall"),
         condition_discourse = fct_relevel(condition_discourse,
                                           "Idiosyncratic\nDiscourse",
                                           "Conventional\nDiscourse",
                                           "Overall\n"))

ggplot(contrasts_framing_draws,
       aes(x = diff_mean_expected_value, y = condition_discourse)) +
  stat_halfeye(aes(fill = after_stat(x < 0)), .width = 0.95, slab_alpha = 0.6, slab_color = "black", slab_linewidth = 0.5) +
  scale_fill_manual(values = c("TRUE" = "grey70", "FALSE" = "grey30"), guide = "none") + 
  geom_text(data = contrasts_framing_summary,
            aes(x = hdi_hi + .03, label = prop_gt_char, hjust = hjust_label),
            size = 4,
            vjust = -3,
            family = "Times New Roman",
            fontface = "bold",
            color = "grey40") +
  geom_vline(xintercept = 0, alpha = 0.3, linetype = 2) +
  labs(x = "Difference in Predicted Informativeness",
       title = "Framing Condition Contrasts",
       subtitle = "Are participants more informative after ''we''-framing or ''you''-framing?",
       y = NULL) +
  scale_x_continuous(labels = function(x) format(x, trim = TRUE),
                     breaks = pretty_breaks(3)) +
  scale_y_discrete(expand = expansion(mult = c(0.05, .75))) +
  theme_classic(base_size = 15, base_family = "Times New Roman") +
  theme(axis.text.y = element_text(vjust = -.75, hjust = .50),
        axis.ticks.y = element_blank(),
        plot.title = element_text(face = "bold", hjust = .5, size = 16),
        plot.subtitle = element_text(face = "italic", hjust = .8, size = 13)) +
  annotate("text", x = -.30, 
           y =  4.2, 
           label = "''You''-Framing", hjust = 0.5, color = "grey70",
           family = "Times New Roman", fontface = "bold", size = 4.5) +
  annotate("text", x = .50, 
           y = 4.2, 
           label = "''We''-Framing", hjust = .5, color = "grey50",
           family = "Times New Roman", fontface = "bold", size = 4.5)
# ggsave("figs/sup_fig3.png", width = 6.5, height = 4)

### 4 ####
# expected probability of forms by framing condition
# step 1: framing * discourse
draws_mean_prob_discourse_framing <- epreds %>%
  mutate(.category_numeric = as.numeric(.category)) %>% # ranges 1-4
  group_by(.draw, condition_discourse, condition_framing, .category) %>% # marginalize across control vars
  summarise(mean_probability = mean(.epred), 
            .groups = "drop")
diff_draws_mean_prob_discourse_framing <- draws_mean_prob_discourse_framing %>% 
  mutate(condition_framing = factor(condition_framing, levels = c("you", "we"))) %>% 
  group_by(condition_discourse, .category) %>% 
  compare_levels(mean_probability, by = condition_framing) %>% 
  ungroup() %>% 
  rename(comparison = condition_framing,
         diff_mean_prob = mean_probability)
summary_diff_draws_mean_prob_discourse_framing <- diff_draws_mean_prob_discourse_framing %>% 
  group_by(comparison, condition_discourse, .category) %>% 
  summarize(hdi_lo = hdi(diff_mean_prob)[ncol(hdi(diff_mean_prob)) - 1], # rightmost HDI LB to safeguard against multimodality
            mean_diff_mean_prob = mean(diff_mean_prob),
            hdi_hi = hdi(diff_mean_prob)[ncol(hdi(diff_mean_prob))], # rightmost HDI UB to safeguard against multimodality
            prop_gt = mean(diff_mean_prob > 0),
            prop_gt_char = sprintf("%.0f%%", 100 * prop_gt),
            .groups = "drop")

# step 2: framing
draws_mean_prob_framing <- epreds %>%
  mutate(.category_numeric = as.numeric(.category)) %>% # ranges 1-4
  group_by(.draw, condition_framing, .category) %>% # marginalize across control vars
  summarise(mean_probability = mean(.epred), 
            .groups = "drop") 
diff_draws_mean_prob_framing <- draws_mean_prob_framing %>% 
  mutate(condition_framing = factor(condition_framing, levels = c("you", "we"))) %>% 
  group_by(condition_framing, .category) %>% 
  compare_levels(mean_probability, by = condition_framing) %>% 
  ungroup() %>% 
  rename(comparison = condition_framing,
         diff_mean_prob = mean_probability)
summary_diff_draws_mean_prob_framing <- diff_draws_mean_prob_framing %>% 
  group_by(comparison, .category) %>% 
  summarize(hdi_lo = hdi(diff_mean_prob)[ncol(hdi(diff_mean_prob)) - 1], # rightmost HDI LB to safeguard against multimodality
            mean_diff_mean_prob = mean(diff_mean_prob),
            hdi_hi = hdi(diff_mean_prob)[ncol(hdi(diff_mean_prob))], # rightmost HDI UB to safeguard against multimodality
            prop_gt = mean(diff_mean_prob > 0),
            prop_gt_char = sprintf("%.0f%%", 100 * prop_gt),
            .groups = "drop")

# step 3: combine df's
contrasts_cat_f <- bind_rows(diff_draws_mean_prob_framing %>% mutate(condition_discourse = "overall"), diff_draws_mean_prob_discourse_framing) %>% 
  mutate(across(starts_with("condition"), ~str_to_title(.)),
         .category = case_when(.category == "nonling_only" ~ "Nonlinguistic",
                               .category == "proform" ~ "Proform",
                               .category == "np_unmod" ~ "Unmod. NP",
                               .category == "np_mod" ~ "Mod. NP"),
         .category = fct_relevel(.category, "Mod. NP", "Unmod. NP", "Proform", "Nonlinguistic"),
         condition_discourse = fct_relevel(condition_discourse, "Overall")) # overall column first
contrasts_cat_summary_f <- bind_rows(summary_diff_draws_mean_prob_framing %>% mutate(condition_discourse = "overall"), summary_diff_draws_mean_prob_discourse_framing) %>% 
  mutate(hjust_label = ifelse(hdi_hi > -0.010 & hdi_hi < 0.010, 0.2, NA),
         across(starts_with("condition"), ~str_to_title(.)),
         .category = case_when(.category == "nonling_only" ~ "Nonlinguistic",
                               .category == "proform" ~ "Proform",
                               .category == "np_unmod" ~ "Unmod. NP",
                               .category == "np_mod" ~ "Mod. NP"),
         .category = fct_relevel(.category, "Mod. NP", "Unmod. NP", "Proform", "Nonlinguistic"),
         condition_discourse = fct_relevel(condition_discourse, "Overall")) # overall column first

# step 4: plot
ggplot() + # first plot: left y-axis label
  geom_richtext(aes(x = 0, y = 0, label = "<span style='color:grey70'>\"You\"-Framing</span><br><br>
                 <span style='color:black'>vs.</span><br><br>
                 <span style='color:grey50'>\"We\"-Framing</span>"),
                fill = NA, label.color = NA,
                size = 5,
                hjust = 0.5, vjust = 0.5, family = "Times New Roman", fontface = "bold") + 
  theme_void() +
  ggplot(contrasts_cat_f, aes(x = diff_mean_prob, y = comparison)) + # second plot: densities
  stat_halfeye(aes(fill = after_stat(x < 0)), .width = 0.95, slab_color = "black",
               slab_alpha = 0.6, slab_linewidth = 0.5) +
  scale_fill_manual(values = c("TRUE" = "grey70", "FALSE" = "grey30"), guide = "none") + 
  facet_grid(.category ~ condition_discourse, scales = "free") +
  geom_text(data = contrasts_cat_summary_f,
            aes(x = ifelse(hdi_hi < .2, hdi_hi + .05, hdi_hi + .06), label = prop_gt_char, hjust = hjust_label),
            vjust = -2.25,
            size = 4.5,
            family = "Times New Roman",
            fontface = "bold",
            color = "grey40") +
  geom_vline(xintercept = 0, alpha = 0.3, linetype = 2) +
  theme_classic(base_size = 15, base_family = "Times New Roman") +
  labs(y = NULL, x = "Difference in Predicted Probability of Use",
       subtitle = "Which referring forms are used more often after ''we''-framing?") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(face = "italic", hjust = .5, size = 14)) +
  scale_x_continuous(breaks = seq(-0.2, 0.2, by = 0.2)) +
  scale_y_discrete(expand = expansion(mult = c(0.1, 0))) +
  plot_annotation(
    title = "Framing Contrasts by Referring Form",
    theme = theme(plot.title = element_text(face = "bold", family = "Times New Roman", hjust = .6, size = 17))) +
  ggplot() + # third plot: right y-axis label
  geom_richtext(aes(x = 0, y = 0, label = "Referring Form"),
                fill = NA, label.color = NA, 
                family = "Times New Roman",
                fontface = "bold",
                size = 5, angle = 270) +
  theme_void() +
  plot_layout(widths = c(1.5, 6, 0.4))
# ggsave("figs/sup_fig4.png", width = 8, height = 6)

### 5 ####
# trial contrasts
contrasts_trial_marginal <- epreds_expected_value %>% # draws of marginal comparison of trial
  group_by(.draw, trial) %>% 
  summarize(mean_expected_value = mean(expected_value),
            .groups = "drop") %>% 
  mutate(trial = factor(trial, levels = c("1", "2"))) %>% 
  group_by(trial) %>%
  compare_levels(mean_expected_value, by = trial) %>%
  ungroup() %>% 
  rename(comparison = trial,
         diff_mean_expected_value = mean_expected_value)

contrasts_trial_marginal_summary <- contrasts_trial_marginal %>% # summary of marginal comparison
  group_by(comparison) %>% 
  summarize(hdi_lo = hdi(diff_mean_expected_value)[ncol(hdi(diff_mean_expected_value)) - 1], # rightmost HDI LB to safeguard against multimodality
            mean_expected_value = mean(diff_mean_expected_value),
            hdi_hi = hdi(diff_mean_expected_value)[ncol(hdi(diff_mean_expected_value))], # rightmost HDI UB to safeguard against multimodality
            prop_gt = mean(diff_mean_expected_value > 0),
            prop_gt_char = sprintf("%.0f%%", 100 * prop_gt),
            .groups = "drop")

ggplot(contrasts_trial_marginal,
       aes(x = diff_mean_expected_value, y = comparison)) +
  stat_halfeye(aes(fill = after_stat(x < 0)), .width = 0.95, slab_alpha = 0.6, slab_linewidth = 0.5, slab_color = "black") +
  scale_fill_manual(values = c("TRUE" = "grey70", "FALSE" = "grey30"), guide = "none") + 
  geom_text(data = contrasts_trial_marginal_summary,
            aes(x = hdi_hi + .02, label = prop_gt_char),
            size = 5,
            vjust = -5,
            family = "Times New Roman",
            fontface = "bold",
            color = "grey40") +
  geom_vline(xintercept = 0, alpha = 0.3, linetype = 2) +
  labs(x = "Difference in Predicted Informativeness",
       title = "Trial Contrast",
       subtitle = "Are participants more informative in trial 1 or trial 2?",
       y = NULL) +
  scale_x_continuous(labels = function(x) format(x, trim = TRUE),
                     breaks = pretty_breaks(3)) +
  scale_y_discrete(expand = expansion(mult = c(0.05, 1.05))) +
  theme_classic(base_size = 16, base_family = "Times New Roman") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        plot.title = element_text(face = "bold", hjust = .5, size = 17),
        plot.subtitle = element_text(face = "italic", hjust = .5, size = 13)) +
  annotate("text", x = -.15, 
           y =  1.9, 
           label = "Trial 1", color = "grey70",
           family = "Times New Roman", fontface = "bold", size = 5) +
  annotate("text", x = .15, 
           y = 1.9, 
           label = "Trial 2", color = "grey50",
           family = "Times New Roman", fontface = "bold", size = 5)
# ggsave("figs/sup_fig5.png", width = 5, height = 2.5)

# reliability ####
reliability_ch <- read_csv("data/reliability_wepd.csv")
reliability_jv <- wepd_data2 %>%
  filter(participant_id %in% reliability_ch$participant_id) %>%
  rename_with(~paste0(., "_jv"), everything())

func_reliability <- function(var){
  var_jv <- paste0(var, "_jv")
  reliability <- tibble(reliability_ch %>% select({{var}}),
                        reliability_jv %>% select(all_of(var_jv))) %>% 
    mutate(across(everything(), ~ifelse(is.na(.), "na", paste0(.))))
  kappa <- kappa2(reliability, weight = "squared") # squared-weighted kappa for ordinal (larger penalties for larger disagreements)
  agree <- agree(reliability)
  return(list(reliability, kappa, agree))
}

reliability_item <- func_reliability("ref_item")
reliability_form <- func_reliability("form")

end <- Sys.time()
end - start
