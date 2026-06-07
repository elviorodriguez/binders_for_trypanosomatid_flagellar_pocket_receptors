
library(tidyverse)

designs_trf <- read.delim(file = "./design_analysis_results_tfr.csv", sep = ",") %>% 
  mutate(Target = "TfR")
designs_isg65 <- read.delim(file = "./design_analysis_results_isg65.csv", sep = ",") %>% 
  mutate(Target = "ISG65")
designs_hphbr1 <- read.delim(file = "./design_analysis_results_hphbr1.csv", sep = ",") %>% 
  mutate(Target = "HpHbR1")
designs_hphbr2 <- read.delim(file = "./design_analysis_results_hphbr2.csv", sep = ",") %>% 
  mutate(Target = "HpHbR2")

designs <- bind_rows(designs_trf, designs_isg65, designs_hphbr1, designs_hphbr2) %>%
  mutate(Target = factor(Target,
                         levels = c("TfR", "ISG65", "HpHbR1", "HpHbR2")))

filtered_designs <- designs %>% 
  filter(rank == 1,
         min_interaction_pae < 5,
         rmsd_to_rfdiff < 2)

# Scatter plot filtering by PAE<5 and RMSD<2
filtered_designs %>% 
  ggplot() +
  aes(x = min_interaction_pae, y = rmsd_to_rfdiff, fill = Target) +
  geom_point(shape = 21,size = 4,alpha = 0.8,color = "black",stroke = 0.5) +
  scale_fill_manual(
    values = c("TfR"="pink","ISG65"="darkorange","HpHbR1"="forestgreen","HpHbR2"="green")) +
  xlab("Minimum interaction PAE (Å)") +
  ylab("RMSD to RFdiffusion (Å)") +
  theme_bw()+ 
  theme(text = element_text(size = 16.5)) 

# No filter scatter plots
designs %>% 
  ggplot() +
  aes(x = min_interaction_pae, y = rmsd_to_rfdiff, fill = Target) +
  geom_point(shape = 21, size = 2, alpha = 0.5) +
  geom_vline(xintercept = 5, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 2, linetype = "dashed", color = "black") +
  scale_fill_manual(
    values = c("TfR"="pink","ISG65"="darkorange","HpHbR1"="forestgreen","HpHbR2"="green")) +
  scale_x_continuous(
    breaks = seq(0, 30, 5),
    labels = c("0", "5", "10", "15", "20", "25", "")
  ) +
  xlab("Minimum interaction PAE (Å)") +
  ylab("RMSD to RFdiffusion (Å)") +
  theme_bw()+ 
  theme(text = element_text(size = 16.5)) +
  facet_wrap(~Target) +
  facet_wrap(~Target, nrow = 1)

# Nº of successful models for each target
nrow(filtered_designs %>% filter(Target == "TfR"))
nrow(filtered_designs %>% filter(Target == "ISG65"))
nrow(filtered_designs %>% filter(Target == "HpHbR1"))
nrow(filtered_designs %>% filter(Target == "HpHbR2"))

# Distribution of raw score
designs %>% filter(raw_score > 0.5) %>% 
  ggplot() +
  aes(x = raw_score, fill = Target) +
  geom_histogram(binwidth = 0.1) +
  scale_fill_manual(
    values = c("TfR"="pink","ISG65"="darkorange","HpHbR1"="forestgreen","HpHbR2"="green")) +
  ylab("Frquency") +
  xlab("Score (S)") +
  theme_bw() +
  theme(text = element_text(size = 16.5))


# > names(designs)
# [1] "design_id"                  "target"                     "rank"                       "mean_binder_plddt"         
# [5] "mean_binder_intrachain_pae" "min_interaction_pae"        "rmsd_to_rfdiff"             "binder_length"             
# [9] "target_length"              "model_file"                 "ptm"                        "iptm"                      
# [13] "pae_rank"                   "raw_score"                 