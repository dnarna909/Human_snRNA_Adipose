# https://kstreet13.github.io/bioc2020trajectories/articles/workshopTrajectories.html
# Differential progression ------
load( paste0(Disk, Project.folder, "/", Rds.folder, "/","Differentiation_PHATE_", groups, ".Peak.df.Rds")) 
Peak.df; Psu.df
# sample.meta <- 
  
# Peak.df <- Peak.df %>% 
#   tibble::rownames_to_column(var = "Dataset") %>% 
#   left_join(sample.meta, by = join_by(Dataset == Dataset) ) %>%
#   distinct( Dataset, .keep_all = TRUE)
# Psu.df <- Psu.df %>% 
#   left_join(sample.meta, by = join_by(Dataset.Group == Dataset) ) 

Compare.Group = "BA_Group"
# groups = levels(Psu.df$Group)[levels(Psu.df$Group) %in% unique(Psu.df$Group)];groups
groups = unique(Psu.df[[Compare.Group]]);groups

p <- ggplot(Psu.df, aes(x=Pt, group=.data[[paste0("Dataset", ".Group")]], color= .data[[Compare.Group]])) +
  geom_density() + # Use semi-transparent fill
  facet_grid(formula(paste(". ~",  "Compare.Group")))+ 
  scale_fill_brewer(palette="Dark2") + theme_minimal()
print(p)

p <- ggplot(Psu.df, aes(x=Pt, group=.data[[Compare.Group]], color= .data[[Compare.Group]])) +
  geom_density() + # Use semi-transparent fill
  facet_grid(formula(paste(". ~",  "Compare.Group")))+ 
  scale_fill_brewer(palette="Dark2") + theme_minimal()
print(p)

# multiple sample Pt density test
library("kSamples")
set.seed(142)
Result <- ad.test(as.vector(Psu.df$Pt) ~ as.vector(Psu.df[[Compare.Group]]) ) 
Result.ad <- mean(as.data.frame(Result[["ad"]])$` asympt. P-value`); Result.ad

#  interpret this question as a univariate analysis of the pseudotime values between the two groups.
########################
### Kolmogorov-Smirnov Test
########################

ks.test(Psu.ls[[1]]$Pt,
        Psu.ls[[2]]$Pt)


# test for peak.df
library(car)
library(rstatix)
Peak.df[[Compare.Group]] <- factor(Peak.df[[Compare.Group]], levels = levels(Peak.df[[Compare.Group]])[levels(Peak.df[[Compare.Group]]) %in% unique(Peak.df[[Compare.Group]])] )
tes
t.result <- Peak.df %>%
  as.data.frame() %>% # nest() %>% pluck("data", 3) %>% 
  rstatix::t_test(formula(paste0("Peak2.y", " ~ ", "Compare.Group"))) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj") %>% 
  mutate(method = "T-test") %>% 
  add_xy_position(x = "Compare.Group")
ano.result <- Peak.df  %>%
  as.data.frame() %>%
  rstatix::anova_test(formula(paste0("Peak2.y", " ~ ", "Compare.Group")), type = 3) %>% 
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj") %>% 
  ungroup() 