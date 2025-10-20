# Sources ----

# Open all libraries

library(mice)
library(readr)
library(stringr)
library(nlme)
library(lmerTest)
library(dplyr)
library(lattice)
library(FactoMineR)
library(MuMIn)
library(vegan)
library(ggplot2)
library(ggpubr)
library(factoextra)
library(effectsize)


# Seed = 195, generated at random.org (min = 1, max = 1000); 
# date/hour = 2023-01-17 09:49:56 UTC


# Data loading ----

# Load base tibble

microcosm_master <- read_csv("data/microcosm_master.csv")


# Activate filter_1 = 1

microcosm_master <- 
  microcosm_master %>%
  filter(
    filter_1 == 1 & pfg_combo != "gc3_legume_forb"
  ) %>%
  mutate(
    sp_combo_dom = str_c(
      sp_combo,
      dom_status,
      sep = "_"
    )
  )


# Create sp_pfg_richness variable

microcosm_master$pfg_number <- as.character(microcosm_master$pfg_number)

microcosm_master$sp_richness <- as.character(microcosm_master$sp_richness)

microcosm_master <-
  microcosm_master %>%
  mutate(
    sp_pfg_richness = str_c(
      sp_richness,
      pfg_number,
      sep = "_"
    )
  )


# Parse columns to factors as appropriate

microcosm_master$sp <- parse_factor(microcosm_master$sp_combo)

microcosm_master$dom_status <- parse_factor(microcosm_master$dom_status)

microcosm_master$organ <- parse_factor(microcosm_master$organ)

microcosm_master$pfg_combo <- parse_factor(microcosm_master$pfg_combo)

microcosm_master$pfg_number <- as.factor(microcosm_master$pfg_number)

microcosm_master$sp_richness <- as.factor(microcosm_master$sp_richness)

microcosm_master$sp_pfg_richness <- as.factor(microcosm_master$sp_pfg_richness)

microcosm_master$sp_pfg_richness <- ordered(microcosm_master$sp_pfg_richness,
                                            levels = c("1_1", "2_1", "2_2", "4_4"))



# H1: Litter quality PCA ----


# Open traits data

i_traits <- read.table(
  "data/i_traits.txt", 
  header = T, sep = "\t", 
  colClasses = c("integer", 
                 rep("factor", 4), 
                 rep("numeric", 10), 
                 "character")
)

i_traits <- as_tibble(i_traits)


# Create DMC and remove useless variables

i_quality <- i_traits %>%
  mutate(
    dmc = dry_mass / humid_mass
  ) %>%
  select(
    -dry_mass, 
    -humid_mass,
    -notes
  )


# Account NA to decide imputation method

ncol(i_quality[, 6:14]) * 96 # Total cells in tibble

sum(is.na(i_quality[, 6:14]))

(68 / 864) * 100 # 7.87% missing data, so we can impute with PMM



## Imputation: PMM ----

set.seed(195)

imputations <- mice(
  i_quality, m = 100, method = "pmm",
  predictorMatrix = quickpred(i_quality, mincor = 0.1),
  include = c("sp", "dom_status", "organ"),
  exclude = c("id", "replicate"),
  seed = 195
)

View(imputations$imp)


# Everything good, so isolate variables in tibble

imp_data <- list(
  in. = as_tibble(t(imputations$imp$in.)),
  ic = as_tibble(t(imputations$imp$ic)),
  ca = as_tibble(t(imputations$imp$ca)),
  mg = as_tibble(t(imputations$imp$mg)),
  mn = as_tibble(t(imputations$imp$mn)),
  p = as_tibble(t(imputations$imp$p)),
  tannins = as_tibble(t(imputations$imp$tannins)),
  dmc = as_tibble(t(imputations$imp$dmc))
)


# Get all medians in a list of named vectors

imp_value <- list(
  in. = c(),
  ic = c(),
  ca = c(),
  mg = c(),
  mn = c(),
  p = c(),
  tannins = c(),
  dmc = c()
)

for (i in seq_along(imp_data)){
  values <- sapply(imp_data[[i]], FUN = median)
  imp_value[[i]] <- values
}

rm(i)


# Make a matrix with the vectors in the list

imp_tibble <- bind_rows(imp_value$in.,
                        imp_value$ic,
                        imp_value$ca,
                        imp_value$mg,
                        imp_value$mn,
                        imp_value$p,
                        imp_value$tannins,
                        imp_value$dmc)

id_imp <- as.integer(colnames(imp_tibble))

names_imp <- names(imp_value)

imp_tibble <- as_tibble(t(imp_tibble))

colnames(imp_tibble) <- names_imp

imp_tibble <-
  imp_tibble %>%
  mutate(
    icn = ic / in.,
    id = id_imp
  ) %>%
  relocate(
    id
  ) %>%
  relocate(
    icn, .after = ic
  )

rm(imp_data, imp_value, imputations)


# Complete original data with the imputations

quality_complete <- 
  left_join(
    i_quality, 
    imp_tibble, 
    by = "id"
  ) %>% 
  mutate(
    in. = ifelse(is.na(in..x), in..y, in..x),
    ic = ifelse(is.na(ic.x), ic.y, ic.x),
    icn = ifelse(is.na(icn.x), icn.y, icn.x),
    ca = ifelse(is.na(ca.x), ca.y, ca.x),
    mg = ifelse(is.na(mg.x), mg.y, mg.x),
    mn = ifelse(is.na(mn.x), mn.y, mn.x),
    p = ifelse(is.na(p.x), p.y, p.x),
    tannins = ifelse(is.na(tannins.x), tannins.y, tannins.x),
    dmc = ifelse(is.na(dmc.x), dmc.y, dmc.x)
  ) %>%
  select(
    -in..x, -in..y,
    -ic.x, -ic.y,
    -icn.x, -icn.y,
    -ca.x, -ca.y,
    -mg.x, -mg.y,
    -mn.x, -mn.y,
    -p.x, -p.y,
    -tannins.x, -tannins.y,
    -dmc.x, -dmc.y
  )


quality_complete$dom_status <- relevel(
  quality_complete$dom_status,
  ref = "w")


# Replace tannins and mn values < 0 with 0

quality_complete$mn[quality_complete$mn < 0] <- 0

quality_complete$tannins[quality_complete$tannins < 0] <- 0

quality_complete <- quality_complete[1:96, ]



## PCA: leaves vs. roots ----

pca_quality <- PCA(quality_complete[, 8:14], 
                   scale.unit = TRUE, 
                   graph = FALSE)

summary(pca_quality) # 30.25%, 22.32% 


# PERMANOVA for validation

scaled_quality <- scale(quality_complete[, 8:14])

dist_quality <- vegdist(scaled_quality, method = "euclidean")

set.seed(195)

permanova_traits <- 
  adonis2(
    dist_quality ~ species * dom_status * organ,
    permutations = 10000,
    data = quality_complete,
    by = "terms"
  )

permanova_traits



## PCA: leaves ----

quality_leaf <-
  quality_complete %>%
  filter(organ == "leaf")

pca_quality_leaf <- PCA(quality_leaf[, 8:14], 
                        scale.unit = TRUE, 
                        graph = FALSE)

summary(pca_quality_leaf) # 34.67%, 23.05%


# PERMANOVA for validation

scaled_quality_leaf <- scale(quality_leaf[, 8:14])

dist_quality_leaf <- vegdist(scaled_quality_leaf, method = "euclidean")

set.seed(195)

permanova_traits_leaf <- 
  adonis2(
    dist_quality_leaf ~ species * dom_status,
    permutations = 10000,
    data = quality_leaf,
    by = "terms"
  )

permanova_traits_leaf



## PCA: roots ----

quality_root <-
  quality_complete %>%
  filter(organ == "root")

pca_quality_root <- PCA(quality_root[, 8:14], 
                        scale.unit = TRUE, 
                        graph = FALSE)

summary(pca_quality_root) # 41.77%, 20.88%


# PERMANOVA for validation

scaled_quality_root <- scale(quality_root[, 8:14])

dist_quality_root <- vegdist(scaled_quality_root, method = "euclidean")

set.seed(195)

permanova_traits_root <- 
  adonis2(
    dist_quality_root ~ species * dom_status,
    permutations = 10000,
    data = quality_root,
    by = "terms"
  )

permanova_traits_root



## Mean values of traits ----

quality_w_leaf <- 
  quality_leaf %>%
  filter(
    dom_status == "w"
  ) %>%
  group_by(
    species,
    dom_status
  ) %>%
  summarise(across(c(icn, ca, mg, mn, p, tannins, dmc), mean, na.rm = TRUE))

quality_c_leaf <- 
  quality_leaf %>%
  filter(
    dom_status == "c"
  ) %>%
  group_by(
    species,
    dom_status
  ) %>%
  summarise(across(c(icn, ca, mg, mn, p, tannins, dmc), mean, na.rm = TRUE))


quality_c_leaf$icn - quality_w_leaf$icn
quality_c_leaf$ca - quality_w_leaf$ca
quality_c_leaf$mg - quality_w_leaf$mg
quality_c_leaf$mn - quality_w_leaf$mn
quality_c_leaf$p - quality_w_leaf$p
quality_c_leaf$tannins - quality_w_leaf$tannins
quality_c_leaf$dmc - quality_w_leaf$dmc


quality_w_root <- 
  quality_root %>%
  filter(
    dom_status == "w"
  ) %>%
  group_by(
    species,
    dom_status
  ) %>%
  summarise(across(c(icn, ca, mg, mn, p, tannins, dmc), mean, na.rm = TRUE))

quality_c_root <- 
  quality_root %>%
  filter(
    dom_status == "c"
  ) %>%
  group_by(
    species,
    dom_status
  ) %>%
  summarise(across(c(icn, ca, mg, mn, p, tannins, dmc), mean, na.rm = TRUE))


quality_c_root$icn - quality_w_root$icn
quality_c_root$ca - quality_w_root$ca
quality_c_root$mg - quality_w_root$mg
quality_c_root$mn - quality_w_root$mn
quality_c_root$p - quality_w_root$p
quality_c_root$tannins - quality_w_root$tannins
quality_c_root$dmc - quality_w_root$dmc



# H2 & 3: Decomposability ----


## Mass loss: leaves + roots ----

lmm_mass_loss <- lmer(
  ml_microcosm ~ dom_status * organ * sp_richness + (1|pfg_combo), 
  data = microcosm_master
)


# Check if the model is correct

plot(lmm_mass_loss, 
     type = c("p", "smooth"), 
     col.line = 1)

plot(lmm_mass_loss,
     sqrt(abs(resid(.))) ~ fitted(.),
     type = c("p", "smooth"), 
     col.line = 1)

qqmath(lmm_mass_loss)

plot(lmm_mass_loss, 
     rstudent(.) ~ hatvalues(.))


# Normality of residuals

par(mfrow = c(1,1))

hist(residuals(lmm_mass_loss)) # Good enough



## ANOVA: leaves + roots ----

anova(lmm_mass_loss)

ranova(lmm_mass_loss)

r.squaredGLMM(lmm_mass_loss)

F_to_eta2(
  f = c(0.0039,
        130.2828,
        0.1474,
        3.2138,
        0.9852,
        1.0077,
        1.4730),
  df = c(1,
         1,
         2,
         1,
         2,
         2,
         2),
  df_error = c(437.18,
               437.19,
               15.25,
               437.19,
               437.18,
               437.19,
               437.20)
)



# H4: Mixing effects ----


# Prepare data

# Calculate mean values of decomposition for leaves

leaf_w_m <-
  microcosm_master %>%
  filter(
    dom_status == "w" & organ == "leaf" & sp_richness == 1
  ) %>%
  pull(
    ml_microcosm
  ) %>%
  mean()


leaf_c_m <-
  microcosm_master %>%
  filter(
    dom_status == "c" & organ == "leaf" & sp_richness == 1
  ) %>%
  pull(
    ml_microcosm
  ) %>%
  mean()


# Calculate mean values of decomposition for roots

root_w_m <-
  microcosm_master %>%
  filter(
    dom_status == "w" & organ == "root" & sp_richness == 1
  ) %>%
  pull(
    ml_microcosm
  ) %>%
  mean()


root_c_m <-
  microcosm_master %>%
  filter(
    dom_status == "c" & organ == "root" & sp_richness == 1
  ) %>%
  pull(
    ml_microcosm
  ) %>%
  mean()

# Create data_tibble to test ME

ME_dummy <-
  microcosm_master %>%
  select(
    id,
    microcosm_random_id,
    dom_status,
    organ,
    sp_richness,
    pfg_combo,
    ml_microcosm
  ) %>%
  filter(
    sp_richness != 1
  )


# Wild leaves

ME_leaf_w <-
  ME_dummy %>%
  filter(
    dom_status == "w" & organ == "leaf"
  ) %>%
  mutate(
    ME = ((ml_microcosm - leaf_w_m) / leaf_w_m) * 100
  )


# Dom. leaves

ME_leaf_c <-
  ME_dummy %>%
  filter(
    dom_status == "c" & organ == "leaf"
  ) %>%
  mutate(
    ME = ((ml_microcosm - leaf_c_m) / leaf_c_m) * 100
  )


# Wild roots

ME_root_w <-
  ME_dummy %>%
  filter(
    dom_status == "w" & organ == "root"
  ) %>%
  mutate(
    ME = ((ml_microcosm - root_w_m) / root_w_m) * 100
  )


# Dom. roots

ME_root_c <-
  ME_dummy %>%
  filter(
    dom_status == "c" & organ == "root"
  ) %>%
  mutate(
    ME = ((ml_microcosm - root_c_m) / root_c_m) * 100
  )


# Merge partial tibbles

ME_all <-
  bind_rows(
    ME_leaf_w,
    ME_leaf_c,
    ME_root_w,
    ME_root_c
  )


t.test(ME_all$ME, mu = 0)



## Factors vs. ME ----


# Raw data model

lmm_me <- lmer(
  ME ~ dom_status * organ * sp_richness + (1|pfg_combo), 
  data = ME_all
)


# Check if the model is correct

plot(lmm_me, 
     type = c("p", "smooth"), 
     col.line = 1)

plot(lmm_me,
     sqrt(abs(resid(.))) ~ fitted(.),
     type = c("p", "smooth"), 
     col.line = 1)

qqmath(lmm_me)

plot(lmm_me, 
     rstudent(.) ~ hatvalues(.))


# Normality of residuals

par(mfrow = c(1,1))

hist(residuals(lmm_me)) # Good enough



### ANOVA ----

anova(lmm_me)

ranova(lmm_me)

r.squaredGLMM(lmm_me)

F_to_eta2(
  f = c(1.1204,
        5.2553,
        0.0944,
        5.3931,
        1.8620,
        0.0097,
        2.3137),
  df = c(1,
         1,
         1,
         1,
         1,
         1,
         1),
  df_error = c(345.53,
               345.55,
               6.15,
               345.57,
               345.53,
               345.55,
               345.57)
)



### Post hoc: pairwise t tests ---- 

# Data preparation

ME_post <-
  ME_all

ME_post$dom_status <- as.character(
  ME_post$dom_status
)

ME_post$organ <- as.character(
  ME_post$organ
)

ME_post <-
  ME_post %>%
  mutate(
    interaction = str_c(
      dom_status,
      organ,
      sep = "_"
    )
  )

#### Actual test ----

pairwise.t.test(
  ME_post$ME,
  ME_post$interaction,
  p.adjust.method = "BH"
)



# Figures ----

## Custom functions for figures ----

# Draws a PCA plot
pca_plot <- function(original_data, 
                     pca_result){
  
  # Create vectors with factor data
  sp_pca <- as.character(original_data$species)
  
  dom_pca <- as.character(original_data$dom_status)
  
  organ_pca <- as.character(original_data$organ)
  
  
  # Create factor based on two factors, so we can see
  # the effect of both at once in the same plot
  dom_organ_pca <- str_c(dom_pca, 
                         organ_pca,
                         sep = "_")
  
  
  # Isolate coordinates for PC1 and PC2
  dim_pca <- pca_result[["ind"]][["coord"]]
  
  dim_1 <- dim_pca[, 1]
  
  dim_2 <- dim_pca[, 2]
  
  
  # Build the tibble we want to plot
  data_ord <- as_tibble(cbind(sp_pca, dom_organ_pca, dim_1, dim_2))
  
  
  # Name elements of the tibble to simplify syntax when programming
  names(data_ord) <- c("sp", "dom_organ", "dim_1", "dim_2")
  
  
  # Parse factors
  data_ord$sp <- parse_factor(data_ord$sp)
  data_ord$dom_organ <- parse_factor(data_ord$dom_organ)
  
  
  # Parse numbers; we use as.numeric() because it coerces non-numeric to NA
  data_ord$dim_1 <- as.numeric(data_ord$dim_1)
  data_ord$dim_2 <- as.numeric(data_ord$dim_2)
  
  
  # Draw the plot
  basic <- ggplot(data_ord, 
                  aes(x = dim_1, 
                      y = dim_2, 
                      colour = sp, 
                      shape = dom_organ))
  
  basic + 
    geom_vline(xintercept = 0, linetype = "dashed") + 
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_point(size = 4) + #do not overlap full screen with 4
    theme_bw() + 
    theme(axis.title.x = element_text(size = 12, face = "plain", vjust = 0),
          axis.title.y = element_text(size = 12, face = "plain", vjust = 2),
          legend.title = element_text(size = 12, face = "plain"),
          legend.text = element_text(size = 10, face = "plain"))
}



# Draws traits boxplot
custom_box <- function(variable,
                       my_data = quality_complete){
  
  # Force to regular text
  variable_y <- sym(variable)
  
  # Get dom_status * organ factor
  my_data$dom_organ <- str_c(
    my_data$dom_status,
    my_data$organ,
    sep = "_"
  )
  
  my_data$dom_organ <- parse_factor(
    my_data$dom_organ,
    levels = c("w_leaf", "c_leaf", "w_root", "c_root"),
    ordered = TRUE
  )
  
  # Actual plot
  ggplot(data = my_data, 
         aes(x = dom_organ, 
             y = !!variable_y,
             fill = dom_organ)) +
    facet_wrap(~species, ncol = 4) +
    geom_boxplot() +
    geom_point(size = 1.5,
               alpha = 1/2) + 
    scale_fill_manual(name = NULL,
                      values = c("#4477AA", "#228833", "#CCBB44", "#EE6677"),
                      labels = c("Wild leaf", 
                                 "Dom. leaf",
                                 "Wild root",
                                 "Dom. root")) +
    scale_x_discrete(labels = c("Wild leaf", 
                                "Dom. leaf",
                                "Wild root",
                                "Dom. root")) +
    xlab(NULL) +
    theme_bw() + 
    theme(
      title = element_text(face = "bold", size = 8),
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
      legend.position = "top",
      panel.grid = element_blank(),
      panel.grid.major.y = element_line(colour = "black",
                                        linewidth = 0.3,
                                        linetype = "dotted")
    )
}



## Figure 2 ----

### PCA plot ----

p2_1 <- pca_plot(quality_complete, 
                 pca_quality,
                 var_colour = "organ",
                 var_shape = "organ") +
  scale_colour_manual(name = "Organ", 
                      values = c("#228833", "#CCBB44")) +
  xlab("Dimension 1 (30.25%)") +
  ylab("Dimension 2 (22.32%)") +
  theme(
    legend.position = "top"
  )

p2_1


### Loadings ----

# Check names
fviz_pca_var(
  pca_quality,
  arrowsize = 1,
  repel = TRUE)

# No names, so we can add them on other software for clarity

p2_2 <- fviz_pca_var(
  pca_quality,
  geom = c("arrow"),
  arrowsize = 1,
  repel = TRUE)

p2_2 <- p2_2 +
  ggtitle(NULL) + 
  xlab("Dimension 1 (31.7%)") +
  ylab("Dimension 2 (23.8%)") +
  theme(
    panel.grid = element_blank()
  )


ggarrange(
  p2_1,
  p2_2,
  ncol = 2,
  nrow = 1,
  common.legend = TRUE,
  legend = "top",
  labels = "auto")


ggsave("figs/beta/fig_2_beta.tiff", 
       width = 18, 
       height = 10,
       unit = "cm", 
       dpi = 300)



## Figure 3 ----


### Leaves' PCA ----


p3_1 <- pca_plot(quality_leaf, 
                  pca_quality_leaf) +
  geom_point(size = 3) + 
  scale_colour_manual(values = c("#E69F00", 
                                 "#56B4E9",
                                 "#009E73",
                                 "#F0E442",
                                 "#0072B2",
                                 "#D55E00",
                                 "#CC79A7",
                                 "#000000")) +
  xlab("Dimension 1 (34.67%)") +
  ylab("Dimension 2 (23.05%)") +
  theme(
    legend.position = "top",
    panel.grid = element_blank()
  )

p3_1



### Leaves' loadings ----

# Check names
fviz_pca_var(
  pca_quality_leaf,
  arrowsize = 1,
  repel = TRUE)

# No names, so we can add them on other software for clarity

p3_2 <- fviz_pca_var(
  pca_quality_leaf,
  geom = c("arrow"),
  arrowsize = 1,
  repel = TRUE)

p3_2 <- p3_2 +
  ggtitle(NULL) + 
  xlab("PC 1 (34.67%)") +
  ylab("PC 2 (23.05%)") +
  theme(
    panel.grid = element_blank()
  )

p3_2



### Roots' PCA ----
p3_3<- pca_plot(quality_root, 
                  pca_quality_root) +
  geom_point(size = 3) + 
  scale_colour_manual(values = c("#E69F00", 
                                 "#56B4E9",
                                 "#009E73",
                                 "#F0E442",
                                 "#0072B2",
                                 "#D55E00",
                                 "#CC79A7",
                                 "#000000")) +
  
  xlab("Dimension 1 (41.77%)") +
  ylab("Dimension 2 (20.88%)") +
  theme(
    legend.position = "top",
    panel.grid = element_blank()
  )

p3_3



### Roots' loadings ----

# Check names
fviz_pca_var(
  pca_quality_root,
  arrowsize = 1,
  repel = TRUE)

# No names, so we can add them on other software for clarity

p3_4 <- fviz_pca_var(
  pca_quality_root,
  geom = c("arrow"),
  arrowsize = 1,
  repel = TRUE)


p3_4 <- p3_4 +
  ggtitle(NULL) + 
  xlab("PC 1 (41.77%)") + 
  ylab("PC 2 (20.88%)") + 
  theme(
    panel.grid = element_blank()
  )

p3_4



### Whole plot ----


ggarrange(
  p3_1,
  p3_2,
  p3_3,
  p3_4,
  ncol = 2,
  nrow = 2,
  common.legend = TRUE,
  legend = "top",
  labels = "auto")

ggsave("figs/beta/fig_3_beta.tiff", width = 18, height = 18, unit = "cm", dpi = 300)



## Figure 4 ----

p4 <- ggplot(microcosm_master_plot, 
             aes(x = dom_organ, 
                 y = ml_microcosm)
)


set.seed(195)
p4 + 
  geom_violin(draw_quantiles = 0.5) +
  facet_wrap(~sp_richness, nrow = 1, 
             ncol = 3
  ) +
  geom_jitter(aes(shape = organ),
              size = 1.2,
              alpha = 1/1.5
  ) +
  scale_y_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80)) +
  xlab(NULL) +
  ylab("Mass loss (%)") +
  geom_hline(yintercept = 47.6, 
             linetype = "dashed",
             linewidth = 1,
             colour = "#228833") + 
  geom_hline(yintercept = 36.7, 
             linetype = "dashed", 
             linewidth = 1,
             colour = "#CCBB44") +
  scale_x_discrete(name = NULL, 
                   labels = c("W_leaf", 
                              "D_leaf", 
                              "W_root", 
                              "D_root")) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, 
                               vjust = 1, 
                               hjust = 1),
    legend.position = "top"
  )


ggsave("figs/beta/fig_4_beta.tiff", width = 14.56, height = 10.5, unit = "cm", dpi = 300)



## Figure 5 ----

p5 <- ggplot(ME_all, 
             aes(x = dom_status, 
                 y = ME)
)


set.seed(195)
p5 + 
  geom_violin(draw_quantiles = 0.5) +
  facet_wrap(~organ,
             labeller = labeller(group = custom_labels),
             nrow = 1, 
             ncol = 3
  ) +
  geom_jitter(aes(shape = dom_status),
              size = 1,
              alpha = 1/1.5) +
  xlab(NULL) +
  ylab("Mixture effects of mass loss (%)") +
  scale_x_discrete(name = NULL, 
                   labels = c("Wild", 
                              "Domesticated")) +
  scale_y_continuous(limits = c(-85, 130),
                     breaks = c(-75, -50, -25, 0, 25, 50, 75, 100, 125)) +
  geom_hline(yintercept = 4.95, 
             linetype = "dashed", 
             linewidth = 0.75,
             colour = "#B8A598") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "top"
  )


ggsave("figs/beta/fig_5_beta.tiff", width = 10.5, height = 10.5, unit = "cm", dpi = 300)



## Figure S1 ----

custom_box("in.") + 
  ylab("Nitrogen") +
  ggtitle("Litter Nitrogen by Dom. Status and Organ")

ggsave("figs/beta/supp/figS1a_beta.tiff", 
       width = 14, 
       height = 13, 
       unit = "cm", 
       dpi = 300)


custom_box("ic") + 
  ylab("Carbon") +
  ggtitle("Litter Carbon by Dom. Status and Organ")


ggsave("figs/beta/supp/figS1b_beta.tiff", 
       width = 14, 
       height = 13, 
       unit = "cm", 
       dpi = 300)


custom_box("icn") + 
  ylab("C/N") +
  ggtitle("Litter C/N by Dom. Status and Organ")


ggsave("figs/beta/supp/figS1c_beta.tiff", 
       width = 14,
       height = 13,
       unit = "cm",
       dpi = 300)


custom_box("ca") + 
  ylab("Calcium") +
  ggtitle("Litter Calcium by Dom. Status and Organ")


ggsave("figs/beta/supp/figS1d_beta.tiff",
       width = 14, 
       height = 13,
       unit = "cm",
       dpi = 300)


custom_box("mg") + 
  ylab("Magnesium") +
  ggtitle("Litter Magnesium by Dom. Status and Organ") 


ggsave("figs/beta/supp/figS1e_beta.tiff", 
       width = 14, 
       height = 13,
       unit = "cm", 
       dpi = 300)


custom_box("mn") + 
  ylab("Manganese") +
  ggtitle("Litter Manganese by Dom. Status and Organ")


ggsave("figs/beta/supp/figS1f_beta.tiff", 
       width = 14, 
       height = 13,
       unit = "cm", 
       dpi = 300)


custom_box("p") + 
  ylab("Phosphorus") +
  ggtitle("Litter Phosphorus by Dom. Status and Organ")


ggsave("figs/beta/supp/figS1g_beta.tiff", 
       width = 14, 
       height = 13, 
       unit = "cm", 
       dpi = 300)


custom_box("tannins") + 
  ylab("Tannins") +
  ggtitle("Litter Tannins by Dom. Status and Organ")


ggsave("figs/beta/supp/figS1h_beta.tiff", 
       width = 14, 
       height = 13, 
       unit = "cm", 
       dpi = 300)


custom_box("dmc") + 
  ylab("DMC (Dry matter content)") +
  ggtitle("Litter DMC by Dom. Status and Organ")


ggsave("figs/beta/supp/figS1i_beta.tiff",
       width = 14, 
       height = 13, 
       unit = "cm", 
       dpi = 300)
