library(yardstick)
library(pROC)
library(PRROC)
library(data.table)
library(tidyverse)
library(igraph)
library(ergm)
library(ergm.tapered)
library(readxl)
library(cvms)
library(viridis)
library(patchwork)
library(umap)
library(pheatmap)

#### JOINT COMPARISONS: 

# Prep graph: 

dryad <- fread(
  'final_predictions_90.tsv',
  sep = '\t',
  header = TRUE,
  quote = ''
)

dryad %>% glimpse() # 17849 edges
length(unique(c(dryad$Protein1, dryad$Protein2))) # 9806 

dryad <- dryad %>%
  filter(Protein1 != Protein2) %>%
  mutate(
    from = pmin(Protein1, Protein2),  # alphabetically smaller ID
    to = pmax(Protein1, Protein2)   # larger ID
  ) %>%
  distinct(from, to, .keep_all = TRUE)

dryad_edgelist <- dryad %>% dplyr::select(Protein1, Protein2)
glimpse(dryad)

length(unique(c(dryad$Protein1, dryad$Protein2))) # 9806

g <- graph_from_edgelist(dryad_edgelist %>% as.matrix(), directed = FALSE)

comp <- components(g)
summary(comp$csize)
hist(comp$csize)

g <- induced_subgraph(g, which(comp$membership == which.max(comp$csize)))
g_net <- intergraph::asNetwork(g)
g_net


## FEATURES 

dryad %>% glimpse()
length(unique(c(dryad$Process1, dryad$Process2)))
unique(c(dryad$Process1, dryad$Process2))

# easiest feature: Known, a feature that accounts for observation bias, each protein has a numerical score that measures how well it is characterized

known <- data.frame(protein = c(dryad$Protein1, dryad$Protein2), known_score = c(dryad$Known1, dryad$Known2))
known <- known %>% 
  distinct(.keep_all = TRUE)

known <- known %>% 
  filter((protein %in% V(g)$name))
ind_miss <- which(known$known_score == 'na')

known <- known %>% 
  mutate(known_score = ifelse(known_score == 'na', min(known_score, na.rm = TRUE), known_score)) %>%  # impute minimum value for caution, the presence of NA is probably biologically significant, poor experimental validation or non-presence in biologically known clusters
  mutate(known_score = as.numeric(known_score))

# subcellular_loc

protein_ids <- unique(names(V(g))) # proteins in the network's largest main component

node_list_1 <- dryad %>% 
  distinct(Locality1, Protein1) # distinct nodes and locations in the from column

node_list_2 <- dryad %>% 
  distinct(Locality2, Protein2) # distinct nodes and locations in the to column

node_list <- rbind(node_list_1, node_list_2, use.names = FALSE) # combine the 2 dataframes
node_list <- node_list %>%  
  distinct(.keep_all = TRUE) # keep only distinct pairs

node_list <- node_list %>% 
  filter(Protein1 %in% protein_ids) # keep only proteins in the main component

dryad_long <- node_list %>%
  separate_rows(Locality1, sep = ',') %>%
  mutate(Locality1 = str_trim(Locality1))   # remove extra spaces

location_mapping <- list(
  Nucleus = c(
    'Nucleus', 'Chromosome', 'Nucleosome core', 'Nuclear pore complex',
    'Telomere', 'Centromere', 'Kinetochore', 'Mitochondrion nucleoid'
  ),
  Cytoplasm = c(
    'Cytoplasm', 'Cytoskeleton', 'Intermediate filament', 'Keratin',
    'Microtubule', 'Cell projection', 'Sarcoplasmic reticulum',
    'Cytoplasmic vesicle', 'Microsome'
  ),
  Membrane = c(
    'Membrane', 'Cell membrane', 'Postsynaptic cell membrane',
    'Plasma membrane', 'Cell junction', 'Tight junction',
    'Gap junction', 'Coated pit', 'Signalosome'
  ),
  Endomembrane_Secretory = c(
    'Endosome', 'Golgi apparatus', 'Endoplasmic reticulum', 'Lysosome',
    'Vacuole', 'Peroxisome', 'Exosome', 'Lipid droplet', 'HDL', 'LDL', 'VLDL',
    'Secreted', 'Extracellular matrix', 'Surface film', 'Amyloid'
  ),
  Mitochondrion = c(
    'Mitochondrion', 'Mitochondrion inner membrane', 'Mitochondrion outer membrane'
  ),
  Complexes = c(
    'Proteasome', 'Spliceosome', 'Primosome', 'Signal recognition particle',
    'Synapse', 'Synaptosome', 'Cilium', 'Flagellum', 'Dynein',
    'MHC I', 'MHC II', 'T cell receptor', 'Immunoglobulin'
  )
)

assign_macro <- function(term) {
  for (categoria in names(location_mapping)) {
    if (term %in% location_mapping[[categoria]]) return(categoria)
  }
  return('Other')
}

dryad_long$macrolocation <- sapply(dryad_long$Locality1, assign_macro)
dryad_wide <- dryad_long %>%
  dplyr::select(Protein1, macrolocation) %>%
  distinct() %>%   # remove duplicates in case a protein appears multiple times in the same macrolocation
  mutate(value = 1) %>%  # create indicator for presence
  pivot_wider(
    names_from = macrolocation,
    values_from = value,
    values_fill = 0   # fill missing combinations with 0
  )

dryad_wide <- dryad_wide %>% 
  mutate(protein = Protein1) %>% 
  dplyr::select(-Protein1)

known <- known %>% 
  left_join(dryad_wide, by = 'protein')

###

vertex_data <- known %>%
  dplyr::filter(protein %in% network.vertex.names(g_net)) %>%
  dplyr::arrange(match(protein, network.vertex.names(g_net)))

write.csv(vertex_data, 'node_features.csv')
pheatmap::pheatmap(cor(vertex_data[,-1]), display_numbers = TRUE)

g_net %v% 'known_score'         <- vertex_data$known_score %>% as.numeric()
g_net %v% 'nucleus'         <- vertex_data$Nucleus %>% as.numeric()
g_net %v% 'cytoplasm'         <- vertex_data$Cytoplasm %>% as.numeric()
g_net %v% 'complexes'         <- vertex_data$Complexes %>% as.numeric()
g_net %v% 'membrane'         <- vertex_data$Membrane %>% as.numeric()
g_net %v% 'endomembrane_secretory'         <- vertex_data$Endomembrane_Secretory %>% as.numeric()
g_net %v% 'mitochondrion'         <- vertex_data$Mitochondrion %>% as.numeric()
g_net %v% 'other'         <- vertex_data$Other %>% as.numeric()

g <- intergraph::asIgraph(g_net, vnames = 'vertex.names')

## ERGM 

test <- read_xlsx('test_data.xlsx')
glimpse(test)

pred_m0 <- read.csv('m0_prediction_ergm.csv')
pred_m1 <- read.csv('m1_prediction_ergm.csv')
pred_m2 <- read.csv('m2_prediction_ergm.csv')
pred_m3 <- read.csv('m3_prediction_ergm.csv')
pred_m4 <- read.csv('m6_prediction_ergm.csv') # predictions from HOPE 

prediction <- tibble(truth = test$label) # CREO TIBBLE 
prediction <- prediction %>% 
  mutate(pred = pred_m0$p_hat_batch) %>%
  mutate(model = 'Model 0')  %>% 
  add_row(truth = test$label, pred = pred_m1$p_hat_batch, model = 'Model 1') %>% 
  add_row(truth = test$label, pred = pred_m2$p_hat_batch, model = 'Model 2') %>% 
  add_row(truth = test$label, pred = pred_m3$p_hat_batch, model = 'Model 3') %>% 
  add_row(truth = test$label, pred = pred_m4$p_hat_batch, model = 'Model 4')


#lift <- prediction %>% mutate(truth = factor(truth)) %>% group_by(model) %>%  ## creo colonna truth e raggruppo tutto per modello 
#  lift_curve(truth = truth, pred, event_level='second') %>% 
#  autoplot() + scale_color_brewer(palette = 'Set1')     
# LIFT not included in the final comparison but can be calculated like this (to view plot, just call lift)

roc <- prediction %>% mutate(truth = factor(truth, levels = c('0', '1'))) %>% group_by(model) %>% 
  roc_curve(truth, pred, event_level='second') %>% 
  autoplot() + scale_color_brewer(palette = 'Set1')

roc

auc <- prediction %>% mutate(truth = factor(truth, levels = c('0', '1'))) %>% group_by(model) %>% 
  roc_auc(truth, pred, event_level='second') %>% arrange(.estimate %>% desc())

auc %>%
  mutate(.estimate = format(.estimate, digits = 6, nsmall = 5)) ## AUC with extra digits



pr <- prediction %>% 
  mutate(truth = factor(truth, levels = c('0', '1'))) %>% 
  group_by(model) %>% 
  pr_curve(truth, pred, event_level = 'second') %>% 
  autoplot() + 
  scale_color_brewer(palette = 'Set1')

pr # precision-recall curve



pr_auc_res <- prediction %>% 
  mutate(truth = factor(truth, levels = c('0', '1'))) %>% 
  group_by(model) %>% 
  pr_auc(truth, pred, event_level = 'second') %>% 
  arrange(desc(.estimate))

pr_auc_res %>% 
  mutate(.estimate = format(.estimate, digits = 6, nsmall = 5)) # PR AUC values 




roc_df <- prediction %>% 
  mutate(truth = factor(truth, levels = c('0','1'))) %>% 
  group_by(model) %>% 
  roc_curve(truth, pred, event_level = 'second')

roc_plot <- ggplot(roc_df,
                   aes(x = 1 - specificity,
                       y = sensitivity,
                       color = model,
                       linewidth = model)) +
  geom_path() +
  scale_color_brewer(palette = 'Set1') +
  scale_linewidth_manual(values = c(
    'Model 0' = 0.4,
    'Model 1' = 0.4,
    'Model 2' = 0.4,
    'Model 3' = 0.4,
    'Model 4' = 1.1
  ), guide = 'none') +
  labs(title = 'ROC Curve',
       x = '1-Specificity',
       y = 'Sensitivity',
       color = 'Model') +
  theme_minimal(base_size = 11)

pr_df <- prediction %>% 
  mutate(truth = factor(truth, levels = c('0','1'))) %>% 
  group_by(model) %>% 
  pr_curve(truth, pred, event_level = 'second')

pr_plot <- ggplot(pr_df,
                  aes(x = recall,
                      y = precision,
                      color = model,
                      linewidth = model)) +
  geom_path() +
  scale_color_brewer(palette = 'Set1') +
  scale_linewidth_manual(values = c(
    'Model 0' = 0.4,
    'Model 1' = 0.4,
    'Model 2' = 0.4,
    'Model 3' = 0.4,
    'Model 4' = 1.1
  ), guide = 'none') +
  labs(title = 'Precision–Recall Curve',
       x = 'Recall',
       y = 'Precision',
       color = 'Model') +
  theme_minimal(base_size = 11) +ylim(0,1)

(roc_plot + pr_plot) +
  plot_layout(guides = 'collect') &
  theme(legend.position = 'bottom') # ROC and PR curves with highlighted best model


# lift (not included but the same format as ROC and PR is implemented here)

lift_df <- prediction %>% 
  mutate(truth = factor(truth, levels = c('0','1'))) %>% 
  group_by(model) %>% 
  lift_curve(truth, pred, event_level = 'second') 

lift_plot <- ggplot(lift_df,
                  aes(x = .percent_tested,
                      y = .lift,
                      color = model,
                      linewidth = model)) +
  geom_path() +
  scale_color_brewer(palette = 'Set1') +
  scale_linewidth_manual(values = c(
    'Model 0' = 0.4,
    'Model 1' = 0.4,
    'Model 2' = 0.4,
    'Model 3' = 0.4,
    'Model 4' = 1.1
  ), guide = 'none') +
  labs(title = 'Lift Curve',
       x = '% Tested',
       y = 'Lift',
       color = 'Model') +
  theme_minimal(base_size = 11) +
  geom_hline(yintercept = 1, linetype = 'dotted', linewidth = 0.3)

lift_plot



heldout <- read_xlsx('test_data.xlsx') # reload everything for ease 
heldout <- heldout %>%
  mutate(pred_m3 = pred_m3$p_hat_batch,
         pred_m4 = pred_m4$p_hat_batch) # add predictions for the 3rd and 4th model (gwesp vs gwnsp) into the dataframe

# top k precision

prediction %>%
  group_by(model) %>%
  arrange(desc(pred), .by_group = TRUE) %>%
  mutate(top_5pct = row_number() <= 0.05 * n()) %>%
  summarise(
    precision_top5 = mean(truth[top_5pct])
  )

prediction %>%
  group_by(model) %>%
  arrange(desc(pred), .by_group = TRUE) %>%
  mutate(top_5pct = row_number() <= 0.1 * n()) %>%
  summarise(
    precision_top5 = mean(truth[top_5pct])
  )

prediction %>%
  group_by(model) %>%
  arrange(desc(pred), .by_group = TRUE) %>%
  mutate(top_5pct = row_number() <= 0.2 * n()) %>%
  summarise(
    precision_top5 = mean(truth[top_5pct])
  )

###########################################################

dryad <- dryad %>%
  mutate(
    key = ifelse(Protein1 < Protein2,
                 paste0(Protein1, '_', Protein2),
                 paste0(Protein2, '_', Protein1))
  )

heldout <- heldout %>%
  mutate(key = ifelse(Protein1 < Protein2, paste0(Protein1, '_', Protein2), paste0(Protein2, '_', Protein1)))

heldout <- heldout %>%
  left_join(dryad %>% select(key, Process1, Process2), by = 'key')

same_process <- function(p1, p2) {
  if (is.na(p1) || is.na(p2)) return(FALSE)
  set1 <- tolower(strsplit(p1,',')[[1]])
  set2 <- tolower(strsplit(p2, ',')[[1]])
  length(intersect(set1, set2)) >0
} # check if 2 nodes in an edge have the same process (similarly to the HOPE implementation, i normalize the strings and check for an intersection)

heldout <- heldout %>%
  mutate( same_process = mapply(same_process, Process1, Process2)) # calculate a same_process variable fro each heldout edge

rate_m3 <- heldout %>%
  arrange(desc(pred_m3)) %>%
  mutate(top5 = row_number() <= 0.2*n()) %>%
  summarise(rate = mean(same_process[top5], na.rm = TRUE))

rate_m3

rate_m4 <- heldout %>%
  arrange(desc(pred_m4)) %>%
  mutate(top5 = row_number() <= 0.2*n()) %>%
  summarise(rate = mean(same_process[top5], na.rm = TRUE))

rate_m4 # calculate the rate at which top 20% of predicitons shared the same protein complex (using data from the original dataset)



################################################
################################################
################################################

##### GNN  #####################################

################################################
################################################
################################################


test <- read_xlsx('test_data.xlsx')
glimpse(test)

pred_m0 <- read.csv('m0_pred_gnn.csv')
pred_m1 <- read.csv('m1_pred_gnn.csv')
pred_m2 <- read.csv('m2_pred_gnn.csv')
pred_m3 <- read.csv('m3_pred_gnn.csv')
pred_m4 <- read.csv('m4_pred_gnn.csv')
pred_m5 <- read.csv('m5_pred_gnn.csv')


prediction <- tibble(truth = test$label) ## CREO TIBBLE 
prediction <- prediction %>% 
  mutate(pred = pred_m0$score) %>%
  mutate(model = 'Model 0')  %>% 
  add_row(truth = test$label, pred = pred_m1$score, model = 'Model 1') %>% 
  add_row(truth = test$label, pred = pred_m2$score, model = 'Model 2') %>% 
  add_row(truth = test$label, pred = pred_m3$score, model = 'Model 3') %>% 
  add_row(truth = test$label, pred = pred_m4$score, model = 'Model 4') %>%
  add_row(truth = test$label, pred = pred_m5$score, model = 'Model 5')
  


lift <- prediction %>% mutate(truth = factor(truth)) %>% group_by(model) %>%  ## creo colonna truth e raggruppo tutto per modello 
  lift_curve(truth = truth, pred, event_level='second') %>% 
  autoplot() + scale_color_brewer(palette = 'Set1')

roc <- prediction %>% mutate(truth = factor(truth, levels = c('0', '1'))) %>% group_by(model) %>% 
  roc_curve(truth, pred, event_level='second') %>% 
  autoplot() + scale_color_brewer(palette = 'Set1')

roc
lift

auc <- prediction %>% mutate(truth = factor(truth, levels = c('0', '1'))) %>% group_by(model) %>% 
  roc_auc(truth, pred, event_level='second') %>% arrange(.estimate %>% desc())

auc %>%
  mutate(.estimate = format(.estimate, digits = 6, nsmall = 5)) ## for extra digits

auc


pr <- prediction %>% 
  mutate(truth = factor(truth, levels = c('0', '1'))) %>% 
  group_by(model) %>% 
  pr_curve(truth, pred, event_level = 'second') %>% 
  autoplot() + 
  scale_color_brewer(palette = 'Set1')


pr_auc_res <- prediction %>% 
  mutate(truth = factor(truth, levels = c('0', '1'))) %>% 
  group_by(model) %>% 
  pr_auc(truth, pred, event_level = 'second') %>% 
  arrange(desc(.estimate))

pr_auc_res %>% 
  mutate(.estimate = format(.estimate, digits = 6, nsmall = 5)) # all of this is exactly the same as the ERGM


roc_df <- prediction %>% 
  mutate(truth = factor(truth, levels = c('0','1'))) %>% 
  group_by(model) %>% 
  roc_curve(truth, pred, event_level = 'second')

roc_plot <- ggplot(roc_df,
                   aes(x = 1 - specificity,
                       y = sensitivity,
                       color = model,
                       linewidth = model)) +
  geom_path() +
  scale_color_brewer(palette = 'Set1') +
  scale_linewidth_manual(values = c(
    'Model 0' = 1,
    'Model 1' = 0.2,
    'Model 2' = 0.2,
    'Model 3' = 0.2,
    'Model 4' = 0.2,
    'Model 5' = 0.2
  ), guide = 'none') +
  labs(title = 'ROC Curve',
       x = '1-Specificity',
       y = 'Sensitivity',
       color = 'Model') +
  theme_minimal(base_size = 11)

pr_df <- prediction %>% 
  mutate(truth = factor(truth, levels = c('0','1'))) %>% 
  group_by(model) %>% 
  pr_curve(truth, pred, event_level = 'second')

pr_plot <- ggplot(pr_df,
                  aes(x = recall,
                      y = precision,
                      color = model,
                      linewidth = model)) +
  geom_path() +
  scale_color_brewer(palette = 'Set1') +
  scale_linewidth_manual(values = c(
    'Model 0' = 1,
    'Model 1' = 0.2,
    'Model 2' = 0.2,
    'Model 3' = 0.2,
    'Model 4' = 0.2,
    'Model 5' = 0.2
  ), guide = 'none') +
  labs(title = 'Precision–Recall Curve',
       x = 'Recall',
       y = 'Precision',
       color = 'Model') +
  theme_minimal(base_size = 11)

(roc_plot + pr_plot) +
  plot_layout(guides = 'collect') &
  theme(legend.position = 'bottom')

lift_df <- prediction %>% 
  mutate(truth = factor(truth, levels = c('0','1'))) %>% 
  group_by(model) %>% 
  lift_curve(truth, pred, event_level = 'second')

lift_plot <- ggplot(lift_df,
                    aes(x = .percent_tested,
                        y = .lift,
                        color = model,
                        linewidth = model)) +
  geom_path() +
  scale_color_brewer(palette = 'Set1') +
  scale_linewidth_manual(values = c(
    'Model 0' = 1,
    'Model 1' = 0.2,
    'Model 2' = 0.2,
    'Model 3' = 0.2,
    'Model 4' = 0.2,
    'Model 5' = 0.2
  ), guide = 'none') +
  labs(title = 'Lift Curve',
       x = '% Tested',
       y = 'Lift',
       color = 'Model') +
  theme_minimal(base_size = 11)

lift_plot

# TOP k PRECISION

prediction %>%
  group_by(model) %>%
  arrange(desc(pred), .by_group = TRUE) %>%
  mutate(top_5pct = row_number() <= 0.05 * n()) %>%
  summarise(
    precision_top5 = mean(truth[top_5pct])
  )

prediction %>%
  group_by(model) %>%
  arrange(desc(pred), .by_group = TRUE) %>%
  mutate(top_5pct = row_number() <= 0.1 * n()) %>%
  summarise(
    precision_top5 = mean(truth[top_5pct])
  )

prediction %>%
  group_by(model) %>%
  arrange(desc(pred), .by_group = TRUE) %>%
  mutate(top_5pct = row_number() <= 0.2 * n()) %>%
  summarise(
    precision_top5 = mean(truth[top_5pct])
  )

pred_unif <- read.csv('m_pred_gnn_unif.csv')


## UMAP embeddings 

message_embeddings <- read.csv('gnn_node_embeddings.csv') # from Python output
glimpse(message_embeddings)

emb_only <- message_embeddings %>%
  select(starts_with('X')) # all embeddings start with X in the csv

cor_mat <- cor(emb_only)

pheatmap(cor_mat,
         show_rownames = FALSE,
         show_colnames = FALSE,
         main = 'Correlation Between Embedding Dimensions') # correlation matrix heatmap for node embeddings

umap_layout <- umap(
  as.matrix(emb_only),
  n_neighbors = 30,
  min_dist = 0.9,
  metric = 'cosine')

top10 <- message_embeddings %>% 
  count(community) %>%
  slice_max(n, n = 10) %>%
  pull(community) # top 10 louvain communities

message_embeddings <- message_embeddings %>%
  mutate(community_top10 = ifelse(community %in% top10, community, 'Other')) # add a column that shows a node's community only if they are in the top 10 (alll others are aggregated into an 'Other' category)

umap_df <- data.frame(
  UMAP1 = umap_layout$layout[,1],
  UMAP2 = umap_layout$layout[,2],
  community = as.factor(message_embeddings$community_top10))  # dataframe with umap results 

umap_df <- umap_df %>%
  mutate(point_size = ifelse(community == 'Other', 0.5, 1.1)) # add a point size variable to accentuate the top 10 in the plot (bigger if in the top 10 communities)

levels <- levels(umap_df$community) # factor levels

color_map <- setNames(rep('grey', length(levels)), levels) # start everything as grey

top_levels <- setdiff(levels, 'Other') # assign colors to top communities (not 'Other')
palette <- viridis(length(top_levels), option = 'turbo') # use the viridis palette (best looking from my perspective)

color_map[top_levels] <- palette 


ggplot(umap_df, aes(UMAP1, UMAP2, color = community, size = point_size)) +
  geom_point(alpha = 0.8) +
  theme_minimal() +
  labs(title = 'UMAP of Test Set Node Embeddings',
       color = 'Community') +
  scale_color_manual(values = color_map, na.value = 'grey') +
  scale_size_identity()  # tells ggplot to use the exact sizes in the column


### 7 different umap for each nodematch term

plot_df <- umap_df %>%
  mutate(Nucleus = vertex_data$Nucleus,
         Cytoplasm = vertex_data$Cytoplasm,
         Complexes = vertex_data$Complexes,
         Membrane = vertex_data$Membrane,
         Endomembrane = vertex_data$Endomembrane_Secretory,
         Mitochondrion = vertex_data$Mitochondrion,
         Other = vertex_data$Other) %>%
  pivot_longer(
    cols = Nucleus:Other,
    names_to = 'Compartment',
    values_to = 'Presence')
ggplot(plot_df, aes(UMAP1, UMAP2, color = as.factor(Presence))) +
  geom_point(size = 0.6, alpha = 0.8) +
  facet_wrap(~Compartment, ncol = 3) +
  scale_color_manual(
    values = c('0' = 'lightblue', '1' = 'red'),
    name = 'Presence') +
  theme_minimal() +
  theme(legend.position = 'right',
        strip.text = element_text(face = 'bold'))


##################  JOINT  ANALYSIS ######################

library(openxlsx)
pred_gnn <- read.csv('m0_pred_gnn.csv')
pred_ergm <- read.csv('m6_prediction_ergm.csv')

predictions <- data.frame(pred_ergm = pred_ergm$p_hat_batch,
                          pred_gnn = pred_gnn$score,
                          label = heldout$label)


cor(predictions$pred_ergm, predictions$pred_gnn, method = 'spearman') # correlation between the 2 prediction vectors

roc_ergm <- roc(heldout$label, pred_ergm$p_hat_batch)
roc_gnn <- roc(predictions$label, predictions$pred_gnn)
auc(roc_ergm)
auc(roc_gnn) # AUC values (as seen before)

youden_ergm <- coords(
  roc_ergm,
  x = 'best',
  best.method = 'youden',
  ret = c('threshold', 'sensitivity', 'specificity', 'youden')) # youden j statistic

youden_gnn <- coords(
  roc_gnn,
  x = 'best',
  best.method = 'youden',
  ret = c('threshold', 'sensitivity', 'specificity', 'youden')) # youden j statistic


predictions <- predictions %>% 
  mutate(prediction_label_ergm = ifelse(pred_ergm < youden_ergm$threshold, 0, 1),
         prediction_label_gnn = ifelse(pred_gnn < youden_gnn$threshold, 0, 1))




predictions_long <- predictions %>%
  select(label, pred_ergm, pred_gnn) %>%
  pivot_longer(-label, names_to = 'model', values_to = 'prob')

ggplot(predictions_long, aes(x = prob, fill = factor(label))) +
  geom_density(alpha = 0.4) +
  facet_wrap(~model, scales = 'free') +   # big boy
  theme_minimal() +
  labs(fill = 'True Label') # predictions stratified by label (ergm is almost illegible so log scale is needed)


predictions_long <- predictions_long %>%
  mutate(prob_log = log(prob)) #log scale


pl1 <- ggplot(predictions) +
  geom_density(aes(x = log(pred_ergm), fill = factor(label)), alpha = 0.4) +
  labs(x = 'ERGM log-predictions',
       y = 'Density',
       fill = 'Label') 

pl2 <- ggplot(predictions) +
  geom_density(aes(x = pred_gnn, fill = factor(label)), alpha = 0.4) +
  labs(x = 'GNN predictions',
       y = 'Density',
       fill = 'Label')

combined <- (pl1 + pl2) +
  plot_layout(guides = 'collect') &
  theme(legend.position = 'right') &
  plot_annotation(title = 'Prediction Seperation via Density') # combine the 2 plots

combined

#### CM 

make_cm <- function(pred, true) {
  as.data.frame(table(Pred = pred, True = true))
} # confusion matrix dataframe

cm_ergm <- make_cm(predictions$prediction_label_ergm, heldout$label)
cm_gnn  <- make_cm(predictions$prediction_label_gnn,  heldout$label)

plot_cm <- function(cm, title){
  ggplot(cm, aes(x = Pred, y = True, fill = Freq)) +
    geom_tile(color = 'white') +
    geom_text(aes(label = Freq,
              color = ifelse(Pred == True, 'white', 'black')),
              size = 4) +
    scale_fill_gradient(low = '#f7fbff', high = '#08306b') + # manually inputted to mimic the Python CM colors
    scale_color_identity() +
    scale_x_discrete(limits = c('0', '1')) +
    scale_y_discrete(limits = c('1', '0')) +
    labs(title = title,
         x = 'Predicted label',
         y = 'True label') +
    theme_classic() +
    theme(panel.grid = element_blank(),
          legend.position = 'none')
} # function to plot the CM

p1 <- plot_cm(cm_ergm, 'ERGM Confusion Matrix')
p2 <- plot_cm(cm_gnn,  'GNN Confusion Matrix')

p1 + p2
