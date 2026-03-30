library(data.table)
library(tidyverse)
library(igraph)
library(ergm)
library(poweRlaw)
library(ggrepel)

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
  mutate(from = pmin(Protein1, Protein2),  # alphabetically smaller ID
        to = pmax(Protein1, Protein2)) %>%   # larger ID
  distinct(from, to, .keep_all = TRUE)

dryad_edgelist <- dryad %>% dplyr::select(Protein1, Protein2)
glimpse(dryad)

length(unique(c(dryad$Protein1, dryad$Protein2))) # 9806

g <- graph_from_data_frame(dryad_edgelist, directed = FALSE)

comp <- components(g)
summary(comp$csize)
hist(comp$csize)

g <- induced_subgraph(g, which(comp$membership == which.max(comp$csize)))
g_net <- intergraph::asNetwork(g)
g_net


## FEATURES 

dryad %>% glimpse()

# easiest feature: Known, a feature that accounts for observation bias, each protein has a numerical score that measures how well it is characterized (unused in final implementations as it didn't prove to provide any useful signal)

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
  pivot_wider(names_from =  macrolocation, values_from = value, values_fill = 0)   # fill missing combinations with 0

dryad_wide <- dryad_wide %>% 
  mutate(protein = Protein1) %>% 
  dplyr::select(-Protein1)

known <- known %>% 
  left_join(dryad_wide, by = 'protein')


vertex_data <- known %>%
  dplyr::filter(protein %in% network.vertex.names(g_net)) %>%
  dplyr::arrange(match(protein, network.vertex.names(g_net)))

pheatmap::pheatmap(cor(vertex_data[,-1]), display_numbers = TRUE) # correlation matrix of features

g_net %v% 'known_score'         <- vertex_data$known_score %>% as.numeric()
g_net %v% 'nucleus'         <- vertex_data$Nucleus %>% as.numeric()
g_net %v% 'cytoplasm'         <- vertex_data$Cytoplasm %>% as.numeric()
g_net %v% 'complexes'         <- vertex_data$Complexes %>% as.numeric()
g_net %v% 'membrane'         <- vertex_data$Membrane %>% as.numeric()
g_net %v% 'endomembrane_secretory'         <- vertex_data$Endomembrane_Secretory %>% as.numeric()
g_net %v% 'mitochondrion'         <- vertex_data$Mitochondrion %>% as.numeric()
g_net %v% 'other'         <- vertex_data$Other %>% as.numeric()

g <- intergraph::asIgraph(g_net)
V(g)$name <- network.vertex.names(g_net)

#### DESCRIPTIVE 

edge_density(g)
deg <- degree(g)
deg %>% summary()
hist(deg, breaks = seq(0,102,2))

ggplot(data.frame(degree = degree(g))) + 
  geom_histogram(
    aes(x = degree), bins = 100,
    fill = ggthemes::economist_pal()(3)[3],
    alpha = 0.7, colour = 'black') +
  xlab('Node Degree') +
  ylab('Absolute Frequency') +
  ggtitle('Degree Distribution') +
  theme(axis.title = element_text(size = 10),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10)),
        title = element_text(size = 10),
        panel.background = element_rect(fill = 'white'),
        plot.background = element_rect(fill = 'white')) # degree distribution histogram


clusters <- cluster_louvain(g) # louvain clustering
clusters
igraph::assortativity(g, directed = FALSE, values = clusters$membership) # assortativity
modularity(g, membership = clusters$membership) # modularity
transitivity(g, type = 'undirected') # transitivity

#

table(dryad_wide$Nucleus) # absolute frequency tables for each location
table(dryad_wide$Cytoplasm)
table(dryad_wide$Complexes)
table(dryad_wide$Membrane)
table(dryad_wide$Endomembrane_Secretory)
table(dryad_wide$Mitochondrion)
table(dryad_wide$Other)



#####

deg <- degree(g)
deg <- deg[deg>0]   # remove isolated nodes(power law defined for k >=1)

pl <- displ$new(deg) # Fit power law

est <- estimate_xmin(pl) # to fit a power law we need to find the minimum degree from which the power law hypothesis holds 
pl$setXmin(est)
pl$setPars(estimate_pars(pl))

cat('xmin =', pl$getXmin(), '\n') # xmin is 1, showing the whole distribution follows a circa power law distribution
cat('gamma =', pl$pars, '\n') # slope gamma parameter

bs <- bootstrap_p(pl, no_of_sims = 1000) # gof test to see how plausible this is 
cat('power law p-value =', bs$p, '\n')


ln   <- dislnorm$new(deg)  # log normal model
expd <- disexp$new(deg)

ln$setXmin(est$xmin)      # same cutoff
ln$setPars(estimate_pars(ln))

expd$setXmin(est$xmin)    # same cutoff
expd$setPars(estimate_pars(expd))

comp_ln  <- compare_distributions(pl, ln) # Vuong test 1
comp_exp <- compare_distributions(pl, expd) # Vuong test 2 
?compare_distributions

comp_ln$test_statistic
comp_ln$p_two_sided

comp_exp$test_statistic
comp_exp$p_two_sided



# Degree distribution
deg <- degree(g)
deg_tab <- as.data.frame(table(deg))
colnames(deg_tab) <- c('k', 'freq')
deg_tab$k <- as.numeric(as.character(deg_tab$k))
deg_tab$freq <- deg_tab$freq / sum(deg_tab$freq)


pl <- displ$new(deg) # fit power law
est <- estimate_xmin(pl)
pl$setPars(estimate_pars(pl))

ln <- dislnorm$new(deg) # fit lognormal
ln$setPars(estimate_pars(ln))

expd <- disexp$new(deg) # fit lognormal
expd$setPars(estimate_pars(expd))


xmin <- 1 # calculated before
deg_tab_fit <- deg_tab %>% filter(k >= xmin)

kvals <- deg_tab_fit$k

fit_df <- data.frame( # fitted curves at all k values 
  k = kvals,
  powerlaw = dist_pdf(pl, kvals),
  lognormal = dist_pdf(ln, kvals),
  exponential = dist_pdf(expd, kvals))

plot_df <- left_join(deg_tab_fit, fit_df, by='k') # merge for final plot

ggplot(plot_df, aes(x = k)) +
  
  # Empirical distribution
  geom_point(
    aes(y = freq),
    color = 'black',
    size = 2,
    alpha = 0.75) + # scatter plot
  geom_line(aes(y = powerlaw,  color = 'Power-law'), linewidth = 0.8) +
  geom_line(aes(y = lognormal, color = 'Log-normal'), linewidth = 0.8) +
  geom_line(aes(y = exponential, color = 'Exponential'), linewidth = 0.8) + # the 3 curves
  scale_x_log10() + # scaling for visibility
  scale_y_log10(breaks = scales::trans_breaks('log10', function(x) 10^x),
                labels = scales::trans_format('log10', scales::math_format(10^.x))) +
  xlab('Degree (k)') +
  ylab('Probability P(k)') +
  ggtitle('Goodness of Fit for Degree Distribution Models (log-log scale)') +
  ggthemes::theme_fivethirtyeight() +
  theme(axis.title = element_text(size = 10),
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.title.y = element_text(margin = margin(r = 10)),
      title = element_text(size = 10),
    panel.background = element_rect(fill = 'white'),
    plot.background = element_rect(fill = 'white'),
    legend.background = element_rect(fill = 'white')
  ) 


#### ANALYSIS FOR CHAPTER

g_net # 8273 nodes, 16638 edges, undirected, no missing edges, unipartite, no self loops

edge_density(g) # edge density incredibly low: 0.0004862477 (sparse)

deg <- data.frame(deg = degree(g), name = V(g)$name)

mean_distance(g, directed = FALSE) # 8.148856
diameter(g, directed = FALSE) # 25, considering the size of the network, it is not very connected, although this was already evident from the density and other stats
?mean_distance

# identification of hubs
which(deg$deg > 50)
jj <- deg[c(which(deg$deg > 20)),] # nodes with highest degrees
jj$deg %>% length()


## 

betweenness_ppi <- betweenness(g, directed = FALSE) 
betweenness_ppi %>% summary()
hist(betweenness_ppi)

which(betweenness_ppi > quantile(betweenness_ppi, 0.999)) # nodes with highest betweenness
betweenness_ppi <- data.frame(betweenness = betweenness_ppi, name = V(g)$name )

ppi_metrics <- betweenness_ppi %>%
  left_join(deg, by = 'name') # merge for plot

ind <- which(ppi_metrics$betweenness > 20)
ppi_metrics[ind,]

hubs <- ppi_metrics[ind, ] # hub node dataframe
hubs


# Optionally, highlight top nodes by degree or betweenness
top_nodes <- ppi_metrics %>%
  dplyr::filter(deg >= quantile(deg, 0.99) | betweenness >= quantile(betweenness, 0.95)) # find the top degree and betweenness nodes

ggplot(hubs, aes(x = deg, y = betweenness)) +
  geom_point(aes(color = betweenness, size = deg), alpha = 0.6) +  # color by betweenness, size by degree
  geom_smooth(method = 'lm', se = FALSE, color = 'green', linetype = 'dashed') +
  geom_text_repel(data = top_nodes, aes(label = name), size = 3, max.overlaps = 15) +
  scale_x_log10(labels = scales::comma_format()) +
  scale_y_log10(labels = scales::comma_format()) +
  scale_color_viridis_c(option = 'C', trans = 'log', guide = guide_colorbar(title = 'Betweenness')) +
  scale_size_continuous(range = c(2, 6), guide = guide_legend(title = 'Degree')) +
  labs(title = 'Relationship Between Degree and Betweenness Centrality',
       subtitle = 'Nodes with highest degree or betweenness are labeled',
       x = 'Degree',
       y = 'Betweenness Centrality') +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = 'bold', size = 14),
        plot.subtitle = element_text(size = 10, color = 'gray30'),
        legend.position = 'right',
        panel.grid.minor = element_blank())

library(ggplot2)
library(ggrepel)
library(dplyr)
library(scales)

cor(hubs$deg, hubs$betweenness, method = 'pearson')
