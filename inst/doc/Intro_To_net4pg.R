## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 10, fig.height = 10)

## ----style, echo = FALSE, results = 'asis'------------------------------------
BiocStyle::markdown()

## ----message = FALSE, include = FALSE-----------------------------------------
library(net4pg)
library(igraph)
library(ggplot2)

## -----------------------------------------------------------------------------
incM_filename <- system.file("extdata"
                        , "incM_example"
                        , package = "net4pg"
                        , mustWork = TRUE)
rownames_filename <- system.file("extdata"
                        , "peptideIDs_incM_example"
                        , package = "net4pg"
                        , mustWork = TRUE)
colnames_filename <- system.file("extdata"
                        , "proteinIDs_incM_example"
                        , package = "net4pg"
                        , mustWork = TRUE)
incM <- read_inc_matrix(incM_filename = incM_filename
                , colnames_filename = colnames_filename
                , rownames_filename = rownames_filename)

## -----------------------------------------------------------------------------
dim(incM)

## -----------------------------------------------------------------------------
incM_reduced <- reduce_inc_matrix(incM)
dim(incM_reduced) # check the size of the reduced incidence matrix

## -----------------------------------------------------------------------------
adjM <- get_adj_matrix(incM_reduced)
dim(adjM) # check the size of the adjacency matrix:

## -----------------------------------------------------------------------------
multProteinCC <- get_cc(adjM)

## -----------------------------------------------------------------------------
cc.multProteins <- multProteinCC$ccs
length(cc.multProteins)

## -----------------------------------------------------------------------------
# Calculate CCs size and percentage of single- vs multi-protein CCs
CCstatsOut <- cc_stats(incM = incM
                       , cc.proteins = cc.multProteins
                       , reducedIncM = TRUE)

# Number of single-protein CCs:
CCstatsOut$N_singleProtCC

# Number of multi-protein CCs
CCstatsOut$N_multiProtCC

# Total number of CCs
totCCs <- CCstatsOut$N_singleProtCC + CCstatsOut$N_multiProtCC
totCCs

# Percentage of single-protein CCs:
PercSingleP <- round(CCstatsOut$N_singleProtCC / totCCs * 100, digits = 2)
PercSingleP

# View table of CC size distribution
CCstatsOut$NproteinsDistribution

# Plot CC size distribution
plot(factor(CCstatsOut$NproteinsDistribution$N_proteins
       , levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", ">10"))
     , as.numeric(as.vector(CCstatsOut$NproteinsDistribution$N_CC))
     , type = "s"
     , xlab = "N_proteins"
     , ylab = "N_CCs")

## -----------------------------------------------------------------------------
peptideStatsOut <- peptide_stats(incM = incM)

# Number of shared peptides
peptideStatsOut$nbShared

# Number of specific peptides
peptideStatsOut$nbSpecific

# Percentage of specific peptides
peptideStatsOut$percSpecific

## -----------------------------------------------------------------------------
cc.peptides.incM <- cc_composition(cc.multProteins, incM = incM)

## ---- fig.height=7------------------------------------------------------------
# Generate the bipartite graph
prot <- "ENSP261"
subgraphCC <- plot_cc(prot = prot
        , cc.proteins = cc.multProteins
        , cc.subincM = cc.peptides.incM$cc.subincM
        , incM = incM
        , tagProt = "ENSP"
        , tagContam = "Contam")
# Plot it
plot.igraph(subgraphCC$g
            , layout = layout_as_bipartite
            , edge.width = 1
            ,  edge.arrow.width = 0.3
            , vertex.size = 35
            , edge.arrow.size = 0.5
            , vertex.size2 = 35
            , vertex.label.cex = 1
            , asp = 0.25
            , margin = -0.1) +
title(paste0("Protein ", prot, " in CC#", subgraphCC$cc_id), line = -1)

## -----------------------------------------------------------------------------
incM_filename <- system.file("extdata"
                        , "incM_example"
                        , package = "net4pg"
                        , mustWork = TRUE)
rownames_filename <- system.file("extdata"
                        , "peptideIDs_incM_example"
                        , package = "net4pg"
                        , mustWork = TRUE)
colnames_filename <- system.file("extdata"
                        , "proteinIDs_incM_example"
                        , package = "net4pg"
                        , mustWork = TRUE)
incM <- read_inc_matrix(incM_filename = incM_filename
                , colnames_filename = colnames_filename
                , rownames_filename = rownames_filename)

## -----------------------------------------------------------------------------
dim(incM)

## -----------------------------------------------------------------------------
# Read input file names
exprTranscriptsFile <- system.file("extdata"
                        , "expressed_transcripts.txt"
                        , package = "net4pg"
                        , mustWork = TRUE)
protein2transcriptFile <- system.file("extdata"
                        , "protein_to_transcript"
                        , package = "net4pg"
                        , mustWork = TRUE)

# Perform filtering
incM_filt <- transcriptome_filter(incM
                            , exprTranscriptsFile = exprTranscriptsFile
                            , proteinToTranscriptFile = protein2transcriptFile
                            , tagContam = "Contam"
                            , remove = "sharedOnly")

# Check size after transcriptome-informed filtering
dim(incM_filt)

## -----------------------------------------------------------------------------
# Reduce incidence matrix size to accelerate downstream computation
incM_filt_reduced <- reduce_inc_matrix(incM_filt)

# Calculate the adjacency matrix describing protein-to-protein connections
adjM_filt <- get_adj_matrix(incM_filt_reduced)

# Generate a graph of protein-to-protein connections by shared peptides and
# calculate its CCs (i.e., sets of proteins connected by shared peptides
multProteinCC_filt <- get_cc(adjM_filt)

# Extract the list of vectors enumerating protein members in each CC 
cc.multProteins_filt <- multProteinCC_filt$ccs

# Calculate CCs size and % of single- vs multi-protein CCs obtained after
# transcriptome-informed filtering
CCstatsOut_filt <- cc_stats(incM = incM_filt
                                , cc.proteins = multProteinCC_filt$ccs
                                , reducedIncM = TRUE)

# Number of single-protein CCs
CCstatsOut_filt$N_singleProtCC

# Number of multi-protein CCs
CCstatsOut_filt$N_multiProtCC

# Total number of CCs
totCCs_filt <- CCstatsOut_filt$N_singleProtCC + CCstatsOut_filt$N_multiProtCC
totCCs_filt

# Percentage of single-protein CCs
PercSingleP_filt <- round(CCstatsOut_filt$N_singleProtCC / totCCs_filt * 100
                          , digits = 2)

# View table of CC size distribution
CCstatsOut_filt$NproteinsDistribution

# Plot CC size distribution
plot(factor(CCstatsOut_filt$NproteinsDistribution$N_proteins
       , levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", ">10"))
     , as.numeric(as.vector(CCstatsOut_filt$NproteinsDistribution$N_CC))
     , type = "s"
     , xlab = "N_proteins"
     , ylab = "N_CCs")

## -----------------------------------------------------------------------------
comp <- as.data.frame(cbind(as.character(as.vector(c("before_filter"
                                                  , "after_filter")))
                         , as.numeric(as.vector(c(PercSingleP
                                                  , PercSingleP_filt)))))
colnames(comp) <- c("Filter", "Perc_SingleP")

ggplot(data = comp
       , aes(x = as.factor(Filter), y = as.numeric(as.vector(Perc_SingleP)))) +
      geom_bar(stat = "identity") +
      theme_classic() +
      xlab("") +
      ylab("% single-prot CCs") +
      ylim(0, 100) +
      coord_flip() +
      geom_text(aes(label = as.numeric(as.vector(Perc_SingleP)))
                , hjust = 1.5, color = "white", size = 4)

## -----------------------------------------------------------------------------
old.par <- par(no.readonly = TRUE) # save default par values

ymax_before <- as.numeric(as.vector(CCstatsOut$NproteinsDistribution$N_CC))
ymax_after <- as.numeric(as.vector(CCstatsOut_filt$NproteinsDistribution$N_CC))

ymax <- max(max(ymax_before), max(ymax_after))

par(mfrow = c(1, 2))
plot(factor(CCstatsOut$NproteinsDistribution$N_proteins
        , levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", ">10"))
      , as.numeric(as.vector(CCstatsOut$NproteinsDistribution$N_CC))
      , type = "s"
      , xlab = "N_proteins"
      , ylab = "N_CCs"
      , ylim = c(0, ymax)
      , main = "before filtering")
plot(factor(CCstatsOut_filt$NproteinsDistribution$N_proteins
      , levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", ">10"))
      , as.numeric(as.vector(CCstatsOut_filt$NproteinsDistribution$N_CC))
      , type = "s"
      , xlab = "N_proteins"
      , ylab = "N_CCs"
      , ylim = c(0, ymax)
      , main = "after filtering")

par(old.par) # restore default par values

## -----------------------------------------------------------------------------
peptideStatsOut_filt <- peptide_stats(incM = incM_filt)

## -----------------------------------------------------------------------------
comp <- as.data.frame(cbind(
                      as.character(as.vector(c("before_filter"
                                            , "after_filter")))
                    , as.numeric(as.vector(c(peptideStatsOut$nbShared
                                        , peptideStatsOut_filt$nbShared)))))
colnames(comp) <- c("Filter", "Perc_sharedPeptides")

ggplot(data = comp
    , aes(x = as.factor(Filter)
          , y = as.numeric(as.vector(Perc_sharedPeptides)))) +
      geom_bar(stat = "identity") +
      theme_classic() +
      xlab("") +
      ylab("% shared peptides") +
      ylim(0, 100) +
      coord_flip() +
      geom_text(aes(label = as.numeric(as.vector(Perc_sharedPeptides)))
                , hjust = 1.5, color = "white", size = 4)

## ---- fig.height=14-----------------------------------------------------------
# Extract peptides and peptide-to-protein mappings for each CC after filtering
cc.peptides.incM_filt <- cc_composition(cc.multProteins_filt
                                            , incM = incM_filt)

# Generate a bipartite graph of the CC which contains the protein of interest,
# before and after transcriptome-informed filtering.
prot <- "ENSP261"
subgraphCC_beforeFilter <- plot_cc(prot = prot
                      , cc.proteins = cc.multProteins
                      , cc.subincM = cc.peptides.incM$cc.subincM
                      , incM = incM
                      , tagProt = "ENSP"
                      , tagContam = "Contam")

subgraphCC_afterFilter <- plot_cc(prot = prot
                     , cc.proteins = cc.multProteins_filt
                     , cc.subincM = cc.peptides.incM_filt$cc.subincM
                     , incM = incM_filt
                     , tagProt = "ENSP"
                     , tagContam = "Contam")

# Plot
old.par <- par(no.readonly = TRUE) # save default par values

par(mfrow = c(2, 1))
plot.igraph(subgraphCC_beforeFilter$g
            , layout = layout_as_bipartite
            , edge.width = 1
            , edge.arrow.width = 0.3
            , vertex.size = 35
            , edge.arrow.size = 0.5
            , vertex.size2 = 35
            , vertex.label.cex = 1
            , asp = 0.45
            , margin = -0.1) +
title(paste0("Protein "
             , prot
             , " in CC #"
             , subgraphCC_beforeFilter$cc_id
             , " before filtering")
      , line = -1)
plot.igraph(subgraphCC_afterFilter$g
            , layout = layout_as_bipartite
            , edge.width = 1
            , edge.arrow.width = 0.3
            , vertex.size = 35
            , edge.arrow.size = 0.5
            , vertex.size2 = 35
            , vertex.label.cex = 1
            , asp = 0.45
            , margin = -0.1) +
title(paste0("Protein "
             , prot, " in CC #"
             , subgraphCC_beforeFilter$cc_id
             , " after filtering")
      , line = -1)

par(old.par) # restore default par values

## -----------------------------------------------------------------------------
sessionInfo()

