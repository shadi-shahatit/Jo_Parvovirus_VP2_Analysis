#### Sanger sequence analysis MSA - v3 FINAL
#### Shadi Shahatit, MBV Lab, JUST, 2026
# libraries ---------------------------------------------------------------

sys_dir <- "/home/shadi/Desktop/ParvoProject/"

library(ips)
library(treeio)
library(tidytree)
library(pheatmap)
library(pegas)
library(openxlsx)
library(ggbreak)
library(systemPipeR)
library(Biostrings)
library(msa)
library(ape)
library(phangorn)
library(ggplot2)
library(ggrepel)
library(ggalt)
library(ggmsa)
library(pegas)
library(haplotypes)
library(ggtree)
library(viridis)
library(patchwork)
library(tidyverse)    # dplyr, tidyr, ggplot2, readr
library(stringr)
library(ggnewscale)
library(readxl)
library(writexl)

# MSA - DNA & Protein analysis --------------------------------------------------------

#### DNA and VP2 protein analysis

## load the files

dna_file <- file.path(sys_dir, "Parvo_sanger_full_seq","combined_sample_v3_vac_ref_ncbi_add_v142.fasta")
protein_file <- file.path(sys_dir,"Parvo_sanger_full_seq/Expasy_translate/Expasy_output/best_orf","combined_best_orf_sample_v3_vac_ref_ncbi_add_v142.fasta")

dna_sequences <- readDNAStringSet(dna_file)
protein_sequences <- readAAStringSet(protein_file)

dna_sequences                     
length(dna_sequences)              
width(dna_sequences)        
names(dna_sequences)       
dna_sequences[[1]]               
as.character(dna_sequences[1])     
letterFrequency(dna_sequences, letters = "GC", as.prob = TRUE)

length(protein_sequences)              

## naming issues

## clean names of FASTA HEADERS
## spaces/commas/slashes -> underscore
## accession_species_country for NCBI viruses
## outgroups do not exist here but will be replaced with add seq

addseq_map <- c(
  "EU659112.1" = "ref_FPV_oldest", 
  "KX434461.1" = "Feline_parvovirus_IZSSI_Italy",
  "KX900570.1" = "Feline_parvovirus_HH6_China",
  "OM638043.1" = "Feline_parvovirus_Egypt",
  "AY742933.1" = "recent_Canine_parvovirus_2a_NewZealand",
  "MG013488.1" = "Canine_parvovirus_2c_China",
  "MF510157.1" = "Canine_parvovirus_2c_Italy"
)

dna_names <- names(dna_sequences)
pipe_hdrs <- dna_names[grepl("\\|", dna_names)]
pipe_hdrs <- pipe_hdrs[!sub("\\s*\\|.*", "", trimws(pipe_hdrs)) %in% names(addseq_map)]
pipe_parts <- strsplit(pipe_hdrs, "\\|")
acc     <- trimws(sapply(pipe_parts, `[`, 1))
species <- trimws(sapply(pipe_parts, `[`, 3))
country <- trimws(sapply(pipe_parts, `[`, 4))
length(unique(acc))
unique(species)
unique(country)
ncbi_viruses_map <- setNames(paste0(gsub("_+", "_", gsub("[ ,/]+", "_", trimws(species))), "_",
                                    gsub("_+", "_", gsub("[ ,/]+", "_", trimws(country)))),
                             acc)

# addseq_map
# outgroup_map

rename_fasta_headers <- function(hdrs) {
  first_tok <- trimws(sub("\\s*\\|.*", "", hdrs))   # strip pipe metadata
  first_tok <- sub("\\s.*", "", first_tok)          # strip anything after first space
  out <- first_tok
  in_ncbi <- first_tok %in% names(ncbi_viruses_map)
  in_addseq_map <- first_tok %in% names(addseq_map)
  out[in_ncbi]     <- paste0(first_tok[in_ncbi], "_", ncbi_viruses_map[first_tok[in_ncbi]])
  out[in_addseq_map] <- paste0(first_tok[in_addseq_map], "_", addseq_map[first_tok[in_addseq_map]])
  out
}

dna_sequences_mod <- dna_sequences
protein_sequences_mod <- protein_sequences

names(dna_sequences_mod) <- rename_fasta_headers(dna_names)
names(protein_sequences_mod) <- rename_fasta_headers(names(protein_sequences))

names(dna_sequences)
names(dna_sequences_mod) 
names(protein_sequences) 
names(protein_sequences_mod) 
dna_sequences_mod                     
length(dna_sequences_mod)              
width(dna_sequences_mod)        
names(dna_sequences_mod)   
protein_sequences_mod                     
length(protein_sequences_mod)              
width(protein_sequences_mod)        
names(protein_sequences_mod) 

setequal(names(dna_sequences_mod), names(protein_sequences_mod))
identical(sort(names(dna_sequences_mod)), sort(names(protein_sequences_mod)))
# both should be T

## Alignment

## DNA

dna_alignment_Muscle <- msa(dna_sequences_mod, method = "Muscle")
dna_phy_Muscle <- phyDat(as.matrix(dna_alignment_Muscle), type = "DNA")

## we will stick with MUSCLE; medium-sized datasets, highly similar sequences, and protein-coding DNA

## trim MSA_Muscle

# M38246.1_ref_FPV VP2 DNA cord = 2655 - 4585

dna_alignment_Muscle_trimmed <- DNAStringSet(dna_alignment_Muscle)
ref_FPV_name <- grep("M38246.1_ref_FPV", names(dna_alignment_Muscle_trimmed), value = TRUE)
ref_FPV_seq <- dna_alignment_Muscle_trimmed[[ref_FPV_name]]
ref_FPV_chars <- strsplit(as.character(ref_FPV_seq), "")[[1]]
non_gap_positions <- which(ref_FPV_chars != "-")

start_col <- non_gap_positions[2655]
end_col <- non_gap_positions[4585]
# ref_FPV_seq[2655]
# ref_FPV_seq[4585]
# ref_FPV_seq[start_col]
# ref_FPV_seq[end_col]

dna_alignment_Muscle_trimmed <- subseq(dna_alignment_Muscle_trimmed, start = start_col, end = end_col)
width(dna_alignment_Muscle_trimmed)[1]
dna_phy_Muscle_trimmed <- phyDat(as.matrix(dna_alignment_Muscle_trimmed), type = "DNA")

## trim MSA_Muscle per FJ222821.1_ref_CPV2C

# FJ222821.1_ref_CPV2C VP2 DNA cord = 1 - 1755

# dna_alignment_Muscle_trimmed_CPV2C <- DNAStringSet(dna_alignment_Muscle_trimmed)
# ref_CPV2C_name <- grep("FJ222821.1_ref_CPV2C", names(dna_alignment_Muscle_trimmed_CPV2C), value = TRUE)
# ref_CPV2C_seq <- dna_alignment_Muscle_trimmed_CPV2C[[ref_CPV2C_name]]
# ref_CPV2C_chars <- strsplit(as.character(ref_CPV2C_seq), "")[[1]]
# non_gap_positions_CPV2C <- which(ref_CPV2C_chars != "-")
# 
# start_col <- non_gap_positions_CPV2C[1]
# end_col <- non_gap_positions_CPV2C[1755-3]
# ref_CPV2C_seq[1]
# ref_CPV2C_seq[1755]
# ref_CPV2C_seq[start_col]
# ref_CPV2C_seq[end_col]
# 
# dna_alignment_Muscle_trimmed_CPV2C <- subseq(dna_alignment_Muscle_trimmed_CPV2C, start = start_col, end = end_col)
# (width(dna_alignment_Muscle_trimmed_CPV2C)[1])/3
# dna_phy_Muscle_trimmed_CPV2C <- phyDat(as.matrix(dna_alignment_Muscle_trimmed_CPV2C), type = "DNA")
# 
# width(dna_alignment_Muscle_trimmed_CPV2C)/3

## trim MSA_Muscle per FJ222821.1_ref_CPV2C - v2

# CPV2C VP2 trim - reference-anchored column selection 4 143-taxaMSA 
dna_alignment_Muscle_trimmed_CPV2C <- DNAStringSet(dna_alignment_Muscle_trimmed)
ref_CPV2C_name <- grep("FJ222821.1_ref_CPV2C", names(dna_alignment_Muscle_trimmed_CPV2C), value = TRUE)
ref_CPV2C_chars <- strsplit(as.character(dna_alignment_Muscle_trimmed_CPV2C[[ref_CPV2C_name]]), "")[[1]]
non_gap_positions_CPV2C <- which(ref_CPV2C_chars != "-")
vp2_cols <- non_gap_positions_CPV2C[1:1752]          # exact ref non-gap columns, VP2 CDS length
length(vp2_cols); length(vp2_cols) / 3               

dna_aln_mat_trimmed_CPV2C <- as.matrix(dna_alignment_Muscle_trimmed_CPV2C)[, vp2_cols]
dna_alignment_Muscle_trimmed_CPV2C <- DNAStringSet(apply(dna_aln_mat_trimmed_CPV2C, 1, paste0, collapse = ""))
width(dna_alignment_Muscle_trimmed_CPV2C)[1] / 3
dna_phy_Muscle_trimmed_CPV2C <- phyDat(as.matrix(dna_alignment_Muscle_trimmed_CPV2C), type = "DNA")

width(dna_alignment_Muscle_trimmed_CPV2C)/3

writeXStringSet(DNAStringSet(dna_alignment_Muscle_trimmed_CPV2C),
                file.path(sys_dir, "Parvo_sanger_full_seq/aln/paper/Parvovirus_DNA_alignment_Muscle_trimmed_CPV2C_FINAL_SVRNA_v3_v142.fasta"))

## Alignment

## Protein 

protein_alignment_Muscle <- msa(protein_sequences_mod, method = "Muscle")
protein_phy_Muscle <- phyDat(as.matrix(protein_alignment_Muscle), type = "AA")

writeXStringSet(DNAStringSet(dna_alignment_Muscle_trimmed),
                file.path(sys_dir, "Parvo_sanger_full_seq/aln/paper/Parvovirus_DNA_alignment_Muscle_trimmed_FINAL_SVRNA_v3_v142.fasta"))
writeXStringSet(AAStringSet(protein_alignment_Muscle),
                file.path(sys_dir, "Parvo_sanger_full_seq/aln/paper/Parvovirus_Protein_alignment_Muscle_FINAL_SVRNA_v3_v142.fasta"))

## define color palette first for groups
group_colors <- c(
  "Jo_Feline"        = "firebrick4",
  "Regional_Feline"  = "firebrick1",
  "Jo_Lion"          = "maroon",
  "Jo_Canine"        = "darkblue",
  "Regional_Canine"  = "royalblue2",
  "Reference"        = "darkcyan",
  "Outgroup"         = "gray20",
  "Vaccine"          = "forestgreen",
  "Other"            = "purple4"
  )

# Neighbor Joining trees - DNA & Protein analysis --------------------------------------------------

## Neighbor Joining trees

## DNA

set.seed(32)
dna_NJ_Muscle_trimmed_JC <- NJ(dist.ml(dna_phy_Muscle_trimmed))
dna_NJ_Muscle_trimmed_F81 <- NJ(dist.ml(dna_phy_Muscle_trimmed,model = "F81"))
dna_NJ_Muscle_trimmed_K80 <- NJ(dist.dna(as.DNAbin(dna_phy_Muscle_trimmed), model = "K80"))

model_test <- modelTest(dna_phy_Muscle_trimmed, model = c("JC","F81","K80"), G = FALSE, I = FALSE)
# Model        df  logLik   AIC      BIC
# JC 281 -5260.684 11083.37 12647.5 
# F81 284 -5174.966 10917.93 12498.77 
# K80 282 -5195.703 10955.41 12525.11 

parsimony(dna_NJ_Muscle_trimmed_JC, dna_phy_Muscle_trimmed)
parsimony(dna_NJ_Muscle_trimmed_F81, dna_phy_Muscle_trimmed)
parsimony(dna_NJ_Muscle_trimmed_K80, dna_phy_Muscle_trimmed)

## F81 is best
dna_NJ_Muscle_trimmed <- dna_NJ_Muscle_trimmed_F81

## plot the trees

tip_df <- data.frame(
  label = dna_NJ_Muscle_trimmed$tip.label,
  group = ifelse(grepl("sample", dna_NJ_Muscle_trimmed$tip.label, ignore.case = TRUE) & grepl("feline", dna_NJ_Muscle_trimmed$tip.label, ignore.case = TRUE), "Jo_Feline",
                 ifelse(grepl("sample", dna_NJ_Muscle_trimmed$tip.label, ignore.case = TRUE) & grepl("canine", dna_NJ_Muscle_trimmed$tip.label, ignore.case = TRUE), "Jo_Canine",
                        ifelse(grepl("sample", dna_NJ_Muscle_trimmed$tip.label, ignore.case = TRUE) & grepl("lion",   dna_NJ_Muscle_trimmed$tip.label, ignore.case = TRUE), "Jo_Lion",
                               ifelse(grepl("NC_001510|NC_001718|NC_029797|NC_038544", dna_NJ_Muscle_trimmed$tip.label), "Outgroup",
                                      ifelse(grepl("_ref_", dna_NJ_Muscle_trimmed$tip.label), "Reference",
                                             ifelse(grepl("vaccine", dna_NJ_Muscle_trimmed$tip.label, ignore.case = TRUE), "Vaccine",
                                                    ifelse(grepl("feline", dna_NJ_Muscle_trimmed$tip.label, ignore.case = TRUE), "Regional_Feline",
                                                           ifelse(grepl("canine", dna_NJ_Muscle_trimmed$tip.label, ignore.case = TRUE), "Regional_Canine",
                                                                  "Other"))))))))
)

## filtering and rooting v2

## remove zero branch
dna_NJ_Muscle_trimmed_collapse <- di2multi(dna_NJ_Muscle_trimmed, tol = 1e-8)
length(dna_NJ_Muscle_trimmed$tip.label)
length(dna_NJ_Muscle_trimmed_collapse$tip.label)
## none were removed; use manual filtering

dna_NJ_Muscle_trimmed_midpoint <- midpoint(dna_NJ_Muscle_trimmed)
dna_NJ_Muscle_trimmed_rooted_oldestFPV <- root(dna_NJ_Muscle_trimmed, outgroup = "EU659112.1_ref_FPV_oldest", resolve.root = TRUE)

parsimony(dna_NJ_Muscle_trimmed, dna_phy_Muscle_trimmed)
parsimony(dna_NJ_Muscle_trimmed_midpoint, dna_phy_Muscle_trimmed)
parsimony(dna_NJ_Muscle_trimmed_rooted_oldestFPV, dna_phy_Muscle_trimmed)

ggtreeplot_dna_NJ_Muscle_trimmed <- ggtree(dna_NJ_Muscle_trimmed, layout = "rectangular",
                                           # branch.length = "none",
                                           open.angle = 10, size = 0.3) %<+% tip_df +
  geom_tiplab(aes(color = group), size = 2, offset = 0.002, align = TRUE, linesize = 0.1) +
  geom_tippoint(aes(color = group), size = 1.2) +
  scale_color_manual(values = group_colors) +
  theme_tree2() +
  labs(title = paste0("Parvovirus DNA - NJ - ", deparse(substitute(dna_NJ_Muscle_trimmed))), color = "Group") +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 8)) +
  coord_cartesian(clip = "off") +                  
  hexpand(0.5, direction = 1)

ggtreeplot_dna_NJ_Muscle_trimmed_midpoint <- ggtree(dna_NJ_Muscle_trimmed_midpoint, layout = "rectangular",
                                                    # branch.length = "none",
                                                    open.angle = 10, size = 0.3) %<+% tip_df +
  geom_tiplab(aes(color = group), size = 2, offset = 0.002, align = TRUE, linesize = 0.1) +
  geom_tippoint(aes(color = group), size = 1.2) +
  scale_color_manual(values = group_colors) +
  theme_tree2() +
  labs(title = paste0("Parvovirus DNA - NJ - ", deparse(substitute(dna_NJ_Muscle_trimmed_midpoint))), color = "Group") +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 8)) +
  coord_cartesian(clip = "off") +                  
  hexpand(0.5, direction = 1)

ggtreeplot_dna_NJ_Muscle_trimmed_rooted_oldestFPV <- ggtree(dna_NJ_Muscle_trimmed_rooted_oldestFPV, layout = "rectangular",
                                                                 # branch.length = "none",
                                                                 open.angle = 10, size = 0.3) %<+% tip_df +
  geom_tiplab(aes(color = group), size = 2, offset = 0.002, align = TRUE, linesize = 0.1) +
  geom_tippoint(aes(color = group), size = 1.2) +
  scale_color_manual(values = group_colors) +
  theme_tree2() +
  labs(title = paste0("Parvovirus DNA - NJ - ", deparse(substitute(dna_NJ_Muscle_trimmed_rooted_oldestFPV))), color = "Group") +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 8)) +
  coord_cartesian(clip = "off") +                  
  hexpand(0.5, direction = 1)

# ## more filtering
# 
# plot_data <- ggtreeplot_dna_NJ_Muscle_trimmed_rooted_oldestFPV$data
# plot_data <- plot_data[plot_data$isTip, ]
# tip_order <- plot_data$label[order(plot_data$y, decreasing = TRUE)]
# 
# tips_to_keep <- c(
#   "MZ056887.1_Canine_parvovirus_2_Egypt",
#   # "MZ056881.1_Canine_parvovirus_2_Egypt",
#   # "MZ056882.1_Canine_parvovirus_2_Egypt",
#   # "MZ056890.1_Canine_parvovirus_2_Egypt",
#   # "MZ056888.1_Canine_parvovirus_2_Egypt",
#   "MW539053.1_Canine_parvovirus_2b_Turkey",
#   # "OM721656.1_Canine_parvovirus_2b_Turkey",
#   # "OM721655.1_Canine_parvovirus_2b_Turkey",
#   "OM100702.1_Canine_parvovirus_2_Egypt",
#   # "OM100700.1_Canine_parvovirus_2_Egypt",
#   # "MZ056892.1_Canine_parvovirus_2_Egypt",
#   # "MZ056883.1_Canine_parvovirus_2_Egypt",
#   # "MZ056891.1_Canine_parvovirus_2_Egypt",
#   "OL330979.1_Canine_parvovirus_2_Iran",
#   # "OL330978.1_Canine_parvovirus_2_Iran",
#   "MW653250.1_Canine_parvovirus_2a_Iran",
#   # "OL330980.1_Canine_parvovirus_2_Iran",
#   # "MW653249.1_Canine_parvovirus_2a_Iran",
#   "OR667800.1_Canine_parvovirus_2_Iraq",
#   # "OR451707.1_Canine_parvovirus_2_Iraq",
#   # "OR667801.1_Canine_parvovirus_2_Iraq",
#   # "OR667804.1_Canine_parvovirus_2_Iraq",
#   # "OR667803.1_Canine_parvovirus_2_Iraq",
#   "OR667802.1_Canine_parvovirus_2_Iraq",
#   "sample18_canine",
#   "sample8_canine_79_rev",
#   "sample68_canine",
#   "Sample47_canine",
#   "OM100699.1_Canine_parvovirus_2_Egypt",
#   "OM100696.1_Canine_parvovirus_2_Egypt",
#   "Sample77_canine",
#   "sample76_canine",
#   "sample57_canine",
#   "sample51_canine",
#   "sample72_canine",
#   "OM100701.1_Canine_parvovirus_2_Egypt",
#   "OM100698.1_Canine_parvovirus_2_Egypt",
#   "MF510157.1_Canine_parvovirus_2c_Italy",
#   "MG013488.1_Canine_parvovirus_2c_China",
#   # "MZ056886.1_Canine_parvovirus_2_Egypt",
#   # "MZ056884.1_Canine_parvovirus_2_Egypt",
#   # "MZ056889.1_Canine_parvovirus_2_Egypt",
#   # "MZ056885.1_Canine_parvovirus_2_Egypt",
#   # "OM100697.1_Canine_parvovirus_2_Egypt",
#   "MW653248.1_Canine_parvovirus_2a_Iran",
#   # "KX268117.1_Canine_parvovirus_2_Turkey",
#   # "KX268105.1_Canine_parvovirus_2_Turkey",
#   "OL330977.1_Canine_parvovirus_2_Iran",
#   # "KX268118.1_Canine_parvovirus_2_Turkey",
#   # "KX268106.1_Canine_parvovirus_2_Turkey",
#   # "KX268107.1_Canine_parvovirus_2_Turkey",
#   "MW653251.1_Canine_parvovirus_2b_Iran",
#   "MW653253.1_Canine_parvovirus_2c_Iran",
#   "KX268109.1_Canine_parvovirus_2_Turkey",
#   "FJ222821.1_ref_CPV2C",
#   "KX268114.1_Canine_parvovirus_2_Turkey",
#   # "KX268110.1_Canine_parvovirus_2_Turkey",
#   # "KX268108.1_Canine_parvovirus_2_Turkey",
#   "MW653256.1_Canine_parvovirus_2a_Iran",
#   "KX268115.1_Canine_parvovirus_2_Turkey",
#   "AY742933.1_recent_Canine_parvovirus_2a_NewZealand",
#   # "KX268116.1_Canine_parvovirus_2_Turkey",
#   # "KX268111.1_Canine_parvovirus_2_Turkey",
#   # "KX268113.1_Canine_parvovirus_2_Turkey",
#   # "KX268112.1_Canine_parvovirus_2_Turkey",
#   "MW653252.1_Canine_parvovirus_2b_Iran",
#   "M74849.1_ref_CPV2B",
#   "M24003.1_ref_CPV2A",
#   "vaccineCanine_Zoetis",
#   "M38245.1_ref_CPV2",
#   "vaccineCanine_MSD",
#   "KX268119.1_Canine_parvovirus_2_Turkey",
#   "M38246.1_ref_FPV",
#   "KP081409.1_Feline_parvovirus_Iran",
#   "vaccineFeline_Zoetis",
#   "vaccineFeline_MSD",
#   "PP663067.2_Feline_parvovirus_Egypt",
#   # "PP663062.2_Feline_parvovirus_Egypt",
#   # "PV808486.1_Feline_parvovirus_Egypt",
#   # "PV808493.1_Feline_parvovirus_Egypt",
#   # "PV808496.1_Feline_parvovirus_Egypt",
#   # "PV521956.1_Feline_parvovirus_Egypt",
#   # "PP663061.2_Feline_parvovirus_Egypt",
#   # "PV808491.1_Feline_parvovirus_Egypt",
#   # "PV521951.1_Feline_parvovirus_Egypt",
#   # "PV808488.1_Feline_parvovirus_Egypt",
#   # "PV521944.1_Feline_parvovirus_Egypt",
#   # "PV521945.1_Feline_parvovirus_Egypt",
#   # "PV521957.1_Feline_parvovirus_Egypt",
#   # "PV521958.1_Feline_parvovirus_Egypt",
#   # "PV521959.1_Feline_parvovirus_Egypt",
#   # "PV521965.1_Feline_parvovirus_Egypt",
#   # "PV521946.1_Feline_parvovirus_Egypt",
#   # "PP663047.2_Feline_parvovirus_Egypt",
#   # "PV808492.1_Feline_parvovirus_Egypt",
#   # "PV521952.1_Feline_parvovirus_Egypt",
#   # "PV521943.1_Feline_parvovirus_Egypt",
#   # "PP663042.2_Feline_parvovirus_Egypt",
#   # "PV521950.1_Feline_parvovirus_Egypt",
#   "PP663052.2_Feline_parvovirus_Egypt",
#   "sample67_feline",
#   "Sample12_Feline",
#   "PV521947.1_Feline_parvovirus_Egypt",
#   "PP663048.2_Feline_parvovirus_Egypt",
#   "sample_lion",
#   "Sample69_Feline",
#   "FPV_JO24_Sample_11",
#   "Sample_44_feline",
#   "Sample70_Feline",
#   "sample7_feline",
#   "Sample25_Feline",
#   "sample38_feline",
#   "sample27_feline",
#   "PP663041.2_Feline_parvovirus_Egypt",
#   # "PV521961.1_Feline_parvovirus_Egypt",
#   # "PP663072.2_Feline_parvovirus_Egypt",
#   # "PP663045.2_Feline_parvovirus_Egypt",
#   # "PV521949.1_Feline_parvovirus_Egypt",
#   # "PP663073.2_Feline_parvovirus_Egypt",
#   # "PP663070.2_Feline_parvovirus_Egypt",
#   # "PP663050.2_Feline_parvovirus_Egypt",
#   # "PV521963.1_Feline_parvovirus_Egypt",
#   # "PV808495.1_Feline_parvovirus_Egypt",
#   # "PV521960.1_Feline_parvovirus_Egypt",
#   # "PV521942.1_Feline_parvovirus_Egypt",
#   # "PV808490.1_Feline_parvovirus_Egypt",
#   # "PV808489.1_Feline_parvovirus_Egypt",
#   # "PV521962.1_Feline_parvovirus_Egypt",
#   # "PV521953.1_Feline_parvovirus_Egypt",
#   # "PV521954.1_Feline_parvovirus_Egypt",
#   # "PP663055.2_Feline_parvovirus_Egypt",
#   # "PP663056.2_Feline_parvovirus_Egypt",
#   "OM638043.1_Feline_parvovirus_Egypt",
#   # "PV521948.1_Feline_parvovirus_Egypt",
#   # "PV521964.1_Feline_parvovirus_Egypt",
#   # "PP663074.2_Feline_parvovirus_Egypt",
#   # "PP663071.2_Feline_parvovirus_Egypt",
#   # "PP663044.2_Feline_parvovirus_Egypt",
#   # "PV521955.1_Feline_parvovirus_Egypt",
#   "PP663059.2_Feline_parvovirus_Egypt",
#   "KX900570.1_Feline_parvovirus_HH6_China",
#   "KX434461.1_Feline_parvovirus_IZSSI_Italy",
#   "EU659112.1_ref_FPV_oldest"
#   )
# 
# tips_to_drop <- setdiff(dna_NJ_Muscle_trimmed$tip.label, tips_to_keep)
# dna_NJ_Muscle_trimmed_filtered <- drop.tip(dna_NJ_Muscle_trimmed, tips_to_drop)
# dna_NJ_Muscle_trimmed_filtered_midpoint <- midpoint(dna_NJ_Muscle_trimmed_filtered)
# dna_NJ_Muscle_trimmed_filtered_rooted_oldestFPV <- root(dna_NJ_Muscle_trimmed_filtered, outgroup = "EU659112.1_ref_FPV_oldest", resolve.root = TRUE)
# 
# length(dna_NJ_Muscle_trimmed$tip.label)
# length(dna_NJ_Muscle_trimmed_filtered_midpoint$tip.label)
# length(dna_NJ_Muscle_trimmed_filtered_rooted_oldestFPV$tip.label)
# 
# ggtreeplot_dna_NJ_Muscle_trimmed_filtered <- ggtree(dna_NJ_Muscle_trimmed_filtered, layout = "rectangular",
#                                                     # branch.length = "none",
#                                                     open.angle = 10, size = 0.3) %<+% tip_df +
#   geom_tiplab(aes(color = group), size = 2, offset = 0.002, align = TRUE, linesize = 0.1) +
#   geom_tippoint(aes(color = group), size = 1.2) +
#   scale_color_manual(values = group_colors) +
#   theme_tree2() +
#   labs(title = paste0("Parvovirus DNA - NJ - ", deparse(substitute(dna_NJ_Muscle_trimmed_filtered))), color = "Group") +
#   theme(legend.position = "bottom",
#         plot.title = element_text(size = 14, face = "bold"),
#         legend.text = element_text(size = 8)) +
#   coord_cartesian(clip = "off") +                  
#   hexpand(0.5, direction = 1)
# 
# ggtreeplot_dna_NJ_Muscle_trimmed_filtered_midpoint <- ggtree(dna_NJ_Muscle_trimmed_filtered_midpoint, layout = "rectangular",
#                                                              # branch.length = "none",
#                                                              open.angle = 10, size = 0.3) %<+% tip_df +
#   geom_tiplab(aes(color = group), size = 2, offset = 0.002, align = TRUE, linesize = 0.1) +
#   geom_tippoint(aes(color = group), size = 1.2) +
#   scale_color_manual(values = group_colors) +
#   theme_tree2() +
#   labs(title = paste0("Parvovirus DNA - NJ - ", deparse(substitute(dna_NJ_Muscle_trimmed_filtered_midpoint))), color = "Group") +
#   theme(legend.position = "bottom",
#         plot.title = element_text(size = 14, face = "bold"),
#         legend.text = element_text(size = 8)) +
#   coord_cartesian(clip = "off") +                  
#   hexpand(0.5, direction = 1)
# 
# ggtreeplot_dna_NJ_Muscle_trimmed_filtered_rooted_oldestFPV <- ggtree(dna_NJ_Muscle_trimmed_filtered_rooted_oldestFPV, layout = "rectangular",
#                                                              # branch.length = "none",
#                                                              open.angle = 10, size = 0.3) %<+% tip_df +
#   geom_tiplab(aes(color = group), size = 2, offset = 0.002, align = TRUE, linesize = 0.1) +
#   geom_tippoint(aes(color = group), size = 1.2) +
#   scale_color_manual(values = group_colors) +
#   theme_tree2() +
#   labs(title = paste0("Parvovirus DNA - NJ - ", deparse(substitute(dna_NJ_Muscle_trimmed_filtered_rooted_oldestFPV))), color = "Group") +
#   theme(legend.position = "bottom",
#         plot.title = element_text(size = 14, face = "bold"),
#         legend.text = element_text(size = 8)) +
#   coord_cartesian(clip = "off") +                  
#   hexpand(0.5, direction = 1)

## Protein

## Neighbor Joining trees

set.seed(32)
protein_NJ_Muscle <- NJ(dist.ml(protein_phy_Muscle))
protein_NJ_Muscle_WAG <- NJ(dist.ml(protein_phy_Muscle,model = "WAG"))
protein_NJ_Muscle_JTT <- NJ(dist.ml(protein_phy_Muscle,model = "JTT"))
protein_NJ_Muscle_LG <- NJ(dist.ml(protein_phy_Muscle,model = "LG"))

model_test <- modelTest(protein_phy_Muscle, model = c("WAG","JTT","LG"), G = FALSE, I = FALSE)
# Model        df  logLik   AIC      BIC
# WAG 281 -2924.178 6410.356 7697.905 
# JTT 281 -2909.251 6380.503 7668.052 
# LG 281 -2948.542 6459.084 7746.633 

parsimony(protein_NJ_Muscle_WAG, protein_phy_Muscle)
parsimony(protein_NJ_Muscle_JTT, protein_phy_Muscle)
parsimony(protein_NJ_Muscle_LG,  protein_phy_Muscle)

## JTT is best
protein_NJ_Muscle <- protein_NJ_Muscle_JTT

## plot the trees

## filtering and rooting v2

## remove zero branch
protein_NJ_Muscle_collapse <- di2multi(protein_NJ_Muscle, tol = 1e-8)
length(protein_NJ_Muscle$tip.label)
length(protein_NJ_Muscle_collapse$tip.label)
## none were removed; use manual filtering

protein_NJ_Muscle_midpoint <- midpoint(protein_NJ_Muscle)
protein_NJ_Muscle_rooted_oldestFPV <- root(protein_NJ_Muscle, outgroup = "EU659112.1_ref_FPV_oldest", resolve.root = TRUE)

ggtreeplot_protein_NJ_Muscle <- ggtree(protein_NJ_Muscle, layout = "rectangular",
                                       # branch.length = "none",
                                       open.angle = 10, size = 0.3) %<+% tip_df +
  geom_tiplab(aes(color = group), size = 2, offset = 0.002, align = TRUE, linesize = 0.1) +
  geom_tippoint(aes(color = group), size = 1.2) +
  scale_color_manual(values = group_colors) +
  theme_tree2() +
  labs(title = paste0("Parvovirus Protein - NJ - ", deparse(substitute(protein_NJ_Muscle))), color = "Group") +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 8)) +
  coord_cartesian(clip = "off") +                  
  hexpand(0.5, direction = 1)

ggtreeplot_protein_NJ_Muscle_midpoint <- ggtree(protein_NJ_Muscle_midpoint, layout = "rectangular",
                                                # branch.length = "none",
                                                open.angle = 10, size = 0.3) %<+% tip_df +
  geom_tiplab(aes(color = group), size = 2, offset = 0.002, align = TRUE, linesize = 0.1) +
  geom_tippoint(aes(color = group), size = 1.2) +
  scale_color_manual(values = group_colors) +
  theme_tree2() +
  labs(title = paste0("Parvovirus Protein - NJ - ", deparse(substitute(protein_NJ_Muscle_midpoint))), color = "Group") +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 8)) +
  coord_cartesian(clip = "off") +                  
  hexpand(0.5, direction = 1)

ggtreeplot_protein_NJ_Muscle_rooted_oldestFPV <- ggtree(protein_NJ_Muscle_rooted_oldestFPV, layout = "rectangular",
                                                # branch.length = "none",
                                                open.angle = 10, size = 0.3) %<+% tip_df +
  geom_tiplab(aes(color = group), size = 2, offset = 0.002, align = TRUE, linesize = 0.1) +
  geom_tippoint(aes(color = group), size = 1.2) +
  scale_color_manual(values = group_colors) +
  theme_tree2() +
  labs(title = paste0("Parvovirus Protein - NJ - ", deparse(substitute(protein_NJ_Muscle_rooted_oldestFPV))), color = "Group") +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 8)) +
  coord_cartesian(clip = "off") +                  
  hexpand(0.5, direction = 1)

# ## more filtering
# 
# protein_NJ_Muscle_filtered <- drop.tip(protein_NJ_Muscle, tips_to_drop)
# protein_NJ_Muscle_filtered_midpoint <- midpoint(protein_NJ_Muscle_filtered)
# protein_NJ_Muscle_filtered_rooted_oldestFPV <- root(protein_NJ_Muscle_filtered, outgroup = "EU659112.1_ref_FPV_oldest", resolve.root = TRUE)
# 
# length(protein_NJ_Muscle$tip.label)
# length(protein_NJ_Muscle_filtered_midpoint$tip.label)
# length(protein_NJ_Muscle_filtered_midpoint$tip.label)
# 
# ggtreeplot_protein_NJ_Muscle_filtered <- ggtree(protein_NJ_Muscle_filtered, layout = "rectangular",
#                                                          open.angle = 10, size = 0.3) %<+% tip_df +
#   geom_tiplab(aes(color = group), size = 3, offset = 0.002, align = TRUE, linesize = 0.1) +
#   geom_tippoint(aes(color = group), size = 1.2) +
#   scale_color_manual(values = group_colors) +
#   theme_tree2() +
#   labs(title = paste0("Parvovirus Protein - NJ - ", deparse(substitute(protein_NJ_Muscle_filtered))), color = "Group") +
#   theme(legend.position = "bottom",
#         plot.title = element_text(size = 14, face = "bold"),
#         legend.text = element_text(size = 8)) +
#   coord_cartesian(clip = "off") +
#   hexpand(0.5, direction = 1)
# 
# ggtreeplot_protein_NJ_Muscle_filtered_midpoint <- ggtree(protein_NJ_Muscle_filtered_midpoint, layout = "rectangular",
#                                                          open.angle = 10, size = 0.3) %<+% tip_df +
#   geom_tiplab(aes(color = group), size = 3, offset = 0.002, align = TRUE, linesize = 0.1) +
#   geom_tippoint(aes(color = group), size = 1.2) +
#   scale_color_manual(values = group_colors) +
#   theme_tree2() +
#   labs(title = paste0("Parvovirus Protein - NJ - ", deparse(substitute(protein_NJ_Muscle_filtered_midpoint))), color = "Group") +
#   theme(legend.position = "bottom",
#         plot.title = element_text(size = 14, face = "bold"),
#         legend.text = element_text(size = 8)) +
#   coord_cartesian(clip = "off") +
#   hexpand(0.5, direction = 1)
# 
# ggtreeplot_protein_NJ_Muscle_filtered_rooted_oldestFPV <- ggtree(protein_NJ_Muscle_filtered_rooted_oldestFPV, layout = "rectangular",
#                                                          open.angle = 10, size = 0.3) %<+% tip_df +
#   geom_tiplab(aes(color = group), size = 3, offset = 0.002, align = TRUE, linesize = 0.1) +
#   geom_tippoint(aes(color = group), size = 1.2) +
#   scale_color_manual(values = group_colors) +
#   theme_tree2() +
#   labs(title = paste0("Parvovirus Protein - NJ - ", deparse(substitute(protein_NJ_Muscle_filtered_rooted_oldestFPV))), color = "Group") +
#   theme(legend.position = "bottom",
#         plot.title = element_text(size = 14, face = "bold"),
#         legend.text = element_text(size = 8)) +
#   coord_cartesian(clip = "off") +
#   hexpand(0.5, direction = 1)

# Maximum Likelihood trees - DNA & Protein analysis -----------------------

## Maximum Likelihood trees from iqtree 

## DNA

dna_ML_Muscle_trimmed_tree <- read.tree(file.path(sys_dir,"Parvo_sanger_full_seq/iqtree_output/paper/Maxll_DNA_aln_Muscle_trimmed_FINAL_SVRNA_v3_v142.contree"))

tip_df_dna <- data.frame(
  label = dna_ML_Muscle_trimmed_tree$tip.label,
  group = ifelse(grepl("sample", dna_ML_Muscle_trimmed_tree$tip.label, ignore.case = TRUE) & grepl("feline", dna_ML_Muscle_trimmed_tree$tip.label, ignore.case = TRUE), "Jo_Feline",
                 ifelse(grepl("sample", dna_ML_Muscle_trimmed_tree$tip.label, ignore.case = TRUE) & grepl("canine", dna_ML_Muscle_trimmed_tree$tip.label, ignore.case = TRUE), "Jo_Canine",
                        ifelse(grepl("sample", dna_ML_Muscle_trimmed_tree$tip.label, ignore.case = TRUE) & grepl("lion",   dna_ML_Muscle_trimmed_tree$tip.label, ignore.case = TRUE), "Jo_Lion",
                               ifelse(grepl("NC_001510|NC_001718|NC_029797|NC_038544", dna_ML_Muscle_trimmed_tree$tip.label), "Outgroup",
                                      ifelse(grepl("_ref_", dna_ML_Muscle_trimmed_tree$tip.label), "Reference",
                                             ifelse(grepl("vaccine", dna_ML_Muscle_trimmed_tree$tip.label, ignore.case = TRUE), "Vaccine",
                                                    ifelse(grepl("feline", dna_ML_Muscle_trimmed_tree$tip.label, ignore.case = TRUE), "Regional_Feline",
                                                           ifelse(grepl("canine", dna_ML_Muscle_trimmed_tree$tip.label, ignore.case = TRUE), "Regional_Canine",
                                                                  "Other"))))))))
)

## rooting v2

dna_ML_Muscle_trimmed_tree_midpoint <- midpoint(dna_ML_Muscle_trimmed_tree)
dna_ML_Muscle_trimmed_tree_rooted_oldestFPV <- root(dna_ML_Muscle_trimmed_tree, outgroup = "EU659112.1_ref_FPV_oldest", resolve.root = TRUE)

parsimony(dna_ML_Muscle_trimmed_tree, dna_phy_Muscle_trimmed)
parsimony(dna_ML_Muscle_trimmed_tree_midpoint, dna_phy_Muscle_trimmed)
parsimony(dna_ML_Muscle_trimmed_tree_rooted_oldestFPV, dna_phy_Muscle_trimmed)

parsimony(dna_NJ_Muscle_trimmed_midpoint, dna_phy_Muscle_trimmed)
parsimony(dna_NJ_Muscle_trimmed_rooted_oldestFPV, dna_phy_Muscle_trimmed)

ggtreeplot_dna_ML_Muscle_trimmed_tree <- ggtree(dna_ML_Muscle_trimmed_tree, layout = "rectangular", size = 0.3) %<+% tip_df_dna +
  geom_tiplab(aes(color = group), size = 2, offset = 0.002, align = TRUE, linesize = 0.1) +
  geom_tippoint(aes(color = group), size = 1.2) +
  geom_text2(aes(label = label, subset = !isTip), size = 1.8, hjust = 1.2, vjust = -0.4) +
  scale_color_manual(values = group_colors) +
  theme_tree2() +
  labs(title = paste0("Parvovirus DNA - ML", deparse(substitute(dna_ML_Muscle_trimmed_tree))), color = "Group") +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 8)) +
  coord_cartesian(clip = "off") +
  hexpand(0.5, direction = 1)

ggtreeplot_dna_ML_Muscle_trimmed_tree_midpoint <- ggtree(dna_ML_Muscle_trimmed_tree_midpoint, layout = "rectangular", size = 0.3) %<+% tip_df_dna +
  geom_tiplab(aes(color = group), size = 2, offset = 0.002, align = TRUE, linesize = 0.1) +
  geom_tippoint(aes(color = group), size = 1.2) +
  geom_text2(aes(label = label, subset = !isTip), size = 1.8, hjust = 1.2, vjust = -0.4) +
  scale_color_manual(values = group_colors) +
  theme_tree2() +
  labs(title = paste0("Parvovirus DNA - ML", deparse(substitute(dna_ML_Muscle_trimmed_tree_midpoint))), color = "Group") +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 8)) +
  coord_cartesian(clip = "off") +
  hexpand(0.5, direction = 1)

ggtreeplot_dna_ML_Muscle_trimmed_tree_rooted_oldestFPV <- ggtree(dna_ML_Muscle_trimmed_tree_rooted_oldestFPV, layout = "rectangular", size = 0.3) %<+% tip_df_dna +
  geom_tiplab(aes(color = group), size = 2, offset = 0.002, align = TRUE, linesize = 0.1) +
  geom_tippoint(aes(color = group), size = 1.2) +
  geom_text2(aes(label = label, subset = !isTip), size = 1.8, hjust = 1.2, vjust = -0.4) +
  scale_color_manual(values = group_colors) +
  theme_tree2() +
  labs(
    # title = paste0("Parvovirus DNA - ML", deparse(substitute(dna_ML_Muscle_trimmed_tree_rooted_oldestFPV))),
    color = "Group") +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 8)) +
  coord_cartesian(clip = "off") +
  hexpand(0.5, direction = 1)

## more filtering

plot_data <- ggtreeplot_dna_ML_Muscle_trimmed_tree_rooted_oldestFPV$data
plot_data <- plot_data[plot_data$isTip, ]
tip_order <- plot_data$label[order(plot_data$y, decreasing = TRUE)]

## FINAL_SVRNA_v3

# tips_to_keep <- c(
#   "OM100701.1_Canine_parvovirus_2_Egypt",
#   # "OM100698.1_Canine_parvovirus_2_Egypt",
#   # "OM100696.1_Canine_parvovirus_2_Egypt",
#   # "OM100697.1_Canine_parvovirus_2_Egypt",
#   "OM100699.1_Canine_parvovirus_2_Egypt",
#   "MG013488.1_Canine_parvovirus_2c_China",
#   "MF510157.1_Canine_parvovirus_2c_Italy",
#   "Sample57_Canine",
#   "Sample68_Canine",
#   "Sample8_Canine",
#   "Sample18_Canine",
#   "Sample77_Canine",
#   "Sample76_Canine",
#   "Sample47_Canine",
#   "MZ056886.1_Canine_parvovirus_2_Egypt",
#   "MZ056884.1_Canine_parvovirus_2_Egypt",
#   "MZ056889.1_Canine_parvovirus_2_Egypt",
#   "sample51_canine",
#   "Sample72_Canine",
#   "MZ056885.1_Canine_parvovirus_2_Egypt",
#   "OR667801.1_Canine_parvovirus_2_Iraq",
#   # "OR667800.1_Canine_parvovirus_2_Iraq",
#   # "OR667802.1_Canine_parvovirus_2_Iraq",
#   "OR451707.1_Canine_parvovirus_2_Iraq",
#   # "OR667804.1_Canine_parvovirus_2_Iraq",
#   # "OR667803.1_Canine_parvovirus_2_Iraq",
#   # "MZ056890.1_Canine_parvovirus_2_Egypt",
#   # "MZ056882.1_Canine_parvovirus_2_Egypt",
#   # "MZ056888.1_Canine_parvovirus_2_Egypt",
#   # "MZ056887.1_Canine_parvovirus_2_Egypt",
#   # "MZ056881.1_Canine_parvovirus_2_Egypt",
#   "MW539053.1_Canine_parvovirus_2b_Turkey",
#   # "OM721656.1_Canine_parvovirus_2b_Turkey",
#   # "OM721655.1_Canine_parvovirus_2b_Turkey",
#   "MZ056892.1_Canine_parvovirus_2_Egypt",
#   # "MZ056883.1_Canine_parvovirus_2_Egypt",
#   # "MZ056891.1_Canine_parvovirus_2_Egypt",
#   "OL330979.1_Canine_parvovirus_2_Iran",
#   # "OM100702.1_Canine_parvovirus_2_Egypt",
#   # "OM100700.1_Canine_parvovirus_2_Egypt",
#   "MW653250.1_Canine_parvovirus_2a_Iran",
#   # "OL330980.1_Canine_parvovirus_2_Iran",
#   # "OL330978.1_Canine_parvovirus_2_Iran",
#   # "MW653249.1_Canine_parvovirus_2a_Iran",
#   # "MW653248.1_Canine_parvovirus_2a_Iran",
#   # "OL330977.1_Canine_parvovirus_2_Iran",
#   "KX268113.1_Canine_parvovirus_2_Turkey",
#   # "KX268112.1_Canine_parvovirus_2_Turkey",
#   "MW653252.1_Canine_parvovirus_2b_Iran",
#   "MW653256.1_Canine_parvovirus_2a_Iran",
#   # "KX268115.1_Canine_parvovirus_2_Turkey",
#   "AY742933.1_recent_Canine_parvovirus_2a_NewZealand",
#   # "KX268116.1_Canine_parvovirus_2_Turkey",
#   # "KX268111.1_Canine_parvovirus_2_Turkey",
#   # "KX268117.1_Canine_parvovirus_2_Turkey",
#   # "KX268105.1_Canine_parvovirus_2_Turkey",
#   # "KX268107.1_Canine_parvovirus_2_Turkey",
#   # "KX268118.1_Canine_parvovirus_2_Turkey",
#   # "KX268106.1_Canine_parvovirus_2_Turkey",
#   # "KX268114.1_Canine_parvovirus_2_Turkey",
#   # "KX268110.1_Canine_parvovirus_2_Turkey",
#   # "KX268108.1_Canine_parvovirus_2_Turkey",
#   "MW653251.1_Canine_parvovirus_2b_Iran",
#   "MW653253.1_Canine_parvovirus_2c_Iran",
#   "KX268109.1_Canine_parvovirus_2_Turkey",
#   "FJ222821.1_ref_CPV2C",
#   "M74849.1_ref_CPV2B",
#   "M24003.1_ref_CPV2A",
#   "M38245.1_ref_CPV2",
#   "vaccineCanine_MSD",
#   "KX268119.1_Canine_parvovirus_2_Turkey",
#   "vaccineCanine_Zoetis",
#   "PV808492.1_Feline_parvovirus_Egypt",
#   # "PV521952.1_Feline_parvovirus_Egypt",
#   # "PV808488.1_Feline_parvovirus_Egypt",
#   # "PV521944.1_Feline_parvovirus_Egypt",
#   # "PV521943.1_Feline_parvovirus_Egypt",
#   # "PP663042.2_Feline_parvovirus_Egypt",
#   # "PV521958.1_Feline_parvovirus_Egypt",
#   # "PV521945.1_Feline_parvovirus_Egypt",
#   # "PV521965.1_Feline_parvovirus_Egypt",
#   # "PV521959.1_Feline_parvovirus_Egypt",
#   "PV521957.1_Feline_parvovirus_Egypt",
#   # "PV808491.1_Feline_parvovirus_Egypt",
#   # "PV521951.1_Feline_parvovirus_Egypt",
#   # "PV521956.1_Feline_parvovirus_Egypt",
#   # "PP663061.2_Feline_parvovirus_Egypt",
#   # "PV808486.1_Feline_parvovirus_Egypt",
#   # "PP663062.2_Feline_parvovirus_Egypt",
#   # "PV808493.1_Feline_parvovirus_Egypt",
#   # "PV808496.1_Feline_parvovirus_Egypt",
#   # "PP663067.2_Feline_parvovirus_Egypt",
#   # "PV521946.1_Feline_parvovirus_Egypt",
#   "PP663047.2_Feline_parvovirus_Egypt",
#   "sample12_feline",
#   "PV521950.1_Feline_parvovirus_Egypt",
#   "PP663052.2_Feline_parvovirus_Egypt",
#   "Sample67_Feline",
#   "Sample70_Feline",
#   "Sample38_Feline",
#   "Sample_Lion",
#   "Sample69_Feline",
#   "Sample25_Feline",
#   "Sample27_Feline",
#   "Sample7_Feline",
#   "Sample44_Feline",
#   # "FPV_JO24_Sample_11",
#   "PV521947.1_Feline_parvovirus_Egypt",
#   # "PP663048.2_Feline_parvovirus_Egypt",
#   "OM638043.1_Feline_parvovirus_Egypt",
#   # "PV521964.1_Feline_parvovirus_Egypt",
#   # "PP663074.2_Feline_parvovirus_Egypt",
#   # "PV521948.1_Feline_parvovirus_Egypt",
#   # "PV808490.1_Feline_parvovirus_Egypt",
#   "PV808489.1_Feline_parvovirus_Egypt",
#   # "PV521962.1_Feline_parvovirus_Egypt",
#   # "PP663055.2_Feline_parvovirus_Egypt",
#   # "PP663056.2_Feline_parvovirus_Egypt",
#   # "PV521954.1_Feline_parvovirus_Egypt",
#   # "PV521953.1_Feline_parvovirus_Egypt",
#   # "PV808495.1_Feline_parvovirus_Egypt",
#   # "PV521960.1_Feline_parvovirus_Egypt",
#   # "PP663073.2_Feline_parvovirus_Egypt",
#   # "PP663041.2_Feline_parvovirus_Egypt",
#   # "PP663050.2_Feline_parvovirus_Egypt",
#   # "PP663070.2_Feline_parvovirus_Egypt",
#   # "PP663072.2_Feline_parvovirus_Egypt",
#   # "PV521961.1_Feline_parvovirus_Egypt",
#   "PV521949.1_Feline_parvovirus_Egypt",
#   # "PV521942.1_Feline_parvovirus_Egypt",
#   # "PP663045.2_Feline_parvovirus_Egypt",
#   # "PV521963.1_Feline_parvovirus_Egypt",
#   # "PP663071.2_Feline_parvovirus_Egypt",
#   # "PP663044.2_Feline_parvovirus_Egypt",
#   # "PV521955.1_Feline_parvovirus_Egypt",
#   "PP663059.2_Feline_parvovirus_Egypt",
#   "KX434461.1_Feline_parvovirus_IZSSI_Italy",
#   "KX900570.1_Feline_parvovirus_HH6_China",
#   "M38246.1_ref_FPV",
#   "KP081409.1_Feline_parvovirus_Iran",
#   "vaccineFeline_Zoetis",
#   "vaccineFeline_MSD",
#   "EU659112.1_ref_FPV_oldest"
# )

## FINAL_SVRNA_v3_v142

tips_to_keep <- c(
  "Sample68_Canine",
  "OM100697.1_Canine_parvovirus_2_Egypt",
  "Sample47_Canine",
  "Sample8_Canine",
  "Sample18_Canine",
  "OM100699.1_Canine_parvovirus_2_Egypt",
  # "OM100701.1_Canine_parvovirus_2_Egypt",
  # "OM100698.1_Canine_parvovirus_2_Egypt",
  "OM100696.1_Canine_parvovirus_2_Egypt",
  "MG013488.1_Canine_parvovirus_2c_China",
  "MF510157.1_Canine_parvovirus_2c_Italy",
  "Sample57_Canine",
  "MZ056886.1_Canine_parvovirus_2_Egypt",
  # "MZ056884.1_Canine_parvovirus_2_Egypt",
  "Sample77_Canine",
  "Sample76_Canine",
  "MZ056885.1_Canine_parvovirus_2_Egypt",
  # "MZ056889.1_Canine_parvovirus_2_Egypt",
  "sample51_canine",
  "Sample72_Canine",
  "OR667801.1_Canine_parvovirus_2_Iraq",
  # "OR667800.1_Canine_parvovirus_2_Iraq",
  # "OR667802.1_Canine_parvovirus_2_Iraq",
  # "OR451707.1_Canine_parvovirus_2_Iraq",
  # "OR667804.1_Canine_parvovirus_2_Iraq",
  # "OR667803.1_Canine_parvovirus_2_Iraq",
  "MW653250.1_Canine_parvovirus_2a_Iran",
  "OL330980.1_Canine_parvovirus_2_Iran",
  # "OL330978.1_Canine_parvovirus_2_Iran",
  "MW653249.1_Canine_parvovirus_2a_Iran",
  "OM100702.1_Canine_parvovirus_2_Egypt",
  # "OM100700.1_Canine_parvovirus_2_Egypt",
  "OL330979.1_Canine_parvovirus_2_Iran",
  "MZ056892.1_Canine_parvovirus_2_Egypt",
  # "MZ056883.1_Canine_parvovirus_2_Egypt",
  # "MZ056891.1_Canine_parvovirus_2_Egypt",
  "OM721655.1_Canine_parvovirus_2b_Turkey",
  "MZ056890.1_Canine_parvovirus_2_Egypt",
  # "MZ056882.1_Canine_parvovirus_2_Egypt",
  # "MZ056887.1_Canine_parvovirus_2_Egypt",
  # "MZ056888.1_Canine_parvovirus_2_Egypt",
  # "MZ056881.1_Canine_parvovirus_2_Egypt",
  "MW539053.1_Canine_parvovirus_2b_Turkey",
  # "OM721656.1_Canine_parvovirus_2b_Turkey",
  "OL330977.1_Canine_parvovirus_2_Iran",
  "KX268117.1_Canine_parvovirus_2_Turkey",
  # "KX268105.1_Canine_parvovirus_2_Turkey",
  # "KX268107.1_Canine_parvovirus_2_Turkey",
  # "KX268118.1_Canine_parvovirus_2_Turkey",
  # "KX268106.1_Canine_parvovirus_2_Turkey",
  # "KX268114.1_Canine_parvovirus_2_Turkey",
  # "KX268110.1_Canine_parvovirus_2_Turkey",
  # "KX268108.1_Canine_parvovirus_2_Turkey",
  "MW653251.1_Canine_parvovirus_2b_Iran",
  "MW653253.1_Canine_parvovirus_2c_Iran",
  # "KX268109.1_Canine_parvovirus_2_Turkey",
  "FJ222821.1_ref_CPV2C",
  "MW653256.1_Canine_parvovirus_2a_Iran",
  "KX268115.1_Canine_parvovirus_2_Turkey",
  "AY742933.1_recent_Canine_parvovirus_2a_NewZealand",
  "KX268116.1_Canine_parvovirus_2_Turkey",
  # "KX268111.1_Canine_parvovirus_2_Turkey",
  # "KX268113.1_Canine_parvovirus_2_Turkey",
  # "KX268112.1_Canine_parvovirus_2_Turkey",
  "MW653252.1_Canine_parvovirus_2b_Iran",
  "MW653248.1_Canine_parvovirus_2a_Iran",
  "M74849.1_ref_CPV2B",
  "M24003.1_ref_CPV2A",
  "M38245.1_ref_CPV2",
  "vaccineCanine_MSD",
  "KX268119.1_Canine_parvovirus_2_Turkey",
  "vaccineCanine_Zoetis",
  "Sample70_Feline",
  "Sample38_Feline",
  "Sample69_Feline",
  "Sample_Lion",
  "Sample25_Feline",
  "Sample27_Feline",
  "Sample7_Feline",
  "Sample44_Feline",
  "PV521947.1_Feline_parvovirus_Egypt",
  # "PP663048.2_Feline_parvovirus_Egypt",
  "OM638043.1_Feline_parvovirus_Egypt",
  # "PV521964.1_Feline_parvovirus_Egypt",
  # "PP663074.2_Feline_parvovirus_Egypt",
  # "PV521948.1_Feline_parvovirus_Egypt",
  # "PV808490.1_Feline_parvovirus_Egypt",
  # "PV808489.1_Feline_parvovirus_Egypt",
  # "PV521962.1_Feline_parvovirus_Egypt",
  # "PV521954.1_Feline_parvovirus_Egypt",
  # "PV521953.1_Feline_parvovirus_Egypt",
  # "PP663056.2_Feline_parvovirus_Egypt",
  # "PP663055.2_Feline_parvovirus_Egypt",
  # "PV808495.1_Feline_parvovirus_Egypt",
  # "PV521960.1_Feline_parvovirus_Egypt",
  # "PP663070.2_Feline_parvovirus_Egypt",
  # "PP663041.2_Feline_parvovirus_Egypt",
  # "PP663072.2_Feline_parvovirus_Egypt",
  # "PP663050.2_Feline_parvovirus_Egypt",
  "PP663073.2_Feline_parvovirus_Egypt",
  # "PV521961.1_Feline_parvovirus_Egypt",
  # "PV521949.1_Feline_parvovirus_Egypt",
  # "PP663045.2_Feline_parvovirus_Egypt",
  # "PV521942.1_Feline_parvovirus_Egypt",
  # "PV521963.1_Feline_parvovirus_Egypt",
  # "PP663071.2_Feline_parvovirus_Egypt",
  # "PP663044.2_Feline_parvovirus_Egypt",
  # "PV521955.1_Feline_parvovirus_Egypt",
  # "PP663059.2_Feline_parvovirus_Egypt",
  # "PV808491.1_Feline_parvovirus_Egypt",
  # "PV521951.1_Feline_parvovirus_Egypt",
  # "PV808488.1_Feline_parvovirus_Egypt",
  # "PV521944.1_Feline_parvovirus_Egypt",
  # "PV521943.1_Feline_parvovirus_Egypt",
  # "PP663042.2_Feline_parvovirus_Egypt",
  # "PV521965.1_Feline_parvovirus_Egypt",
  # "PV521945.1_Feline_parvovirus_Egypt",
  # "PV521959.1_Feline_parvovirus_Egypt",
  # "PV521958.1_Feline_parvovirus_Egypt",
  # "PV521957.1_Feline_parvovirus_Egypt",
  # "PV808492.1_Feline_parvovirus_Egypt",
  # "PV521952.1_Feline_parvovirus_Egypt",
  # "PV521956.1_Feline_parvovirus_Egypt",
  # "PP663061.2_Feline_parvovirus_Egypt",
  # "PV808496.1_Feline_parvovirus_Egypt",
  # "PP663062.2_Feline_parvovirus_Egypt",
  # "PV808493.1_Feline_parvovirus_Egypt",
  # "PV808486.1_Feline_parvovirus_Egypt",
  # "PP663067.2_Feline_parvovirus_Egypt",
  # "PV521946.1_Feline_parvovirus_Egypt",
  "PP663047.2_Feline_parvovirus_Egypt",
  "sample12_feline",
  "PV521950.1_Feline_parvovirus_Egypt",
  # "PP663052.2_Feline_parvovirus_Egypt",
  "Sample67_Feline",
  "KX434461.1_Feline_parvovirus_IZSSI_Italy",
  "KX900570.1_Feline_parvovirus_HH6_China",
  "M38246.1_ref_FPV",
  "KP081409.1_Feline_parvovirus_Iran",
  "vaccineFeline_Zoetis",
  "vaccineFeline_MSD",
  "EU659112.1_ref_FPV_oldest"
  )

length(tips_to_keep) ## 65

tips_to_drop <- setdiff(dna_ML_Muscle_trimmed_tree$tip.label, tips_to_keep)
dna_ML_Muscle_trimmed_tree_filtered <- drop.tip(dna_ML_Muscle_trimmed_tree, tips_to_drop)
dna_ML_Muscle_trimmed_tree_filtered_midpoint <- midpoint(dna_ML_Muscle_trimmed_tree_filtered)
dna_ML_Muscle_trimmed_tree_filtered_rooted_oldestFPV <- root(dna_ML_Muscle_trimmed_tree_filtered, outgroup = "EU659112.1_ref_FPV_oldest", resolve.root = TRUE)

length(dna_ML_Muscle_trimmed_tree$tip.label)
length(dna_ML_Muscle_trimmed_tree_filtered_midpoint$tip.label)
length(dna_ML_Muscle_trimmed_tree_filtered_rooted_oldestFPV$tip.label)

ggtreeplot_dna_ML_Muscle_trimmed_tree_filtered <- ggtree(dna_ML_Muscle_trimmed_tree_filtered, layout = "rectangular", size = 0.3) %<+% tip_df_dna +
  geom_tiplab(aes(color = group), size = 2, offset = 0.002, align = TRUE, linesize = 0.1) +
  geom_tippoint(aes(color = group), size = 1.2) +
  geom_text2(aes(label = label, subset = !isTip), size = 1.8, hjust = 1.2, vjust = -0.4) +
  scale_color_manual(values = group_colors) +
  theme_tree2() +
  labs(title = paste0("Parvovirus DNA - ML", deparse(substitute(dna_ML_Muscle_trimmed_tree_filtered))), color = "Group") +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 8)) +
  coord_cartesian(clip = "off") +                  
  hexpand(0.5, direction = 1)

ggtreeplot_dna_ML_Muscle_trimmed_tree_filtered_midpoint <- ggtree(dna_ML_Muscle_trimmed_tree_filtered_midpoint, layout = "rectangular", size = 0.3) %<+% tip_df_dna +
  geom_tiplab(aes(color = group), size = 2, offset = 0.002, align = TRUE, linesize = 0.1) +
  geom_tippoint(aes(color = group), size = 1.2) +
  geom_text2(aes(label = label, subset = !isTip), size = 1.8, hjust = 1.2, vjust = -0.4) +
  scale_color_manual(values = group_colors) +
  theme_tree2() +
  labs(title = paste0("Parvovirus DNA - ML", deparse(substitute(dna_ML_Muscle_trimmed_tree_filtered_midpoint))), color = "Group") +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 8)) +
  coord_cartesian(clip = "off") +                  
  hexpand(0.5, direction = 1)

ggtreeplot_dna_ML_Muscle_trimmed_tree_filtered_rooted_oldestFPV <- ggtree(dna_ML_Muscle_trimmed_tree_filtered_rooted_oldestFPV, layout = "rectangular", size = 0.3) %<+% tip_df_dna +
  geom_tiplab(aes(color = group), size = 2, offset = 0.002, align = TRUE, linesize = 0.1) +
  geom_tippoint(aes(color = group), size = 1.2) +
  geom_text2(aes(label = label, subset = !isTip), size = 1.8, hjust = 1.2, vjust = -0.4) +
  scale_color_manual(values = group_colors) +
  theme_tree2() +
  labs(
    # title = paste0("Parvovirus DNA - ML", deparse(substitute(dna_ML_Muscle_trimmed_tree_filtered_rooted_oldestFPV))),
    color = "Group") +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 8)) +
  coord_cartesian(clip = "off") +                  
  hexpand(0.5, direction = 1)

## NO Bootstrap values 

ggtree(dna_ML_Muscle_trimmed_tree_filtered_rooted_oldestFPV, layout = "rectangular", size = 0.3) %<+% tip_df_dna + 
  geom_tiplab(aes(color = group), size = 2, offset = 0.002, align = TRUE, linesize = 0.1) + 
  geom_tippoint(aes(color = group), size = 1.2) + 
  geom_point2(aes(subset = !isTip & !is.na(as.numeric(label)), 
                  size = cut(as.numeric(label), 
                             breaks = c(-Inf, 50, 85, Inf), 
                             labels = c("<50", "50-84", ">85"))), 
              shape = 16, colour = "gray35") + 
  scale_size_manual(values = c("<50" = 1, 
                               "50-84" = 1.8, 
                               ">85" = 2.8 
  ), name = "Bootstrap", drop = FALSE, na.translate = FALSE) + 
  scale_color_manual(values = group_colors) + 
  theme_tree2() + 
  labs(title = paste0("Parvovirus DNA - ML", deparse(substitute(dna_ML_Muscle_trimmed_tree_filtered_rooted_oldestFPV))), color = "Group") + 
  theme(legend.position = "bottom", 
        plot.title = element_text(size = 14, face = "bold"), 
        legend.text = element_text(size = 8), 
        legend.title = element_text(size = 9)) + 
  coord_cartesian(clip = "off") + 
  hexpand(0.5, direction = 1)

## Protein

protein_ML_Muscle_tree <- read.tree(file.path(sys_dir,"Parvo_sanger_full_seq/iqtree_output/paper/Maxll_Protein_aln_Muscle_FINAL_SVRNA_v3_v142.contree"))

tip_df_protein <- data.frame(
  label = protein_ML_Muscle_tree$tip.label,
  group = ifelse(grepl("sample", protein_ML_Muscle_tree$tip.label, ignore.case = TRUE) & grepl("feline", protein_ML_Muscle_tree$tip.label, ignore.case = TRUE), "Jo_Feline",
                 ifelse(grepl("sample", protein_ML_Muscle_tree$tip.label, ignore.case = TRUE) & grepl("canine", protein_ML_Muscle_tree$tip.label, ignore.case = TRUE), "Jo_Canine",
                        ifelse(grepl("sample", protein_ML_Muscle_tree$tip.label, ignore.case = TRUE) & grepl("lion",   protein_ML_Muscle_tree$tip.label, ignore.case = TRUE), "Jo_Lion",
                               ifelse(grepl("NC_001510|NC_001718|NC_029797|NC_038544", protein_ML_Muscle_tree$tip.label), "Outgroup",
                                      ifelse(grepl("_ref_", protein_ML_Muscle_tree$tip.label), "Reference",
                                             ifelse(grepl("vaccine", protein_ML_Muscle_tree$tip.label, ignore.case = TRUE), "Vaccine",
                                                    ifelse(grepl("feline", protein_ML_Muscle_tree$tip.label, ignore.case = TRUE), "Regional_Feline",
                                                           ifelse(grepl("canine", protein_ML_Muscle_tree$tip.label, ignore.case = TRUE), "Regional_Canine",
                                                                  "Other"))))))))
)

## rooting v2 

protein_ML_Muscle_tree_midpoint <- midpoint(protein_ML_Muscle_tree)
protein_ML_Muscle_tree_rooted_oldestFPV <- root(protein_ML_Muscle_tree, outgroup = "EU659112.1_ref_FPV_oldest", resolve.root = TRUE)

parsimony(protein_ML_Muscle_tree, protein_phy_Muscle)
parsimony(protein_ML_Muscle_tree_midpoint,  protein_phy_Muscle)
parsimony(protein_ML_Muscle_tree_rooted_oldestFPV, protein_phy_Muscle)

parsimony(protein_NJ_Muscle_midpoint, protein_phy_Muscle)
parsimony(protein_NJ_Muscle_rooted_oldestFPV, protein_phy_Muscle)

ggtreeplot_protein_ML_Muscle_tree <- ggtree(protein_ML_Muscle_tree, layout = "rectangular", size = 0.3) %<+% tip_df_protein +
  geom_tiplab(aes(color = group), size = 2, offset = 0.002, align = TRUE, linesize = 0.1) +
  geom_tippoint(aes(color = group), size = 1.2) +
  geom_text2(aes(label = label, subset = !isTip), size = 1.8, hjust = 1.2, vjust = -0.4) +
  scale_color_manual(values = group_colors) +
  theme_tree2() +
  labs(title = paste0("Parvovirus Protein - ML", deparse(substitute(protein_ML_Muscle_tree))), color = "Group") +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 8)) +
  coord_cartesian(clip = "off") +
  hexpand(0.5, direction = 1)

ggtreeplot_protein_ML_Muscle_tree_midpoint <- ggtree(protein_ML_Muscle_tree_midpoint, layout = "rectangular", size = 0.3) %<+% tip_df_protein +
  geom_tiplab(aes(color = group), size = 2, offset = 0.002, align = TRUE, linesize = 0.1) +
  geom_tippoint(aes(color = group), size = 1.2) +
  geom_text2(aes(label = label, subset = !isTip), size = 1.8, hjust = 1.2, vjust = -0.4) +
  scale_color_manual(values = group_colors) +
  theme_tree2() +
  labs(title = paste0("Parvovirus Protein - ML", deparse(substitute(protein_ML_Muscle_tree_midpoint))), color = "Group") +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 8)) +
  coord_cartesian(clip = "off") +
  hexpand(0.5, direction = 1)

ggtreeplot_protein_ML_Muscle_tree_rooted_oldestFPV <- ggtree(protein_ML_Muscle_tree_rooted_oldestFPV, layout = "rectangular", size = 0.3) %<+% tip_df_protein +
  geom_tiplab(aes(color = group), size = 2, offset = 0.002, align = TRUE, linesize = 0.1) +
  geom_tippoint(aes(color = group), size = 1.2) +
  geom_text2(aes(label = label, subset = !isTip), size = 1.8, hjust = 1.2, vjust = -0.4) +
  scale_color_manual(values = group_colors) +
  theme_tree2() +
  labs(
    # title = paste0("Parvovirus Protein - ML", deparse(substitute(protein_ML_Muscle_tree_rooted_oldestFPV))),
    color = "Group") +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 8)) +
  coord_cartesian(clip = "off") +
  hexpand(0.5, direction = 1)

## more filtering

protein_ML_Muscle_tree_filtered <- drop.tip(protein_ML_Muscle_tree, tips_to_drop)
protein_ML_Muscle_tree_filtered_midpoint <- midpoint(protein_ML_Muscle_tree_filtered)
protein_ML_Muscle_tree_filtered_rooted_oldestFPV <- root(protein_ML_Muscle_tree_filtered, outgroup = "EU659112.1_ref_FPV_oldest", resolve.root = TRUE)
protein_ML_Muscle_tree_filtered_rooted_FPV <- root(protein_ML_Muscle_tree_filtered, outgroup = "M38246.1_ref_FPV", resolve.root = TRUE)

length(protein_ML_Muscle_tree$tip.label)
length(protein_ML_Muscle_tree_filtered_midpoint$tip.label)
length(protein_ML_Muscle_tree_filtered_rooted_oldestFPV$tip.label)

ggtreeplot_protein_ML_Muscle_tree_filtered <- ggtree(protein_ML_Muscle_tree_filtered, layout = "rectangular", size = 0.3) %<+% tip_df_protein +
  geom_tiplab(aes(color = group), size = 2, offset = 0.002, align = TRUE, linesize = 0.1) +
  geom_tippoint(aes(color = group), size = 1.2) +
  geom_text2(aes(label = label, subset = !isTip), size = 1.8, hjust = 1.2, vjust = -0.4) +
  scale_color_manual(values = group_colors) +
  theme_tree2() +
  labs(title = paste0("Parvovirus Protein - ML", deparse(substitute(protein_ML_Muscle_tree_filtered))), color = "Group") +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 8)) +
  coord_cartesian(clip = "off") +                  
  hexpand(0.5, direction = 1)

ggtreeplot_protein_ML_Muscle_tree_filtered_midpoint <- ggtree(protein_ML_Muscle_tree_filtered_midpoint, layout = "rectangular", size = 0.3) %<+% tip_df_protein +
  geom_tiplab(aes(color = group), size = 2, offset = 0.002, align = TRUE, linesize = 0.1) +
  geom_tippoint(aes(color = group), size = 1.2) +
  geom_text2(aes(label = label, subset = !isTip), size = 1.8, hjust = 1.2, vjust = -0.4) +
  scale_color_manual(values = group_colors) +
  theme_tree2() +
  labs(title = paste0("Parvovirus Protein - ML", deparse(substitute(protein_ML_Muscle_tree_filtered_midpoint))), color = "Group") +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 8)) +
  coord_cartesian(clip = "off") +                  
  hexpand(0.5, direction = 1)

ggtreeplot_protein_ML_Muscle_tree_filtered_rooted_oldestFPV <- ggtree(protein_ML_Muscle_tree_filtered_rooted_oldestFPV, layout = "rectangular", size = 0.3) %<+% tip_df_protein +
  geom_tiplab(aes(color = group), size = 2, offset = 0.002, align = TRUE, linesize = 0.1) +
  geom_tippoint(aes(color = group), size = 1.2) +
  geom_text2(aes(label = label, subset = !isTip), size = 1.8, hjust = 1.2, vjust = -0.4) +
  scale_color_manual(values = group_colors) +
  theme_tree2() +
  labs(
    # title = paste0("Parvovirus Protein - ML", deparse(substitute(protein_ML_Muscle_tree_filtered_rooted_oldestFPV))),
    color = "Group") +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 8)) +
  coord_cartesian(clip = "off") +                  
  hexpand(0.5, direction = 1)

## filter by bootstrap value

bootstrap_cutoff <- 50

node_support <- as.numeric(dna_ML_Muscle_trimmed_tree$node.label)
Ntip <- length(dna_ML_Muscle_trimmed_tree$tip.label)
low_support_nodes <- which(node_support < bootstrap_cutoff) + Ntip

tips_to_drop <- unique(unlist(lapply(low_support_nodes, function(n) {
  extract.clade(dna_ML_Muscle_trimmed_tree, n)$tip.label
})))

dna_ML_Muscle_trimmed_tree_filtered <- drop.tip(dna_ML_Muscle_trimmed_tree, tips_to_drop)

tip_df_dna_filtered <- data.frame(
  label = dna_ML_Muscle_trimmed_tree_filtered$tip.label,
  group = ifelse(grepl("sample", dna_ML_Muscle_trimmed_tree_filtered$tip.label, ignore.case = TRUE) & grepl("feline", dna_ML_Muscle_trimmed_tree_filtered$tip.label, ignore.case = TRUE), "Jo_Feline",
                 ifelse(grepl("sample", dna_ML_Muscle_trimmed_tree_filtered$tip.label, ignore.case = TRUE) & grepl("canine", dna_ML_Muscle_trimmed_tree_filtered$tip.label, ignore.case = TRUE), "Jo_Canine",
                        ifelse(grepl("sample", dna_ML_Muscle_trimmed_tree_filtered$tip.label, ignore.case = TRUE) & grepl("lion",   dna_ML_Muscle_trimmed_tree_filtered$tip.label, ignore.case = TRUE), "Jo_Lion",
                               ifelse(grepl("NC_001510|NC_001718|NC_029797|NC_038544", dna_ML_Muscle_trimmed_tree_filtered$tip.label), "Outgroup",
                                      ifelse(grepl("_ref_", dna_ML_Muscle_trimmed_tree_filtered$tip.label), "Reference",
                                             ifelse(grepl("vaccine", dna_ML_Muscle_trimmed_tree_filtered$tip.label, ignore.case = TRUE), "Vaccine",
                                                    ifelse(grepl("feline", dna_ML_Muscle_trimmed_tree_filtered$tip.label, ignore.case = TRUE), "Regional_Feline",
                                                           ifelse(grepl("canine", dna_ML_Muscle_trimmed_tree_filtered$tip.label, ignore.case = TRUE), "Regional_Canine",
                                                                  "Other"))))))))
)

ggtree(dna_ML_Muscle_trimmed_tree_filtered, layout = "rectangular", size = 0.3) %<+% tip_df_dna_filtered +
  geom_tiplab(aes(color = group), size = 2, offset = 0.002, align = TRUE, linesize = 0.1) +
  geom_tippoint(aes(color = group), size = 1.2) +
  geom_text2(aes(label = label, subset = !isTip), size = 1.8, hjust = 1.2, vjust = -0.4) +
  scale_color_manual(values = group_colors) +
  theme_tree2() +
  labs(title = paste0("Parvovirus DNA - ML filtered (BS >= ", bootstrap_cutoff, ")"), color = "Group") +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 8)) +
  coord_cartesian(clip = "off") +
  hexpand(0.5, direction = 1)

ggtreeplot_dna_ML_Muscle_trimmed_tree_rooted_oldestFPV

dna_ML_Muscle_trimmed_tree_rooted_oldestFPV

bootstrap_cutoff <- 50

## get bootstrap support values for internal nodes (node.label order matches internal node order)
node_support_dna <- (as.numeric(dna_ML_Muscle_trimmed_tree_filtered_rooted_oldestFPV$node.label))
node_support_dna[is.na(node_support_dna)] <- 0

## internal nodes are numbered after tips: Ntip+1 ... Ntip+Nnode
Ntip_dna <- Ntip(dna_ML_Muscle_trimmed_tree_filtered_rooted_oldestFPV)
node_ids_dna <- seq_along(node_support_dna) + Ntip_dna

## which internal nodes fall below cutoff
low_support_node_ids_dna <- node_ids_dna[node_support_dna < bootstrap_cutoff]
low_support_node_ids_dna

## find edges leading into those nodes (edge[,2] = child node of each edge)
edges_to_zero_dna <- dna_ML_Muscle_trimmed_tree_filtered_rooted_oldestFPV$edge[,2] %in% low_support_node_ids_dna
sum(edges_to_zero_dna)

## set those edge lengths to zero
dna_ML_Muscle_trimmed_tree_rooted_oldestFPV_collapsed <- dna_ML_Muscle_trimmed_tree_filtered_rooted_oldestFPV
dna_ML_Muscle_trimmed_tree_rooted_oldestFPV_collapsed$edge.length[edges_to_zero_dna] <- 0

## collapse zero-length branches into polytomies
dna_ML_Muscle_trimmed_tree_rooted_oldestFPV_collapsed <- di2multi(dna_ML_Muscle_trimmed_tree_rooted_oldestFPV_collapsed, tol = 1e-8)

## check result
Ntip(dna_ML_Muscle_trimmed_tree_filtered_rooted_oldestFPV)
Ntip(dna_ML_Muscle_trimmed_tree_rooted_oldestFPV_collapsed)
Nnode(dna_ML_Muscle_trimmed_tree_filtered_rooted_oldestFPV)
Nnode(dna_ML_Muscle_trimmed_tree_rooted_oldestFPV_collapsed)

## plot
ggtree(dna_ML_Muscle_trimmed_tree_rooted_oldestFPV_collapsed, layout = "rectangular", size = 0.3) %<+% tip_df_dna +
  geom_tiplab(aes(color = group), size = 2, offset = 0.002, align = TRUE, linesize = 0.1) +
  geom_tippoint(aes(color = group), size = 1.2) +
  scale_color_manual(values = group_colors) +
  theme_tree2() +
  labs(title = paste0("Parvovirus DNA - ML collapsed (BS < ", bootstrap_cutoff, ")"), color = "Group") +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 8)) +
  coord_cartesian(clip = "off") +
  hexpand(0.5, direction = 1)

bootstrap_cutoff <- 30

node_support_dna <- as.numeric(dna_ML_Muscle_trimmed_tree_rooted_oldestFPV$node.label)
Ntip_dna <- length(dna_ML_Muscle_trimmed_tree_rooted_oldestFPV$tip.label)
Nnode_dna <- Nnode(dna_ML_Muscle_trimmed_tree_rooted_oldestFPV)

low_support_nodes_dna <- which(node_support_dna < bootstrap_cutoff) + Ntip_dna
low_support_nodes_dna

tips_to_drop_dna <- unique(unlist(lapply(low_support_nodes_dna, function(n) {
  extract.clade(dna_ML_Muscle_trimmed_tree_rooted_oldestFPV, n)$tip.label
})))

length(tips_to_drop_dna)
tips_to_drop_dna

dna_ML_Muscle_trimmed_tree_rooted_oldestFPV_filtered <- drop.tip(dna_ML_Muscle_trimmed_tree_rooted_oldestFPV, tips_to_drop_dna)

length(dna_ML_Muscle_trimmed_tree_rooted_oldestFPV$tip.label)
length(dna_ML_Muscle_trimmed_tree_rooted_oldestFPV_filtered$tip.label)

ggtree(dna_ML_Muscle_trimmed_tree_rooted_oldestFPV_filtered, layout = "rectangular", size = 0.3) %<+% tip_df_dna +
  geom_tiplab(aes(color = group), size = 2, offset = 0.002, align = TRUE, linesize = 0.1) +
  geom_tippoint(aes(color = group), size = 1.2) +
  scale_color_manual(values = group_colors) +
  theme_tree2() +
  labs(title = paste0("Parvovirus DNA - ML collapsed (BS < ", bootstrap_cutoff, ")"), color = "Group") +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 8)) +
  coord_cartesian(clip = "off") +
  hexpand(0.5, direction = 1)

## TEST - replacing bootstrap values with SH_aLRT and/or site coverage

support_values <- do.call(rbind, strsplit(dna_treefile$node.label, "/"))
SH_aLRT <- as.numeric(support_values[,1])
UFBoot <- as.numeric(support_values[,2])

dna_treefile$node.label <- paste0(round(SH_aLRT, 1), "/", UFBoot)

ggtree(dna_treefile, layout = "rectangular", size = 0.3) %<+% tip_df_dna +
  geom_tiplab(aes(color = group), size = 2, offset = 0.002, align = TRUE, linesize = 0.1) +
  geom_tippoint(aes(color = group), size = 1.2) +
  geom_text2(aes(label = label, subset = !isTip), size = 1.8, hjust = 1.2, vjust = -0.4) +
  scale_color_manual(values = group_colors) +
  theme_tree2() +
  labs(title = paste0("Parvovirus DNA - ML", deparse(substitute(dna_treefile))),
       color = "Group") +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 8)) +
  coord_cartesian(clip = "off") +
  hexpand(0.5, direction = 1)

## calculate per-site coverage across alignment
## assign site coverage to internal nodes

tree <- dna_ML_Muscle_trimmed_tree_filtered_rooted_oldestFPV
node_coverage <- sapply((Ntip(tree) + 1):(Ntip(tree) + Nnode(tree)), function(node) {
  tips <- extract.clade(tree, node)$tip.label
  cols <- which(colSums(aln_mat[tips, , drop = FALSE] != "-") > 0)
  mean(site_coverage[cols], na.rm = TRUE)
})
tree$node.label <- round(node_coverage, 1)

ggtree(tree, layout = "rectangular", size = 0.3) %<+% tip_df_dna +
  geom_tiplab(aes(color = group), 
              size = 2, offset = 0.002, align = TRUE, linesize = 0.1) +
  geom_tippoint(aes(color = group), size = 1.2) +
  geom_text2(aes(label = label, subset = !isTip),
             size = 1.8, hjust = 1.2, vjust = -0.4) +
  scale_color_manual(values = group_colors) +
  theme_tree2() +
  labs(title = "Parvovirus DNA - ML (site coverage %)",
       color = "Group") +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 8)) +
  coord_cartesian(clip = "off") +
  hexpand(0.5, direction = 1)

# Bayesian Analysis of Phylogeny trees - DNA & Protein analysis -----------

#### DNA

dna_list_for_nexus <- strsplit(as.character(dna_alignment_Muscle_trimmed), "")
names(dna_list_for_nexus) <- names(dna_alignment_Muscle_trimmed)

## sanity check before writing
length(dna_list_for_nexus)              ## should be 143
length(dna_list_for_nexus[[1]])         ## should be 1932, not 1

# write.nexus.data(dna_list_for_nexus,
#                  file = file.path(sys_dir, "Parvo_sanger_full_seq/test_dna_msa_trimmed_mrbayes.nex"),
#                  format = "dna",
#                  interleaved = FALSE)
# mrbayes_block <- "
# begin mrbayes;
#   lset nst=6 rates=gamma;
#   mcmc ngen=2000000 samplefreq=1000 printfreq=1000 diagnfreq=10000 nchains=4;
#   sumt burnin=25;
#   sump burnin=25;
# end;
# "
# cat(mrbayes_block,
#     file = "/home/shadi/Desktop/ParvoProject/Parvo_sanger_full_seq/test_dna_msa_trimmed_mrbayes.nex",
#     append = TRUE)

write.nexus.data(dna_list_for_nexus,
                 file = file.path(sys_dir, "Parvo_sanger_full_seq/mrbayes_output/mrbayes_dna_MSA_trimmed_FINAL_SVRNA_v3.nex"),
                 format = "dna",
                 interleaved = FALSE)
mrbayes_block <- "
begin mrbayes;
  set autoclose=yes nowarn=yes;
  lset nst=2 rates=gamma;
  prset statefreqpr=fixed(empirical);
  mcmc ngen=2000000 samplefreq=1000 printfreq=1000 diagnfreq=10000 nchains=4;
  sump burnin=500;
  sumt burnin=500;
end;
"
cat(mrbayes_block,
    file = "/home/shadi/Desktop/ParvoProject/Parvo_sanger_full_seq/mrbayes_output/mrbayes_dna_MSA_trimmed_FINAL_SVRNA_v3.nex",
    append = TRUE)

## read the tree

mrbayes_tree_dna <- read.nexus("/home/shadi/Desktop/ParvoProject/Parvo_sanger_full_seq/mrbayes_output/mrbayes_dna_MSA_trimmed_FINAL_SVRNA_v3.nex.con.tre")
mrbayes_tree_dna_rooted <- root(mrbayes_tree_dna, outgroup = "EU659112.1_ref_FPV_oldest", resolve.root = TRUE)

tip_df_mrbayes <- data.frame(
  label = mrbayes_tree_dna_rooted$tip.label,
  group = ifelse(grepl("sample", mrbayes_tree_dna_rooted$tip.label, ignore.case = TRUE) & grepl("feline", mrbayes_tree_dna_rooted$tip.label, ignore.case = TRUE), "Jo_Feline",
                 ifelse(grepl("sample", mrbayes_tree_dna_rooted$tip.label, ignore.case = TRUE) & grepl("canine", mrbayes_tree_dna_rooted$tip.label, ignore.case = TRUE), "Jo_Canine",
                        ifelse(grepl("sample", mrbayes_tree_dna_rooted$tip.label, ignore.case = TRUE) & grepl("lion", mrbayes_tree_dna_rooted$tip.label, ignore.case = TRUE), "Jo_Lion",
                               ifelse(grepl("NC_001510|NC_001718|NC_029797|NC_038544", mrbayes_tree_dna_rooted$tip.label), "Outgroup",
                                      ifelse(grepl("_ref_", mrbayes_tree_dna_rooted$tip.label), "Reference",
                                             ifelse(grepl("vaccine", mrbayes_tree_dna_rooted$tip.label, ignore.case = TRUE), "Vaccine",
                                                    ifelse(grepl("feline", mrbayes_tree_dna_rooted$tip.label, ignore.case = TRUE), "Regional_Feline",
                                                           ifelse(grepl("canine", mrbayes_tree_dna_rooted$tip.label, ignore.case = TRUE), "Regional_Canine",
                                                                  "Other"))))))))
)

ggtree(mrbayes_tree_dna_rooted, layout = "rectangular", size = 0.3) %<+% tip_df_mrbayes +
  geom_tiplab(aes(color = group), size = 2, offset = 0.002, align = TRUE, linesize = 0.1) +
  geom_tippoint(aes(color = group), size = 1.2) +
  geom_text2(aes(label = round(as.numeric(label), 2), subset = !isTip), size = 1.8, hjust = 1.2, vjust = -0.4) +
  scale_color_manual(values = group_colors) +
  theme_tree2() +
  labs(title = "Parvovirus DNA - MrBayes consensus DNA", color = "Group") +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 8)) +
  coord_cartesian(clip = "off") +
  hexpand(0.5, direction = 1)

## with probability

mrbayes_tree_dna <- read.beast("/home/shadi/Desktop/ParvoProject/Parvo_sanger_full_seq/mrbayes_output/mrbayes_dna_MSA_trimmed_FINAL_SVRNA_v3.nex.con.tre")
mrbayes_data <- as_tibble(mrbayes_tree_dna)[, c("node", "prob")]
class(mrbayes_data) <- c("tbl_df", "tbl", "data.frame")
n_tips <- length(mrbayes_tree_dna@phylo$tip.label)
print(mrbayes_data[mrbayes_data$prob > 0.95, ], n = 30)
nrow(mrbayes_data[mrbayes_data$prob > 0.95, ])
nrow(mrbayes_data[, ])

## MrBayes trees is now built on FINAL_SVRNA_v3 with 143 seq, so lets remove FPV_JO24_Sample_11 manually 

seq_to_remove <- "FPV_JO24_Sample_11"
mrbayes_tree_dna_noJO24 <- drop.tip(mrbayes_tree_dna, seq_to_remove)

length(mrbayes_tree_dna_noJO24@phylo$tip.label)
length(mrbayes_tree_dna_noJO24@phylo$tip.label)
seq_to_remove %in% mrbayes_tree_dna_noJO24@phylo$tip.label

mrbayes_tree_dna_noJO24_rooted_oldestFPV <- root(as.phylo(mrbayes_tree_dna_noJO24), outgroup = "EU659112.1_ref_FPV_oldest", resolve.root = TRUE)
mrbayes_tree_dna_noJO24_rooted_oldestFPV <- treedata(phylo = mrbayes_tree_dna_noJO24_rooted_oldestFPV, data = mrbayes_data)

ggtreeplot_mrbayes_tree_dna_noJO24_rooted_oldestFPV <- ggtree(mrbayes_tree_dna_noJO24_rooted_oldestFPV, layout = "rectangular", size = 0.3) %<+% tip_df_mrbayes +
  geom_tiplab(aes(color = group), size = 2, offset = 0.002, align = TRUE, linesize = 0.1) +
  geom_tippoint(aes(color = group), size = 1.2) +
  geom_text2(aes(label = ifelse(!isTip & as.numeric(prob) >= 0, round(as.numeric(prob), 2), NA)), size = 1.8, hjust = 1.2, vjust = -0.4) +
  scale_color_manual(values = group_colors) +
  theme_tree2() +
  labs(
    # title = "Parvovirus DNA - MrBayes consensus DNA",
    color = "Group") +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 8)) +
  coord_cartesian(clip = "off") +
  hexpand(0.5, direction = 1)

parsimony(dna_NJ_Muscle_trimmed_rooted_oldestFPV, dna_phy_Muscle_trimmed)
parsimony(dna_ML_Muscle_trimmed_tree_rooted_oldestFPV, dna_phy_Muscle_trimmed)
# parsimony(mrbayes_tree_dna_noJO24_rooted_oldestFPV, dna_phy_Muscle_trimmed) 
parsimony(as.phylo(mrbayes_tree_dna_noJO24_rooted_oldestFPV), dna_phy_Muscle_trimmed)
# parsimony(mrbayes_tree_dna_noJO24_rooted_oldestFPV@phylo, dna_phy_Muscle_trimmed)

## filtration

tips_to_keep

# tips_to_drop <- setdiff(mrbayes_tree_dna_noJO24$tip.label, tips_to_keep)
tips_to_drop <- setdiff(as.phylo(mrbayes_tree_dna_noJO24)$tip.label, tips_to_keep)
mrbayes_tree_dna_noJO24_filtered <- drop.tip(mrbayes_tree_dna_noJO24, tips_to_drop)
mrbayes_tree_dna_noJO24_filtered_rooted_oldestFPV <- root(mrbayes_tree_dna_noJO24_filtered, outgroup = "EU659112.1_ref_FPV_oldest", resolve.root = TRUE)

# length(mrbayes_tree_dna_noJO24$tip.label)
# length(mrbayes_tree_dna_noJO24_filtered_rooted_oldestFPV$tip.label)
length(as.phylo(mrbayes_tree_dna_noJO24)$tip.label)
length(as.phylo(mrbayes_tree_dna_noJO24_filtered)$tip.label)

ggtree(mrbayes_tree_dna_noJO24_filtered_rooted_oldestFPV, layout = "rectangular", size = 0.3) %<+% tip_df_mrbayes +
  geom_tiplab(aes(color = group), size = 2, offset = 0.002, align = TRUE, linesize = 0.1) +
  geom_tippoint(aes(color = group), size = 1.2) +
  geom_text2(aes(label = label, subset = !isTip), size = 1.8, hjust = 1.2, vjust = -0.4) +
  scale_color_manual(values = group_colors) +
  theme_tree2() +
  labs(title = "Parvovirus DNA - MrBayes consensus", color = "Group") +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 8)) +
  coord_cartesian(clip = "off") +                  
  hexpand(0.5, direction = 1)

#### protein

protein_list_for_nexus <- strsplit(as.character(protein_alignment_Muscle), "")
names(protein_list_for_nexus) <- names(protein_alignment_Muscle)

## sanity check before writing
length(protein_list_for_nexus)              ## should be 143
length(protein_list_for_nexus[[1]])         ## should not be 1

# write.nexus.data(protein_list_for_nexus,
#                  file = file.path(sys_dir, "Parvo_sanger_full_seq/test_protein_msa_mrbayes.nex"),
#                  format = "protein",
#                  interleaved = FALSE)
# 
# mrbayes_block <- "
# begin mrbayes;
#   lset nst=6 rates=gamma;
#   mcmc ngen=2000000 samplefreq=1000 printfreq=1000 diagnfreq=10000 nchains=4;
#   sumt burnin=25;
#   sump burnin=25;
# end;
# "
# cat(mrbayes_block,
#     file = "/home/shadi/Desktop/ParvoProject/Parvo_sanger_full_seq/test_protein_msa_mrbayes.nex",
#     append = TRUE)

write.nexus.data(protein_list_for_nexus,
                 file = file.path(sys_dir, "Parvo_sanger_full_seq/mrbayes_output/mrbayes_protein_MSA_FINAL_SVRNA_v3.nex"),
                 format = "protein",
                 interleaved = FALSE)
mrbayes_block <- "
begin mrbayes;
  set autoclose=yes nowarn=yes;
  prset aamodelpr=fixed(jones);
  lset rates=kmixture nmixtcat=2;
  mcmc ngen=2000000 samplefreq=1000 printfreq=1000 diagnfreq=10000 nchains=4;
  sump burnin=500;
  sumt burnin=500;
end;
"
cat(mrbayes_block,
    file = "/home/shadi/Desktop/ParvoProject/Parvo_sanger_full_seq/mrbayes_output/mrbayes_protein_MSA_FINAL_SVRNA_v3.nex",
    append = TRUE)

## read the tree

mrbayes_tree_protein <- read.nexus("/home/shadi/Desktop/ParvoProject/Parvo_sanger_full_seq/mrbayes_output/mrbayes_protein_MSA_FINAL_SVRNA_v3.nex.con.tre")
mrbayes_tree_protein_rooted <- root(mrbayes_tree_protein, outgroup = "EU659112.1_ref_FPV_oldest", resolve.root = TRUE)

tip_df_mrbayes <- data.frame(
  label = mrbayes_tree_protein_rooted$tip.label,
  group = ifelse(grepl("sample", mrbayes_tree_protein_rooted$tip.label, ignore.case = TRUE) & grepl("feline", mrbayes_tree_protein_rooted$tip.label, ignore.case = TRUE), "Jo_Feline",
                 ifelse(grepl("sample", mrbayes_tree_protein_rooted$tip.label, ignore.case = TRUE) & grepl("canine", mrbayes_tree_protein_rooted$tip.label, ignore.case = TRUE), "Jo_Canine",
                        ifelse(grepl("sample", mrbayes_tree_protein_rooted$tip.label, ignore.case = TRUE) & grepl("lion", mrbayes_tree_protein_rooted$tip.label, ignore.case = TRUE), "Jo_Lion",
                               ifelse(grepl("NC_001510|NC_001718|NC_029797|NC_038544", mrbayes_tree_protein_rooted$tip.label), "Outgroup",
                                      ifelse(grepl("_ref_", mrbayes_tree_protein_rooted$tip.label), "Reference",
                                             ifelse(grepl("vaccine", mrbayes_tree_protein_rooted$tip.label, ignore.case = TRUE), "Vaccine",
                                                    ifelse(grepl("feline", mrbayes_tree_protein_rooted$tip.label, ignore.case = TRUE), "Regional_Feline",
                                                           ifelse(grepl("canine", mrbayes_tree_protein_rooted$tip.label, ignore.case = TRUE), "Regional_Canine",
                                                                  "Other"))))))))
)

ggtree(mrbayes_tree_protein_rooted, layout = "rectangular", size = 0.3) %<+% tip_df_mrbayes +
  geom_tiplab(aes(color = group), size = 2, offset = 0.002, align = TRUE, linesize = 0.1) +
  geom_tippoint(aes(color = group), size = 1.2) +
  geom_text2(aes(label = round(as.numeric(label), 2), subset = !isTip), size = 1.8, hjust = 1.2, vjust = -0.4) +
  scale_color_manual(values = group_colors) +
  theme_tree2() +
  labs(title = "Parvovirus DNA - MrBayes consensus protein", color = "Group") +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 8)) +
  coord_cartesian(clip = "off") +
  hexpand(0.5, direction = 1)

## with probability

mrbayes_tree_protein <- read.beast("/home/shadi/Desktop/ParvoProject/Parvo_sanger_full_seq/mrbayes_output/mrbayes_protein_MSA_FINAL_SVRNA_v3.nex.con.tre")
mrbayes_data <- as_tibble(mrbayes_tree_protein)[, c("node", "prob")]
class(mrbayes_data) <- c("tbl_df", "tbl", "data.frame")
n_tips <- length(mrbayes_tree_protein@phylo$tip.label)
print(mrbayes_data[mrbayes_data$prob > 0.95, ], n = 30)
nrow(mrbayes_data[mrbayes_data$prob > 0.95, ])
nrow(mrbayes_data[, ])

## MrBayes trees is now built on FINAL_SVRNA_v3 with 143 seq, so lets remove FPV_JO24_Sample_11 manually 

seq_to_remove <- "FPV_JO24_Sample_11"
mrbayes_tree_protein_noJO24 <- drop.tip(mrbayes_tree_protein, seq_to_remove)

length(mrbayes_tree_protein@phylo$tip.label)
length(mrbayes_tree_protein_noJO24@phylo$tip.label)
seq_to_remove %in% mrbayes_tree_protein_noJO24@phylo$tip.label

mrbayes_tree_protein_noJO24_rooted_oldestFPV <- root(as.phylo(mrbayes_tree_protein_noJO24), outgroup = "EU659112.1_ref_FPV_oldest", resolve.root = TRUE)
mrbayes_tree_protein_noJO24_rooted_oldestFPV <- treedata(phylo = mrbayes_tree_protein_noJO24_rooted_oldestFPV, data = mrbayes_data)

ggtreeplot_mrbayes_tree_protein_noJO24_rooted_oldestFPV <- ggtree(mrbayes_tree_protein_noJO24_rooted_oldestFPV, layout = "rectangular", size = 0.3) %<+% tip_df_mrbayes +
  geom_tiplab(aes(color = group), size = 2, offset = 0.002, align = TRUE, linesize = 0.1) +
  geom_tippoint(aes(color = group), size = 1.2) +
  geom_text2(aes(label = ifelse(!isTip & as.numeric(prob) >= 0, round(as.numeric(prob), 2), NA)), size = 1.8, hjust = 1.2, vjust = -0.4) +
  scale_color_manual(values = group_colors) +
  theme_tree2() +
  labs(
    # title = "Parvovirus protein - MrBayes consensus protein",
    color = "Group") +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 8)) +
  coord_cartesian(clip = "off") +
  hexpand(0.5, direction = 1)

parsimony(protein_NJ_Muscle_rooted_oldestFPV, protein_phy_Muscle)
parsimony(protein_ML_Muscle_tree_rooted_oldestFPV, protein_phy_Muscle)
# parsimony(mrbayes_tree_protein_noJO24_rooted_oldestFPV, protein_phy_Muscle)
parsimony(as.phylo(mrbayes_tree_protein_noJO24_rooted_oldestFPV), protein_phy_Muscle)

## filtration

tips_to_keep

# tips_to_drop <- setdiff(mrbayes_tree_protein_noJO24$tip.label, tips_to_keep)
tips_to_drop <- setdiff(as.phylo(mrbayes_tree_protein_noJO24)$tip.label, tips_to_keep)
mrbayes_tree_protein_noJO24_filtered <- drop.tip(mrbayes_tree_protein_noJO24, tips_to_drop)
mrbayes_tree_protein_noJO24_filtered_rooted_oldestFPV <- root(mrbayes_tree_protein_noJO24_filtered, outgroup = "EU659112.1_ref_FPV_oldest", resolve.root = TRUE)

# length(mrbayes_tree_protein_noJO24$tip.label)
# length(mrbayes_tree_protein_noJO24_filtered_rooted_oldestFPV$tip.label)
length(as.phylo(mrbayes_tree_protein_noJO24)$tip.label)
length(as.phylo(mrbayes_tree_protein_noJO24_filtered)$tip.label)

ggtree(mrbayes_tree_protein_noJO24_filtered_rooted_oldestFPV, layout = "rectangular", size = 0.3) %<+% tip_df_mrbayes +
  geom_tiplab(aes(color = group), size = 2, offset = 0.002, align = TRUE, linesize = 0.1) +
  geom_tippoint(aes(color = group), size = 1.2) +
  geom_text2(aes(label = label, subset = !isTip), size = 1.8, hjust = 1.2, vjust = -0.4) +
  scale_color_manual(values = group_colors) +
  theme_tree2() +
  labs(title = "Parvovirus DNA - MrBayes consensus", color = "Group") +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 8)) +
  coord_cartesian(clip = "off") +                  
  hexpand(0.5, direction = 1)

## Model selection for MrBayes - DNA & Protein

## DNA
dna_model_test <- modelTest(dna_phy_Muscle_trimmed,
                            model = c("JC","F81","K80","HKY","SYM","GTR"),
                            G = TRUE, I = TRUE)
dna_model_test[order(dna_model_test$BIC), ]   # lowest BIC = best fit

## Protein
protein_model_test <- modelTest(protein_phy_Muscle,
                                model = c("JTT","WAG","LG","Dayhoff"),
                                G = TRUE, I = TRUE)
protein_model_test[order(protein_model_test$BIC), ]

# save ggtree plots ------------------------------------------------------------

## save all ggtree plots

plot_dir <- file.path(sys_dir, "Parvo_sanger_full_seq", "Rplots/paper")
plot_names <- c(
  "ggtreeplot_dna_ML_Muscle_trimmed_tree_filtered_rooted_oldestFPV",
  "ggtreeplot_dna_ML_Muscle_trimmed_tree_rooted_oldestFPV",    
  "ggtreeplot_protein_ML_Muscle_tree_filtered_rooted_oldestFPV",                  
  "ggtreeplot_protein_ML_Muscle_tree_rooted_oldestFPV",
  "ggtreeplot_mrbayes_tree_dna_noJO24_rooted_oldestFPV",            
  "ggtreeplot_mrbayes_tree_protein_noJO24_rooted_oldestFPV" 
  )

# plot_dir <- file.path(sys_dir, "Parvo_sanger_full_seq", "Rplots/FINAL_SVRNA_v3")
# plot_names <- ls(pattern = "^ggtreeplot_")

# plot_dir <- file.path(sys_dir, "Parvo_sanger_full_seq", "Rplots/FINAL_SVRNA_v2/MrBayes_HapNet")
# plot_names <- ls(pattern = "^ggtreeplot_mrbayes")

for (i in plot_names) {
  
  ggsave(
    filename = file.path(plot_dir, paste0(sub("^ggtreeplot_", "", i), ".png")),
    plot = get(i),
    width = 12,
    height = 10,
    units = "in",
    dpi = 600,
    bg = "white",
    limitsize = FALSE)
  
  ggsave(
    filename = file.path(plot_dir, paste0(sub("^ggtreeplot_", "", i), ".pdf")),
    plot = get(i),
    width = 12,
    height = 10,
    units = "in",
    dpi = 600,
    bg = "white",
    limitsize = FALSE)
  
}

# ggsave(
#   filename = file.path(plot_dir, paste0(sub("^ggtreeplot_", "", i), "noJoFPV_v2.png")),
#   plot = ggtreeplot_protein_ML_Muscle_tree_filtered_rooted_oldestFPV,
#   width = 12,
#   height = 10,
#   units = "in",
#   dpi = 600,
#   bg = "white",
#   limitsize = FALSE)

# PCA + π - DNA & Protein analysis --------------------------------------------

## MSA with the outgroup

## PCA

## DNA

length(dna_alignment_Muscle_trimmed)

dna_distance <- dist.dna(as.DNAbin(dna_alignment_Muscle_trimmed), model = "F81") ## check model
dna_pcoa <- cmdscale(dna_distance, k = 2, eig = TRUE)
dna_var_explained <- round(dna_pcoa$eig / sum(dna_pcoa$eig[dna_pcoa$eig > 0]) * 100, 1)
dna_pcoa_df <- data.frame(
  PC1 = dna_pcoa$points[,1],
  PC2 = dna_pcoa$points[,2],
  Sample = rownames(dna_pcoa$points))
dna_pcoa_df <- merge(dna_pcoa_df, tip_df, by.x = "Sample", by.y = "label", all.x = TRUE)

pca_dna_v142 <- ggplot(dna_pcoa_df, aes(PC1, PC2, color = group, label = Sample)) +
  geom_point(size = 3) +
  # geom_text_repel(size = 3, max.overlaps = Inf) +
  scale_color_manual(values = group_colors) +
  labs(
    # title = "DNA PCoA",
    x = paste0("PC1 (", dna_var_explained[1], "%)"),
    y = paste0("PC2 (", dna_var_explained[2], "%)"),
    color = "Group") +
  theme_bw()

## Protein

protein_alignment_Muscle <- AAStringSet(protein_alignment_Muscle)
length(protein_alignment_Muscle)

protein_phy_Muscle <- phyDat(as.matrix(protein_alignment_Muscle), type = "AA")
protein_distance <- dist.ml(protein_phy_Muscle,model = "JTT")                    ## check model
protein_pcoa <- cmdscale(protein_distance, k = 2, eig = TRUE)
protein_var_explained <- round(protein_pcoa$eig / sum(protein_pcoa$eig[protein_pcoa$eig > 0]) * 100, 1)
protein_pcoa_df <- data.frame(
  PC1 = protein_pcoa$points[,1],
  PC2 = protein_pcoa$points[,2],
  Sample = rownames(protein_pcoa$points))
protein_pcoa_df <- merge(protein_pcoa_df, tip_df, by.x = "Sample", by.y = "label", all.x = TRUE)

pca_protein_v142 <- ggplot(protein_pcoa_df, aes(PC1, PC2, color = group, label = Sample)) +
  geom_point(size = 3) +
  # geom_text_repel(size = 3, max.overlaps = Inf) +
  scale_color_manual(values = group_colors) +
  labs(
    # title = "Protein PCoA",
    x = paste0("PC1 (", protein_var_explained[1], "%)"),
    y = paste0("PC2 (", protein_var_explained[2], "%)"),
    color = "Group") +
  theme_bw()

ggsave(filename = file.path(sys_dir, "Parvo_sanger_full_seq/Rplots/paper", "pca_dna_v142.png"),
       plot = pca_dna_v142, width = 8, height = 6, units = "in", dpi = 600, bg = "white", limitsize = FALSE)
ggsave(filename = file.path(sys_dir, "Parvo_sanger_full_seq/Rplots/paper", "pca_protein_v142.png"),
       plot = pca_protein_v142, width = 8, height = 6, units = "in", dpi = 600, bg = "white", limitsize = FALSE)

## Nucleotide diversity (π)

## π is the avg # of nucleotide differences per site between all pairs of sequences per group

dna_alignment_Muscle_trimmed

dna_aln_bins <- as.DNAbin(dna_alignment_Muscle_trimmed)
group_df <- tip_df[match(names(dna_aln_bins), tip_df$label), ]
group <- group_df$group

dna_aln_bins <- as.matrix(as.DNAbin(dna_alignment_Muscle_trimmed))
is.matrix(dna_aln_bins) 
dim(dna_aln_bins)

nuc.div(dna_aln_bins)

nucleotide_diversity <- data.frame(
  Group = names(table(group)),
  N = as.numeric(table(group)))

nucleotide_diversity$Pi <- sapply(nucleotide_diversity$Group,
                                  function(x) {round(nuc.div(dna_aln_bins[group == x, , drop = FALSE]), 6)})
nucleotide_diversity
# Group  N       Pi
# 1       Jo_Canine  9 0.002495
# 2       Jo_Feline  9 0.002975
# 3         Jo_Lion  1      NaN
# 4       Reference  6 0.007521
# 5 Regional_Canine 57 0.008239
# 6 Regional_Feline 56 0.004777
# 7         Vaccine  4 0.006921

# Variant calling - Protein analysis --------------------------------------

## calling variants

## using 
# protein_sequences_mod
# protein_alignment_Muscle

dna_names_mod <- names(protein_sequences_mod)
SamplesVacRef_seq <- dna_names_mod[grepl("sample|vaccine|ref", dna_names_mod, ignore.case = TRUE)]

SamplesVacRef_seq_order <- c(
  "FJ222821.1_ref_CPV2C" ,
  "M38246.1_ref_FPV", "EU659112.1_ref_FPV_oldest",
  "M38245.1_ref_CPV2", "M24003.1_ref_CPV2A", "M74849.1_ref_CPV2B",
  "Sample47_Canine"      ,"Sample68_Canine", "Sample72_Canine",
  "Sample76_Canine"   ,"Sample77_Canine",
  "Sample8_Canine"  , "sample51_canine","Sample57_Canine",  "Sample18_Canine",  
  "sample12_feline","Sample25_Feline","Sample27_Feline",  
  "Sample38_Feline","Sample44_Feline"        ,  "Sample69_Feline" , "Sample70_Feline" ,
  "Sample7_Feline", "Sample67_Feline",
  "Sample_Lion",
  # "FPV_JO24_Sample_11",
  "vaccineCanine_MSD" , "vaccineCanine_Zoetis" , "vaccineFeline_MSD"  ,  "vaccineFeline_Zoetis"
  )

protein_alignment_Muscle_ss <- AAStringSet(protein_alignment_Muscle)
protein_phy_Muscle_ss <- phyDat(as.matrix(protein_alignment_Muscle_ss), type = "AA")

protein_sequences_SamplesVacRef <- protein_sequences_mod[names(protein_sequences_mod) %in% SamplesVacRef_seq]
protein_alignment_Muscle_ss_SamplesVacRef <- protein_alignment_Muscle_ss[
  names(protein_alignment_Muscle_ss) %in% SamplesVacRef_seq,]

## AA mapping 

aa_3to1 <- c(Ala="A", Arg="R", Asn="N", Asp="D", Cys="C", Glu="E", Gln="Q", Gly="G",
             His="H", Ile="I", Leu="L", Lys="K", Met="M", Phe="F", Pro="P", Ser="S",
             Thr="T", Trp="W", Tyr="Y", Val="V")
aa_1to3 <- setNames(names(aa_3to1), aa_3to1)

# ref_FPV <- "M38246.1_ref_FPV"           # set reference
## you cannot use FPV ref as it is the full VP1 and VP2 region, so the protein pos are for VP1; Use CPV2c as it is shorter

ref_CPV2C <- "FJ222821.1_ref_CPV2C"     # set reference

names(protein_alignment_Muscle_ss_SamplesVacRef)
protein_alignment_Muscle_ss_SamplesVacRef_mat <- as.matrix(as(protein_alignment_Muscle_ss_SamplesVacRef, "AAStringSet"))
# rows = seqs, cols = alignment positions

## map alignment columns -> ungapped reference residue numbers for publication
## since FJ222821.1_ref_CPV2C is shorter than other refs (VP2 vs VP1), count VP1 as - pos

ref_seq_CPV2C <- protein_alignment_Muscle_ss_SamplesVacRef_mat[ref_CPV2C, ]
ref_nongap_cols_CPV2C <- which(ref_seq_CPV2C != "-")
ref_start_col_CPV2C <- min(ref_nongap_cols_CPV2C)
# positions after ref start pos
pos_cols <- ref_nongap_cols_CPV2C
pos_pos  <- seq_along(pos_cols)
# positions before ref start pos
n_before <- ref_start_col_CPV2C - 1
if (n_before > 0) {
  neg_cols <- (ref_start_col_CPV2C - 1):1
  neg_pos  <- -seq_len(n_before)
} else {
  neg_cols <- integer(0)
  neg_pos  <- integer(0)
}
col_to_refpos_CPV2C <- setNames(c(neg_pos, pos_pos), c(neg_cols, pos_cols))

## variant calling vs a reference (aka, CPV2c)

ProMSA_SamplesVacRef_Var_ref_CPV2C <- lapply(rownames(protein_alignment_Muscle_ss_SamplesVacRef_mat), function(sn) {
  q <- protein_alignment_Muscle_ss_SamplesVacRef_mat[sn, ]
  diffs <- which(q != ref_seq_CPV2C & q != "-" & ref_seq_CPV2C != "-")
  if (length(diffs) == 0) return(NULL)
  refpos <- col_to_refpos_CPV2C[as.character(diffs)]
  data.frame(seq = sn, pos = refpos, ref_aa = ref_seq_CPV2C[diffs], alt_aa = q[diffs],
             mutation = paste0(ref_seq_CPV2C[diffs], refpos, q[diffs]))
})

ProMSA_SamplesVacRef_Var_ref_CPV2C_df <- do.call(rbind, ProMSA_SamplesVacRef_Var_ref_CPV2C)

## loop through all 6 refs and get all divergence in a variant table

ref_names <- grep("_ref_", rownames(protein_alignment_Muscle_ss_SamplesVacRef_mat), value = TRUE)
# SamplesVac_seq <- setdiff(rownames(protein_alignment_Muscle_ss_SamplesVacRef_mat), ref_names)
SamplesVacRef_seq

all_cols_CPV2C <- as.numeric(names(col_to_refpos_CPV2C))

ProMSA_SamplesVacRef_Var_allref_CPV2C_list <- lapply(SamplesVacRef_seq, function(sn) {
  sample_row <- protein_alignment_Muscle_ss_SamplesVacRef_mat[sn, all_cols_CPV2C]
  ref_matrix <- protein_alignment_Muscle_ss_SamplesVacRef_mat[ref_names, all_cols_CPV2C, drop = FALSE]
  
  disagrees <- sapply(seq_along(all_cols_CPV2C), function(i) {
    sv <- sample_row[i]
    rv <- ref_matrix[, i]
    sv != "-" && any(sv != rv & rv != "-")
  })
  
  if (!any(disagrees)) return(NULL)
  
  hit_cols <- all_cols_CPV2C[disagrees]
  hit_refpos <- col_to_refpos_CPV2C[as.character(hit_cols)]
  
  out <- data.frame(seq = sn, pos = hit_refpos, sample_aa = sample_row[disagrees])
  ref_vals <- t(ref_matrix[, disagrees, drop = FALSE])
  colnames(ref_vals) <- paste0("aa_", ref_names)
  cbind(out, ref_vals)
})

ProMSA_SamplesVacRef_Var_allref_CPV2C_df <- do.call(rbind, ProMSA_SamplesVacRef_Var_allref_CPV2C_list)
ProMSA_SamplesVacRef_Var_allref_CPV2C_df <- ProMSA_SamplesVacRef_Var_allref_CPV2C_df[order(ProMSA_SamplesVacRef_Var_allref_CPV2C_df$seq, ProMSA_SamplesVacRef_Var_allref_CPV2C_df$pos), ]
rownames(ProMSA_SamplesVacRef_Var_allref_CPV2C_df) <- NULL

ProMSA_SamplesVacRef_Var_allref_CPV2C_df

## pivoting to a table

all_cols_CPV2C <- as.numeric(names(col_to_refpos_CPV2C))
all_refpos_CPV2C <- col_to_refpos_CPV2C[as.character(all_cols_CPV2C)]

var_cols_final_CPV2C <- sort(unique(ProMSA_SamplesVacRef_Var_allref_CPV2C_df$pos))
var_cols_alncol_CPV2C <- all_cols_CPV2C[match(var_cols_final_CPV2C, all_refpos_CPV2C)]

var_table_mat_CPV2C <- sapply(var_cols_alncol_CPV2C, function(c) {
  fpv_aa <- ref_seq_CPV2C[c]
  col_vals <- protein_alignment_Muscle_ss_SamplesVacRef_mat[, c]
  out <- ifelse(col_vals == "-", "-", ifelse(col_vals == fpv_aa, ".", col_vals))
  out[ref_CPV2C] <- fpv_aa
  out
})

colnames(var_table_mat_CPV2C) <- var_cols_final_CPV2C
var_table_df_CPV2C <- data.frame(Strain = rownames(protein_alignment_Muscle_ss_SamplesVacRef_mat), var_table_mat_CPV2C, check.names = FALSE)
rownames(var_table_df_CPV2C) <- NULL

var_table_df_CPV2C <- var_table_df_CPV2C[match(SamplesVacRef_seq_order, var_table_df_CPV2C$Strain), ]
rownames(var_table_df_CPV2C) <- NULL
var_table_df_CPV2C

# "." = matches CPV2c at this position
# "-" = no sequence here (aka, missing data)

## AA mapping 

ProMSA_SamplesVacRef_Var_ref_CPV2C_df$ref_aa_3 <- aa_1to3[ProMSA_SamplesVacRef_Var_ref_CPV2C_df$ref_aa]
ProMSA_SamplesVacRef_Var_ref_CPV2C_df$alt_aa_3 <- aa_1to3[ProMSA_SamplesVacRef_Var_ref_CPV2C_df$alt_aa]
ProMSA_SamplesVacRef_Var_ref_CPV2C_df$mutation_3 <- paste0(
  ProMSA_SamplesVacRef_Var_ref_CPV2C_df$ref_aa_3,
  ProMSA_SamplesVacRef_Var_ref_CPV2C_df$pos,
  ProMSA_SamplesVacRef_Var_ref_CPV2C_df$alt_aa_3)

var_table_df_CPV2C_3letter <- var_table_df_CPV2C
cols_to_convert <- setdiff(names(var_table_df_CPV2C_3letter), "Strain")

var_table_df_CPV2C_3letter[cols_to_convert] <- lapply(var_table_df_CPV2C_3letter[cols_to_convert], function(col) {
  ifelse(col %in% c(".", "-"), col, aa_1to3[col])
})

write.xlsx(ProMSA_SamplesVacRef_Var_allref_CPV2C_df,
           file.path(sys_dir, "Parvo_sanger_full_seq/ProMSA_SamplesVacRef_Var_ref_CPV2C_v3_v142.xlsx"))
write.xlsx(var_table_df_CPV2C_3letter,
           file.path(sys_dir, "Parvo_sanger_full_seq/VP2_aa_var_table_3letter_ref_CPV2C_v3_v142.xlsx"))

# Haplotype network -------------------------------------------------------------

# ## haplotypes from your trimmed VP2 alignment
# dna_bin_hap <- as.DNAbin(dna_alignment_Muscle_modified_trimmed)
# hap <- pegas::haplotype(dna_bin_hap)
# hap
# 
# ## map each sample to its haplotype (tutorial's stack/setNames trick)
# hapInfo <- stack(setNames(attr(hap, "index"), rownames(hap)))
# names(hapInfo) <- c("index", "haplotype")
# hapInfo$sample <- labels(dna_bin_hap)[hapInfo$index]
# 
# ## merge with your existing group metadata
# merged <- merge(hapInfo, tip_df, by.x = "sample", by.y = "label")
# head(merged)
# table(merged$haplotype)
# 
# ## build network
# net <- pegas::haploNet(hap)
# net
# 
# ## plain plot first
# plot(net, size = attr(net, "freq"), fast = TRUE)
# 
# ## pie matrix: haplotype x group
# pie <- table(merged$haplotype, merged$group)
# pie
# 
# ## colored plot using your existing group_colors palette
# plot(net, size = attr(net, "freq"), pie = pie,
#      bg = group_colors[colnames(pie)],
#      labels = FALSE,
#      fast = TRUE,
#      cex = 0.6)
# legend("topright", colnames(pie), col = group_colors[colnames(pie)], pch = 19, ncol = 2)
# 
# dev.off()
# 
# plot(net, size = attr(net, "freq"), pie = pie,
#      bg = group_colors[colnames(pie)],
#      labels = FALSE,
#      fast = TRUE,
#      cex = 0.6,
#      show.mutation = 2)
# legend("topright", colnames(pie), col = group_colors[colnames(pie)], pch = 19, ncol = 2)

## build a hap net for your samples only

mysamp_vec_ref_to_keep <- c(
  "Sample57_Canine",
  "Sample68_Canine",
  "Sample47_Canine",
  "Sample77_Canine",
  "Sample76_Canine",
  "Sample8_Canine",
  "Sample18_Canine",
  "Sample72_Canine",
  "sample51_canine",
  "FJ222821.1_ref_CPV2C",
  "M74849.1_ref_CPV2B",
  "M24003.1_ref_CPV2A",
  "M38245.1_ref_CPV2",
  "vaccineCanine_MSD",
  "vaccineCanine_Zoetis",
  "Sample70_Feline",
  "Sample38_Feline",
  "Sample_Lion",
  "Sample69_Feline",
  "Sample27_Feline",
  "Sample25_Feline",
  "Sample7_Feline",
  "Sample44_Feline",
  # "FPV_JO24_Sample_11",
  "sample12_feline",
  "Sample67_Feline",
  "M38246.1_ref_FPV",
  "vaccineFeline_Zoetis",
  "vaccineFeline_MSD",
  "EU659112.1_ref_FPV_oldest"
)

length(mysamp_vec_ref_to_keep)

dna_sequences_modified <- dna_sequences_mod[names(dna_sequences_mod) %in% mysamp_vec_ref_to_keep]
protein_sequences_modified <- protein_sequences_mod[names(protein_sequences_mod) %in% mysamp_vec_ref_to_keep]

length(dna_sequences_modified)

dna_sequences_modified                     
length(dna_sequences_modified)              
width(dna_sequences_modified)        
names(dna_sequences_modified)       
dna_sequences_modified[[1]]               
as.character(dna_sequences_modified[1])     
letterFrequency(dna_sequences_modified, letters = "GC", as.prob = TRUE)

dna_alignment_Muscle_modified <- msa(dna_sequences_modified, method = "Muscle")
dna_phy_Muscle_modified <- phyDat(as.matrix(dna_alignment_Muscle_modified), type = "DNA")

dna_alignment_Muscle_modified_trimmed <- DNAStringSet(dna_alignment_Muscle_modified)
ref_FPV_name <- grep("M38246.1_ref_FPV", names(dna_alignment_Muscle_modified_trimmed), value = TRUE)
ref_FPV_seq <- dna_alignment_Muscle_modified_trimmed[[ref_FPV_name]]
ref_FPV_chars <- strsplit(as.character(ref_FPV_seq), "")[[1]]
non_gap_positions <- which(ref_FPV_chars != "-")

start_col <- non_gap_positions[2655]
end_col <- non_gap_positions[4585]

dna_alignment_Muscle_modified_trimmed <- subseq(dna_alignment_Muscle_modified_trimmed, start = start_col, end = end_col)
width(dna_alignment_Muscle_modified_trimmed)[1]
dna_phy_Muscle_modified_trimmed <- phyDat(as.matrix(dna_alignment_Muscle_modified_trimmed), type = "DNA")

## haplotypes from your trimmed VP2 alignment
dna_bin_hap <- as.DNAbin(dna_alignment_Muscle_modified_trimmed)
hap <- pegas::haplotype(dna_bin_hap)
hap
# Number of haplotypes: 23 
# Sequence length: 1931 

## map each sample to its haplotype (tutorial's stack/setNames trick)
hapInfo <- stack(setNames(attr(hap, "index"), rownames(hap)))
names(hapInfo) <- c("index", "haplotype")
hapInfo$sample <- labels(dna_bin_hap)[hapInfo$index]

## merge with your existing group metadata
merged <- merge(hapInfo, tip_df, by.x = "sample", by.y = "label")
head(merged)
table(merged$haplotype)
table(merged$group)
length(unique(merged$haplotype))

## build network
net <- pegas::haploNet(hap)
net

## plain plot first
plot(net, size = attr(net, "freq"), fast = TRUE)

## pie matrix: haplotype x group
pie <- table(merged$haplotype, merged$group)
pie

## colored plot using your existing group_colors palette
plot(net, size = attr(net, "freq"), pie = pie,
     bg = group_colors[colnames(pie)],
     labels = FALSE,
     fast = TRUE,
     cex = 0.6)
legend("topright", colnames(pie), col = group_colors[colnames(pie)], pch = 19, ncol = 2)

# save it
dev.off()
# png("/home/shadi/Desktop/ParvoProject/Parvo_sanger_full_seq/Rplots/FINAL_SVRNA_v3/MrBayes_HapNet/haplotype_network.png",
#     width = 8, height = 8, units = "in", res = 300)
png("/home/shadi/Desktop/ParvoProject/Parvo_sanger_full_seq/Rplots/paper/haplotype_network.png",
    width = 8, height = 6, units = "in", res = 600)
plot(net, size = attr(net, "freq"), pie = pie,
     bg = group_colors[colnames(pie)],
     labels = FALSE,
     fast = TRUE,
     cex = 0.6)
legend("topleft", colnames(pie), col = group_colors[colnames(pie)], pch = 19, ncol = 2)
dev.off()

# Positive selection - DNA ------------------------------------------------

## we need the trimmed coding seq of VP2 with no stop codons - dna_alignment_Muscle_trimmed_CPV2C

width(dna_alignment_Muscle_trimmed_CPV2C)/3 ## aka Parvovirus_DNA_alignment_Muscle_trimmed_CPV2C_FINAL_SVRNA_v3_v142.fasta

## filter for the 65 rep seq not the 142 to run sel scan using Data Monkey

## FINAL_SVRNA_v3_v142

seq_to_keep <- c(
  "Sample68_Canine",
  "OM100697.1_Canine_parvovirus_2_Egypt",
  "Sample47_Canine",
  "Sample8_Canine",
  "Sample18_Canine",
  "OM100699.1_Canine_parvovirus_2_Egypt",
  # "OM100701.1_Canine_parvovirus_2_Egypt",
  # "OM100698.1_Canine_parvovirus_2_Egypt",
  "OM100696.1_Canine_parvovirus_2_Egypt",
  "MG013488.1_Canine_parvovirus_2c_China",
  "MF510157.1_Canine_parvovirus_2c_Italy",
  "Sample57_Canine",
  "MZ056886.1_Canine_parvovirus_2_Egypt",
  # "MZ056884.1_Canine_parvovirus_2_Egypt",
  "Sample77_Canine",
  "Sample76_Canine",
  "MZ056885.1_Canine_parvovirus_2_Egypt",
  # "MZ056889.1_Canine_parvovirus_2_Egypt",
  "sample51_canine",
  "Sample72_Canine",
  "OR667801.1_Canine_parvovirus_2_Iraq",
  # "OR667800.1_Canine_parvovirus_2_Iraq",
  # "OR667802.1_Canine_parvovirus_2_Iraq",
  # "OR451707.1_Canine_parvovirus_2_Iraq",
  # "OR667804.1_Canine_parvovirus_2_Iraq",
  # "OR667803.1_Canine_parvovirus_2_Iraq",
  "MW653250.1_Canine_parvovirus_2a_Iran",
  "OL330980.1_Canine_parvovirus_2_Iran",
  # "OL330978.1_Canine_parvovirus_2_Iran",
  "MW653249.1_Canine_parvovirus_2a_Iran",
  "OM100702.1_Canine_parvovirus_2_Egypt",
  # "OM100700.1_Canine_parvovirus_2_Egypt",
  "OL330979.1_Canine_parvovirus_2_Iran",
  "MZ056892.1_Canine_parvovirus_2_Egypt",
  # "MZ056883.1_Canine_parvovirus_2_Egypt",
  # "MZ056891.1_Canine_parvovirus_2_Egypt",
  "OM721655.1_Canine_parvovirus_2b_Turkey",
  "MZ056890.1_Canine_parvovirus_2_Egypt",
  # "MZ056882.1_Canine_parvovirus_2_Egypt",
  # "MZ056887.1_Canine_parvovirus_2_Egypt",
  # "MZ056888.1_Canine_parvovirus_2_Egypt",
  # "MZ056881.1_Canine_parvovirus_2_Egypt",
  "MW539053.1_Canine_parvovirus_2b_Turkey",
  # "OM721656.1_Canine_parvovirus_2b_Turkey",
  "OL330977.1_Canine_parvovirus_2_Iran",
  "KX268117.1_Canine_parvovirus_2_Turkey",
  # "KX268105.1_Canine_parvovirus_2_Turkey",
  # "KX268107.1_Canine_parvovirus_2_Turkey",
  # "KX268118.1_Canine_parvovirus_2_Turkey",
  # "KX268106.1_Canine_parvovirus_2_Turkey",
  # "KX268114.1_Canine_parvovirus_2_Turkey",
  # "KX268110.1_Canine_parvovirus_2_Turkey",
  # "KX268108.1_Canine_parvovirus_2_Turkey",
  "MW653251.1_Canine_parvovirus_2b_Iran",
  "MW653253.1_Canine_parvovirus_2c_Iran",
  # "KX268109.1_Canine_parvovirus_2_Turkey",
  "FJ222821.1_ref_CPV2C",
  "MW653256.1_Canine_parvovirus_2a_Iran",
  "KX268115.1_Canine_parvovirus_2_Turkey",
  "AY742933.1_recent_Canine_parvovirus_2a_NewZealand",
  "KX268116.1_Canine_parvovirus_2_Turkey",
  # "KX268111.1_Canine_parvovirus_2_Turkey",
  # "KX268113.1_Canine_parvovirus_2_Turkey",
  # "KX268112.1_Canine_parvovirus_2_Turkey",
  "MW653252.1_Canine_parvovirus_2b_Iran",
  "MW653248.1_Canine_parvovirus_2a_Iran",
  "M74849.1_ref_CPV2B",
  "M24003.1_ref_CPV2A",
  "M38245.1_ref_CPV2",
  "vaccineCanine_MSD",
  "KX268119.1_Canine_parvovirus_2_Turkey",
  "vaccineCanine_Zoetis",
  "Sample70_Feline",
  "Sample38_Feline",
  "Sample69_Feline",
  "Sample_Lion",
  "Sample25_Feline",
  "Sample27_Feline",
  "Sample7_Feline",
  "Sample44_Feline",
  "PV521947.1_Feline_parvovirus_Egypt",
  # "PP663048.2_Feline_parvovirus_Egypt",
  "OM638043.1_Feline_parvovirus_Egypt",
  # "PV521964.1_Feline_parvovirus_Egypt",
  # "PP663074.2_Feline_parvovirus_Egypt",
  # "PV521948.1_Feline_parvovirus_Egypt",
  # "PV808490.1_Feline_parvovirus_Egypt",
  # "PV808489.1_Feline_parvovirus_Egypt",
  # "PV521962.1_Feline_parvovirus_Egypt",
  # "PV521954.1_Feline_parvovirus_Egypt",
  # "PV521953.1_Feline_parvovirus_Egypt",
  # "PP663056.2_Feline_parvovirus_Egypt",
  # "PP663055.2_Feline_parvovirus_Egypt",
  # "PV808495.1_Feline_parvovirus_Egypt",
  # "PV521960.1_Feline_parvovirus_Egypt",
  # "PP663070.2_Feline_parvovirus_Egypt",
  # "PP663041.2_Feline_parvovirus_Egypt",
  # "PP663072.2_Feline_parvovirus_Egypt",
  # "PP663050.2_Feline_parvovirus_Egypt",
  "PP663073.2_Feline_parvovirus_Egypt",
  # "PV521961.1_Feline_parvovirus_Egypt",
  # "PV521949.1_Feline_parvovirus_Egypt",
  # "PP663045.2_Feline_parvovirus_Egypt",
  # "PV521942.1_Feline_parvovirus_Egypt",
  # "PV521963.1_Feline_parvovirus_Egypt",
  # "PP663071.2_Feline_parvovirus_Egypt",
  # "PP663044.2_Feline_parvovirus_Egypt",
  # "PV521955.1_Feline_parvovirus_Egypt",
  # "PP663059.2_Feline_parvovirus_Egypt",
  # "PV808491.1_Feline_parvovirus_Egypt",
  # "PV521951.1_Feline_parvovirus_Egypt",
  # "PV808488.1_Feline_parvovirus_Egypt",
  # "PV521944.1_Feline_parvovirus_Egypt",
  # "PV521943.1_Feline_parvovirus_Egypt",
  # "PP663042.2_Feline_parvovirus_Egypt",
  # "PV521965.1_Feline_parvovirus_Egypt",
  # "PV521945.1_Feline_parvovirus_Egypt",
  # "PV521959.1_Feline_parvovirus_Egypt",
  # "PV521958.1_Feline_parvovirus_Egypt",
  # "PV521957.1_Feline_parvovirus_Egypt",
  # "PV808492.1_Feline_parvovirus_Egypt",
  # "PV521952.1_Feline_parvovirus_Egypt",
  # "PV521956.1_Feline_parvovirus_Egypt",
  # "PP663061.2_Feline_parvovirus_Egypt",
  # "PV808496.1_Feline_parvovirus_Egypt",
  # "PP663062.2_Feline_parvovirus_Egypt",
  # "PV808493.1_Feline_parvovirus_Egypt",
  # "PV808486.1_Feline_parvovirus_Egypt",
  # "PP663067.2_Feline_parvovirus_Egypt",
  # "PV521946.1_Feline_parvovirus_Egypt",
  "PP663047.2_Feline_parvovirus_Egypt",
  "sample12_feline",
  "PV521950.1_Feline_parvovirus_Egypt",
  # "PP663052.2_Feline_parvovirus_Egypt",
  "Sample67_Feline",
  "KX434461.1_Feline_parvovirus_IZSSI_Italy",
  "KX900570.1_Feline_parvovirus_HH6_China",
  "M38246.1_ref_FPV",
  "KP081409.1_Feline_parvovirus_Iran",
  "vaccineFeline_Zoetis",
  "vaccineFeline_MSD",
  "EU659112.1_ref_FPV_oldest"
)

width(dna_alignment_Muscle_trimmed_CPV2C)/3
length(dna_alignment_Muscle_trimmed_CPV2C)

dna_alignment_Muscle_trimmed_CPV2C_repseq <- dna_alignment_Muscle_trimmed_CPV2C[names(dna_alignment_Muscle_trimmed_CPV2C) %in% seq_to_keep]

width(dna_alignment_Muscle_trimmed_CPV2C_repseq)/3
length(dna_alignment_Muscle_trimmed_CPV2C_repseq)
length(seq_to_keep)

writeXStringSet(DNAStringSet(dna_alignment_Muscle_trimmed_CPV2C_repseq),
                file.path(sys_dir, "Parvo_sanger_full_seq/aln/paper/Parvovirus_DNA_alignment_Muscle_trimmed_CPV2C_FINAL_SVRNA_v3_v142_repseq.fasta"))

## run Data Monkey on Parvovirus_DNA_alignment_Muscle_trimmed_CPV2C_FINAL_SVRNA_v3_v142_repseq.fasta

## SLAC

slac_path <- file.path(sys_dir, "Parvo_sanger_full_seq/Evo_Analysis_Datamonkey/SLAC","datamonkey-table.csv")

slac_res <- read.csv(slac_path, na.strings = "null") 
colnames (slac_res) <- c(
   "Partition","Site","ES","EN","S","N","P[S]","dS","dN","dN_dS",
   "P [dN/dS > 1]","P [dN/dS < 1]","Total branch length")
slac_res <- slac_res %>%
  mutate(
    Selection = case_when(
      `P [dN/dS > 1]` < 0.05 ~ "Positive_selection",
      `P [dN/dS < 1]` < 0.05 ~ "Purifying_selection",
      TRUE ~ "Not_significant"))

summary(slac_res$dN_dS)

slac_res_summary <- slac_res %>%
  count(Selection) %>%
  mutate(Percentage = 100 * n / sum(n))

slac_sum_plot <- ggplot(slac_res_summary, aes(x = Selection, y = n)) +
  geom_col() +
  geom_text(aes(label = paste0(n, " (", round(Percentage, 1), "%)")), vjust = -0.4) +
  labs(x = "Selection category", y = "Number of codon sites") +
  theme_classic()

ggplot(slac_res, aes(x = Site, y = dN_dS)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_segment(aes(xend = Site, y = 1, yend = dN_dS, color = Selection)) +
  geom_point(aes(color = Selection), size = 1.5) +
  scale_color_manual(values = c("Purifying_selection" = "firebrick", "Positive_selection" = "darkgreen", "Not_significant" = "grey")) +
  scale_x_continuous(breaks = seq(0, max(slac_res$Site, na.rm = TRUE), by = 50)) +
  labs(x = "Codon site", y = "Nonsynonymous-to-synonymous substitution rate ratio (dN/dS)") +
  theme_classic()

## FEL

fel_path <- file.path(sys_dir, "Parvo_sanger_full_seq/Evo_Analysis_Datamonkey/FEL","fel_res_table.txt")

fel_res <- read.csv(fel_path, sep = "\t")
fel_res <- fel_res %>%
  select(-X) %>%
  mutate(
    # Selection = case_when(
    #   p.value < 0.05 & beta > alpha ~ "Positive_selection",
    #   p.value < 0.05 & beta < alpha ~ "Purifying_selection",
    #   TRUE ~ "Not_significant"),
    beta_alpha = ifelse(alpha > 0, beta / alpha, NA))

summary(fel_res$beta_alpha)

table(fel_res$class)

fel_res_purifying_codons <- fel_res[fel_res$class == "Purifying",]$codon

fel_res_summary <- fel_res %>%
  count(class) %>%
  mutate(Percentage = 100 * n / sum(n))

fel_sum_plot <- ggplot(fel_res_summary, aes(x = class, y = n)) +
  geom_col() +
  geom_text(aes(label = paste0(n, " (", round(Percentage, 1), "%)")), vjust = -0.4) +
  labs(x = "Selection category", y = "Number of codon sites") +
  theme_classic()

ggplot(fel_res, aes(x = codon, y = alpha.beta)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_segment(aes(xend = codon, y = 1, yend = alpha.beta, color = class)) +
  geom_point(aes(color = class), size = 1.5) +
  scale_color_manual(values = c("Purifying" = "firebrick","Diversifying" = "darkgreen","Neutral" = "gray60", "Invariable" = "grey")) +
  scale_x_continuous(breaks = seq(0, max(fel_res$codon, na.rm = TRUE), by = 50)) +
  labs(x = "Codon site", y = "Nonsynonymous-to-synonymous substitution rate ratio (β/α)") +
  theme_classic()

ggplot((fel_res), aes(x = codon)) +
  geom_col(aes(y = beta, fill = class), width = 1) +
  geom_col(aes(y = -alpha, fill = class), width = 1) +
  geom_hline(yintercept = 0) +
  scale_x_continuous(breaks = seq(0, max(fel_res$codon, na.rm = TRUE), by = 50)) +
  scale_y_continuous(limits = c(-10, 10), breaks = seq(-10, 10, by = 2), labels = abs) +
  scale_fill_manual(values = c("Purifying" = "firebrick","Diversifying" = "darkgreen","Neutral" = "gray60","Invariable" = "grey")) +
  labs(x = "Codon site", y = "Rate estimate") +
  theme_classic()

## MEME

meme_path <- file.path(sys_dir, "Parvo_sanger_full_seq/Evo_Analysis_Datamonkey/MEME","meme_res_table.txt")

meme_res <- read.csv(meme_path, sep = "\t")
meme_res <- meme_res %>%
  select(-X) %>%
  mutate(
    `β.` = as.numeric(gsub(",", "", `β.`)),
    `β1` = as.numeric(gsub(",", "", `β1`)),
    `α` = as.numeric(gsub(",", "", `α`)))

sum(is.na(meme_res$α))

summary(meme_res$p.value)

table(meme_res$Class)

meme_res_summary <- meme_res %>%
  count(Class) %>%
  mutate(Percentage = 100 * n / sum(n))

meme_sum_plot <- ggplot(meme_res_summary, aes(x = Class, y = n)) +
  geom_col() +
  geom_text(aes(label = paste0(n, " (", round(Percentage, 1), "%)")), vjust = -0.4) +
  labs(x = "Selection category", y = "Number of codon sites") +
  theme_classic()

ggplot(meme_res, aes(x = Codon, y = `β.`)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_segment(aes(xend = Codon, y = 1, yend = `β.`, color = Class)) +
  geom_point(aes(color = Class), size = 1.5) +
  scale_color_manual(values = c("Purifying" = "firebrick","Diversifying" = "darkgreen","Neutral" = "gray60","Invariable" = "grey")) +
  scale_x_continuous(breaks = seq(0, max(meme_res$Codon, na.rm = TRUE), by = 50)) +
  labs(x = "Codon site", y = "Positive-selection rate (β+)") +
  theme_classic()

## FUBAR

fubar_path <- file.path(sys_dir, "Parvo_sanger_full_seq/Evo_Analysis_Datamonkey/FUBAR","datamonkey-table.csv")

fubar_res <- read.csv(fubar_path) %>%
  # select(-Partition) %>%
  rename(
    alpha = X.alpha.,
    beta = X.beta.,
    beta_alpha = X.beta...alpha.,
    prob_alpha_beta = Prob..alpha...beta..,
    prob_beta_alpha = Prob..alpha...beta...1,
    BayesFactor = BayesFactor..alpha...beta..) %>%
  mutate(
    Selection = case_when(
      prob_beta_alpha > 0.9 ~ "Diversifying",
      prob_alpha_beta > 0.9 ~ "Purifying",
      TRUE ~ "Neutral"))

summary(fubar_res$beta_alpha)

table(fubar_res$Selection)

fubar_res_summary <- fubar_res %>%
  count(Selection) %>%
  mutate(Percentage = 100 * n / sum(n))

fubar_sum_plot <- ggplot(fubar_res_summary, aes(x = Selection, y = n)) +
  geom_col() +
  geom_text(aes(label = paste0(n, " (", round(Percentage, 1), "%)")), vjust = -0.4) +
  labs(x = "Selection category", y = "Number of codon sites") +
  theme_classic()

ggplot(fubar_res, aes(x = Site, y = beta_alpha)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_segment(aes(xend = Site, y = 1, yend = beta_alpha, color = Selection)) +
  geom_point(aes(color = Selection), size = 1.5) +
  scale_color_manual(values = c("Purifying" = "firebrick","Diversifying" = "darkgreen","Neutral" = "gray60")) +
  scale_x_continuous(breaks = seq(0, max(fubar_res$Site, na.rm = TRUE), by = 50)) +
  labs(x = "Codon site", y = "Nonsynonymous-to-synonymous substitution rate ratio (β/α)") +
  theme_classic()

## combine slac + fel + meme + fubar

combined_sel_sum_plot <- (meme_sum_plot / fubar_sum_plot | fel_sum_plot / slac_sum_plot )

ggsave(filename = file.path(sys_dir, "Parvo_sanger_full_seq/Rplots/paper", "combined_selection_summary_plot.png"),
       plot = combined_sel_sum_plot, width = 12, height = 10, units = "in", dpi = 600, bg = "white", limitsize = FALSE)

slac_res_combine <- slac_res %>%
  select(Site, Selection) %>%
  rename(SLAC_sel = Selection)

fel_res_combine <- fel_res %>%
  select(codon, class) %>%
  rename(Site = codon, FEL_sel = class)

meme_res_combine <- meme_res %>%
  select(Codon, Class) %>%
  rename(Site = Codon, MEME_sel = Class)

fubar_res_combine <- fubar_res %>%
  select(Site, Selection) %>%
  rename(FUBAR_sel = Selection)

table(meme_res$Class)
table(fubar_res$Selection)
table(slac_res$Selection)
table(fel_res$class)

# meme_res$dN_dS <- meme_res$"β1"/meme_res$"α"
# 
# selection_res_combine <- meme_res %>%
#   left_join(
#     selection_res_mapping %>%
#       select(Site,SLAC_sel,FEL_sel,MEME_sel,combined_sel), 
#     by = c("Codon" = "Site"))
# 
# ggplot(selection_res_combine, aes(x = Codon, y = dN_dS)) +
#   geom_hline(yintercept = 1, linetype = "dashed") +
#   geom_segment(aes(xend = Codon, y = 1, yend = dN_dS, color = combined_sel)) +
#   geom_point(aes(color = combined_sel), size = 1.8) +
#   scale_color_manual(values = c(
#     "Purifying_FEL" = "firebrick2",
#     "Purifying_SLAC_FEL" = "firebrick4",
#     "Diversifying_MEME" = "darkgreen",
#     "Not_significant_SLAC_FEL_MEME" = "gray60")) +
#   scale_x_continuous(breaks = seq(0, max(selection_res_combine$Codon, na.rm = TRUE), by = 50)) +
#   labs(x = "Codon site", y = "Substitution rate ratio (dN/dS)", color = "Combined selection signal") +
#   theme_classic()

selection_res_mapping <- slac_res_combine %>%
  left_join(fel_res_combine, by = "Site") %>%
  left_join(meme_res_combine, by = "Site") %>%
  left_join(fubar_res_combine, by = "Site") %>%
  mutate(
    SLAC_sel = recode(
      SLAC_sel,
      "Positive_selection" = "Diversifying",
      "Purifying_selection" = "Purifying",
      "Not_significant" = "Not_significant"),
    FEL_sel = recode(
      FEL_sel,
      "Diversifying" = "Diversifying",
      "Purifying" = "Purifying",
      "Neutral" = "Not_significant",
      "Invariable" = "Not_significant"),
    MEME_sel = recode(
      MEME_sel,
      "Diversifying" = "Diversifying",
      "Purifying" = "Purifying",
      "Neutral" = "Not_significant",
      "Invariable" = "Not_significant"),
    FUBAR_sel = recode(
      FUBAR_sel,
      "Diversifying" = "Diversifying",
      "Purifying" = "Purifying",
      "Neutral" = "Not_significant"),
    combined_sel = case_when(
      SLAC_sel == "Not_significant" & FEL_sel == "Not_significant" & MEME_sel == "Not_significant" & FUBAR_sel != "Not_significant" ~ paste0(FUBAR_sel, "_FUBAR"),
      SLAC_sel == "Not_significant" & FEL_sel == "Not_significant" & MEME_sel != "Not_significant" & FUBAR_sel == "Not_significant" ~ paste0(MEME_sel, "_MEME"),
      SLAC_sel == "Not_significant" & FEL_sel != "Not_significant" & MEME_sel == "Not_significant" & FUBAR_sel == "Not_significant" ~ paste0(FEL_sel, "_FEL"),
      SLAC_sel != "Not_significant" & FEL_sel == "Not_significant" & MEME_sel == "Not_significant" & FUBAR_sel == "Not_significant" ~ paste0(SLAC_sel, "_SLAC"),
      SLAC_sel == FEL_sel & FEL_sel == MEME_sel & MEME_sel == FUBAR_sel ~ paste0(SLAC_sel, "_SLAC_FEL_MEME_FUBAR"),
      SLAC_sel == FEL_sel & FEL_sel == MEME_sel ~ paste0(SLAC_sel, "_SLAC_FEL_MEME"),
      SLAC_sel == FEL_sel & FEL_sel == FUBAR_sel ~ paste0(SLAC_sel, "_SLAC_FEL_FUBAR"),
      SLAC_sel == MEME_sel & MEME_sel == FUBAR_sel ~ paste0(SLAC_sel, "_SLAC_MEME_FUBAR"),
      FEL_sel == MEME_sel & MEME_sel == FUBAR_sel ~ paste0(FEL_sel, "_FEL_MEME_FUBAR"),
      MEME_sel == FUBAR_sel ~ paste0(MEME_sel, "_MEME_FUBAR"),
      SLAC_sel == FEL_sel ~ paste0(SLAC_sel, "_SLAC_FEL"),
      FEL_sel == FUBAR_sel ~ paste0(FEL_sel, "_FEL_FUBAR"),
      SLAC_sel == FUBAR_sel ~ paste0(SLAC_sel, "_SLAC_FUBAR"),
      SLAC_sel == MEME_sel ~ paste0(SLAC_sel, "_SLAC_MEME"),
      FEL_sel == MEME_sel ~ paste0(FEL_sel, "_FEL_MEME"),
      TRUE ~ "Not_significant"))

table(selection_res_mapping$combined_sel)
table(selection_res_mapping$SLAC_sel)
table(selection_res_mapping$FEL_sel)
table(selection_res_mapping$MEME_sel)
table(selection_res_mapping$FUBAR_sel)

length(selection_res_mapping[selection_res_mapping$combined_sel == "Purifying_SLAC_FEL_FUBAR", ]$Site)

selection_res_combine <- slac_res %>%
  left_join(
    selection_res_mapping %>%
      select(Site, SLAC_sel, FEL_sel, MEME_sel, FUBAR_sel, combined_sel),
    by = "Site") %>%
  select(-Selection)

table(selection_res_combine$combined_sel)

Diversifying           <- 2+2+1
Purifying              <- 36+18+7
Not_significant        <- 518
codon_number <- sum(Diversifying,Not_significant,Purifying)
100-((codon_number-Diversifying)/codon_number)*100 
100-((codon_number-Purifying)/codon_number)*100 
100-((codon_number-Not_significant)/codon_number)*100 

dN_dS_selection_VP2_plot <- ggplot(selection_res_combine, aes(x = Site, y = dN_dS)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_segment(aes(xend = Site, y = 1, yend = dN_dS, color = combined_sel)) +
  geom_point(aes(color = combined_sel), size = 1.8) +
  scale_color_manual(
    values = c(
    "Not_significant_SLAC_FEL_MEME_FUBAR" = "gray40",
    "Purifying_FUBAR" = "firebrick1",
    "Purifying_FEL_FUBAR" = "firebrick", 
    "Purifying_SLAC_FEL_FUBAR" = "firebrick",
    "Diversifying_MEME_FUBAR" = "darkgreen",
    "Diversifying_FUBAR" = "#4DAF4A",
    "Diversifying_MEME" = "#4DAF4A"),
    labels = c(
      "Not_significant_SLAC_FEL_MEME_FUBAR" = "Not significant",
      "Purifying_FUBAR" = "Purifying (FUBAR)",
      "Purifying_FEL_FUBAR" = "Purifying (FUBAR + FEL)",
      "Purifying_SLAC_FEL_FUBAR" = "Purifying (FUBAR + FEL + SLAC)",
      "Diversifying_MEME_FUBAR" = "Diversifying (MEME + FUBAR)",
      "Diversifying_FUBAR" = "",
      "Diversifying_MEME" = "Diversifying (MEME or FUBAR)"),
    breaks = c(
      "Diversifying_MEME_FUBAR",
      "Diversifying_MEME",
      "Purifying_SLAC_FEL_FUBAR",
      "Purifying_FEL_FUBAR",
      "Purifying_FUBAR",
      "Not_significant_SLAC_FEL_MEME_FUBAR",
      "Diversifying_FUBAR")) +
  scale_x_continuous(breaks = seq(0, max(slac_res$Site, na.rm = TRUE), by = 50)) +
  labs(x = "Codon site", y = "Substitution rate ratio (dN/dS)", color = "") +
  theme_classic()+
  theme(
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.position = "bottom")

ggsave(filename = file.path(sys_dir, "Parvo_sanger_full_seq/Rplots/paper", "dN_dS_selection_VP2_plot.png"),
       plot = dN_dS_selection_VP2_plot, width = 14, height = 8, units = "in", dpi = 600, bg = "white", limitsize = FALSE)

## write supp table to export selection scan results to Excel

overview_sel_res_df <- data.frame(
  Sheet = c("SLAC", "FEL", "MEME", "FUBAR", "Combined"),
  Description = c(
    "SLAC (Single-Likelihood Ancestor Counting): counting-based method estimating dN/dS per codon site",
    "FEL (Fixed Effects Likelihood): maximum-likelihood method estimating a single dN/dS ratio per site",
    "MEME (Mixed Effects Model of Evolution): detects episodic diversifying selection at individual sites, allowing selection to vary across branches",
    "FUBAR (Fast Unconstrained Bayesian AppRoximation): Bayesian method estimating selection per site with posterior probabilities",
    "Combined: SLAC, FEL, MEME and FUBAR results merged per codon site used to plot the paper figure"),
  stringsAsFactors = FALSE)

selection_sheets <- list(
  "Overview" = overview_sel_res_df,
  "SLAC"     = slac_res,
  "FEL"      = fel_res,
  "MEME"     = meme_res,
  "FUBAR"    = fubar_res,
  "Combined" = selection_res_combine)

wb <- createWorkbook()

for (sheet_name in names(selection_sheets)) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, selection_sheets[[sheet_name]], headerStyle = createStyle(textDecoration = "bold"))
  freezePane(wb, sheet_name, firstRow = TRUE)
  setColWidths(wb, sheet_name, cols = seq_len(ncol(selection_sheets[[sheet_name]])), widths = "auto")
}

setColWidths(wb, "Overview", cols = 1:2, widths = c(15, 100))
addStyle(wb, "Overview",
         style = createStyle(wrapText = TRUE, valign = "top"),
         rows = 2:(nrow(overview_sel_res_df) + 1), cols = 2, gridExpand = TRUE)

saveWorkbook(
  wb,
  file.path(sys_dir, "Parvo_sanger_full_seq/Evo_Analysis_Datamonkey", "supptable_selection_combined.xlsx"),
  overwrite = TRUE)

# MODIFIED MSA & Maximum Likelihood trees - DNA & Protein analysis ---------------------------------------------------------------------

#### DNA and VP2 protein analysis
#### re do the MSA and ML tree for best rep seqs

# rep_seq_to_keep <- c(
#   "OM100701.1_Canine_parvovirus_2_Egypt",
#   # "OM100698.1_Canine_parvovirus_2_Egypt",
#   # "OM100696.1_Canine_parvovirus_2_Egypt",
#   # "OM100697.1_Canine_parvovirus_2_Egypt",
#   "OM100699.1_Canine_parvovirus_2_Egypt",
#   "MG013488.1_Canine_parvovirus_2c_China",
#   "MF510157.1_Canine_parvovirus_2c_Italy",
#   "Sample57_Canine",
#   "Sample68_Canine",
#   "Sample47_Canine",
#   "Sample77_Canine",
#   "Sample76_Canine",
#   "Sample8_Canine",
#   "Sample18_Canine",
#   "MZ056886.1_Canine_parvovirus_2_Egypt",
#   "MZ056884.1_Canine_parvovirus_2_Egypt",
#   "MZ056889.1_Canine_parvovirus_2_Egypt",
#   "Sample72_Canine",
#   "sample51_canine",
#   "MZ056885.1_Canine_parvovirus_2_Egypt",
#   "OR667801.1_Canine_parvovirus_2_Iraq",
#   # "OR667800.1_Canine_parvovirus_2_Iraq",
#   # "OR667802.1_Canine_parvovirus_2_Iraq",
#   # "OR451707.1_Canine_parvovirus_2_Iraq",
#   # "OR667804.1_Canine_parvovirus_2_Iraq",
#   # "OR667803.1_Canine_parvovirus_2_Iraq",
#   "MW653250.1_Canine_parvovirus_2a_Iran",
#   "OL330980.1_Canine_parvovirus_2_Iran",
#   # "OL330978.1_Canine_parvovirus_2_Iran",
#   # "MW653249.1_Canine_parvovirus_2a_Iran",
#   "OM100702.1_Canine_parvovirus_2_Egypt",
#   # "OM100700.1_Canine_parvovirus_2_Egypt",
#   # "OL330979.1_Canine_parvovirus_2_Iran",
#   # "MZ056892.1_Canine_parvovirus_2_Egypt",
#   # "MZ056883.1_Canine_parvovirus_2_Egypt",
#   # "MZ056891.1_Canine_parvovirus_2_Egypt",
#   "OM721655.1_Canine_parvovirus_2b_Turkey",
#   "MZ056890.1_Canine_parvovirus_2_Egypt",
#   # "MZ056882.1_Canine_parvovirus_2_Egypt",
#   # "MZ056887.1_Canine_parvovirus_2_Egypt",
#   # "MZ056888.1_Canine_parvovirus_2_Egypt",
#   # "MZ056881.1_Canine_parvovirus_2_Egypt",
#   # "MW539053.1_Canine_parvovirus_2b_Turkey",
#   # "OM721656.1_Canine_parvovirus_2b_Turkey",
#   "OL330977.1_Canine_parvovirus_2_Iran",
#   "KX268117.1_Canine_parvovirus_2_Turkey",
#   # "KX268105.1_Canine_parvovirus_2_Turkey",
#   # "KX268107.1_Canine_parvovirus_2_Turkey",
#   # "KX268118.1_Canine_parvovirus_2_Turkey",
#   # "KX268106.1_Canine_parvovirus_2_Turkey",
#   # "KX268114.1_Canine_parvovirus_2_Turkey",
#   # "KX268110.1_Canine_parvovirus_2_Turkey",
#   # "KX268108.1_Canine_parvovirus_2_Turkey",
#   "MW653251.1_Canine_parvovirus_2b_Iran",
#   "MW653253.1_Canine_parvovirus_2c_Iran",
#   "KX268109.1_Canine_parvovirus_2_Turkey",
#   "FJ222821.1_ref_CPV2C",
#   "MW653256.1_Canine_parvovirus_2a_Iran",
#   "KX268115.1_Canine_parvovirus_2_Turkey",
#   "AY742933.1_recent_Canine_parvovirus_2a_NewZealand",
#   # "KX268116.1_Canine_parvovirus_2_Turkey",
#   # "KX268111.1_Canine_parvovirus_2_Turkey",
#   # "KX268113.1_Canine_parvovirus_2_Turkey",
#   # "KX268112.1_Canine_parvovirus_2_Turkey",
#   "MW653252.1_Canine_parvovirus_2b_Iran",
#   "MW653248.1_Canine_parvovirus_2a_Iran",
#   "M74849.1_ref_CPV2B",
#   "M24003.1_ref_CPV2A",
#   "M38245.1_ref_CPV2",
#   "vaccineCanine_MSD",
#   "KX268119.1_Canine_parvovirus_2_Turkey",
#   "vaccineCanine_Zoetis",
#   "PP663072.2_Feline_parvovirus_Egypt",
#   # "PP663041.2_Feline_parvovirus_Egypt",
#   # "PP663073.2_Feline_parvovirus_Egypt",
#   # "PP663070.2_Feline_parvovirus_Egypt",
#   # "PP663050.2_Feline_parvovirus_Egypt",
#   # "PV521961.1_Feline_parvovirus_Egypt",
#   # "PV521949.1_Feline_parvovirus_Egypt",
#   # "PP663045.2_Feline_parvovirus_Egypt",
#   # "PV808495.1_Feline_parvovirus_Egypt",
#   # "PV521960.1_Feline_parvovirus_Egypt",
#   # "PV521963.1_Feline_parvovirus_Egypt",
#   # "PV521942.1_Feline_parvovirus_Egypt",
#   "OM638043.1_Feline_parvovirus_Egypt",
#   # "PV521964.1_Feline_parvovirus_Egypt",
#   # "PP663074.2_Feline_parvovirus_Egypt",
#   # "PV521948.1_Feline_parvovirus_Egypt",
#   # "PV808490.1_Feline_parvovirus_Egypt",
#   # "PV808489.1_Feline_parvovirus_Egypt",
#   # "PV521962.1_Feline_parvovirus_Egypt",
#   # "PV521954.1_Feline_parvovirus_Egypt",
#   # "PV521953.1_Feline_parvovirus_Egypt",
#   # "PP663056.2_Feline_parvovirus_Egypt",
#   # "PP663055.2_Feline_parvovirus_Egypt",
#   # "PP663071.2_Feline_parvovirus_Egypt",
#   # "PP663044.2_Feline_parvovirus_Egypt",
#   # "PV521955.1_Feline_parvovirus_Egypt",
#   "PP663059.2_Feline_parvovirus_Egypt",
#   "Sample70_Feline",
#   "Sample38_Feline",
#   "Sample_Lion",
#   "Sample69_Feline",
#   "Sample27_Feline",
#   "Sample25_Feline",
#   "Sample7_Feline",
#   "Sample44_Feline",
#   # "FPV_JO24_Sample_11",
#   "PV521947.1_Feline_parvovirus_Egypt",
#   # "PP663048.2_Feline_parvovirus_Egypt",
#   # "PV808492.1_Feline_parvovirus_Egypt",
#   # "PV521952.1_Feline_parvovirus_Egypt",
#   # "PV808488.1_Feline_parvovirus_Egypt",
#   # "PV521944.1_Feline_parvovirus_Egypt",
#   # "PV808491.1_Feline_parvovirus_Egypt",
#   # "PV521951.1_Feline_parvovirus_Egypt",
#   # "PV521956.1_Feline_parvovirus_Egypt",
#   # "PP663061.2_Feline_parvovirus_Egypt",
#   # "PV521965.1_Feline_parvovirus_Egypt",
#   # "PV521945.1_Feline_parvovirus_Egypt",
#   # "PV521959.1_Feline_parvovirus_Egypt",
#   # "PV521958.1_Feline_parvovirus_Egypt",
#   # "PV521957.1_Feline_parvovirus_Egypt",
#   # "PV521943.1_Feline_parvovirus_Egypt",
#   # "PP663042.2_Feline_parvovirus_Egypt",
#   # "PV808493.1_Feline_parvovirus_Egypt",
#   # "PP663062.2_Feline_parvovirus_Egypt",
#   # "PV808496.1_Feline_parvovirus_Egypt",
#   # "PV808486.1_Feline_parvovirus_Egypt",
#   # "PP663067.2_Feline_parvovirus_Egypt",
#   # "PV521946.1_Feline_parvovirus_Egypt",
#   "PP663047.2_Feline_parvovirus_Egypt",
#   "sample12_feline",
#   "PV521950.1_Feline_parvovirus_Egypt",
#   "PP663052.2_Feline_parvovirus_Egypt",
#   "Sample67_Feline",
#   "KX434461.1_Feline_parvovirus_IZSSI_Italy",
#   "KX900570.1_Feline_parvovirus_HH6_China",
#   "M38246.1_ref_FPV",
#   "KP081409.1_Feline_parvovirus_Iran",
#   "vaccineFeline_Zoetis",
#   "vaccineFeline_MSD",
#   "EU659112.1_ref_FPV_oldest"
# )
# 
# identical(sort(rep_seq_to_keep), sort(tips_to_keep))
# 
# dna_sequences_modified <- dna_sequences_mod[names(dna_sequences_mod) %in% rep_seq_to_keep]
# protein_sequences_modified <- protein_sequences_mod[names(protein_sequences_mod) %in% rep_seq_to_keep]
# 
# length(dna_sequences_modified)
# length(protein_sequences_modified)
# 
# setdiff(rep_seq_to_keep, names(dna_sequences_modified))
# setdiff(rep_seq_to_keep, names(protein_sequences_modified))
# 
# dna_sequences_modified                     
# length(dna_sequences_modified)              
# width(dna_sequences_modified)        
# names(dna_sequences_modified)       
# dna_sequences_modified[[1]]               
# as.character(dna_sequences_modified[1])     
# letterFrequency(dna_sequences_modified, letters = "GC", as.prob = TRUE)
# 
# setequal(names(dna_sequences_modified), names(protein_sequences_modified))
# 
# #### MSA Alignment
# 
# ## DNA
# 
# dna_alignment_Muscle_modified <- msa(dna_sequences_modified, method = "Muscle")
# dna_phy_Muscle_modified <- phyDat(as.matrix(dna_alignment_Muscle_modified), type = "DNA")
# 
# ## we will stick with Muscle_modified; medium-sized datasets, highly similar sequences, and protein-coding DNA
# 
# ## trim MSA_Muscle_modified
# 
# # M38246.1_ref_FPV VP2 DNA cord = 2655 - 4585
# 
# dna_alignment_Muscle_modified_trimmed <- DNAStringSet(dna_alignment_Muscle_modified)
# ref_FPV_name <- grep("M38246.1_ref_FPV", names(dna_alignment_Muscle_modified_trimmed), value = TRUE)
# ref_FPV_seq <- dna_alignment_Muscle_modified_trimmed[[ref_FPV_name]]
# ref_FPV_chars <- strsplit(as.character(ref_FPV_seq), "")[[1]]
# non_gap_positions <- which(ref_FPV_chars != "-")
# 
# start_col <- non_gap_positions[2655]
# end_col <- non_gap_positions[4585]
# # ref_FPV_seq[2655]
# # ref_FPV_seq[4585]
# # ref_FPV_seq[start_col]
# # ref_FPV_seq[end_col]
# 
# dna_alignment_Muscle_modified_trimmed <- subseq(dna_alignment_Muscle_modified_trimmed, start = start_col, end = end_col)
# width(dna_alignment_Muscle_modified_trimmed)[1]
# dna_phy_Muscle_modified_trimmed <- phyDat(as.matrix(dna_alignment_Muscle_modified_trimmed), type = "DNA")
# 
# ## trim MSA_Muscle_modified per FJ222821.1_ref_CPV2C
# 
# # FJ222821.1_ref_CPV2C VP2 DNA cord = 1 - 1755
# 
# dna_alignment_Muscle_modified_trimmed_CPV2C <- DNAStringSet(dna_alignment_Muscle_modified)
# ref_CPV2C_name <- grep("FJ222821.1_ref_CPV2C", names(dna_alignment_Muscle_modified_trimmed_CPV2C), value = TRUE)
# ref_CPV2C_seq <- dna_alignment_Muscle_modified_trimmed_CPV2C[[ref_CPV2C_name]]
# ref_CPV2C_chars <- strsplit(as.character(ref_CPV2C_seq), "")[[1]]
# non_gap_positions_CPV2C <- which(ref_CPV2C_chars != "-")
# 
# start_col <- non_gap_positions_CPV2C[1]
# end_col <- non_gap_positions_CPV2C[1755-3]
# ref_CPV2C_seq[1]
# ref_CPV2C_seq[1755]
# ref_CPV2C_seq[start_col]
# ref_CPV2C_seq[end_col]
# 
# dna_alignment_Muscle_modified_trimmed_CPV2C <- subseq(dna_alignment_Muscle_modified_trimmed_CPV2C, start = start_col, end = end_col)
# (width(dna_alignment_Muscle_modified_trimmed_CPV2C)[1])/3
# dna_phy_Muscle_modified_trimmed_CPV2C <- phyDat(as.matrix(dna_alignment_Muscle_modified_trimmed_CPV2C), type = "DNA")
# 
# width(dna_alignment_Muscle_modified_trimmed_CPV2C)/3
# 
# # writeXStringSet(DNAStringSet(dna_alignment_Muscle_modified_trimmed_CPV2C),
# #                 file.path(sys_dir, "Parvo_sanger_full_seq/aln/Parvovirus_DNA_alignment_Muscle_modified_trimmed_CPV2C_FINAL_SVRNA_repseq_v2.fasta"))
# 
# ## Alignment
# 
# ## Protein 
# 
# protein_alignment_Muscle_modified <- msa(protein_sequences_modified, method = "Muscle")
# protein_phy_Muscle_modified <- phyDat(as.matrix(protein_alignment_Muscle_modified), type = "AA")
# 
# # writeXStringSet(DNAStringSet(dna_alignment_Muscle_modified_trimmed),
# #                 file.path(sys_dir, "Parvo_sanger_full_seq/aln/Parvovirus_DNA_alignment_Muscle_modified_trimmed_FINAL_SVRNA_repseq_v2.fasta"))
# # writeXStringSet(AAStringSet(protein_alignment_Muscle_modified),
# #                 file.path(sys_dir, "Parvo_sanger_full_seq/aln/Parvovirus_Protein_alignment_Muscle_modified_FINAL_SVRNA_repseq_v2.fasta"))
# 
# #### Maximum Likelihood trees
# 
# ## Maximum Likelihood trees from iqtree 
# 
# ## DNA
# 
# dna_ML_Muscle_modified_trimmed_tree <- read.tree(file.path(sys_dir,"Parvo_sanger_full_seq/iqtree_output/Maxll_DNA_aln_Muscle_trimmed_FINAL_SVRNA_repseq_v2.contree"))
# 
# tip_df_dna_modified <- data.frame(
#   label = dna_ML_Muscle_modified_trimmed_tree$tip.label,
#   group = ifelse(grepl("sample", dna_ML_Muscle_modified_trimmed_tree$tip.label, ignore.case = TRUE) & grepl("feline", dna_ML_Muscle_modified_trimmed_tree$tip.label, ignore.case = TRUE), "Jo_Feline",
#                  ifelse(grepl("sample", dna_ML_Muscle_modified_trimmed_tree$tip.label, ignore.case = TRUE) & grepl("canine", dna_ML_Muscle_modified_trimmed_tree$tip.label, ignore.case = TRUE), "Jo_Canine",
#                         ifelse(grepl("sample", dna_ML_Muscle_modified_trimmed_tree$tip.label, ignore.case = TRUE) & grepl("lion",   dna_ML_Muscle_modified_trimmed_tree$tip.label, ignore.case = TRUE), "Jo_Lion",
#                                ifelse(grepl("NC_001510|NC_001718|NC_029797|NC_038544", dna_ML_Muscle_modified_trimmed_tree$tip.label), "Outgroup",
#                                       ifelse(grepl("_ref_", dna_ML_Muscle_modified_trimmed_tree$tip.label), "Reference",
#                                              ifelse(grepl("vaccine", dna_ML_Muscle_modified_trimmed_tree$tip.label, ignore.case = TRUE), "Vaccine",
#                                                     ifelse(grepl("feline", dna_ML_Muscle_modified_trimmed_tree$tip.label, ignore.case = TRUE), "Regional_Feline",
#                                                            ifelse(grepl("canine", dna_ML_Muscle_modified_trimmed_tree$tip.label, ignore.case = TRUE), "Regional_Canine",
#                                                                   "Other"))))))))
# )
# 
# ## rooting v2
# 
# dna_ML_Muscle_modified_trimmed_tree_midpoint <- midpoint(dna_ML_Muscle_modified_trimmed_tree)
# dna_ML_Muscle_modified_trimmed_tree_rooted_oldestFPV <- root(dna_ML_Muscle_modified_trimmed_tree, outgroup = "EU659112.1_ref_FPV_oldest", resolve.root = TRUE)
# 
# parsimony(dna_ML_Muscle_modified_trimmed_tree, dna_phy_Muscle_modified_trimmed)
# parsimony(dna_ML_Muscle_modified_trimmed_tree_midpoint, dna_phy_Muscle_modified_trimmed)
# parsimony(dna_ML_Muscle_modified_trimmed_tree_rooted_oldestFPV, dna_phy_Muscle_modified_trimmed)
# 
# # parsimony(dna_NJ_Muscle_modified_trimmed_midpoint, dna_phy_Muscle_modified_trimmed)
# # parsimony(dna_NJ_Muscle_modified_trimmed_tree_rooted_oldestFPV, dna_phy_Muscle_modified_trimmed)
# 
# ggtreeplot_dna_ML_Muscle_modified_trimmed_tree <- ggtree(dna_ML_Muscle_modified_trimmed_tree, layout = "rectangular", size = 0.3) %<+% tip_df_dna_modified +
#   geom_tiplab(aes(color = group), size = 2, offset = 0.002, align = TRUE, linesize = 0.1) +
#   geom_tippoint(aes(color = group), size = 1.2) +
#   geom_text2(aes(label = label, subset = !isTip), size = 1.8, hjust = 1.2, vjust = -0.4) +
#   scale_color_manual(values = group_colors) +
#   theme_tree2() +
#   labs(title = paste0("Parvovirus DNA - ML", deparse(substitute(dna_ML_Muscle_modified_trimmed_tree))), color = "Group") +
#   theme(legend.position = "bottom",
#         plot.title = element_text(size = 14, face = "bold"),
#         legend.text = element_text(size = 8)) +
#   coord_cartesian(clip = "off") +
#   hexpand(0.5, direction = 1)
# 
# ggtreeplot_dna_ML_Muscle_modified_trimmed_tree_midpoint <- ggtree(dna_ML_Muscle_modified_trimmed_tree_midpoint, layout = "rectangular", size = 0.3) %<+% tip_df_dna_modified +
#   geom_tiplab(aes(color = group), size = 2, offset = 0.002, align = TRUE, linesize = 0.1) +
#   geom_tippoint(aes(color = group), size = 1.2) +
#   geom_text2(aes(label = label, subset = !isTip), size = 1.8, hjust = 1.2, vjust = -0.4) +
#   scale_color_manual(values = group_colors) +
#   theme_tree2() +
#   labs(title = paste0("Parvovirus DNA - ML", deparse(substitute(dna_ML_Muscle_modified_trimmed_tree_midpoint))), color = "Group") +
#   theme(legend.position = "bottom",
#         plot.title = element_text(size = 14, face = "bold"),
#         legend.text = element_text(size = 8)) +
#   coord_cartesian(clip = "off") +
#   hexpand(0.5, direction = 1)
# 
# ggtreeplot_dna_ML_Muscle_modified_trimmed_tree_rooted_oldestFPV <- ggtree(dna_ML_Muscle_modified_trimmed_tree_rooted_oldestFPV, layout = "rectangular", size = 0.3) %<+% tip_df_dna_modified +
#   geom_tiplab(aes(color = group), size = 2, offset = 0.002, align = TRUE, linesize = 0.1) +
#   geom_tippoint(aes(color = group), size = 1.2) +
#   geom_text2(aes(label = label, subset = !isTip), size = 1.8, hjust = 1.2, vjust = -0.4) +
#   scale_color_manual(values = group_colors) +
#   theme_tree2() +
#   labs(title = paste0("Parvovirus DNA - ML", deparse(substitute(dna_ML_Muscle_modified_trimmed_tree_rooted_oldestFPV))), color = "Group") +
#   theme(legend.position = "bottom",
#         plot.title = element_text(size = 14, face = "bold"),
#         legend.text = element_text(size = 8)) +
#   coord_cartesian(clip = "off") +
#   hexpand(0.5, direction = 1)
# 
# ggtreeplot_dna_ML_Muscle_trimmed_tree_filtered
# ggtreeplot_dna_ML_Muscle_trimmed_tree_filtered_midpoint
# ggtreeplot_dna_ML_Muscle_trimmed_tree_filtered_rooted_oldestFPV
# # ggtreeplot_dna_ML_Muscle_trimmed_tree
# # ggtreeplot_dna_ML_Muscle_trimmed_tree_midpoint
# # ggtreeplot_dna_ML_Muscle_trimmed_tree_rooted_oldestFPV 
# 
# ## Protein
# 
# protein_ML_Muscle_modified_tree <- read.tree(file.path(sys_dir,"Parvo_sanger_full_seq/iqtree_output/Maxll_Protein_aln_Muscle_FINAL_SVRNA_repseq_v2.contree"))
# 
# tip_df_protein_modified <- data.frame(
#   label = protein_ML_Muscle_modified_tree$tip.label,
#   group = ifelse(grepl("sample", protein_ML_Muscle_modified_tree$tip.label, ignore.case = TRUE) & grepl("feline", protein_ML_Muscle_modified_tree$tip.label, ignore.case = TRUE), "Jo_Feline",
#                  ifelse(grepl("sample", protein_ML_Muscle_modified_tree$tip.label, ignore.case = TRUE) & grepl("canine", protein_ML_Muscle_modified_tree$tip.label, ignore.case = TRUE), "Jo_Canine",
#                         ifelse(grepl("sample", protein_ML_Muscle_modified_tree$tip.label, ignore.case = TRUE) & grepl("lion",   protein_ML_Muscle_modified_tree$tip.label, ignore.case = TRUE), "Jo_Lion",
#                                ifelse(grepl("NC_001510|NC_001718|NC_029797|NC_038544", protein_ML_Muscle_modified_tree$tip.label), "Outgroup",
#                                       ifelse(grepl("_ref_", protein_ML_Muscle_modified_tree$tip.label), "Reference",
#                                              ifelse(grepl("vaccine", protein_ML_Muscle_modified_tree$tip.label, ignore.case = TRUE), "Vaccine",
#                                                     ifelse(grepl("feline", protein_ML_Muscle_modified_tree$tip.label, ignore.case = TRUE), "Regional_Feline",
#                                                            ifelse(grepl("canine", protein_ML_Muscle_modified_tree$tip.label, ignore.case = TRUE), "Regional_Canine",
#                                                                   "Other"))))))))
# )
# 
# ## rooting v2 
# 
# protein_ML_Muscle_modified_tree_midpoint <- midpoint(protein_ML_Muscle_modified_tree)
# protein_ML_Muscle_modified_tree_rooted_oldestFPV <- root(protein_ML_Muscle_modified_tree, outgroup = "EU659112.1_ref_FPV_oldest", resolve.root = TRUE)
# 
# parsimony(protein_ML_Muscle_modified_tree, protein_phy_Muscle_modified)
# parsimony(protein_ML_Muscle_modified_tree_midpoint,  protein_phy_Muscle_modified)
# parsimony(protein_ML_Muscle_modified_tree_rooted_oldestFPV, protein_phy_Muscle_modified)
# 
# # parsimony(protein_NJ_Muscle_modified_midpoint, protein_phy_Muscle_modified)
# # parsimony(protein_NJ_Muscle_modified_rooted_oldestFPV, protein_phy_Muscle_modified)
# 
# ggtreeplot_protein_ML_Muscle_modified_tree <- ggtree(protein_ML_Muscle_modified_tree, layout = "rectangular", size = 0.3) %<+% tip_df_protein_modified +
#   geom_tiplab(aes(color = group), size = 2, offset = 0.002, align = TRUE, linesize = 0.1) +
#   geom_tippoint(aes(color = group), size = 1.2) +
#   geom_text2(aes(label = label, subset = !isTip), size = 1.8, hjust = 1.2, vjust = -0.4) +
#   scale_color_manual(values = group_colors) +
#   theme_tree2() +
#   labs(title = paste0("Parvovirus Protein - ML", deparse(substitute(protein_ML_Muscle_modified_tree))), color = "Group") +
#   theme(legend.position = "bottom",
#         plot.title = element_text(size = 14, face = "bold"),
#         legend.text = element_text(size = 8)) +
#   coord_cartesian(clip = "off") +
#   hexpand(0.5, direction = 1)
# 
# ggtreeplot_protein_ML_Muscle_modified_tree_midpoint <- ggtree(protein_ML_Muscle_modified_tree_midpoint, layout = "rectangular", size = 0.3) %<+% tip_df_protein_modified +
#   geom_tiplab(aes(color = group), size = 2, offset = 0.002, align = TRUE, linesize = 0.1) +
#   geom_tippoint(aes(color = group), size = 1.2) +
#   geom_text2(aes(label = label, subset = !isTip), size = 1.8, hjust = 1.2, vjust = -0.4) +
#   scale_color_manual(values = group_colors) +
#   theme_tree2() +
#   labs(title = paste0("Parvovirus Protein - ML", deparse(substitute(protein_ML_Muscle_modified_tree_midpoint))), color = "Group") +
#   theme(legend.position = "bottom",
#         plot.title = element_text(size = 14, face = "bold"),
#         legend.text = element_text(size = 8)) +
#   coord_cartesian(clip = "off") +
#   hexpand(0.5, direction = 1)
# 
# ggtreeplot_protein_ML_Muscle_modified_tree_rooted_oldestFPV <- ggtree(protein_ML_Muscle_modified_tree_rooted_oldestFPV, layout = "rectangular", size = 0.3) %<+% tip_df_protein_modified +
#   geom_tiplab(aes(color = group), size = 2, offset = 0.002, align = TRUE, linesize = 0.1) +
#   geom_tippoint(aes(color = group), size = 1.2) +
#   geom_text2(aes(label = label, subset = !isTip), size = 1.8, hjust = 1.2, vjust = -0.4) +
#   scale_color_manual(values = group_colors) +
#   theme_tree2() +
#   labs(title = paste0("Parvovirus Protein - ML", deparse(substitute(protein_ML_Muscle_modified_tree_rooted_oldestFPV))), color = "Group") +
#   theme(legend.position = "bottom",
#         plot.title = element_text(size = 14, face = "bold"),
#         legend.text = element_text(size = 8)) +
#   coord_cartesian(clip = "off") +
#   hexpand(0.5, direction = 1)
# 
# ggtreeplot_protein_ML_Muscle_tree_filtered
# ggtreeplot_protein_ML_Muscle_tree_filtered_midpoint
# ggtreeplot_protein_ML_Muscle_tree_filtered_rooted_oldestFPV
# # ggtreeplot_protein_ML_Muscle_tree
# # ggtreeplot_protein_ML_Muscle_tree_midpoint
# # ggtreeplot_protein_ML_Muscle_tree_rooted_oldestFPV

## compare bootstrap values

## extract_bootstrap_info: extract tip-to-tip bootstrap support/label info from a given tree
extract_bootstrap_info <- function(tree_list) {
  get_tip_support <- function(tree) {
    node_support <- suppressWarnings(as.numeric(tree$node.label))
    Ntip <- length(tree$tip.label)
    parent_node <- sapply(seq_len(Ntip), function(tip) tree$edge[tree$edge[, 2] == tip, 1])
    setNames(node_support[parent_node - Ntip], tree$tip.label)
  }
  
  classify_support <- function(x) {
    cut(x,
        breaks = c(-Inf, 50, 85, Inf),
        labels = c("Weak", "Moderate", "Strong"),
        right = FALSE)
  }
  
  all_tips <- unique(unlist(lapply(tree_list, function(t) t$tip.label)))
  
  support_df <- data.frame(
    tip = all_tips,
    sapply(tree_list, function(tree) get_tip_support(tree)[all_tips]),
    row.names = NULL
  )
  
  tree_cols <- setdiff(names(support_df), "tip")
  label_df <- data.frame(
    tip = support_df$tip,
    lapply(support_df[tree_cols], classify_support)
  )
  names(label_df)[-1] <- paste0(tree_cols, "_label")
  
  cbind(support_df, label_df[, -1, drop = FALSE])
}

## filtered is 143 MSA and tree then filtered the tips
## modified is 60 MSA and tree directly

ML_tree_SVRNA_list <- list(
  tree_dna_FINAL_SVRNA            = dna_ML_Muscle_trimmed_tree_filtered_rooted_oldestFPV,
  tree_dna_FINAL_SVRNA_repseq     = dna_ML_Muscle_modified_trimmed_tree_rooted_oldestFPV,
  tree_protein_FINAL_SVRNA        = protein_ML_Muscle_tree_filtered_rooted_oldestFPV,
  tree_protein_FINAL_SVRNA_repseq = protein_ML_Muscle_modified_tree_rooted_oldestFPV)

# ML_tree_SVRNA_list <- list(
#   tree_dna_FINAL_SVRNA            = dna_ML_Muscle_trimmed_tree_filtered,
#   tree_dna_FINAL_SVRNA_repseq     = dna_ML_Muscle_modified_trimmed_tree,
#   tree_protein_FINAL_SVRNA        = protein_ML_Muscle_tree_filtered,
#   tree_protein_FINAL_SVRNA_repseq = protein_ML_Muscle_modified_tree)

ML_tree_SVRNA_bootstrap_comparison <- extract_bootstrap_info(ML_tree_SVRNA_list)

ML_tree_SVRNA_bootstrap_summary <- do.call(rbind, lapply(names(ML_tree_SVRNA_list), function(nm) {
  col <- paste0(nm, "_label")
  tab <- table(ML_tree_SVRNA_bootstrap_comparison[[col]], useNA = "ifany")
  data.frame(
    tree     = nm,
    Strong   = as.numeric(tab["Strong"]),
    Moderate = as.numeric(tab["Moderate"]),
    Weak     = as.numeric(tab["Weak"]),
    NA_count = sum(is.na(ML_tree_SVRNA_bootstrap_comparison[[col]])))}))

ML_tree_SVRNA_bootstrap_summary[is.na(ML_tree_SVRNA_bootstrap_summary)] <- 0

ML_tree_SVRNA_bootstrap_summary

# Pairwise identity heatmap - DNA ------------------------------------------

sample_seq_to_keep <- c(
  "Sample57_Canine",
  "Sample68_Canine",
  "Sample47_Canine",
  "Sample77_Canine",
  "Sample76_Canine",
  "Sample8_Canine",
  "Sample18_Canine",
  "Sample72_Canine",
  "sample51_canine",
  "FJ222821.1_ref_CPV2C",
  "M74849.1_ref_CPV2B",
  "M24003.1_ref_CPV2A",
  "M38245.1_ref_CPV2",
  "vaccineCanine_MSD",
  "vaccineCanine_Zoetis",
  
  "Sample70_Feline",
  "Sample38_Feline",
  "Sample_Lion",
  "Sample69_Feline",
  "Sample27_Feline",
  "Sample25_Feline",
  "Sample7_Feline",
  "Sample44_Feline",
  "sample12_feline",
  "Sample67_Feline",
  "vaccineFeline_Zoetis",
  "vaccineFeline_MSD",
  # "FPV_JO24_Sample_11",
  "M38246.1_ref_FPV",
  "EU659112.1_ref_FPV_oldest"
)

dna_sequences_Josamples <- dna_sequences_mod[names(dna_sequences_mod) %in% sample_seq_to_keep]
length(dna_sequences_Josamples)
dna_sequences_Josamples                     

#### MSA Alignment

## DNA

dna_alignment_Muscle_Josamples <- msa(dna_sequences_Josamples, method = "Muscle")
dna_phy_Muscle_Josamples <- phyDat(as.matrix(dna_alignment_Muscle_Josamples), type = "DNA")
dna_alignment_Muscle_Josamples_trimmed <- DNAStringSet(dna_alignment_Muscle_Josamples)
ref_FPV_name <- grep("M38246.1_ref_FPV", names(dna_alignment_Muscle_Josamples_trimmed), value = TRUE)
ref_FPV_seq <- dna_alignment_Muscle_Josamples_trimmed[[ref_FPV_name]]
ref_FPV_chars <- strsplit(as.character(ref_FPV_seq), "")[[1]]
non_gap_positions <- which(ref_FPV_chars != "-")
start_col <- non_gap_positions[2655]
end_col <- non_gap_positions[4585]
dna_alignment_Muscle_Josamples_trimmed <- subseq(dna_alignment_Muscle_Josamples_trimmed, start = start_col, end = end_col)
width(dna_alignment_Muscle_Josamples_trimmed)[1]
dna_phy_Muscle_Josamples_trimmed <- phyDat(as.matrix(dna_alignment_Muscle_Josamples_trimmed), type = "DNA")

dna_dist_Muscle_trimmed_Josamples <- dist.dna(as.DNAbin(dna_phy_Muscle_Josamples_trimmed), model = "raw", pairwise.deletion = TRUE)
dna_identity_Muscle_trimmed_Josamples <- 1 - as.matrix(dna_dist_Muscle_trimmed_Josamples)
dna_hclust_Muscle_trimmed_Josamples <- hclust(dna_dist_Muscle_trimmed_Josamples, method = "average")
dna_order_rank <- match(seq_along(dna_hclust_Muscle_trimmed_Josamples$order), dna_hclust_Muscle_trimmed_Josamples$order)
dna_identity_Muscle_trimmed_masked <- dna_identity_Muscle_trimmed_Josamples
dna_identity_Muscle_trimmed_masked[outer(dna_order_rank, dna_order_rank, ">")] <- NA

pheatmap(dna_identity_Muscle_trimmed_masked,
         cluster_rows = dna_hclust_Muscle_trimmed_Josamples,
         cluster_cols = dna_hclust_Muscle_trimmed_Josamples,
         color = colorRampPalette(c("darkgreen", "white", "firebrick3"))(100),
         na_col = "white",
         border_color = NA,
         fontsize = 4,
         fontsize_row = 3,
         fontsize_col = 3,
         main = "Parvovirus DNA - Pairwise Identity")

dna_dist_Muscle_trimmed_Josamples <- dist.dna(as.DNAbin(dna_phy_Muscle_modified_trimmed), model = "raw", pairwise.deletion = TRUE)
dna_identity_Muscle_trimmed_Josamples <- 1 - as.matrix(dna_dist_Muscle_trimmed_Josamples)
dna_hclust_Muscle_trimmed_Josamples <- hclust(dna_dist_Muscle_trimmed_Josamples, method = "average")
dna_order_rank <- match(seq_along(dna_hclust_Muscle_trimmed_Josamples$order), dna_hclust_Muscle_trimmed_Josamples$order)
dna_identity_Muscle_trimmed_masked <- dna_identity_Muscle_trimmed_Josamples
dna_identity_Muscle_trimmed_masked[outer(dna_order_rank, dna_order_rank, ">")] <- NA

pheatmap(dna_identity_Muscle_trimmed_masked,
         cluster_rows = dna_hclust_Muscle_trimmed_Josamples,
         cluster_cols = dna_hclust_Muscle_trimmed_Josamples,
         color = colorRampPalette(c("darkgreen", "white", "firebrick3"))(100),
         na_col = "white",
         border_color = NA,
         fontsize = 4,
         fontsize_row = 3,
         fontsize_col = 3,
         main = "Parvovirus DNA - Pairwise Identity")

p <- pheatmap(dna_identity_Muscle_trimmed_masked,
              cluster_rows = dna_hclust_Muscle_trimmed_Josamples,
              cluster_cols = dna_hclust_Muscle_trimmed_Josamples,
              color = colorRampPalette(c("darkgreen", "white", "firebrick3"))(100),
              na_col = "white",
              border_color = NA,
              fontsize = 8,
              fontsize_row = 7,
              fontsize_col = 7,
              main = "Parvovirus DNA - Pairwise Identity",
              silent = TRUE)

p$gtable$grobs[[which(p$gtable$layout$name == "row_names")]]$gp$fontface <- 2
p$gtable$grobs[[which(p$gtable$layout$name == "col_names")]]$gp$fontface <- 2

grid::grid.newpage()
grid::grid.draw(p$gtable)


