#### Sanger Sequencing Analysis – Parvovirus Project
#### Shadi Shahatit, RA - JUST, 2025
# Libraries ---------------------------------------------------------------

## install sangeranalyseR package
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("sangeranalyseR")

## load the package
library(sangeranalyseR)
library(Biostrings)
library(DECIPHER)

## documentation: https://sangeranalyser.readthedocs.io/en/latest/#

# sangeranalyseR - run 1 -------------------------------------------------

## define path for your seq fasta files

# sanger_ab1_parent_directory <- "C:/Users/Shadi Shahatit/OneDrive/Desktop/ParvoProject/sanger_seq_analysis/seq_run_1/"
# 
# sanger_ab1_file_names <- c(
#   "H1_12_R.81_20251007_132231.ab1",
#   "G1_12_F.81_20251007_132230.ab1",
#   "F1_12_R.79_20251007_132229.ab1",
#   "E1_12_F.79_20251007_132228.ab1",
#   "D1_12_R.555_20251007_125421.ab1",
#   "C1_12_F.555_20251007_125420.ab1")

sanger_ab1_directory <- "C:/Users/Shadi Shahatit/OneDrive/Desktop/ParvoProject/sanger_seq_analysis/seq_run_1/"

## SangerRead

sangerReadF <- SangerRead(readFeature           = "Forward Read",
                          readFileName          = paste0(sanger_ab1_directory,"sample12_81_1_F.ab1"),
                          geneticCode           = GENETIC_CODE,
                          TrimmingMethod        = "M1",                      ## can be M1 or M2 trimming methods
                          M1TrimmingCutoff      = 0.0001,
                          M2CutoffQualityScore  = NULL,
                          M2SlidingWindowSize   = NULL,
                          baseNumPerRow         = 100,
                          heightPerRow          = 200,
                          signalRatioCutoff     = 0.33,
                          showTrimmed           = TRUE)
qualityBasePlot(sangerReadF)

# writeFasta(sangerReadF,
#            outputDir         = sanger_ab1_directory,
#            compress          = FALSE,
#            compression_level = NA)
generateReport(sangerReadF,
               outputDir = sanger_ab1_directory)

## SangerAlignment

my_aligned_contigs <- SangerAlignment(ABIF_Directory      = sanger_ab1_directory,
                                      REGEX_SuffixForward = "_[0-9]*_F.ab1$",
                                      REGEX_SuffixReverse = "_[0-9]*_F.ab1$")

launchApp(my_aligned_contigs)

my_sangerAlignment <- SangerAlignment(inputSource          = "ABIF",
                                      processMethod        = "REGEX",
                                      ABIF_Directory       = sanger_ab1_directory,
                                      REGEX_SuffixForward  = "_[0-9]*_F.ab1$",
                                      REGEX_SuffixReverse  = "_[0-9]*_R.ab1$",
                                      TrimmingMethod       = "M1",
                                      M1TrimmingCutoff     = 0.0001,
                                      M2CutoffQualityScore = NULL,
                                      M2SlidingWindowSize  = NULL,
                                      baseNumPerRow        = 100,
                                      heightPerRow         = 200,
                                      signalRatioCutoff    = 0.33,
                                      showTrimmed          = TRUE,
                                      refAminoAcidSeq      = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN",
                                      minReadsNum          = 2,
                                      minReadLength        = 20,
                                      minFractionCall      = 0.5,
                                      maxFractionLost      = 0.5,
                                      geneticCode          = GENETIC_CODE,
                                      acceptStopCodons     = TRUE,
                                      readingFrame         = 1,
                                      processorsNum        = 2)

my_sangerAlignment@objectResults@readResultTable

generateReport(my_sangerAlignment,
               outputDir           = sanger_ab1_directory,
               includeSangerRead   = FALSE,
               includeSangerContig = FALSE)

# sangeranalyseR - run 2 -------------------------------------------------

## define path for your seq fasta files

sanger_ab1_directory <- "C:/Users/Shadi Shahatit/OneDrive/Desktop/sanger_seq_analysis/seq_run_2/HC01461648/"

## SangerRead

sangerReadF <- SangerRead(readFeature           = "Forward Read",
                          readFileName          = paste0(sanger_ab1_directory,"Sample18C_81_1_F.ab1"),
                          geneticCode           = GENETIC_CODE,
                          TrimmingMethod        = "M1",                      ## can be M1 or M2 trimming methods
                          M1TrimmingCutoff      = 0.0001,
                          M2CutoffQualityScore  = NULL,
                          M2SlidingWindowSize   = NULL,
                          baseNumPerRow         = 100,
                          heightPerRow          = 200,
                          signalRatioCutoff     = 0.33,
                          showTrimmed           = TRUE)
qualityBasePlot(sangerReadF)

# writeFasta(sangerReadF,
#            outputDir         = sanger_ab1_directory,
#            compress          = FALSE,
#            compression_level = NA)
generateReport(sangerReadF,
               outputDir = sanger_ab1_directory)

## SangerContig

my_sangerContig <- SangerContig(inputSource           = "ABIF",
                                processMethod         = "REGEX",
                                ABIF_Directory        = paste0(sanger_ab1_directory),
                                contigName            = "sample77_canine",
                                REGEX_SuffixForward   = "_[0-9]*_F.ab1$",
                                REGEX_SuffixReverse   = "_[0-9]*_R.ab1$",
                                # REGEX_SuffixForward  = "_F.ab1$",
                                # REGEX_SuffixReverse  = "_R.ab1$",
                                TrimmingMethod        = "M1",
                                M1TrimmingCutoff      = 0.0001,
                                M2CutoffQualityScore  = NULL,
                                M2SlidingWindowSize   = NULL,
                                baseNumPerRow         = 100,
                                heightPerRow          = 200,
                                signalRatioCutoff     = 0.33,
                                showTrimmed           = TRUE,
                                refAminoAcidSeq       = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN",
                                minReadsNum           = 2,
                                minReadLength         = 20,
                                minFractionCall       = 0.5,
                                maxFractionLost       = 0.5,
                                geneticCode           = GENETIC_CODE,
                                acceptStopCodons      = TRUE,
                                readingFrame          = 1,
                                processorsNum         = 2)

# sangeranalyseR - run 3 -------------------------------------------------

## define path for your seq fasta files

sanger_ab1_directory <- "/home/shadi/Desktop/ParvoProject/sanger_res_final/sangeranalyseR_samples_v3/"

## SangerRead

## for one tace file

sangerReadF <- SangerRead(readFeature           = "Forward Read",
                          readFileName          = paste0(sanger_ab1_directory,"sample18_canine/sample18_canine_55_F.ab1"),
                          geneticCode           = GENETIC_CODE,
                          ## can be M1 or M2 trimming methods
                          # TrimmingMethod        = "M1",                      
                          # M1TrimmingCutoff      = 0.0001,
                          # M2CutoffQualityScore  = NULL,
                          # M2SlidingWindowSize   = NULL,
                          TrimmingMethod        = "M2",
                          M1TrimmingCutoff      = NULL,
                          M2CutoffQualityScore  = 35,
                          M2SlidingWindowSize   = 4,
                          baseNumPerRow         = 100,
                          heightPerRow          = 200,
                          signalRatioCutoff     = 0.33,
                          showTrimmed           = TRUE)

qualityBasePlot(sangerReadF)
# generateReport(sangerReadF, outputDir = sanger_ab1_directory)

## SangerAlignment

## loop through the samples

sanger_ab1_directory <- "/home/shadi/Desktop/ParvoProject/sanger_res_final/sangeranalyseR_samples_v3/"

parvo_samples <- c(
  "sample18_canine",
  "sample27_feline",
  "sample38_canine",
  "sample44_feline",
  "sample47_canine",
  "sample51_canine",
  "sample56_feline",
  "sample57_canine",
  "sample67_feline",
  "sample69_feline",
  "sample70_feline",
  "sample72_canine",
  "sample76_canine",
  "sample77_canine",
  "sample_lion"
)

for (sample in parvo_samples) {
  sanger_ab1_directory_sample <- paste0(sanger_ab1_directory, sample)
  cat("Processing:", sample, "\n")
  
  my_sangerAlignment <- SangerAlignment(inputSource          = "ABIF",
                                        processMethod        = "REGEX",
                                        ABIF_Directory       = paste0(sanger_ab1_directory_sample, "/"),
                                        REGEX_SuffixForward  = "_F.ab1$",
                                        REGEX_SuffixReverse  = "_R.ab1$",
                                        ## can be M1 or M2 trimming methods
                                        # TrimmingMethod        = "M1",                      
                                        # M1TrimmingCutoff      = 0.0001,
                                        # M2CutoffQualityScore  = NULL,
                                        # M2SlidingWindowSize   = NULL,
                                        TrimmingMethod        = "M2",
                                        M1TrimmingCutoff      = NULL,
                                        M2CutoffQualityScore  = 40,             # or 29
                                        M2SlidingWindowSize   = 4,              # or 15
                                        baseNumPerRow        = 100,
                                        heightPerRow         = 200,
                                        signalRatioCutoff    = 0.20,            # could be 0.33 or 0.25
                                        showTrimmed          = TRUE,
                                        refAminoAcidSeq      = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN",
                                        minReadsNum          = 2,
                                        minReadLength        = 20,
                                        minFractionCall      = 0.5,
                                        maxFractionLost      = 0.5,
                                        geneticCode          = GENETIC_CODE,
                                        acceptStopCodons     = TRUE,
                                        readingFrame         = 1,
                                        processorsNum        = 2)
  
  fail_table <- my_sangerAlignment@objectResults@readResultTable
  if (nrow(fail_table[fail_table$creationResult == FALSE, ]) == 0) {cat("all is good, proceed :-)\n")}
  
  assign(paste0("sangerAlignment_", sample), my_sangerAlignment)
  
  contig_df <- do.call(rbind, lapply(seq_along(my_sangerAlignment@contigList), function(i) {
    data.frame(
      sample               = sample,
      contig               = names(my_sangerAlignment@contigList)[i],
      contigSeq            = as.character(my_sangerAlignment@contigList[[i]]@contigSeq),
      contigSeq_length     = nchar(as.character(my_sangerAlignment@contigList[[i]]@contigSeq)),
      
      # Fwd_PrimarySeq       = as.character(my_sangerAlignment@contigList[[i]]@forwardReadList[[1]]@primarySeq),
      # Rev_PrimarySeq       = as.character(my_sangerAlignment@contigList[[i]]@reverseReadList[[1]]@primarySeq),
      # Fwd_PriSeq_length    = nchar(as.character(my_sangerAlignment@contigList[[i]]@forwardReadList[[1]]@primarySeq)),
      # Rev_PriSeq_length    = nchar(as.character(my_sangerAlignment@contigList[[i]]@reverseReadList[[1]]@primarySeq)),
      
      Fwd_start_trim       = (my_sangerAlignment@contigList[[i]]@forwardReadList[[1]]@QualityReport@trimmedStartPos),
      Fwd_end_trim         = (my_sangerAlignment@contigList[[i]]@forwardReadList[[1]]@QualityReport@trimmedFinishPos),
      Rev_start_trim       = (my_sangerAlignment@contigList[[i]]@reverseReadList[[1]]@QualityReport@trimmedStartPos),
      Rev_end_trim         = (my_sangerAlignment@contigList[[i]]@reverseReadList[[1]]@QualityReport@trimmedFinishPos),
      
      Fwd_PrimSeq_trimmed  = as.character(subseq(my_sangerAlignment@contigList[[i]]@forwardReadList[[1]]@primarySeq,
                                                 start = (my_sangerAlignment@contigList[[i]]@forwardReadList[[1]]@QualityReport@trimmedStartPos)+1,
                                                 end = (my_sangerAlignment@contigList[[i]]@forwardReadList[[1]]@QualityReport@trimmedFinishPos))),
      Rev_PrimSeq_trimmed  = as.character(subseq(my_sangerAlignment@contigList[[i]]@reverseReadList[[1]]@primarySeq,
                                                 start = (my_sangerAlignment@contigList[[i]]@reverseReadList[[1]]@QualityReport@trimmedStartPos)+1,
                                                 end = (my_sangerAlignment@contigList[[i]]@reverseReadList[[1]]@QualityReport@trimmedFinishPos))),

      Fwd_PriSeq_trim_len  = nchar(as.character(subseq(my_sangerAlignment@contigList[[i]]@forwardReadList[[1]]@primarySeq,
                                                 start = (my_sangerAlignment@contigList[[i]]@forwardReadList[[1]]@QualityReport@trimmedStartPos)+1,
                                                 end = (my_sangerAlignment@contigList[[i]]@forwardReadList[[1]]@QualityReport@trimmedFinishPos)))),
      Rev_PriSeq_trim_len  = nchar(as.character(subseq(my_sangerAlignment@contigList[[i]]@reverseReadList[[1]]@primarySeq,
                                                 start = (my_sangerAlignment@contigList[[i]]@reverseReadList[[1]]@QualityReport@trimmedStartPos)+1,
                                                 end = (my_sangerAlignment@contigList[[i]]@reverseReadList[[1]]@QualityReport@trimmedFinishPos)))),
      
      Fwd_MeanQualityScore = my_sangerAlignment@contigList[[i]]@forwardReadList[[1]]@QualityReport@trimmedMeanQualityScore,
      Rev_MeanQualityScore = my_sangerAlignment@contigList[[i]]@reverseReadList[[1]]@QualityReport@trimmedMeanQualityScore)
  }))
  
  assign(paste0("contig_df_", sample), contig_df)
  
  fwd_rev_seqs <- c(
    setNames(DNAStringSet(contig_df$Fwd_PrimSeq_trimmed), paste0(contig_df$contig, "_fwd")),
    setNames(DNAStringSet(contig_df$Rev_PrimSeq_trimmed), paste0(contig_df$contig, "_rev")))
  
  writeXStringSet(fwd_rev_seqs, filepath = paste0(sanger_ab1_directory, "cap3_analysis_v3/", sample, ".fasta"))
  
  cat("FASTA written for:", sample, "\n")
  
  }

## explore my_sangerAlignment & contig_df

sangerAlignment_sample72_canine@objectResults@readResultTable
contig_df_sample72_canine

sangerAlignment_sample72_canine@contigList[[1]]@forwardReadList[[1]]@primaryAASeqS1
qualityBasePlot(sangerAlignment_sample72_canine@contigList$sample72_canine_81@forwardReadList[[1]])

## then do CAP3 for assembly after activate tracetrack-local
# for fasta in cap3_analysis_v3/*.fasta; do
# cap3 "$fasta"
# done
## inspect with head or cat cap3_analysis_v3/sample*.fasta.cap.contigs

writeFasta(sangerAlignment_sample18_canine,
           outputDir         = sanger_ab1_directory,
           compress          = FALSE,
           compression_level = NA,
           selection         = "all")

launchApp(sangerAlignment_sample38_canine)

generateReport(sangerAlignment_sample18_canine,
               outputDir           = sanger_ab1_directory,
               includeSangerRead   = FALSE,
               includeSangerContig = FALSE)

writeFasta(sangerAlignment_sample18_canine,
           outputDir         = sanger_ab1_directory,
           compress          = FALSE,
           compression_level = NA,
           selection         = "all")

## for one sample

my_sangerAlignment <- SangerAlignment(inputSource          = "ABIF",
                                      processMethod        = "REGEX",
                                      ABIF_Directory       = paste0(sanger_ab1_directory_sample,"/"),
                                      REGEX_SuffixForward  = "_F.ab1$",
                                      REGEX_SuffixReverse  = "_R.ab1$",
                                      TrimmingMethod       = "M2",
                                      # M1TrimmingCutoff     = 0.0001,
                                      M1TrimmingCutoff     = NULL,
                                      M2CutoffQualityScore = 29,
                                      M2SlidingWindowSize  = 15,
                                      baseNumPerRow        = 100,
                                      heightPerRow         = 200,
                                      signalRatioCutoff    = 0.33,
                                      showTrimmed          = TRUE,
                                      refAminoAcidSeq      = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN",
                                      minReadsNum          = 2,
                                      minReadLength        = 20,
                                      minFractionCall      = 0.5,
                                      maxFractionLost      = 0.5,
                                      geneticCode          = GENETIC_CODE,
                                      acceptStopCodons     = TRUE,
                                      readingFrame         = 1,
                                      processorsNum        = 2)

my_sangerAlignment@objectResults@readResultTable

if (nrow(my_sangerAlignment@objectResults@readResultTable[my_sangerAlignment@objectResults@readResultTable$creationResult == F,]) == 0) {
  print("all is good, proceed :-)")
}

contigSeq_55 <- as.character(my_sangerAlignment@contigList[[1]]@contigSeq)
contigSeq_79 <- as.character(my_sangerAlignment@contigList[[2]]@contigSeq)
contigSeq_81 <- as.character(my_sangerAlignment@contigList[[3]]@contigSeq)
nchar(contigSeq_55)
nchar(contigSeq_79)
nchar(contigSeq_81)

df <- do.call(rbind, lapply(seq_along(my_sangerAlignment@contigList), function(i)
  data.frame(
    contig = names(my_sangerAlignment@contigList)[i],
    Fwd_MeanQualityScore = my_sangerAlignment@contigList[[i]]@forwardReadList[[1]]@QualityReport@trimmedMeanQualityScore,
    Rev_MeanQualityScore = my_sangerAlignment@contigList[[i]]@reverseReadList[[1]]@QualityReport@trimmedMeanQualityScore))
)

writeFasta(my_sangerAlignment,
           outputDir         = sanger_ab1_directory,
           compress          = FALSE,
           compression_level = NA,
           selection         = "all")

generateReport(my_sangerAlignment,
               outputDir           = sanger_ab1_directory,
               includeSangerRead   = FALSE,
               includeSangerContig = FALSE)

launchApp(my_sangerAlignment)

## define path for your seq fasta files

sanger_ab1_directory <- "/home/shadi/Desktop/ParvoProject/sanger_res_final/sample18/primer_all"

my_aligned_contigs <- SangerAlignment(ABIF_Directory      = sanger_ab1_directory,
                                      REGEX_SuffixForward = "_[0-9]*_F.ab1$",
                                      REGEX_SuffixReverse = "_[0-9]*_F.ab1$")

launchApp(my_aligned_contigs)

sanger_ab1_directory <- "/home/shadi/Desktop/ParvoProject/sanger_res_final/sangeranalyseR_samples"

my_sangerAlignment <- SangerAlignment(inputSource          = "ABIF",
                                      processMethod        = "REGEX",
                                      ABIF_Directory       = sanger_ab1_directory,
                                      # REGEX_SuffixForward  = "_[0-9]*_F.ab1$",
                                      # REGEX_SuffixReverse  = "_[0-9]*_R.ab1$",
                                      REGEX_SuffixForward  = "_F.ab1$",
                                      REGEX_SuffixReverse  = "_R.ab1$",
                                      TrimmingMethod       = "M1",
                                      M1TrimmingCutoff     = 0.0001,
                                      M2CutoffQualityScore = NULL,
                                      M2SlidingWindowSize  = NULL,
                                      baseNumPerRow        = 100,
                                      heightPerRow         = 200,
                                      signalRatioCutoff    = 0.33,
                                      showTrimmed          = TRUE,
                                      refAminoAcidSeq      = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN",
                                      minReadsNum          = 2,
                                      minReadLength        = 20,
                                      minFractionCall      = 0.5,
                                      maxFractionLost      = 0.5,
                                      geneticCode          = GENETIC_CODE,
                                      acceptStopCodons     = TRUE,
                                      readingFrame         = 1,
                                      processorsNum        = 2)

my_sangerAlignment@objectResults@readResultTable
nrow(my_sangerAlignment@objectResults@readResultTable[my_sangerAlignment@objectResults@readResultTable$creationResult==F,])
head(my_sangerAlignment@objectResults@readResultTable[my_sangerAlignment@objectResults@readResultTable$creationResult==F,])
my_sangerAlignment@contigsTree$tip.label

my_sangerAlignment@contigList[[1]]@alignment
my_sangerAlignment@contigList[[1]]@contigSeq

as.character(my_sangerAlignment@contigList$sample27_feline_55@contigSeq)
my_sangerAlignment@contigList[[1]]@forwardReadList[[1]]@QualityReport
my_sangerAlignment@contigList[[1]]@reverseReadList[[1]]@QualityReport

sangerAlignment_sample18_canine@contigList[[1]]@forwardReadList[[1]]
sangerAlignment_sample18_canine@contigList[[1]]@reverseReadList[[1]]

sangerAlignment_sample18_canine@contigList[["sample18_canine_81"]]@forwardReadList[[1]]@primarySeq

seqs <- DNAStringSet(contig_df_sample18_canine$contigSeq)
names(seqs) <- contig_df_sample18_canine$contig

writeXStringSet(seqs, filepath = paste0(sanger_ab1_directory, "sample18_canine.fasta"))

fwd <- as.character(sangerAlignment_sample18_canine@contigList[["sample18_canine_81"]]@alignment$"1_Read_sample18_canine_81_F.ab1")
rev <- as.character(sangerAlignment_sample18_canine@contigList[["sample18_canine_81"]]@alignment$"2_Read_sample18_canine_81_R.ab1")
con <- as.character(sangerAlignment_sample18_canine@contigList[["sample18_canine_81"]]@alignment$Consensus)

contig_df_sample18_canine$contigSeq[3]

seqs <- DNAStringSet(contig_df_sample18_canine$contigSeq)
names(seqs) <- contig_df_sample18_canine$contig

writeXStringSet(seqs, filepath = paste0(sanger_ab1_directory, "sample18_canine.fasta"))

# sangeranalyseR - run 4 -------------------------------------------------

## define path for your seq fasta files

sanger_ab1_directory <- "/home/shadi/Desktop/ParvoProject/sanger_res_final/sangeranalyseR_samples_v4/"

## SangerRead

## TO ANALYZE one sample trace file

sangerReadF <- SangerRead(readFeature           = "Forward Read",
                          readFileName          = paste0(sanger_ab1_directory,"sample18_canine/sample18_canine_55_F.ab1"),
                          geneticCode           = GENETIC_CODE,
                          ## can be M1 or M2 trimming methods
                          TrimmingMethod        = "M1",
                          M1TrimmingCutoff      = 0.001,
                          M2CutoffQualityScore  = NULL,
                          M2SlidingWindowSize   = NULL,
                          # TrimmingMethod        = "M2",
                          # M1TrimmingCutoff      = NULL,
                          # M2CutoffQualityScore  = 20,
                          # M2SlidingWindowSize   = 10,
                          baseNumPerRow         = 100,
                          heightPerRow          = 200,
                          signalRatioCutoff     = 0.33,
                          showTrimmed           = TRUE)

qualityBasePlot(sangerReadF)

## SangerAlignment

## TO ANALYZE all samples res via a loop

sanger_ab1_directory <- "/home/shadi/Desktop/ParvoProject/sanger_res_final/sangeranalyseR_samples_v4/"

parvo_samples <- c(
  "sample18_canine",
  "sample25_feline",
  "sample27_feline",
  "sample38_canine",
  "sample44_feline",
  "sample47_canine",
  "sample51_canine",
  "sample56_feline",
  "sample57_canine",
  "sample67_feline",
  "sample68_canine",
  "sample69_feline",
  "sample70_feline",
  "sample72_canine",
  "sample76_canine",
  "sample77_canine",
  "sample8_canine",
  "sample_lion"
)

for (sample in parvo_samples) {
  sanger_ab1_directory_sample <- paste0(sanger_ab1_directory, sample)
  cat("Processing:", sample, "\n")
  
  my_sangerAlignment <- SangerAlignment(inputSource          = "ABIF",
                                        processMethod        = "REGEX",
                                        ABIF_Directory       = paste0(sanger_ab1_directory_sample, "/"),
                                        REGEX_SuffixForward  = "_F.ab1$",
                                        REGEX_SuffixReverse  = "_R.ab1$",
                                        ## can be M1 or M2 trimming methods
                                        TrimmingMethod        = "M1",
                                        M1TrimmingCutoff      = 0.001,         # 0.0001 or 0.01 or 0.001
                                        M2CutoffQualityScore  = NULL,
                                        M2SlidingWindowSize   = NULL,
                                        # TrimmingMethod        = "M2",
                                        # M1TrimmingCutoff      = NULL,
                                        # M2CutoffQualityScore  = 20,             # or 29 or 20 or 40
                                        # M2SlidingWindowSize   = 10,             # or 15 or 10 or 4
                                        baseNumPerRow        = 100,
                                        heightPerRow         = 200,
                                        signalRatioCutoff    = 0.20,            # could be 0.33 or 0.25
                                        showTrimmed          = TRUE,
                                        refAminoAcidSeq      = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN",
                                        minReadsNum          = 2,
                                        minReadLength        = 20,
                                        minFractionCall      = 0.5,
                                        maxFractionLost      = 0.5,
                                        geneticCode          = GENETIC_CODE,
                                        acceptStopCodons     = TRUE,
                                        readingFrame         = 1,
                                        processorsNum        = 2)
  
  fail_table <- my_sangerAlignment@objectResults@readResultTable
  if (nrow(fail_table[fail_table$creationResult == FALSE, ]) == 0) {cat("all is good, proceed :-)\n")}
  
  assign(paste0("sangerAlignment_", sample), my_sangerAlignment)
  
  contig_df <- do.call(rbind, lapply(seq_along(my_sangerAlignment@contigList), function(i) {
    data.frame(
      sample               = sample,
      contig               = names(my_sangerAlignment@contigList)[i],
      contigSeq            = as.character(my_sangerAlignment@contigList[[i]]@contigSeq),
      contigSeq_length     = nchar(as.character(my_sangerAlignment@contigList[[i]]@contigSeq)),
      
      # Fwd_PrimarySeq       = as.character(my_sangerAlignment@contigList[[i]]@forwardReadList[[1]]@primarySeq),
      # Rev_PrimarySeq       = as.character(my_sangerAlignment@contigList[[i]]@reverseReadList[[1]]@primarySeq),
      # Fwd_PriSeq_length    = nchar(as.character(my_sangerAlignment@contigList[[i]]@forwardReadList[[1]]@primarySeq)),
      # Rev_PriSeq_length    = nchar(as.character(my_sangerAlignment@contigList[[i]]@reverseReadList[[1]]@primarySeq)),
      
      Fwd_start_trim       = (my_sangerAlignment@contigList[[i]]@forwardReadList[[1]]@QualityReport@trimmedStartPos),
      Fwd_end_trim         = (my_sangerAlignment@contigList[[i]]@forwardReadList[[1]]@QualityReport@trimmedFinishPos),
      Rev_start_trim       = (my_sangerAlignment@contigList[[i]]@reverseReadList[[1]]@QualityReport@trimmedStartPos),
      Rev_end_trim         = (my_sangerAlignment@contigList[[i]]@reverseReadList[[1]]@QualityReport@trimmedFinishPos),
      
      Fwd_PrimSeq_trimmed  = as.character(subseq(my_sangerAlignment@contigList[[i]]@forwardReadList[[1]]@primarySeq,
                                                 start = (my_sangerAlignment@contigList[[i]]@forwardReadList[[1]]@QualityReport@trimmedStartPos)+1,
                                                 end = (my_sangerAlignment@contigList[[i]]@forwardReadList[[1]]@QualityReport@trimmedFinishPos))),
      Rev_PrimSeq_trimmed  = as.character(subseq(my_sangerAlignment@contigList[[i]]@reverseReadList[[1]]@primarySeq,
                                                 start = (my_sangerAlignment@contigList[[i]]@reverseReadList[[1]]@QualityReport@trimmedStartPos)+1,
                                                 end = (my_sangerAlignment@contigList[[i]]@reverseReadList[[1]]@QualityReport@trimmedFinishPos))),
      
      Fwd_PriSeq_trim_len  = nchar(as.character(subseq(my_sangerAlignment@contigList[[i]]@forwardReadList[[1]]@primarySeq,
                                                       start = (my_sangerAlignment@contigList[[i]]@forwardReadList[[1]]@QualityReport@trimmedStartPos)+1,
                                                       end = (my_sangerAlignment@contigList[[i]]@forwardReadList[[1]]@QualityReport@trimmedFinishPos)))),
      Rev_PriSeq_trim_len  = nchar(as.character(subseq(my_sangerAlignment@contigList[[i]]@reverseReadList[[1]]@primarySeq,
                                                       start = (my_sangerAlignment@contigList[[i]]@reverseReadList[[1]]@QualityReport@trimmedStartPos)+1,
                                                       end = (my_sangerAlignment@contigList[[i]]@reverseReadList[[1]]@QualityReport@trimmedFinishPos)))),
      
      Fwd_MeanQualityScore = my_sangerAlignment@contigList[[i]]@forwardReadList[[1]]@QualityReport@trimmedMeanQualityScore,
      Rev_MeanQualityScore = my_sangerAlignment@contigList[[i]]@reverseReadList[[1]]@QualityReport@trimmedMeanQualityScore)
  }))
  
  assign(paste0("contig_df_", sample), contig_df)
  
  fwd_rev_seqs <- c(
    setNames(DNAStringSet(contig_df$Fwd_PrimSeq_trimmed), paste0(contig_df$contig, "_fwd")),
    setNames(DNAStringSet(contig_df$Rev_PrimSeq_trimmed), paste0(contig_df$contig, "_rev")))
  
  writeXStringSet(fwd_rev_seqs, filepath = paste0(sanger_ab1_directory, "cap3_analysis_v4/", sample, "_trimmed.fasta"))
  
  cat("FASTA written for:", sample, "\n")
  
}

## explore my_sangerAlignment & contig_df

# sangerAlignment_sample72_canine@objectResults@readResultTable
# contig_df_sample72_canine
# qualityBasePlot(sangerAlignment_sample72_canine@contigList$sample72_canine_81@forwardReadList[[1]])

## passed with no errors
sangerAlignment_sample18_canine@objectResults@readResultTable
sangerAlignment_sample25_feline@objectResults@readResultTable
sangerAlignment_sample27_feline@objectResults@readResultTable
sangerAlignment_sample38_canine@objectResults@readResultTable
sangerAlignment_sample44_feline@objectResults@readResultTable
sangerAlignment_sample51_canine@objectResults@readResultTable
sangerAlignment_sample57_canine@objectResults@readResultTable
sangerAlignment_sample67_feline@objectResults@readResultTable
sangerAlignment_sample68_canine@objectResults@readResultTable
sangerAlignment_sample70_feline@objectResults@readResultTable
sangerAlignment_sample72_canine@objectResults@readResultTable
sangerAlignment_sample76_canine@objectResults@readResultTable
sangerAlignment_sample77_canine@objectResults@readResultTable
sangerAlignment_sample_lion@objectResults@readResultTable

## have errors
sangerAlignment_sample47_canine@objectResults@readResultTable
# 3 sample47_canine_79_F.ab1          FALSE MIN_READ_LENGTH_ERROR
# 4 sample47_canine_79_R.ab1          FALSE     READ_NUMBER_ERROR
sangerAlignment_sample56_feline@objectResults@readResultTable
# 1 sample56_feline_55_R.ab1          FALSE MIN_READ_LENGTH_ERROR
# 2 sample56_feline_55_F.ab1          FALSE     READ_NUMBER_ERROR
# 3 sample56_feline_79_R.ab1          FALSE MIN_READ_LENGTH_ERROR
# 4 sample56_feline_79_F.ab1          FALSE     READ_NUMBER_ERROR
sangerAlignment_sample69_feline@objectResults@readResultTable
# 3 sample69_feline_79_R.ab1          FALSE MIN_READ_LENGTH_ERROR
# 4 sample69_feline_79_F.ab1          FALSE     READ_NUMBER_ERROR
sangerAlignment_sample8_canine@objectResults@readResultTable
# 5 sample8_canine_81_R.ab1          FALSE MIN_READ_LENGTH_ERROR
# 6 sample8_canine_81_F.ab1          FALSE     READ_NUMBER_ERROR

## inspect errors visually and play around with M1 and M2

sangerReadF <- SangerRead(readFeature           = "Reverse Read",
                          readFileName          = paste0(sanger_ab1_directory,"sample56_feline/sample56_feline_79_R.ab1"),
                          geneticCode           = GENETIC_CODE,
                          ## can be M1 or M2 trimming methods
                          TrimmingMethod        = "M1",
                          M1TrimmingCutoff      = 0.01,
                          M2CutoffQualityScore  = NULL,
                          M2SlidingWindowSize   = NULL,
                          # TrimmingMethod        = "M2",
                          # M1TrimmingCutoff      = NULL,
                          # M2CutoffQualityScore  = 30,
                          # M2SlidingWindowSize   = 10,
                          baseNumPerRow         = 100,
                          heightPerRow          = 200,
                          signalRatioCutoff     = 0.33,
                          showTrimmed           = TRUE)

qualityBasePlot(sangerReadF)

## M2 is not the best option for Sanger
## all MIN_READ_LENGTH_ERROR are valid per chromatograms (expect sample56_feline_79_R.ab1 ?)
## all READ_NUMBER_ERROR are not valid as their chromatograms are good so analyze alone

## TO ANALYZE the in-house Sanger res

sanger_ab1_directory <- "/home/shadi/Desktop/ParvoProject/sanger_res_final/sangeranalyseR_samples_v4/inhouse_sanger_SQ/"

parvo_samples <- c(
  "sample7_feline",
  "sample12_feline"
)

for (sample in parvo_samples) {
  sanger_ab1_directory_sample <- paste0(sanger_ab1_directory, sample)
  cat("Processing:", sample, "\n")
  
  my_sangerAlignment <- SangerAlignment(inputSource          = "ABIF",
                                        processMethod        = "REGEX",
                                        ABIF_Directory       = paste0(sanger_ab1_directory_sample, "/"),
                                        REGEX_SuffixForward  = "_F.ab1$",
                                        REGEX_SuffixReverse  = "_R.ab1$",
                                        ## can be M1 or M2 trimming methods
                                        TrimmingMethod        = "M1",
                                        M1TrimmingCutoff      = 0.001,         # 0.0001 or 0.01 or 0.001
                                        M2CutoffQualityScore  = NULL,
                                        M2SlidingWindowSize   = NULL,
                                        # TrimmingMethod        = "M2",
                                        # M1TrimmingCutoff      = NULL,
                                        # M2CutoffQualityScore  = 20,             # or 29 or 20 or 40
                                        # M2SlidingWindowSize   = 10,             # or 15 or 10 or 4
                                        baseNumPerRow        = 100,
                                        heightPerRow         = 200,
                                        signalRatioCutoff    = 0.20,            # could be 0.33 or 0.25
                                        showTrimmed          = TRUE,
                                        refAminoAcidSeq      = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN",
                                        minReadsNum          = 2,
                                        minReadLength        = 20,
                                        minFractionCall      = 0.5,
                                        maxFractionLost      = 0.5,
                                        geneticCode          = GENETIC_CODE,
                                        acceptStopCodons     = TRUE,
                                        readingFrame         = 1,
                                        processorsNum        = 2)
  
  fail_table <- my_sangerAlignment@objectResults@readResultTable
  if (nrow(fail_table[fail_table$creationResult == FALSE, ]) == 0) {cat("all is good, proceed :-)\n")}
  
  assign(paste0("sangerAlignment_", sample), my_sangerAlignment)
  
  contig_df <- do.call(rbind, lapply(seq_along(my_sangerAlignment@contigList), function(i) {
    data.frame(
      sample               = sample,
      contig               = names(my_sangerAlignment@contigList)[i],
      contigSeq            = as.character(my_sangerAlignment@contigList[[i]]@contigSeq),
      contigSeq_length     = nchar(as.character(my_sangerAlignment@contigList[[i]]@contigSeq)),
      
      # Fwd_PrimarySeq       = as.character(my_sangerAlignment@contigList[[i]]@forwardReadList[[1]]@primarySeq),
      # Rev_PrimarySeq       = as.character(my_sangerAlignment@contigList[[i]]@reverseReadList[[1]]@primarySeq),
      # Fwd_PriSeq_length    = nchar(as.character(my_sangerAlignment@contigList[[i]]@forwardReadList[[1]]@primarySeq)),
      # Rev_PriSeq_length    = nchar(as.character(my_sangerAlignment@contigList[[i]]@reverseReadList[[1]]@primarySeq)),
      
      Fwd_start_trim       = (my_sangerAlignment@contigList[[i]]@forwardReadList[[1]]@QualityReport@trimmedStartPos),
      Fwd_end_trim         = (my_sangerAlignment@contigList[[i]]@forwardReadList[[1]]@QualityReport@trimmedFinishPos),
      Rev_start_trim       = (my_sangerAlignment@contigList[[i]]@reverseReadList[[1]]@QualityReport@trimmedStartPos),
      Rev_end_trim         = (my_sangerAlignment@contigList[[i]]@reverseReadList[[1]]@QualityReport@trimmedFinishPos),
      
      Fwd_PrimSeq_trimmed  = as.character(subseq(my_sangerAlignment@contigList[[i]]@forwardReadList[[1]]@primarySeq,
                                                 start = (my_sangerAlignment@contigList[[i]]@forwardReadList[[1]]@QualityReport@trimmedStartPos)+1,
                                                 end = (my_sangerAlignment@contigList[[i]]@forwardReadList[[1]]@QualityReport@trimmedFinishPos))),
      Rev_PrimSeq_trimmed  = as.character(subseq(my_sangerAlignment@contigList[[i]]@reverseReadList[[1]]@primarySeq,
                                                 start = (my_sangerAlignment@contigList[[i]]@reverseReadList[[1]]@QualityReport@trimmedStartPos)+1,
                                                 end = (my_sangerAlignment@contigList[[i]]@reverseReadList[[1]]@QualityReport@trimmedFinishPos))),
      
      Fwd_PriSeq_trim_len  = nchar(as.character(subseq(my_sangerAlignment@contigList[[i]]@forwardReadList[[1]]@primarySeq,
                                                       start = (my_sangerAlignment@contigList[[i]]@forwardReadList[[1]]@QualityReport@trimmedStartPos)+1,
                                                       end = (my_sangerAlignment@contigList[[i]]@forwardReadList[[1]]@QualityReport@trimmedFinishPos)))),
      Rev_PriSeq_trim_len  = nchar(as.character(subseq(my_sangerAlignment@contigList[[i]]@reverseReadList[[1]]@primarySeq,
                                                       start = (my_sangerAlignment@contigList[[i]]@reverseReadList[[1]]@QualityReport@trimmedStartPos)+1,
                                                       end = (my_sangerAlignment@contigList[[i]]@reverseReadList[[1]]@QualityReport@trimmedFinishPos)))),
      
      Fwd_MeanQualityScore = my_sangerAlignment@contigList[[i]]@forwardReadList[[1]]@QualityReport@trimmedMeanQualityScore,
      Rev_MeanQualityScore = my_sangerAlignment@contigList[[i]]@reverseReadList[[1]]@QualityReport@trimmedMeanQualityScore)
  }))
  
  assign(paste0("contig_df_", sample), contig_df)
  
  fwd_rev_seqs <- c(
    setNames(DNAStringSet(contig_df$Fwd_PrimSeq_trimmed), paste0(contig_df$contig, "_fwd")),
    setNames(DNAStringSet(contig_df$Rev_PrimSeq_trimmed), paste0(contig_df$contig, "_rev")))
  
  writeXStringSet(fwd_rev_seqs, filepath = paste0("/home/shadi/Desktop/ParvoProject/sanger_res_final/sangeranalyseR_samples_v4/", 
                                                  "cap3_analysis_v4/", sample, "_trimmed.fasta"))
  
  cat("FASTA written for:", sample, "\n")
  
}

## explore my_sangerAlignment & contig_df

## passed with no errors
sangerAlignment_sample7_feline@objectResults@readResultTable
sangerAlignment_sample12_feline@objectResults@readResultTable

## TO ANALYZE the single primer res (READ_NUMBER_ERROR + samples 12F_79R & 7F_79R)

sanger_ab1_directory <- "/home/shadi/Desktop/ParvoProject/sanger_res_final/sangeranalyseR_samples_v4/singleprimer_res/"

parvo_samples <- c(
  "sample12_feline_79_R.ab1",
  "sample7_feline_79_R.ab1",
  
  # "sample47_canine_79_R.ab1"
  # "sample56_feline_79_F.ab1",
  "sample56_feline_55_F.ab1",
  # "sample69_feline_79_F.ab1",
  "sample8_canine_81_F.ab1"
)

single_read_df <- do.call(rbind, lapply(parvo_samples, function(sample_file) {
  
  ## determine read direction from filename
  read_direction <- ifelse(grepl("_F\\.ab1$", sample_file), "Forward Read", "Reverse Read")
  sample_name    <- sub("\\.ab1$", "", sample_file)
  
  cat("Processing:", sample_file, "as", read_direction, "\n")
  
  sangerRead <- SangerRead(readFeature           = read_direction,
                           readFileName          = paste0(sanger_ab1_directory, sample_file),
                           geneticCode           = GENETIC_CODE,
                           ## can be M1 or M2 trimming methods
                           TrimmingMethod        = "M1",
                           M1TrimmingCutoff      = 0.001,         # 0.0001 or 0.01 or 0.001
                           M2CutoffQualityScore  = NULL,
                           M2SlidingWindowSize   = NULL,
                           # TrimmingMethod        = "M2",
                           # M1TrimmingCutoff      = NULL,
                           # M2CutoffQualityScore  = 20,             # or 29 or 20 or 40
                           # M2SlidingWindowSize   = 10,             # or 15 or 10 or 4
                           baseNumPerRow         = 100,
                           heightPerRow          = 200,
                           signalRatioCutoff     = 0.20,
                           showTrimmed           = TRUE)
  
  assign(paste0("sangerRead_", sample_name), sangerRead, envir = .GlobalEnv)
  
  trimmed_seq <- as.character(subseq(sangerRead@primarySeq,
                                     start = sangerRead@QualityReport@trimmedStartPos + 1,
                                     end   = sangerRead@QualityReport@trimmedFinishPos))
  
  data.frame(
    sample             = sample_name,
    read_direction     = read_direction,
    start_trim         = sangerRead@QualityReport@trimmedStartPos,
    end_trim           = sangerRead@QualityReport@trimmedFinishPos,
    PrimSeq_trimmed    = trimmed_seq,
    PriSeq_trim_len    = nchar(trimmed_seq),
    MeanQualityScore   = sangerRead@QualityReport@trimmedMeanQualityScore,
    stringsAsFactors   = FALSE)
  
}))

## write trimmed seqs to FASTA

single_fasta_seqs <- setNames(DNAStringSet(single_read_df$PrimSeq_trimmed),
                              single_read_df$sample)
writeXStringSet(single_fasta_seqs,
                filepath = paste0("/home/shadi/Desktop/ParvoProject/sanger_res_final/sangeranalyseR_samples_v4/",
                                  "cap3_analysis_v4/", "single_primer_reads_trimmed.fasta"))

## explore the res

single_read_df

sangerRead_sample12_feline_79_R
sangerRead_sample56_feline_79_F
sangerRead_sample8_canine_81_F
sangerRead_sample47_canine_79_R
sangerRead_sample69_feline_79_F
sangerRead_sample56_feline_55_F
sangerRead_sample7_feline_79_R

qualityBasePlot(sangerRead_sample12_feline_79_R)
qualityBasePlot(sangerRead_sample56_feline_79_F)
qualityBasePlot(sangerRead_sample8_canine_81_F)
qualityBasePlot(sangerRead_sample47_canine_79_R)
qualityBasePlot(sangerRead_sample69_feline_79_F)
qualityBasePlot(sangerRead_sample56_feline_55_F)
qualityBasePlot(sangerRead_sample7_feline_79_R)

## contigs

## then do CAP3 for assembly after activate tracetrack-local for cap3_analysis_v4 fasta files content
# for fasta in cap3_analysis_v3/*.fasta; do
# cap3 "$fasta"
# done
## inspect with head or cat cap3_analysis_v3/sample*.fasta.cap.contigs

## OR 

## then do CAP3 for assembly after activate tracetrack-local for cap3_analysis_v4/fasta_files_sangeranalyseR_trimmed_M1_0.0001/fasta_files_cap_contigs_M1_0.0001/ fasta files content
# for fasta in *.fasta; do
# cap3 "$fasta"
# done
## inspect with head or cat cap3_analysis_v4/fasta_files_sangeranalyseR_trimmed_M1_0.0001/fasta_files_cap_contigs_M1_0.0001/sample*.fasta.cap.contigs
## rm *.fasta *.fasta.cap.ace *.fasta.cap.contigs.links *.fasta.cap.contigs.qual *.fasta.cap.info *.fasta.cap.singlets

# sangeranalyseR - vaccine run --------------------------------------------

## TO ANALYZE the vaccine_v1 sanger in lab seq for feline and canine

sanger_ab1_directory <- "/home/shadi/Desktop/ParvoProject/sanger_res_final/sangeranalyseR_samples_v4/"

parvo_samples <- c(
  "vaccine_v1_canine",
  "vaccine_v1_feline"
)

for (sample in parvo_samples) {
  sanger_ab1_directory_sample <- paste0(sanger_ab1_directory, sample)
  cat("Processing:", sample, "\n")
  
  my_sangerAlignment <- SangerAlignment(inputSource          = "ABIF",
                                        processMethod        = "REGEX",
                                        ABIF_Directory       = paste0(sanger_ab1_directory_sample, "/"),
                                        REGEX_SuffixForward  = "_F.ab1$",
                                        REGEX_SuffixReverse  = "_R.ab1$",
                                        ## can be M1 or M2 trimming methods
                                        TrimmingMethod        = "M1",
                                        M1TrimmingCutoff      = 0.01,         # 0.0001 or 0.01 or 0.001
                                        M2CutoffQualityScore  = NULL,
                                        M2SlidingWindowSize   = NULL,
                                        # TrimmingMethod        = "M2",
                                        # M1TrimmingCutoff      = NULL,
                                        # M2CutoffQualityScore  = 20,             # or 29 or 20 or 40
                                        # M2SlidingWindowSize   = 10,             # or 15 or 10 or 4
                                        baseNumPerRow        = 100,
                                        heightPerRow         = 200,
                                        signalRatioCutoff    = 0.20,            # could be 0.33 or 0.25
                                        showTrimmed          = TRUE,
                                        refAminoAcidSeq      = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN",
                                        minReadsNum          = 2,
                                        minReadLength        = 20,
                                        minFractionCall      = 0.5,
                                        maxFractionLost      = 0.5,
                                        geneticCode          = GENETIC_CODE,
                                        acceptStopCodons     = TRUE,
                                        readingFrame         = 1,
                                        processorsNum        = 2)
  
  fail_table <- my_sangerAlignment@objectResults@readResultTable
  if (nrow(fail_table[fail_table$creationResult == FALSE, ]) == 0) {cat("all is good, proceed :-)\n")}
  
  assign(paste0("sangerAlignment_", sample), my_sangerAlignment)
  
  contig_df <- do.call(rbind, lapply(seq_along(my_sangerAlignment@contigList), function(i) {
    data.frame(
      sample               = sample,
      contig               = names(my_sangerAlignment@contigList)[i],
      contigSeq            = as.character(my_sangerAlignment@contigList[[i]]@contigSeq),
      contigSeq_length     = nchar(as.character(my_sangerAlignment@contigList[[i]]@contigSeq)),
      
      # Fwd_PrimarySeq       = as.character(my_sangerAlignment@contigList[[i]]@forwardReadList[[1]]@primarySeq),
      # Rev_PrimarySeq       = as.character(my_sangerAlignment@contigList[[i]]@reverseReadList[[1]]@primarySeq),
      # Fwd_PriSeq_length    = nchar(as.character(my_sangerAlignment@contigList[[i]]@forwardReadList[[1]]@primarySeq)),
      # Rev_PriSeq_length    = nchar(as.character(my_sangerAlignment@contigList[[i]]@reverseReadList[[1]]@primarySeq)),
      
      Fwd_start_trim       = (my_sangerAlignment@contigList[[i]]@forwardReadList[[1]]@QualityReport@trimmedStartPos),
      Fwd_end_trim         = (my_sangerAlignment@contigList[[i]]@forwardReadList[[1]]@QualityReport@trimmedFinishPos),
      Rev_start_trim       = (my_sangerAlignment@contigList[[i]]@reverseReadList[[1]]@QualityReport@trimmedStartPos),
      Rev_end_trim         = (my_sangerAlignment@contigList[[i]]@reverseReadList[[1]]@QualityReport@trimmedFinishPos),
      
      Fwd_PrimSeq_trimmed  = as.character(subseq(my_sangerAlignment@contigList[[i]]@forwardReadList[[1]]@primarySeq,
                                                 start = (my_sangerAlignment@contigList[[i]]@forwardReadList[[1]]@QualityReport@trimmedStartPos)+1,
                                                 end = (my_sangerAlignment@contigList[[i]]@forwardReadList[[1]]@QualityReport@trimmedFinishPos))),
      Rev_PrimSeq_trimmed  = as.character(subseq(my_sangerAlignment@contigList[[i]]@reverseReadList[[1]]@primarySeq,
                                                 start = (my_sangerAlignment@contigList[[i]]@reverseReadList[[1]]@QualityReport@trimmedStartPos)+1,
                                                 end = (my_sangerAlignment@contigList[[i]]@reverseReadList[[1]]@QualityReport@trimmedFinishPos))),
      
      Fwd_PriSeq_trim_len  = nchar(as.character(subseq(my_sangerAlignment@contigList[[i]]@forwardReadList[[1]]@primarySeq,
                                                       start = (my_sangerAlignment@contigList[[i]]@forwardReadList[[1]]@QualityReport@trimmedStartPos)+1,
                                                       end = (my_sangerAlignment@contigList[[i]]@forwardReadList[[1]]@QualityReport@trimmedFinishPos)))),
      Rev_PriSeq_trim_len  = nchar(as.character(subseq(my_sangerAlignment@contigList[[i]]@reverseReadList[[1]]@primarySeq,
                                                       start = (my_sangerAlignment@contigList[[i]]@reverseReadList[[1]]@QualityReport@trimmedStartPos)+1,
                                                       end = (my_sangerAlignment@contigList[[i]]@reverseReadList[[1]]@QualityReport@trimmedFinishPos)))),
      
      Fwd_MeanQualityScore = my_sangerAlignment@contigList[[i]]@forwardReadList[[1]]@QualityReport@trimmedMeanQualityScore,
      Rev_MeanQualityScore = my_sangerAlignment@contigList[[i]]@reverseReadList[[1]]@QualityReport@trimmedMeanQualityScore)
  }))
  
  assign(paste0("contig_df_", sample), contig_df)
  
  fwd_rev_seqs <- c(
    setNames(DNAStringSet(contig_df$Fwd_PrimSeq_trimmed), paste0(contig_df$contig, "_fwd")),
    setNames(DNAStringSet(contig_df$Rev_PrimSeq_trimmed), paste0(contig_df$contig, "_rev")))
  
  writeXStringSet(fwd_rev_seqs, filepath = paste0(sanger_ab1_directory, "cap3_analysis_v4/", sample, "_trimmed.fasta"))
  
  cat("FASTA written for:", sample, "\n")
  
}

## explore my_sangerAlignment & contig_df

# sangerAlignment_vaccine_v1_canine@objectResults@readResultTable
# contig_df_vaccine_v1_canine
# qualityBasePlot(sangerAlignment_vaccine_v1_canine@contigList$vaccine_v1_canine_81@forwardReadList[[1]])

## passed with no errors
sangerAlignment_vaccine_v1_canine@objectResults@readResultTable
sangerAlignment_vaccine_v1_feline@objectResults@readResultTable

## have errors
sangerAlignment_vaccine_v1_feline@objectResults@readResultTable
# 5 vaccine_v1_feline_81_R.ab1          FALSE MIN_READ_LENGTH_ERROR
# 6 vaccine_v1_feline_81_F.ab1          FALSE     READ_NUMBER_ERROR

## TO ANALYZE the single primer res for vaccine_v1 seq that failed above (READ_NUMBER_ERROR)

sanger_ab1_directory <- "/home/shadi/Desktop/ParvoProject/sanger_res_final/sangeranalyseR_samples_v4/singleprimer_res/"

parvo_samples <- c(
  "vaccine_v1_feline_81_F.ab1"
)

single_read_df <- do.call(rbind, lapply(parvo_samples, function(sample_file) {
  
  ## determine read direction from filename
  read_direction <- ifelse(grepl("_F\\.ab1$", sample_file), "Forward Read", "Reverse Read")
  sample_name    <- sub("\\.ab1$", "", sample_file)
  
  cat("Processing:", sample_file, "as", read_direction, "\n")
  
  sangerRead <- SangerRead(readFeature           = read_direction,
                           readFileName          = paste0(sanger_ab1_directory, sample_file),
                           geneticCode           = GENETIC_CODE,
                           ## can be M1 or M2 trimming methods
                           TrimmingMethod        = "M1",
                           M1TrimmingCutoff      = 0.01,         # 0.0001 or 0.01 or 0.001
                           M2CutoffQualityScore  = NULL,
                           M2SlidingWindowSize   = NULL,
                           # TrimmingMethod        = "M2",
                           # M1TrimmingCutoff      = NULL,
                           # M2CutoffQualityScore  = 20,             # or 29 or 20 or 40
                           # M2SlidingWindowSize   = 10,             # or 15 or 10 or 4
                           baseNumPerRow         = 100,
                           heightPerRow          = 200,
                           signalRatioCutoff     = 0.20,
                           showTrimmed           = TRUE)
  
  assign(paste0("sangerRead_", sample_name), sangerRead, envir = .GlobalEnv)
  
  trimmed_seq <- as.character(subseq(sangerRead@primarySeq,
                                     start = sangerRead@QualityReport@trimmedStartPos + 1,
                                     end   = sangerRead@QualityReport@trimmedFinishPos))
  
  data.frame(
    sample             = sample_name,
    read_direction     = read_direction,
    start_trim         = sangerRead@QualityReport@trimmedStartPos,
    end_trim           = sangerRead@QualityReport@trimmedFinishPos,
    PrimSeq_trimmed    = trimmed_seq,
    PriSeq_trim_len    = nchar(trimmed_seq),
    MeanQualityScore   = sangerRead@QualityReport@trimmedMeanQualityScore,
    stringsAsFactors   = FALSE)
  
}))

## write trimmed seqs to FASTA

single_fasta_seqs <- setNames(DNAStringSet(single_read_df$PrimSeq_trimmed),
                              single_read_df$sample)
writeXStringSet(single_fasta_seqs,
                filepath = paste0("/home/shadi/Desktop/ParvoProject/sanger_res_final/sangeranalyseR_samples_v4/",
                                  "cap3_analysis_v4/", "single_primer_reads_trimmed.fasta"))

## explore the res

single_read_df
sangerRead_vaccine_v1_feline_81_F
qualityBasePlot(sangerRead_vaccine_v1_feline_81_F)

## TO ANALYZE the vaccines ZOETIS and MSD sanger res from Macrogene for feline and canine

sanger_ab1_directory <- "/home/shadi/Desktop/ParvoProject/sanger_res_final/sangeranalyseR_samples_v4/"

parvo_samples <- c(
  # "vaccine_ZOETIS_canine",
  # "vaccine_ZOETIS_feline",
  "vaccine_MSD_canine"
  # "vaccine_MSD_feline"
)

for (sample in parvo_samples) {
  sanger_ab1_directory_sample <- paste0(sanger_ab1_directory, sample)
  cat("Processing:", sample, "\n")
  
  my_sangerAlignment <- SangerAlignment(inputSource          = "ABIF",
                                        processMethod        = "REGEX",
                                        ABIF_Directory       = paste0(sanger_ab1_directory_sample, "/"),
                                        REGEX_SuffixForward  = "_F.ab1$",
                                        REGEX_SuffixReverse  = "_R.ab1$",
                                        ## can be M1 or M2 trimming methods
                                        TrimmingMethod        = "M1",
                                        M1TrimmingCutoff      = 0.0001,         # 0.0001 or 0.01 or 0.001
                                        M2CutoffQualityScore  = NULL,
                                        M2SlidingWindowSize   = NULL,
                                        # TrimmingMethod        = "M2",
                                        # M1TrimmingCutoff      = NULL,
                                        # M2CutoffQualityScore  = 20,             # or 29 or 20 or 40
                                        # M2SlidingWindowSize   = 10,             # or 15 or 10 or 4
                                        baseNumPerRow        = 100,
                                        heightPerRow         = 200,
                                        signalRatioCutoff    = 0.20,            # could be 0.33 or 0.25
                                        showTrimmed          = TRUE,
                                        refAminoAcidSeq      = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN",
                                        minReadsNum          = 2,
                                        minReadLength        = 20,
                                        minFractionCall      = 0.5,
                                        maxFractionLost      = 0.5,
                                        geneticCode          = GENETIC_CODE,
                                        acceptStopCodons     = TRUE,
                                        readingFrame         = 1,
                                        processorsNum        = 2)
  
  fail_table <- my_sangerAlignment@objectResults@readResultTable
  if (nrow(fail_table[fail_table$creationResult == FALSE, ]) == 0) {cat("all is good, proceed :-)\n")}
  
  assign(paste0("sangerAlignment_", sample), my_sangerAlignment)
  
  contig_df <- do.call(rbind, lapply(seq_along(my_sangerAlignment@contigList), function(i) {
    data.frame(
      sample               = sample,
      contig               = names(my_sangerAlignment@contigList)[i],
      contigSeq            = as.character(my_sangerAlignment@contigList[[i]]@contigSeq),
      contigSeq_length     = nchar(as.character(my_sangerAlignment@contigList[[i]]@contigSeq)),
      
      # Fwd_PrimarySeq       = as.character(my_sangerAlignment@contigList[[i]]@forwardReadList[[1]]@primarySeq),
      # Rev_PrimarySeq       = as.character(my_sangerAlignment@contigList[[i]]@reverseReadList[[1]]@primarySeq),
      # Fwd_PriSeq_length    = nchar(as.character(my_sangerAlignment@contigList[[i]]@forwardReadList[[1]]@primarySeq)),
      # Rev_PriSeq_length    = nchar(as.character(my_sangerAlignment@contigList[[i]]@reverseReadList[[1]]@primarySeq)),
      
      Fwd_start_trim       = (my_sangerAlignment@contigList[[i]]@forwardReadList[[1]]@QualityReport@trimmedStartPos),
      Fwd_end_trim         = (my_sangerAlignment@contigList[[i]]@forwardReadList[[1]]@QualityReport@trimmedFinishPos),
      Rev_start_trim       = (my_sangerAlignment@contigList[[i]]@reverseReadList[[1]]@QualityReport@trimmedStartPos),
      Rev_end_trim         = (my_sangerAlignment@contigList[[i]]@reverseReadList[[1]]@QualityReport@trimmedFinishPos),
      
      Fwd_PrimSeq_trimmed  = as.character(subseq(my_sangerAlignment@contigList[[i]]@forwardReadList[[1]]@primarySeq,
                                                 start = (my_sangerAlignment@contigList[[i]]@forwardReadList[[1]]@QualityReport@trimmedStartPos)+1,
                                                 end = (my_sangerAlignment@contigList[[i]]@forwardReadList[[1]]@QualityReport@trimmedFinishPos))),
      Rev_PrimSeq_trimmed  = as.character(subseq(my_sangerAlignment@contigList[[i]]@reverseReadList[[1]]@primarySeq,
                                                 start = (my_sangerAlignment@contigList[[i]]@reverseReadList[[1]]@QualityReport@trimmedStartPos)+1,
                                                 end = (my_sangerAlignment@contigList[[i]]@reverseReadList[[1]]@QualityReport@trimmedFinishPos))),
      
      Fwd_PriSeq_trim_len  = nchar(as.character(subseq(my_sangerAlignment@contigList[[i]]@forwardReadList[[1]]@primarySeq,
                                                       start = (my_sangerAlignment@contigList[[i]]@forwardReadList[[1]]@QualityReport@trimmedStartPos)+1,
                                                       end = (my_sangerAlignment@contigList[[i]]@forwardReadList[[1]]@QualityReport@trimmedFinishPos)))),
      Rev_PriSeq_trim_len  = nchar(as.character(subseq(my_sangerAlignment@contigList[[i]]@reverseReadList[[1]]@primarySeq,
                                                       start = (my_sangerAlignment@contigList[[i]]@reverseReadList[[1]]@QualityReport@trimmedStartPos)+1,
                                                       end = (my_sangerAlignment@contigList[[i]]@reverseReadList[[1]]@QualityReport@trimmedFinishPos)))),
      
      Fwd_MeanQualityScore = my_sangerAlignment@contigList[[i]]@forwardReadList[[1]]@QualityReport@trimmedMeanQualityScore,
      Rev_MeanQualityScore = my_sangerAlignment@contigList[[i]]@reverseReadList[[1]]@QualityReport@trimmedMeanQualityScore)
  }))
  
  assign(paste0("contig_df_", sample), contig_df)
  
  fwd_rev_seqs <- c(
    setNames(DNAStringSet(contig_df$Fwd_PrimSeq_trimmed), paste0(contig_df$contig, "_fwd")),
    setNames(DNAStringSet(contig_df$Rev_PrimSeq_trimmed), paste0(contig_df$contig, "_rev")))
  
  writeXStringSet(fwd_rev_seqs, filepath = paste0(sanger_ab1_directory, "cap3_analysis_v4/", sample, "_trimmed.fasta"))
  
  cat("FASTA written for:", sample, "\n")
  
}

## explore my_sangerAlignment & contig_df

# sangerAlignment_vaccine_ZOETIS_canine@objectResults@readResultTable
# contig_df_vaccine_ZOETIS_canine
# qualityBasePlot(sangerAlignment_vaccine_ZOETIS_canine@contigList$vaccine_ZOETIS_canine_81@forwardReadList[[1]])
# sangerAlignment_vaccine_MSD_canine@objectResults@readResultTable
# contig_df_vaccine_MSD_canine
# qualityBasePlot(sangerAlignment_vaccine_MSD_canine@contigList$vaccine_MSD_canine_55@forwardReadList[[1]])

## passed with no errors
sangerAlignment_vaccine_ZOETIS_canine@objectResults@readResultTable
sangerAlignment_vaccine_ZOETIS_feline@objectResults@readResultTable

sangerAlignment_vaccine_MSD_canine@objectResults@readResultTable
sangerAlignment_vaccine_MSD_feline@objectResults@readResultTable

## have errors
## NONE for 0.0001
## NONE for 0.001
## NONE for 0.01

## TO ANALYZE the single primer res for vaccine_ZOETIS and vaccine_MSD seq that failed above (READ_NUMBER_ERROR)

sanger_ab1_directory <- "/home/shadi/Desktop/ParvoProject/sanger_res_final/sangeranalyseR_samples_v4/singleprimer_res/"

parvo_samples <- c(
  "NONE"
)


              
