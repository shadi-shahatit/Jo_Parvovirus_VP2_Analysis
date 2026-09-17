# Parvo Project - Sanger Analysis

## Overview

This workflow describes the analysis of Sanger sequences for Parvovirus, including trimming sequences, comparing to reference sequences, VP2 ORF identification, multiple sequence alignment, phylogenetic analysis, and downstream sequence variation analyses.

The workflow includes:

1. Sanger sequence trimming using sangeranalyseR
2. Contig assembly and full sequence preparation
3. Addition of references and regional NCBI Viruses
4. DNA translation and VP2 ORF selection
5. DNA and protein multiple sequence alignment
6. Phylogenetic analyses with Neighbor-Joining, Maximum Likelihood, and Bayesian Inference 
7. Additional sequence variation analyses

---

# Tools

Required software:

- sangeranalyseR
- SnapGene
- SeqKit
- ExPASy Translate
- BLAST+
- IQ-TREE2
- MrBayes
- R packages for sequence alignment, phylogenetics, and visualization

---

# Workflow

## 1.1 Sanger trimming using sangeranalyseR

Our samples included 18 qPCR-positive samples of 9 canine, 9 feline, and 1 lion, together with 4 vaccine-associated strains: vaccineCanine_MSD, vaccineCanine_Zoetis, vaccineFeline_MSD, and vaccineFeline_Zoetis. Raw Sanger chromatograms were processed using `sangeranalyseR` (v4) to trim sequences and generate cleaned sequence reads. See sangeranalyseR.R

---

## 1.2 Contig assembly and full sequence preparation

Trimmed Sanger reads were assembled into contigs and full sequences using SnapGene (alternatively, with Contig). Sequence names were manually checked and adjusted to FASTA formatting.

---

## 1.3 Addition of references and regional NCBI Viruses

## Reference sequences

Reference sequences included:

| Accession | Description |
|-----------|-------------|
| EU659112.1 | FPV oldest isolate |
| M38246.1 | FPV |
| M38245.1 | CPV2 |
| M24003.1 | CPV2A |
| M74849.1 | CPV2B |
| FJ222821.1 | CPV2c |

## NCBI Virus sequences

Sequences were downloaded from NCBI Virus: https://www.ncbi.nlm.nih.gov/labs/virus/vssi/ (accessed on 20/07/2026). 

Filtering criteria:

- Taxonomy:
  - Feline parvovirus (taxid:10785)
  - Canine parvovirus 2 (taxid:246878)

- Minimum sequence length:
  - 1000 bp

- Hosts:
  - Felis catus, domestic cat (taxid:9685)
  - Felidae cat family (taxid:9681)
  - Canis lupus familiaris dog (taxid:9615)
  - Canidae dog, coyote, wolf, fox (taxid:9608)

- Geographic regions:
  - Egypt
  - Iran
  - Turkey
  - Iraq

A total of 107 NCBI Virus sequences were included. We also added key global sequences for comparison, including KX434461.1 (Italy), KX900570.1 (China), OM638043 (Egypt), AY742933 (recent CPV2a), MG013488 (CPV2c, China), and MF510157 (CPV2c, Italy).

## Outgroup sequences

We tested several outgroup sequences, including NC_001510 (Protoparvovirus rodent1), NC_001718.1 (Protoparvovirus ungulate1), NC_029797 (Megabat bufavirus), and NC_038544 (Primate protoparvovirus1), but divergence was too high, so we eventually rooted the tree with the oldest FPV sequence instead (EU659112.1).

## Combine sequences

All sample, vaccine, reference, and NCBI Virus sequences were combined:

```bash
seqkit seq *.fa refs/*.fa NCBIvirus_seq_filtered_final.fasta add_seq/*.fasta > combined_sample_v3_vac_ref_ncbi_add_v142.fasta
seqkit stats combined_sample_v3_vac_ref_ncbi_add_v142.fasta 
```

## 1.4 DNA translation and VP2 ORF selection

## ExPASy translation

DNA sequences were translated using **ExPASy six-frame translation**.

```bash
## run Expasy on bash for all samples from sample_vac_ref_ncbi_out_combined

f="combined_sample_vac_ref_ncbi_out.fasta" 

awk -v RS=">" 'NR>1 {
    n=split($0, lines, "\n")
    header=lines[1]
    seq=""
    for(i=2;i<=n;i++) seq=seq lines[i]
    print header "\t" seq
}' "$f" | while IFS=$'\t' read -r header seq; do
    id=$(echo "$header" | cut -d'/' -f1 | cut -d' ' -f1)
    curl -s -d "dna_sequence=${seq}&output_format=fasta" \
        https://web.expasy.org/cgi-bin/translate/dna2aa.cgi \
        > "Expasy_translate/Expasy_output/${id}_translated.fasta"
    sleep 1
done

## pick best ORF from Expasy output files

## use pickbestORF_v2.py
## this scirpt parses ExPASy 6-frame translation FASTA, collects candidate ORFs (M-stop/M-end) >=200aa across all 6 frames, BLASTs against VP2 protein ref, keeps hits >=75% identity,
## picks the longest passing candidate per sample, writes best-ORF record + full log. Run as python3 pickbestORF_v2.py input.fasta output.fasta vp2_prot_ref.fasta [min_len] [min_pident]

# build the BLAST db once

makeblastdb -in refs/M38246.1_ref_FPV_VP2_aa_Expasy.txt -dbtype prot -out refs/blastdb/vp2_db
ls -la refs/blastdb/

vp2_db="refs/blastdb/vp2_db"
logfile="Expasy_translate/Expasy_output/best_orf/best_orf_log_final_21072026.txt"
> "$logfile"

for out in Expasy_translate/Expasy_output/*_translated.fasta; do
    id=$(basename "$out" _translated.fasta)
    python3 Expasy_translate/pickbestORF_v2.py "$out" "Expasy_translate/Expasy_output/best_orf/${id}_best_orf.fasta" "$vp2_db" 150 30 2>&1 | tee -a "$logfile"
done

# combine all best-ORF fasta

seqkit seq Expasy_translate/Expasy_output/best_orf/*_best_orf.fasta > Expasy_translate/Expasy_output/best_orf/combined_best_orf.fasta
seqkit stats Expasy_translate/Expasy_output/best_orf/combined_best_orf.fasta
```

## 5. DNA and protein multiple sequence alignment

DNA and VP2 protein sequences were aligned using **MUSCLE** and **ClustalOmega** via R `msa` package.
The final downstream analyses were performed using MUSCLE alignments.
DNA alignments were generated for the full nucleotide dataset (dna) and the VP2 coding region (protein), while DNA alignments were then trimmed for the VP2 region per M38246.1_ref_FPV: 2655-4585 bp.

MSA and downstream analyses were performed using `ParvoProject_MSA_Tree_v2.R`.

---

## 6. Neighbor-Joining and Maximum Likelihood phylogenetic analysis

Phylogenetic trees were reconstructed from DNA and protein alignments using **Neighbor-Joining (NJ)** and **Maximum Likelihood (ML)** approaches.

Neighbor-Joining trees were generated in R using the `ape` and `phangorn` packages, with nucleotide and amino acid substitution models evaluated before selecting the best-performing models.

Maximum Likelihood trees were generated using `IQ-TREE2` via

```bash
iqtree2
	-s aln/MSA_file.fasta \
	-m MFP \
	-B 1000 \
	-alrt 1000 \
	-bnni \
	-T 6 \
	-seed 1001 \
	-vv \
	-pre iqtree_output/Maxll_DNA_aln_Muscle_trimmed \
	-redo
```

---

## 7. Additional sequence variation analyses

Sequence variation analyses included:

- Principal Coordinate Analysis (PCoA) based on DNA and VP2 protein distance matrices.
- Nucleotide diversity (π) estimation across defined sequence groups.
- VP2 amino acid variant identification from protein alignments relative to all reference sequences.
- Generation of amino acid substitution tables comparing samples and reference sequences.

## Output

The analysis generates the following final outputs:

| Output | Description |
|--------|-------------|
| Combined DNA FASTA file | Final sequences for our samples, vaccines, references, NCBI viruses, and outgroup sequences |
| Combined protein FASTA file | Final Expasy translated amino acid sequences selected from six-frame translations after ORF filtertion |
| DNA multiple sequence alignment | MUSCLE and ClustalOmega-generated alignments of DNA and protein sequences |
| Neighbor-Joining Phylogenetic trees | Neighbor-Joining and Maximum Likelihood trees generated from DNA and protein (rooted and unrooted) |
| Amino acid variants from FPV | Amino acid substitutions identified relative to the FPV reference sequence |
| VP2 amino acid variants table | Three-letter formated amino acid substitution table to compare between samples, vaccine, and references |


