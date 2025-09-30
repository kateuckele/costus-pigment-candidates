# Costus Pigment Candidates

This repository documents my work identifying candidate genes involved in floral pigment variation in Costus.

## Purpose

In previous work, we identified a major QTL underlying the presence/absence of nectar guides. The pigment involved was identified as cyanidin, a product of the anthocyanin biosynthetic pathway. The aim of this project is to identify and evaluate candidate genes in the anthocyanin biosynthetic pathway that are located within the QTL region. The analyses here combine genome annotation, curated gene lists, BLAST searches, ortholog inference, and phylogenetic trees to generate a working set of pigment candidates.

## Repository Contents

costus-pigment-candidates/
├── data/               # Input datasets (genome annotation, candidate lists, references)
├── scripts/            # Python, bash, and R scripts used in analyses
├── results/            # Outputs (see below)
├── config/             # Configuration files (paths, tool settings)
└── README.md           # Project documentation (this file)

## Results directory
- pigment_hits/ - Filtered candidate gene hits from BLAST.
- blast_tables/ - Raw BLAST results in tabular format.
- blastdb/ - Local BLAST databases built from reference sequences.
- orthologs/ - Ortholog groupings of candidate genes.
- alignments+trees/ - Multiple sequence alignments and phylogenetic trees.
- transcripts/ - Extracted transcript sequences for candidate genes.

## Outputs
The typical end products are:
- A curated set of Costus pigment candidate genes (pigment_hits/)
- Supporting orthologs for downstream evolutionary analyses

## Notes
Analyses were run with Python 3.12 and standard genomics tools (e.g., BLAST+, gffread, AGAT).

Not all intermediate files are versioned here; some large raw data are external.

## Workflow Diagram

```mermaid
flowchart TD

classDef scripted fill:#e6f4ea,stroke:#2e7d32,color:#1b5e20
classDef manual   fill:#ffefea,stroke:#c62828,color:#8e0000
classDef output   fill:#e3f2fd,stroke:#1565c0,color:#0d47a1
classDef legend   fill:#f5f5f5,stroke:#9e9e9e,color:#424242

A["Step 1 (script): extract_transcripts.sh<br/>Translated transcripts from QTL region → <code>results/transcripts/chr6_QTL_aa_transcripts.fa</code>"]:::scripted
B["Step 2 (manual): BLASTp at NCBI (protein DB)<br/>Download XML → <code>results/blast_tables/*Alignment.xml</code>"]:::manual
C["Step 3 (script): blast_xml_to_table.py<br/>XML → TSV → <code>results/blast_tables/lasius_blast_table.tsv</code>"]:::scripted
D["Step 4 (script): summarize_blast_hits_pigments.R<br/>Filter/classify → <code>..._pigment_related.csv</code>, <code>..._pigment_candidates.tsv</code>"]:::scripted
E["Step 5 (manual): curate candidates (C1, CHI1, F3primeH, FLS)<br/>Create FASTAs → <code>data/candidates/*.fa</code>"]:::manual
F["Step 6 (script): ortholog_finder.py (via configs/orthologs.yaml)<br/>Find orthologs in other Costus species → <code>results/orthologs/*/</code>"]:::scripted
G["Step 7 (manual): MAFFT (alignments) + IQ-TREE (gene trees)<br/>→ <code>results/alignments+trees/</code>"]:::manual
H["Final outputs:<br/><code>pigment_hits/</code>, <code>orthologs/</code>, <code>alignments+trees/</code>"]:::output

A --> B --> C --> D --> E --> F --> G --> H

subgraph Legend
L1["Scripted step"]:::scripted
L2["Manual step"]:::manual
L3["Output collection"]:::output
end
```



