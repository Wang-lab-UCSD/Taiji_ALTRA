[website layout](https://github.com/Wang-lab-UCSD/Taiji_ALTRA/blob/main/website/RA%20website.pdf)

# Homepage

- [hypothesis image](https://github.com/Wang-lab-UCSD/Taiji_ALTRA/blob/main/figures/hypothesis.png)
- text: "Proposed hypothesis to RA onset. Under the influence of risk factors such as genetics and environmental exposures, epigenetic remodeling took place in multiple cell types involving signature pathways like SUMOylation, RUNX2, YAP1, NOTCH3, and β-Catenin Pathways. The signature TFs drive a characteristic set of pro-inflammatory genes in receiver cells that can, in turn, contribute to the onset and perpetuation of RA. Diverse cell types and pathogenic mechanisms can drive a common clinical phenotype known as RA and could explain the wide variation in clinical response to agents that target individual cytokines or cell types"

# PageRank page

- [heatmap image](https://github.com/Wang-lab-UCSD/Taiji_ALTRA/blob/main/website/hp_normed2_allTFs_partial_labels.pdf)
- text: "PageRank scores heatmap of all TFs across all pseudo-bulk clusters. Labeled TFs are top 10 TFs from each Kmeans group. Clusters in columns are ordered by Kmeans group and their disease states are labeled with two color palettes. The first palette is red for Early RA, yellow for At-Risk, and green for controls. The second is green for CON and orange for At-Risk/ERA. Color of the cell indicates the normalized PageRank scores with red displaying high scores. Each Kmeans group displayed distinct dynamic patterns of TF activity. G2 has the largest number of specific TFs."
- The kmeans group information of TFs can be found [here](https://github.com/Wang-lab-UCSD/Taiji_ALTRA/blob/main/tables/Table_S5_Kmeans_groups.xlsx). For the TFs not appearing in the spreadsheet, assign them as 'non-specific'.

# Umap page

- images for each individual can be found: https://github.com/Wang-lab-UCSD/Taiji_ALTRA/tree/main/website/umap 
- Example image: [umap for participant 1](https://github.com/Wang-lab-UCSD/Taiji_ALTRA/blob/main/website/umap/p1_UMAPs.pdf).
- Each pdf includes five images:
   1. RNA-seq only
   2. ATAC-seq only
   3. Co-embedding colored by cell types
   4. Co-embedding colored by assays
   5. color palatte for cell types

# Cellular network page

- Overview: Understanding how the signature cell inflammation signals are transmitted
<!-- -  [heatmap image](https://github.com/Wang-lab-UCSD/Taiji_ALTRA/blob/main/website/hp_patients_across_celltypes_n_sig_clusters.pdf)
-  text: "Heatmap of all participants in G2 across cell types. The horizontal axis shows the individual participants and the vertical axis shows each cell type. Side bar represents the disease states of participants. Color represents the number of clusters per cell type for each participant. All the At-Risk and ERA participants had the signature in at least one cell type but the combination and distribution of cell types are highly variable. For the participants with more than one signature clusters, cellular network plots can be visualized using the search box. Try searching: p3, p10, etc" -->
<!-- -  summary images can be found: https://github.com/Wang-lab-UCSD/Taiji_ALTRA/tree/main/website/cellchat/ind_pseudobulk_cluster_network
    - we only need left hand side of each image
    - text: "Color represents the cell type and thickness of edge weight is proportional to the interaction strength. Thicker edge line indicates stronger signal. Solid and open circles represent source and target respectively. The edge thickness was normalized and comparable across different participants." -->
- specific pathway images can be found: https://github.com/Wang-lab-UCSD/Taiji_ALTRA/tree/main/website/cellchat/pathway_network
  - we only need **left** hand side of each image
  - example image title: "participant 7 - TGFB1"
  - text: "Edge thickness: interaction strength"

- outgoing/signaling figures: https://github.com/Wang-lab-UCSD/Taiji_ALTRA/tree/main/website/cellchat/sender_receiver_heatmap
  - example image title: "signaling pathway: TGFB1_TGFBR1_TGFBR2"
  - text: Red color represents the total outgoing signal strength while blue represents the total incoming signal strength
  - note: a txt file with all the pairs that user can choose. In total 69 pairs.

- expression levels of induced genes: https://github.com/Wang-lab-UCSD/Taiji_ALTRA/tree/main/website/cellchat/receptor_expression_heatmap
  - example image title: "TGFB1 induced genes expression in receptors"
  - note: a txt file with all the pairs that user can choose. In total 71 pairs. Two more pairs: `COL4A4_CD44`, `TNFSF13B_TNFRSF13C`


# Pathogenic genes page

- [heatmap image](https://github.com/Wang-lab-UCSD/Taiji_ALTRA/blob/main/website/hp_top63_feature_rna_max.pdf)
- text: "Normalized gene expression of 63 pathogenic genes in At-Risk/ERA and control participants. For each gene, the maximum gene expression across clusters was taken within each Kmeans group and each individual. 
  
  Rows represent mediators ordered by the descending predictive power. For example, MMP23B is the most predictive gene to distinguish At-Risk/ERA from controls
  
  Identified pathogenic genes are uniformly more active in At-Risk/ERAs compared to controls. Example genes include *MMP23B*, *TGFB1*, *IFNL1*, *IL15*, *NOTCH1*, and *CCL5*."

- images for each gene can be found: https://github.com/Wang-lab-UCSD/Taiji_ALTRA/tree/main/website/pathogenic 

# additional notes

1. update pathogenic page to include all 63 genes, place search bar in the center, remove the "Gene Name" text because each plot already has title
2. paper link and github link should be larger and closer
3. cellular network: 
    - shrink the image.
    - update the overview
    - center the search bar and title 

4. additional search bar: how to design the layout
    - instead of searching mediator, let user search signaling pair instead. For example: "TGFB1_TGFBR1_TGFBR2": mediator is TGFB1 and receptors are TGFBR1 and TGFBR2
    
