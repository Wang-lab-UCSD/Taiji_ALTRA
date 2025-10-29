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
<!-- -  [heatmap image](https://github.com/Wang-lab-UCSD/Taiji_ALTRA/blob/main/website/hp_patients_across_celltypes_n_sig_clusters.pdf)
-  text: "Heatmap of all participants in G2 across cell types. The horizontal axis shows the individual participants and the vertical axis shows each cell type. Side bar represents the disease states of participants. Color represents the number of clusters per cell type for each participant. All the At-Risk and ERA participants had the signature in at least one cell type but the combination and distribution of cell types are highly variable. For the participants with more than one signature clusters, cellular network plots can be visualized using the search box. Try searching: p3, p10, etc" -->
<!-- -  summary images can be found: https://github.com/Wang-lab-UCSD/Taiji_ALTRA/tree/main/website/cellchat/ind_pseudobulk_cluster_network
    - we only need left hand side of each image
    - text: "Color represents the cell type and thickness of edge weight is proportional to the interaction strength. Thicker edge line indicates stronger signal. Solid and open circles represent source and target respectively. The edge thickness was normalized and comparable across different participants." -->
- specific pathway images can be found: https://github.com/Wang-lab-UCSD/Taiji_ALTRA/tree/main/website/cellchat/pathway_network
  - we only need **left** hand side of each image
  - example image title: "participant 7 - TGFB1"
  - text: "Edge thickness: interaction strength"

# Pathogenic genes page

- [heatmap image](https://github.com/Wang-lab-UCSD/Taiji_ALTRA/blob/main/website/hp_top30_feature_rna_max.pdf)
- text: "Normalized gene expression of top 30 predictors for At-Risk/ERA and control participants in G2 and G4 clusters respectively. For each gene, the maximum gene expression across clusters was taken within each Kmeans group and each individual. Rows represent mediators while columns represent patients. Red cell represents a higher expression level of the cytokine in the patient. Top 30 predictor cytokines are uniformly more active in At-Risk/ERAs compared to controls. Example genes include *MMP23B*, *TGFB1*, *IFNL1*, *IL15*, *NOTCH1*, and *CCL5*."
- images for each gene can be found: https://github.com/Wang-lab-UCSD/Taiji_ALTRA/tree/main/website/pathogenic 

# additional notes

1. text don't take the total width, should leave some space from both sides
2. paper link and github link should be larger and closer
3. Example Umap should take the whole screen.
4. Provide toggle to select participant in Umap and toggle to select gene in pathogenic section
5. affected samples: KT00006 (p34), KT00012 (p35), KT00015 (p36), KT00055 (p10), KT00056 (p11), KT00057 (p12), KT00060 (p13)

