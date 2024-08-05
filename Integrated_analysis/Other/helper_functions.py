
####################
# HELPER FUNCTIONS #
####################

# Function to save gene markers in an excel file where each sheet contains the markers for a specific cluster
def save_markers(adata, 
                 var_name, 
                 save_path, 
                 add_custom_score = True,  
                 order_var = "CustomScore",
                 pval_threshold = 0.01
                ):

    # import requires packages
    import pandas as pd
    import numpy as np
    
    # dictionary to hold DataFrames for each cluster
    dfs = {}  
    
    # iterate over clusters 
    for cluster in adata.obs[var_name].unique().sort_values():
        
        # get marker genes for the current cluster
        marker_genes = adata.uns["dea_"+var_name]["names"][str(cluster)].tolist()  
        scores = adata.uns["dea_"+var_name]["scores"][str(cluster)].tolist()  
        pvals = adata.uns["dea_"+var_name]["pvals"][str(cluster)].tolist() 
        pvals_adj = adata.uns["dea_"+var_name]["pvals_adj"][str(cluster)].tolist()
        logfoldchanges = adata.uns["dea_"+var_name]["logfoldchanges"][str(cluster)].tolist() 
    
        
        # create a df for the current cluster's marker genes
        df = pd.DataFrame({"Gene": marker_genes, "Score": scores, "Pval": pvals, "PvalAdj": pvals_adj, "Log2Fold": logfoldchanges})
    
        # remove non signficant genes
        df = df[df["PvalAdj"] < 0.01]
        df = df[df["Pval"] < 0.01]
        
        # include just positive log fold changes 
        df = df[df["Log2Fold"] > 0.1]
        
        # scale score 
        scores = df["Score"]
        scaled_scores = (np.array(scores) - np.min(scores)) / (np.max(scores) - np.min(scores)) 
    
        # scale log2fold
        logfoldchanges = df["Log2Fold"]
        scaled_logfoldchanges = (np.array(logfoldchanges) - np.min(logfoldchanges)) / (np.max(logfoldchanges) - np.min(logfoldchanges)) 
        
        # add scaled metrics as new df columns
        df["ScaledScore"] = scaled_scores
        df["ScaledLog2Fold"] = scaled_logfoldchanges
        
        # compute the final score
        df["CustomScore"] = (df["ScaledLog2Fold"] * df["ScaledScore"]) 
    
        # Sort df based on the final score in descending order
        df = df.sort_values(by="CustomScore", ascending=False)
        
        # add the DataFrame to the dictionary with the cluster as key
        dfs[cluster] = df
    
    # save to excel
    with pd.ExcelWriter(save_path) as writer:
        # save each cluster markers into a separate sheet
        for cluster, df in dfs.items():
            df.to_excel(writer, sheet_name=f'Cluster_{cluster}', index=False)


# function to calculate variances on *sparse* matrix
def vars(a, axis=None):
    """ Variance of sparse matrix a
    var = mean(a**2) - mean(a)**2
    """
    a_squared = a.copy()
    a_squared.data **= 2
    return a_squared.mean(axis) - np.square(a.mean(axis))




# function to compute gene signature scores and plot them
def compute_signature_score(adata, gene_set, score_name, palette, plot=True, drop=False):


    # import requires packages
    import pandas as pd
    import numpy as np   
    import scanpy as sc
    import matplotlib.pyplot as plt
    import seaborn as sns
    import matplotlib.colors as mcolors
    
    # Inital setting for plot size
    from matplotlib import rcParams
    
    FIGSIZE = (3, 3)
    rcParams["figure.figsize"] = FIGSIZE
    
    # Centered non-symmetric palette
    
    # Make mock column for plotting, here we use B cell score
    sc.tl.score_genes(adata, gene_list=gene_set, score_name=score_name, use_raw=False)
    
    
    # Palette normalization with centering and adapted dynamic range to correspond to
    # the distance of vmin and vmax from the cenetr
    # Adapted from https://stackoverflow.com/a/50003503
    class MidpointNormalize(mcolors.Normalize):
        def __init__(self, vmin=None, vmax=None, midpoint=0, clip=False):
            self.midpoint = midpoint
            mcolors.Normalize.__init__(self, vmin, vmax, clip)
    
        def __call__(self, value, clip=None):
            value = np.array(value).astype(float)
            normalized_min = max(
                0.0,
                0.5
                * (1.0 - abs((self.midpoint - self.vmin) / (self.midpoint - self.vmax))),
            )
            normalized_max = min(
                1.0,
                0.5
                * (1.0 + abs((self.vmax - self.midpoint) / (self.midpoint - self.vmin))),
            )
            normalized_mid = 0.5
            x, y = (
                [self.vmin, self.midpoint, self.vmax],
                [normalized_min, normalized_mid, normalized_max],
            )
            return np.ma.masked_array(np.interp(value, x, y))
    
    
    # Add padding arround vmin and vmax as Colorbar sets value limits to round numbers below and
    # above the vmin and vmax, respectively, which means that they can not be assigned the correct
    # color with our nomalisation function that is limited to vmin and vmax
    # However, this padding reduces the dynamic range as we set a broad padding and
    # then later discard values that are not needed for the rounding up and down
    # of the vmin and vmax on the Colorbar, respectively
    vmin = adata.obs[score_name].min()
    vmax = adata.obs[score_name].max()
    vpadding = (vmax - vmin) * 0.2
    norm = MidpointNormalize(vmin=vmin - vpadding, vmax=vmax + vpadding, midpoint=0)

    if plot == True:
        
        # Plot umap
        fig = sc.pl.umap(
            adata,
            color=score_name,
            cmap=palette,
            s=20,
            norm=norm,
            return_fig=True,
            show=False,
            use_raw=False
        )
        # Adjust Colorbar ylim to be just outside of vmin,vmax and not far outside of this range
        # as the padding we set initially may be too broad
        cmap_yticklabels = np.array([t._y for t in fig.axes[1].get_yticklabels()])
        fig.axes[1].set_ylim(
            max(cmap_yticklabels[cmap_yticklabels < vmin]),
            min(cmap_yticklabels[cmap_yticklabels > vmax]),
        )
        
    if drop == True: 
        
        adata.obs.drop(score_name, axis=1, inplace=True)


def split_umap(adata, split_by, ncol=2, nrow=None, **kwargs):
    categories = adata.obs[split_by].cat.categories
    if nrow is None:
        nrow = int(np.ceil(len(categories) / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(5*ncol, 4*nrow))
    axs = axs.flatten()
    for i, cat in enumerate(categories):
        ax = axs[i]
        sc.pl.umap(adata[adata.obs[split_by] == cat], ax=ax, show=False, title=cat, **kwargs)
    plt.tight_layout()