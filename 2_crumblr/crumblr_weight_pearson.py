# plotting
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
import seaborn as sns

# data
import numpy as np
import pandas as pd
from scipy import stats
from scipy import sparse
import h5py
from sklearn.linear_model import LinearRegression

pal = {'FRMD4A': '#279e68', # dark green
       'PICALM': '#98df8a', # light green
       'CD163': '#ff7f0e', # dark orange
       'GPNMB': '#d62728', # red ***
       'MKI67': '#dbdb8d', # light yellow ***
       'TMEM163': '#1f77b4', # blue
       'AIF1': '#c5b0d5', # light purple
       'IFI44L': '#17becf', # cyan
       'HIF1A': '#aec7e8', # light blue
       'CCL3': '#aa40fc', # purple
       'HSPA1A': '#9edae5', # light blue
       'HIST': '#867CB9', # violet ***
       'ERN1': '#ff9896'} # light pink

# Set random seed for reproducibility
np.random.seed(42)
    
def scatter_crumblr(crumblr_x, crumblr_y, name_x='logFC_x', name_y='logFC_y'):
    
    crumblr_x.index = crumblr_x.assay
    crumblr_x = crumblr_x[['estimate','std.error']]
    crumblr_x.columns = [name_x,'SE_x']

    crumblr_y.index = crumblr_y.assay
    crumblr_y = crumblr_y[['estimate','std.error']]
    crumblr_y.columns = [name_y,'SE_y']
    
    # merge
    merged = pd.merge(crumblr_x, crumblr_y, left_index=True, right_index=True)
    merged['weight'] = 2/(merged.SE_x+merged.SE_y)
    merged['color'] = [pal[x] for x in merged.index]
    merged = merged.reset_index()
    
    # calc weighted pearson corr
    corcoef=corr(merged[name_x], merged[name_y], merged['weight']).round(2)
    
    # weighted linear regression model
    weighted_model = LinearRegression()
    weighted_model.fit(merged[name_x].values.reshape(-1, 1), merged[name_y].values.reshape(-1, 1), sample_weight=merged['weight'].values)
    X_new = np.linspace(min(merged[name_x]), max(merged[name_x]), 100).reshape(-1, 1)
    y_pred = weighted_model.predict(X_new)

    plt.rcParams["figure.figsize"] = (4, 4)
    sns.set_style("white")

    ax=sns.scatterplot(data=merged,
                       x=name_x, 
                       y=name_y,
                       size='weight',
                       sizes=(20, 200), color=merged['color'])
    # ax.set(xlabel=name_x, ylabel=name_y)
    plt.axhline(y=0, color="black", linestyle="--", lw=0.5)
    plt.axvline(x=0, color="black", linestyle="--", lw=0.5)

    for i, subtype in enumerate (merged.assay):
        plt.annotate(subtype.split('_')[-1], (merged[name_x][i]+0.002, merged[name_y][i]) ,fontsize=7)

    ax = plt.gca() # Get a matplotlib's axes instance
    plt.title("Weighted Pearson's r ={}". format(corcoef), fontsize=10)

    # The following code block adds the correlation line:
    plt.plot(X_new, y_pred, '-')

    sns.move_legend(ax, "upper left", bbox_to_anchor=(0, 1), title='Weight (1/SE)', frameon=False, alignment='left', labelspacing=0.4)

def m(x, w):
    """Weighted Mean"""
    return np.sum(x * w) / np.sum(w)

def cov(x, y, w):
    """Weighted Covariance"""
    return np.sum(w * (x - m(x, w)) * (y - m(y, w))) / np.sum(w)

def corr(x, y, w):
    """Weighted Correlation"""
    return cov(x, y, w) / np.sqrt(cov(x, x, w) * cov(y, y, w))

scatter_crumblr(crumblr_x = pd.read_csv("meta_crumblr_controls_age.csv", index_col=0),
                crumblr_y = pd.read_csv("meta_crumblr_dx.csv", index_col=0),
                name_x = 'logFC_Aging',
                name_y = 'logFC_dxAD')