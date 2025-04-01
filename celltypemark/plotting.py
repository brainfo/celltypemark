import seaborn as sns
import scanpy as sc
import matplotlib.pyplot as plt

def heatmap(mark_score, save: str=None):
    sns.clustermap(mark_score, annot=False)
    if save is not None:
        plt.savefig(save, bbox_inches = 'tight')
    else:
        plt.show()