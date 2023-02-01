import os
from scipy.sparse import data
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from sklearn.metrics import roc_curve, auc, f1_score

def draw(work_dir, dataset):
    df = pd.read_csv(open(work_dir + '/performance.csv', 'r'))
    df_val = df[['AUC_val']]
    df_val.columns = ['AUC']
    df_val['Epoch'] = df_val.index + 1
    df_val['Dataset'] = 'Validation'

    df_test = df[['AUC_test']]
    df_test.columns = ['AUC']
    df_test['Dataset'] = 'Test'
    df_test['Epoch'] = df_test.index + 1

    df_infer = df[['AUC_infer']]
    df_infer.columns = ['AUC']
    df_infer['Dataset'] = 'Gencode'
    df_infer['Epoch'] = df_infer.index + 1

    df = pd.concat([df_val, df_test, df_infer])
    ax = sns.lineplot(data=df, x='Epoch', y='AUC', hue='Dataset', palette='Set1')
    ax.set_ylim((0,1))
    ax.set_xticks(np.arange(0, 21, 1))
    plt.savefig(work_dir + '/performance.csv')
    plt.close()


if __name__ == '__main__':
    draw()