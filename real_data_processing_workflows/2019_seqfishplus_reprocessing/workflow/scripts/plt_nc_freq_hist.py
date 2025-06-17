# -*- coding: utf-8 -*-
"""
Created on Wed Feb 23 19:53:07 2022

@author: jonat
"""

import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import nbinom, poisson
import scipy.stats as stats
import numpy as np

bc_files = snakemake.input

barcodes = pd.read_csv(bc_files[0])

def neg_cntrl_saturated_ccs(df):
    n_nctrls = np.sum(df.gene == "negative_control")
    remove_if_true = len(df) > 5 and n_nctrls/len(df) > 0.3
    #print("remove_if_true: ",  remove_if_true)
    return not remove_if_true

barcodes_filtered = barcodes.groupby("cc").filter(neg_cntrl_saturated_ccs)

barcodes.set_index("gene", inplace=True)
barcodes_filtered.set_index("gene", inplace=True)

ch = bc_files[0].split('_')[2]
try:
    nc_bcs = pd.DataFrame(barcodes.loc['negative_control'])
except:
    nc_bcs = pd.DataFrame(columns = barcodes.columns)

try:
    nc_bcs_filtered = pd.DataFrame(barcodes_filtered.loc['negative_control'])
except:
    nc_bcs_filtered = pd.DataFrame(columns = barcodes_filtered.columns)


nc_bcs.loc[:,"ch"] = ch

for bcs_fname in bc_files[1:]:
    #print(len(nc_counts))

    barcodes = pd.read_csv(bcs_fname)
    barcodes_filtered = barcodes.groupby("cc").filter(neg_cntrl_saturated_ccs)


    barcodes.set_index("gene", inplace=True)

    try:
        _nc_bcs = pd.DataFrame(barcodes.loc['negative_control'])
        ch = bcs_fname.split('_')[2]

        _nc_bcs.loc[:,"ch"] = ch
        nc_bcs = nc_bcs.append(_nc_bcs)
    except:
        pass

    try:
        _nc_bcs = pd.DataFrame(barcodes_filtered.loc['negative_control'])
        ch = bcs_fname.split('_')[2]

        _nc_bcs.loc[:,"ch"] = ch
        nc_bcs_filtered = nc_bcs_filtered.append(_nc_bcs)
    except:
        pass

#print("nc_bcs")
#print(nc_bcs)

ncells = len(list(nc_bcs.groupby(['pos','cellid'])))

def get_mean_sd_poisson(srs):
    vec = np.array(srs)
    nzeros = ncells - len(vec)
    vec = np.append(vec, [0]*nzeros)
    #df_w_zeros = df.append([0]*nzeros)

    return(np.mean(vec), np.var(vec))

#nc_bcs_grpd = nc_bcs.groupby(["gene_number","ch"])

nc_counts = nc_bcs.groupby(["gene_number","ch"]).size()

cell_nc_counts = nc_bcs.groupby(["gene_number","pos","cellid","ch"]).size()

var_mus = cell_nc_counts.groupby(["gene_number","ch"]).apply(get_mean_sd_poisson)

#print(var_mus)

#print(nc_counts)
#print("mu: ")
#print(nc_counts.mean())
#print("var: ")
#print(nc_counts.var())

max_barcodes = int(nc_counts.max())

barcode_cnts = np.array(range(max_barcodes+1))

mu = nc_counts.mean()
sigsq = nc_counts.var()
nhat_mm = mu**2/(sigsq - mu)
phat_mm = mu/sigsq

pm = nbinom.pmf(barcode_cnts, nhat_mm, phat_mm)

mu = nc_counts.mean()


#print("max(nc_counts): ", max(nc_counts))
#print("min(nc_counts): ", min(nc_counts))

counts, bins, bc_ob = plt.hist(nc_counts, bins=max(nc_counts)-min(nc_counts), label="observed")

u_false_counts, ncws = np.unique(nc_counts, return_counts=True)

false_counts = list(range(int(max(u_false_counts))+1))
n_neg_ctrl_cws = np.zeros(len(false_counts))
n_neg_ctrl_cws[u_false_counts] = ncws

#expected = poisson.pmf(false_counts, mu)*n_neg_ctrl_cws.sum()
#expected = poisson.pmf(false_counts, mu)*len(nc_counts)
nbinom_expected = nbinom.pmf(barcode_cnts, nhat_mm, phat_mm)*len(nc_counts)
poisson_expected = poisson.pmf(barcode_cnts, mu)*len(nc_counts)

def poisson_dispersion_test(cnts): #(df):
    #cnts = df["0"]
    cnts = np.array(cnts)
    cnts = np.append(cnts, [0]*(ncells-len(cnts)))
    chi2_score = ncells*np.var(cnts)/np.mean(cnts)
    p = 1 - stats.chi2.cdf(chi2_score, ncells-1)
    return p

pdt_p = poisson_dispersion_test(nc_counts)

#plt.plot(false_counts, n_neg_ctrl_cws, label="observed")
plt.plot(false_counts, nbinom_expected, '.', label="Negative Binomial")
plt.plot(false_counts, poisson_expected, '*', label="Poisson")

plt.title("phat = %.2f, nhat = %.2f, mu = %.2f, pdt_p = %.2f" % (phat_mm, nhat_mm, mu, pdt_p))
plt.legend()
plt.xlabel("Number of Barcodes Found")
plt.ylabel("Number of Negative Control Codewords")
plt.savefig(snakemake.output[0])

nc_counts.to_csv(snakemake.output[1], index=False)
cell_nc_counts.to_csv(snakemake.output[2], index=True)


#### Filtered
#print(nc_bcs_filtered)
nc_bcs_filtered = pd.DataFrame(nc_bcs_filtered)
#print(type(nc_bcs_filtered))

#print(nc_bcs_filtered.columns)

if len(nc_bcs_filtered) == 0:
    nc_counts = np.array([])
    cell_nc_counts = pd.DataFrame(columns=cell_nc_counts.columns)
else:
    nc_counts = nc_bcs_filtered.groupby(["gene_number","ch"]).size()
    nc_counts = np.array(nc_counts)
    cell_nc_counts = nc_bcs_filtered.groupby(["gene_number","pos","cellid","ch"]).size()

nc_counts = np.append(nc_counts, [0]*(14000-len(nc_counts)))


var_mus = cell_nc_counts.groupby(["gene_number","ch"]).apply(get_mean_sd_poisson)

max_barcodes = int(nc_counts.max())

barcode_cnts = np.array(range(max_barcodes+1))

mu = nc_counts.mean()
sigsq = nc_counts.var()
nhat_mm = mu**2/(sigsq - mu)
phat_mm = mu/sigsq

pm = nbinom.pmf(barcode_cnts, nhat_mm, phat_mm)

mu = nc_counts.mean()


#print("max(nc_counts): ", max(nc_counts))
#print("min(nc_counts): ", min(nc_counts))
nbins = max(nc_counts)-min(nc_counts)
print("nbins: ", nbins)
if nbins == 0:
    nbins =1

plt.figure()

counts, bins, bc_ob = plt.hist(nc_counts, bins=nbins, label="observed")

u_false_counts, ncws = np.unique(nc_counts, return_counts=True)

false_counts = list(range(int(max(u_false_counts))+1))
n_neg_ctrl_cws = np.zeros(len(false_counts))
n_neg_ctrl_cws[u_false_counts] = ncws

expected = nbinom.pmf(barcode_cnts, nhat_mm, phat_mm)*len(nc_counts)
poisson_expected = poisson.pmf(barcode_cnts, mu)*len(nc_counts)

pdt_p = poisson_dispersion_test(nc_counts)


#plt.plot(false_counts, n_neg_ctrl_cws, label="observed")
plt.plot(false_counts, expected, '.', label="Negative Binomial")
plt.plot(false_counts, poisson_expected, '.', label="Poisson")

plt.title("phat = %.2f, nhat = %.2f, mu = %.2f, pdt_p = %.2f" % (phat_mm, nhat_mm, mu, pdt_p))
plt.legend()
plt.xlabel("Number of Barcodes Found")
plt.ylabel("Number of Negative Control Codewords")
plt.savefig(snakemake.output[3])

#nc_counts.to_csv(snakemake.output[4], index=False)
#cell_nc_counts.to_csv(snakemake.output[5], index=True)
