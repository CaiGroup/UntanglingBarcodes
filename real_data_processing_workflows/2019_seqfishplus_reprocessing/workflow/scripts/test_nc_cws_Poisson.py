# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 13:15:24 2022

@author: jonat
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

cell_nc_cws = pd.read_csv(snakemake.input[0])
cell_nc_cws.set_index(["gene_number", "ch"], inplace=True)

ncells = len(cell_nc_cws.groupby(["pos", "cellid"]))

"""
accepts a dataframe of subdataframe containing the number of counts
that a negative control codeword was found in each each cell for cells
in which there were non-zero counts. adds zeros, and returns as numpy array.
"""
def get_cell_nc_cnts_array(df):
    cnts = df["0"]
    cnts = np.array(cnts)
    cnts = np.append(cnts, [0]*(ncells-len(cnts)))
    return cnts

"""
accepts an array of counts, and performs a Poisson dispersion test.
returns a p value where the null hypothesis is that the data is Poisson distributed.
"""
def poisson_dispersion_test(cnts):
    chi2_score = ncells*np.var(cnts)/np.mean(cnts)
    return 1 - stats.chi2.cdf(chi2_score, ncells-1)

"""
accepts a dataframe of subdataframe containing the number of counts
that a negative control codeword was found in each each cell for cells
in which there were non-zero counts. Performs a poisson dispersion test to
check if the distribution is Poisson.
"""
def nccw_poisson_dispersion_test(df):
    cnts = get_cell_nc_cnts_array(df)
    return poisson_dispersion_test(cnts)

#def cell_poisson_
    #return

def get_cw_mu_var(df):
    cnts = get_cell_nc_cnts_array(df)
    return pd.DataFrame({'mu': [cnts.mean()], 'var':[cnts.var()]})

cw_mus_vars = cell_nc_cws.groupby(["gene_number", "ch"]).apply(get_cw_mu_var)
cw_mus_vars.sort_values(by='mu', inplace=True)

pvals = cell_nc_cws.groupby(["gene_number", "ch"]).apply(nccw_poisson_dispersion_test)
nccw_cnts = cell_nc_cws.groupby(["gene_number", "ch"]).apply(get_cell_nc_cnts_array)


cnts = np.array(cell_nc_cws["0"])
cnts = np.append(cnts, [0]*(14000*ncells-len(cnts)))
plt.hist(cnts,bins=max(cnts))#+1)
#plt.rcParams.update({"text.usetex": True})
plt.plot(list(range(max(cnts)+1)), stats.poisson.pmf(list(range(max(cnts)+1)), np.mean(cnts))*len(cnts), '.',label="Poisson")
#means = cell_nc_cws.groupby(["gene_number", "ch",]).mean()
plt.title("Poisson Dispersion Test: p = %.2f, mu = %.2f, var = %.2f" % (poisson_dispersion_test(cnts), cnts.mean(), cnts.var()))
plt.ylabel("Number of negative control codewords-cells")
plt.xlabel("Number of Barcodes")
plt.legend()
plt.savefig(snakemake.output[0])

fig = plt.figure()
plt.hist(pvals,bins=20)
plt.title("Poisson Dispersion Test for negctrl bc count in cells")
plt.ylabel("Number of Negative Control Codewords")
plt.xlabel("p-value")
plt.savefig(snakemake.output[1])

fig = plt.figure()
print(cw_mus_vars)
print(cw_mus_vars['mu'])
plt.semilogy(list(cw_mus_vars['mu']),'.',label="mean", markersize=0.1)
plt.semilogy(list(cw_mus_vars['var']),'.',label="variance", markersize=0.1)
plt.title("Codeword cell barcode count")
plt.xlabel("Codeword")
plt.ylabel("mean/varaince of barcode counts in cells")
plt.legend()
plt.savefig(snakemake.output[2])


"""
cnts = np.array(cell_nc_cws["0"])
cnts = np.append(cnts, [0]*(14000*ncells-len(cnts)))
plt.hist(cnts,bins=max(cnts)+1)
plt.plot(list(range(max(cnts)+1)), stats.poisson.pmf(list(range(max(cnts)+1)), np.mean(cnts))*len(cnts),label="Poisson")
#means = cell_nc_cws.groupby(["gene_number", "ch",]).mean()
plt.title("Poisson Dispersion Test : p = %.2f" % p)
plt.ylabel("Number of cellsnegative control codewords-cells")
plt.xlabel("Number of Barcodes")
plt.legend()
plt.savefig(snakemake.output[0])
"""
