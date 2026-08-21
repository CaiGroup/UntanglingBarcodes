import pandas as pd
import numpy as np


def load_codebook(cb_path):
    cb = pd.read_csv(cb_path)
    gene_col = cb.columns[0]
    block_cols = [c for c in cb.columns if c != gene_col]
    matrix = cb[block_cols].values
    n_rounds = len(block_cols)
    q = len(np.unique(matrix))
    return cb, matrix, n_rounds, q, gene_col


def to_bardensr(cb_path, zeros_probed=False):
    """Returns codebook as (N_images, J) indicator matrix for BarDenSR.

    zeros_probed: pass True for codebook families (e.g. seqFISH) where every
    round always carries a real dot and symbol 0 is a valid pseudocolor (the
    Z_n group identity), not "no dot". Symbol 0 is then rendered/matched as
    pseudocolor q (see gen_sim_data.jl and get_cand_cpaths.jl)."""
    cb, matrix, n_rounds, q, _ = load_codebook(cb_path)
    n_codewords = matrix.shape[0]
    if q == 2:
        return matrix.astype(float).T
    n_channels = q if zeros_probed else q - 1
    n_images = n_rounds * n_channels
    indicator = np.zeros((n_images, n_codewords), dtype=float)
    for j in range(n_codewords):
        for r in range(n_rounds):
            symbol = int(matrix[j, r])
            if zeros_probed:
                pc = symbol if symbol != 0 else q
                frame_idx = r * n_channels + (pc - 1)
                indicator[frame_idx, j] = 1.0
            elif symbol != 0:
                frame_idx = r * n_channels + (symbol - 1)
                indicator[frame_idx, j] = 1.0
    return indicator


def to_istdeco(cb_path, zeros_probed=False):
    cb, matrix, n_rounds, q, _ = load_codebook(cb_path)
    n_codewords = matrix.shape[0]
    if q == 2:
        n_channels = 1
    else:
        n_channels = q if zeros_probed else q - 1
    D = np.zeros((n_codewords, n_rounds, n_channels), dtype=float)
    for j in range(n_codewords):
        for r in range(n_rounds):
            symbol = int(matrix[j, r])
            if q == 2:
                if symbol != 0:
                    D[j, r, 0] = 1.0
            elif zeros_probed:
                pc = symbol if symbol != 0 else q
                D[j, r, pc - 1] = 1.0
            else:
                if symbol != 0:
                    D[j, r, symbol - 1] = 1.0
    norms = D.sum(axis=(1, 2), keepdims=True)
    norms[norms == 0] = 1.0
    D = D / norms
    return D


def to_graph_iss_taglist(cb_path):
    cb, matrix, n_rounds, q, gene_col = load_codebook(cb_path)
    taglist = {}
    for j in range(matrix.shape[0]):
        barcode = tuple(int(matrix[j, r]) for r in range(n_rounds))
        gene = int(cb.iloc[j][gene_col])
        taglist[barcode] = gene
    return taglist
