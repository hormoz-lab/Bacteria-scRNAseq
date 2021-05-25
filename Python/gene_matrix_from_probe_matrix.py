import os
import numpy as np

from utils import *


def gene_matrix_from_probe_matrix(input_path, outdir, map_method):

    assert not os.path.isdir(outdir), 'Output path already exists'

    if map_method not in ['max_prob_in_cell', 'median_probe_in_bulk']:
        raise Exception(f'Unrecorgnized map_method: {map_method}')

    counts, CBs, probe_name = load_10x_h5_matrix(
        'filtered_feature_bc_matrix_10X.h5')
    N_genes = len(probe_name)
    N_cells = len(CBs)
    print(f'Analysis for {len(CBs)} common cells using method={map_method}')

    data = parse_molecule_info('molecule_info_10X.h5', CBs, N_genes)
    assert np.all(data['umi_counts'].todense(), counts.todense())
    del counts

    gene_name = [str.split(b'_')[0].lower() for str in probe_name]
    probe_id = [int(str.split(b'_')[1]) for str in probe_name]
    gene_name, which_gene = np.unique(gene_name, return_inverse=True)

    reads_per_cell = np.sum(data['read_counts'], axis=1)
    umis_per_cell = np.sum(data['umi_counts'], axis=1)
    probes_per_cell = np.count_nonzero(data['umi_counts'].toarray(), axis=1)
    _genes_per_cell = np.zeros(
        (max(data['barcode_id'])+1, max(data['gene_id'])+1))
    np.add.at(_genes_per_cell,
              (data['barcode_id'], which_gene[data['gene_id']]), 1)
    genes_per_cell = np.count_nonzero(_genes_per_cell, axis=1)

    gene_idx = [np.argwhere(which_gene == i).flatten()
                for i in range(len(gene_name))]
    sum_probe = np.array(
        [data['umi_counts'][:, idx].sum(axis=1) for idx in gene_idx]
    ).squeeze().T

    counts = {
        'by_probe': data['umi_counts'],
        'by_gene_total': sum_probe,
    }

    if map_method == 'max_probe_in_cell':

        which_probe_max = np.array(
            [data['umi_counts'][:, idx].argmax(axis=1) for idx in gene_idx]
        ).squeeze().T

        max_umi_reads = np.array(
            [data['umi_counts'][:, idx][np.arange(N_cells), max_idx]
             for (idx, max_idx) in zip(gene_idx, which_probe_max.T)]
        ).squeeze().T

        max_probe_reads = np.array(
            [data['read_counts'][:, idx][np.arange(N_cells), max_idx]
             for (idx, max_idx) in zip(gene_idx, which_probe_max.T)]
        ).squeeze().T

        counts['by_gene_max_probe'] = max_probe_umis,
        counts['which_probe_max'] = which_probe_max,
        counts['max_probe_reads'] = max_probe_reads

    else:
        bulk_gene = [
            np.array(data['umi_counts'][:, idx].sum(axis=0)).flatten()
            for idx in gene_idx
        ]
        which_probe_bulk_median = [
            y[np.argmax(np.abs(x-np.median(x)) == min(np.abs(x-np.median(x))))]
            for (x, y) in zip(bulk_gene, gene_idx)
        ]
        bulk_median_probe_reads = data['read_counts'][:,
                                                      which_probe_bulk_median]

        counts['by_median_bulk_probe'] = data['umi_counts'][:,
                                                            which_proobe_bulk_median],
        counts['which_probe_bulk_median'] = which_probe_bulk_median,
        counts['bulk_median_reads'] = bulk_median_probe_reads

    return counts, data, CBs, gene_name
