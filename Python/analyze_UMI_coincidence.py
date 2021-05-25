import numpy as np
import os

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

from utils import *


def analyze_UMI_coincidence(path_10X, path_probe, outdir):

    _, CBs_10X, probe_name_10X = load_10x_h5_matrix(
        f'{path_10X}/filtered_feature_bc_matrix.h5')
    _, CBs_Probe, probe_name_Probe = load_10x_h5_matrix(
        f'{path_probe}/filtered_feature_bc_matrix.h5')
    assert np.all(probe_name_10X == probe_name_Probe)
    N_genes = len(probe_name_10X)

    CBs = CBs_10X[np.in1d(CBs_10X, CBs_Probe)]
    print(f'Analysis for {len(CBs)} common cells')

    tenX = parse_molecule_info(f'{path_10X}/molecule_info.h5', CBs, N_genes)
    probe = parse_molecule_info(f'{path_probe}molecule_info.h5', CBs, N_genes)

    joint_ids = np.array(sorted(
        set(zip(tenX['barcode_id'], tenX['gene_id'])).union(
            set(zip(probe['barcode_id'], probe['gene_id'])))
    ))

    tenX = populate_common(tenX, joint_ids)
    probe = populate_common(probe, joint_ids)

    coincidence_by_umis = np.zeros(
        shape=(max(probe['umi_pop'])+1, max(tenX['umi_pop'])+1),
        dtype=int,
    )
    coincidence_by_reads = np.zeros(
        shape=(max(probe['umi_pop'])+1, max(tenX['umi_pop'])+1),
    )
    indices = (probe['umi_pop'], tenX['umi_pop'])
    np.add.at(coincidence_by_umis, indices, 1)
    np.add.at(coincidence_by_reads, indices,
              probe['read_pop']+tenX['read_pop'])

    os.makedirs(outdir, exist_ok=True)
    np.save(f'{outdir}/coincidence_by_umis.npy', coincidence_by_umis)
    np.save(f'{outdir}/coincidence_by_reads.npy', coincidence_by_reads)

    N_trunc = min(11, min(coincidence_by_umis.shape))

    save_bargraph(coincidence_by_umis[:N_trunc, :N_trunc], outdir)


def populate_common(out, joint_ids):

    r, c = out['umi_counts'].nonzero()
    v = out['umi_counts'][(r, c)]
    indices = list(zip(r, c))
    _, where = _in2d(indices, joint_ids)
    umi_pop = np.zeros(len(joint_ids), dtype=int)
    umi_pop[where] = v
    out['umi_pop'] = umi_pop

    r, c = out['read_counts'].nonzero()
    v = out['read_counts'][(r, c)]
    indices = list(zip(r, c))
    _, where = _in2d(indices, joint_ids)
    read_pop = np.zeros(len(joint_ids))
    read_pop[where] = v / out['umi_pop'][where]
    out['read_pop'] = read_pop

    return out
