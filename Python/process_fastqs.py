import os
import gzip
import shutil
import time
import pickle
from subprocess import check_output
from datetime import timedelta

import numpy as np

import parasail

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from utils import *


def process_fastqs(input_fastq, outdir, tenX_version):

    start_time = time.time()

    assert os.path.isfile(input_fastq), 'Invalid path to 10x counts file'

    assert not os.path.isdir(outdir), 'Output path already exists'

    os.makedirs(outdir)
    os.makedirs(f'{outdir}/10X')
    os.makedirs(f'{outdir}/Probe')
    os.makedirs('meta')

    print(f'Processing {input_fastq}')

    complement = {
        'A': 'T',
        'C': 'G',
        'G': 'C',
        'T': 'A',
    }
    Extender = ''.join(
        [complement.get(nuc) for nuc in reversed('CATAGTTTCTTCGAGCAA')]
    )
    L_ext = len(Extender)
    L_CB = 16
    if tenX_version == 2:
        L_UMI = 10
    elif tenX_version == 3:
        L_UMI = 12
    else:
        raise ValueError('Unrecognized 10X version')

    if input_fastq.endswith('gz'):
        print('Unzipping FASTQ file')
        with gzip.open(input_fastq, 'rt') as f_in:
            with open(input_fastq[:-3], 'wt') as f_out:
                shutil.copyfileobj(f_in, f_out)
        input_fastq = input_fastq[:-3]

    command = f'wc -l {input_fastq}'
    N_reads = int(int(check_output(command, shell=True).split()[0])/4)
    block_size = 1e7
    N_blocks = int(np.ceil(N_reads/block_size))
    print(f'Splitting {N_reads} over {N_blocks}')

    record_iter = SeqIO.parse(input_fastq, 'fastq')

    L_orig = []
    masks = {
        'full_seq': [],
        'extender_found': [],
        'valid_probe': [],
        'valid_seq': [],
    }
    probe_umi_end_pos = []
    sc = []

    for records in batch_iterator(record_iter, block_size):
        _L, _ma, _pr, _sc = \
            process_per_block(records, Extender, L_CB, L_UMI, L_ext, outdir)
        L_orig.append(_L)
        for k in masks.keys():
            masks[k].append(_ma[k])
        probe_umi_end_pos.append(_pr)
        sc.append(_sc)

    L_orig = np.hstack(L_orig)
    for k in masks.keys():
        masks[k] = np.hstack(masks[k])
    probe_umi_end_pos = np.hstack(probe_umi_end_pos)
    sc = np.hstack(probe_umi_end_pos)

    command = f'wc -l {outdir}/10X/R2_001.fastq'
    N_reads_processed = int(
        int(check_output(command, shell=True).split()[0])/4)

    print('Removing unzipped FASTQ file')
    os.remove(input_fastq)

    assert len(masks['valid_seq']) == N_reads_processed

    print('Output and compress R2...')

    with open(f'{outdir}/10X/R2_001.fastq', 'rt') as f_in:
        with gzip.open(f'{outdir}/Probe/R2_001.fastq.gz', 'wt') as f_out:
            shutil.copyfileobj(f_in, f_out)
    with open(f'{outdir}/10X/R2_001.fastq', 'rt') as f_in:
        with gzip.open(f'{outdir}/10X/R2_001.fastq.gz', 'wt') as f_out:
            shutil.copyfileobj(f_in, f_out)
    os.remove(f'{outdir}/10X/R2_001.fastq')

    print('Compress 10X R1...')

    with open(f'{outdir}/10X/R1_001.fastq', 'rt') as f_in:
        with gzip.open(f'{outdir}/10X/R1_001.fastq.gz', 'wt') as f_out:
            shutil.copyfileobj(f_in, f_out)
    os.remove(f'{outdir}/10X/R1_001.fastq')

    print('Compress Probe R1...')

    with open(f'{outdir}/Probe/R1_001.fastq', 'rt') as f_in:
        with gzip.open(f'{outdir}/Probe/R1_001.fastq.gz', 'wt') as f_out:
            shutil.copyfileobj(f_in, f_out)
    os.remove(f'{outdir}/Probe/R1_001.fastq')

    print('Output provenance...')

    metric = {
        'L_orig': L_orig,
        'masks': masks,
        'probe_umi_end_pos': probe_umi_end_pos,
        'sc': sc,
    }

    pickle.dump(metric, open('meta/meta.pkl', 'wb'))

    time_elapsed = time.time() - start_time
    print(f'FASTQ processed in: {timedelta(seconds=time_elapsed)}')


# def batch_iterator(iterator, batch_size):
#    entry = True  # Make sure we loop once
#    while entry:
#        batch = []
#        while len(batch) < batch_size:
#            try:
#                entry = next(iterator)
#            except StopIteration:
#                entry = None
#            if entry is None:
#                # End of file
#                break
#            batch.append(entry)
#        if batch:
#            yield batch


def process_per_block(records, Extender, L_CB, L_UMI, L_ext, outdir):
    S = np.array([rec.seq for rec in records], dtype=object)

    print('Trimming by length...')
    L_orig = np.array([len(seq) for seq in S])
    full_seq = np.where(L_orig > 110)[0]
    S = S[full_seq]

    print('Aligning...')
    sc, al = [], []
    for s in S:
        result = parasail.sg_trace(s._data, Extender, 8, 8, parasail.nuc44)
        sc.append(0.277316 * result.score)
        al.append(result.cigar.decode.split(b'I')[0])
    extender_found = np.where(np.array(sc) >= 17)[0]
    al = np.array(al)[extender_found]
    S = S[extender_found]
    extender_found = full_seq[extender_found]

    print('Trimming by extender...')
    probe_umi_end_pos = np.array([int(match) for match in al])
    # Should be a polyT region separating 10X (CB/UMI) and Probe UMI
    probe_mask = np.where(probe_umi_end_pos >= L_CB + 2*L_UMI)[0]
    S = S[probe_mask]
    valid_probe = extender_found[probe_mask]

    # Keep the extender since the sequences are ~30 bps without it which
    # may affect alignment, but remove sequences consisting of just the
    # extender
    seq_mask = np.where(L_orig[valid_probe] >
                        probe_umi_end_pos[probe_mask] + L_ext)[0]
    S = S[seq_mask]
    valid_seq = valid_probe[seq_mask]

    print('Splitting S...')
    S_CB = np.array([seq[:L_CB]._data for seq in S])
    S_10X_UMI = np.array([seq[L_CB:L_CB+L_UMI]._data for seq in S])
    S_probe_UMI = np.array(
        [seq[e-L_UMI:e]._data
         for (seq, e) in zip(S, probe_umi_end_pos[probe_mask[seq_mask]])]
    )
    S = [seq[s:]._data
         for (seq, s) in zip(S, probe_umi_end_pos[probe_mask[seq_mask]])]

    print('Reading Q...')
    Q = np.array(
        [rec.letter_annotations['phred_quality'] for rec in records],
        dtype=object,
    )
    Q = Q[valid_seq]

    print('Splitting Q...')
    Q_CB = np.array([x[:L_CB] for x in Q])
    Q_10X_UMI = np.array([x[L_CB:L_CB+L_UMI] for x in Q])
    Q_probe_UMI = np.array(
        [x[e-L_UMI:e]
         for (x, e) in zip(Q, probe_umi_end_pos[probe_mask[seq_mask]])]
    )
    Q = [x[s:] for (x, s) in zip(Q, probe_umi_end_pos[probe_mask[seq_mask]])]

    print('Reading H...')
    H = np.array([rec.description for rec in records])
    H = H[valid_seq]

    print('Output R2...')

    output_fastq2 = []
    for (s, h, q) in zip(S, H, Q):
        record = SeqRecord(Seq(s), id=h, name=h, description='')
        record.letter_annotations['phred_quality'] = q
        output_fastq2.append(record)
    with open(f'{outdir}/10X/R2_001.fastq', 'at+') as f:
        SeqIO.write(output_fastq2, f, 'fastq')

    print('Output 10X R1...')
    S = np.char.add(S_CB, S_10X_UMI)
    Q = np.hstack((Q_CB, Q_10X_UMI))
    output_fastq1 = []
    for (s, h, q) in zip(S, H, Q):
        record = SeqRecord(Seq(s), id=h, name=h, description='')
        record.letter_annotations['phred_quality'] = q
        output_fastq1.append(record)
    with open(f'{outdir}/10X/R1_001.fastq', 'at+') as f:
        SeqIO.write(output_fastq1, f, 'fastq')

    print('Output Probe R1...')
    S = np.char.add(S_CB, S_probe_UMI)
    Q = np.hstack((Q_CB, Q_probe_UMI))
    output_fastq1 = []
    for (s, h, q) in zip(S, H, Q):
        record = SeqRecord(Seq(s), id=h, name=h, description='')
        record.letter_annotations['phred_quality'] = q
        output_fastq1.append(record)
    with open(f'{outdir}/Probe/R1_001.fastq', 'at+') as f:
        SeqIO.write(output_fastq1, f, 'fastq')

    masks = {
        'full_seq': full_seq,
        'extender_found': extender_found,
        'valid_probe': valid_probe,
        'valid_seq': valid_seq,
    }

    return L_orig, masks, probe_umi_end_pos, sc
