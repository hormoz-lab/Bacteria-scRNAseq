from collections import Counter, defaultdict

import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import jensenshannon
from scipy.optimize import minimize


## Helper functions
# for sampling 
def _increase_via_binomial(table, n_pcr_cycle, rng):
    """simulate imperfect pcr cycles with 0.8 success rate
    count/read of each UMI is increased via binomial sampling
    """
    reads = np.array(list(table.values()))
    for cycle in range(n_pcr_cycle):
        reads += rng.binomial(reads, p=0.8)
    return {k: v for (k, v) in zip(table.keys(), reads)} 

def _decrease_via_poisson(table, bottleneck, rng):
    """simulate dowmsampling/bottleneck process
    count/read of each UMI is sampled from poisson
    """
    reads = np.array(list(table.values()))
    bottleneck_samples = rng.poisson(reads * bottleneck)
    return {k: v for (k, v) in zip(table.keys(), bottleneck_samples)}

# for computing jensen-shannon divergence 
def _equalize_length(dist1, dist2):
    "make the number of bins to be the same for both distribution"
    if len(dist1) > len(dist2):
        _density = dist2.copy()
        for _ in range(len(dist1) - len(dist2)):
            _density = np.append(_density, 0)
        return dist1, _density
    elif len(dist2) > len(dist1):
        _density = dist1.copy()
        for _ in range(len(dist2) - len(dist1)):
            _density = np.append(_density, 0)  
        return _density, dist2
    else:
        return dist1, dist2

# computing normalization factor
def _calculate_n_total_umi_per_transcript(table):
    hit_per_transcript = defaultdict(int)
    for umi_type, read in table.items():
        transcript_index = int(umi_type.split('-')[0])
        if read > 0:
            hit_per_transcript[transcript_index] += 1
        else:
            hit_per_transcript[transcript_index] += 0
    n_total_umi_type = []
    return list(hit_per_transcript.values())


## Simulation suite
class Simulation:
    def __init__(self, 
        density_exp,
        n_transcript, 
        n_pcr_trial,
        n_subsampling_trial,
        p, 
        rng=None,
    ):
        self.density_exp = density_exp # density histogram from experiment
        self.n_transcript = n_transcript
        self.n_pcr_trial = n_pcr_trial
        self.n_subsampling_trial = n_subsampling_trial
        self.p = p
        self.rng = np.random.default_rng() if rng is None else rng
        self.result = {}

    def simulate_ingem_pcr(self, cycles=6):
        """simulate 6 rounds of ingem PCR with copying success probability 
        given by p. Whenever a reverse strand is copied by the forward 10X
        primer, a new UMI is introduced. Each read is represented as a tuple 
        of UMI-type and strand-direction. 
        """
        print(f'Starting {cycles} cycles of ingem (in-droplet) PCR')
        counters = []   # histogram of UMI counts
        for _ in range(self.n_transcript):
            umis_pool = list(range(2,2000))     # pool of UMI-type 
            reads = [(1, 'reverse')] # probe with poly-A bound to a transcript
            for i in range(cycles):
                reads_next = []
                for read in reads:
                    if self.rng.binomial(1, p=self.p): 
                        if len(reads) == 1:  # first round
                            new_read = (read[0], 'forward')
                            reads_next += [read, new_read]
                        else:   # rest 
                            if read[1] == 'reverse':
                                umi_chosen = self.rng.choice(umis_pool)
                                umis_pool.remove(umi_chosen)
                                new_read = (umi_chosen, 'forward')
                                reads_next += [read, new_read]
                            else:
                                new_read = (read[0], 'reverse')
                                reads_next += [read, new_read]
                    else:
                        reads_next += [read]
                reads = reads_next  # update existing reads with PCR results
            umi_type = [read[0] for read in reads]
            counter = Counter(umi_type)
            bins = [1, 2, 3, 4, 5, 6, 7, 8]
            hist, bin_edges = np.histogram(list(counter.values()), bins)
            histogram = {k: v for k, v in zip(bin_edges[:-1], hist)}
            counters.append(histogram)
            table = self._populate_table(counters)
        print('ingem (in-droplet) PCR ended')
        self.result['post_ingem'] = table

    def simulate_normal_PCR(self, cycles=16):
        print(f'Starting {self.n_pcr_trial} of {cycles} cycles of out-of-droplet PCR')
        tables = [
            _increase_via_binomial(self.result['post_ingem'], cycles, self.rng)
            for _ in range(self.n_pcr_trial)
        ]
        print('out-of-droplet PCR ended')
        self.result['post_pcr'] = tables

    def _populate_table(self, counters):
        types = []
        mRNA_index = np.arange(1, self.n_transcript+1).astype(str)
        for mRNA, counter in zip(mRNA_index, counters):
            reads_index = np.arange(1, sum(counter.values())+1).astype(str)
            types.append([f'{mRNA}-{read}' for read in reads_index])
                    
        n_reads = []
        for counter in counters:
            n_reads.append([1]*counter[1] + 
                           [2]*counter[2] + 
                           [3]*counter[3] + 
                           [4]*counter[4] + 
                           [5]*counter[5] + 
                           [6]*counter[6] + 
                           [7]*counter[7])
            
        table = {}
        for _type, _read in zip(types, n_reads):
            for _type_, _read_ in zip(_type, _read):
                table[_type_] = _read_
        return table

    def plot(self, exp):
        fig, ax = plt.subplots()
        if exp == 'post_ingem':
            bins = [1, 2, 3, 4, 5, 6, 7, 8]
            read_count = list(self.result[exp].values())
            hist, bin_edges = np.histogram(read_count, bins, density=True)
            ax.bar(bin_edges[:-1], hist)
            ax.set_title('after in-droplet PCR')
        elif exp == 'post_pcr':
            read_count = []
            for table in self.result[exp]:
                read_count += list(filter(lambda x: x > 0, table.values()))
            ax.hist(read_count, bins='auto', density=True)
            ax.set_title('after PCR #2')
        ax.set_xlabel('number of reads with the same UMI')
        ax.set_ylabel('normalized count')
        fig.show()

    def _subsample(self, r, n_trial):
        post_bottleneck = [
            _decrease_via_poisson(table, r, self.rng)
            for _ in range(n_trial)
            for table in self.result['post_pcr']
        ]
        return post_bottleneck


    def js_objective(self, r, n_trial=None):
        n_trial = n_trial if n_trial else self.n_subsampling_trial
        # 1. simulate subsampling
        post_bottleneck = self._subsample(r, n_trial)

        # 2. make density histogram of reads
        read_count = []
        for table in post_bottleneck:
            read_count += list(filter(lambda x: x > 0, table.values()))
        bins = [i for i in range(1, max(read_count)+4)]  # will get cutoff...
        density_sim, bin_edges = np.histogram(read_count, bins=bins, density=True)

        # 2.5. preprocess simulation and experimental density histogram
        density_sim = density_sim[density_sim>0]
        density_sim = density_sim[:len(self.density_exp)] # ...here
        density_sim, density_exp = _equalize_length(density_sim, self.density_exp)
        density_sim /= density_sim.sum()
        _density_exp = self.density_exp.copy()
        _density_exp /= _density_exp.sum()

        # 3. compute Jensen-Shannon divergence
        dist = jensenshannon(density_sim, _density_exp)

        return dist

    def fit_r(self, initial_r, n_opt_trial):
        self.opt_result = {'dist': [], 'r': []}
        for _ in range(n_opt_trial):
            res = minimize(
                self.js_objective, 
                np.array([initial_r]),
                method='nelder-mead',
                options={'xatol': 1e-8, 'disp': False})
            self.opt_result['dist'].append(res.fun)
            self.opt_result['r'].append(res.x[0])
        self.r = np.mean(self.opt_result['r'])

    def sweep_r(self, low_r, high_r, n_eval=10, n_trial=2):
        n_trial = n_trial if n_trial else self.n_subsampling_trial
        r_array = np.linspace(low_r, high_r, n_eval)
        dists = [self.js_objective(r, n_trial) for r in r_array]
        fig, ax = plt.subplots()
        ax.plot(r_array, dists, marker='o')
        ax.set_ylabel('J.S. distance')
        ax.set_xlabel('subsampling factor') 
        fig.show()

    def calculate_normalization_factor(self, n_trial=None):
        n_trial = n_trial if n_trial else self.n_subsampling_trial
        post_bottleneck = self._subsample(self.r, n_trial)
        n_total_umi = []
        for table in post_bottleneck:
            n_total_umi += _calculate_n_total_umi_per_transcript(table)
        self.normalization_factor = np.mean(n_total_umi)
