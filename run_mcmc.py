import sys
import numpy as np
import pandas as pd
import pymc3 as pm
import tqdm
from multiprocessing import Pool
import warnings
import logging
import json

warnings.filterwarnings("ignore")
logger = logging.getLogger("pymc3")
logger.setLevel(logging.ERROR)

gbc_path = sys.argv[1]
pr_path = sys.argv[2]
reps = int(sys.argv[3])
nchain = int(sys.argv[4])
bins_pool = int(sys.argv[5])
seed = int(sys.argv[6])
ncores = int(sys.argv[7])
outdir = sys.argv[8]

post_outpath = outdir + "/posteriors.tmp"

with open(pr_path) as f:
        prior = f.read()
prior = json.loads(prior)

grouped_bin_counts = pd.read_table(gbc_path, sep="\t")

window_nums = grouped_bin_counts['window'].tolist()
window_nums_uniq = list(set(window_nums))

def single_post(window):
    dat = grouped_bin_counts.loc[grouped_bin_counts['window'] == window, :]

    if dat.shape[0] < bins_pool:
        dat = dat.sample(bins_pool, replace=True)

    with pm.Model() as model_cnv:
        mu = pm.Lognormal('mu', mu=prior.get('mu')[0], sigma=prior.get('mu')[1])
        alpha = pm.Weibull('alpha', alpha=prior.get('size')[0], beta=prior.get('size')[1])
        y = pm.NegativeBinomial('y', mu=mu, alpha=alpha, observed=dat['RC_GCcorr_count'].tolist())
        trace = pm.sample(reps, tune=reps/2, chains=nchain, cores=1, random_seed=seed, progressbar=False, target_accept=0.99)

    ppc = pm.fast_sample_posterior_predictive(trace, samples=1, model=model_cnv, random_seed=seed)
    return ppc.get('y')

if __name__ == "__main__":
    p = Pool(ncores)
    posteriors = list(tqdm.tqdm(p.imap(single_post, window_nums_uniq), total=len(window_nums_uniq)))
    posteriors = np.matrix(np.squeeze(np.array(posteriors))).astype(int)
print(posteriors.shape)

with open(post_outpath, 'w') as f:
    for line in posteriors:
        np.savetxt(f, line, fmt='%.0f')