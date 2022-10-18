import numpy as np
import pymc3 as pm
import tqdm
import multiprocess as mp
import warnings
import logging

warnings.filterwarnings("ignore")
logger = logging.getLogger("pymc3")
logger.setLevel(logging.ERROR)

def generate_posterior(grouped_bin_counts, prior, reps, nchain, bins_pool, sd, ncores):
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
            trace = pm.sample(reps, tune=reps/2, chains=nchain, cores=1, random_seed=sd, progressbar=False, target_accept=0.99)

        ppc = pm.fast_sample_posterior_predictive(trace, samples=reps, model=model_cnv, random_seed=sd)
        return ppc.get('y')

    if __name__ == "__main__":
        with mp.Pool(ncores) as p:
            posteriors = np.empty((len(window_nums_uniq), 1000))
            for _ in tqdm.tqdm(p.imap(single_post, window_nums_uniq), total=len(window_nums_uniq)):
                posteriors = np.vstack([posteriors, _])
                pass

    return posteriors