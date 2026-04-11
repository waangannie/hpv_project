import numpy as np
import sciris as sc
import hpvsim as hpv
import pandas as pd
import pylab as plt

pars = dict(
    n_agents  = 10e3,             # population size
    genotypes = [16, 18, 'hr'],   # HPV16, HPV18, plus all other high-risk genotypes pooled
    start     = 1980,             # simulation start year
    end       = 2030,             # simulation end year
    dt        = 0.25,             # timestep (quarterly)
    location  = 'nigeria',        # country for demographic data
    verbose   = 0,                # suppress print output
)

sim_baseline = hpv.Sim(pars, label='Baseline') ##fix - label does not show
sim_baseline.run()
sim_baseline.plot()


##change to beta sensitivity sweep ##

# beta = per-contact transmission probability
beta_values = np.linspace(0.1, 0.8, 5)  # sweep 5 values from 0.1 to 0.8
sims_beta = []
 
for b in beta_values:
    sim = hpv.Sim(pars, beta=b, label=f'beta={b:.2f}') #label is for formatting
    sims_beta.append(sim)
 
msim_beta = hpv.MultiSim(sims_beta)
msim_beta.run()
msim_beta.plot(['infections', 'hpv_prevalence', 'cancers', 'cancer_incidence'])


##age-stratified disease burden ##

az = hpv.age_results(
    result_args=sc.objdict(
        hpv_prevalence=sc.objdict(
            years=2020,
            edges=np.array([0., 15., 20., 25., 30., 40., 45., 50., 55., 65., 100.]),
        ),
        hpv_incidence=sc.objdict(
            years=2020,
            edges=np.array([0., 15., 20., 25., 30., 40., 45., 50., 55., 65., 100.]),
        ),
        cancer_incidence=sc.objdict(
            years=2020,
            edges=np.array([0., 20., 25., 30., 40., 45., 50., 55., 65., 100.]),
        ),
        cancer_mortality=sc.objdict(
            years=2020,
            edges=np.array([0., 20., 25., 30., 40., 45., 50., 55., 65., 100.]),
        ),
    )
)
 
sim_age = hpv.Sim(pars, analyzers=[az], label='Age-stratified analysis')##fix - label does not show
sim_age.run()
 
# Plot age-stratified results
a = sim_age.get_analyzer()
a.plot()


##uncertainty analysis ##
sim_uncertainty = hpv.Sim(pars, label='Uncertainty analysis')##fix - label does not show
msim_uncertainty = hpv.MultiSim(sim_uncertainty)
msim_uncertainty.run(n_runs=6)  # 6 runs with different random seeds
 
# Plot individual runs
msim_uncertainty.plot('hpv_prevalence')
plt.title('HPV Prevalence — Stochastic Uncertainty (6 Seeds)') ##fix - title does not show

# Plot mean with uncertainty band
msim_uncertainty.mean()
msim_uncertainty.plot('hpv_prevalence')
plt.title('HPV Prevalence — Mean ± Uncertainty') ##fix - title shows for empty graph
 
msim_uncertainty.mean()
msim_uncertainty.plot('cancer_incidence')
plt.title('Cancer Incidence — Mean ± Uncertainty') ##fix - title shows for empty graph