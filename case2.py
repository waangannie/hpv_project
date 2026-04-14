import hpvsim as hpv
import pandas as pd
import numpy as np

pars = dict(
    n_agents=10e3,
    start=1950,
    end=2060,
    dt=0.25, #quarterly
    location='nigeria',
)

def make_vaccine(coverage, start_year):
    vaccine = hpv.routine_vx(
        prob=coverage,           # vaccination probability per eligible person
        start_year=start_year,   # rollout start
        product='bivalent',      # bivalent vaccine covers HPV16 + HPV18
        age_range=(9, 14),       # target age group
        label=f'{int(coverage*100)}% coverage from {start_year}'
    )
    return vaccine

scenarios = {
    'No vaccination':       None,
    '30% coverage (2010)':  make_vaccine(0.30, 2010),
    '70% coverage (2010)': make_vaccine(0.7,2010)
}

sims = []
for label, intervention in scenarios.items():
    sim = hpv.Sim(
        pars,
        interventions=intervention,
        label=label
    )
    sim.run()
    sims.append(sim)

msim = hpv.MultiSim(sims)
msim.plot(['cancers', 'hpv_prevalence', 'cins'])