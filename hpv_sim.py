import hpvsim as hpv

sim = hpv.Sim()
sim.run()
sim.plot()


#running this:
pars = dict(
    n_agents = 10e3,
    genotypes = [16, 18, 'hr'], # Simulate genotypes 16 and 18, plus all other high-risk HPV genotypes pooled together
    start = 1980,
    end = 2030,
)

sim = hpv.Sim(pars)
sim.run()

#is the same as this:
sim = hpv.Sim(n_agents=10e3, start=1980, end=2030)
sim.run()

#can also mix and match creating a paramater and sim:
sim = hpv.Sim(pars, end=2050) # Use parameters defined above, except set the end data to 2050
sim.run()

#plot by:
fig = sim.plot()


# Custom vaccination intervention - "In the year 2000, find every agent aged 10–13 and set their immunity against disease 0 to maximum — as if they were just vaccinated."
def custom_vx(sim): #sim.t is the current time step index, and sim.yearvec is the array of all years in the simulation.
    if sim.yearvec[sim.t] == 2000: 
        target_group = (sim.people.age>9) * (sim.people.age<14) #multiplying them together acts like an AND 
        sim.people.peak_imm[0, target_group] = 1 #0 refers to the first disease/pathogen in the simulation. Setting it to 1 for everyone in target_group means those agents instantly get 100% peak immunity

pars = dict(
    location = 'tanzania', # Use population characteristics for Japan
    n_agents = 10e3, # Have 50,000 people total in the population
    start = 1980, # Start the simulation in 1980
    n_years = 50, # Run the simulation for 50 years
    burnin = 10, # Discard the first 20 years as burnin period
    verbose = 0, # Do not print any output
)

# Running with multisims
s1 = hpv.Sim(pars, label='Default') #baseline/control - no interventions
s2 = hpv.Sim(pars, interventions=custom_vx, label='Custom vaccination')
msim = hpv.MultiSim([s1, s2])
msim.run()
fig = msim.plot(['cancers', 'cins'])


#checking objects on level of detail in the order of ascending
sim.brief()
sim.summarize()
sim.disp()

#default plotting
sim.plot()
#Plotting w/ Matplotlib
import pylab as plt # Shortcut for import matplotlib.pyplot as plt
plt.plot(sim.results['year'], sim.results['infections']);
#plotting more than 1 figure - to_plot
sim.plot(to_plot=['infections', 'hpv_incidence']);

#save figures
hpv.savefig('my-fig.png')

#change style of plots - default style of hpvsim is Matplotlib style 
sim.plot(style='ggplot');
sim.plot(style='simple', legend_args={'frameon':True}, style_args={'ytick.direction':'in'});

#controlling plot font size
import numpy as np
with hpv.options.with_style(fontsize=6):
    sim.plot() # This will have 6 point font
    plt.figure(); plt.plot(np.random.rand(20), 'o') # So will this

#saving sims
sim.save('my-sim.sim')

#export the results and parameters to an excel file
import pandas as pd

sim.to_excel('my-sim.xlsx')
df = pd.read_excel('my-sim.xlsx')
print(df)


#Making and running a multisim based on a single sim (running one sim with uncertainty)
sim = hpv.Sim()
msim = hpv.MultiSim(sim)
msim.run(n_runs=5)
msim.plot(); #by default the multisim simply plots each simulation.

#calculate either the mean or the median of the results across all the sims
msim.mean()
msim.plot('infections');

msim.median()
msim.plot('infections');

#can also treat each of the individual sims as part of a larger single sim, and “combine” the results into one sim
msim.combine()
msim.plot('infections');

#run a set of different sims - do a sweep across relative transmissibility of HPV
import numpy as np

rel_trans_vals = np.linspace(0.1, 0.8, 5) # Sweep from 0.5 to 1.5 with 5 values
sims = []
for rel_trans in rel_trans_vals:
    sim = hpv.Sim(beta=rel_trans, label=f'Rel trans HPV = {rel_trans}')
    sims.append(sim)
msim = hpv.MultiSim(sims)
msim.run()
msim.plot('infections');

#can use multisims to do very compact scenario explorations
def custom_vx(sim):
    if sim.yearvec[sim.t] == 2000:
        target_group = (sim.people.age>9) * (sim.people.age<14)
        sim.people.peak_imm[0, target_group] = 1

pars = dict(
    location = 'tanzania', # Use population characteristics for Japan
    n_agents = 10e3, # Have 50,000 people total in the population
    start = 1980, # Start the simulation in 1980
    n_years = 50, # Run the simulation for 50 years
    burnin = 10, # Discard the first 20 years as burnin period
    verbose = 0, # Do not print any output
)

s1 = hpv.Sim(pars, label='Default')
s2 = hpv.Sim(pars, interventions=custom_vx, label='Custom vaccination')
hpv.parallel(s1, s2).plot(['hpv_incidence', 'cancer_incidence']); #the command hpv.parallel() is an alias for hpv.MultiSim().run()

#merge or split different multisims together
n_sims = 3
rel_trans_vals = [0.25, 0.5, 0.75]

msims = []
for rel_trans in rel_trans_vals:
    sims = []
    for s in range(n_sims):
        sim = hpv.Sim(n_agents=10e3, beta=rel_trans, rand_seed=s, label=f'Rel trans = {rel_trans}')
        sims.append(sim)
    msim = hpv.MultiSim(sims)
    msim.run()
    msim.mean()
    msims.append(msim)

merged = hpv.MultiSim.merge(msims, base=True)
merged.plot(color_by_sim=True);


#running w/ scenario objects - starts from the same base sim, then modify the parameters as you specify, and finally add uncertainty if desired
# Set base parameters -- these will be shared across all scenarios
basepars = {'n_agents':10e3}

# Configure the settings for each scenario
scenarios = {'baseline': {
              'name':'Baseline',
              'pars': {}
              },
            'high_rel_trans': {
              'name':'High rel trans (0.75)',
              'pars': {
                  'beta': 0.75,
                  }
              },
            'low_rel_trans': {
              'name':'Low rel trans(0.25)',
              'pars': {
                  'beta': 0.25,
                  }
              },
             }

# Run and plot the scenarios
scens = hpv.Scenarios(basepars=basepars, scenarios=scenarios)
scens.run()
scens.plot();


#Screening and treatment interventions
# Define a series of interventions to screen, triage, assign treatment, and administer treatment
prob = 0.6
screen      = hpv.routine_screening(start_year=2015, prob=prob, product='via', label='screen') # Routine screening
to_triage   = lambda sim: sim.get_intervention('screen').outcomes['positive'] # Define who's eligible for triage
triage      = hpv.routine_triage(eligibility=to_triage, prob=prob, product='hpv', label='triage') # Triage people
to_treat    = lambda sim: sim.get_intervention('triage').outcomes['positive'] # Define who's eligible to be assigned treatment
assign_tx   = hpv.routine_triage(eligibility=to_treat, prob=prob, product='tx_assigner', label='assign_tx') # Assign treatment
to_ablate   = lambda sim: sim.get_intervention('assign_tx').outcomes['ablation'] # Define who's eligible for ablation treatment
ablation    = hpv.treat_num(eligibility=to_ablate, prob=prob, product='ablation') # Administer ablation
to_excise   = lambda sim: sim.get_intervention('assign_tx').outcomes['excision'] # Define who's eligible for excision
excision    = hpv.treat_delay(eligibility=to_excise, prob=prob, product='excision') # Administer excision

# Define the parameters
pars = dict(
    n_agents      = 20e3,       # Population size
    n_years       = 35,         # Number of years to simulate
    verbose       = 0,          # Don't print details of the run
    rand_seed     = 2,          # Set a non-default seed
    genotypes     = [16, 18],   # Include the two genotypes of greatest general interest
)

# Create the sim with and without interventions
orig_sim = hpv.Sim(pars, label='Baseline')
sim = hpv.Sim(pars, interventions = [screen, triage, assign_tx, ablation, excision], label='With screen & treat')

# Run and plot
msim = hpv.parallel(orig_sim, sim)
msim.plot();


# making a custom product - Define a new treatment called new_tx that works on all HPV genotypes, with increasing effectiveness as the disease gets more severe — then package it into an HPVsim intervention object ready to use in a simulation
import hpvsim as hpv
import pandas as pd #pandas library for creating and managing data tables
my_treatment_data = pd.DataFrame({'name':'new_tx', 'state':['precin','cin1','cin2','cin3','cancerous'],'genotype':'all','efficacy':[.2,.3,.3,.4,.4]}) #creates spreadsheet like structure 
my_treatment = hpv.tx(df=my_treatment_data) #passes dataframe into hpvsim tx


#Prophylactic vaccination - targets a vaccine product towards a subset of the population
vx = hpv.routine_vx(prob=prob, start_year=2015, age_range=[9,10], product='bivalent')

# Create the sim with and without interventions
orig_sim = hpv.Sim(pars, label='Baseline')
sim = hpv.Sim(pars, interventions = vx, label='With vaccination')

# Run and plot
msim = hpv.parallel(orig_sim, sim)
msim.plot();


# therapeutic vaccination
import numpy as np #library for numerical computing
import hpvsim as hpv

# Define mass therapeutic vaccination:
campaign_txvx_dose1 = hpv.campaign_txvx(prob = 0.9, years = 2015, age_range = [30,50], product = 'txvx1', label = 'campaign txvx')
second_dose_eligible = lambda sim: (sim.people.txvx_doses == 1) | (sim.t > (sim.people.date_tx_vaccinated + 0.5 / sim['dt']))
campaign_txvx_dose2 = hpv.campaign_txvx(prob = 0.7, years=[2015,2016], age_range=[30, 70], product = 'txvx2', eligibility = second_dose_eligible, label = 'campaign txvx 2nd dose')
routine_txvx_dose1 = hpv.routine_txvx(prob = 0.9, start_year = 2016, age_range = [30,31], product = 'txvx2',label = 'routine txvx')
second_dose_eligible = lambda sim: (sim.people.txvx_doses == 1) | (sim.t > (sim.people.date_tx_vaccinated + 0.5 / sim['dt']))
routine_txvx_dose2 = hpv.routine_txvx(prob = 0.8, start_year = 2016, age_range = [30,31], product = 'txvx1', eligibility=second_dose_eligible, label = 'routine txvx 2nd dose')
mass_vx_intvs = [campaign_txvx_dose1, campaign_txvx_dose2, routine_txvx_dose1, routine_txvx_dose2]
for intv in mass_vx_intvs: intv.do_plot=False

# Define therapeutic vaccination within screen and treat
campaign_txvx_dose1 = hpv.campaign_txvx(prob = 0.9, years = 2015, age_range = [30,50], product = 'txvx1', label = 'campaign txvx')
second_dose_eligible = lambda sim: (sim.people.txvx_doses == 1) | (sim.t > (sim.people.date_tx_vaccinated + 0.5 / sim['dt']))
campaign_txvx_dose2 = hpv.campaign_txvx(prob = 0.7, years=[2015,2016], age_range=[30, 70], product = 'txvx2', eligibility = second_dose_eligible, label = 'campaign txvx 2nd dose')
routine_txvx_dose1 = hpv.routine_txvx(prob = 0.9, start_year = 2016, age_range = [30,31], product = 'txvx2',label = 'routine txvx')
second_dose_eligible = lambda sim: (sim.people.txvx_doses == 1) | (sim.t > (sim.people.date_tx_vaccinated + 0.5 / sim['dt']))
routine_txvx_dose2 = hpv.routine_txvx(prob = 0.8, start_year = 2016, age_range = [30,31], product = 'txvx1', eligibility=second_dose_eligible, label = 'routine txvx 2nd dose')
mass_vx_intvs = [campaign_txvx_dose1, campaign_txvx_dose2, routine_txvx_dose1, routine_txvx_dose2]
for intv in mass_vx_intvs: intv.do_plot=False


# Screen, triage, assign treatment, treat
screen_eligible = lambda sim: np.isnan(sim.people.date_screened) | (sim.t > (sim.people.date_screened + 5 / sim['dt']))
routine_screen = hpv.routine_screening(start_year=2016, product='hpv', prob=0.1, eligibility=screen_eligible, age_range=[30, 50], label='routine screening')
screened_pos = lambda sim: sim.get_intervention('routine screening').outcomes['positive'] # Get those who screen positive
pos_screen_assesser = hpv.routine_triage(start_year=2016, product = 'txvx_assigner', prob = 1.0, annual_prob=False, eligibility = screened_pos, label = 'txvx assigner') # Offer TxVx or refer them for further testing
txvx_eligible = lambda sim: sim.get_intervention('txvx assigner').outcomes['txvx'] # Get people who've been classified as txvx eligible based on the positive screen assessment
deliver_txvx = hpv.linked_txvx(prob = 0.8, product = 'txvx1', eligibility = txvx_eligible, label = 'txvx') # Deliver txvx to them

screen_vx_intv = [routine_screen, pos_screen_assesser, deliver_txvx]
for intv in screen_vx_intv: intv.do_plot=False

sim0 = hpv.Sim(pars=pars, label='Baseline')
sim1 = hpv.Sim(pars=pars, interventions=mass_vx_intvs, label='Mass therapeutic vaccination')
sim2 = hpv.Sim(pars=pars, interventions=screen_vx_intv, label='Therapeutic vaccination through screening')

# Run and plot
msim = hpv.parallel(sim0, sim1, sim2)
msim.plot();


#analyzers - objects that do not change the behavior of a simulation, but just report on its internal state
#results by age
import numpy as np
import sciris as sc #acts as a "library of the gaps"
import hpvsim as hpv

# Create some parameters, setting beta (per-contact transmission probability) higher
# to create more cancers for illutration
pars = dict(beta=0.5, n_agents=50e3, start=1970, n_years=50, dt=1., location='tanzania')

# Also set initial HPV prevalence to be high, again to generate more cancers
pars['init_hpv_prev'] = {
    'age_brackets'  : np.array([  12,   17,   24,   34,  44,   64,    80, 150]),
    'm'             : np.array([ 0.0, 0.75, 0.9, 0.45, 0.1, 0.05, 0.005, 0]),
    'f'             : np.array([ 0.0, 0.75, 0.9, 0.45, 0.1, 0.05, 0.005, 0]),
}

# Create the age analyzers.
az1 = hpv.age_results(
    result_args=sc.objdict(
        hpv_prevalence=sc.objdict( # The keys of this dictionary are any results you want by age, and can be any key of sim.results
            years=2019, # List the years that you want to generate results for
            edges=np.array([0., 15., 20., 25., 30., 40., 45., 50., 55., 65., 100.]),
        ),
        hpv_incidence=sc.objdict(
            years=2019,
            edges=np.array([0., 15., 20., 25., 30., 40., 45., 50., 55., 65., 100.]),
        ),
        cancer_incidence=sc.objdict(
            years=2019,
            edges=np.array([0.,20.,25.,30.,40.,45.,50.,55.,65.,100.]),
        ),
        cancer_mortality=sc.objdict(
            years=2019,
            edges=np.array([0., 20., 25., 30., 40., 45., 50., 55., 65., 100.]),
        )
    )
)

sim = hpv.Sim(pars, genotypes=[16, 18], analyzers=[az1])
sim.run()
a = sim.get_analyzer()
a.plot();
#also possible to plot these results alongside data
az2 = hpv.age_results(
    result_args=sc.objdict(
        cancers=sc.objdict(
            datafile='example_cancer_cases.csv',
        ),
    )
)
sim = hpv.Sim(pars, genotypes=[16, 18], analyzers=[az2])
sim.run()
a = sim.get_analyzer()
a.plot();


#snapshots - take “pictures” of the sim.people object at specified points in time
snap = hpv.snapshot(timepoints=['2020'])
sim = hpv.Sim(pars, analyzers=snap)
sim.run()

a = sim.get_analyzer()
people = a.snapshots[0]

# Plot age mixing
import pylab as pl
import matplotlib as mpl
fig, ax = pl.subplots(nrows=1, ncols=1, figsize=(5, 4))

fc = people.contacts['m']['age_f'] # Get the age of female contacts in marital partnership
mc = people.contacts['m']['age_m'] # Get the age of male contacts in marital partnership
h = ax.hist2d(fc, mc, bins=np.linspace(0, 75, 16), density=True, norm=mpl.colors.LogNorm())
ax.set_xlabel('Age of female partner')
ax.set_ylabel('Age of male partner')
fig.colorbar(h[3], ax=ax)
ax.set_title('Marital age mixing')
pl.show();


#age pyramids (like snapshots but in pyramid figures)
# Create some parameters
pars = dict(n_agents=50e3, start=2000, n_years=30, dt=0.5)

# Make the age pyramid analyzer
age_pyr = hpv.age_pyramid(
    timepoints=['2010', '2020'],
    datafile='south_africa_age_pyramid.csv',
    edges=np.linspace(0, 100, 21))

# Make the sim, run, get the analyzer, and plot
sim = hpv.Sim(pars, location='south africa', analyzers=age_pyr)
sim.run()
a = sim.get_analyzer()
fig = a.plot(percentages=True);


#calibration runs your simulation hundreds of times, each time trying different parameter values, and scores each run on how closely it matches your real-world data
#calibration object contains:
    #- an hpv.Sim() instance with details of the model configuration
    #- two lists of parameters to vary, one for parameters that vary by genotype and one for those that don’t
    #- dataframes that hold the calibration targets, which are typically added as csv files
    #- a list of any additional results to plot
    #- settings that are passed to the Optuna package (an open source hyperparameter optimization framework that automates calibration for HPVsim)
import hpvsim as hpv

# Configure a simulation with some parameters
pars = dict(n_agents=10e3, start=1980, end=2020, dt=0.25, location='nigeria')
sim = hpv.Sim(pars)

# Specify some parameters to adjust during calibration.
# The parameters in the calib_pars dictionary don't vary by genotype,
# whereas those in the genotype_pars dictionary do. Both kinds are
# given in the order [best, lower_bound, upper_bound].
calib_pars = dict(
        beta=[0.05, 0.010, 0.20],
    )

genotype_pars = dict(
    hpv16=dict(
        cin_fn=dict(k=[0.5, 0.2, 1.0]),
        dur_cin=dict(par1=[6, 4, 12])
    ),
    hpv18=dict(
        cin_fn=dict(k=[0.5, 0.2, 1.0]),
        dur_cin=dict(par1=[6, 4, 12])
    )
)

# List the datafiles that contain data that we wish to compare the model to:
datafiles=['nigeria_cancer_cases.csv',
           'nigeria_cancer_types.csv']

# List extra results that we don't have data on, but wish to include in the
# calibration object so we can plot them.
results_to_plot = ['cancer_incidence', 'asr_cancer_incidence']

# Create the calibration object, run it, and plot the results
calib = hpv.Calibration(
    sim,
    calib_pars=calib_pars,
    genotype_pars=genotype_pars,
    extra_sim_result_keys=results_to_plot,
    datafiles=datafiles,
    total_trials=3, n_workers=1
)
calib.calibrate(die=True)
calib.plot(res_to_plot=4);

