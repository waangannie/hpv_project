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