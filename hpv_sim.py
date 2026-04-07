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