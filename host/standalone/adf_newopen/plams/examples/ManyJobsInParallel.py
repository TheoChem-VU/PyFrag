import multiprocessing

# Run as many job simultaneously as there are cpu on the system:
maxjobs = multiprocessing.cpu_count()
print('Running up to {} jobs in parallel simultaneously'.format(maxjobs))

# Set the default jobrunner to be parallel:
config.default_jobrunner = JobRunner(parallel=True, maxjobs=maxjobs)

# Number of cores for each job:
config.job.runscript.nproc = 1

# Load all molecules in the folder 'molecules':
molecules = read_molecules('molecules')

settings = Settings()
settings.input.ams.Task = 'GeometryOptimization'
settings.input.dftb.Model = 'GFN1-xTB'

results = []
for name, molecule in sorted(molecules.items()):
    job = AMSJob(molecule=molecule, settings=settings, name=name)
    results.append(job.run())

# Only print the results of the succesful caluclations:
for result in [r for r in results if r.ok()]:
    print('Energy for {:<12}: {:>10.3f} kcal/mol'.format(result.name, result.get_energy(unit='kcal/mol')))
