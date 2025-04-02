
import os
import numpy as np
import pandas as pd
import json
import notion_df as ndf
from typing import Iterable
from tqdm import tqdm
from itertools import product
from datetime import datetime
from shutil import copy, rmtree
from subprocess import run, PIPE


ndf.pandas()


def to_iter(x):
    return x if isinstance(x, Iterable) and not isinstance(x, str) else [x]


def to_val(value: str):
    if value.lower() == 'false':
        return False
    if value.lower() == 'true':
        return True
    try:
        return int(value)
    except ValueError:
        try:
            return float(value)
        except ValueError:
            return value


def list_sum(*arglist):
    start = []
    for l in arglist:
        start += l
    return start


"""
Slurm Job Array creation script (Please pardon our dust, work in progress)
"""


# Parameters for your sweep, you only need to enter values that
#   differ from the base
coeff = int(round(2 * 6.2e-6 / (31e4 * 20e-15 * 1e-4 )))
config_list = [
   
    {
       "pts": ["H1"],
        'temp0': [0.3],
        'temp1': [2.0],
        'preprocess': [True]
    }
]


# Specify the graph class and graph names that you want to run,
#   and how many iterations of each
_iter = np.arange(20)
nvals = np.arange(260, 501, 10)
sweep_job_list = [     
        {
    'graph_class': ['set'],
    'graph': list_sum(  
      [f'G0{i}' for i in range(22, 30)], 
      [f'G0{i}' for i in range(35, 42)], 
    ),
    'iter': range(100)
    },
    
]
# Define the simulation name and short justification
# (folder will be sim{sim_no:03d}_{name})
project = 'PT+ICM'
sim_no = 0
name = "pt_sat_qubo"
simname = f'sim{sim_no:03d}_{name}'
justif = """
    Testing the efficicacy of the Nusslein formulation
"""
timeval = datetime.now().isoformat(timespec='seconds')



# Simulation system parameters
time = '5-00:00:00'  # runtime in days-hours:minutes:seconds
partition = 'ising'
memory = 8  # req memory in GB
concurrent = 200  # max number of jobs running concurrently
makejobs = 8

# Get path information from system environment variables
#  (should only need to configure once)
# SET THIS TO YOUR DATA PATH

thispath = os.path.dirname(os.path.abspath(__file__))
queue_path = f'{thispath}/queue_host.py'
directory_path = '/scratch/mburns13/data/ISING_RESULTS'
# automatically get sampler version e.g. BRIM-fast
solver_name = os.path.abspath(__file__).split(os.sep)[-2]
bluehive = os.path.exists('/scratch/mhuang_lab')
basepath = './base_config.json'  # Change as desired
base_config = json.loads(open(basepath).read())
base_config = dict([(u, to_val(v))
                   for u, v in base_config.items()])
cmdlineparams = set(list(base_config.keys()))
# Replace with your path to GSET folder (or keep GSET)
graph_dir = os.environ['GSET']
# Replace with your path to random inits (or don't)
random_dir = os.environ['RND']

config_list = [dict([(u, to_iter(v))
                     for u, v in i.items()]) for i in config_list]
sweep_job_list = [dict([(u, to_iter(v))
                   for u, v in sweep_jobs.items()]) for sweep_jobs in sweep_job_list]
# __________SHOULD NOT NEED TO RECONFIGURE ANYTHING BELOW THIS LINE__________
config_list = [{**i, **sw} for i, sw in product(config_list, sweep_job_list)]
njobs = 0
for i in config_list:
    njobs += np.prod([len(j) for j in i.values()])
print(njobs)
if njobs > 20000:
    print(f'Estimated {njobs} jobs. Continue? (y/n)')
    token = input()
    while (token not in {"y", "n"}):
            print('Please enter either y or n')
            token = input().lower()
    if token == 'n':
        quit()

colnames = list_sum(*[list(
    list_sum(*[[i] if isinstance(i, str) else list(i) for i in i.keys()])
) for i in config_list])
colname = []
paramset = list(set(colnames))


def main():
    print(random_dir)
    result_dir = f'{directory_path}/{solver_name}/{simname}'
    if os.path.exists(result_dir):
        print('Directory already exists, ' +
              'do you want to rebuild (b), replace job files (r), or exit (q)')
        token = input().lower()
        while (token not in {"r", "b", "q"}):
            print('Please enter either r, b or q')
            token = input().lower()
        if token == 'b':
            print('Removing previous directory...')
            rmtree(result_dir)
        elif token == 'q':
            print('Exiting...')
            exit(0)
    
    if not os.path.exists(result_dir + '/build'):
        print('Making/building directory...')
        os.makedirs(result_dir)
        os.mkdir(result_dir+'/build')
        copy_and_compile(result_dir+'/build')
    parampath = result_dir+'/parameters'
    if not os.path.exists(parampath):
        os.mkdir(parampath)
    # send a copy of this script and the base config to the parameter subdir
    copy(__file__, parampath)
    copy(basepath, parampath)
    jobpath = result_dir+'/jobs'
    if not os.path.exists(jobpath):
        os.mkdir(jobpath)
    jobrecords = []
    execpath = result_dir+'/build/pt_icm'
    
    print('Enumerating jobs and writing launch files...')
    for ind, c in tqdm(enumerate(configurations(config_list))):
        jobno = ind + 1
        jobrecords.append(tuple([jobno] + [c[i] for i in paramset]))
        make_job_file(execpath=execpath,
                      job_no=jobno,
                      config=c,
                      target_directory=result_dir)
    numjobs = len(jobrecords)
    record_df = pd.DataFrame.from_records(jobrecords, 
                                          columns=['job']+paramset)
    record_df.to_csv(parampath+'/records.csv', index=False)
    print(f'Wrote record of {numjobs} jobs to {parampath+"/records.csv"}')
    print('Writing the sbatch launch script...')
    if not os.path.exists(result_dir+'/data'):
        os.mkdir(result_dir+'/data')
    if not os.path.exists(result_dir+'/err'):
        os.mkdir(result_dir+'/err')
    make_launch_file(result_dir)
    print(f'To run, enter: sbatch {result_dir}/parameters/launch_queue.sh')

    
    newrow = pd.DataFrame(data=[('In Progress',
                                 name,
                                 timeval,
                                 solver_name,
                                 project,
                                 justif,
                                 result_dir,
                                 parampath)],
                      columns=[
                          'Status',
                          'Name',
                          'Date',
                          'Solver',
                          'Project',
                          'Justification',
                          'Result Path',
                          'File Path'
                      ])
    
    
# define configurations as values within an enumerated variable space
def configurations(config_list: list):
    for parameters in config_list:
        param_names = list(parameters.keys())
        param_space = product(*list(parameters.values()))
        for config in param_space:
            newdict = {}
            for name, val in zip(param_names, config):
                if not isinstance(name, str):
                    for name1, val1 in zip(name, val):
                        newdict[name1] = val1
                else:
                    newdict[name] = val
            if not os.path.exists(
                f'{graph_dir}/{newdict["graph_class"]}/{newdict["graph"]}'
            ):
                print(os.path.join(
                    graph_dir,
                    newdict["graph_class"],
                    newdict["graph"]
                ))
                continue
            yield {**base_config, **newdict}


def copy_and_compile(target_directory: str):
    """Copy the source files to `target_directory` and compile the binary

    Args:
        target_directory (str): Directory containing the compiled/source 
            code for the experiment
    """
    if bluehive:
        os.system('module unload gcc')
        os.system('module load gcc/11.2.0')
        os.system('module load cmake')
    output = run(['cmake', f'-B{target_directory}', f'-S{thispath}'], stderr=PIPE)
    try:
        output.check_returncode()
    except:
        print("CMake error with code:\n{}".format(output.stderr.decode()))
        exit(0)
    output = run(['make', f'--directory={target_directory}', f'-j{makejobs}'], stderr=PIPE)
    try:
        output.check_returncode()
    except:
        print("Make error with code:\n{}".format(output.stderr.decode()))
        exit(0)

def make_launch_file(target_directory: str):
    script_contents = f"""#!/bin/bash
#SBATCH --job-name QUEUE_{simname}
#SBATCH --output={target_directory}/parameters/queue.run
#SBATCH --error={target_directory}/err/queue.err
#SBATCH -p {partition}
#SBATCH --mem=1G
#SBATCH --ntasks=1
#SBATCH --time={time}

python3 {queue_path} {simname}
"""
    with open(target_directory+'/parameters/launch_queue.sh', mode='w') as launch:
        launch.write(script_contents)
        
def make_job_file(execpath: str, 
                  target_directory: str, 
                  config: dict, 
                  job_no: int):
    """Generate launch files for each PT experiment

    Args:
        target_directory (str): Path to the directory containing the job    
            scripts
        config (dict)
    """
    header = f"""#!/bin/bash
#SBATCH --job-name {job_no}_pt
#SBATCH --output={target_directory}/data/pt_output_{job_no}.run
#SBATCH --error={target_directory}/err/pt_output_{job_no}.err
#SBATCH -p {partition}
#SBATCH --mem={memory}G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time={time}
    """
    filepath = f'{target_directory}/jobs/pt_job_{job_no}.sh'
    with open(filepath, mode='w') as jobfile:
        
        jobfile.write(f'{header}\n')
        jobfile.write(f'{execpath}\\\n')
        for parameter, value in config.items():
            if parameter not in cmdlineparams:
                continue
            if isinstance(value, bool):
                if value:
                    jobfile.write(f'\t--{parameter}\\\n')
            else:
                jobfile.write(f'\t--{parameter} {value}\\\n')
        jobfile.write(
            f'\t--problem {graph_dir}/' +
            f'{config["graph_class"]}/{config["graph"]}\\\n')
            
        jobfile.write(f'\t--seed {config["iter"]}\\\n')


if __name__ == '__main__':
    main()
