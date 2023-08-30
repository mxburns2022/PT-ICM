from collections import defaultdict
from itertools import product
from simple_slurm import Slurm
from pathlib import Path
from datetime import datetime, timedelta
from utils import query_yes_no
import os

ROOT = Path(__file__).parent.parent # Escape from scripts folder
GRAPH_LOCATION = Path(ROOT, 'example_graphs')
get_all_files: callable = lambda d : [os.path.join(d, f) for f in os.listdir(d) if os.path.isfile(os.path.join(d, f))]

CONFIG = defaultdict(list)
CONFIG['-p'] = ['1.gset']
CONFIG['-T0'] = [15]
CONFIG['-T1'] = [2]
CONFIG['-pts'] = [300]
CONFIG['--inverse'] = ['true', 'false']
CONFIG['--preprocess'] = ['true', 'false']
CONFIG['--icm'] = ['true', 'false']
TOTAL_HOST = 8

def generate_commands() -> list[str]:
    cmds = []
    for config in product(*CONFIG.values()):
        cmd = ['mpirun', '-np', '8', './pt_icm']
        for k, v in zip(CONFIG.keys(), config):
            cmd.append(k)
            cmd.append(str(v))
        cmds.append(cmd)
    return cmds

def run() -> None:
    now: str = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    slurm = Slurm(
        cpus_per_task=8,
        job_name=f"PT_ICM_{now}",
        time=timedelta(days=1, hours=2, minutes=3, seconds=4),
        partition='ising',
    )
    join_cmd = lambda cmds : ' && '.join([' '.join(cmd) for cmd in cmds])

    # Split the commands into TOTAL_HOST groups and then join them together
    cmds = generate_commands()
    slurm_cmds_group = [[] for _ in range(TOTAL_HOST)]
    for i in range(len(cmds)):
        slurm_cmds_group[i % TOTAL_HOST].append(cmds[i])
    for i in range(TOTAL_HOST):
        slurm_cmds_group[i] = join_cmd(slurm_cmds_group[i])

    # Instructions
    print(f"You have {len(cmds)} commands to run. Divided into {TOTAL_HOST} groups.\n")
    print(f"The header of each slurm script is as follows:\n")
    print(slurm)
    
    # Submit the jobs
    if query_yes_no("Do you want to submit the jobs?"):
        print("Submitting jobs...")
        for cmd in slurm_cmds_group:
            print(cmd)
            # slurm.sbatch(cmd)
    else:
        print("Bye")

if __name__ == '__main__':
    run()