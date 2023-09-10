from collections import defaultdict
from itertools import product
from pathlib import Path
from datetime import datetime, timedelta
from utils import query_yes_no
import os

ROOT = Path(__file__).parent.parent  # Escape from scripts folder
GRAPH_LOCATION = Path(ROOT, 'example_graphs')
get_all_files: callable = lambda d: [os.path.join(d, f) for f in os.listdir(d) if os.path.isfile(os.path.join(d, f))]

CONFIG = defaultdict(list)
CONFIG['-p'] = list(map(lambda x: os.path.join(GRAPH_LOCATION, x), ['G022', 'G039']))
CONFIG['-T0'] = [2, 3, 4]
CONFIG['-T1'] = [0.25, 0.5, 0.75]
CONFIG['--pts'] = [20000, 30000, 40000]
CONFIG['--inverse'] = ['false']
CONFIG['--preprocess'] = ['false']
CONFIG['--icm'] = ['true', 'false']
CONFIG['--seed'] = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
TOTAL_HOST = 64
ISODATE = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')


def generate_commands() -> list[str]:
    cmds = []
    for config in product(*CONFIG.values()):
        cmd = ['mpirun', '-np', '1', './pt_icm']
        for k, v in zip(CONFIG.keys(), config):
            if v == 'true':
                cmd.append(k)
            elif v == 'false':
                continue
            else:
                cmd.append(k)
                cmd.append(str(v))
        cmds.append(cmd)
    return cmds


def run():

    header = f"""#!/bin/bash
#SBATCH -p ising
#SBATCH -J PT_{ISODATE}
#SBATCH -e PT_{ISODATE}.err.%j
#SBATCH -o PT_{ISODATE}.run.%j
#SBATCH --mem-per-cpu=8GB
#SBATCH -t 48:00:00
#SBATCH --cpus-per-task=8

module load mpich;
module unload gcc;
module load gcc/9.1.0;
"""
    cmds = generate_commands()
    slurm_cmds_group = [[] for _ in range(TOTAL_HOST)]
    for i in range(len(cmds)):
        slurm_cmds_group[i % TOTAL_HOST].append(cmds[i])
    script_cnt = 0
    for i in range(TOTAL_HOST):
        filename = f"slurm_{script_cnt}.sh"
        with open(os.path.join(Path(__file__).parent, filename), 'w') as f:
            f.write(header)
            for cmd in slurm_cmds_group[i]:
                f.write(' '.join(cmd) + '\n')
        script_cnt += 1
    with open('final.sh', 'w') as f:
        for i in range(script_cnt):
            f.write(f'sbatch slurm_{i}.sh\n')
    print(f"Total {script_cnt} scripts generated.")
    # cnt = 0
    # sh_cnt = 0
    # cmd = ""
    # UPPER_BOUND = 150
    # for i in range(UPPER_BOUND):
    #     graph_index = int(i)
    #     graph_name = f"H{graph_index:03d}"
    #     cmd += f"matlab -nodisplay -nodesktop -r 'clear; cd /home/lsun16/Desktop/OPO; load tiny_all.mat; OPO_XY_nondim_BO_ep({graph_name}); exit' | tee -a {graph_name}.txt"
    #     cmd += "\n"
    #     if cnt == 2 or i == UPPER_BOUND - 1:
    #         with open(f"slurm_{sh_cnt}.sh", 'w') as f:
    #             f.write(header)
    #             f.write(cmd)
    #         cmd = ""
    #         cnt = 0
    #         sh_cnt += 1
    #     cnt += 1

    # with open("final.sh", 'w') as f:
    #     for i in range(sh_cnt):
    #         f.write(f"sbatch slurm_{i}.sh\n")

# def run() -> None:
#     now: str = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
#     slurm = Slurm(
#         cpus_per_task=8,
#         job_name=f"PT_ICM_{now}",
#         time=timedelta(days=1, hours=2, minutes=3, seconds=4),
#         partition='ising',
#     )

#     join_cmd = lambda cmds : ' && '.join([' '.join(cmd) for cmd in cmds])

#     # Split the commands into TOTAL_HOST groups and then join them together
#     cmds = generate_commands()
#     slurm_cmds_group = [[] for _ in range(TOTAL_HOST)]
#     for i in range(len(cmds)):
#         slurm_cmds_group[i % TOTAL_HOST].append(cmds[i])
#     for i in range(TOTAL_HOST):
#         slurm_cmds_group[i] = join_cmd(slurm_cmds_group[i])

#     # Instructions
#     print(f"You have {len(cmds)} commands to run. Divided into {TOTAL_HOST} groups.\n")
#     print(f"The header of each slurm script is as follows:\n")
#     print(slurm)

#     # Submit the jobs
#     if query_yes_no("Do you want to submit the jobs?"):
#         print("Submitting jobs...")
#         for cmd in slurm_cmds_group:
#             print(cmd)
#             # slurm.sbatch(cmd)
#     else:
#         print("Bye")


if __name__ == '__main__':
    run()
