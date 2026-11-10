# Slurm templates

This directory contains slurm job scripts for the heavier steps of the tutorial. With the small simulated data set you can run everything directly in the terminal, but with real data you should submit steps like these as jobs, so they run on the compute nodes of the cluster instead of the login node.

## The scripts

Each script belongs to a section of the main README. Submit it **from that section's directory** (the scripts use the same relative paths as the README, e.g. `../data/...`).

| Script | Submit from | README section | Type |
|---|---|---|---|
| `1.1_plink_pca.sh` | `1.1_pca_with_plink` | Using Plink | single job |
| `1.2_cram_to_bam.sh` | `1.2_pca_with_angsd` | Using angsd (CRAM to BAM) | array, one task per CRAM |
| `1.2_angsd_pcangsd.sh` | `1.2_pca_with_angsd` | Using angsd (ANGSD + PCAngsd) | single job |
| `2.1_admixture_array.sh` | `2.1_admixture` | Run ADMIXTURE for several values of K | array, one task per K |
| `2.2_angsd_ngsadmix_input.sh` | `2.2_ngsadmix` | With NGSadmix (likelihoods + thinning) | single job |
| `2.2_ngsadmix_array.sh` | `2.2_ngsadmix` | With NGSadmix (run NGSadmix) | array, one task per K |
| `3_dofasta_and_mask_array.sh` | `3_phylogeny` | Make masked sequences with ANGSD | array, one task per sample |
| `3_extract_filter_windows.sh` | `3_phylogeny` | Make masked sequences with ANGSD (windows) | single job |
| `3_iqtree_array.sh` | `3_phylogeny` | Estimate a tree for each window | array, one task per window |
| `3_collapse_astral.sh` | `3_phylogeny` | Combine the window trees with ASTRAL-IV | single job |

Quick steps such as plotting with R, making population files with `cut`/`awk` and `bedtools makewindows` take seconds, so just run those directly.

Every script has a short usage note at the top that shows the exact command to submit it.

## Before you submit

1. **Fill in the placeholders.** At the top of every script you'll find:

        #SBATCH --account=<SLURM_ACCOUNT>
        #SBATCH --partition=<PARTITION>

    Replace `<SLURM_ACCOUNT>` with your compute project (e.g. `naiss2026-x-xxx`) and `<PARTITION>` with a partition on your cluster (e.g. `shared` on Dardel). You can also leave the script untouched and give these on the command line instead, which overrides the values in the script:

        sbatch --account=naiss2026-x-xxx --partition=shared ../slurmscripts/1.1_plink_pca.sh

2. **Create the `logs` directory.** The scripts write their output and error messages to `logs/`. Slurm does not create this directory for you, and if it's missing the job fails without writing any log. So always run this first, in the directory you submit from:

        mkdir -p logs

3. **Make sure the software is found.** If you use the container, `export POPGEN_SIF=...` (and `APPTAINER_BIND`) first; if you use conda, activate the environment. See the next section.

## Software: container or conda

The scripts can get their software from the tutorial container (recommended) or from a conda environment. See the "Software" section of the main README for how to get either of them. A slurm job runs on a different computer (a compute node) than the one you're typing on, but by default it inherits the environment variables of the terminal you submitted from, which is what both options below rely on.

### With the container (recommended)

Tell the scripts where the container file is, by setting `POPGEN_SIF` in your terminal before you submit:

```bash
export POPGEN_SIF=/path/to/popgen_tutorial.sif
export APPTAINER_BIND=/path/to/popgen_tutorial   # the cloned repo (which contains data/), see below
sbatch ../slurmscripts/2.1_admixture_array.sh
```

You need to export these once per terminal session (or put the two `export` lines in your `~/.bashrc`). Instead of exporting `POPGEN_SIF`, you can also write the path directly in the script, on the line `sif=${POPGEN_SIF:-}`.

Inside the scripts you'll see `$run` in front of every program, e.g.:

```bash
$run angsd -b bamlist.txt ...
```

When a container is given, `$run` means `apptainer exec popgen_tutorial.sif`, which runs that program inside the container. Without a container, `$run` is empty and the line is just `angsd -b bamlist.txt ...`. In two scripts (`3_extract_filter_windows.sh` and `3_collapse_astral.sh`) a whole loop is run inside the container at once with `$run bash <<'EOF' ... EOF`, since starting the container separately for thousands of small commands would be slow.

**About `APPTAINER_BIND`:** a container can only see the folders that are *bound* (made visible) inside it. By default these are your home directory, `/tmp` and the directory you're in. The scripts are run from a section directory, e.g. `2.1_admixture`, but read files from `../data`, which is outside it. If your tutorial folder is in your home directory this works anyway, but if it's in project storage (e.g. `/cfs/klemming/projects/...` on Dardel), you need to bind it, otherwise you'll get "No such file or directory" errors for files that clearly exist.

### With conda

Leave `POPGEN_SIF` unset. **If `popgen_tutorial_env` is activated when you run `sbatch`, the job will usually find all programs.** If that doesn't work on your cluster, or you prefer the job to set things up itself, each script has a commented-out block you can enable by removing the `#` at the start of the lines:

```bash
# ml PDC miniconda3                    # Dardel only: load conda as a module
# source activate popgen_tutorial_env  # or: conda activate popgen_tutorial_env
```

### The software check

Before doing any work, every script checks that the programs it needs can be found (in the container, or in your environment), and stops right away with a message like this if one is missing:

    ERROR: 'plink2' was not found. Check the software settings at the top of this script.

This way you find out after a few seconds rather than after waiting in the queue and running half of the job. If you use the container and see this, look in the `.err` log for a message from Apptainer just above it, e.g. that the `.sif` file can't be found. The scripts never install software themselves.

## Single jobs and array jobs

A **single job** runs one script once, e.g. `1.2_angsd_pcangsd.sh`.

An **array job** runs the same script many times in parallel. Each copy (a *task*) gets its own number in the variable `$SLURM_ARRAY_TASK_ID`, and the script uses that number to decide what to work on:

- In `2.1_admixture_array.sh` and `2.2_ngsadmix_array.sh`, the line `#SBATCH --array=2-6` starts five tasks numbered 2-6, and the task number is used directly as K.
- In the other array scripts, the task number picks one file from a list, e.g. task 3 converts the 3rd CRAM file:

        cram=$(ls "$bamdir"/*.dedup.cram | sed -n "${SLURM_ARRAY_TASK_ID}p")

    Here you need to tell slurm how many tasks to start, which is the number of files, when you submit:

        n=$(ls ../data/map_and_call_out_standard/02_bamfiles/*.dedup.cram | wc -l)
        sbatch --array=1-"$n" ../slurmscripts/1.2_cram_to_bam.sh

If you have many tasks (e.g. thousands of IQ-TREE windows), add `%` and a number to limit how many run at the same time, e.g. `--array=1-"$n"%50`. Clusters also have a maximum array size. If you get an error about this, ask your support team or split the array into several submissions.

## Keeping track of jobs

```bash
squeue -u $USER              # your queued and running jobs
sacct -j <jobid>             # status of a job (also after it has finished)
scancel <jobid>              # cancel a job (or a whole array)
```

The output of each job ends up in `logs/`:

- single jobs: `logs/<jobname>-<jobid>.out` and `.err`
- array jobs: `logs/<jobname>-<arrayjobid>_<tasknumber>.out` and `.err`

If a job fails, look in the `.err` file first. Each step needs the output of the previous one, so wait until a job has finished successfully (state `COMPLETED` in `sacct`) before submitting the next one.

## Adjusting resources

The `--time` and `--cpus-per-task` values are set for the tutorial data. For real data, start with a test on a few samples or windows, check how long it took with `sacct -j <jobid> --format=JobID,Elapsed,MaxRSS`, and adjust. If a job runs out of memory, request more CPUs (on many clusters memory scales with the number of CPUs) or add an `#SBATCH --mem=...` line if your cluster allows it.
