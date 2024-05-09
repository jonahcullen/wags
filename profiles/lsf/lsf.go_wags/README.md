## LSF profile

The profile was collected from: https://github.com/Snakemake-Profiles/lsf. Cookiecutter parameters may need to be modified to fit your HPC, read instructions to create a custom lsf profile and:


```
cookiecutter --output-dir /path/to/wags/profiles/lsf /path/to/lsf
```

### The cookiecutter parameters used to make this profile (most are default):

- LSF_UNIT_FOR_LIMITS: 3 (Gb)
- UNKNWN_behaviour: 1 (wait)
- ZOMBI_behaviour: 1 (ignore)
- latency_wait: 5
- use_conda: no
- use_singularity: no
- restart_times: 2
- print_shell_commmands: no 
- jobs: 300
- default_mem_mb: 1024
- default_cluster_logdir: logs/cluster
- default queue: None
- default_project: None
- max_status_checks_per_second: 10
- max_jobs_per_second: 10
- max_status_checks: 1
- wait_between_tries: 0.001
- jobscript_timeout: 10
- profile_name: lsf.go_wags
