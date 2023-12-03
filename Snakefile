configfile: 'config.yaml'


rule all:
    pass

include: "rules/generate_simulations.smk"
include: "rules/run_methods.smk"