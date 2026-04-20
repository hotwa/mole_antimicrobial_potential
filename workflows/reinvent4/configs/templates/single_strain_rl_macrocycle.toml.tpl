run_type = "staged_learning"
device = ${device}
tb_logdir = ${tb_logdir}

[parameters]
prior_file = ${prior_file}
agent_file = ${agent_file}
smiles_file = ${smiles_file}
summary_csv_prefix = ${summary_csv_prefix}
batch_size = ${batch_size}
randomize_smiles = true
unique_sequences = true

[learning_strategy]
type = "dap"
sigma = ${sigma}
rate = ${learning_rate}

[diversity_filter]
type = "IdenticalMurckoScaffold"
bucket_size = ${bucket_size}
minscore = ${diversity_minscore}

[[stage]]
chkpt_file = ${chkpt_file}
termination = "simple"
max_score = ${max_score}
min_steps = ${min_steps}
max_steps = ${max_steps}

[stage.scoring]
type = "geometric_mean"

[[stage.scoring.component]]
[stage.scoring.component.custom_alerts]
[[stage.scoring.component.custom_alerts.endpoint]]
name = "Unwanted SMARTS"
params.smarts = [
    "[#8][#8]",
    "[#6;+]",
    "[#16][#16]"
]

[[stage.scoring.component]]
[stage.scoring.component.ExternalProcess]
[[stage.scoring.component.ExternalProcess.endpoint]]
name = ${score_component_name}
weight = ${score_component_weight}
params.executable = ${score_bridge_executable}
params.args = ${score_bridge_args}
params.property = ${score_property}
