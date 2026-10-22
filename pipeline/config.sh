# epel-af2 에서 한 번만 고치는 곳. 각 스크립트 실행 전에 `source pipeline/config.sh`.
# 모든 값은 pipeline/lib/common.py 의 기본값을 덮어쓴다.

export CHCODH2_MNT=/mnt/af2results/mingyu
export CHCODH2_DOCK="$CHCODH2_MNT/workspace/docking_simulation"
export CHCODH2_SPLIT="$CHCODH2_DOCK/result/homolog_split"
export CHCODH2_OLDR="$CHCODH2_DOCK/result/homolog_selection"

export CHCODH2_IPR_FASTA="$CHCODH2_MNT/chcodh2/protein-matching-IPR010047.fasta"
export CHCODH2_MASTER="$CHCODH2_SPLIT/master_1827.tsv"
export CHCODH2_ALIGN="$CHCODH2_OLDR/chcodh2_alignment.tsv"

# 이번 작업의 산출물이 쌓이는 곳
export CHCODH2_OUT="$CHCODH2_SPLIT/cluster_sweep"

# 인계문서 6절 — Jupyter 의 PATH 를 믿지 않는다. 00_env_probe.py 가 찾아준 절대경로를 여기 박아라.
# export CHCODH2_MMSEQS=/절대/경로/mmseqs
