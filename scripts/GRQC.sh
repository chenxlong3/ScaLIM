cd ../build

EPS=0.1
NUM_CAND_EDGES=0
K_EDGES=50
RAND_SEED=2024
K_SEED=50

# ./format_graph -dataset ../data/GRQC -filename edgelist_ic.txt

# SEED_MODE="RAND"
SEED_MODE="OUTDEG"
NUM_SAMPLES=0
K_EDGES=100
./ScaLIM -dataset ../data/GRQC -epsilon $EPS -delta 0.001 -k_seed $K_SEED -k_edges $K_EDGES -rand_seed $RAND_SEED -num_cand_edges $NUM_CAND_EDGES -fast True -probability WC -seed_mode $SEED_MODE -method ScaLIM_minus -num_samples $NUM_SAMPLES
./ScaLIM -dataset ../data/GRQC -epsilon $EPS -delta 0.001 -k_seed $K_SEED -k_edges $K_EDGES -rand_seed $RAND_SEED -num_cand_edges $NUM_CAND_EDGES -fast True -probability WC -seed_mode $SEED_MODE -method ScaLIM -num_samples $NUM_SAMPLES
./ScaLIM -dataset ../data/GRQC -epsilon $EPS -delta 0.001 -k_seed $K_SEED -k_edges $K_EDGES -rand_seed $RAND_SEED -num_cand_edges $NUM_CAND_EDGES -fast True -probability WC -seed_mode $SEED_MODE -method Greedy2 -num_samples $NUM_SAMPLES
