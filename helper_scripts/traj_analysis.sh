#!/bin/bash


files=()

for file in /work/CAND/shared/Chens/AD_fibrils/AD_redesign/7YMN_5_min/top_N_fibrils/*; do
    files+=("$(basename "$file")")
done

for file in "${files[@]}"; do 

    pdb_id=$file
    sim_id=Simulation1

    n_rep=10

    work_dir=/work/CAND/shared/Chens/AD_fibrils/AD_redesign/7YMN_5_min/Simulations
    input_dir=$work_dir/inputs
    output_dir=$work_dir/outputs
    run_dir=$output_dir/$sim_id

    if [ ! -d "$work_dir/results/" ]; then
        mkdir $work_dir/results/
    fi
    if [ ! -d "$work_dir/results/$sim_id" ]; then
        mkdir $work_dir/results/$sim_id
    fi
    for i in `seq 0 $((n_rep-1))` 
    
    do
        traj=$output_dir/$sim_id/$pdb_id.run.$i.up
        echo $traj
        python $UPSIDE_HOME/py/get_info_from_upside_traj.py  $traj $work_dir/results/$sim_id/${pdb_id}_${sim_id}_$i
        #python $UPSIDE_HOME/py/extract_vtf.py                $traj $work_dir/results/$sim_id/${pdb_id}_${sim_id}_$i.vtf
    done
done

