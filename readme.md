This is the rerun of the entire process but building in the entire design pipeline for easness of use
**Overall Proces**
1. **ProteinMPNN**
    1. git submodule init
    2. git submodule update
    3. Edit the mpnn_config.yaml based on instructions
        - This config will be in charge from fetch the pdb and cleaning them to actually running mpnn.
    4. python run_mpnn.py
        - This will do everything then launch the job for you
