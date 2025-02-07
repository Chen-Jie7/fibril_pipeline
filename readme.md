Allows control for fibril design pipeline through single config file.

Fetch pdb -> chain renaming based on 3D coordinates(since only takes one fibril), res_num renumbering, pyrosetta fast-relax -> mpnn formatting and config file -> submit

**Overall Proces**
'git clone (https://github.com/Chen-Jie7/fibril_pipeline.git)'

0. **If you don't have a conda installed, would recommend getting mamba(bascially same thing as conda)**
    - Follow instructions at (https://github.com/conda-forge/miniforge)
    
1. **Runnig ProteinMPNN**
    *Once time Setup Instructions*
    0. cd to fibril_pipline
    1. git submodule init
    2. git submodule update
    3. mamba create env -n environment.yml
    4. mamba install pytorch torchvision torchaudio cudatoolkit=11.8 -c pytorch or go to (https://pytorch.org/) for newest cuda version
    5. install pyrosetta
        a. pip install pyrosetta-installer=0.1.1
        b. python -c 'import pyrosetta_installer; pyrosetta_installer.install_pyrosetta(mirror=1)'

    *Verify run*:
    python run_mpnn.py
    - This will run the current 7ymn sample based on mpnn_config file spcifications.
    - Aka pdb is in cleaned up, minimized, then launch the job mpnn for you
