Allows control for fibril design pipeline through single config file.

Fetch pdb -> chain renaming based on 3D coordinates(since only takes one fibril), res_num renumbering, pyrosetta fast-relax -> mpnn formatting and config file -> submit

**Overall Proces**
'git clone (https://github.com/Chen-Jie7/fibril_pipeline.git)'

0. **If you don't have a conda installed, would recommend getting mamba(bascially same thing as conda)**
    - Follow instructions at (https://github.com/conda-forge/miniforge)

1. **Running ProteinMPNN**
    *Once time Setup Instructions*
    0. cd to fibril_pipline
    1. 'git submodule init'
    2. 'git submodule update'
    3. 'mamba env create -n environment.yml'
    4. 'mamba install pytorch torchvision torchaudio cudatoolkit=11.8 -c pytorch' or go to (https://pytorch.org/) for newest cuda version
    5. install pyrosetta
        a. pip install pyrosetta-installer==0.1.1
        b. python -c 'import pyrosetta_installer; pyrosetta_installer.install_pyrosetta(mirror=1)'

    *Verify run*:
    python 1_run_mpnn.py
    - This will run the current 7ymn sample based on mpnn_config file spcifications.
    - Aka pdb is in cleaned up, minimized, then launch the job mpnn for you
2. **Energetic Calculations**
    Config setting in Phase 3 to parallelize energetic calculations and saving those pdbs
    
3. **Running Upside2**
    1. git clone https://github.com/sosnicklab/upside2-md.git
    2. cd upside2-md
    3. ./install.sh
        - Here you might get a little upside install bug problem where it says can't find eigen3
            0. mamba activate upside2
            1. $CONDA_PREFIX -name "Macros.h" 2>/dev/null
            2. Copy paste the full path that ends with 'include/Eigen/src/Core/util/Macros.h'
            3. Go to upside2-md/cmake/FindEigen3.cmake
            4. Ctrl f to find 'file(READ "${EIGEN3_INCLUDE_DIR}/Eigen/src/Core/util/Macros.h" _eigen3_version_header)'
            5. Replace it with 'file(READ "[Your path above]" _eigen3_version_header)'
                Example: file(READ "/work/CAND/shared/Chens/miniforge3/envs/upside2/include/Eigen/src/Core/util/Macros.h" _eigen3_version_header)

            

            
