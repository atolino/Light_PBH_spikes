# N-body simulations of gravitationally interacting PBHs in a fixed potential

This repository provides the configuration and initial conditions files to run N-body simulations in GADGET4 of light primordial black hole (PBH) spikes around a central solar mass PBH, represented by a fixed central potential.

We first present a quick-start guide that assumes GADGET4, conda, and screen are already installed. Further below, we provide more detailed instructions on how to install these tools and any additional packages required to run the simulations.


## THE LONG STORY SHORT (QUICK START)

Assuming GADGET4, conda, and screen are already installed and working:

1. Compile GADGET4 
   - Go to the GADGET4 source directory.
   - Copy the provided configuration file `Config.sh` into this directory.
   - In `Template-Make-system`, uncomment the compiler you want to use
     (e.g. GCC or Intel).
   - Rename the file:
         mv Template-Make-system Make-system
   - Compile:
         make

   This creates the `Gadget4` executable.

2. Prepare an example run
   - Create a directory inside `gadget4/examples/`, e.g. `myExample`.
   - Copy into it:
       - the compiled `Gadget4` executable
       - `run_pipeline.sh`
       - `param.txt`
       - `initial_conditions.py`
      - any analysis scripts (i.e. `analysis.py` and `units.py`)

3. Make the pipeline executable: <br>
         chmod +x run_pipeline.sh

4. Run the simulation <br>
         screen <br>
         conda activate gadget-env <br>
         ./run_pipeline.sh <run_name> <mass_ratio> <num_particles> <br>
   The simulation will keep running even if you disconnect with *Ctrl + A, then D*.
   Re-attach later with:
         screen -r

You can also give a name to the screen process with *screen -S name* and reattach it with *screen -r name*.

5. If you don't want to use screen and conda (ay!) <br>
   ./run_pipeline.sh <run_name> <mass_ratio> <num_particles> <br>


## PROJECT STRUCTURE


### GADGET INITIAL CONDITIONS AND ANALYSIS


Folder: gadget_fixed_potential

This directory contains all files needed to generate initial conditions,
run GADGET simulations, and analyze the results.

Main files:
- run_pipeline.sh <br>
    This script automates the full workflow for running a GADGET4 simulation of PBH systems.
    Usage: <br>
    ./run_pipeline.sh <run_name> <mass_ratio> <num_particles> <br>
   Requirements:<br>
   python3<br>
   mpirun (MPI)<br>
   Compiled Gadget4 executable in the working directory <br>
   Input files: <br>
   initial_conditions.py <br>
   param.txt <br>

- initial_conditions_ta.py <br>
  Generates and analyzes the initial conditions for a system with N−1 light
  PBHs and one heavy central PBH. The initial density is given at turnaround.

- param.txt <br>
  Parameter file used by GADGET.

- Config.sh <br>
  Configuration file used by GADGET.

- analysis/units.py <br>
  Utility files for units and snapshot loading.



## USEFUL TOOLS: CONDA, SCREEN, AND GADGET4

This project relies on three key tools that help managing the simulations:

1) screen   – session management
   `screen` is used to run long simulations safely on remote or shared
    machines. It allows the user to detach from a terminal session without
    stopping the running job, which is useful when you run time-consuming simulations.

2) conda    – Python environment management
    `conda` is used to manage a reproducible Python environment.
3) GADGET4  – simulation engine
   `GADGET4` is the N-body simulation code used to evolve the PBH system.
    It must be compiled after all system dependencies (MPI, HDF5, compiler)
    are available.



### TYPICAL WORKFLOW

The recommended workflow is:

1. Install SCREEN <br>

2. Install CONDA <br>

3. Compile GADGET4 <br>

On a remote machine, a typical session proceeds as follows:

1. Start a screen session: <br>
       screen

2. Activate the conda environment: <br>
       conda activate <env_name>

3. Compile or run GADGET4: <br>
       make
       ./run_pipeline.sh <run_name> <mass_ratio> <num_particles>

This workflow ensures that simulations continue running even if the
connection is lost, while maintaining a clean and reproducible
software environment.

### INSTALLING CONDA


Option A: Install Miniconda (recommended, no sudo required)

1. Download installer: <br>
   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

2. Run installer: <br>
   bash Miniconda3-latest-Linux-x86_64.sh

3. Follow prompts: <br>
   - Accept license
   - Install to default location (~/miniconda3)
   - Choose YES to initialize conda

4. Reload shell: <br>
   source ~/.bashrc

5. Verify installation: <br>
   conda --version

6. Create environment <br>
    
   conda create -n gadget-env -c conda-forge \
      python=3.11 \
      openmpi \
      compilers \
      gsl \
      fftw \
      hdf5 \
      zlib \
      openblas \
      cmake \
      make


7. Activate the environment: <br>

   conda activate gadget-env

It might be useful to install some packages:

    conda install numpy scipy matplotlib h5py

Optional (recommended for notebooks):

    conda install jupyter ipykernel

Option B: System-wide Conda (requires sudo)

   sudo apt install conda
   OR
   sudo yum install conda



### INSTALLING SCREEN (WITH SUDO)


Ubuntu / Debian:
   sudo apt update
   sudo apt install screen

RHEL / CentOS / Fedora:
   sudo dnf install screen
   (older systems: sudo yum install screen)

Verify:
   screen --version



### INSTALLING SCREEN (WITHOUT SUDO)


A) Check if screen is already available:
   screen
   screen --version

B) Local installation from source:

1. Download source: <br>
   wget https://ftp.gnu.org/gnu/screen/screen-4.9.1.tar.gz <br>
   tar xzf screen-4.9.1.tar.gz <br>
   cd screen-4.9.1 <br>

2. Configure (use old C standard for modern compilers): <br>
   ./configure --prefix=$HOME/.local CC="gcc -std=gnu89" <br>

3. Build and install: <br>
   make <br>
   make install <br>

4. Add screen to PATH: <br>
   export PATH="$HOME/.local/bin:$PATH" <br>
   (add this line to ~/.bashrc)

5. Verify: <br>
   ~/.local/bin/screen --version


C) Install terminfo locally (important, no sudo):

   mkdir -p ~/.terminfo
   tic -x -o ~/.terminfo terminfo/screencap

   export TERMINFO=$HOME/.terminfo
   export TERM=screen-256color

(Add the above exports to ~/.bashrc)



### BASIC SCREEN USAGE


Start a session:
   screen

Detach:
   Ctrl + A, then D

List sessions:
   screen -ls

Reattach:
   screen -r



## GADGET DEPENDENCY CHECKS
Activate your conda environment. Then check the following are installed:
### MPI CHECKS


Check MPI compiler:
   which mpicc
   mpicc --version

Check MPI runtime:
   which mpirun
   mpirun --version

Test MPI:
   mpirun -np 2 hostname

If this fails, MPI is not usable.

### HDF5 CHECKS


Check HDF5 compiler wrapper:
   which h5cc
   h5cc -showconfig

Check HDF5 libraries (optional, may not work without sudo):
   ldconfig -p | grep hdf5

Test compilation:
   h5cc test.c -o test_hdf5

If h5cc is missing, HDF5 is not available.

### MPI AND HDF5 INSTALLATION WITH CONDA

MPI and HDF5 can be installed directly inside a conda environment.
This is the recommended approach if you do not have sudo access or are
working on shared machines.

All commands below assume that a conda environment is already created
and activated:

    conda activate gadget-env


Install MPI using conda-forge:

    conda install -c conda-forge openmpi mpi4py

This provides:
- mpicc
- mpirun
- MPI runtime libraries
- Python MPI bindings (optional, via mpi4py)

Verify installation:

    which mpicc
    mpicc --version

    which mpirun
    mpirun --version

Test MPI:

    mpirun -np 2 hostname


Install HDF5 and its compiler wrapper:

    conda install -c conda-forge hdf5 h5py

This provides:
- HDF5 libraries
- h5cc compiler wrapper
- Python HDF5 bindings (via h5py)

Verify installation:

    which h5cc
    h5cc -showconfig

Test compilation:

    echo '#include <hdf5.h>' > test.c
    h5cc test.c -o test_hdf5

Important notes: 
- Conda-installed MPI and HDF5 are isolated inside the conda environment
  and do not require sudo privileges.
- When compiling GADGET4, make sure that the conda environment is
  activated so that `mpicc` and `h5cc` from conda are used.
- If your system provides MPI/HDF5 via environment modules, you should
  use either modules OR conda, but not both at the same time.

Also compilers should be installed.

Note that if you are still getting errors, try

      conda install -c conda-forge compilers mpi
      conda install -c conda-forge gsl
as the compilers  `x86_64-conda-linux-gnu-c++, gsl` might not be installed.
## GETTING GADGET4 SOURCE CODE


Clone the Gadget4 repository:

      git clone http://gitlab.mpcdf.mpg.de/vrs/gadget4

Enter the source directory:

      cd gadget4

You are now ready to configure and compile Gadget4.


## COMPILING AND RUNNING GADGET4



### COMPILING GADGET4


Inside the GADGET4 source directory, replace the default configuration
file with the provided one:

   cp /path/to/gadget_fixed_potential/Config.sh .

Compile the code:

    make



### SETTING UP A RUN EXAMPLE


Create a new example directory inside:

    gadget4/examples/

For instance:

    mkdir gadget4/examples/myExample

Copy the required analysis and run files into this directory:

   analysis.py
   units.py
   initial_conditions.py
   param.txt
   run_pipeline.sh



### RUNNING THE PIPELINE


Make the pipeline script executable:

    chmod +x run_pipeline.sh

Before running, open `run_pipeline.sh` and adjust the number of MPI
processes in the following line:

    nice -n 19 mpirun -np N_CORES --use-hwthread-cpus ./Gadget4 param.txt

where `N_CORES` should be set according to the available resources.

Notes:
- `nice -n 19` runs GADGET4 at low priority and is recommended when using
  shared machines.
- `--use-hwthread-cpus` allows the use of hardware threads (typically
  two threads per physical core).

You can then start the simulation by running:

    ./run_pipeline.sh <run_name> <mass_ratio> <num_particles>
