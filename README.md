# Quantum-Circuits-with-ITensor.jl

This repository comprises of multi-qubit circuits that were implemented and can be simulated through Julia's ITensors Library which uses matrix products states to represent qubits. The quantum circuits provided multi-qubit adder, variational ansantz circuit, and a time evolving block decimation (TEBD) circuit, that when run, will give you information such as the entanglement of a range of qubits, the measured Sz value of a qubit, or the sum of two binary numbers.

## Downloading Julia

### For Windows
Open command prompt and type the function:
```
winget install julia -s msstore
```
You will be prompted to agree to source agreement terms, when this occurs just type the letter 'y' as show below and Julia will begin downloading.

![Installation-of-Julia](Images/installing%20julia.png) <br><br>
Once download is complete, type 'julia' in the command prompt and that will start the julia terminal!<br><br>
![Installation-of-Julia](Images/complete%20installation%20julia.png) <br>

### For Mac/Linux 
Open terminal and type the function:
```
curl -fsSL https://install.julialang.org | sh
```
If prompted, agree to the terms presented and afterwards Julia will begin downloading. <br>

### Other operating systems or having trouble?

For more detailed instructions to install Julia click [here.](https://julialang.org/downloads/) <br>


## Installing Libraries
Before you can start running simualations you must download the some of the libraries. Simply open the Julia Terminal type the code below. For libraries such as Plots.jl and ITensors.jl that are not already downloaded when Julia is first installed, it will prompt you to confirm installation and if that occurs please follow the on-screen instructions. For example, windows to confirm on Windows, all you have to do is type the letter 'y' when prompted.
```
using ITensors
using Plots
using Statistics
using Printf
```

## Setting Up Simulation
Clone this Github repository, if you do not know click [here.](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository). Once that is complete, open your Julia terminal and go to the directory where this repository is located by using the cd command. A sample is shown below (path represents the file path to get to the repository's folder). An alternative option is to change the directory before you start julia terminal.
```
cd("path/Quantum-Circuits-with-ITensor.jl")
```
For example:<br><br>
![Going-to-directory](Images/going%20into%20repository%20julia.png) <br><br>

## Running Initial Simulation
Now we are ready to start simulating circuits! Go into the Simulatin Circuits folder and copy the code in the SimulationFile.jl folder and paste it into the julia terminal and press enter. The simulation will start and there will printed output to show the progress of the simulation. Once the simulation is complete, the outputs for each circuit will be in the IO Files/Output Files folder with the corresponding file names. Additionally, scatterplots for the Sz measurements v.s. time for the TEBD and variational ansantz circuit are located in IO Files/Scatterplots folder. <br>

Sample of when Simulation Code is ran:<br>
![Sample-Of-Running-Simulation-Code](Images/Complete%20Run%20of%20Simulation%20Code.png) <br><br>

To rerun the simulation, you must copy and paste the code from SimulationFile.jl and paste it into the Julia Terminal. We must do this since there is no way to run julia files in the terminal unlike java and c files. On the bright side, when you rerun the code is it much faster.

## Changing Input Files for Each Circuit in Simulation
Input for each circuit in the simulation can be changed in the IO FIles/Input files folder. Below, will describe the format for the input and what each input represents. Note, the larger you make the input the longer the simulation is going to take as it takes more resources to setup and run the circuit

- Multi-Qubit Adder Circuit(IO Files/Input Files/bitAdderInput.txt):
    - First Line: Binary Integer of any length
    - Second Line: Binary Integer of any length
    - Third Line: Site to measure the entanglement relative to the first site (From site 2 to (4*length of longer bit))
- Variational Ansantz Circuit(IO Files/Input Files/variationalAnsantzInput.txt):
    - First Line: Number of sites
    - Second Line: Duration of each time step
    - Third Line: Total number of time(Note the total number of time divided by the duration for each time step is the amount of time steps that the circuit will undergo)
    - Fourth Line: Site to measure the entanglement relative to the first site (From site 2 to the total amount of sites)
- TEBD Circuit(IO Files/Input Files/tebdInput.txt):
    - First Line: Number of sites
    - Second Line: Duration of each time step
    - Third Line: Total number of time(Note the total number of time divided by the duration for each time step is the amount of time steps that the circuit will undergo)
    - Fourth Line: Site to measure the entanglement relative to the first site (From site 2 to the total amount of sites)
