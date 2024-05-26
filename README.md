# Welcome to NcorpiON

NcorpiON is an N-body software dedicated to the simulation of any gravitational system, and in particular collisional and fragmenting systems. It is very fast and can realistically handle fragmentations.
NcorpiON has been specifically designed to be both fast and precise in systems where the number $N$ of bodies is large.

## Features

- Written in C
- Four built-in modules for mutual interactions, namely
- Brute-force $\mathcal{O}(N^2)$ method
- Barnes-Hut $\mathcal{O}(N \ln N)$ tree code
- FalcON $\mathcal{O}(N)$ Fast Multipole Method
- Mesh $\mathcal{O}(N)$ algorithm
- A built-in fragmentation model that can realistically handle violent collisions
- A python add-on to produce 2D animations of the simulations
- 3D visualization of the simulations in real-time (using REBOUND)
- Completely open source. All features are available in the public github repository.
- Runs on Linux and MacOS.


## Installation

The installation is straightforward. You can install NcorpiON with the ```git``` command

	git clone https://github.com/Jeremycouturier/NcorpiON.git
	
or if you prefer, you can simply download this repository to your computer. Using ```git``` is recommended so you can later on get updates by running ```git pull```


## Setting up a simulation

In order to set up your simulation, all you have to do is to update the file ```src/parameters.h``` with your preferences. This file is used to define the behavior of NcorpiON. From there, you can decide

- Which module you want to use for mutual interactions (collisions and mutual gravity)
- How collisions should be resolved
- What physical effects to take into account
- What are the initial conditions (random or read from file)
- The input/output path
- Whether or not NcorpiON should : output the simulation to files, produces 2D animations, or launch a real-time 3D visualization in a browser tab
- The system of units of your simulation
- And much more ...

Every item of the parameter file is commented and explained within the file, so you might already have a pretty decent idea of the capabilities of NcorpiON after having a look at it. A more
complete documentation is also available on the NcorpiON official website at <https://ncorpion.com/#setup>.


## Running a simulation

Once the parameters file is updated and your simulation set up, you can compile NcorpiON and run the simulation from the ```src``` directory all with the unique command

	./ncorpion.sh
	

## Exploiting a simulation

A numerical simulation is only useful if you have ways to exploit the data it computes. NcorpiON gives you three ways to do so

- By outputing to files the state of the simulation (coordinates, masses and radii of the bodies) at a frequency you determined
- By producing, at the end of the simulation, a 2D animation (gif and mp4).
- By launching a 3D visualization in your web browser while the simulation is running.

These options are compatible together and you can choose one or several of them in ```src/parameters.h```.


## 3D real-time visualization

The 3D visualization of simulations in real time is made possible by the REBOUND software. REBOUND uses a web browser for the visualization in order to keep to a minimum dependancies on external libraries.
In order to make real-time 3D visualization possible, NcorpiON has to communicate with REBOUND. At every timestep of the simulation, NcorpiON computes the new state of the system and sends that information to REBOUND.
REBOUND then communicates with your web browser for the 3D rendering.

The communication between NcorpiON and REBOUND uses MPI. If you are a Linux user on a Debian-based distro, you can download MPI on your system with the command

      sudo apt install mpich
      
If you are a MacOS user, you can install MPI with brew

      brew install openmpi

Note that the parts of REBOUND needed for the 3D visualization are already embedded in NcorpiON, so you don't need to download REBOUND to enjoy the 3D visualization (although you are still encouraged to check out the REBOUND project at <https://rebound.readthedocs.io>).
Depending on your needs, REBOUND could be more adapted to your N-body problem than NcorpiON. Once MPI is installed, NcorpiON takes care of everything and you only need to run the command

      ./ncorpion.sh
      
from within the ```src``` directory. Note that the 3D visualization will make the simulation vastly slower.


## 2D animations

The production of images for the 2D animations uses python (with libraries numpy and matplotlib). The assembly of these images into animations requires ```ffmpeg``` to be installed on your system. Run

      sudo apt install ffmpeg
	
on Debian based distros or

      brew install ffmpeg
      
on MacOS to install it. Note that you will still be able to produce the images even if ```ffmpeg``` is not installed. You can then use another tool to assemble them.


## Non-random initial conditions

If you decide to use non-random initial conditions, then you also need to give NcorpiON a file ```init.txt``` containing the initial conditions. This file must be put at the input/output location before launching the simulation.

It must contain one line per body and 8 columns, namely (x, y, z, vx, vy, vz, mass, radius) if you decide to give your initial conditions in cartesian coordinates or (a, e, i, nu, omega, Omega, mass, radius) if you decide to give your initial conditions in elliptic elements.

If your simulation has a central mass, then its initial conditions must not be given in ```init.txt``` because they are determined in the parameter file when you set up your system of units. The initial conditions must be given in units consistent with what you decided in the parameter file.


## Why NcorpiON

Other N-body softwares exist, like REBOUND, and you might wonder what NcorpiON has that they don't. NcorpiON can use the algorithm FalcON of Walter Dehnen to compute mutual gravitational interactions and to detect collisions between bodies in $\mathcal{O}(N)$.
This FFM algorithm is a considerable improvement of the Barnes-Hut tree code, significantly faster for the same precision, and that conserves the total momentum of the system to machine precision (Barnes-Hut doesn't).
Using NcorpiON for a N-body numerical integration is thus very relevant for a system with large N where precision on the computation of mutual gravity is still needed.

In a N-body simulation with large N, collisions can be numerous and can lead to the fragmentation of the bodies. NcorpiON features a fragmentation module based on the litterature on crater scaling and ejecta models
to realistically render impacts between the bodies. Using NcorpiON is therefore relevant in a system with numerous collisions, like a protoplanetary disk for example.


## Official website

The official website of NcorpiON is available at <https://ncorpion.com>


## Acknowledgments

If you use this code in a scientific publication, please cite NcorpiON.


## Contributors

- Jérémy Couturier, University of Rochester, <https://jeremycouturier.com>

NcorpiON is open source and you are encouraged to contribute to it 


## License

NcorpiON is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

NcorpiON is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with NcorpiON.  If not, see <http://www.gnu.org/licenses/>.
