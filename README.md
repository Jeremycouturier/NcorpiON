# Welcome to Ncorpi$\mathcal{O}$N

Ncorpi$\mathcal{O}$N is an N-body software dedicated to the simulation of collisional and fragmenting systems. It is very fast and can realistically handle fragmentations


## Features

- Written in C
- Four built-in modules for mutual interactions, namely
- Brute-force $\mathcal{O}(N^2)$ method
- Barnes-Hut $\mathcal{O}(N ln N)$ tree code
- Mesh $\mathcal{O}(N)$ algorithm
- FalcON $\mathcal{O}(N)$ Fast Multipole Method
- A built-in fragmentation model that can realistically handle violent collisions
- A python add-on to produce animations of the simulations
- Requires only a C compiler (e.g. gcc) to run
- Completely open source. All features are available in the public github repository.


## Installation

You can install Ncorpi$\mathcal{O}$N with the ```git``` command

	git clone git@github.com:Jeremycouturier/Ncorpi$\mathcal{O}$N.git
	
or if you prefer, you can simply download this repository to your computer. Using ```git``` is recommended so you can later on get updates
by running ```git pull```

Note that the production of animations requires ```ffmpeg``` to be installed on your system. Run

	sudo apt-get update && sudo apt-get install ffmpeg
	
on Debian based distros or

	sudo dnf update && sudo dnf install ffmpeg
	
on Red Hat based distros to install it. You will also need python with matplotlib and numpy for the animations.


## Simulation

To run a simulation, you need to update the file ```src/parameters.h``` (see <https://Ncorpi$\mathcal{O}$N.com/#setup> for details)

In order to simulate and to produce an animation, first run

	./Ncorpi$\mathcal{O}$N.sh number_of_images some_path
	
where the path ```some_path``` is the same as given in ```src/parameters.h```. This will run the simulation and create the images of the animation.
The number of images depends on the length of the simulation and on the
frequency of outputs chosen in the file ```src/parameters.h```. Then, once all the images have been produced, assemble the animation with

	./Ncorpi$\mathcal{O}$N_animation.sh some_path
	
This will return two files ```Ncorpi$\mathcal{O}$N.mp4``` and ```Ncorpi$\mathcal{O}$N.gif``` in ```some_path/gif```

In order to simulate without producing any animation, run

	make clean && make && make clean && ./Ncorpi$\mathcal{O}$N
	
The files of the numerical simulation can then be found in ```some_path```.


## Documentation

The full documentation of Ncorpi$\mathcal{O}$N is available at <https://Ncorpi$\mathcal{O}$N.com>


## Acknowledgments

If you use this code in a scientific publication, please cite Ncorpi$\mathcal{O}$N.


## Contributors

- Jérémy Couturier, University of Rochester, <https://jeremycouturier.com>

Ncorpi$\mathcal{O}$N is open source and you are encouraged to contribute to it 


## License

Ncorpi$\mathcal{O}$N is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Ncorpi$\mathcal{O}$N is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Ncorpi$\mathcal{O}$N.  If not, see <http://www.gnu.org/licenses/>.
