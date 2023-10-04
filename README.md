# Welcome to NcorpiON

NcorpiON is an N-body software dedicated to the simulation of collisional and fragmenting systems.


## Features

- Written in C
- Four built_in modules for mutual interactions :
      - Brute-force O(N^2) method
      - Barnes-Hut O(N ln N) tree code
      - Mesh O(N) algorithm
      - FalcON O(N) Fast Multipole Method
- A built-in fragmentation model that can realistically handle violent collisions
- A python add-on to produce animated gifs of the simulations
- Requires only a C compiler (e.g. gcc) to run
- Completely open source. All features are available in the public github repository.


## Installation

You can install NcorpiON with the git command

      - git clone
    
Later on, you can update NcorpiON by running

      - git pull

To run a simulation, first update the file src/parameters.h and then :

      - Without producing an animation :
      - Producing an animation :

 
## Documentation

The full documentation of NcorpiON is available at <https://ncorpion.com>


## Acknowledgments

If you use this code or parts of this code for results presented in a scientific publication, please cite NcorpiON.


## Contributors

- Jérémy Couturier, University of Rochester, <jeremycouturier@rochester.edu>


NcorpiON is open source and you are encouraged to contribute to it 


## License

NcorpiON is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

NcorpiON is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with NcorpiON.  If not, see <http://www.gnu.org/licenses/>.
