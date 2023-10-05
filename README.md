# Welcome to NcorpiON

NcorpiON is an N-body software dedicated to the simulation of collisional and fragmenting systems.


## Features

- Written in C
- Four built-in modules for mutual interactions, namely
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

      git clone git@github.com:Jeremycouturier/NcorpiON_code.git
    
Later on, you can update NcorpiON by running

      git pull

To run a simulation, first update the file src/parameters.h

Then to run a simulation and to produce an animation

      ./ncorpion.sh number_of_images_to_be_produced path_where_to_produce_them
      cd path_where_to_produce_them/gif && ffmpeg -i %d.png -vf palettegen palette.png && ffmpeg -framerate 25 -i "%d.png" -i palette.png -lavfi paletteuse ncorpion.gif && ffmpeg -framerate 25 -i "%d.png" -i palette.png -lavfi paletteuse ncorpion.mp4 

To run the simulation without producing an animation

      make clean && make && make clean && ./ncorpion


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
