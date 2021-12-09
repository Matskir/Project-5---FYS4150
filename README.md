# Project-5---FYS4150


The file header_proj5.hpp is the header for this project, and contains all needed functions to set up and run the simulation of the probability function over 2D space and time.


The C++ code file proj5_fys4150.cpp included the main function of this project, where the needed functions from the header are utilised to simulate the probability functions and save them as slices in a data cube. The .cpp file should be run as:

g++ -O3 proj5_fys4150.cpp -o proj5_fys4150.exe -larmadillo -Wno-narrowing

./proj5_fys4150.exe


The python script read_plot_cube.py reads the data cube and makes all the necessary plots used in the report for the specified potential as well as produces an animation of the probability function over time. The .bin file needed to run the script is not here included because of its size, but could easily be made by just running proj5_fys4150.cpp first. The python script should be run as:

python read_plot_cube.py


The animations folder contains the animations of the probability function for the different potentials; no and double slit. 
