# virion-spike-model

Generate HIV and SIV virion models with a given number of spikes and measure distances between spikes.

Author: Wonmuk Hwang ( hwm@tamu.edu ). see LICENSE.

How to cite and for details, see:  JOURNAL REF TO BE ADDED UPON PUBLICATION

-------------------------------------------------------

## Directories:
- [./inp1] Contains input files for execution of the code.
- [./example1,./data1]: Directories where outputs are written. Initially empty.
- [vmd]: VMD files for visualization.
-------------------------------------------------------

## How to run:

The following operations will render models in VMD and generate graphs used in Fig. 6 of the reference cited above.


### 0) Dependent package: Visual Molecular Dynamics (VMD).

VMD can be downloaded from: https://www.ks.uiuc.edu/Research/vmd/

### 1) To create a model virion to be visualized by VMD:

- In virion_main.cpp, uncomment line 13:

#define VMD // uncomment to write vmd file

- Compile:

g++ vector_stl.cpp ftn_virion.cpp virion_main.cpp -O3 -o a.out

- Run: The following will write outputs to ./example1

./a.out inp1/demo14g1.dat

./a.out inp1/demo14g3.dat

./a.out inp1/demo73g1.dat

./a.out inp1/demo73g3.dat

- Visualize using VMD:

cd vmd

vmd -size 1000 1000 -e view_demo_14g1-1.vmd

vmd -size 1000 1000 -e view_demo_14g3-1.vmd

vmd -size 1000 1000 -e view_demo_73g1-0.vmd

vmd -size 1000 1000 -e view_demo_73g3-0.vmd

### 2) To measure inter-spike distance distribution:

- In virion_main.cpp, comment out line 13:

//#define VMD // uncomment to write vmd file

- Compile:

g++ vector_stl.cpp ftn_virion.cpp virion_main.cpp -O3 -o a.out

- Run: The following will write outputs to ./data1

./a.out inp1/meas14g1.dat &

./a.out inp1/meas14g3.dat &

./a.out inp1/meas73g1.dat &

./a.out inp1/meas73g3.dat

- Create histograms:

python plot_histo_rcut1r.py

output: histo_rcut1r.pdf, ./data1/histo*_rcut.dat

python plot_histo_nn1.py

output: histo_nn1.pdf, ./data1/histo{14,73}.dat
