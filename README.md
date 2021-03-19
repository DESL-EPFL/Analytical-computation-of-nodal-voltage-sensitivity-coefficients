# Analytical-computation-of-nodal-voltage-sensitivity-coefficients
We provide a MATLAB code for the analytical computation of nodal voltage sensitivity coefficients based on the proposed method in [1]. The numerical simulation presented in [1] can be ran by running the main.m (see Files for details). Other subfunctions are available and can work independently (see Files for details).

## Files
* [main.m] function that runs the simulation presented in [1]
* [SC_Voltage.m] function that computes complex, magnitude and angle nodal voltage sensitivity coefficients using the method presented in [1]
* [3ph Load Flow - IEEE34/Ymatrix_3ph.m] function that computes the compound admittance matrix. It takes as an input an excel file with a specific format.
* [3ph Load Flow - IEEE34/IEEE34Feeder.xlsx] Excel file formatted to fit the function Ymatrix_3ph.m. It contains the line parameters and base values values of the IEEE34 feeder [2] without the presence of the voltage regulators.
* [3ph Load Flow - IEEE34/NR_VoltageDepJacob.m] function that performs the Newton-Raphson (N-R) algorithm and outputs the network state (i.e. phase-to-ground nodal voltages at each node). It supports multiphases automitically depending on the inputted compound admittance matrix.
* [3ph Load Flow - IEEE34/NR_VoltageDepJacob_iteration.m] function that performs a sub-iteration of the Netwon-Raphson algorithm where the Jacobian of the k-th iteration is computed and the network state is updated using the inverted Jacobian and power mismatches of the k-th iteration.
* [LoadFlow input examples/E_star_IEEE34_example_noPQnodes.mat] setpoints for the magnitudes of nodal voltages for pv nodes for a simulation with no pv nodes. Input needed for N-R algorithm.
* [LoadFlow input examples/E_star_IEEE34_example_withPQnodes.mat] setpoints for the magnitudes of nodal voltages for pv nodes for a simulation with pv nodes. Input needed for N-R algorithm.
* [LoadFlow input examples/S_star_IEEE34_example.mat] setpoints for the apparent nodal power injections for pq nodes. Input needed for N-R algorithm.

## Simulation notes
* Running the function main.m should directly output all the results of the work in [1]. 
* All throughout the main.m function there are configurable inputs. These configurations can be set in the code line following the commented line: % (!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!) INPUT (!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!)

## Software 
The following software is required to run the model:
* MATLAB ( in principle any version however tested on MATLAB version 2020a Update 1 (9.8.0.1359463) )

## References 
[1] reference will be published soon. For more info contact the author at sherif.fahmy@epfl.ch
[2] Kersting, W. H. (1991). Radial distribution test feeders. IEEE Transactions on Power Systems, 6(3), 975-985.
