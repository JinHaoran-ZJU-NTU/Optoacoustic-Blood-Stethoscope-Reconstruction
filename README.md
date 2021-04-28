# Optoacoustic Blood Stethoscope Reconstruction
 This project reports the phase shift non-uniform Fast Fourier transform (PS-NUFFT) method for optoacoustic image reconstruction

### Software dependence

1. MATLAB>= R2017. 
2. k-wave toolbox (Version 1.3) --- http://www.k-wave.org/

### Installation guide of k-wave toolbox

k-wave toolbox provides a framework for generate the simulation data.

1. Download k-wave MATLAB toolbox.zip
2. unzip k-wave MATLAB toolbox.zip
3. In MATLAB, add k-Wave folder (with subfolders) into the '**set path**'

### Instructions to run on the data

1. In MATLAB, please run  **Simulation.m** to generate data (it may take several minutes), or directly use **Data.mat**.
2. run **PS-NUFFT .m** to reconstruction the image.



There are detailed description of the code's functionality in the source code (**Simulation.m**, **PS-NUFFT .m**, **Omega_K_3D_NUFFT.m**, **InterpNUFFT.m**, **InterpSinc.m**)

In  **Simulation.m**, you can change parameters **source.p0 (optoacoustic sources) ** to generate yourself simulation data. 



