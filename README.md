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

1. In MATLAB, please run  **Simulation.m** to generate data (it may take several minutes), or directly use **Data.mat** (No need k-wave toolbox).
2. run **Run.m** to reconstruction the image.

### Note

1. There are detailed description of the code's functionality in the source code (**Simulation.m**, **Run .m**, **PS_3D_NUFFT_Fast.m** **PS_3D_NUFFT.m**, **InterpNUFFT.m**, **InterpSinc.m**)
2. Using **PS_3D_NUFFT_Fast** for this two-layer medium demo. For other multi-layered medium, please use **PS_3D_NUFFT**

2. In  **Simulation.m**, you can change parameters **source.p0 (optoacoustic sources)  and medium.sound_speed ** to generate yourself simulation data.  Different **medium. sound_speed** should set corresponding parameters **c and layer** in **Run.m**



