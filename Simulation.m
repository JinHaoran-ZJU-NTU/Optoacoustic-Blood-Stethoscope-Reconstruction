% Simulations In Three Dimensions Example
%
% This example provides a simple demonstration of using k-Wave for the
% simulation and detection of the pressure field generated by an initial
% pressure distribution within a three-dimensional heterogeneous
% propagation medium. It builds on the Homogeneous Propagation Medium and
% Heterogeneous Propagation Medium examples.    
%
% author: Bradley Treeby
% date: 1st July 2009
% last update: 13th April 2017
%  
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2017 Bradley Treeby

% This file is part of k-Wave. k-Wave is free software: you can
% redistribute it and/or modify it under the terms of the GNU Lesser
% General Public License as published by the Free Software Foundation,
% either version 3 of the License, or (at your option) any later version.
% 
% k-Wave is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
% more details. 
% 
% You should have received a copy of the GNU Lesser General Public License
% along with k-Wave. If not, see <http://www.gnu.org/licenses/>. 

clearvars;

% =========================================================================
% SIMULATION
% =========================================================================

% create the computational grid
PML_size = 20;          % size of the PML in grid points
Nx = 256;            % number of grid points in the x direction
Ny = 256;            % number of grid points in the y direction
Nz = 64;            % number of grid points in the z direction
dx = 0.1e-3;        % grid point spacing in the x direction [m]
dy = 0.1e-3;        % grid point spacing in the y direction [m]
dz = 0.1e-3;        % grid point spacing in the z direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

% define the properties of the propagation medium
medium.sound_speed = 1450 * ones(Nx, Ny, Nz);	% [m/s]
medium.sound_speed(:, :, 1:Nz/2) = 1750;        % [m/s]
medium.density = 950 * ones(Nx, Ny, Nz);       % [kg/m^3]
medium.density(:, :, 1:Nz/2) = 1150;          % [kg/m^3]

% % create initial pressure distribution using makeBall
% ball_magnitude = 10;    % [Pa]
% ball_x_pos = 38;        % [grid points]
% ball_y_pos = 32;        % [grid points]
% ball_z_pos = 32;        % [grid points]
% ball_radius = 5;        % [grid points]
% ball_1 = ball_magnitude * makeBall(Nx, Ny, Nz, ball_x_pos, ball_y_pos, ball_z_pos, ball_radius);
% 
% ball_magnitude = 10;    % [Pa]
% ball_x_pos = 20;        % [grid points]
% ball_y_pos = 20;        % [grid points]
% ball_z_pos = 20;        % [grid points]
% ball_radius = 3;        % [grid points]
% ball_2 = ball_magnitude * makeBall(Nx, Ny, Nz, ball_x_pos, ball_y_pos, ball_z_pos, ball_radius);
% 
% source.p0 = ball_1 + ball_2;

% load the initial pressure distribution from an image and scale
p0_magnitude = 2;
p0 = p0_magnitude * loadImage('EXAMPLE_source_two.bmp');



% resize the input image to the desired number of grid points
p0 = resize(p0, [Nx, Ny]);
p0(p0>0.01) = 1;
p0(p0<=0.01) = 0;
% smooth the initial pressure distribution and restore the magnitude
% p0 = smooth(kgrid, p0, true);

% assign to the source structure
source.p0 = zeros(Nx,Ny,Nz);
source.p0(:,:,60) = p0;
% define a series of Cartesian points to collect the data
sensor.mask = zeros(size(medium.sound_speed));
sensor.mask(:,:,10) = 1;

% input arguments
input_args = {'PlotLayout', true, 'PlotPML', false,'PMLSize', PML_size, ...
    'DataCast', 'single', 'CartInterp', 'nearest'};

% run the simulation
sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});

sensor_data = reshape(sensor_data,[Nx,Ny,size(sensor_data,2)]);
save('Data.mat');

