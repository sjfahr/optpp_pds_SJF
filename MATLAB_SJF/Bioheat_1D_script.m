%This is the updated Bioheat_script that should be used with DF's DAKOTA
%run.

cd '/FUS4/data2/sjfahrenholtz/non_MATLAB_git/optpp_pds/optpp_pds'

load('optpp_pds.in.1.mat');

%Write every string as a number

g=str2num(anfact_healthy);      % Optical anisotropy
k=str2num(k_0_healthy);         % Thermal conductivity
mu_a=str2num(mu_a_healthy);     % Optical absorption
mu_s=str2num(mu_s_healthy);     % Optical scattering
probe_u=str2num(probe_init);    % Initial probe temperature
robin_co=str2num(robin_coeff);  % Value of Robin boundary coefficient
w=str2num(w_0_healthy);         % Blood perfusion
x_disp=str2num(x_displace);     % Horizontal displacement
x_rot=str2num(x_rotate);        % Rotation about horizontal axis
y_disp=str2num(y_displace);     % Vertical displacement
y_rot=str2num(y_rotate);        % Rotation about vertical axis
z_disp=str2num(z_displace);     % Depth displacement
z_rot=str2num(z_rotate);        % Rotation about depth axis

cd /FUS4/data2/BioTex/BrainNonMDA/processed/Patient0002/010/laser_log
load 'power_log.txt';

[P,~,delta_P]=power_parser(power_log);

%Define the domain and scaling
modeled_points = [71 71 1]; %[x y z] but remember MATLAB has origin in upper left and x-positive points down and y-positive points right
image_matrix = [256 256 1];
scaling = [2 2 1];
FOV = [0.25 0.25 0.007];

%Build the domain
[dom,MRTI_pix,mod_pix]=modeled_domain_arrays(FOV,image_matrix,scaling,modeled_points);
%Display the domain
dom
MRTI_pix
mod_pix

%Define the source; source.n must be greater than 1 to make sense. Odd is
%better than even
source.n=5;
source.length=0.035;  %~0.033 is when n=5 is visible
source.laser=linspace((-source.length/2),(source.length/2),source.n);

%Run the Bioheat model with the unique powers
tic
[tmap_unique]=Bioheat1D(P,dom,source);
toc