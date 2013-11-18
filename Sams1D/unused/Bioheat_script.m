% This is the updated Bioheat_script that should be used with DF's DAKOTA
% run.

cd '/FUS4/data2/sjfahrenholtz/non_MATLAB_git/optpp_pds/optpp_pds'

load('optpp_pds.in.1.mat');

% Write every string as a number

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

% This is a run script for the Bioheat1D.m file. First load and parse the
% power. The convention for x,y,z arrays is [y x z] with y being vertical in
% the image (row in array); x is horizontal (column in array); z is depth
cd /FUS4/data2/BioTex/BrainNonMDA/processed/Patient0002/000/laser_log
load 'power_log.txt';
[P,~,delta_P] = power_parser(power_log);
%P=[1 -3; 5 0;15 3; 20 6; 25 9; 30 12; 35 15];

% Get the laser registration
cd /FUS4/data2/BioTex/BrainNonMDA/processed/Patient0002/000/laser
reg_laser = csvimport ( 'reg_laser.csv' );
reg = cell2mat( reg_laser(2,:) );           % Grabs the numeric portion of the registration cell array

% 2D vs 3D flag; 0=2D, 1=3D
flag3D=0;

% Define the domain and scaling
mod_point.x=171;  % x position on image
mod_point.y=171;
mod_point.z=1;

mod_center=[floor(mod_point.y/2) floor(mod_point.x/2) 1]; % [y x z]

matrix.x=256;
matrix.y=256;
matrix.z=1;

% For now, x and y scaling must be equal; z = 1
scaling.x=1;
scaling.y=1;
scaling.z=1;

FOV.x=0.25; % change for image/dataset
FOV.y=0.25;
FOV.z=0.007; 

% Build the domain
[dom,MRTI_pix,mod_pix]=modeled_domain(FOV,matrix,scaling,mod_point);
% Display the domain
dom
MRTI_pix
mod_pix

% Define the source; source.n must be greater than 1 to make sense. Odd is
% better than even
source.n=5;
source.length=0.01;  %~0.033 is when n=5 is visible
source.laser=linspace((-source.length/2),(source.length/2),source.n);

% Run the Bioheat model with the unique powers
tic
[tmap_unique]=Bioheat1D(P,dom,source,w,k,g,mu_a,mu_s,probe_u,robin_co);
toc

tmap_unique=tmap_unique+37;
tmap_unique(:,:,:,1)=37;
% Make the full temperature history with the unique tmaps
[tmap]=Build_tmap_history(tmap_unique,delta_P);

clear FOV matrix power_log delta_P mod_point unique_P source mod_pix MRTI_pix dom

% Set delta time 'dt' for dose and run it
dt = 2.5;
[w,Iso]=ArrDose_for_1D(tmap,dt);
%[w_unique,Iso_unique]=ArrDose_for_1D(tmap_unique);


figure(1);imagesc(tmap_unique(:,:,1,1,2));
figure(2);imagesc(tmap_unique(:,:,1,1,4));
figure(3);imagesc( Iso );

w ( w > 4) = 4;  %Put a cap on the dose


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% This section does rotations,
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% which are irrelevant because
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% the model lesion is huge
% start = [50 144]; %[x y] from MRTI image
% stop = [116 125];
% center = [112 126];

%Find the angle that links the start and stop; For this dataset, it's only an angle about the z-positive axis
%diff = stop - start;  %For 2-D [y_diff x_diff]; for 3-D [y_diff x_diff z_diff]

% angle = [0 0 0]; %Initialize angle [y_angle x_angle z_angle]
% if flag3D==0  %For 2-D
%     angle(3)=atan(diff(1)/diff(2)); % z_angle
% else     %For 3-D
%     angle(1)=atan(diff(3)/diff(2)); % x_angle
%     angle(2)=atan(diff(3)/diff(1)); % y_angle
%     angle(3)=atan(diff(1)/diff(2)); % z_angle
% end

% x_angle=0;
% y_angle=0;
% z_angle=angle;

%[aa] = threeDrotate(w,angle(2),angle(1),angle(3));  %Rotated dose from model

%aa=imresize(aa,1/scaling.x); %For now, scaling.x must = scaling.y

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

aa = imresize (w , 1/scaling.x);


cd /FUS4/data2/BioTex/BrainNonMDA/processed/Patient0002/010/matlab
load 'arrheniusDose.mat'
MRTI_dose_size=size(arrheniusDose.mean);

% Need to make the model timing equal the MRTI timing
% 
% jj=0;
% for ii=1:MRTI_dose_size(4);
%     jj=jj+2;
%     model_dose_rot_time(ii)=aa(jj);
% end

%Make the giant matrix; use aa_size because the rotated matrix is bigger
%than the original w matrix
aa_size=size(aa);
size_diff=[(MRTI_dose_size(1)-aa_size(1)) (MRTI_dose_size(2)-aa_size(2))];
upper_left_mod = zeros((size(aa,1)+size_diff(1)),(size(aa,2)+size_diff(2)));
upper_left_mod(1:size(aa,1),1:size(aa,2)) = aa;
%upper_left_mod = padarray(aa, [size_diff(1) size_diff(2)]);  %Dose


%Define intervals that will be written
% y_range = [ (center(1)-floor(aa_size(1)/2)) (center(1)+floor(aa_size(1)/2))];
% x_range = [ (center(2)-floor(aa_size(2)/2)) (center(2)+floor(aa_size(2)/2))];

% Uses registration instead of center
y_range = [ (reg(1)-floor(aa_size(1)/2)) (reg(1)+floor(aa_size(1)/2))];
x_range = [ (reg(2)-floor(aa_size(2)/2)) (reg(2)+floor(aa_size(2)/2))];

%matched_mod = zeros(arrheniusDose.mean(1),arrheniusDose.mean(2),);
matched_mod = zeros (MRTI_dose_size(1), MRTI_dose_size(2));
matched_mod ( y_range(1):y_range(2), x_range(1):x_range(2) ,:,:) = upper_left_mod( 1:aa_size(1) , 1:aa_size(2) );

figure(4); imagesc(matched_mod,[0 1]);
figure(5); imagesc( arrheniusDose.mean , [0 1] );