% This is the updated Bioheat_script that should be used with DF's DAKOTA
% run.

function [metric,dom,MRTI_pix,mod_pix] = sandbox_obj_fxn ( path22, iteration );



%Write every string as a number
%g=0.9;      % Optical anisotropy
k= 0.5   ;    % Thermal conductivity
mu_a= 50;     % Optical absorption
mu_s= 8000;     % Optical scattering
probe_u= 21;    % Initial probe temperature
robin_co=0;  % Value of Robin boundary coefficient
perfusion=6;         % Blood perfusion
x_disp=0.096445;    % Horizontal displacement
x_rot=0;        % Rotation about horizontal axis
y_disp=0.1196;     % Vertical displacement
y_rot=0;        % Rotation about vertical axis
z_disp=0;     % Depth displacement
z_rot=0;        % Rotation about depth axis

g=0.7;      % Optical anisotropy
%g=0.2;

cd (path22);
cd '/FUS4/data2/BioTex/BrainNonMDA/processed/Patient0002/000/laser_log'
load 'power_log.txt';
[P,~,delta_P] = power_parser(power_log);

% Define the domain and scaling
mod_point.x=171;  % x position on image
mod_point.y=171;
mod_point.z=1;

% Import the VTK header info
cd (path22);
cd '/FUS4/data2/BioTex/BrainNonMDA/processed/Patient0002/000/laser'
FOV_import = csvimport ( 'FOV.csv' );
fov = FOV_import ( 2, : );
fov = cell2mat (fov);

matrix.x = fov(1);
matrix.y = fov(2);
matrix.z = fov(3);

spacing.x = fov(4);
spacing.y = fov(5);
spacing.z = fov(6);

FOV.x = matrix.x * spacing.x;
FOV.y = matrix.y * spacing.y;
FOV.z = matrix.z * spacing.z;

% For now, x and y scaling must be equal; z = 1
scaling.x=1;
scaling.y=1;
scaling.z=1;

% Build the domain
[dom,MRTI_pix,mod_pix]=modeled_domain(FOV,matrix,scaling,mod_point);



% Define the source; source.n must be greater than 1 to make sense. Odd is
% better than even
source.n=5;
source.length=0.01;  %~0.033 is when n=5 is visible
source.laser=linspace((-source.length/2),(source.length/2),source.n);

% Run the Bioheat model with the unique powers
[tmap_unique]=Bioheat1D(P,dom,source,perfusion,k,g,mu_a,mu_s,probe_u,robin_co);

tmap_unique=tmap_unique+37;
tmap_unique(:,:,:,1)=37;
% Make the full temperature history with the unique tmaps
% [tmap]=Build_tmap_history(tmap_unique,delta_P);
% 
% % Set delta time 'dt' for dose and run it
% dt = 2.5;
% [www,Iso]=ArrDose_for_1D(tmap,dt);
% www ( www > 4 ) = 4;  %Put a cap on the dose

figure(1);imagesc(tmap_unique(:,:,1,1,2));
%figure(2);imagesc(tmap_unique(:,:,1,1,4));
%figure(2);imagesc( Iso );

% Calculate the number of pixels that must be shifted
pixel_reg.x = round (x_disp / mod_pix.x);
pixel_reg.y = round (y_disp / mod_pix.y);
pixel_reg.z = round (z_disp / mod_pix.z);

%aa = imresize (www , 1/scaling.x);
%aa = imresize (tmap , 1/scaling.x);
aa = imresize (tmap_unique , 1/scaling.x);

bb = aa(round(matrix.x/2),round(matrix.y/2),1,:);

[cc,dd] = max (bb);

cd '/FUS4/data2/BioTex/BrainNonMDA/processed/Patient0002/000/matlab/'
load 'arrheniusDose.mat'
MRTI_dose_size=size(arrheniusDose.mean);

%aa_size=size(aa);
aa_size = size ( aa(:,:,1,dd) );
size_diff=[(MRTI_dose_size(1)-aa_size(1)) (MRTI_dose_size(2)-aa_size(2))];
upper_left_mod = zeros((size(aa,1)+size_diff(1)),(size(aa,2)+size_diff(2)));
upper_left_mod(1:size(aa,1),1:size(aa,2)) = aa(:,:,1,dd);

%Define intervals that will be written
x_range = [ (pixel_reg.x - floor(aa_size(1)/2)) (pixel_reg.x + floor(aa_size(1)/2))];
y_range = [ (pixel_reg.y - floor(aa_size(2)/2)) (pixel_reg.y + floor(aa_size(2)/2))];

roi_x   = [ (pixel_reg.x - 50) (pixel_reg.x + 50) ];  %For model
roi_y   = [ (pixel_reg.y - 50) (pixel_reg.y + 50) ];

roi_x_MRTI   = [ (pixel_reg.x - 20) (pixel_reg.x + 20) ];  %For MRTI
roi_y_MRTI   = [ (pixel_reg.y - 20) (pixel_reg.y + 20) ];

%matched_mod = zeros(arrheniusDose.mean(1),arrheniusDose.mean(2),);
matched_mod = zeros (MRTI_dose_size(1), MRTI_dose_size(2));
matched_mod ( x_range(1):x_range(2), y_range(1):y_range(2) ,:,:) = upper_left_mod( 1:aa_size(1) , 1:aa_size(2) );

model_Iso = (matched_mod >57) ;

model_Iso ( 1:roi_x(1), : ) = 0;
model_Iso ( roi_x(2):end,: )  = 0;
model_Iso ( :, 1:roi_y(1) ) = 0;
model_Iso ( :,roi_y(2):end) = 0;

MRTI_Iso = (arrheniusDose.mean >1) ;

MRTI_Iso ( 1:roi_x_MRTI(1), : ) = 0;
MRTI_Iso ( roi_x_MRTI(2):end, : )  = 0;
MRTI_Iso ( :, 1:roi_y_MRTI(1) ) = 0;
MRTI_Iso ( :,roi_y_MRTI(2):end) = 0;


sum_Iso = model_Iso + MRTI_Iso;
diff_Iso= model_Iso - MRTI_Iso;

figure(3); imagesc(matched_mod);
figure(4); imagesc(model_Iso);
figure(5); imagesc(MRTI_Iso);
figure(6); imagesc(sum_Iso);
figure(7); imagesc(diff_Iso);

metric = abs(sum(sum(diff_Iso)));

cd (path22);
%[metric] = Iso_obj ( matched_mod, MRTI_Iso);

output_param = 'optpp_pds.out.';
index = num2str(iteration);
output_filename = strcat( output_param, index);

csvwrite ( output_filename, metric);

end