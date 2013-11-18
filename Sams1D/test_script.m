% This is the superscript that has the paths for all of the patient data.

% Paths.
tic
setenv ( 'PATH22' , pwd);
path22 = getenv ( 'PATH22' );

% cd (path22)
% load index.txt
index=1;
[metric] = test_obj_fxn ( path22, index );

% index = index + 1;
% csvwrite ('index.txt' , index);

%metric
% dom
% MRTI_pix
% model_pix
toc