% This will be the objective function for the ISMRM 2014 abstract w/ DF

function [metric] = dose_obj ( model_D, MRTI_D);

% Set the dose limits to look at stuff 2 and below.
model_D (model_D > 2) = 2;
MRTI_D (MRTI_D > 2) = 2;

% The pixel-wise metric in the sqrt; then summed to make one number
%metric = sum( sum( sqrt ( abs( model_D.^2 - MRTI_D.^2 ))));

% Alternative metric
metric = sum( sum( abs ( model_D - MRTI_D )));

end