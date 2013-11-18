function [metric] = Iso_obj ( model_Iso, MRTI_D);

MRTI_Iso = (MRTI_D > 1);

metric = sum( sum( abs ( model_Iso - MRTI_Iso )));
end