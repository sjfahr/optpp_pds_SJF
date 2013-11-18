%This function accepts the FOV, matrix, and scaling from an image and returns a
%suggested domain for the Bioheat1D.m function. The FOV is a structure
%describing the number of x and y pixels and the number of slices. The the output domain has
%tighter spacing as the scaling increases. I.e. the number of points is
%linear with the scaling. Scaling is a 3-field structure too (for x,y,z).
%'mod_dom' is the size of the modeled volume in meters.

function [dom,MRTI_pix,mod_pix]=Pixel_size_rounding(FOV,matrix,scaling,mod_dom)

%Fraction of FOV that is modeled
fraction_mod.x=mod_dom.x/FOV.x;
fraction_mod.y=mod_dom.y/FOV.y;
fraction_mod.z=mod_dom.z/FOV.z;

%Define domain
dom.x=mod_dom.x;
dom.y=mod_dom.y;
dom.z=mod_dom.z;
dom.pointx=matrix.x*scaling.x*fraction_mod.x;
dom.pointy=matrix.y*scaling.y*fraction_mod.y;
dom.pointz=matrix.z*scaling.z*fraction_mod.z;

%Find the MRTI pixel size
MRTI_pix.x=FOV.x/matrix.x;
MRTI_pix.y=FOV.y/matrix.y;
MRTI_pix.z=FOV.z/matrix.z;

%Modeled pixel size
mod_pix.x=dom.x/dom.pointx;
mod_pix.y=dom.y/dom.pointy;
mod_pix.z=dom.z/dom.pointz;

end