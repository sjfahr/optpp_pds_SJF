%This script is meant to be run in SI units.
%As a function of power, P, this makes a prediction of the temperature
%within the domain, dom. 'dom' is a structure representing the simulated domain with
%three spacial dimensions in meters and the number of points in each
%dimension (a total of 6 fields). The number of points in 'dom' should be
%odd in order for the centroid to be on the origin and a point to be on the
%origin.

%'source' is a structure with the details of the source. The 'n' field is
%the number of sources. The 'length' is the diffusing tip length, commonly
%0.01 m.

function [tmap]=Bioheat1D_array(unique_P,dom,dom_point,source);

%List of space and time details
R1=0.0015/2; % (m) R1 is the distance from the isotropic laser source point and the edge of the fiber
R2=1; % (m) R2 is the maximum edge of the domain;
time=length(P); %Returns how many timesteps will be calculated
P=P/source.n;  %Scales the power to the number of source points

%List of constants
mua=500; % 1/m
mus=14000; %1/m
g=0.88; % Unity
k=0.527; % W/(m * K)
w=6; % kg / (m^3 * s)
u0=37+273.15; % K
ua=37+273.15; % K

%Points structure
points=linspace(-dom/2,dom/2,dom_point);

%Initialize tmap, t_sample, and r
tmap=zeros(dom_point);
t_sample=zeros(source.n,time);
r=zeros(source.n,1);

%Spatial locations of the sources; My convention is that the long axis of
%the laser is parallel to the y-axis.
laser=linspace(-source.length/2,source.length/2,source.n);

%Giant for loop vector for each source, calculate the t_sample
for i=1:dom_point(1)
    
    for ii=1:dom_point(2)
        
        for iii=1:dom_point(3)
            
            for j=1:source.n
                r(j)=sqrt(points.x(i)^2+(laser(j)-points.y(ii))^2+points.z(iii)^2);
                
                for jj=1:max(size(unique_P))
                    [t_sample(j,jj)]=sammodel1D(u0,ua,k,w,unique_P(jj),r(j),mua,mus,R1,R2,g);
                    tmap(i,ii,iii,jj)=sum(t_sample(j,jj,1));
                end
            end
        end
    end
end


end