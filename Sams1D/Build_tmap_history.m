%This function accepts the unique tmaps and returns the tmaps that
%correspond to MRTI.

function [tmap]=Build_tmap_history(tmap_unique,delta_P);

tmap=zeros(size(tmap_unique,1),size(tmap_unique,2),size(tmap_unique,3),size(delta_P,1));

j=1;
for i=1:length(delta_P)

    if delta_P(i) ~= 0
        j=j+1;
    end
    
    tmap(:,:,:,i)=tmap_unique(:,:,:,j);
    
end
    
end
% 
% j=1;
% i=1;
% 
% 
% 
% while i <= size(delta_P,1)
%     if delta_P(i,1)==0
%         tmap(:,:,:,i)=tmap_unique(:,:,:,j);
%     else
%         
%         while j<= size(tmap_unique,4)
% 
%             tmap(:,:,:,i)=tmap_unique(:,:,:,j);
%             i=i+1;
%             j=j+1;
%         end
%     end
%     i=i+1;
% end


% for i=1:size(delta_P,1)
%     if delta_P(i,1)==1
%         tmap(:,:
%         j=j+1;
% 
% 
% P_percent=unique_P*100/15;
% 
% for j=size(intervals)
% for i=size(power_log,1)
% while i<intervals
% tmap(:,:,:,1:P_percent(i)=tmap_unique
% 
% time=size(power_log,1);

