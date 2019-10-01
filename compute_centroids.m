function [centroidsscaled,height,prob_det]=compute_centroids(X,Y,Z,thresh,debug)
% function compute centroids takes the three 2d vector inputs used on for
% the matlab function contourf(X,Y,Z) and compute the centroid for them
% threshold is the threshold for probability of detection
% debug is a boolean of whether or not you want debugging output
% note add the following lines uncommented below the contourf() function
% call and before the save figure line to add an overlay plot to any
% whatever type of plot these are:
%%% add it dimitri stuff here for centroid plotting
%[centroidsscaled3,height3,prob_det3]=compute_centroids(X,Y,rdmap_compansated',threshold,debug);
%hold on;
%plot(centroidsscaled3(1),centroidsscaled3(2),'b*');
%hold off;
%%% end dimitri stuff for plotting
%note you will need to change the output and input variables appropriately


%calculate the centroid of the whole space
cent=regionprops(true(size(Z)), Z,  'WeightedCentroid');

%get data out of struct and into an array
centroids = cat(1,cent.WeightedCentroid);

%round to nearest integer so we can use the index as a lookup value
centroidsidx(1)=round(centroids(1));
centroidsidx(2)=round(centroids(2));

%look up index and get out actual value from scale
centroidsscaled(1)=X(1,centroidsidx(1));
centroidsscaled(2)=Y(centroidsidx(2),1);

%compute centroid height
height=mean(mean(Z));

%calculate probabality of detection based on height and threshold
if height>=thresh
   prob_det=true;
else
   prob_det=false;
end

%debugging output ifblock
if debug==true
   prob_det
   height
end

end