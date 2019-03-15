function [location, velocity] = propagate(location, velocity, eField, timeStep)
%propagate Returns a new location vector based Newton's laws
%   propagate(location, velocity, system, timeStep)
%   Inpus:
%       location - vector marking original locations
%       velocity - vector marking speed of particles
%       timeStep - Amount of time particle is travelling before update
%   Outputs:
%       newLocation - Vector containing new locations

currentSize = size(location, 1);

%Add a row of blank zeros to get correct size for new locations
location(currentSize + 1,:) = zeros(1,size(location,2));
velocity(currentSize + 1,:) = zeros(1,size(location,2));

%Add the new location values into the new row, using old currentSize
try    
    location(currentSize + 1,:) = location(currentSize,:) + velocity(currentSize,:)*timeStep;
    velocity(currentSize + 1,:) = velocity(currentSize,:);
catch
    fprintf("busted \n");
end

%newLocation = location;

end

