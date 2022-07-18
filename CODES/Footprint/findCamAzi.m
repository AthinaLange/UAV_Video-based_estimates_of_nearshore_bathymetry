function [az] = findCamAzi(numCam,heading,hfov,olap)

% Function to find the azimuth angle for each camera (in degrees) in a camera
% array based on the total number of cameras, the overall heading of the
% field of view, and the horizontal field of view of the camera lens.
%
% Syntax: [az, phi] = findAzimuth(hfov,gam,N)
% 
% Inputs:  
% numCam   the number of cameras to be used
% heading  the desired heading of the center of the full camera view (deg.)
% hfov     horizontal field of view for the camera lens (deg.)
%          **it is assumed that all cameras in the array have the same hfov
%
% Output:
% az       vector of azimuths for the cameras
%

overlap = olap; % overlap between camera footprints, default is 5 degrees

azExtent = numCam*hfov - (numCam-1)*overlap; %including overlap

leftAz = heading - 90 + (180-azExtent)/2;        % in argus coords, clockwise from y=0 "north"

currentLeft = leftAz;

az = NaN(size(numCam));

for ii = 1:numCam
    az(ii) = currentLeft + hfov/2;
    currentLeft = currentLeft + hfov - overlap;
end


%%
%    Developed at Naval Research Laboratory, Stennis Space Center (2019)
%    Public Release Number: 19-1231-0203
%
%    This program is free software: you can redistribute it and/or  
%    modify it under the terms of the GNU General Public License as 
%    published by the Free Software Foundation, version 3 of the 
%    License.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.

%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see
%                                <http://www.gnu.org/licenses/>.
%