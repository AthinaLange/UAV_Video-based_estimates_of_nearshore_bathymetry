% Argus Design Demo.
% 
% This routine is an example of how you can use resolution maps to design a
% possible new Argus station.  As a user, you need to enter information on
% the set of desired cameras and the region of interest.  Cameras are
% defined by their number of pixels, NV by NU, their horizontal azimuth
% (viewing) angles and vertical tilt, and the selected camera position.
% All locations must be in local coordinates (whatever you chose).  The
% routine will then compute an R structure with resolutions in the range
% and cross-range directions, but more importantly in the x and y
% directions.  See help for makeResMapLocal.  Finally the resolution maps
% are displayed.
%
% Note that a small overlap is built in using variable overlap, a
% percentage.

clear all;

% Define the domain in terms of xmin xmax ymin ymax dx dy zLevel.
% dx and dy are only for the presentation of the res map and have no
% lasting design value.  The last number is the mapped vertical level,
% usually sea level.  You also define the left edge azimuth.

XYZ = [0 700 -650 1250 10 10 0];
leftAz = 0;        % in argus coords, clockwise from y=0 "north"
leftAz = leftAz/180*pi; % convert to radians

% Now choose the camera location 
cx = 0;             % wherever you want to put the camera in your
cy = 0;             % coordinates.
cz = 120/3.2808;    % cam height converted from ft to m.

% Now define the cameras including the number of pixels NU and NV and the 
% horizontal fields of view of each N camera atnd 
NU = 2448;  NV = 2048;          % horizontal and vertical num pixels

% list the horizontal field of view and the top of view limits.  Listed in
% degrees.  topOfView is taken to just include the horizon.
% These are examples from a real installation, but they may not be
% appropriate to lenses that you have available.
hfov(1) = 34.4; topOfView(1) = 91;   % listed left to right looking seaward.
hfov(2) = 43.2; topOfView(2) = 91; 
hfov(3) = 43.2; topOfView(3) = 91;
hfov(4) = 43.2; topOfView(4) = 91;
hfov(5) = 34.4; topOfView(5) = 91;
overlap = 0.05;         % percentage of overlap between adjacent cam views
azExtent = sum(hfov);   % excluding overlap

% check if you span enough azimuth
disp(['Sum of azimuths is ' num2str(azExtent,'%.1f') ' degrees.'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%% End of User Input %%%%%%%%%%%%%%%%%%%%%%

% Now go through the cameras from left to right.  Build an array of
% azimuths, tilts, field of views, etc, that will be fed to the routine to
% find the resolution map.

currentLeft = leftAz;
ca = [];
disp( ['Cam  Azimuth     Tilt      FOV      VFOV']) ;
for i = 1:length(hfov)
    fov = deg2rad(hfov(i));     % fov in radians
    vfov = fov * NV/NU;         % assumes roughly square pixels
    c = [ currentLeft+((1-overlap)*fov/2), deg2rad(topOfView(i))-vfov/2, fov, NU, NV, cx, cy, cz];
    disp( [ sprintf( '%3d', i ) ...
            sprintf( '%9.2f', rad2deg(c(1)) ) ...
            sprintf( '%9.2f', rad2deg(c(2)) ) ...
            sprintf( '%9.2f', rad2deg(c(3)) ) ...
            sprintf( '%9.2f', rad2deg(vfov) )] );
    currentLeft = currentLeft + .9 * fov;
    ca = [ca; c]; 
end

% make the resolution map structure
R = makeResMapLocal( ca, XYZ );

% Plot the resolution maps, showing camera boundaries
figure(1); clf
plotResMaps(ca, R)


%   Copyright (C) 2017  Coastal Imaging Research Network
%                       and Oregon State University

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

% CIRN: https://coastal-imaging-research-network.github.io/
% CIL:  http://cil-www.coas.oregonstate.edu
