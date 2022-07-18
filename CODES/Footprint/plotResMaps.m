function plotResMaps(ca, R)
%   plotResMaps(camArray, R)
%
%  plot maps of the cross-shore and alongshore resolutions for any camera
%  installation.  The cameraArray is described by vector camArray, with a
%  row for each camera with content
%    Azimuth Tilt hFOV NU NV cx cy cz
% 
%  where Azimuth, Tilt and hFOV are in radians, representing the camera
%  azimuth (compass convention, CW from y-axis North), tilt is from nadir
%  (horizontal is pi/2), and hFOV is the horizontal field of view.  Roll is
%  assumed to be 0.  NU and NV and the image size in columns and rows while
%  cx, cy and cz are the camera location coordinates.  
%  R is the resolution map structure, created by makeResMapLocal.  We graph
%  the resolution in x and y components
%  Plots will be made in the current figure

clf
subplot(121); 
imagesc(R.x, R.y, R.dcProj);        % cross-shore res.
axis xy; axis equal; axis tight;
colorbar
title('Cross-shore resolution map for demo')

hold on
az = pi/2-ca(:,1);      % convert to math convention
cmap = jet(size(ca,1));
for i = 1: size(ca,1)
    cx = ca(i,6); cy = ca(i,7);
    azLeft = tan(az(i)+ca(i,3)/2);
    azRight = tan(az(i)-ca(i,3)/2);
    if azLeft>0
        my = max(R.y);
    else
        my = min(R.y);
    end
    x = cx+(my-cy)/(azLeft+eps);
    h=plot([cx x],[cy my], '--');
    set(h,'color', cmap(i,:))
    if azRight>0
        my = max(R.y);
    else
        my = min(R.y);
    end
    x = cx+(my-cy)/(azRight+eps);
    h=plot([cx x],[cy my], '--');
    set(h,'color', cmap(i,:))
end
axis([min(R.x) max(R.x) min(R.y) max(R.y)])

subplot(122); 
imagesc(R.x, R.y, R.daProj);        % cross-shore res.
axis xy; axis equal; axis tight;
colorbar
title('Alongshore resolution map for demo')

hold on
az = pi/2-ca(:,1);      % convert to math convention
cmap = jet(size(ca,1));
for i = 1: size(ca,1)
    cx = ca(i,6); cy = ca(i,7);
    azLeft = tan(az(i)+ca(i,3)/2);
    azRight = tan(az(i)-ca(i,3)/2);
    if azLeft>0
        my = max(R.y);
    else
        my = min(R.y);
    end
    x = cx+(my-cy)/(azLeft+eps);
    h=plot([cx x],[cy my], '--');
    set(h,'color', cmap(i,:))
    if azRight>0
        my = max(R.y);
    else
        my = min(R.y);
    end
    x = cx+(my-cy)/(azRight+eps);
    h=plot([cx x],[cy my], '--');
    set(h,'color', cmap(i,:))
end
axis([min(R.x) max(R.x) min(R.y) max(R.y)])


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
