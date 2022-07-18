function myFOV = fov( width, height, fl )

% function FOV = FOV( width, height, fl )
%
%  calculate the field of view (FOV) for a lens/camera pair 
%  using the width and height of the sensor and the focal length
%  of the lens. This will be a "book" value; normal cameras and 
%  lenses have physical differences that will determine final values.
%
%  height:   height of the sensor in mm
%  width:    width of the sensor in mm
%  fl:       focal length of the lens in mm
%
% FOV is a struct containing the following fields:
%  horizFOV:  horizontal FOV in degrees
%  vertFOV:   vertical FOV in degrees
%  diagFOV:   diagonal FOV in degrees
%  hdist:     horizontal view in feet at 1000 feet. Or meters at 1000m. 
%  vdist:     vertical view in feet at 1000 feet. Or microns at 1000 microns.
%  ddist:     diagonal view in feet at 1000 feet. Or angstroms at ...
%
% Note: width and height are not the "sensor size" that is often used
% in the camera literature, it is the physical sensor size. You can usually
% determine this by searching out the sensor data sheet and looking for
% the pixel width and heigh. Multiply that by the number of pixels.
%

	% calculate the diagonal FOV, too.
	diag = sqrt(height*height+width*width);

	myFOV.horizFOV = dofov( width, fl );
	myFOV.vertFOV = dofov( height, fl );
	myFOV.diagFOV = dofov( diag, fl );

	myFOV.hdist = fovft( myFOV.horizFOV );
	myFOV.vdist = fovft( myFOV.vertFOV );
	myFOV.ddist = fovft( myFOV.diagFOV );

function f = dofov( l, fl )

%  l = length (width or height)
%  fl = focal length

	% from "Applied Photographic Optics" by Ray, p.112
	f = 2 * atan2( l, 2*fl );

	% convert to degrees
	f = f * 180 / pi; 

    return;
    
function ft = fovft( fov )

% calcualte ft@1000ft for a given fov

	% convert fov to radians, use only half of FOV to calculate leg
	fovr = (fov/2) * pi / 180;

	% calculate the tangent (perl gots no tan())
	fovinft = sin( fovr ) / cos( fovr );

	% and double to make full FOV, and at 1000 ft
	ft = 2 * 1000 * fovinft;

	return;


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

