function [az,tilt,roll] = R2Angles(R)

% function [az,tilt,roll] = R2Angles(R)
%
% Opposite of angles2R
%

tilt = acos(R(3,3)*-1);
roll = asin(R(1,3) / sin(tilt));
az1 = acos(R(3,2) / sin(tilt));
az2 = asin(R(3,1) / sin(tilt));

if (sign(az1) ~= sign(az2))
  az = az1*-1;
else
  az = az1;
end

az = mod(az, 2*pi);

%

%
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
%
%key support routines Argus CIL CIRN
%

