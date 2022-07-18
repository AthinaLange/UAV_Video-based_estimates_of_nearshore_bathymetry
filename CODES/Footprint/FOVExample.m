% example for use of function 'fov' with typical values for
% different kinds of cameras.
%
%  set the value of 'fl' for the lens you want information about

fl = 9;  % 9mm lens

% a table of several cameras, some obsolete, some not so much

cams = { ...
        { 36, 24, 'Standard 35mm film' }, ...
        { 8.8, 6.6, 'Sony XC-77 (uncal)' }, ...
        { 6.2, 4.9, 'Sony XC-75, Grasshopper. other 1/2"' }, ...
        { 8.3, 6.2, 'Grasshopper2 and other 2/3"' }, ...
        { 4.7, 3.7, 'MicroPix and other 1/3"' }, ...
        { 9.9, 8.7, '2/3" 5MP (Flea, Blackfly, etc)' }, ...
        { 11, 7, '1/1.2" 2.3MP BFLY-23S6C' }, ...
        { 7.2, 5.4, '1/1.8" 4.5um sensor' } };
        
fprintf(' Hfov  ft@  Vfov  ft@  Dfov  ft@ Sensor\n');
fprintf('      1000       1000       1000\n');
fprintf('---------------------------------------------------\n');

for ii = [1:length(cams)]
    
    cdata = cams{ii};
    
    myFOV = fov( cdata{1}, cdata{2}, fl );
    
    fprintf('%5.1f %4.0f %5.1f %4.0f %5.1f %4.0f %s\n', ...
        myFOV.horizFOV, ...
        myFOV.hdist, ...
        myFOV.vertFOV, ...
        myFOV.vdist, ...
        myFOV.diagFOV, ...
        myFOV.ddist, ...
        cdata{3} );
    
end


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

