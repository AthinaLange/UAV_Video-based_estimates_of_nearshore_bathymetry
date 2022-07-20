%% UAV video-based estimates of nearshore bathymetry
% 
%% Input
%   survey (structure)
%       date - YYYYMMDD
%       moplines - [579, 580, 581, 582] - should be the same naming
%       convention as timestack images and cBathy grid x
%       z - [moplines, timestack cross-shore] - should be interpolated onto
%       timestack cross-shore grid (0.1m spacing)
%
%   cBathy (structure)
%       date - YYYYMMDD
%       hover - hover number of the day (used to match with timestack name)
%       x - x coordinates from localX [alongshore coordinates, cross-shore coordinates](e.g., [775:5:-85])
%       y - y coordinates from localY (MOP numbers) [alongshore coordinates, cross-shore coordinates] (e.g., [576:0.2:588.8])
%       z - depth values (-bathy.fCombined.h), with tidal correction [alongshore coordinates, cross-shore coordinates]
%       zerr - 95per confidence interval (bathy.fCombined.hErr) [alongshore coordinates, cross-shore coordinates]
%       tide - tide level during hover (pulled from cBathy stucture array)
%
%   timestacks - located in [date_dir '/timestacks/data]
%
%% Installation
%
%       <https://github.com/Coastal-Imaging-Research-Network/Support-Routines>
%       <https://github.com/Coastal-Imaging-Research-Network/Station-Design-Toolbox>
%
%% Dependencies
%   bathy_inversion
%       pixel_res - change parameters based on camera
%       getNOAAtide
%       breakpt_calculator
%       create_composite_bathys
%       calc_errors
%
%% Output
%   
%   Video_bathy (structure)
%       date - YYYYMMDD
%       location - ie. 'Torrey' (pulled from timestack name)
%       flight - flight number
%       mop - mop number
%       x10 - [0:0.1:500]m
%       survey - z - Depth pulled from survey on MOP line
%       cbathy - 
%           z - cBathy direct output on given transect
%           zerr - cBathy hErr 
%           cbathy_hErr - cBathy with hErr > 0.5m region interpolated over
%           cbathy_gamma - cBathy with breaking region removed for given variable gamma(x)
%       tide - tide level (pulled from cbathy)
%       crests - 
%           t - time for wave tracks (sec)
%           x - cross-shore location for wave tracks (m) (0 is offshore, 500 is onshore)
%           c - phase speed of wave tracks - interp to x10
%       bp - breakpoint index of wave tracks
%       xshift - cross-shore shift index to match subaerial survey with subaqueous bathymetry
%       h_avg - interped to x10
%           lin - linear crest-tracking c = sqrt(gh)
%           nlin - nonlinear crest-tracking c = sqrt(gh(1+0.42))
%           bp - breakpoint transition crest-tracking c = sqrt(gh(1+gamma(x)))
%       gamma - transition between 0 and 0.42 for wave tracks (interped to x10)
%       gamma_mean - gamma(x) - mean of step function gamma for wave tracks
%       composite - constructed bathymetry with subaerial survey with ...
%           cbathy_hErr - cBathy with hErr < 0.5m 
%           cbathy_gamma - cBathy with breaking region removed
%           cbathy_nlin - cBathy_gamma with gamma(x) correction
%           cbathyCT - breakpoint transition crest-tracking surfzone bathymetry and cBathy offshore
%       lims - index of [1st BP valid onshore point, onshore cutoff of breaking, offshore cutoff of breaking]
%       Error - 
%           RMSE - root-mean-square error
%           Skill - skill of estimated versus observed survey bathymetry
%           Bias - bias of estimated bathymetry
%
%           insz - inner surfzone (between shoreline and end of active wave breaking - gamma(x)=0.42)
%           break - breaking region (active wave breaking region - where gamma changes)
%           sz - full surfzone region (between shoreline and wave breaking region)
%           full - [0 500] region
%           offshore - offshore region (between beginning of wave breaking and 500m - gamma(x) = 0)
%           
%           cb - default cBathy (no region removed)
%           cb_hErr - cBathy with hErr > 0.5m removed     
%           cb_gamma - cBathy with breaking region removed
%           lin - linear crest-tracking
%           nlin - nonlinear crest-tracking
%           bp - breakpoint transition crest-tracking
%           comp_hErr - composite cBathy_hErr
%           comp_gamma - composite cBathy_gamma
%           comp_nlin - composite cBathy_nlin
%           comp_CT - composite BP + cBathy
%%
%% Copyright
% (c) Lange, Athina M.Z., Fiedler, Julia W., Merrifield, Mark A., Guza, R.T. (2022)
% UAV video-based estimates of nearshore bathymetry, submitted to Coastal Engineering

%% Define path and import data
close all
clear all
% define local path
local_dir = uigetdir(pwd,'Local Directory');
cd(local_dir)

% define location of data and codes
code_dir = fullfile(local_dir, 'CODES');
addpath(genpath(code_dir))
data_dir = fullfile(local_dir, 'DATA');
addpath(genpath(data_dir))

%%
%%% check that 'survey' is in correct format
load(fullfile(data_dir, 'survey.mat'))
if ~isfield(survey, 'date')
    fprintf('Date of survey required\n')
end
if ~isfield(survey, 'moplines')
    fprintf('x-location of survey required\n')
end
if ~isfield(survey, 'z')
    fprintf('Survey z data required\n')
end

% check for empty fields in survey
for ii = 1:length(survey)
    if isempty(survey(ii).date)
        survey(ii).date = input(sprintf('Date of survey #%i required:\n Format YYYYMMDD: \n', ii))
    end
    if isempty(survey(ii).moplines)
        survey(ii).moplines = input(sprintf('MOP numbers of survey #%i required\n Format: [580, 581, 582]:\n', ii))
    end
    if isempty(survey(ii).z)
        sprintf('Survey data of survey #%i required\n Please exist and reinpurt ''survey''.\n', ii)
        break
    end
end

%%% check that 'cbathy' is in correct format
load(fullfile(data_dir, 'cbathy.mat'))
if ~isfield(cbathy, 'date')
    fprintf('Date of cBathy survey required\n')
end
if ~isfield(cbathy, 'hover')
    fprintf('Hover # of cBathy survey required\n')
end
if ~isfield(cbathy, 'x')
    fprintf('x-values of cBathy survey required\n')
end
if ~isfield(cbathy, 'y')
    fprintf('y-values of cBathy survey required\n')
end
if ~isfield(cbathy, 'z')
    fprintf('z-values of cBathy survey required\n')
end
if ~isfield(cbathy, 'zerr')
    fprintf('hErr of cBathy survey required\n')
end
if ~isfield(cbathy, 'tide')
    fprintf('Tide level of cBathy survey required\n')
end

% check for empty fields in cbathy
for ii = 1:length(cbathy)
    if isempty(cbathy(ii).date) || isempty(cbathy(ii).hover) ||isempty(cbathy(ii).x) || isempty(cbathy(ii).y) || isempty(cbathy(ii).z) || isempty(cbathy(ii).zerr) || isempty(cbathy(ii).tide)
        sprintf('cBathy data required\n Please exist and reinpurt ''cbathy''.\n')
        break
    end
end
clear ii

%% run WaveCrestDetection

% Wave Crest Detector - <https://github.com/AthinaLange/WaveCrestDetection>
% This must be executed through the terminal following the steps in
% 'Predict'. Training data included in [data_dir
% '/timestacks/training_data']. The currect version of the model is given
% in the GitHub Repository. We recommend running the prediction step on a
% GPU for time efficiency. 

% execute the Wave Crest Detector on the timestacks in [date_dir '/timestacks/data]
% collect all processed wave tracks in [data_dir '/timestacks/processed/']


%% run bathy_inversion

bathy_inversion

close all

%% bathy example plots

bathy_plots
