function dynareroot = dynare_config(path_to_dynare)
%function dynareroot = dynare_config(path_to_dynare)
% This function tests the existence of valid mex files (for qz
% decomposition, solution to sylvester equation and kronecker
% products...) and, if needed, add paths to the matlab versions
% of these routines.
% Also adds other directories to the path.
%
% INPUTS
%   none
%             
% OUTPUTS
%   none
%        
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2001-2011 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

if nargin
    addpath(path_to_dynare);
end
dynareroot = strrep(which('dynare'),'dynare.m','');

addpath([dynareroot '/distributions/'])
addpath([dynareroot '/kalman/'])
addpath([dynareroot '/kalman/likelihood'])
addpath([dynareroot '/AIM/'])
addpath([dynareroot '/partial_information/'])
addpath([dynareroot '/parallel/'])

% For functions that exist only under some Octave versions
% or some MATLAB versions, and for which we provide some replacement functions

if ~exist('OCTAVE_VERSION')
    % Replacements for rows() and columns() (inexistent under MATLAB)
    addpath([dynareroot '/missing/rows_columns'])
    % Replacement for vec() (inexistent under MATLAB)
    addpath([dynareroot '/missing/vec'])
    if isempty(ver('stats'))
        % Replacements for functions of the stats toolbox
        addpath([dynareroot '/missing/stats/'])
    end
end

% ordeig() was introducted in MATLAB 7.0.1, and doesn't exist in Octave
if exist('OCTAVE_VERSION') || matlab_ver_less_than('7.0.1')
    addpath([dynareroot '/missing/ordeig'])
end

% rcond() and bicgstable() were introduced in Octave 3.2.0
if exist('OCTAVE_VERSION') && octave_ver_less_than('3.2.0')
    addpath([dynareroot '/missing/rcond'])
    addpath([dynareroot '/missing/bicgstab'])
end

% bsxfun is missing in old versions of MATLAB (and exists in Octave)
if ~exist('OCTAVE_VERSION') && matlab_ver_less_than('7.4')
    addpath([dynareroot '/missing/bsxfun'])
end


% Add path to MEX files
if exist('OCTAVE_VERSION')
    addpath([dynareroot '../mex/octave/']);
else
    % Add win32 specific paths for Dynare Windows package
    if strcmp(computer, 'PCWIN')
        if matlab_ver_less_than('7.5')
            mexpath = [dynareroot '../mex/matlab/win32-7.0-7.4'];
            if exist(mexpath, 'dir')
                addpath(mexpath)
            end
        else
            mexpath = [dynareroot '../mex/matlab/win32-7.5-7.12'];
            if exist(mexpath, 'dir')
                addpath(mexpath)
            end
        end
    end

    % Add win64 specific paths for Dynare Windows package
    if strcmp(computer, 'PCWIN64')
        if matlab_ver_less_than('7.3')
            mexpath = [dynareroot '../mex/matlab/win64-7.2'];
            if exist(mexpath, 'dir')
                addpath(mexpath)
            end
        elseif matlab_ver_less_than('7.5')
            mexpath = [dynareroot '../mex/matlab/win64-7.3-7.4'];
            if exist(mexpath, 'dir')
                addpath(mexpath)
            end
        elseif matlab_ver_less_than('7.8')
            mexpath = [dynareroot '../mex/matlab/win64-7.5-7.7'];
            if exist(mexpath, 'dir')
                addpath(mexpath)
            end
        else
            mexpath = [dynareroot '../mex/matlab/win64-7.8-7.12'];
            if exist(mexpath, 'dir')
                addpath(mexpath)
            end
        end
    end

    % Add generic MATLAB path (with higher priority than the previous ones)
    addpath([dynareroot '../mex/matlab/']);
end

%% Set mex routine names
mex_status = cell(1,3);
mex_status(1,1) = {'mjdgges'};
mex_status(1,2) = {'qz'};
mex_status(1,3) = {'Generalized QZ'};
mex_status(2,1) = {'gensylv'};
mex_status(2,2) = {'gensylv'};
mex_status(2,3) = {'Sylvester equation solution'};
mex_status(3,1) = {'A_times_B_kronecker_C'};
mex_status(3,2) = {'kronecker'};
mex_status(3,3) = {'Kronecker products'};
mex_status(4,1) = {'sparse_hessian_times_B_kronecker_C'};
mex_status(4,2) = {'kronecker'};
mex_status(4,3) = {'Sparse kronecker products'};
number_of_mex_files = size(mex_status,1);
%% Remove some directories from matlab's path. This is necessary if the user has
%% added dynare_v4/matlab with the subfolders. Matlab has to ignore these
%% subfolders if valid mex files exist.
matlab_path = path;
for i=1:number_of_mex_files
    test = strfind(matlab_path,[dynareroot mex_status{i,2}]);
    action = length(test);
    if action
        rmpath([dynareroot mex_status{i,2}]);
        matlab_path = path;
    end
end
%% Test if valid mex files are available, if a mex file is not available
%% a matlab version of the routine is included in the path.
disp(' ')
disp('Configuring Dynare ...')

for i=1:number_of_mex_files
    test = (exist(mex_status{i,1},'file') == 3);
    if ~test
        addpath([dynareroot mex_status{i,2}]);
        message = '[m]   ';
    else
        message = '[mex] ';
    end
    disp([ message mex_status{i,3} '.' ])
end

% Test if bytecode DLL is present
if exist('bytecode', 'file') == 3
    message = '[mex] ';
else
    message = '[no]  ';
end
disp([ message 'Bytecode evaluation.' ])

% Test if k-order perturbation DLL is present
if exist('k_order_perturbation', 'file') == 3
    message = '[mex] ';
else
    message = '[no]  ';
end
disp([ message 'k-order perturbation solver.' ])

% Test if dynare_simul_ DLL is present
if exist('dynare_simul_', 'file') == 3
    message = '[mex] ';
else
    message = '[no]  ';
end
disp([ message 'k-order solution simulation.' ])

disp(' ')
