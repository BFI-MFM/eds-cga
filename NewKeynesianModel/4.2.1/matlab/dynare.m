function dynare(fname, varargin)
%       This command runs dynare with specified model file in argument
%       Filename.
%       The name of model file begins with an alphabetic character, 
%       and has a filename extension of .mod or .dyn.
%       When extension is omitted, a model file with .mod extension
%       is processed.
%
% INPUTS
%   fname:      file name
%   varargin:   list of arguments following fname
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

if strcmpi(fname,'help')
    disp(' ')
    disp(['This is dynare version ' dynare_version() '.'])
    disp(' ')
    disp('USAGE: dynare FILENAME[.mod,.dyn] [OPTIONS]')
    disp(' ')
    disp('dynare executes instruction included in FILENAME.mod.')
    disp(' ')
    disp('OPTIONS:')
    disp(' o noclearall:  By default, dynare  will issue a clear all command to Matlab or Octave,')
    disp('                thereby deleting all workspace variables; this options instructs dynare') 
    disp('                not to clear the workspace.')
    disp(' o debug:       Instructs the preprocessor to write some debugging informations about the') 
    disp('                scanning and parsing of the .mod file.')
    disp(' o notmpterms:  Do not include temporary terms in the generated static and dynamic files')
    disp(' o savemacro[=MACROFILE]: Instructs dynare  to save the intermediary file which is obtained')
    disp('                after macro-processing.')
    disp('                By default, saved output will go in FILENAME-macroexp.mod (where FILENAME is')
    disp('                taken from the mod file), but you can specify another filename.')
    disp(' ')
    return
end

warning_config()

if exist('OCTAVE_VERSION')
    if octave_ver_less_than('3.0.0')
        warning('This version of Dynare has only been tested on Octave 3.0.0 and above. Since your Octave version is older than that, Dynare may fail to run, or give unexpected results. Consider upgrading your Octave installation.');
    end
else
    if matlab_ver_less_than('7.0')
        warning('This version of Dynare has only been tested on MATLAB 7.0 (R14) and above. Since your MATLAB version is older than that, Dynare may fail to run, or give unexpected results. Consider upgrading your MATLAB installation, or switch to Octave.');
    end
end

% disable output paging (it is on by default on Octave)
more off

% sets default format for save() command
if exist('OCTAVE_VERSION')
    default_save_options('-mat')
end

% detect if MEX files are present; if not, use alternative M-files
dynareroot = dynare_config;

if nargin < 1
    error('DYNARE: you must provide the name of the MOD file in argument')
end

if ~ischar(fname)
    error('DYNARE: argument of dynare must be a text string')
end

% Testing if file have extension
% If no extension default .mod is added
if isempty(strfind(fname,'.'))
    fname1 = [fname '.dyn'];
    d = dir(fname1);
    if length(d) == 0
        fname1 = [fname '.mod'];
    end
    fname = fname1;
    % Checking file extension
else
    if ~strcmp(upper(fname(size(fname,2)-3:size(fname,2))),'.MOD') ...
            && ~strcmp(upper(fname(size(fname,2)-3:size(fname,2))),'.DYN')
        error('DYNARE: argument must be a filename with .mod or .dyn extension')
    end;
end;
d = dir(fname);
if length(d) == 0
    error(['DYNARE: can''t open ' fname])
end

command = ['"' dynareroot 'dynare_m" ' fname] ;
for i=2:nargin
    command = [command ' ' varargin{i-1}];
end

% Workaround for bug in Octave >= 3.2
% See http://bugs.debian.org/cgi-bin/bugreport.cgi?bug=550823
if exist('OCTAVE_VERSION') && ~octave_ver_less_than('3.2.0')
    sleep(2)
end

[status, result] = system(command);
disp(result)
if status
    % Should not use "error(result)" since message will be truncated if too long
    error('DYNARE: preprocessing failed')
end

if ~ isempty(find(abs(fname) == 46))
    fname = fname(:,1:find(abs(fname) == 46)-1) ;
end
evalin('base',fname) ;
