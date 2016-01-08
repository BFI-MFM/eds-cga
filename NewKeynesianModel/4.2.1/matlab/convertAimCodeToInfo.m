function [info] = convertAimCodeToInfo(aimCode)
% function [info] = convertAimCodeToInfo(aimCode)
% Returns an appropriate code for print_info
%
% INPUTS
%   aimCode     [integer]    code returned by AIM
%      (aimCode==1)  e='Aim: unique solution.';
%      (aimCode==2)  e='Aim: roots not correctly computed by real_schur.';
%      (aimCode==3)  e='Aim: too many big roots.';
%      (aimCode==35) e='Aim: too many big roots, and q(:,right) is singular.';
%      (aimCode==4)  e='Aim: too few big roots.';
%      (aimCode==45) e='Aim: too few big roots, and q(:,right) is singular.';
%      (aimCode==5)  e='Aim: q(:,right) is singular.';
%      (aimCode==61) e='Aim: too many exact shiftrights.';
%      (aimCode==62) e='Aim: too many numeric shiftrights.';
%
% OUTPUTS
%   info        [integer]    Code to be used to print error in print_info.m

% Copyright (C) 2011 Dynare Team
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

switch aimCode
  case 1
    info = 0; % no problem encountered
  case 2
    info = 102;
  case 3
    info = 103;
  case 35
    info = 135;
  case 4
    info = 104;
  case 45
    info = 145;
  case 5
    info = 105;
  case 61
    info = 161;
  case 62
    info = 162;
  otherwise
    info = 1;
end