function [x,info] = dynare_solve(func,x,jacobian_flag,varargin)
% function [x,info] = dynare_solve(func,x,jacobian_flag,varargin)
% proposes different solvers
%
% INPUTS
%    func:             name of the function to be solved
%    x:                guess values
%    jacobian_flag=1:  jacobian given by the 'func' function
%    jacobian_flag=0:  jacobian obtained numerically
%    varargin:         list of arguments following jacobian_flag
%    
% OUTPUTS
%    x:                solution
%    info=1:           the model can not be solved
%
% SPECIAL REQUIREMENTS
%    none

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

global options_

options_ = set_default_option(options_,'solve_algo',2);
info = 0;
if options_.solve_algo == 0
    if ~exist('OCTAVE_VERSION') && isempty(ver('optim'))
        error('You can''t use solve_algo=0 since you don''t have MATLAB''s Optimization Toolbox')
    end
    options=optimset('fsolve');
    options.MaxFunEvals = 50000;
    options.MaxIter = 2000;
    options.TolFun=1e-8;
    options.Display = 'iter';
    if jacobian_flag
        options.Jacobian = 'on';
    else
        options.Jacobian = 'off';
    end
    if ~exist('OCTAVE_VERSION')
        [x,fval,exitval,output] = fsolve(func,x,options,varargin{:});
    else
        % Under Octave, use a wrapper, since fsolve() does not have a 4th arg
        func2 = str2func(func);
        func = @(x) func2(x, varargin{:});
        % The Octave version of fsolve do not converge when it starts from the solution
        [fvec,fjac] = feval(func,x,varargin{:});
        if max(abs(fvec)) >= options_.solve_tolf
            [x,fval,exitval,output] = fsolve(func,x,options);
        else
            exitval = 3;
        end;
    end
    
    if exitval > 0
        info = 0;
    else
        info = 1;
    end
elseif options_.solve_algo == 1
    nn = size(x,1);
    [x,info]=solve1(func,x,1:nn,1:nn,jacobian_flag,1,varargin{:});
elseif options_.solve_algo == 2 || options_.solve_algo == 4
    nn = size(x,1) ;
    tolf = options_.solve_tolf ;

    if jacobian_flag
        [fvec,fjac] = feval(func,x,varargin{:});
    else
        fvec = feval(func,x,varargin{:});
        fjac = zeros(nn,nn) ;
    end

    i = find(~isfinite(fvec));
    
    if ~isempty(i)
        disp(['STEADY:  numerical initial values incompatible with the following' ...
              ' equations'])
        disp(i')
        error('exiting ...')
    end
    
    if max(abs(fvec)) < tolf
        return ;
    end

    if ~jacobian_flag
        fjac = zeros(nn,nn) ;
        dh = max(abs(x),options_.gstep*ones(nn,1))*eps^(1/3);
        for j = 1:nn
            xdh = x ;
            xdh(j) = xdh(j)+dh(j) ;
            fjac(:,j) = (feval(func,xdh,varargin{:}) - fvec)./dh(j) ;
        end
    end

    [j1,j2,r,s] = dmperm(fjac);
    
    if options_.debug
        disp(['DYNARE_SOLVE (solve_algo=2|4): number of blocks = ' num2str(length(r))]);
    end

    % Activate bad conditioning flag for solve_algo = 2, but not for solve_algo = 4
    bad_cond_flag = (options_.solve_algo == 2);
    
    for i=length(r)-1:-1:1
        if options_.debug
            disp(['DYNARE_SOLVE (solve_algo=2|4): solving block ' num2str(i) ', of size ' num2str(r(i+1)-r(i)) ]);
        end
        [x,info]=solve1(func,x,j1(r(i):r(i+1)-1),j2(r(i):r(i+1)-1),jacobian_flag, bad_cond_flag, varargin{:});
        if info
            return
        end
    end
    fvec = feval(func,x,varargin{:});
    if max(abs(fvec)) > tolf
        [x,info]=solve1(func,x,1:nn,1:nn,jacobian_flag, bad_cond_flag, varargin{:});
    end
elseif options_.solve_algo == 3
    if jacobian_flag
        [x,info] = csolve(func,x,func,1e-6,500,varargin{:});
    else
        [x,info] = csolve(func,x,[],1e-6,500,varargin{:});
    end
else
    error('DYNARE_SOLVE: option solve_algo must be one of [0,1,2,3,4]')
end
