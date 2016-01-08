function [x,check] = solve1(func,x,j1,j2,jacobian_flag,bad_cond_flag,varargin)
% function [x,check] = solve1(func,x,j1,j2,jacobian_flag,bad_cond_flag,varargin)
% Solves systems of non linear equations of several variables
%
% INPUTS
%    func:            name of the function to be solved
%    x:               guess values
%    j1:              equations index for which the model is solved
%    j2:              unknown variables index
%    jacobian_flag=1: jacobian given by the 'func' function
%    jacobian_flag=0: jacobian obtained numerically
%    bad_cond_flag=1: when Jacobian is badly conditionned, use an
%                     alternative formula to Newton step
%    varargin:        list of arguments following bad_cond_flag
%    
% OUTPUTS
%    x:               results
%    check=1:         the model can not be solved
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2001-2010 Dynare Team
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

global M_ options_ fjac  

nn = length(j1);
fjac = zeros(nn,nn) ;
g = zeros(nn,1) ;

tolf = options_.solve_tolf ;
tolx = options_.solve_tolx;
tolmin = tolx ;

stpmx = 100 ;
maxit = options_.solve_maxit ;

check = 0 ;

fvec = feval(func,x,varargin{:});
fvec = fvec(j1);

i = find(~isfinite(fvec));

if ~isempty(i)
    disp(['STEADY:  numerical initial values incompatible with the following' ...
          ' equations'])
    disp(j1(i)')
end

f = 0.5*fvec'*fvec ;

if max(abs(fvec)) < tolf
    return ;
end

stpmax = stpmx*max([sqrt(x'*x);nn]) ;
first_time = 1;
for its = 1:maxit
    if jacobian_flag
        [fvec,fjac] = feval(func,x,varargin{:});
        fvec = fvec(j1);
        fjac = fjac(j1,j2);
    else
        dh = max(abs(x(j2)),options_.gstep*ones(nn,1))*eps^(1/3);
        
        for j = 1:nn
            xdh = x ;
            xdh(j2(j)) = xdh(j2(j))+dh(j) ;
            t = feval(func,xdh,varargin{:});
            fjac(:,j) = (t(j1) - fvec)./dh(j) ;
            g(j) = fvec'*fjac(:,j) ;
        end
    end

    g = (fvec'*fjac)';
    if options_.debug
        disp(['cond(fjac) ' num2str(cond(fjac))])
    end
    M_.unit_root = 0;
    if M_.unit_root
        if first_time
            first_time = 0;
            [q,r,e]=qr(fjac);
            n = sum(abs(diag(r)) < 1e-12);
            fvec = q'*fvec;
            p = e*[-r(1:end-n,1:end-n)\fvec(1:end-n);zeros(n,1)];
            disp(' ')
            disp('STEADY with unit roots:')
            disp(' ')
            if n > 0
                disp(['   The following variable(s) kept their value given in INITVAL' ...
                      ' or ENDVAL'])
                disp(char(e(:,end-n+1:end)'*M_.endo_names))
            else
                disp('   STEADY can''t find any unit root!')
            end
        else
            [q,r]=qr(fjac*e);
            fvec = q'*fvec;
            p = e*[-r(1:end-n,1:end-n)\fvec(1:end-n);zeros(n,1)];
        end     
    elseif bad_cond_flag && cond(fjac) > 1/sqrt(eps)
        fjac2=fjac'*fjac;
        p=-(fjac2+sqrt(nn*eps)*max(sum(abs(fjac2)))*eye(nn))\(fjac'*fvec);
    else
        p = -fjac\fvec ;
    end
    xold = x ;
    fold = f ;

    [x,f,fvec,check]=lnsrch1(xold,fold,g,p,stpmax,func,j1,j2,varargin{:});

    if options_.debug
        disp([its f])
        disp([xold x])
    end
    
    if check > 0
        den = max([f;0.5*nn]) ;
        if max(abs(g).*max([abs(x(j2)') ones(1,nn)])')/den < tolmin
            return
        else
            disp (' ')
            disp (['SOLVE: Iteration ' num2str(its)])
            disp (['Spurious convergence.'])
            disp (x)
            return
        end

        if max(abs(x(j2)-xold(j2))./max([abs(x(j2)') ones(1,nn)])') < tolx
            disp (' ')
            disp (['SOLVE: Iteration ' num2str(its)])
            disp (['Convergence on dX.'])
            disp (x)
            return
        end
    elseif max(abs(fvec)) < tolf
        return
    end
end

check = 1;
disp(' ')
disp('SOLVE: maxit has been reached')

% 01/14/01 MJ lnsearch is now a separate function
% 01/16/01 MJ added varargin to function evaluation
% 04/13/01 MJ added test  f < tolf !!
% 05/11/01 MJ changed tests for 'check' so as to remove 'continue' which is
%             an instruction which appears only in version 6










