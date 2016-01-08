function [hessian_mat, gg, htol1, ihh, hh_mat0, hh1] = mr_hessian(init,x,func,hflag,htol0,varargin)
%  [hessian_mat, gg, htol1, ihh, hh_mat0, hh1] = mr_hessian(init,x,func,hflag,htol0,varargin)
%
%  numerical gradient and Hessian, with 'automatic' check of numerical
%  error 
%
% adapted from Michel Juillard original rutine hessian.m
%
%  func =  name of the function: func must give two outputs: 
%    - the log-likelihood AND the single contributions at times t=1,...,T
%    of the log-likelihood to compute outer product gradient
%  x = parameter values
%  hflag = 0, Hessian computed with outer product gradient, one point
%  increments for partial derivatives in  gradients
%  hflag = 1, 'mixed' Hessian: diagonal elements computed with numerical second order derivatives
%             with correlation structure as from outer product gradient;
%             two point evaluation of derivatives for partial derivatives
%             in gradients
%  hflag = 2, full numerical Hessian, computes second order partial derivatives
%          uses Abramowitz and Stegun (1965) formulas 25.3.24 and 25.3.27
%          p. 884.
%  htol0 = 'precision' of increment of function values for numerical
%  derivatives
%
%  varargin: other parameters of func

% Copyright (C) 2004-2011 Dynare Team
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

global options_ bayestopt_
persistent h1 htol

n=size(x,1);
if init,
    gstep_=options_.gstep;
    htol = 1.e-4;
    %h1=max(abs(x),sqrt(gstep_)*ones(n,1))*eps^(1/4);
    h1=options_.gradient_epsilon*ones(n,1);
    return,
end
func = str2func(func);
[f0, ff0]=feval(func,x,varargin{:});
h2=bayestopt_.ub-bayestopt_.lb;
hmax=bayestopt_.ub-x;
hmax=min(hmax,x-bayestopt_.lb);

h1 = min(h1,0.5.*hmax);

if htol0<htol,
    htol=htol0;
end
xh1=x;
f1=zeros(size(f0,1),n);
f_1=f1;
ff1=zeros(size(ff0));
ff_1=ff1;
ggh=zeros(size(ff0,1),n);

i=0;
while i<n,
    i=i+1;
    h10=h1(i);
    hcheck=0;
    xh1(i)=x(i)+h1(i);
    try
        [fx, ffx]=feval(func,xh1,varargin{:});
    catch
        fx=1.e8;
    end
    it=1;
    dx=(fx-f0);
    ic=0;
    
    icount = 0;
    h0=h1(i);
    while (abs(dx(it))<0.5*htol || abs(dx(it))>(3*htol)) && icount<10 && ic==0,
        %while abs(dx(it))<0.5*htol && icount< 10 && ic==0,
        icount=icount+1;
        if abs(dx(it))<0.5*htol
            if abs(dx(it)) ~= 0,
                h1(i)=min(max(1.e-10,0.3*abs(x(i))), 0.9*htol/abs(dx(it))*h1(i));
            else
                h1(i)=2.1*h1(i);
            end
            h1(i) = min(h1(i),0.5*hmax(i));
            xh1(i)=x(i)+h1(i);
            try
                [fx, ffx]=feval(func,xh1,varargin{:});
            catch
                fx=1.e8;
            end
        end
        if abs(dx(it))>(3*htol),
            h1(i)= htol/abs(dx(it))*h1(i);
            xh1(i)=x(i)+h1(i);
            try
                [fx, ffx]=feval(func,xh1,varargin{:});
            catch
                fx=1.e8;
            end
            while (fx-f0)==0,
                h1(i)= h1(i)*2;
                xh1(i)=x(i)+h1(i);
                [fx, ffx]=feval(func,xh1,varargin{:});
                ic=1;
            end
        end
        it=it+1;
        dx(it)=(fx-f0);
        h0(it)=h1(i);
        if (h1(i)<1.e-12*min(1,h2(i)) && h1(i)<0.5*hmax(i)),% || (icount==10 &&  abs(dx(it))>(3*htol)),
            ic=1;
            hcheck=1;
        end
    end
    f1(:,i)=fx;
    if any(isnan(ffx)),
        ff1=ones(size(ff0)).*fx/length(ff0);
    else
        ff1=ffx;
    end
    xh1(i)=x(i)-h1(i);
    [fx, ffx]=feval(func,xh1,varargin{:});
    f_1(:,i)=fx;
    if any(isnan(ffx)),
        ff_1=ones(size(ff0)).*fx/length(ff0);
    else
        ff_1=ffx;
    end
    ggh(:,i)=(ff1-ff_1)./(2.*h1(i));
    xh1(i)=x(i);
    if hcheck & htol<1,
        htol=min(1,max(min(abs(dx))*2,htol*10));
        h1(i)=h10;
        i=0;
    end
    save hess.mat
end

h_1=h1;
xh1=x;
xh_1=xh1;

gg=(f1'-f_1')./(2.*h1);

if hflag==2,
    gg=(f1'-f_1')./(2.*h1);
    hessian_mat = zeros(size(f0,1),n*n);
    for i=1:n
        if i > 1
            k=[i:n:n*(i-1)];
            hessian_mat(:,(i-1)*n+1:(i-1)*n+i-1)=hessian_mat(:,k);
        end 
        hessian_mat(:,(i-1)*n+i)=(f1(:,i)+f_1(:,i)-2*f0)./(h1(i)*h_1(i));
        temp=f1+f_1-f0*ones(1,n);
        for j=i+1:n
            xh1(i)=x(i)+h1(i);
            xh1(j)=x(j)+h_1(j);
            xh_1(i)=x(i)-h1(i);
            xh_1(j)=x(j)-h_1(j);
            temp1 = feval(func,xh1,varargin{:});
            
            temp2 = feval(func,xh_1,varargin{:});
            
            hessian_mat(:,(i-1)*n+j)=-(-temp1 -temp2+temp(:,i)+temp(:,j))./(2*h1(i)*h_1(j));
            xh1(i)=x(i);
            xh1(j)=x(j);
            xh_1(i)=x(i);
            xh_1(j)=x(j);
            j=j+1;
            save hess.mat
        end
        i=i+1;
    end
    
elseif hflag==1,
    hessian_mat = zeros(size(f0,1),n*n);
    for i=1:n,
        dum = (f1(:,i)+f_1(:,i)-2*f0)./(h1(i)*h_1(i));
        if dum>eps,
            hessian_mat(:,(i-1)*n+i)=dum;
        else
            hessian_mat(:,(i-1)*n+i)=max(eps, gg(i)^2);
        end                        
    end
    %hessian_mat2=hh_mat(:)';
end

gga=ggh.*kron(ones(size(ff1)),2.*h1');  % re-scaled gradient
hh_mat=gga'*gga;  % rescaled outer product hessian 
hh_mat0=ggh'*ggh;  % outer product hessian
A=diag(2.*h1);  % rescaling matrix
% igg=inv(hh_mat);  % inverted rescaled outer product hessian
ihh=A'*(hh_mat\A);  % inverted outer product hessian
if hflag>0 & min(eig(reshape(hessian_mat,n,n)))>0,
    hh0 = A*reshape(hessian_mat,n,n)*A';  %rescaled second order derivatives
    hh = reshape(hessian_mat,n,n);  %rescaled second order derivatives
    sd0=sqrt(diag(hh0));   %rescaled 'standard errors' using second order derivatives
    sd=sqrt(diag(hh_mat));  %rescaled 'standard errors' using outer product
    hh_mat=hh_mat./(sd*sd').*(sd0*sd0');  %rescaled inverse outer product with 'true' std's
    igg=inv(hh_mat);   % rescaled outer product hessian with 'true' std's
    ihh=A'*(hh_mat\A);  % inverted outer product hessian
    hh_mat0=inv(A)'*hh_mat*inv(A);  % outer product hessian with 'true' std's
    sd=sqrt(diag(ihh));   %standard errors
    sdh=sqrt(1./diag(hh));   %diagonal standard errors
    for j=1:length(sd),
        sd0(j,1)=min(bayestopt_.p2(j), sd(j));  %prior std
        sd0(j,1)=10^(0.5*(log10(sd0(j,1))+log10(sdh(j,1))));
        %sd0(j,1)=0.5*(sd0(j,1)+sdh(j,1));
    end
    ihh=ihh./(sd*sd').*(sd0*sd0');  %inverse outer product with modified std's
    igg=inv(A)'*ihh*inv(A);  % inverted rescaled outer product hessian with modified std's
    hh_mat=inv(igg);   % outer product rescaled hessian with modified std's
    hh_mat0=inv(A)'*hh_mat*inv(A);  % outer product hessian with modified std's
                                    %     sd0=sqrt(1./diag(hh0));   %rescaled 'standard errors' using second order derivatives
                                    %     sd=sqrt(diag(igg));  %rescaled 'standard errors' using outer product
                                    %     igg=igg./(sd*sd').*(sd0*sd0');  %rescaled inverse outer product with 'true' std's
                                    %     hh_mat=inv(igg);   % rescaled outer product hessian with 'true' std's
                                    %     ihh=A'*igg*A;  % inverted outer product hessian
                                    %     hh_mat0=inv(A)'*hh_mat*inv(A);  % outer product hessian with 'true' std's
end
if hflag<2,
    hessian_mat=hh_mat0(:);
end

if any(isnan(hessian_mat)),
    hh_mat0=eye(length(hh_mat0));
    ihh=hh_mat0;
    hessian_mat=hh_mat0(:);    
end
hh1=h1;
htol1=htol;
save hess.mat
% 11/25/03 SA Created from Hessian_sparse (removed sparse)


