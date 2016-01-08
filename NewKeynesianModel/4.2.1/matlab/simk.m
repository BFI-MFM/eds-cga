function simk
% function simk
% performs deterministic simulations with lead or lag on more than one
% period
%
% Currently used only for purely forward models.
%
% INPUTS
%   ...
% OUTPUTS
%   ...
% ALGORITHM
%   Laffargue, Boucekkine, Juillard (LBJ)
%   see Juillard (1996) Dynare: A program for the resolution and
%   simulation of dynamic models with forward variables through the use
%   of a relaxation algorithm. CEPREMAP. Couverture Orange. 9602.
%
% SPECIAL REQUIREMENTS
%   None.
%  

% Copyright (C) 1996-2010 Dynare Team
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

global M_ options_ oo_

nk = M_.maximum_endo_lag + M_.maximum_endo_lead + 1 ;
ny = size(M_.lead_lag_incidence,2) ;
icc1 = M_.lead_lag_incidence(nk,:) > 0;

for i = 1:M_.maximum_lead -1
    icc1 = [M_.lead_lag_incidence(nk-i,:) | icc1(1,:); icc1] ;
end

icc1 = find(icc1') ;
iy = M_.lead_lag_incidence > 0 ;
isc = cumsum(sum(iy',1))' ;
iyr0 = find(M_.lead_lag_incidence') ;
ncc1 = size(icc1,1) ;
ncc = ncc1 + 1 ;
ncs = size(iyr0,1) ;

ky = zeros(ny,nk) ;            % indices of variables at each lead or lag
lky = zeros(nk,1) ;
for i = 1:nk
    j = find(M_.lead_lag_incidence(i,:))' ;
    if isempty(j)
        lky(i) = 0;
    else
        lky(i) = size(j,1) ;
        ky(1:lky(i),i) = j ;
    end
end

jwc = find(iy(2:M_.maximum_endo_lead+1,:)') ; % indices of columns for
                                              % triangularization
                                              % as many rows as lags in model

if isempty(jwc)
    jwc = 0 ;
    ljwc = 0 ;
    temp = icc1 ;
else
    ljwc = size(jwc,1) ;          % length of each row in jwc
    temp = union(jwc,icc1,'rows') ;      % prepares next iteration
end

j1 = ky(1:lky(1),1) ;
lj1 = lky(1) ;

for i = 2:M_.maximum_endo_lag
    [j1,lj1] = ffill(j1,lj1,selif(temp+(i-1)*ny,temp <= ny)) ;
    if M_.maximum_endo_lead == 1
        if lky(i+M_.maximum_endo_lead) > 0
            [jwc,ljwc] = ffill(jwc,ljwc, ky(1:lky(i+M_.maximum_endo_lead),i+M_.maximum_endo_lead)+(M_.maximum_endo_lead-1)*ny) ;
            if ljwc(i) == 0
                temp = icc1;
            else
                temp = union(jwc(1:ljwc(i),i),icc1,'rows') ;
            end
        else
            [jwc,ljwc] = ffill(jwc,ljwc,[]) ;
            temp = icc1 ;
        end
    else
        temp = temp(lj1(i)+1:size(temp,1),:) - ny ;
        if lky(i+M_.maximum_endo_lead) > 0
            [jwc,ljwc] = ffill(jwc,ljwc,[temp;ky(1:lky(i+M_.maximum_endo_lead),i+M_.maximum_endo_lead)+(M_.maximum_endo_lead-1)*ny]);
        else
            [jwc,ljwc] = ffill(jwc,ljwc,temp) ;
        end
        temp = union(jwc(1:ljwc(i),i),icc1,'rows') ;
    end
end

[j1,lj1] = ffill(j1,lj1,selif(temp+M_.maximum_endo_lag*ny, temp <= ny)) ;
ltemp = zeros(M_.maximum_endo_lag,1) ;
jwc1 = zeros(ncc1,M_.maximum_endo_lag) ;

for i = 1:M_.maximum_endo_lag
    temp = union(jwc(1:ljwc(i),i),icc1,'rows') ;
    ltemp(i) = size(temp,1) ;
    if ljwc(i) > 0
        jwc(1:ljwc(i),i) = indnv(jwc(1:ljwc(i),i),temp) ;
    end
    jwc1(:,i) = indnv(icc1,temp) ;
end

h1 = clock ;

disp (['-----------------------------------------------------']) ;
disp ('MODEL SIMULATION') ;
fprintf ('\n') ;

for iter = 1:options_.maxit_
    h2 = clock ;
    oo_.endo_simul = oo_.endo_simul(:);
    err_f = 0;
    
    fid = fopen([M_.fname '.swp'],'w+') ;

    it_ = 1+M_.maximum_lag ;
    ic = [1:ny] ;
    iyr = iyr0 ;
    i = M_.maximum_endo_lag+1 ;
    while (i>1) & (it_<=options_.periods+M_.maximum_endo_lag)
        h3 = clock ;
        [d1,jacobian] = feval([M_.fname '_dynamic'],oo_.endo_simul(iyr),oo_.exo_simul, M_.params, it_);
        d1 = -d1 ;
        err_f = max(err_f,max(abs(d1)));
        if lky(i) ~= 0
            j1i = ky(1:lky(i),i) ;
            w0 = jacobian(:,isc(i-1)+1:isc(i)) ;
        else
            w0 = [];
        end
        ttemp = iy(i+1:i+M_.maximum_endo_lead,:)' ;
        jwci = find(ttemp) ;
        if ~ isempty(jwci)
            w = jacobian(:,isc(i)+1:isc(i+M_.maximum_endo_lead)) ;
        end
        j = i ;
        while j <= M_.maximum_endo_lag
            if ~isempty(w0)

                ofs = ((it_-M_.maximum_lag-M_.maximum_endo_lag+j-2)*ny)*ncc*8 ;
                junk = fseek(fid,ofs,-1) ;
                c = fread(fid,[ncc,ny],'float64')';

                if isempty(jwci)
                    w = -w0*c(j1i,1:ncc1) ;
                    jwci = icc1 ;
                else
                    iz = union(jwci,icc1,'rows') ;
                    ix = indnv(jwci,iz) ;
                    iy__ = indnv(icc1,iz) ;
                    temp = zeros(size(w,1),size(iz,1)) ;
                    temp(:,ix) = w;
                    temp(:,iy__) = temp(:,iy__)-w0*c(j1i,1:ncc1) ;
                    w = temp ;
                    jwci = iz ;
                    clear temp iz ix iy__ ;
                end
                d1 = d1-w0*c(j1i,ncc) ;
                clear c ;
            end
            j = j + 1 ;
            if isempty(jwci)
                j1i = [];
                if lky(j+M_.maximum_endo_lead) ~= 0
                    jwci = ky(1:lky(j+M_.maximum_endo_lead),j+M_.maximum_endo_lead) + (M_.maximum_endo_lead-1)*ny ;
                    w = jacobian(:,isc(j+M_.maximum_endo_lead-1)+1:isc(j+M_.maximum_endo_lead)) ;
                else
                    jwci = [] ;
                end
            else
                j1i = selif(jwci,jwci<(ny+1)) ;
                w0 = w(:,1:size(j1i,1)) ;
                if size(jwci,1) == size(j1i,1)
                    if lky(j+M_.maximum_endo_lead) ~= 0
                        jwci = ky(1:lky(j+M_.maximum_endo_lead),j+M_.maximum_endo_lead)+(M_.maximum_endo_lead-1)*ny ;
                        w = jacobian(:,isc(j+M_.maximum_endo_lead-1)+1:isc(j+M_.maximum_endo_lead)) ;
                    else
                        jwci = [] ;
                    end
                else
                    jwci = jwci(size(j1i,1)+1:size(jwci,1),:)-ny ;
                    w = w(:,size(j1i,1)+1:size(w,2)) ; 
                    if lky(j+M_.maximum_endo_lead) ~= 0
                        jwci = [ jwci; ky(1: lky(j+M_.maximum_endo_lead),j+M_.maximum_endo_lead)+(M_.maximum_endo_lead-1)*ny] ;
                        w = [w jacobian(:,isc(j+M_.maximum_endo_lead-1)+1:isc(j+M_.maximum_endo_lead))] ;
                        %         else
                        %           jwci = [] ;
                    end
                end
            end
        end
        jwci = [indnv(jwci,icc1);ncc] ;
        w = [w d1] ;
        c = zeros(ny,ncc) ;
        c(:,jwci) = w0\w ;
        clear w w0 ;

        junk = fseek(fid,0,1) ;
        fwrite(fid,c','float64') ;
        clear c ;

        it_ = it_ + 1;
        ic = ic + ny ;
        iyr = iyr + ny ;
        i = i - 1 ;
    end
    icr0 = (it_-M_.maximum_lag-M_.maximum_endo_lag -1)*ny ;
    while it_ <= options_.periods+M_.maximum_lag
        [d1,jacobian] = feval([M_.fname '_dynamic'],oo_.endo_simul(iyr),oo_.exo_simul, M_.params, it_);
        d1 = -d1 ;
        err_f = max(err_f,max(abs(d1)));
        w0 = jacobian(:,1:isc(1)) ;
        w = jacobian(:,isc(1)+1:isc(1+M_.maximum_endo_lead)) ;
        j = 1 ;
        while j <= M_.maximum_endo_lag
            icr = j1(1:lj1(j),j)-(j-1)*ny ;

            ofs = ((icr0+(j-1)*ny+1)-1)*ncc*8 ;
            junk = fseek(fid,ofs,-1) ;
            c = fread(fid,[ncc,ny],'float64')' ;

            temp = zeros(ny,ltemp(j)) ;
            if ljwc(j) > 0
                temp(:,jwc(1:ljwc(j),j)) = w ;
            end
            temp(:,jwc1(:,j))=temp(:,jwc1(:,j))-w0*c(icr,1:ncc1) ;
            w = temp ;
            clear temp ;
            d1 = d1-w0*c(icr,ncc) ;
            clear c ;
            j = j + 1 ;
            w0 = w(:,1:lj1(j)) ;
            if M_.maximum_endo_lead == 1
                w = jacobian(:,isc(j+M_.maximum_endo_lead-1)+1:isc(j+M_.maximum_endo_lead)) ;
            else
                w = w(:,lj1(j)+1:size(w,2)) ;

                if lky(j+M_.maximum_endo_lead) > 0
                    w = [w jacobian(:,isc(j+M_.maximum_endo_lead-1)+1:isc(j+M_.maximum_endo_lead))] ;
                end
            end
        end
        c = w0\[w d1] ;
        d1 = [] ;
        clear w w0 ;
        junk = fseek(fid,0,1) ;
        fwrite(fid,c','float64') ;
        clear c ;
        it_ = it_ + 1 ;
        ic = ic + ny ;
        iyr = iyr + ny ;
        icr0 = icr0 + ny ;
    end
    if options_.terminal_condition == 1

        ofs = (((it_-M_.maximum_lag-2)*ny+1)-1)*ncc*8 ;
        junk = fseek(fid,ofs,-1) ;
        c = fread(fid,[ncc,ny],'float64')';

        for i = 1:M_.maximum_endo_lead
            w = tril(triu(ones(ny,ny+ncc1))) ;
            w(:,jwc1(:,M_.maximum_endo_lag)) = w(:,jwc1(:,M_.maximum_endo_lag))+c(:,1:ncc1) ;
            c = [w(:,ny+1:size(w,2))' c(:,ncc)]/w(:,1:ny) ;

            junk = fseek(fid,0,1) ;
            fwrite(fid,c','float64') ;

            it_ = it_+1 ;
            ic = ic + ny ;

        end
    end
    oo_.endo_simul = reshape(oo_.endo_simul,ny,options_.periods+M_.maximum_lag+M_.maximum_endo_lead) ;
    if options_.terminal_condition == 1
        hbacsup = clock ;
        c = bksupk(ny,fid,ncc,icc1) ;
        hbacsup = etime(clock,hbacsup) ;
        c = reshape(c,ny,options_.periods+M_.maximum_endo_lead)' ;
        y(:,1+M_.maximum_endo_lag:(options_.periods+M_.maximum_endo_lead+M_.maximum_endo_lag)) = y(:,1+M_.maximum_endo_lag:(options_.periods+M_.maximum_endo_lead+M_.maximum_endo_lag))+options_.slowc*c' ;
    else
        hbacsup = clock ;
        c = bksupk(ny,fid,ncc,icc1) ;
        hbacsup = etime(clock,hbacsup) ;
        c = reshape(c,ny,options_.periods)' ;
        oo_.endo_simul(:,1+M_.maximum_endo_lag:(options_.periods+M_.maximum_endo_lag)) = oo_.endo_simul(:,1+M_.maximum_endo_lag:(options_.periods+M_.maximum_endo_lag))+options_.slowc*c' ;
    end

    fclose(fid) ;

    h2 = etime(clock,h2) ;
    [junk,i1] = max(abs(c));
    [junk,i2] = max(junk);
    disp(['variable ' M_.endo_names(i2,:) ' period ' num2str(i1(i2))])
    err = max(max(abs(c./options_.scalv'))) ;
    disp ([num2str(iter) '-     err = ' num2str(err)]) ;
    disp (['err_f = ' num2str(err_f)])
    disp (['    Time of this iteration  : ' num2str(h2)]) ;
    if options_.timing
        disp (['        Back substitution               : ' num2str(hbacsup)]) ;
    end
    if err < options_.dynatol
        h1 = etime(clock,h1) ;
        fprintf ('\n') ;
        disp (['        Total time of simulation        : ' num2str(h1)]) ;
        fprintf ('\n') ;
        disp (['        Convergence achieved.']) ;
        disp (['-----------------------------------------------------']) ;
        fprintf ('\n') ;
        return ;
    end
end
disp(['WARNING : the maximum number of iterations is reached.']) ;
fprintf ('\n') ;
disp (['-----------------------------------------------------']) ;
