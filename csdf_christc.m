function chidegf = csdf_christc(var,chivp,chivs1,chivs2,vp_e,vs1_e,vs2_e,err_vp,err_vs1,err_vs2,rho)

% Chi square for degree of freedom for the experimental data against the
% Christoffel equation solutions. 'Var' are the variables to be minimized
% stored in a vector [om0 chi0 th0 c11 c12 c44]. 'chivp': vector containing
% the degrees of rotation at which were detected the longitudinal
% velocities. 'chivs1': degrees of rotation at which were detected the
% shear waves (first polarization). 'chivs2': degrees of rotation at which
% were detected the shear waves (second polarization). 'vp_e': vector with
% the experimental longitudinal velocities. 'vs1_e' and 'vs2_e':
% experimental shear waves (first and second polarization). 'err_vp',
% 'err_vs1' and 'err_vs2' standard deviations associated with the
% experimental velocities. 'rho': density (g/cm3) of the crystal.

%
% ========________Enrico Bruschini - Roma 20/04/2012________========
% ========________    Last modified on 09/09/2012   ________========



orient = var(1:3);          cij = var(4:6);

model_vp = chriscub(orient,chivp,cij,rho);
vp_c  = model_vp(:,3);

model_vs1 = chriscub(orient,chivs1,cij,rho);
vs1_c  = model_vs1(:,1);

model_vs2 = chriscub(orient,chivs2,cij,rho);
vs2_c  = model_vs2(:,2);



% Calculate the residuals between experimental and calculated data

r1 = vp_e - vp_c;       % Residuals of experimental and calculated Vp
r2 = vs1_e - vs1_c;     % Residuals of experimental Vs1 and calculated Vs1
r3 = vs2_e - vs2_c;     % Residuals of experimental Vs2 and calculated Vs2
%r4 = vs1_e - vs2_c;     % Residuals of experimental Vs1 and calculated Vs2
%r5 = vs2_e - vs1_c;     % Residuals of experimental Vs2 and calculated Vs1

% Identifies the two shear polarization modes

%rr2 = sum(r2.^2);
%rr3 = sum(r3.^2);
%rr4 = sum(r4.^2);
%rr5 = sum(r5.^2);

%if rr2 == rr4
%    r2 = r2;
%elseif rr2 < rr4
%    r2 = r2;
%elseif rr4 < rr2
%    r2 = r4;
%end

%if rr3 == rr5
%    r3 = r3;
%elseif rr3 < rr5
%    r3 = r3;
%elseif rr5 < rr3
%    r3 = r5;
%end

Nvp = length(chivp);        Nvs1 = length(chivs1);       Nvs2 = length(chivs2);
N = Nvp + Nvs1 + Nvs2;      p = length(var);    % parameters of the degree of freedom

sq_vp  = sum((r1./err_vp).^2);
sq_vs1 = sum((r2./err_vs1).^2);
sq_vs2 = sum((r3./err_vs2).^2);

chi_vp  = sq_vp ./(N-p);          % Chi square for degree of freedom of Vp
chi_vs1 = sq_vs1./(N-p);         % Chi square for degree of freedom of Vs1
chi_vs2 = sq_vs2./(N-p);         % Chi square for degree of freedom of Vs2

chidegf = chi_vp + chi_vs1 + chi_vs2;





