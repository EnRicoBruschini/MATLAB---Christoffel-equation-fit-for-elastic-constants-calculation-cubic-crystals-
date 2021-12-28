% Christoffel cubic - Fit of velocity data to find elastic constants and
% crystal orientation (cubic crystals)


% Read data
% The input .txt file has to be prepared in the following way:
% column1 -- azimuthal angles (in degree)
% column2 -- measured velocities (km/s)
% column3 -- standard deviations of the velocities in column2
% column4 -- wave mode (0 = Vp; 1 = Vs_first transverse polarization; 2 = Vs_second transverse polarization)


close all; clear all;

fid = fopen('test_file.txt');
M  = textscan(fid,'%f %f %f %f' ,'Headerlines',1,'CommentStyle','%');
MM = [M{1} M{2} M{3} M{4}];

% Assign variables
mode = M{4};                        % Assign the mode to each velocity
vpi  = find(mode==0);               
vs1i = find(mode==1);               
vs2i = find(mode==2);               

for i = 1:length(vpi)               % Create a matrix with the Vps     
    vp_m(i,:) = MM(vpi(i),:);
end

for i = 1:length(vs1i)              % Create a matrix with the Vs1s
    vs1_m(i,:) = MM(vs1i(i),:);
end

for i = 1:length(vs2i)              % Create a matrix with the Vs2s
    vs2_m(i,:) = MM(vs2i(i),:);     
end

chivp    = vp_m(:,1);                % Rotation angle of the platelet
vp_e     = vp_m(:,2);                % Vp's vector
err_vp = vp_m(:,3);                  % Standard deviation of Vp's

chivs1    = vs1_m(:,1);              % Rotation angle of the platelet
vs1_e     = vs1_m(:,2);              % Vs1's vector
err_vs1   = vs1_m(:,3);              % Standard deviation of Vs1's

chivs2    = vs2_m(:,1);              % Rotation angle of the platelet
vs2_e     = vs2_m(:,2);              % Vs2's vector
err_vs2   = vs2_m(:,3);              % Standard deviation of Vs2's


rho = 3.6244;              % Density (g/cm^3)

vel = [vp_e; vs1_e; vs2_e];               % Creates a velocities matix to be used in the least square minimization routine
err = [err_vp; err_vs1; err_vs2];         % Creates an errors matix to be used in the least square minimization routine

c110 = 311.14;            % Initial guess for c11 (GPa)
c120 = 186;               % Initial guess for c12 (GPa)
c440 = 156;               % Initial guess for c44 (GPa)
omega0 = 49;              % Initial guess for omega (degrees)
chi0   = 42;              % Initial guess for chi (degrees)
theta0 = -7;              % Initial guess for theta (degrees)

init_val = [omega0 chi0 theta0 c110 c120 c440];
options = optimset('TolFun',1e-8,'TolX',1e-0,'MaxFunEvals',4000,'MaxIter',2000);
[coeff fval exfl] = fminsearch(@(var) csdf_christc(var,chivp,chivs1,chivs2,vp_e,vs1_e,vs2_e,err_vp,err_vs1,err_vs2,rho),init_val,options);


% Create a calculated model to plot the
% experimental velocities against the calculated ones.

orient   = coeff(1:3);          cij = coeff(4:6);

model_vp = chriscub(orient,chivp,cij,rho);
vp_c  = model_vp(:,3);

model_vs1 = chriscub(orient,chivs1,cij,rho);
vs1_c  = model_vs1(:,1);

model_vs2 = chriscub(orient,chivs2,cij,rho);
vs2_c  = model_vs2(:,2);

chideg   = 0:180;
model    = chriscub(orient,chideg,cij,rho);
vp_calc  = model(:,3);
vs1_calc = model(:,1);
vs2_calc = model(:,2);

figure(1)
err_vp_05 = err_vp./2; 
errorbar(chivp,vp_e,err_vp_05,'r*')
hold on
plot(chideg,vp_calc,'b-')
hold off
xlim([0 180])

figure(2)
err_vs1_05 = err_vs1./2;
err_vs2_05 = err_vs2./2;
plot(chideg,vs1_calc,'b-',chideg,vs2_calc,'g-')

hold on
errorbar(chivs1,vs1_e,err_vs1_05,'r*')
errorbar(chivs2,vs2_e,err_vs2_05,'b*')
hold off
xlim([0 180])





