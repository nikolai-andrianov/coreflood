%% History match of the data for Greensand Phase 1 Experiment 3

function exp3_hm

mrstModule add ad-core ad-props ad-blackoil

plot_orig = false;
plot_orig = true;
plot_calib = true;

%do_hm = false;
do_hm = true;

% Select a representation for relative permeabilities and capillary pressure
krpc = 'Corey';
%krpc = 'Genuchten';
%krpc = 'Pointwise';


%% Fluid properties

% Densities and viscosities for water and CO2
rhow = 1.0630*gram/(centi*meter)^3;
rhon = 0.7237*gram/(centi*meter)^3;
muw = 0.565*centi*poise;
mun = 0.06004*centi*poise;

%% Read the reference solution
T = readtable('Exp3.csv', 'ReadVariableNames', false, 'HeaderLines', 2);

T.Properties.VariableNames = {'Date', 'Time', 'DateTime', 'TotRate', ...
    'TotCVol','density', 'InletP', 'OutletP', 'FlowStatus', 'FluidID', ...
    'BrinePerm', 'CO2Perm', 'PV', 'dP'};

%% Convert the data entries
Date = T.Date;
Time = T.Time;
q_exp = T.TotRate;
density_exp = T.density;
Kbrine_exp = str2double(T.BrinePerm);
KCO2_exp = str2double(T.CO2Perm);
pv_exp = str2double(T.PV);
dp_exp = str2double(T.dP);

%% Get the time instants

% Apparently T.Day + T.Time <> T.DateTime 
% Also, the seconds in DateTime are always zero
% Therefore creating a new DateTime = T.Day + T.Time
DateTime = strcat(T.Date, {' '}, T.Time);
dn = datenum(DateTime, 'dd-mm-yy HH:MM:SS');
ddt = diff(dn) * day;
t_exp = cumsum(ddt);


%% Cut the initial data points when flow = 0
cut_zero_flow = true;

% Cut all FPs with zero flow
if cut_zero_flow
    iflow = find(q_exp);
    q_exp = q_exp(iflow);
    density_exp = density_exp(iflow);
    Kbrine_exp = Kbrine_exp(iflow);
    KCO2_exp = KCO2_exp(iflow);
    pv_exp = pv_exp(iflow);
    dp_exp = dp_exp(iflow);
    DateTime = DateTime(iflow);
    %dn = dn(iflow);
    ddt = ddt(iflow);
    t_exp = cumsum(ddt);
    %dn = t_exp;
end

% Keep these among the parameters of the objective function
param.t_exp = t_exp;
param.pv_exp = pv_exp;
param.dp_exp = dp_exp;
param.q_exp = q_exp;

%% Detect the brine permeability higher than Klim
Klim = 1200;
ia = find(Kbrine_exp > 1200);
t1 = DateTime{ia(1)};
t2 = DateTime{ia(end)};
disp(' ')
disp(['Brine permeability > ' num2str(Klim) ' detected between' t1 ' and ' t2])
%disp([t1 ' and ' t2])
disp(['   ... setting these values to ' num2str(Klim) ' for plotting'])
Kbrine_exp(ia) = Klim;

%% Idenify different flow periods
disp(' ')
disp('Original flow periods:');
eps = 0.1;
rho_thresh = 0.5*(rhow + rhon);
n = 1;
qinj(n) = q_exp(n);
inj(n) = 1;
qtmp = q_exp;
while 1
    % Define a flow period with rates close to the first value in that period
    ind = find(abs(qtmp - qinj(n)) < eps);    
    
    % The flow period starts at the current time index
    ind = ind(ind >= inj(n));
    
    % The flow period has to be continuous
    gap = diff(ind);
    ig = find(gap > 1);
    if ~isempty(ig)
        ind = ind(1:ig(1));
    end
    
    % Idenify FPs of water or CO2 injection by comparing the density readings
    % with a density threshold
    rho = mean(density_exp(ind)) * gram/(milli*litre);
	if rho < rho_thresh
        wat_inj(n) = 0;
        phase = 'CO2';
    else
        wat_inj(n) = 1;
        phase = 'brine';
    end
        
    inj(n + 1) = ind(end) + 1;
    disp(['FP #' num2str(n) ' is between ' ...
        DateTime{inj(n)} ' and ' DateTime{ind(end)} ...
        ' with ' phase ' rate of ' num2str(qinj(n)) ])
    if ind(end) + 1 < numel(q_exp)
        qinj(n + 1) = q_exp(ind(end) + 1);
    else
        %qinj(n + 1) = 0;
        % Set the proper length of flow periods
        %inj = inj(1:end - 1);
        % Idenify the end of the last FP 
        inj(n + 1) = ind(end);
        qinj(n + 1) = 0;
        wat_inj(n + 1) = wat_inj(n);
        break;
    end
    qtmp(ind) = qinj(n);
    n = n + 1;    

end



%% Plotting

if plot_orig
    figure(1)
    % -------------------------------------------------------------------------
    subplot(2, 1, 1)
    [ax1, h1, h2] = plotyy(pv_exp, q_exp, [pv_exp, pv_exp], [KCO2_exp, Kbrine_exp]);
    set(ax1(1), 'xlim', [0 max(pv_exp)])
    set(ax1(2), 'xlim', [0 max(pv_exp)])
    %set(ax(2), 'ylim', [0 Klim])

    set(get(ax1(1), 'ylabel'), 'string', 'Rate (ml/h)')
    set(get(ax1(2), 'ylabel'), 'string', 'Permeability (mD)')

    % set(h1, 'Color', [1 0.5 0])
    % set(h2(1), 'Color', 'r')
    % set(h2(2), 'Color', 'b')

    xlabel('PV')
    legend('Rate', 'K_{CO2}', 'K_{brine}')
    
    % Mark the flow periods
    hold(ax1(1))
    plot(ax1(1), pv_exp(inj(1:end-1)),  qinj(1:end-1), 'o')
    ifp = [1:numel(qinj(1:end-1))];
    sv = strcat(' FP#', num2str(ifp(:)));
    text(pv_exp(inj(1:end-1)),  qinj(1:end-1), sv)
    hold(ax1(1))

    % -------------------------------------------------------------------------
    subplot(2, 1, 2)
    [ax2, h1, h2] = plotyy(pv_exp, dp_exp, pv_exp, density_exp);
    set(ax2(1), 'xlim', [0 max(pv_exp)])
    set(ax2(2), 'xlim', [0 max(pv_exp)])

    set(get(ax2(1), 'ylabel'), 'string', '\Deltap (bar)')
    set(get(ax2(2), 'ylabel'), 'string', 'Density (g/ml)')
    % set(h2, 'Marker', 'o')
    % set(h1, 'Color', 'g')
    % set(h2, 'Color', 'm')

    xlabel('PV')
    legend('\Deltap', '\rho')

end

%% Sample data
Lx = 15.053*centi*meter; 
D = 3.794*centi*meter;
A = pi*D^2/4;

phi = mean([35.4, 34.5, 34.7])/100;             % Average of porosities for samples 100, 101, 102
Kg = harmmean([1140, 1230, 1240])*milli*darcy;  % Average gas permeability
pv = phi * A * Lx;

% Keep the PV among the parameters of the objective function
param.pv = pv;

%% Define the simulation geometry

N = 500;    % Number of gridblocks in the flow direction (x)
M = 1;      % Number of gridblocks in the transverse direction (y & z)

% Represent the core as a rectangular parallelepiped 
Ly = sqrt(A);
Lz = Ly;

%G = cartGrid([N M M],[Lx Ly Lz]);

% Refine the cells towards the edges of the core
dxi = Lx / (N - 5);
dx = [0.1, 0.2, 0.4, 0.8, ones(1, N - 8), 0.8, 0.4, 0.2, 0.1]' * dxi;
x = [0; cumsum(dx)];
y = [0: Ly/M: Ly]';
z = [0: Lz/M: Lz]';

% Add two ghost cells at the inflow and outlfow
dx = [0.1*dxi; dx; 0.1*dxi];
x = [0; cumsum(dx)] - 0.1*dxi;

G = computeGeometry( tensorGrid(x, y, z) );

% Keep the geometry among the parameters of the objective function
param.G = G;

%% Rock properties

% Create the rock structure
rock = makeRock(G, Kg, phi);

% Set rock compressibility
pref = 1 * barsa;
cr = 0;
rock.cr = cr;
rock.pref = pref;

% plotCellData(G, rock.perm/milli/darcy), view(3), colorbar

% Keep rock among the parameters of the objective function
param.rock = rock;


%% Fluid model

fluid = initSimpleADIFluid('phases','WO',    ... % Fluid phases: water and oil
                           'mu',  [muw, mun],      ... % Viscosities
                           'rho', [rhow,rhon],     ... % Surface densities [kg/m^3]
                           'c',   [0, 0],          ... % Fluid compressibility[Cm, Cm],              
                           'pRef', 0,              ... % Reference pressure 
                           'cR',  0                ... % Rock compressibility
                           );

% Keep fluid among the parameters of the objective function
param.fluid = fluid;


%% Initial guess for the optimization variables

if strcmp(krpc, 'Corey')
    % Assign the initial guess for the Corey parameters
    % The critical saturations and the end-point permeabilities are taken from
    % "AF790 Kw Ko Rt Ro".xls
    % The initial parameters for the Corey capillary pressure are from
    % MICP_100M.xlsx
    % termed ini0
    var0 = [292/1140, 0.75, 3, 3, 0.376, 0.16, 0.04*barsa, 2]; % [krw, kro, nw, no, Swr, Sor, pd, lambda]

    % Use initial guesses such that krw > kro as supported from the
    % comparison with the imbibition data, Bennion & Bachu SPE99326, and
    % Ott et al SCA2014-019.
    % termed ini1
    var0 = [0.7, 0.3, 3, 3, 0.376, 0.16, 0.04*barsa, 2]; % [krw, kro, nw, no, Swr, Sor, pd, lambda]
    
    % Test various initial guesses.
    % termed ini2
    var0 = [0.5, 0.5, 2, 2, 0.376, 0.16, 0.04*barsa, 2]; % [krw, kro, nw, no, Swr, Sor, pd, lambda]
    
    % Test various initial guesses.
    % termed ini3
    var0 = [0.5, 0.5, 3, 3, 0.376, 0.16, 0.04*barsa, 2]; % [krw, kro, nw, no, Swr, Sor, pd, lambda]        
    
%     disp('test pc = 0');
%     var0 = [292/1140, 0.75, 3, 3, 0.376, 0.16, 0, 2]; % [krw, kro, n1w, no, Swr, Sor, pd, lambda]    
    
    % Assign the Corey parameters
    var0c = num2cell(var0);
    [krw, kro, nw, no, Swr, Sor, pd, lambda] = deal(var0c{:});
    
elseif strcmp(krpc, 'Genuchten') 
    
    % Assign the van Genuchten parameters  46
    % The critical saturations and the end-point permeabilities are taken from
    % "AF790 Kw Ko Rt Ro".xls
    % eps, gamma from Helmig p. 46
    % aplha, m, n from Helmig p. 37
    var0 = [292/1140, 0.75, 0.5, 1/3, 0.376, 0.16, 0.37, 0.77, 4.37]; % [krw, kro, eps, gamma, Swr, Sor, alpha, m, n]

    % Testing
    %var0 = [1, 1, 0.5, 1/3, 0.1, 0.1, 0.37, 0.77, 4.37]; % [krw, kro, eps, gamma, Swr, Sor, alpha, m, n]    
    %var0 = [292/1140, 0.75, 3, 3, 0.376, 0.16, 1/(0.04*barsa), 2, 1]; % [krw, kro, n1w, no, Swr, Sor, pd, lambda]
    %var0 = [292/1140, 0.75, 0.5, 1/3, 0.376, 0.16, 0, 0.77, 4.37]; % [krw, kro, eps, gamma, Swr, Sor, alpha, m, n]
    
    var0c = num2cell(var0);
    [krw, kro, eps, gamma, Swr, Sor, alpha, m, n] = deal(var0c{:});    

elseif strcmp(krpc, 'Pointwise')
    
    % Pointwise representation of kr and pc leads to convergence issues.
    % The number of points in the tabular representation has to be limited
    % in order to reduce the degrees of freedom during the optimization.
    % If linear interpolation is used to create kr and pc, the curves are
    % non-smooth at the interpolation points, leading to poor convergence
    % and at some point the forward solution breaks, even if started from
    % the optimal Corey kr and pc shapes.
    % Using spline or cubic interpolation in kr and pc curves creates weird 
    % curve shapes (e.g. with under- or overshoots), which again leads to 
    % non-convergence.     
    
    % Initial guess for the tabular representation of kr - linear functions
    param.nsw = 5;
    Swr = 0.1; 
    Sor = 0.1;
    sw = linspace(Swr, 1 - Sor, param.nsw);
    krw = linspace(0, 1, parim.nsw);
    kro = flip(krw);
    % Use Corey capillary pressure as an initial guess for the tabular representation of pc  
    % Adjust pc(Swr) by eps to avoid Inf 
    [pd, lambda] = deal(0.04*barsa, 2);
    eps = 1e-8;
    pc_corey = @(sw) pd*max(min(((sw - (Swr - eps))./(1 - Swr - Sor)), 1), 0).^(-1./lambda);
    pc = pc_corey(sw);
    
%     % Initial guess from the Corey optimal solution
%     results = dlmread('exp3_iterations.csv', ',', 1, 0);
%     err = results(:, 1);
%     par = results(:, 2:end);      
%     [err_min, ie] = min(err);
%     par_opt = par(ie, :);
%      
%     var0 = num2cell(par_opt);
%     [krw, kro, nw, no, Swr, Sor, pd, lambda] = deal(var0{:});  
%     sw = linspace(Swr, 1 - Sor, param.nsw);
%     
%     krW  = coreyPhaseRelpermAD(nw, Swr, krw, Swr + Sor);
%     krO = coreyPhaseRelpermAD(no, Sor, kro, Swr + Sor);  
%     pc_corey = @(sw) pd*max(min(((sw - (Swr - eps))./(1 - Swr - Sor)), 1), 0).^(-1./lambda);
%     krw = krW(sw);
%     kro = krO(1 - sw);
%     pc = pc_corey(sw);
    
    % Exclude points at Swr and 1-Sor because krw(Swr) = 0 and kro(1-Sor) = 0 
    var0 = [Swr, Sor, krw(2:end), kro(1:end-1), pc];    
else
    disp('The representation for relative permeabilities and capillary pressure is not defined!'); 
    return;
end 


%% Initial conditions

% Set initial pressure to 200 bar
p0 = 200*barsa;
state.pressure = ones(G.cells.num,1) * p0;  

% Initially saturated with brine
state.s = repmat([1 0], G.cells.num,1);


%% Boundary conditions

% Keep only the FPs longer than 1 minute 
dt_fp = diff(t_exp(inj));
idx = find(dt_fp > 60);
dt_fp = dt_fp(idx);
qinj = qinj(idx);
wat_inj = wat_inj(idx);

disp(' ')
disp('Keeping only the FPs longer than 1 minute:');
phase = {'CO2', 'brine'};
for n = 1:numel(idx)
    disp(['FP #' num2str(idx(n)) ' with ' phase{wat_inj(n) + 1} ' rate of ' num2str(qinj(n)) ])
end

% Simulate only the selected FPs
 sel_fp = [1, 2];
%  sel_fp = [3];
%  sel_fp = [1: numel(idx)];
str = 'Simulating only the FPs: ';
phase = {'CO2', 'brine'};
for n = 1:numel(sel_fp)
    str = [str num2str(idx(sel_fp(n))) ' ']; 
end

% The initial time instant
if sel_fp(1) > 1
    t0 = sum(dt_fp(1:sel_fp(1) - 1));
else
    t0 = 0;
end

disp(' ')
disp(str)
idx = idx(sel_fp);
dt_fp = dt_fp(sel_fp);
qinj = qinj(sel_fp);
wat_inj = wat_inj(sel_fp);

% Overwrite the initial conditions 
if sel_fp(1) > 1
   fac = 1.1; % 1.2221
   Sw_ini = 0.4595; % The value of average Sw after CO2 flooding
   disp(' ')
   disp(['Setting the initial conditions to ' num2str(fac) ' * Swr..'])
   state.s = repmat([fac*Swr 1 - fac*Swr], G.cells.num,1);
%    disp(['Setting the initial conditions to ' num2str(Sw_ini)])
%    state.s = repmat([Sw_ini 1 - Sw_ini], G.cells.num,1);
end

% Keep the initial state among the parameters of the objective function
param.state = state;


% The volumetric rates change at specified time instants
ti = t0 + [0; cumsum(dt_fp)];
qw = qinj .* wat_inj       * (milli*litre)/hour;
qn = qinj .* (1 - wat_inj) * (milli*litre)/hour;


qt = qw + qn;

% Volumetric composition of the injected fluid  
sat = [qw; qn] ./ [qt; qt];

% Time steps, corresponding to different flow periods
dt = diff(ti);

ndt = numel(dt);

% Dirichlet BC at the right, no-flow conditions otherwise
% The initial pressure is interpreted as the right BC pressure 
bc = pside([], G, 'East', p0, 'sat', [1 0]);  
%bc = pside([], G, 'East', p0, 'sat', []);  

% A template for the Neumann BC at the left
bc = fluxside(bc, G, 'West', 1* (centi*meter)^3/hour, 'sat', [1 0]);

% Create the schedule with the template BC
schedule = simpleSchedule(dt, 'bc', bc);

% Replicate the control structure to accomodate for changing BC
schedule.control = repmat(schedule.control, ndt, 1);

% Clear the steps in the schedule as they will be overwritten with the
% small time steps within each flow period 
schedule = rmfield(schedule, 'step');

% Max time step within a FP
maxdts = 10*minute;
%maxdts = 60*minute;

% Number and length of min time steps at a start of a FP
nmin = 30;
dtmin = 1; 

% % Too coarse time steps to provoke non-convergence
% nmin = 5;
% dtmin = 10; 

% Number of linearly increasing time steps in a ramp-up
Nr = 10;

for n = 1:ndt
    
    % Adjust the injected fluid rate and composition to the current flow period
    schedule.control(n).bc.value = [p0, qt(n)]';
    schedule.control(n).bc.sat = [[1 0]; sat(:, n)'];
    
    % Introduce smaller timesteps for each flow period
    %dts = rampupTimesteps(dt(n), maxdts);
    
    % Partitition the first time interval within a FP into nmin timesteps 
    % of length dtmin, followed by Nr linearly increasing time steps 
    dtmax = min(dt(n), maxdts);
    dtn = 2*(dtmax - dtmin*nmin) / Nr - dtmin;
    dd = (dtn - dtmin) / (Nr - 1);
    dtinc = ones(Nr, 1) * dd;
    dtt = (cumsum(dtinc) - dd) + dtmin;
    dts1 = [ones(nmin, 1)*dtmin; dtt];
    dts1(end) = dtmax - sum(dts1(1:end-1));
    
    % Even steps
    dt_rem = repmat(maxdts, floor((dt(n) - maxdts)/maxdts), 1);
    % Final ministep if present
    dt_final = dt(n) - dtmax - sum(dt_rem);
    if dt_final <= 0
        dt_final = [];
    end
    % Combined timesteps
    dts = [dts1; dt_rem; dt_final];    

    % Save the small time steps in the schedule for each flow period and
    % provide an index to the corresponding entry in schedule.control
    if ~isfield(schedule, 'step')
        schedule.step.val = dts;
        schedule.step.control = ones(numel(dts), 1) * n;
    else
        schedule.step.val = [schedule.step.val; dts];
        schedule.step.control = [schedule.step.control; ones(numel(dts), 1) * n];
    end
        
end

% Keep schedule among the parameters of the objective function
param.schedule = schedule;
param.t0 = t0;

%% Select the calibration points

% Simulation tim instants
t_sim = t0 + cumsum(schedule.step.val);

% Choose nfr points within each flow period for evaluating the residual
nfr = 30;
tfp0 = t0;
ind_exp = zeros(nfr * ndt, 1);
ind_sim = zeros(nfr * ndt, 1);
for n = 1:ndt
    
    % Split each FP in a sequence of nr intervals with linearly increasing
    % lengths
    dtmin = dt(n) / 100;    
    dtn = 2*dt(n) / nfr - dtmin;
    dd = (dtn - dtmin) / (nfr - 1);
    dtinc = ones(nfr, 1) * dd;
    dtt = (cumsum(dtinc) - dd) + dtmin;
    dts1 = [dtmin; dtt];
    dts1(end) = dt(n) - sum(dts1(1:end-1));  
    
    % Get the indices of the corresponding time instants in the
    % experimental data and in simulation results
    tt = tfp0 + cumsum(dts1(1:end-1));
    for i = 1:nfr
        [~, ind_exp((n - 1)*nfr + i)] = min(abs(t_exp - tt(i)));
        [~, ind_sim((n - 1)*nfr + i)] = min(abs(t_sim - tt(i)));
    end
    
    % Start of the next FP
    tfp0 = tfp0 + sum(dt(n));
    
end

% Visualize the calibration points
if plot_calib
    figure(100)
    plot(t_exp/hour, dp_exp, 'k')
    hold on
    plot(t_exp(ind_exp)/hour, dp_exp(ind_exp), 'o')
end

% Pass the calibration points to the objective function
param.ind_exp = ind_exp;
param.dp_exp = dp_exp;
param.ind_sim = ind_sim;

% Save the experimental data vs time
expname = 'exp3_tdp.csv';
exp_tdp = [t_exp'/hour; dp_exp'];
dlmwrite(expname, exp_tdp);
disp(' ')
disp(['Experimental data vs time is saved in ' expname])
disp(' ')


%% Fitting the Corey parameters to match the experimental dp curve

if do_hm
    param.plot_snapshots = false;
    param.plot_fluid_props = true;
    param.plot_results_t = true;
    param.plot_results_pv = false;
    
    param.fname = 'iterations.csv';  
    param.dpname = 'iterations_dp.csv';    
    param.tname = 'iterations_time.csv';    
    
    disp(' ')
    if strcmp(krpc, 'Corey')
        
        disp('Using Corey representation for relative permeabilities and capillary pressure..');        
        % The names of optimization parameters
        param.var_names = {'krw', 'kro', 'nw' , 'no' , 'Swr', 'Sor' , 'pd', 'lambda'};        
        fun_exp3_obj = @(var0) exp3_obj_corey(var0, param); 
        
    elseif strcmp(krpc, 'Genuchten') 
        
        disp('Using van Genuchten representation for relative permeabilities and capillary pressure..');        
        % The names of optimization parameters
        param.var_names = {'krw', 'kro', 'eps', 'gamma' , 'Swr', 'Sor' , 'aplha', 'm', 'n'}; 
        fun_exp3_obj = @(var0) exp3_obj_genuchten(var0, param);         
            
    elseif strcmp(krpc, 'Pointwise')
        disp('Using pointwise representation for relative permeabilities and capillary pressure..');
        
        % The names of optimization parameters
        ind = [1:param.nsw];    
        cind = cellstr(num2str(ind'));
        krw_ind = strcat('krw', cind(1:end-1));
        krn_ind = strcat('krn', cind(1:end-1));
        pc_ind = strcat('pc', cind);
        param.var_names = [{'Swr', 'Sor'}, krw_ind', krn_ind', pc_ind']; 
    
        fun_exp3_obj = @(var0) exp3_obj_pointwise(var0, param);    
    else
        disp('The representation for relative permeabilities and capillary pressure is not defined!'); 
        return;
    end    
    disp(' ')      

    % Save the residuals and the optimization parameters in fname 
    header = [{'residual'}, param.var_names];
    header = strjoin(header, repmat({','}, 1, numel(param.var_names)));
    fileID = fopen(param.fname, 'w');
    fprintf(fileID,'%s\n', header);
    fclose(fileID);

    % Save the calculated time instants (1st line) and dp's (other lines) in dpname 
    dlmwrite(param.dpname, t_sim'/hour, 'delimiter', ',');

    % Save the elapsed time per iteration in tname 
    tit = 'Elapsed time (sec)';
    fileID = fopen(param.tname, 'w');
    fprintf(fileID,'%s\n', tit);
    fclose(fileID);  
    
    % Unconstrained minimization using Nelder-Mead
    [var_opt, fval] = fminsearch(fun_exp3_obj, var0);

else

    disp(' ')
    
    % Residuals and Corey parameters per iteration
    results = dlmread('exp3_corey_drain_iterations.csv', ',', 1, 0);
    err = results(:, 1);
    par = results(:, 2:end);      
    [err_min, ie] = min(err);
    par_opt = par(ie, :);
    par_optc = num2cell(par_opt);
    
    if strcmp(krpc, 'Corey')
        % Setup the Corey relative permeabilities and the capillary pressure function
        [krw, kro, nw, no, Swr, Sor, pd, lambda] = deal(par_optc{:});
        f.krW  = coreyPhaseRelpermAD(nw, Swr, krw, Swr + Sor);
        f.krO = coreyPhaseRelpermAD(no, Sor, kro, Swr + Sor);  
        f.pcOW = @(sw) pd*max(min(((sw - Swr)./(1 - Swr - Sor)), 1), 0).^(-1./lambda); 
    elseif strcmp(krpc, 'Genuchten')
        disp('To be done..');  
        return;        
    elseif strcmp(krpc, 'Pointwise')
        disp('To be done..');  
        return;
    else
        disp('The representation for relative permeabilities and capillary pressure is not defined!'); 
        return;
    end              
    
    % Plot the simulation results for the optimal set of Corey parameters
    plot_saved = false;
    plot_saved = true;
    if plot_saved
        
        figure(98)
        sw = [0:0.01:1];
        plot(sw, f.krW(sw), 'LineStyle', '-', 'Color', 'b')
        hold on
        plot(sw, f.krO(1-sw), 'LineStyle', '-', 'Color', 'r')
        title('Relative permeabilities')
        legend('k_{rw}', 'k_{ro}')
        xlabel('S_w')
        figure(99)
        plot(sw, f.pcOW(sw)/barsa, 'b')
        hold on
        title('Capillary pressure')
        xlabel('S_w')
        ylabel('p_c (bar)')  
    
        disp('Plotting the saved best history matched parameters..')        
        % Simulation time instants and dp per iteration
        tdp = dlmread('exp3_iterations_dp.csv', ',');
        t_sim = tdp(1, :);
        % Apparently dlmread reads extra empty entries; cut them  
        [~, nsim] = max(t_sim);
        t_sim = t_sim(1:nsim);    
        dp_sim = tdp(2:end, 1:nsim);

        dp_opt = dp_sim(ie, :);

        figure
        plot(t_exp/hour, dp_exp, 'k')
        hold on
        plot(t_sim, dp_opt, 'Color', 'b')
        xlabel('Time from start (hours)')
        ylabel('Differential pressure (bar)')
        legend('Experimental \DeltaP', 'Simulation')  
    else
        disp('Run the simulation with the optimal set of Corey parameters..')
        disp(' ')
        % Run the simulation with the optimal set of Corey parameters
        param.plot_snapshots = false;
        param.plot_fluid_props = true;
        param.plot_results = true;   
        
        % Save the residuals and the Corey parameters in fname 
        param.fname = 'iterations_opt.csv';
        var_names = 'residual,krw,kro,nw,no,Swr,Sor,pd,lambda';
        %dlmwrite(param.fname, var_names, 'delimiter', ',');
        fileID = fopen(param.fname, 'w');
        fprintf(fileID,'%s\n', var_names);
        fclose(fileID);

        % Save the calculated time instants (1st line) and dp's (other lines) in dpname 
        param.dpname = 'iterations_dp_opt.csv';
        dlmwrite(param.dpname, t_sim'/hour, 'delimiter', ',');

        % Save the elapsed time per iteration in tname 
        param.tname = 'iterations_time_opt.csv';
        tit = 'Elapsed time (sec)';
        fileID = fopen(param.tname, 'w');
        fprintf(fileID,'%s\n', tit);
        fclose(fileID);    
    
        param.color = 'b';
        
        if strcmp(krpc, 'Corey')
            exp3_obj_corey(par_opt, param);    
        elseif strcmp(krpc, 'Pointwise')
            disp('To be done..');  
            return;
        else
            disp('The representation for relative permeabilities and capillary pressure is not defined!'); 
            return;
        end          
    end
    
    % Elapsed simulation time per iteration
    t_elapsed = dlmread('exp3_iterations_time.csv', ',', 1, 0);
    
    
    
    


end
