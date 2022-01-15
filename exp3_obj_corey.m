
%% Objective function using Corey representation of rel perms and pc
function residual = exp3_obj_corey(var0, param)

tic

% % Check whether the parameters lie within their physical bounds
% if any(cell2mat(var0) <= 0)
%     var_names = {'krw', 'kro', 'nw' , 'no' , 'Swr', 'Sor' , 'pd', 'lambda'};
%     vals = cell2mat(var0);
%     ind_neg = vals <= 0;
%     disp(['Replacing negative ' strjoin(var_names(ind_neg), {' '}) ' with a small value..']);
% end

% Check whether the parameters lie within their physical bounds
if any(var0 <= 0)
    ind_neg = var0 <= 0;
    if nnz(ind_neg) > 1
        par_neg = strjoin(param.var_names(ind_neg), {' and '});
    else
        par_neg = param.var_names{ind_neg};
    end
    disp(' ')
    disp(['Replacing negative ' par_neg ' with a small value..']);
    disp(' ')
    var0(ind_neg) = 1e-8;
end

% Assign the Corey parameters
var0 = num2cell(var0);
[krw, kro, nw, no, Swr, Sor, pd, lambda] = deal(var0{:});

% Get the simulation setup
G = param.G;
rock = param.rock;
fluid = param.fluid;
state = param.state;
schedule = param.schedule;
t0 = param.t0;
pv = param.pv;
t_exp = param.t_exp;
pv_exp = param.pv_exp;
dp_exp = param.dp_exp;
q_exp = param.q_exp;
ind_exp = param.ind_exp;
ind_sim = param.ind_sim;

plot_snapshots = param.plot_snapshots;
plot_fluid_props = param.plot_fluid_props;
plot_results_t = param.plot_results_t;
plot_results_pv = param.plot_results_pv;


% % Setup the Corey relative permeabilities and the capillary pressure function
% f.krW  = coreyPhaseRelpermAD(nw, Swr, krw, Swr + Sor);
% f.krO = coreyPhaseRelpermAD(no, Sor, kro, Swr + Sor);  

f.krW  = myCoreyPhaseRelpermAD(nw, Swr, krw, Swr + Sor);
f.krO = myCoreyPhaseRelpermAD(no, Sor, kro, Swr + Sor);  
f.pcOW = coreyCapPressureAD(pd, lambda, Swr, Sor);
% eps = 1e-8;
% f.pcOW = @(sw) pd*max(min(((sw - (Swr - eps))./(1 - Swr - Sor)), 1), 0).^(-1./lambda);

% Unless explicitly provided in param, infer the line color from the 
% number of entries in the residuals file
if isfield(param, 'color')
    color = param.color;    
else
    ncol = 10;
    colors = colorcube(ncol + 1);
    % Remove yellow
    colors = colors(2:end, :);
    try
        % Get the current number of iterations (excluding the header)
        res = dlmread('iterations.csv', ',', 1, 0);
        ires = size(res, 1);
    catch
        % dlmread produces error when the file only contains the header
        ires = 0;
    end
    ic = rem(ires, ncol) + 1;
    color = colors(ic, :);
end

if plot_fluid_props
    figure(98)
    sw = [0:0.01:1];
    plot(sw, f.krW(sw), 'LineStyle', '-', 'Color', color)
    hold on
    plot(sw, f.krO(1-sw), 'LineStyle', '--', 'Color', color)
    title('Relative permeabilities')
    legend('k_{rw}', 'k_{ro}')
    xlabel('S_w')
    figure(99)
    plot(sw, f.pcOW(sw)/barsa, 'Color', color)
    hold on
    title('Capillary pressure')
    xlabel('S_w')
    ylabel('p_c (bar)')
    drawnow
end

%% Introduce a second saturation region for the boundary cells

% The 1st region corresponds to the interior cells, the second one - to the
% boundary cells with the same relative permeabilities as for the interior,
% but with pc = 0
fluid.krW = {f.krW, f.krW};
fluid.krO = {f.krO, f.krO};
fluid.pcOW = {f.pcOW, @(sw)0};

% Get the indices of the inflow and outflow boundaries
ib_inj = boundaryFaceIndices(G, 'West');
ib_prod = boundaryFaceIndices(G, 'East');

% Indices of cells adjacent to the inflow and outflow boundaries
ic_inj = max(G.faces.neighbors(ib_inj, :));  
ic_prod = max(G.faces.neighbors(ib_prod, :)); 

% Mark cells with the corresponding saturation region
reg = ones(G.cells.num, 1);
reg(ic_inj) = 2;
reg(ic_prod) = 2;

% Set up the different regions 
rock.regions = struct('saturation', reg);


%% Set up the two-phase oil-water model

model = TwoPhaseOilWaterModel(G, rock, fluid);

%% Run the simulation

% Set a max number of iterations greater than default (=6)
maxTimestepCuts = 10;
%maxTimestepCuts = 3;
solver = NonLinearSolver('maxTimestepCuts', maxTimestepCuts);

% The simulation can fail due to too coarse time steps.
% Refine the time steps until convergence.
converged = false;
while ~converged
    try
        [~, states] = simulateScheduleAD(state, model, schedule, ...
                                         'NonLinearSolver', solver);
        converged = true;
    catch
        disp(' ')
        disp('Simulation failed due to too coarse time steps...')
        disp('Refining all time steps by factor 2..')
        ndt = numel(schedule.step.control);
        cc = repmat(schedule.step.control, 1, 2);
        schedule.step.control = reshape(cc',  ndt * 2, 1);
        ccval = repmat(schedule.step.val / 2, 1, 2);
        schedule.step.val = reshape(ccval', ndt * 2, 1);
    end
end

% Get the simulation results
ns = numel(states);
t = t0 + cumsum(schedule.step.val);

% Get the time-dependent parameters
[pw_inj, pw_prod, pn_inj, pn_prod, ...
 qw_inj, qn_inj, qw_prod, qn_prod, sw_aver, pvi] = deal(zeros(ns, 1));

if plot_snapshots
    fig = figure(100);
end

t_prev = 0;
for n = 1:ns
    
    % Influxes and outfluxes
    qw_inj(n) = sum(states{n}.flux(ib_inj, 1));
    qn_inj(n) = sum(states{n}.flux(ib_inj, 2));
    qw_prod(n) = sum(states{n}.flux(ib_prod, 1));
    qn_prod(n) = sum(states{n}.flux(ib_prod, 2));
    
    % Pore volumes injected
    pvi(n) = (qw_inj(n) + qn_inj(n)) * (states{n}.time - t_prev) / pv;
    t_prev = states{n}.time;
    
    % Simulation time
    %t(n) = states{n}.time;
        
    % Pressures
    pn_inj(n) = mean(states{n}.pressure(ic_inj));
    sw = mean(states{n}.s(ic_inj, 1));
    if isfield(fluid, 'pcOW')
        pc = fluid.pcOW{2}(sw);
    else
        pc = 0;
    end
    pw_inj(n) = pn_inj(n) - pc;

    pn_prod(n) = mean(states{n}.pressure(ic_prod));
    sw = mean(states{n}.s(ic_prod, 1));
    if isfield(fluid, 'pcOW')
        pc = fluid.pcOW{2}(sw);
    else
        pc = 0;
    end
    pw_prod(n) = pn_prod(n) - pc;
    
    % Average water saturation
    sw_aver(n) = mean(states{n}.s(:, 1));
    
    if plot_snapshots
        clf
        hp = uipanel('Parent', fig, 'BorderType', 'none', ...
                    'TitlePosition', 'centertop', 'FontSize', 10); 
        set(hp, 'Title', ['At ' num2str(t(n)/hour) ' hours']);

        subplot(1, 2, 1, 'Parent', hp)
        plotCellData(G, states{n}.pressure / barsa);
        view(3), colorbar
        axis equal
        title('Pressure [bar]')

        subplot(1, 2, 2, 'Parent', hp)
        plotCellData(G, states{n}.s(:,1));
        view(3), colorbar
        axis equal
        title('Water saturation')   
        
        drawnow
    end

end

pvi = cumsum(pvi);

%% Calculate the residual between the experimental and simulation data
% The time instants at which dp is calculate are slightly different for the
% simulation and experimental data because of different time intervals;
% ignore this difference for the moment.
dp_sim = (pw_inj - pw_prod)/barsa;
residual = norm(dp_sim(ind_sim) - dp_exp(ind_exp), 2)^2;


%% Save the iterations results 
dlmwrite(param.fname, [residual, var0], 'delimiter', ',', '-append');
dlmwrite(param.dpname, dp_sim', 'delimiter', ',', '-append');

el_time = toc;
dlmwrite(param.tname, el_time, 'delimiter', ',', '-append');


%% Plotting

if plot_results_pv
    figure(2)
    plot(pv_exp, dp_exp, 'k')
    hold on
    plot(pvi, (pw_inj - pw_prod)/barsa, 'Color', color)
    xlabel('PV')
    ylabel('Differential pressure (bar)')
    legend('Experimental \DeltaP', 'Simulation')

    figure(3)
    plot(pv_exp, q_exp)
    hold on
    plot(pvi, (qw_inj + qn_inj)/(milli*litre/hour), 'Color', color)
    xlabel('PV')
    ylabel('q_{t, inj} (ml/hour)')
    legend('Experimental', 'Simulation')

    figure(4)
    plot(pvi, sw_aver, 'Color', color)
    xlabel('PV')
    ylabel('Average S_w')
end

if plot_results_t
    figure(5)
    plot(t_exp/hour, dp_exp, 'k')
    hold on
    plot(t/hour, (pw_inj - pw_prod)/barsa, 'Color', color)
    xlabel('Time from start (hours)')
    ylabel('Differential pressure (bar)')
    legend('Experimental \DeltaP', 'Simulation')

    figure(6)
    plot(t_exp/hour, q_exp, 'k')
    hold on
    plot(t/hour, (qw_inj + qn_inj)/(milli*litre/hour), 'Color', color)
    xlabel('Time from start (hours)')
    ylabel('q_{t, inj} (ml/hour)')
    legend('Experimental', 'Simulation')

    figure(7)
    plot(t/hour, sw_aver, 'Color', color)
    xlabel('Time from start (hours)')
    ylabel('Average S_w')
end   

end











