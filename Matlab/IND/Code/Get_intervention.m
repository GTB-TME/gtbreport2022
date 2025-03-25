clear all; load('estims.mat');

% Initial condition as of 2016
[out, aux] = get_objective_wRNTCP(x1, prm, ref, sel, agg, gps, data);

init = aux.soln(end,1:end-1);


x = x1;
y = [0.85 0.25];     % Proportion of private sector engaged, and proportion reduction in patient delay


r.beta = x(1);
r.careseeking = x(2);


opt = 2;

% if opt == 1
% 
% % --- Current state
% M0 = make_model(p, r, i, s, gps);
% 
% % -- Plugging the private sector gaps



if opt == 2

% --- Current state
M0 = make_model(p, r, i, s, gps);
% --- PSE only
p1 = p; r1 = r; 
p1.pu = p.pu + y(1)*(1-p.pu);              % PSE
M1 = make_model(p1, r1, i, s, gps);
% --- Plugging public sector gaps
p2 = p1; r2 = r1;
p2.Dx(1) = 0.95; p2.Tx_init(1) = 0.95;
r2.default(1) = r2.Tx*0.05./(1-0.05);
M2 = make_model(p2, r2, i, s, gps);        % Plugging public sector gaps
% --- Increase public sector access
p3 = p2; r3 = r2;
r3.access = 0.5;
M3 = make_model(p3, r3, i, s, gps);        % Increasing public sector access

% Baseline
[t0, soln0] = ode15s(@(t,in) goveqs_basis2(t, in, M0, i, s, p, sel, agg), [2016:2036], init(end,:), odeset('NonNegative',[1:i.nstates]));

% PSE only
[t1, soln1] = ode15s(@(t,in) goveqs_scaleup(t, in, M0, M1, [2016 2024], i, s, p, sel, agg), [2016:2036], init(end,:), odeset('NonNegative',[1:i.nstates]));

% PSE + public sector cascades
[t2, soln2] = ode15s(@(t,in) goveqs_scaleup(t, in, M0, M2, [2016 2024], i, s, p, sel, agg), [2016:2036], init(end,:), odeset('NonNegative',[1:i.nstates]));

% PSE + public sector cascades + public sector access
[t3, soln3] = ode15s(@(t,in) goveqs_scaleup(t, in, M0, M3, [2016 2024], i, s, p, sel, agg), [2016:2036], init(end,:), odeset('NonNegative',[1:i.nstates]));

getinc = @(t,sol)   diff(interp1(t,sol(:,i.aux.inc),t(1):t(end)),1);
getpop = @(col,sol) (sum(sol(1:end-1,col),2) + sum(sol(2:end,col),2))/2;

inc0 = getinc(t0, soln0)'./getpop(1:i.nstates, soln0);
inc1 = getinc(t1, soln1)'./getpop(1:i.nstates, soln1);
inc2 = getinc(t2, soln2)'./getpop(1:i.nstates, soln2);
inc3 = getinc(t3, soln3)'./getpop(1:i.nstates, soln3);
inc  = [inc0 inc1 inc2 inc3]*1e5;

figure; hold on; years = [2016:2035]; 

lw = 1.5; fs = 14;
p2 = plot(years, inc, 'linewidth', lw);

yl = ylim; yl(1) = 0; ylim(yl); xlim([years(1),2030]);
legend(p2,'Baseline','PPIA','PPIA + public sector cascades (TB care)','PPIA + public sector cascades + public sector access','Location','SouthEast');
ylabel('Incidence per lakh population','fontsize',fs);
set(gca,'fontsize',fs);

cinc0 = sum(inc0); cinc1 = sum(inc1);
pca = 1-cinc1/cinc0

end