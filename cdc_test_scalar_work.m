% testing cdc code - scalar equation
% error versus work (# function evaluations)

clear;

plot_str{1} = 'ko-'; % 1st order
plot_str{2} = 'b^-'; % 2nd
plot_str{3} = 'r*-'; % 3rd
plot_str{4} = 'cx-'; % 4th
plot_str{5} = 'gd-'; % 5th
plot_str{6} = 'ys-'; % 6th
plot_str{7} = 'k+--'; % 7th
plot_str{8} = 'b<--'; %8th


%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('CDC-FE, uniform grid');
opts.dc = 1; % CDC
opts.grid = 1; % uniform grid
opts.pred = 1; % FE
opts.corr = 1; % FE

tspan = [0 1]; % initial and final time
y0 = 1; % initial condition
opts.exact = 0.5; %exact solution if it exists

ti = 0;
tf = 1;

figure(1), clf

for k = 1:5
    opts.ncorr = k-1; % number of correctors
    opts.levels = k;  % total number of levels (predictor + correctors)
    opts.nquad = k+1;

    % convergence study
    approx_func_evals = [20:20:320]';
    Nk = length(approx_func_evals);

    % size of system
    m = length(y0);
    y_store = zeros(m,Nk);
    N_store = zeros(1,Nk);
    work_store = zeros(1,Nk);
    
    feval_per_interval = k*(opts.nquad - 1);
             
    for kk = 1:Nk
        %interavls 
        N = floor(approx_func_evals(kk)/feval_per_interval); % number of large time intervals
        dt = (tf - ti)/N; % grid spacing, large intervals
        tspan = [0:N]*dt;
    
        sol = deferred_correction(@scalar_ode,tspan,y0,opts);
        y_store(:,kk) = sol(:,end);
        N_store(1,kk) = N;
        work_store(1,kk) = N*feval_per_interval;
    end
    
    
    if isfield(opts,'exact')
        % exact solution exists    
        err_store = zeros(1,Nk);
        for kk = 1:Nk 
            err_store(kk) = max(abs(y_store(:,kk)-opts.exact));
        end
        loglog(work_store,err_store,plot_str{k});
        ind = (work_store>0);
        p = polyfit(log(work_store(ind)),log(err_store(ind)),1);
        rate = p(1);
        legend_str{k} = sprintf('slope=%4.1f',abs(rate));                
    else
        err_store = zeros(1,Nk-1);
        for kk = 1:(Nk - 1)
            err_store(kk) = max(abs(y_store(kk+1)-y_store(kk)));
        end
        loglog(work_store(1:end-1),err_store,plot_str{k});
        p = polyfit(log(work_store(1:end-1)),log(err_store),1);
        rate = p(1);
        legend_str{k} = sprintf('slope=%4.1f',abs(rate));
    end
    xlabel('number of function evaluations');
    ylabel('absolute error');
    title('CDC, FE, uniform grid')
    set(gca,'FontSize',24)
    hold on
    
    
        
    % print data for tikz
    % for kk = 1:Nk
    %     fprintf('(%d,%g)\n',work_store(kk),err_store(kk));
    % end
    % fprintf('\n')
        
    
end
legend(legend_str,'Location','NorthEastOutside');

%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('CDC-FE, rk3 predictor, rk2 corrector');
opts.dc = 1; % CDC
opts.grid = 1; % uniform grid
opts.pred = 2; % RK2
opts.corr = 2; % RK2

tspan = [0 1]; % initial and final time
y0 = 1; % initial condition
opts.exact = 0.5; %exact solution if it exists

ti = 0;
tf = 1;

legend_str = {};

figure(2), clf

for k = 1:3
    opts.ncorr = k-1; % number of correctors
    opts.levels = k;  % total number of levels (predictor + correctors)
                      
    % TODO FIX
    if k == 1
        opts.nquad = 2;
    else
        opts.nquad = 2 + 2*opts.ncorr + 1;
    end

    % convergence study
    approx_func_evals = [20:20:640]';
    Nk = length(approx_func_evals);

    % size of system
    m = length(y0);
    y_store = zeros(m,Nk);
    N_store = zeros(1,Nk);
    work_store = zeros(1,Nk);
    
    % predictor is rk2, 2 func evals
    % corrector is rk2, 2 func evals per corrector
    feval_per_interval = (2 + (k-1)*2) *(opts.nquad - 1);
    

             
    for kk = 1:Nk
        %interavls 
        N = floor(approx_func_evals(kk)/feval_per_interval); % number of large time intervals
        dt = (tf - ti)/N; % grid spacing, large intervals
        tspan = [0:N]*dt;
    
        sol = deferred_correction(@scalar_ode,tspan,y0,opts);
        y_store(:,kk) = sol(:,end);
        N_store(1,kk) = N;
        work_store(1,kk) = N*feval_per_interval;
    end
    
    
    if isfield(opts,'exact')
        % exact solution exists    
        err_store = zeros(1,Nk);
        for kk = 1:Nk 
            err_store(kk) = max(abs(y_store(:,kk)-opts.exact));
        end
        loglog(work_store,err_store,plot_str{k});
        ind = (work_store>0);
        p = polyfit(log(work_store(ind)),log(err_store(ind)),1);
        %p = polyfit(log(work_store),log(err_store),1);
        rate = p(1);
        legend_str{k} = sprintf('slope=%4.1f',abs(rate));                
    else
        err_store = zeros(1,Nk-1);
        for kk = 1:(Nk - 1)
            err_store(kk) = max(abs(y_store(kk+1)-y_store(kk)));
        end
        loglog(work_store(1:end-1),err_store,plot_str{k});
        p = polyfit(log(work_store(1:end-1)),log(err_store),1);
        rate = p(1);
        legend_str{k} = sprintf('slope=%4.1f',abs(rate));
    end
    xlabel('number of function evaluations');
    ylabel('absolute error');
    title('CDC, RK3 predictor, RK2 corrector, uniform grid')
    set(gca,'FontSize',24)
    hold on
    
    % print data for tikz
    for kk = 1:Nk
        fprintf('(%d,%g)\n',work_store(kk),err_store(kk));
    end
    fprintf('\n')
        
end
legend(legend_str,'Location','NorthEastOutside');
error('break')



%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('CDC, FE, Gauss-Legendre');

opts.dc = 1; % CDC
opts.grid = 2; % gauss legendre nodes
opts.pred = 1; % FE
opts.corr = 1; % FE

tspan = [0 1]; % initial and final time
y0 = 1; % initial condition
opts.exact = 0.5; %exact solution if it exists

ti = 0;
tf = 1;

figure(2), clf
opts.nquad = 4; % number of quadrature points
for k = 1:6
    opts.levels = k; % number of levels (predictor + (levels-1)
                     % correctors

    % convergence study
    Nk = 10;
    m = length(y0);
    y_store = zeros(m,Nk);
    N_store = zeros(1,Nk);
    
    for kk = 1:Nk
        N = 2*kk+10;  % number of large (equi-spaced) time interavls 
        dt = (tf - ti)/N; % grid spacing, large intervals
        tspan = [0:N]*dt;
    
        sol = deferred_correction(@scalar_ode,tspan,y0,opts);
        y_store(:,kk) = sol(:,end);
        N_store(1,kk) = N;
    end
    
    
    if isfield(opts,'exact')
        % exact solution exists    
        err_store = zeros(1,Nk);
        for kk = 1:Nk 
            err_store(kk) = max(abs(y_store(:,kk)-opts.exact));
        end
        loglog(N_store,err_store,plot_str{k});
        p = polyfit(log(N_store),log(err_store),1);
        rate = p(1);
        legend_str{k} = sprintf('slope=%4.1f',abs(rate));
    else
        err_store = zeros(1,Nk-1);
        for kk = 1:(Nk - 1)
            err_store(kk) = max(abs(y_store(kk+1)-y_store(kk)));
        end
        loglog(N_store(1:end-1),err_store,plot_str{k});
        p = polyfit(log(N_store(1:end-1)),log(err_store),1);
        rate = p(1);
        legend_str{k} = sprintf('slope=%4.1f',abs(rate));
    end
    xlabel('number of intervals');
    ylabel('absolute error');
    title('CDC, FE, Gauss--Legendre')
    set(gca,'FontSize',24)
    hold on
end
legend(legend_str,'Location','NorthEastOutside');


%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('CDC, FE, chebychev');

opts.dc = 1; % CDC
opts.grid = 3; % chebychev grid
opts.pred = 1; % FE
opts.corr = 1; % FE

tspan = [0 1]; % initial and final time
y0 = 1; % initial condition
opts.exact = 0.5; %exact solution if it exists

ti = 0;
tf = 1;

figure(3), clf
opts.nquad = 6; % number of quadrature points
for k = 1:6
    opts.levels = k; % number of levels (predictor + (levels-1)
                     % correctors

    % convergence study
    Nk = 10;
    m = length(y0);
    y_store = zeros(m,Nk);
    N_store = zeros(1,Nk);
    
    for kk = 1:Nk
        N = 2*kk+10;  % number of large (equi-spaced) time interavls 
        dt = (tf - ti)/N; % grid spacing, large intervals
        tspan = [0:N]*dt;
    
        sol = deferred_correction(@scalar_ode,tspan,y0,opts);
        y_store(:,kk) = sol(:,end);
        N_store(1,kk) = N;
    end
    
    
    if isfield(opts,'exact')
        % exact solution exists    
        err_store = zeros(1,Nk);
        for kk = 1:Nk 
            err_store(kk) = max(abs(y_store(:,kk)-opts.exact));
        end
        loglog(N_store,err_store,plot_str{k});
        p = polyfit(log(N_store),log(err_store),1);
        rate = p(1);
        legend_str{k} = sprintf('slope=%4.1f',abs(rate));
    else
        err_store = zeros(1,Nk-1);
        for kk = 1:(Nk - 1)
            err_store(kk) = max(abs(y_store(kk+1)-y_store(kk)));
        end
        loglog(N_store(1:end-1),err_store,plot_str{k});
        p = polyfit(log(N_store(1:end-1)),log(err_store),1);
        rate = p(1);
        legend_str{k} = sprintf('slope=%4.1f',abs(rate));
    end
    xlabel('number of intervals');
    ylabel('absolute error');
    title('CDC, FE, chebychev grid')
    set(gca,'FontSize',24)
    hold on
end
legend(legend_str,'Location','NorthEastOutside');


%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('CDC, FE, gauss lobatto');

opts.dc = 1; % CDC
opts.grid = 4; % gauss lobatto
opts.pred = 1; % FE
opts.corr = 1; % FE

tspan = [0 1]; % initial and final time
y0 = 1; % initial condition
opts.exact = 0.5; %exact solution if it exists

ti = 0;
tf = 1;

figure(4), clf
opts.nquad = 4; % number of quadrature points
for k = 1:6
    opts.levels = k; % number of levels (predictor + (levels-1)
                     % correctors

    % convergence study
    Nk = 10;
    m = length(y0);
    y_store = zeros(m,Nk);
    N_store = zeros(1,Nk);
    
    for kk = 1:Nk
        N = 2*kk+10;  % number of large (equi-spaced) time interavls 
        dt = (tf - ti)/N; % grid spacing, large intervals
        tspan = [0:N]*dt;
    
        sol = deferred_correction(@scalar_ode,tspan,y0,opts);
        y_store(:,kk) = sol(:,end);
        N_store(1,kk) = N;
    end
    
    
    if isfield(opts,'exact')
        % exact solution exists    
        err_store = zeros(1,Nk);
        for kk = 1:Nk 
            err_store(kk) = max(abs(y_store(:,kk)-opts.exact));
        end
        loglog(N_store,err_store,plot_str{k});
        p = polyfit(log(N_store),log(err_store),1);
        rate = p(1);
        legend_str{k} = sprintf('slope=%4.1f',abs(rate));
    else
        err_store = zeros(1,Nk-1);
        for kk = 1:(Nk - 1)
            err_store(kk) = max(abs(y_store(kk+1)-y_store(kk)));
        end
        loglog(N_store(1:end-1),err_store,plot_str{k});
        p = polyfit(log(N_store(1:end-1)),log(err_store),1);
        rate = p(1);
        legend_str{k} = sprintf('slope=%4.1f',abs(rate));
    end
    xlabel('number of intervals');
    ylabel('absolute error');
    title('CDC, FE, Gauss Lobatto')
    set(gca,'FontSize',24)
    hold on
end
legend(legend_str,'Location','NorthEastOutside');



%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('CDC-Midpoint, uniform grid');
opts.dc = 1; % CDC
opts.grid = 1; % uniform grid

opts.pred = 2; % RK2
opts.corr = 2; % RK2

tspan = [0 1]; % initial and final time
y0 = 1; % initial condition
opts.exact = 0.5; %exact solution if it exists

ti = 0;
tf = 1;

figure(5), clf
opts.nquad = 7; % number of quadrature points
legend_str={};
for k = 1:3
    opts.levels = k; % number of levels (predictor + (levels-1)
                     % correctors

    % convergence study
    Nk = 5;
    m = length(y0);
    y_store = zeros(m,Nk);
    N_store = zeros(1,Nk);
    
    for kk = 1:Nk
        N = 2*kk+2;  % number of large (equi-spaced) time interavls 
        dt = (tf - ti)/N; % grid spacing, large intervals
        tspan = [0:N]*dt;
    
        sol = deferred_correction(@scalar_ode,tspan,y0,opts);
        y_store(:,kk) = sol(:,end);
        N_store(1,kk) = N;
    end
    
    
    if isfield(opts,'exact')
        % exact solution exists    
        err_store = zeros(1,Nk);
        for kk = 1:Nk 
            err_store(kk) = max(abs(y_store(:,kk)-opts.exact));
        end
        loglog(N_store,err_store,plot_str{k});
        p = polyfit(log(N_store),log(err_store),1);
        rate = p(1);
        legend_str{k} = sprintf('slope=%4.1f',abs(rate));                
    else
        err_store = zeros(1,Nk-1);
        for kk = 1:(Nk - 1)
            err_store(kk) = max(abs(y_store(kk+1)-y_store(kk)));
        end
        loglog(N_store(1:end-1),err_store,plot_str{k});
        p = polyfit(log(N_store(1:end-1)),log(err_store),1);
        rate = p(1);
        legend_str{k} = sprintf('slope=%4.1f',abs(rate));
    end
    xlabel('number of intervals');
    ylabel('absolute error');
    title('CDC, Midpoint, uniform grid')
    set(gca,'FontSize',24)
    hold on
end
legend(legend_str,'Location','NorthEastOutside');


%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('CDC, midpoint, gauss lobatto');

opts.dc = 1; % CDC
opts.grid = 4; % gauss lobatto
opts.pred = 2; % RK2
opts.corr = 2; % RK@

tspan = [0 1]; % initial and final time
y0 = 1; % initial condition
opts.exact = 0.5; %exact solution if it exists

ti = 0;
tf = 1;

figure(6), clf
opts.nquad = 4; % number of quadrature points
legend_str={};
for k = 1:3
    opts.levels = k; % number of levels (predictor + (levels-1)
                     % correctors

    % convergence study
    Nk = 5;
    m = length(y0);
    y_store = zeros(m,Nk);
    N_store = zeros(1,Nk);
    
    for kk = 1:Nk
        N = 2*kk+2;  % number of large (equi-spaced) time interavls 
        dt = (tf - ti)/N; % grid spacing, large intervals
        tspan = [0:N]*dt;
    
        sol = deferred_correction(@scalar_ode,tspan,y0,opts);
        y_store(:,kk) = sol(:,end);
        N_store(1,kk) = N;
    end
    
    
    if isfield(opts,'exact')
        % exact solution exists    
        err_store = zeros(1,Nk);
        for kk = 1:Nk 
            err_store(kk) = max(abs(y_store(:,kk)-opts.exact));
        end
        loglog(N_store,err_store,plot_str{k});
        p = polyfit(log(N_store),log(err_store),1);
        rate = p(1);
        legend_str{k} = sprintf('slope=%4.1f',abs(rate));
    else
        err_store = zeros(1,Nk-1);
        for kk = 1:(Nk - 1)
            err_store(kk) = max(abs(y_store(kk+1)-y_store(kk)));
        end
        loglog(N_store(1:end-1),err_store,plot_str{k});
        p = polyfit(log(N_store(1:end-1)),log(err_store),1);
        rate = p(1);
        legend_str{k} = sprintf('slope=%4.1f',abs(rate));
    end
    xlabel('number of intervals');
    ylabel('absolute error');
    title('CDC, midpoint, Gauss-Lobatto')
    set(gca,'FontSize',24)
    hold on
end
legend(legend_str,'Location','NorthEastOutside');




%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('CDC-BE, uniform grid');
opts.dc = 1; % CDC
opts.grid = 1; % uniform grid
opts.pred = 6; % BE
opts.corr = 6; % BE

tspan = [0 1]; % initial and final time
y0 = 1; % initial condition
opts.exact = 0.5; %exact solution if it exists

ti = 0;
tf = 1;

figure(7), clf
for k = 1:6
    opts.levels = k; % number of levels (predictor + (levels-1)
                     % correctors
    opts.nquad = k+1; % number of quadrature points 

    % convergence study
    Nk = 8;
    m = length(y0);
    y_store = zeros(m,Nk);
    N_store = zeros(1,Nk);
    
    for kk = 1:Nk
        N = kk+10;  % number of large (equi-spaced) time interavls 
        dt = (tf - ti)/N; % grid spacing, large intervals
        tspan = [0:N]*dt;
    
        sol = deferred_correction(@scalar_ode,tspan,y0,opts);
        y_store(:,kk) = sol(:,end);
        N_store(1,kk) = N;
    end
    
    
    if isfield(opts,'exact')
        % exact solution exists    
        err_store = zeros(1,Nk);
        for kk = 1:Nk 
            err_store(kk) = max(abs(y_store(:,kk)-opts.exact));
        end
        loglog(N_store,err_store,plot_str{k});
        p = polyfit(log(N_store),log(err_store),1);
        rate = p(1);
        legend_str{k} = sprintf('slope=%4.1f',abs(rate));                
    else
        err_store = zeros(1,Nk-1);
        for kk = 1:(Nk - 1)
            err_store(kk) = max(abs(y_store(kk+1)-y_store(kk)));
        end
        loglog(N_store(1:end-1),err_store,plot_str{k});
        p = polyfit(log(N_store(1:end-1)),log(err_store),1);
        rate = p(1);
        legend_str{k} = sprintf('slope=%4.1f',abs(rate));
    end
    xlabel('number of intervals');
    ylabel('absolute error');
    title('CDC, BE, uniform grid')
    set(gca,'FontSize',24)
    hold on
end
legend(legend_str,'Location','NorthEastOutside');


%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('CDC-Trapezoid, uniform grid');
opts.dc = 1; % CDC
opts.grid = 1; % uniform grid
opts.pred = 3; % Trapezoid
opts.corr = 3; % Trapezoid

tspan = [0 1]; % initial and final time
y0 = 1; % initial condition
opts.exact = 0.5; %exact solution if it exists

ti = 0;
tf = 1;

legend_str={}
figure(8), clf
for k = 1:3
    opts.levels = k; % number of levels (predictor + (levels-1)
                     % correctors
    opts.nquad = 2*k+1; % number of quadrature points 

    % convergence study
    Nk = 8;
    m = length(y0);
    y_store = zeros(m,Nk);
    N_store = zeros(1,Nk);
    
    for kk = 1:Nk
        N = kk+5;  % number of large (equi-spaced) time interavls 
        dt = (tf - ti)/N; % grid spacing, large intervals
        tspan = [0:N]*dt;
    
        sol = deferred_correction(@scalar_ode,tspan,y0,opts);
        y_store(:,kk) = sol(:,end);
        N_store(1,kk) = N;
    end
    
    
    if isfield(opts,'exact')
        % exact solution exists    
        err_store = zeros(1,Nk);
        for kk = 1:Nk 
            err_store(kk) = max(abs(y_store(:,kk)-opts.exact));
        end
        loglog(N_store,err_store,plot_str{k});
        p = polyfit(log(N_store),log(err_store),1);
        rate = p(1);
        legend_str{k} = sprintf('slope=%4.1f',abs(rate));                
    else
        err_store = zeros(1,Nk-1);
        for kk = 1:(Nk - 1)
            err_store(kk) = max(abs(y_store(kk+1)-y_store(kk)));
        end
        loglog(N_store(1:end-1),err_store,plot_str{k});
        p = polyfit(log(N_store(1:end-1)),log(err_store),1);
        rate = p(1);
        legend_str{k} = sprintf('slope=%4.1f',abs(rate));
    end
    xlabel('number of intervals');
    ylabel('absolute error');
    title('CDC, Trapezoid, uniform grid')
    set(gca,'FontSize',24)
    hold on
end
legend(legend_str,'Location','NorthEastOutside');


%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('CDC-Trapezoid, uniform grid');
opts.dc = 1; % CDC
opts.grid = 4; % gauss lobatto
opts.pred = 3; % Trapezoid
opts.corr = 3; % Trapezoid

tspan = [0 1]; % initial and final time
y0 = 1; % initial condition
opts.exact = 0.5; %exact solution if it exists

ti = 0;
tf = 1;

legend_str={}
figure(9), clf
for k = 1:3
    opts.levels = k; % number of levels (predictor + (levels-1)
                     % correctors
    opts.nquad = 5; % number of quadrature points 

    % convergence study
    Nk = 8;
    m = length(y0);
    y_store = zeros(m,Nk);
    N_store = zeros(1,Nk);
    
    for kk = 1:Nk
        N = kk+5;  % number of large (equi-spaced) time interavls 
        dt = (tf - ti)/N; % grid spacing, large intervals
        tspan = [0:N]*dt;
    
        sol = deferred_correction(@scalar_ode,tspan,y0,opts);
        y_store(:,kk) = sol(:,end);
        N_store(1,kk) = N;
    end
    
    
    if isfield(opts,'exact')
        % exact solution exists    
        err_store = zeros(1,Nk);
        for kk = 1:Nk 
            err_store(kk) = max(abs(y_store(:,kk)-opts.exact));
        end
        loglog(N_store,err_store,plot_str{k});
        p = polyfit(log(N_store),log(err_store),1);
        rate = p(1);
        legend_str{k} = sprintf('slope=%4.1f',abs(rate));                
    else
        err_store = zeros(1,Nk-1);
        for kk = 1:(Nk - 1)
            err_store(kk) = max(abs(y_store(kk+1)-y_store(kk)));
        end
        loglog(N_store(1:end-1),err_store,plot_str{k});
        p = polyfit(log(N_store(1:end-1)),log(err_store),1);
        rate = p(1);
        legend_str{k} = sprintf('slope=%4.1f',abs(rate));
    end
    xlabel('number of intervals');
    ylabel('absolute error');
    title('CDC, Trapezoid, gauss lobatto nodes')
    set(gca,'FontSize',24)
    hold on
end
legend(legend_str,'Location','NorthEastOutside');
