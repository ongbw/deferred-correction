% testing idc in code on nbody problem
clear;

plot_str{1} = 'ko-'; % 1st order
plot_str{2} = 'b^-'; % 2nd
plot_str{3} = 'r*-'; % 3rd
plot_str{4} = 'cx-'; % 4th
plot_str{5} = 'gd-'; % 5th
plot_str{6} = 'ys-'; % 6th
plot_str{7} = 'k+--'; % 7th
plot_str{8} = 'b<--'; %8th

ode = @nbody ;
tspan = [0 1]; % initial and final time

% initial condition
nplus = 200;
nminus = 200;

xip = ([1:nplus]'-0.5)/nplus;
xim = ([1:nminus]'-0.5)/nminus;

y0 = [ xip;
       zeros(nplus,1);...
       xim;
       sin(6*pi*xim)];
       

ti = tspan(1);
tf = tspan(2);

% testing sdc/idc
opts.dc = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('IDC-FE, uniform grid');
% performs as expected

opts.grid = 1; % uniform grid
opts.pred = 1; % FE
opts.corr = 1; % FE

figure(1), clf
for k = 1:1
    opts.nquad = max(2,k); % number of quadrature points

    opts.levels = k; % number of levels (predictor + (levels-1)
                     % correctors

    % convergence study
    Nk = 4;
    m = length(y0);
    y_store = zeros(m,Nk);
    N_store = zeros(1,Nk);
    
    for kk = 1:Nk
        %N = 2*kk+10;  % number of large (equi-spaced) time
        %interavls 
        N = 1000*2^(kk-1);
        dt = (tf - ti)/N; % grid spacing, large intervals
        tspan = [0:N]*dt;
    
        sol = deferred_correction(ode,tspan,y0,opts);
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
    title('IDC, FE, uniform grid')
    set(gca,'FontSize',24)
    hold on
end
legend(legend_str,'Location','NorthEastOutside');


error('break');

%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('SDC, Gauss-Legendre');
% order of accuracy limited by interpolation, even though
% quadrature sufficiently accurate

opts.grid = 2; % gaussian grid
opts.pred = 1; % FE
opts.corr = 1; % FE


figure(2), clf

%order = [ 1 2 3 4 5 6];
%quad = [ 2 2 2 
for k = 1:6
    opts.nquad = 7; % number of quadrature points needed
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
    
        sol = deferred_correction(ode,tspan,y0,opts);
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
    title('SDC, Gauss--Legendre')
    set(gca,'FontSize',24)
    hold on
end
legend(legend_str,'Location','NorthEastOutside');


%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('IDC-FE, chebychev');
% we lose one order because the polynomial has to be interpolated

opts.grid = 3; % chebychev grid
opts.pred = 1; % FE
opts.corr = 1; % FE

figure(3), clf

for k = 1:6
    opts.nquad = max(2,k) + 1; % number of quadrature points

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
    
        sol = deferred_correction(ode,tspan,y0,opts);
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
    title('IDC, FE, chebychev grid')
    set(gca,'FontSize',24)
    hold on
end
legend(legend_str,'Location','NorthEastOutside');


%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('SDC, FE, gauss lobatto');
% one order increase each time observed as expected

opts.grid = 4; % gauss lobatto
opts.pred = 1; % FE
opts.corr = 1; % FE


figure(4), clf

for k = 1:6
    opts.nquad = max(2,ceil((k+3)/2)); % number of quadrature points
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
    
        sol = deferred_correction(ode,tspan,y0,opts);
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
    title('SDC, FE, gauss lobatto')
    set(gca,'FontSize',24)
    hold on
end
legend(legend_str,'Location','NorthEastOutside');



%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('IDC-Midpoint, uniform grid');
% two order increase each correction observed.  behaves as expected

opts.grid = 1; % uniform grid
opts.pred = 2; % RK2
opts.corr = 2; % RK2


tspan = [0 1]; % initial and final time

figure(5), clf
legend_str={};
for k = 1:3
    opts.nquad = 2*k; % number of quadrature points
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
    
        sol = deferred_correction(ode,tspan,y0,opts);
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
    title('IDC-Midpoint, uniform grid')
    set(gca,'FontSize',24)
    hold on
end
legend(legend_str,'Location','NorthEastOutside');


%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('SDC, midpoint, gauss lobatto');
% order 2 increase seen for each corrector.  not expected.
% previously, took incorrect number of quadrature nodes? fixed at
% n=4.  that should have provided at (2n-3) 5th order accurate interpolant

opts.grid = 4; % gauss lobatto
opts.pred = 2; % RK2
opts.corr = 2; % RK2

figure(6), clf
legend_str={};
for k = 1:3
    opts.nquad = max(2,ceil((2*k+3)/2)); % number of quadrature
                                         % points
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
    
        sol = deferred_correction(ode,tspan,y0,opts);
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
    title('SDC, midpoint, Gauss-Lobatto')
    set(gca,'FontSize',24)
    hold on
end
legend(legend_str,'Location','NorthEastOutside');



%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('IDC, RK3 predictor, RK2 corrector');

opts.grid = 1; % uniform
opts.pred = 4; % RK3
opts.corr = 2; % RK2

figure(7), clf

legend_str={};
expected_order = [3 5 7];
for k = 1:3
    
    opts.nquad = expected_order(k); % number of quadrature points
    
    opts.levels = k; % number of levels (predictor + (levels-1)
                     % correctors

    % convergence study
    Nk = 6;
    m = length(y0);
    y_store = zeros(m,Nk);
    N_store = zeros(1,Nk);
    
    for kk = 1:Nk
        N = kk + 2;  % number of large (equi-spaced) time interavls 
        dt = (tf - ti)/N; % grid spacing, large intervals
        tspan = [0:N]*dt;
    
        sol = deferred_correction(ode,tspan,y0,opts);
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
    title('IDC, RK3 predictor, RK2 corrector')
    set(gca,'FontSize',24)
    hold on
end
legend(legend_str,'Location','NorthEastOutside');


%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('IDC, IRK2 ');

opts.grid = 1; % uniform
opts.pred = 3; % Implicit RK2
opts.corr = 3; % implicit RK2

figure(8), clf

legend_str={};
expected_order = [2 4 6];
for k = 1:3
    
    opts.nquad = expected_order(k); % number of quadrature points
    
    opts.levels = k; % number of levels (predictor + (levels-1)
                     % correctors

    % convergence study
    Nk = 6;
    m = length(y0);
    y_store = zeros(m,Nk);
    N_store = zeros(1,Nk);
    
    for kk = 1:Nk
        N = kk + 2;  % number of large (equi-spaced) time interavls 
        dt = (tf - ti)/N; % grid spacing, large intervals
        tspan = [0:N]*dt;
    
        sol = deferred_correction(ode,tspan,y0,opts);
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
    title('IDC, Implicit RK2')
    set(gca,'FontSize',24)
    hold on
end
legend(legend_str,'Location','NorthEastOutside');


%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('SDC, trapezoid, gauss lobatto');
% order 2 increase seen for each corrector.  not expected.
% previously, took incorrect number of quadrature nodes? fixed at
% n=4.  that should have provided at (2n-3) 5th order accurate interpolant

opts.grid = 4; % gauss lobatto
opts.pred = 3; % RK2
opts.corr = 3; % RK2

figure(9), clf
legend_str={};
for k = 1:3
    opts.nquad = k+2;
    %opts.nquad = max(2,ceil((2*k+3)/2)); % number of quadrature
    %                                     % points
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
    
        sol = deferred_correction(ode,tspan,y0,opts);
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
    title('SDC, trapezoid, Gauss-Lobatto')
    set(gca,'FontSize',24)
    hold on
end
legend(legend_str,'Location','NorthEastOutside');
