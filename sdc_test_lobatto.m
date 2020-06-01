% testing idc in code on vector ode
clear;

plot_str{1} = 'ko-'; % 1st order
plot_str{2} = 'b^-'; % 2nd
plot_str{3} = 'r*-'; % 3rd
plot_str{4} = 'cx-'; % 4th
plot_str{5} = 'gd-'; % 5th
plot_str{6} = 'ys-'; % 6th
plot_str{7} = 'k+--'; % 7th
plot_str{8} = 'b<--'; %8th

ode = @vector_ode ;
tspan = [0 1]; % initial and final time
y0 = [1;1]; % initial condition
ti = tspan(1);
tf = tspan(2);
opts.exact = [exp(tf)*(cos(tf^2/2)+sin(tf^2/2));
              exp(tf)*(cos(tf^2/2)-sin(tf^2/2))]; %exact solution if it exists

% testing sdc/idc
opts.dc = 2;

PRINTOUTPUT=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('SDC, Gauss-Legendre, RK2');
% order of accuracy limited by interpolation, even though
% quadrature sufficiently accurate

opts.grid = 2; % gauss--legendre

figure(7), clf

opts.pred = 4;
opts.corr=4;

%opts.pred = 1; % FE
%opts.corr = 1; % FE

%opts.pred = 4; % RK3
%opts.corr = 4; % RK3

%opts.pred = 5; % RK4
%opts.corr = 5; % RK4


%opts.pred = 2; % RK2 - Midpoint
%opts.corr = 2; % RK2

expected_order = [ 2 4 6];

quad = [3];

%quad = [4 4 4 4 4 ];
%quad = [3 3 3 3 3 3];
%quad = [ 4 4 4 4 4 4 4 4]

%quad = [2 2 3];% order = 2*quad 
%quad = [2 4 6];

% quad = [3 5 7];
%quad = [4 6 8 ];

for k = 1:1
    opts.levels = k; % number of quadrature points needed
    opts.nquad = quad(k); % number of levels (predictor + (levels-1)
                     % correctors
    % convergence study
    Nk = 4;
    m = length(y0);
    y_store = zeros(m,Nk);
    N_store = zeros(1,Nk);
    
    for kk = 1:Nk
        N = 2^(kk+1);  % number of large (equi-spaced) time interavls 
        dt = (tf - ti)/N; % grid spacing, large intervals
        tspan = [0:N]*dt;
    
        sol = deferred_correction(ode,tspan,y0,opts);
        y_store(:,kk) = sol(:,end);
        N_store(1,kk) = N;
    end
    
    if PRINTOUTPUT,  fprintf('\n');    end
    if isfield(opts,'exact')
        % exact solution exists    
        err_store = zeros(1,Nk);
        for kk = 1:Nk 
            err_store(kk) = max(abs(y_store(:,kk)-opts.exact));
            if PRINTOUTPUT
                fprintf('(%d,%g)\n',N_store(kk),err_store(kk));
            end        
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
    title('SDC, Gauss--Legendre, RK2')
    set(gca,'FontSize',24)
    hold on
end
legend(legend_str,'Location','NorthEastOutside');

