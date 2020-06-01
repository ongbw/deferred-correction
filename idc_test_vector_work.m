% testing idc in code
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
              exp(tf)*(cos(tf^2/2)-sin(tf^2/2))]; 


%ode = @scalar_ode ;
%tspan = [0 1]; % initial and final time
%y0 = 1; % initial condition
%opts.exact = 0.5; %exact solution if it exists

ti = tspan(1);
tf = tspan(2);

WRITETOFILE = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this generates figure 6 (old SIREV figure 4.4)

disp('IDC-FE, uniform grid');
% performs as expected

opts.dc = 2; % SDC/IDC
opts.grid = 1; % uniform grid
opts.pred = 1; % FE
opts.corr = 1; % FE

figure(1), clf
for k = 1:5
    opts.nquad = max(2,k); % number of quadrature points
    opts.ncorr = k-1; % number of correctors
    opts.levels = k; % number of levels (predictor + (levels-1)
                     % correctors

    % convergence study
    approx_func_evals = [20:40:320]';
    Nk = length(approx_func_evals);
    
    % size of system
    m = length(y0);
    y_store = zeros(m,Nk);
    N_store = zeros(1,Nk);
    work_store = zeros(1,Nk);
    
    feval_per_interval = k*(opts.nquad - 1);
    
     for kk = 1:Nk
        %interavls 
        N = floor(approx_func_evals(kk)/feval_per_interval);
        dt = (tf - ti)/N; % grid spacing, large intervals         
        tspan = [0:N]*dt;
    
        
        sol = deferred_correction(ode,tspan,y0,opts);
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
        ind = (work_store>0);
        p = polyfit(log(work_store(ind)),log(err_store(ind)),1);
        loglog(work_store,err_store,plot_str{k});
        rate = abs(p(1));
        legend_str{k} = sprintf('slope=%4.1f',abs(rate));                
    else
        err_store = zeros(1,Nk-1);
        for kk = 1:(Nk - 1)
            err_store(kk) = max(abs(y_store(kk+1)-y_store(kk)));
        end
        loglog(work_store(1:end-1),err_store,plot_str{k});
        p = polyfit(log(work_store(1:end-1)),log(err_store),1);
        rate = abs(p(1));
        legend_str{k} = sprintf('slope=%4.1f',abs(rate));
    end
    xlabel('number of function evaluations');
    ylabel('absolute error');
    title('IDC, FE, uniform grid')
    set(gca,'FontSize',24)
    hold on

    if WRITETOFILE
        % output data to file
        filename=['idc-L',num2str(k),'-unif-fe-fe-work.dat']
        fid = fopen(filename,'w');
        fprintf(fid,'%% Levels = %d\n',k);
        fprintf(fid,'\\addplot coordinates {\n');
        for kk = 1:(Nk-1)
            fprintf(fid,'  (%d,%g)\n',work_store(kk),err_store(kk));
        end
        if isfield(opts,'exact')
            fprintf(fid,'  (%d,%g)\n',work_store(Nk),err_store(Nk));
        end
        fprintf(fid,'};\n');
        if (k==1) % provisional solution
            fprintf(fid,'\\addlegendentry{Provisional \\\\ slope = %3.1f};\n\n',rate);
        elseif (k==2)
            fprintf(fid,'\\addlegendentry{1 correction \\\\ slope = %3.1f};\n\n',rate);
        else
            fprintf(fid,'\\addlegendentry{%d corrections \\\\ slope = %3.1f};\n\n',k-1,rate);    
        end
        fclose(fid);
    end
    
end
legend(legend_str,'Location','NorthEastOutside');



%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this generates figure 7 (old SIREV figure 4.5)

disp('IDC-RK2 predictor, RK2 corrector, uniform grid');
% performs as expected

% fig 4.5
opts.dc = 2; % SDC/IDC 
opts.grid = 1; % uniform grid
opts.pred = 2; % RK2
opts.corr = 2; % RK2

figure(2), clf
legend_str = {};

for k = 1:3
    opts.ncorr = k-1; % number of correctors
    opts.levels = k; % number of levels (predictor + (levels-1)
                     % correctors

    if k == 1
        opts.nquad = 2;
    else
        opts.nquad = 2 + 2*opts.ncorr;
    end

    % convergence study
    approx_func_evals = [20:20:320]';
    Nk = length(approx_func_evals);
    
    % size of system
    m = length(y0);
    y_store = zeros(m,Nk);
    N_store = zeros(1,Nk);
    work_store = zeros(1,Nk);
    
    feval_per_interval = 2*(opts.nquad-1) + ...% predictor
        opts.ncorr * 2 * (opts.nquad-1);
    
    
     for kk = 1:Nk
        %interavls 
        N = floor(approx_func_evals(kk)/feval_per_interval);
        dt = (tf - ti)/N; % grid spacing, large intervals         
        tspan = [0:N]*dt;
    
        
        sol = deferred_correction(ode,tspan,y0,opts);
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
        ind = (work_store>0);
        p = polyfit(log(work_store(ind)),log(err_store(ind)),1);
        loglog(work_store,err_store,plot_str{k});
        rate = abs(p(1));
        legend_str{k} = sprintf('slope=%4.1f',abs(rate));                
    else
        err_store = zeros(1,Nk-1);
        for kk = 1:(Nk - 1)
            err_store(kk) = max(abs(y_store(kk+1)-y_store(kk)));
        end
        loglog(work_store(1:end-1),err_store,plot_str{k});
        p = polyfit(log(work_store(1:end-1)),log(err_store),1);
        rate = abs(p(1));
        legend_str{k} = sprintf('slope=%4.1f',abs(rate));
    end
    xlabel('number of function evaluations');
    ylabel('absolute error');
    title('IDC, RK2 predictor RK2 correctors, uniform grid')
    set(gca,'FontSize',24)
    hold on

    if WRITETOFILE
        % output data to file
        filename=['idc-L',num2str(k),'-unif-rk2-rk2-work.dat']
        fid = fopen(filename,'w');
        fprintf(fid,'%% Levels = %d\n',k);
        fprintf(fid,'\\addplot coordinates {\n');
        for kk = 1:(Nk-1)
            fprintf(fid,'  (%d,%g)\n',work_store(kk),err_store(kk));
        end
        if isfield(opts,'exact')
            fprintf(fid,'  (%d,%g)\n',work_store(Nk),err_store(Nk));
        end
        fprintf(fid,'};\n');
        if (k==1) % provisional solution
            fprintf(fid,'\\addlegendentry{Provisional \\\\ slope = %3.1f};\n\n',rate);
        elseif (k==2)
            fprintf(fid,'\\addlegendentry{1 correction \\\\ slope = %3.1f};\n\n',rate);
        else
            fprintf(fid,'\\addlegendentry{%d corrections \\\\ slope = %3.1f};\n\n',k-1,rate);    
        end
        fclose(fid);
    end

end
legend(legend_str,'Location','NorthEastOutside');


%%%%%%%%%%%%%%%%%%%%%%%%%%
% this generates figure 10 (old SIREV figure 4.8)

disp('SDC-FE predictor, FE corrector, lobatto');
% performs as expected

opts.dc = 2; % SDC/IDC 
opts.grid = 4; % lobatto
opts.pred = 1; % FE
opts.corr = 1; % FE

figure(3), clf
legend_str = {};

for k = 1:6
    opts.ncorr = k-1; % number of correctors
    opts.levels = k; % number of levels (predictor + (levels-1)
                     % correctors

    if k == 1
        opts.nquad = 2;
    else
        opts.nquad = ceil( (k+2)/2);
    end

    % convergence study
    approx_func_evals = [20:20:320]';
    Nk = length(approx_func_evals);
    
    % size of system
    m = length(y0);
    y_store = zeros(m,Nk);
    N_store = zeros(1,Nk);
    work_store = zeros(1,Nk);
    
    feval_per_interval = (opts.nquad-1) + ...% predictor
        opts.ncorr *  (opts.nquad-1);
    
    
     for kk = 1:Nk
        %interavls 
        N = floor(approx_func_evals(kk)/feval_per_interval);
        dt = (tf - ti)/N; % grid spacing, large intervals         
        tspan = [0:N]*dt;
    
        
        sol = deferred_correction(ode,tspan,y0,opts);
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
        ind = (work_store>0);
        p = polyfit(log(work_store(ind)),log(err_store(ind)),1);
        loglog(work_store,err_store,plot_str{k});
        rate = abs(p(1));
        legend_str{k} = sprintf('slope=%4.1f',abs(rate));                
    else
        err_store = zeros(1,Nk-1);
        for kk = 1:(Nk - 1)
            err_store(kk) = max(abs(y_store(kk+1)-y_store(kk)));
        end
        loglog(work_store(1:end-1),err_store,plot_str{k});
        p = polyfit(log(work_store(1:end-1)),log(err_store),1);
        rate = abs(p(1));
        legend_str{k} = sprintf('slope=%4.1f',abs(rate));
    end
    xlabel('number of function evaluations');
    ylabel('absolute error');
    title('SDC, FE predictor FE correctors, lobatto')
    set(gca,'FontSize',24)
    hold on

    if WRITETOFILE
        % output data to file
        filename=['sdc-L',num2str(k),'-lob-fe-fe-work.dat']
        fid = fopen(filename,'w');
        fprintf(fid,'%% Levels = %d\n',k);
        fprintf(fid,'\\addplot coordinates {\n');
        for kk = 1:(Nk-1)
            fprintf(fid,'  (%d,%g)\n',work_store(kk),err_store(kk));
        end
        if isfield(opts,'exact')
            fprintf(fid,'  (%d,%g)\n',work_store(Nk),err_store(Nk));
        end
        fprintf(fid,'};\n');
        if (k==1) % provisional solution
            fprintf(fid,'\\addlegendentry{Provisional \\\\ slope = %3.1f};\n\n',rate);
        elseif (k==2)
            fprintf(fid,'\\addlegendentry{1 correction \\\\ slope = %3.1f};\n\n',rate);
        else
            fprintf(fid,'\\addlegendentry{%d corrections \\\\ slope = %3.1f};\n\n',k-1,rate);    
        end
        fclose(fid);
    end

end
legend(legend_str,'Location','NorthEastOutside');




%%%%%%%%%%%%%%%%%%%%%%%%%%
% this generate figure 11 (old SIREV figure 4.9)

disp('SDC-FE predictor, FE corrector, legendre');
% performs as expected

opts.dc = 3; % SDC
opts.grid = 2; % legendre nodes


opts.pred = 1; 
opts.corr = 1; 

figure(4), clf
legend_str = {};



expected_order = [2 4 6 8];
quad = [1 2 3 4];
nlevels = [2 4 6 8];

opts.interp = 0; % generate collocation solution

for counter = 1:length(nlevels)
    k = nlevels(counter);
    opts.ncorr = k-1; % number of correctors
    opts.levels = k; % number of levels (predictor + (levels-1)
                     % correctors

    opts.nquad = quad(counter);

    % convergence study
    approx_func_evals = [40:20:320]';
    Nk = length(approx_func_evals);
        
    
    % size of system
    m = length(y0);
    y_store = zeros(m,Nk);
    N_store = zeros(1,Nk);
    work_store = zeros(1,Nk);

    feval_per_interval = (opts.levels+1)*(opts.nquad);
    % note: there is a +1 in opts.level because we are doing one
    % final collocation solve.  
    
    for kk = 1:Nk
        %interavls 
        N = floor(approx_func_evals(kk)/feval_per_interval);
        dt = (tf - ti)/N; % grid spacing, large intervals         
        tspan = [0:N]*dt;
            
        sol = deferred_correction(ode,tspan,y0,opts);
        y_store(:,kk) = sol(:,end);
        N_store(1,kk) = N;
        work_store(1,kk) =  N*feval_per_interval;
     end
    
    if isfield(opts,'exact')
        % exact solution exists    
        err_store = zeros(1,Nk);
        for kk = 1:Nk 
            err_store(kk) = max(abs(y_store(:,kk)-opts.exact));
        end
        p = polyfit(log(work_store),log(err_store),1);
        loglog(work_store,err_store,plot_str{k});
        rate = abs(p(1));
        legend_str{counter} = sprintf('slope=%4.1f',abs(rate));                
    else
        err_store = zeros(1,Nk-1);
        for kk = 1:(Nk - 1)
            err_store(kk) = max(abs(y_store(kk+1)-y_store(kk)));
        end
        loglog(work_store(1:end-1),err_store,plot_str{k});
        p = polyfit(log(work_store(1:end-1)),log(err_store),1);
        rate = abs(p(1));
        legend_str{counter} = sprintf('slope=%4.1f',abs(rate));
    end
    xlabel('number of function evaluations');
    ylabel('absolute error');
    title('SDC, FE predictor FE correctors, legendre')
    set(gca,'FontSize',24)
    hold on

    if WRITETOFILE
        % output data to file
        filename=['sdc-L',num2str(k),'-leg-fe-fe-work.dat']
        fid = fopen(filename,'w');
        fprintf(fid,'%% Levels = %d\n',k);
        fprintf(fid,'\\addplot coordinates {\n');
        for kk = 1:(Nk-1)
            fprintf(fid,'  (%d,%g)\n',work_store(kk),err_store(kk));
        end
        if isfield(opts,'exact')
            fprintf(fid,'  (%d,%g)\n',work_store(Nk),err_store(Nk));
        end
        fprintf(fid,'};\n');
        if (k==1) % provisional solution
            fprintf(fid,'\\addlegendentry{Provisional \\\\ slope = %3.1f};\n\n',rate);
        elseif (k==2)
            fprintf(fid,'\\addlegendentry{1 correction \\\\ slope = %3.1f};\n\n',rate);
        else
            fprintf(fid,'\\addlegendentry{%d corrections \\\\ slope = %3.1f};\n\n',k-1,rate);    
        end
        fclose(fid);
    end

end
legend(legend_str,'Location','NorthEastOutside');


                


error('break')

%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('SDC, Gauss-Legendre');
% order of accuracy limited by interpolation, even though
% quadrature sufficiently accurate

opts.dc = 2; % SDC/IDC 
opts.grid = 2; % gaussian grid
opts.pred = 1; % FE
opts.corr = 1; % FE


figure(2), clf

for k = 1:6
    opts.nquad = 2*ceil(k/2); % number of quadrature points needed
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

opts.dc = 2; % SDC/IDC 
opts.grid = 3; % chebychev grid
opts.pred = 1; % FE
opts.corr = 1; % FE

figure(3), clf

for k = 1:6
    opts.nquad = max(2,k); % number of quadrature points

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

opts.dc = 2; % SDC/IDC 
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

opts.dc = 2; % SDC/IDC 
opts.grid = 1; % uniform grid
opts.pred = 2; % RK2
opts.corr = 2; % RK2


tspan = [0 1]; % initial and final time
y0 = 1; % initial condition
opts.exact = 0.5; %exact solution if it exists

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

opts.dc = 2; % SDC/IDC 
opts.grid = 4; % gauss lobatto
opts.pred = 3; % RK2
opts.corr = 3; % RK2

figure(6), clf
legend_str={};


expected_order = [2 4 6];
%quad = [ 2 3 5]; % max order = 2*quad-3
quad = [ 2 3 7]; % max order = 2*quad-3

for k = 1:3
    opts.nquad = quad(k);
    %    opts.nquad = max(2,ceil((2*k+3)/2)); % number of quadrature
                                         % points
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

opts.dc = 2; % SDC/IDC 
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

opts.dc = 2; % SDC/IDC 
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

