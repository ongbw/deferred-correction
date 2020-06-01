% testing zadunaisky code
clear;

plot_str{1} = 'ko-'; % 1st order
plot_str{2} = 'b^-'; % 2nd
plot_str{3} = 'r*-'; % 3rd
plot_str{4} = 'cx-'; % 4th
plot_str{5} = 'gd-'; % 5th
plot_str{6} = 'ys-'; % 6th
plot_str{7} = 'k+--'; % 7th
plot_str{8} = 'b<--'; %8th

% ode = @vode2;
%tspan = [0 2]; % initial and final time
%y0 = vode2_exact(tspan(1)); % initial condition
%opts.exact = vode2_exact(tspan(end)); %exact solution if it exists

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

WRITETOFILE = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this generates figure 3 (old SIREV figure 4.1)

disp('Zadunaisky-FE, uniform grid');
opts.dc = 4; % zadunaisky
opts.grid = 1; % uniform grid
opts.pred = 1; % FE
opts.corr = 1; % FE

figure(1), clf
legend_str={};

for k = 1:5
    opts.nquad = max(2,k+1); % number of quadrature points

    opts.levels = k; % number of levels (predictor + (levels-1)
                     % correctors

    % convergence study
    Nk = 5;
    m = length(y0);
    y_store = zeros(m,Nk);
    N_store = zeros(1,Nk);
    
    for kk = 1:Nk

        N = 5*2^(kk);
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
        rate = abs(p(1));
        legend_str{k} = sprintf('slope=%4.1f',abs(rate));                
    else
        err_store = zeros(1,Nk-1);
        for kk = 1:(Nk - 1)
            err_store(kk) = max(abs(y_store(kk+1)-y_store(kk)));
        end
        loglog(N_store(1:end-1),err_store,plot_str{k});
        p = polyfit(log(N_store(1:end-1)),log(err_store),1);
        rate = abs(p(1));
        legend_str{k} = sprintf('slope=%4.1f',abs(rate));
    end
    xlabel('number of intervals');
    ylabel('absolute error');
    title('Zadunaisky, FE, uniform grid')
    set(gca,'FontSize',24)
    hold on
    
    if WRITETOFILE
        % output data to file
        filename=['zd-L',num2str(k),'-unif-fe-fe.dat']
        fid = fopen(filename,'w');
        fprintf(fid,'%% Levels = %d\n',k);
        fprintf(fid,'\\addplot coordinates {\n');
        for kk = 1:(Nk-1)
            fprintf(fid,'  (%d,%g)\n',N_store(kk),err_store(kk));
        end
        if isfield(opts,'exact')
            fprintf(fid,'  (%d,%g)\n',N_store(Nk),err_store(Nk));
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
% this generates figure 4 (old SIREV figure 4.2)

disp('Zadunaisky-FE, uniform grid');
opts.dc = 4; % zadunaisky
opts.grid = 1; % uniform grid
opts.pred = 1; % FE
opts.corr = 1; % FE

figure(2), clf
legend_str={};

for k = 1:5
    opts.nquad = max(2,k+1); % number of quadrature points
    opts.ncorr = k-1; % number of correctors
    opts.levels = k; % number of levels (predictor + (levels-1)
                     % correctors

    % convergence study
    approx_func_evals = [40:40:320]';
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
        loglog(work_store,err_store,plot_str{k});
        p = polyfit(log(work_store),log(err_store),1);
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
    title('Zadunaisky, FE, uniform grid')
    set(gca,'FontSize',24)
    hold on
    
    if WRITETOFILE
        % output data to file
        filename=['zd-L',num2str(k),'-unif-fe-fe-work.dat']
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
% this generates figure 5 (old SIREV figure 4.3)

disp('Zadunaisky-RK3-RK2, uniform grid');
opts.dc = 4; % zadunaisky
opts.grid = 1; % uniform grid
opts.pred = 4; % RK3
opts.corr = 2; % RK2

figure(3), clf
legend_str={};

    
for k = 1:3

    opts.ncorr = k-1; % number of correctors
    opts.levels = k; % number of levels (predictor + (levels-1)
                     % correctors
    
    if k == 1
        opts.nquad = 2;
    else
        opts.nquad = 4 + 2*opts.ncorr;
    end
    %opts.nquad = max(2,k); % number of quadrature points



    % convergence study
    approx_func_evals = [50:50:300]';
    Nk = length(approx_func_evals);
    
    % size of system
    m = length(y0);
    y_store = zeros(m,Nk);
    N_store = zeros(1,Nk);
    work_store = zeros(1,Nk);
    
    feval_per_interval = 3*(opts.nquad-1) + ...% predictor
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
        loglog(work_store,err_store,plot_str{k});
        p = polyfit(log(work_store),log(err_store),1);
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
    title('Zadunaisky, RK3-RK2, uniform grid')
    set(gca,'FontSize',24)
    hold on
    
    if WRITETOFILE
        % output data to file
        filename=['zd-L',num2str(k),'-unif-rk3-rk2-work.dat']
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
% this generates figure 8 (old SIREV figure 4.6)

disp('Zadunaisky-FE, gauss-lobatto grid');
opts.dc = 4; % zadunaisky
opts.grid = 4; % gauss-lobatto
opts.pred = 1; % FE
opts.corr = 1; % FE

figure(4), clf
legend_str={};

for k = 1:5
    opts.nquad = max(2,k+1); % number of quadrature points

    opts.levels = k; % number of levels (predictor + (levels-1)
                     % correctors

    % convergence study
    Nk = 5;
    m = length(y0);
    y_store = zeros(m,Nk);
    N_store = zeros(1,Nk);
    
    for kk = 1:Nk

        N = 5*2^(kk);
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
        rate = abs(p(1));
        legend_str{k} = sprintf('slope=%4.1f',abs(rate));                
    else
        err_store = zeros(1,Nk-1);
        for kk = 1:(Nk - 1)
            err_store(kk) = max(abs(y_store(kk+1)-y_store(kk)));
        end
        loglog(N_store(1:end-1),err_store,plot_str{k});
        p = polyfit(log(N_store(1:end-1)),log(err_store),1);
        rate = abs(p(1));
        legend_str{k} = sprintf('slope=%4.1f',abs(rate));
    end
    xlabel('number of intervals');
    ylabel('absolute error');
    title('Zadunaisky, FE, gauss-lobatto grid')
    set(gca,'FontSize',24)
    hold on
    
    if WRITETOFILE
        % output data to file
        filename=['zd-L',num2str(k),'-lob-fe-fe.dat']
        fid = fopen(filename,'w');
        fprintf(fid,'%% Levels = %d\n',k);
        fprintf(fid,'\\addplot coordinates {\n');
        for kk = 1:(Nk-1)
            fprintf(fid,'  (%d,%g)\n',N_store(kk),err_store(kk));
        end
        if isfield(opts,'exact')
            fprintf(fid,'  (%d,%g)\n',N_store(Nk),err_store(Nk));
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
% this generates figure 9 (old SIREV figure 4.7)

disp('Zadunaisky-RK2, gauss-lobatto grid');
opts.dc = 4; % zadunaisky
opts.grid = 4; % gauss-lobatto
opts.pred = 2; % RK2
opts.corr = 2; % RK2

figure(5), clf
legend_str={};

for k = 1:3
    opts.ncorr = k - 1;
    opts.levels = k; % number of levels (predictor + (levels-1)
                     % correctors
    if k == 1
        opts.nquad = 2;
    else
        opts.nquad = 3 + 2*opts.ncorr;
    end
    
    % convergence study
    Nk = 5;
    m = length(y0);
    y_store = zeros(m,Nk);
    N_store = zeros(1,Nk);
    
    for kk = 1:Nk

        N = 5*2^(kk);
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
        rate = abs(p(1));
        legend_str{k} = sprintf('slope=%4.1f',abs(rate));                
    else
        err_store = zeros(1,Nk-1);
        for kk = 1:(Nk - 1)
            err_store(kk) = max(abs(y_store(kk+1)-y_store(kk)));
        end
        loglog(N_store(1:end-1),err_store,plot_str{k});
        p = polyfit(log(N_store(1:end-1)),log(err_store),1);
        rate = abs(p(1));
        legend_str{k} = sprintf('slope=%4.1f',abs(rate));
    end
    xlabel('number of intervals');
    ylabel('absolute error');
    title('Zadunaisky, FE, gauss-lobatto grid')
    set(gca,'FontSize',24)
    hold on
    
    if WRITETOFILE
        % output data to file
        filename=['zd-L',num2str(k),'-lob-rk2-rk2.dat']
        fid = fopen(filename,'w');
        fprintf(fid,'%% Levels = %d\n',k);
        fprintf(fid,'\\addplot coordinates {\n');
        for kk = 1:(Nk-1)
            fprintf(fid,'  (%d,%g)\n',N_store(kk),err_store(kk));
        end
        if isfield(opts,'exact')
            fprintf(fid,'  (%d,%g)\n',N_store(Nk),err_store(Nk));
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




%% rk4
% disp('Zadunaisky-RK4, uniform grid');
% opts.dc = 4; % zadunaisky
% opts.grid = 1; % uniform grid
% opts.pred = 5; % RK4
% opts.corr = 5; % RK4


% figure(2), clf
% opts.nquad = 9; % number of quadrature points 
% for k = 1:2
%     opts.levels = k; % number of levels (predictor + (levels-1)
%                      % correctors

%     % convergence study
%     Nk = 5;
%     m = length(y0);
%     y_store = zeros(m,Nk);
%     N_store = zeros(1,Nk);
    
%     for kk = 1:Nk
%         N = kk+2;  % number of large (equi-spaced) time interavls 
%         dt = (tf - ti)/N; % grid spacing, large intervals
%         tspan = [0:N]*dt;
    
%         sol = deferred_correction(ode,tspan,y0,opts);
%         y_store(:,kk) = sol(:,end);
%         N_store(1,kk) = N;
%     end
    
    
%     if isfield(opts,'exact')
%         % exact solution exists    
%         err_store = zeros(1,Nk);
%         for kk = 1:Nk 
%             err_store(kk) = max(abs(y_store(:,kk)-opts.exact));
%         end
%         loglog(N_store,err_store,plot_str{k});
%         p = polyfit(log(N_store),log(err_store),1);
%         rate = p(1);
%         legend_str{k} = sprintf('slope=%4.1f',abs(rate));                
%     else
%         err_store = zeros(1,Nk-1);
%         for kk = 1:(Nk - 1)
%             err_store(kk) = max(abs(y_store(kk+1)-y_store(kk)));
%         end
%         loglog(N_store(1:end-1),err_store,plot_str{k});
%         p = polyfit(log(N_store(1:end-1)),log(err_store),1);
%         rate = p(1);
%         legend_str{k} = sprintf('slope=%4.1f',abs(rate));
%     end
%     xlabel('number of intervals');
%     ylabel('absolute error');
%     title('Zadunaisky, RK4, uniform grid')
%     set(gca,'FontSize',24)
%     hold on
% end
% legend(legend_str,'Location','NorthEastOutside');

