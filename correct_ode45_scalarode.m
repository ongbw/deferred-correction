% testing idc corrector over top of ode45
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
disp('ode45 - FE correctors, adaptive mesh');

% initial conditions
y0 = [1];
   
tspan = [0 50]; % initial and final time

% ode 45 as predictor
predictor = 2;

if predictor == 1

    % use default options
    sol = ode45(@scalar_ode,tspan,y0);
    
    % adaptive time steps chosen by solver
    t_full = sol.x;
    % provisional solution at each of these time steps
    y_full = sol.y;
else

    % FE as predictor
    t_full = linspace(tspan(1),tspan(2),201);

    h = t_full(2:end)-t_full(1:end-1);
    
    y_full(:,1) = y0;
    for n = 2:length(t_full)
        y_full(:,n) = y_full(:,n-1) + ...
            h(n-1)*scalar_ode(t_full(n-1),y_full(:,n-1));
    end
end


%Etot = sum(energy);

figure(1); clf
plot(t_full,y_full,'k-');
xlabel('time');
ylabel('solution');


figure(2); clf
semilogy(t_full,abs(y_full - 1./(1+t_full.^2)),'k-');
xlabel('time');
ylabel('error');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perform a set number of idc FE correction steps
ncorr = 3; % 
nquad = 1 + ncorr; % RK4 + ncorr correctors.  Note because of
                   % adaptive mesh, no reason to use higher order
                   % integrator, unless possibly stability.

if length(t_full) < nquad
    error(['not enough quadrature points to obtain desired ' ...
           'accuracy']);
end

% precompute dydt from the predictor
dydt_full = zeros( size(y_full) );

for nt = 1:length(t_full)
    dydt_full(:,nt) = scalar_ode(t_full(nt),y_full(:,nt));
end
    
n_subintervals = nquad - 1;

% N: number of large intervals
N = ceil( (length(t_full)-1) / n_subintervals); 

ind_start = 1;
ind_end = nquad;
yi = y0;

% todo: maybe splitting equally like this isn't the right thing to
% do.  maybe it is better to do a la ridc, and one correction at a
% time, picking the most reasonable mesh


for nintervals = 1:(N-1)
    f{1} = dydt_full(:,ind_start:ind_end);
    t = t_full(ind_start:ind_end);
    h = (t(2:end)-t(1:end-1));
    
    % this rescaling should fix the roundoff errors in the
    % quadrature weights
    dt = t(end)-t(1);
    tt = (t - t(1))/dt;
    S = generate_integration_matrix(tt,tt);
    S = dt*S;
    
    for corr = 1:ncorr
        level = 1 + corr;
        y{level}(:,1) = yi;
        for n = 2:nquad
            y{level}(:,n) = y{level}(:,n-1) + ...
                h(n-1)*scalar_ode(t(n-1),y{level}(:,n-1)) - ...
                h(n-1)*f{level-1}(:,n-1) + ...
                f{level-1}*(S(n,:)' - S(n-1,:)'); 
            %scalar_ode(t(n-1),y{level}(:,n-1))
            %f{level-1}(:,n-1) 
            %f{level-1}
        end
        
        if corr < ncorr
            for n = 1:nquad
                f{level}(:,n) = scalar_ode(t(n),y{level}(:,n));
            end
        end

        figure(1);
        hold on
        plot(t,y{level},plot_str{mod(nintervals,8)+1})
        
        %[y{level};
        % abs(y{level}-1./(1+t.^2))]

        % pause
        
        figure(2); hold on
        semilogy(t,abs(y{level} - 1./(1+t.^2)),plot_str{mod(nintervals,8)+1})
    

    end
    %pause
    
    ind_start = ind_end;
    ind_end = ind_start + nquad - 1;
    yi = y{end}(:,end);

end

