% testing idc corrector over top of ode45
clear;


%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('ode45 - FE correctors, adaptive mesh');

% initial conditions
q0 = [1 0 0 1/100 0 0];
p0 = [1 0 0 1 0 0];
y0 = [q0'; p0'];
   
%tspan = [0 500]; % initial and final time
tspan = [0 1]; % initial and final time

% use default options
sol = ode45(@fpu_ode,tspan,y0);

% adaptive time steps chosen by solver
t_full = sol.x;
% provisional solution at each of these time steps
y_full = sol.y;

% show non-conservation of energy
energy = 0.5 * (y_full(10:12,:).^2 +  100^2 * y_full(4:6,:).^2);

Etot = sum(energy);

figure(1);
plot(t_full,Etot,'c-',...
     t_full,energy(1,:),'b-',...
     t_full,energy(2,:),'g-',...
     t_full,energy(3,:),'r-')

drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perform a set number of idc FE correction steps
ncorr = 1; % 
nquad = 4 + ncorr; % RK4 + ncorr correctors.  Note because of
                   % adaptive mesh, no reason to use higher order
                   % integrator, unless possibly stability.

if length(t_full) < nquad
    error(['not enough quadrature points to obtain desired ' ...
           'accuracy']);
end

% precompute dydt from the predictor
dydt_full = zeros( size(y_full) );

for nt = 1:length(t_full)
    dydt_full(:,nt) = fpu_ode(t_full(nt),y_full(:,nt));
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

figure(2), clf;

for nintervals = 1:(N-1)
    f{1} = dydt_full(:,ind_start:ind_end);
    t = t_full(ind_start:ind_end);
    h = (t(2:end)-t(1:end-1));
    
    S = generate_integration_matrix(t,t);
    
    for corr = 1:ncorr
        level = 1 + corr;
        y{level}(:,1) = yi;
        for n = 2:nquad
            y{level}(:,n) = y{level}(:,n-1) + ...
                h(n-1)*fpu_ode(t(n-1),y{level}(:,n-1)) - ...
                h(n-1)*f{level-1}(:,n-1) + ...
                f{level-1}*(S(n,:)' - S(n-1,:)'); 
        end
        
        if corr < ncorr
            for n = 1:nquad
                f{level}(:,n) = fpu_ode(t(n),y{level}(:,n));
            end
        end
    end

    % show non-conservation of energy
    energy = 0.5 * (y{end}(10:12,:).^2 +  100^2 * y{end}(4:6,:).^2);

    Etot = sum(energy);
    
    figure(2);
    plot(t,Etot,'c-',...
         t,energy(1,:),'b-',...
         t,energy(2,:),'g-',...
         t,energy(3,:),'r-')
    hold on
    
    ind_start = ind_end;
    ind_end = ind_start + nquad - 1;
    yi = y{end}(:,n);

    %pause
end

% [todo]: figure out last large interval.


    
    
    

% we'll use the unstable approach to compute quadrature weights
% first...  might need to be more careful with scaling.

