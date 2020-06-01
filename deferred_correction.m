function [sol] = deferred_correction(ode,tspan,y0,opts)


%% specify quadrature grid
% opts.grid = 1; % uniform
% opts.grid = 2; % gaussian nodes (Gauss Legendre)
% opts.grid = 3; % Chebysehv nodes (Gauss Legendre)

%% integrator building blocks.  
% opts.integrator = 1; 

[param.pred.A, param.pred.b, param.pred.c,...
 param.pred.s, param.pred.type, param.pred.order] = ...
    init_integrator(opts.pred);

[param.corr.A, param.corr.b, param.corr.c,...
 param.corr.s, param.corr.type, param.corr.order] = ...
    init_integrator(opts.corr);



N = length(tspan);
m = length(y0);
sol = zeros(m,N);
sol(:,1) = y0;

for n = 1:(N-1)

    param.dt = tspan(n+1) - tspan(n);    
    
    switch opts.dc
      case{1} % classical deferred correction
        sol(:,n+1) = cdc_step(ode,tspan(n),sol(:,n),opts,param);            
      
      case{2}
        % idc
        sol(:,n+1) = idc_step(ode,tspan(n),sol(:,n),opts,param);
     
        
      case{3}
        % sdc (minion style)
        sol(:,n+1) = sdc_step(ode,tspan(n),sol(:,n),opts,param);
     
      case{4}
        % zadunisky
        sol(:,n+1) = z1976_step(ode,tspan(n),sol(:,n),opts,param);
      
      otherwise
        error('opts.dc not valid');
    end
    

end 







