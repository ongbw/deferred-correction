function [sol,t] = pipeline_deferred_correction(ode,tspan,y0,opts)


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

h = tspan(2:end)-tspan(1:end-1);

% predictor
level = 1;
y{1}(:,1) = y0;

switch param.pred.type
  case{0}
    % explicit

    for n = 2:N
        % stages for RK method
        K = zeros(m,param.pred.s);
        for s = 1:param.pred.s
            temp = y{1}(:,n-1) + K(:,1:s-1)*param.pred.A(s,1:s-1)';
            K(:,s) = h(n-1)*ode(tspan(n-1)+param.pred.c(s)*h(n-1),temp);
        end
        
        y{1}(:,n) = y{1}(:,n-1) + K*param.pred.b';
    end
    
  case{1}
    % diagonally implicit
    
    for n = 2:N
        % stages for RK method
        K = zeros(m,param.pred.s);
        tk = tspan(n-1);
        for s = 1:param.pred.s
            temp = y{level}(:,n-1) + K(:,1:s-1)*param.pred.A(s,1:s-1)';
            func_hand = @(Ki) Ki - h(n-1)*ode(tk+param.pred.c(s)*h(n-1),...
                                              temp+param.pred.A(s,s)*Ki);
            k_init = h(n-1)*ode(tk+param.pred.c(s)*h(n-1),...
                                  y{level}(:,n-1)+...
                                  param.pred.c(s)*h(n-1)*ode(tk,y{level}(:,n-1)));
            [K(:,s),blah,flag] = ...
                fsolve(func_hand, k_init, optimset('Display','off',...
                                                   'TolFun',1e-9));

            if flag > 1
                error('did not converge to root')
            end
        end
        
        y{1}(:,n) = y{1}(:,n-1) + K*param.pred.b';
        tk = tk + h(n-1);
    end
    
  case{2}
    % built in integrator
    
    % use ode45 to get provisional solution
    sol = ode45(ode,[tspan(1),tspan(end)],y0);
    tspan = sol.x;
    y{1} = sol.y;
    
    % redefine variables based on tspan.
    N = length(tspan);
    m = length(y0);
    h = tspan(2:end)-tspan(1:end-1);
    
  otherwise
    error('fully implicit integrators not coded yet')
    
end

if opts.levels > 1
    
    switch opts.dc
      case{1} % classical deferred correction
        for level = 2:opts.levels
            e{level} = zeros(m,N);

            % compute number of interpolation points to use.
            ngrid = param.pred.order + (level-1)*param.corr.order + 1;

            if length(tspan) < ngrid
                error(['not enough interpolation points to obtain desired ' ...
                       'accuracy']);
            end

            for n = 2:N
                
                ind_start = min(n-1,N-ngrid + 1);
                ind_end =  min(n - 1 + ngrid - 1,N);

                tk = tspan(n-1);
                grid = tspan(ind_start:ind_end);
                lgrid = grid(end)-grid(1);                
                scaled_grid = (grid - grid(1)) / lgrid;

                
                % RK correctors 
                staget = tspan(n-1) + param.corr.c*h(n-1);
                staget = (staget-grid(1))/lgrid;


                Dstage = generate_differentiation_matrix(scaled_grid,staget);
                Dstage = Dstage / lgrid;
                Lstage = generate_interpolation_matrix(scaled_grid,staget);

                
                switch param.corr.type
                  case{0} % explicit
                    K = zeros(m,param.corr.s);
                    for s = 1:param.corr.s
                        temp = e{level}(:,n-1) + K(:,1:s-1)*param.corr.A(s,1:s-1)' + ...
                               y{level-1}(:,ind_start:ind_end) * Lstage(s,:)';
                        
                        %Lstage(s,:)
                        %pause
                        
                        K(:,s) = h(n-1)*ode(tk+param.corr.c(s)*h(n-1),temp) - ...
                                 h(n-1)*y{level-1}(:,ind_start:ind_end)*Dstage(s,:)';            
                        
                    end

                    
                  case{1} % diagonally implicit
                    
                    K = zeros(m,param.corr.s);
                    for s = 1:param.corr.s
                        temp = e{level}(:,n-1) + ...
                               K(:,1:s-1)*param.corr.A(s,1:s-1)' + ...
                               y{level-1}(:,ind_start:ind_end)*Lstage(s,:)';
                
                        func_hand = @(Ki) Ki - ...
                            h(n-1)*ode(tk+param.corr.c(s)*h(n-1),...
                                       temp+param.corr.A(s,s)*Ki) + ...
                            h(n-1)*y{level-1}(:,ind_start:ind_end)*Dstage(s,:)'; 
                
                        k_init = h(n-1)*ode(tk+param.corr.c(s)*h(n-1),temp) - ...
                                 y{level-1}(:,ind_start:ind_end)*Dstage(s,:)'; 

                        [K(:,s),blah,flag] = ...
                            fsolve(func_hand, k_init, optimset('Display','off',...
                                                               'TolFun',1e-12));
                        if flag > 1
                            error('did not converge to root')
                        end
                    end
                                    
                  otherwise
                    error('param.corr.type not valid');
                    end % switch
                    
                e{level}(:,n) = e{level}(:,n-1) + K*param.corr.b' ;
                                
            end
            y{level} = y{level-1} + e{level};
        end
        
        
      case{2} % sdc and idc
        
        for level = 2:opts.levels
            % compute f(t,y) at quadrature points
            dydt = zeros(m,N);
            for n = 1:N
                dydt(:,n) = ode(tspan(n),y{level-1}(:,n));
            end

            y{level} = zeros(m,N);
            y{level}(:,1) = y0;            

            % compute number of quadrature points to use.
            nquad = param.pred.order + (level-1)*param.corr.order;

            
            if length(tspan) < nquad
                error(['not enough quadrature points to obtain desired ' ...
                       'accuracy']);
            end

            for n = 2:N
                ind_start = min(n-1,N-nquad + 1);
                ind_end =  min(n - 1 + nquad - 1,N);

                tk = tspan(n-1);
                f = dydt(:,ind_start:ind_end);
                t = tspan(ind_start:ind_end);
                
                quad_interval = tspan(n-1:n);
                dt = t(end)-t(1);
                tt = (t - t(1))/dt;
                quad_interval = (quad_interval-t(1))/dt;
                S = generate_integration_matrix(tt,quad_interval);
                S = dt*S;
                

                % RK correctors 
                quad_interval = tspan(n-1) + param.corr.c*h(n-1);
                quad_interval = (quad_interval-t(1))/dt;

                Sstage = generate_integration_matrix(tt,quad_interval);
                Sstage = dt*Sstage;
                Lstage = generate_interpolation_matrix(tt,quad_interval);

                K = zeros(m,param.corr.s);
                switch param.corr.type
                  case{0} % explicit
                
                    for s = 1:param.corr.s
                        temp = y{level}(:,n-1) + K(:,1:s-1)*param.corr.A(s,1:s-1)' + ...
                               f * (Sstage(s,:)' - S(1,:)');
                        
                        K(:,s) = h(n-1)*ode(tk+param.corr.c(s)*h(n-1),temp) - ...
                                 h(n-1)*f*Lstage(s,:)';            
                        
                    end                
                
                    y{level}(:,n) = y{level}(:,n-1) + K*param.corr.b' + f*(S(2,:)'-S(1,:)');
                  
                  case{1} % diagonally implicit

                    for s = 1:param.corr.s
                        temp = y{level}(:,n-1) + K(:,1:s-1)*param.corr.A(s,1:s-1)' + ...
                               f * (Sstage(s,:)' - S(1,:)');


                        func_hand = @(Ki) Ki - ...
                            h(n-1)*ode(tk+param.corr.c(s)*h(n-1),...
                                       temp+param.corr.A(s,s)*Ki) + ...
                            h(n-1)*f*Lstage(s,:)'; 
                
                        k_init = h(n-1)*ode(tk+param.corr.c(s)*h(n-1),...
                                            y{level}(:,n-1)+param.corr.c(s)*h(n-1)*ode(tk,y{level}(:,n-1)));
                        
                        [K(:,s),blah,flag] = ...
                            fsolve(func_hand, k_init, optimset('Display','off',...
                                                               'TolFun',1e-9));                
                        
                        if flag > 1
                            error('did not converge to root')
                        end
                    end

                    y{level}(:,n) = y{level}(:,n-1) + K*param.corr.b' + f*(S(2,:)'-S(1,:)');
                    
                  otherwise
                    error('fully implicit integrators not coded');
                end % switch coor type
            end
                
        end
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      case{3} % classical deferred correction, interpolated grid
        
        if opts.hardcode == 1
            N = 20000;
        end
        
        % specify new uniform mesh, to match Erin's setup
        tnew = linspace(tspan(1),tspan(end),N);
        dt = (tspan(end)-tspan(1))/(N-1);
        
        warning('resetting opts.level = 2 ...')
        
        sp = spline(tspan,y{1});
        dsp = splinederiv(tspan,y{1});
        P = @(x) ppval(sp, x);
        dP = @(x) ppval(dsp, x);
        
        level = 2;
        e{level} = zeros(m,N);

        
        for n = 2:N
                
            tk = tnew(n-1);
                
            switch param.corr.type
              case{0} % explicit
                K = zeros(m,param.corr.s);
                for s = 1:param.corr.s
                    temp = e{level}(:,n-1) + K(:,1:s-1)*param.corr.A(s,1:s-1)' + ...
                           P(tk + param.corr.c(s)*dt);
                    
                    K(:,s) = dt*ode(tk+param.corr.c(s)*dt,temp) - ...
                             dt*dP(tk + param.corr.c(s)*dt);
                end
                
                
              case{1} % diagonally implicit
                
                K = zeros(m,param.corr.s);
                for s = 1:param.corr.s
                    temp = e{level}(:,n-1) + ...
                           K(:,1:s-1)*param.corr.A(s,1:s-1)' + ...
                           P(tk + param.corr.c(s)*dt);
                    
                    func_hand = @(Ki) Ki - ...
                        dt*ode(tk+param.corr.c(s)*dt, ...
                                   temp+param.corr.A(s,s)*Ki) + ...
                        dt*dP(tk+para.corr.c(s)*dt); 
                    
                    k_init = dt*ode(tk+param.corr.c(s)*dt,temp) - ...
                             dP(tk+param.corr.c(s)*dt); 
                    
                    [K(:,s),blah,flag] = ...
                        fsolve(func_hand, k_init, optimset('Display','off',...
                                                           'TolFun',1e-12));
                    if flag > 1
                        error('did not converge to root')
                    end
                end
                
              otherwise
                error('param.corr.type not valid');
            end % switch
            
            e{level}(:,n) = e{level}(:,n-1) + K*param.corr.b' ;
            
        end
        y{level} = P(tnew) + e{level};
        
        tspan = tnew;
        
      otherwise
        error('opts.dc not valid');
    end % switch opts.dc
end






sol = y{level};

if nargout > 1
    t = tspan;
end