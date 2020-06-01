function sol = idc_step(ode,ti,yi,opts,param)

% get "unscaled" quadrature nodes
switch opts.grid
  case 1 
    % uniform nodes
    t = linspace(0,1,opts.nquad);

  case 2
    % Gaussian (Gauss Legendre) nodes
    t = gaussnodes(0,1,opts.nquad);
    
  case 3
    % Chebychev 
    t = chebyshevnodes(0,1,opts.nquad);
        
  case 4
    % Gauss-Lobatto nodes
    t = gauss_lobatto(0,1,opts.nquad);
    
  otherwise
    error('unknown option in opts.grid')
    
end % switch opts.grid

m = length(yi);

% storing solution at quadrature nodes
m = length(yi);
y{1} = zeros(m,opts.nquad);

h = (t(2:end)-t(1:end-1))*param.dt;


        
% predictor
level = 1;
switch param.pred.type
  case{0}
    % explicit
    
    % treat first substep specially since we don't want to pollute
    % structure for values at quadrature nodes
    K = zeros(m,param.pred.s);
    tk = ti;
    
    switch opts.grid
      case {1,4} 
        % uniform node/gauss lobattos (includes left endpoint)
        y{1}(:,1) = yi;
        
      otherwise
        h0 = param.dt*(t(1)-0);
        for s = 1:param.pred.s
            temp = yi + K(:,1:s-1)*param.pred.A(s,1:s-1)';
            K(:,s) = h0*ode(tk+param.pred.c(s)*h0,temp);
        end
        y{1}(:,1) = yi + K*param.pred.b';
        tk = tk+h0;
    end
    
    for n = 2:opts.nquad
        K = zeros(m,param.pred.s);
        for s = 1:param.pred.s
            temp = y{1}(:,n-1) + K(:,1:s-1)*param.pred.A(s,1:s-1)';
            K(:,s) = h(n-1)*ode(tk+param.pred.c(s)*h(n-1),temp);
        end
        
        y{1}(:,n) = y{1}(:,n-1) + K*param.pred.b';
        tk = tk + h(n-1);
    end
    
  case{1}
    % diagonally implicit
    
    % structure for values at quadrature nodes
    K = zeros(m,param.pred.s);
    tk = ti;
    
    switch opts.grid
      case {1,4} 
        % uniform node/gauss lobattos (includes left endpoint)
        y{1}(:,1) = yi;
        
      otherwise
        h0 = param.dt*(t(1)-0);
        for s = 1:param.pred.s
            temp = yi + K(:,1:s-1)*param.pred.A(s,1:s-1)';            
            func_hand = @(Ki) Ki - h0*ode(tk+param.pred.c(s)*h0,...
                                          temp+param.pred.A(s,s)*Ki);
            k_init = param.dt*ode(tk+param.pred.c(s)*h0,...
                                  yi+param.pred.c(s)*param.dt*ode(tk,yi));
            [K(:,s),~,flag] = ...
                fsolve(func_hand, k_init, optimset('Display','off',...
                                                   'TolFun',1e-9));
        
            if flag > 1
                error('did not converge to root')
            end
            
        end
        y{1}(:,1) = yi + K*param.pred.b';
        tk = tk+h0;
    end
    
    for n = 2:opts.nquad
        K = zeros(m,param.pred.s);
        for s = 1:param.pred.s
            temp = y{level}(:,n-1) + K(:,1:s-1)*param.pred.A(s,1:s-1)';
            func_hand = @(Ki) Ki - h(n-1)*ode(tk+param.pred.c(s)*h(n-1),...
                                              temp+param.pred.A(s,s)*Ki);
            k_init = param.dt*ode(tk+param.pred.c(s)*h(n-1),...
                                  y{level}(:,n-1)+...
                                  param.pred.c(s)*param.dt*ode(tk,y{level}(:,n-1)));
            [K(:,s),~,flag] = ...
                fsolve(func_hand, k_init, optimset('Display','off',...
                                                   'TolFun',1e-9));

            if flag > 1
                error('did not converge to root')
            end
        end
        
        y{1}(:,n) = y{1}(:,n-1) + K*param.pred.b';
        tk = tk + h(n-1);
    end
  otherwise
    error('fully implicit integrators not coded yet')
    
end

        
if opts.levels > 1
    % there are correction loops!
    % in practice, can compute S{} and L{} once and reuse for all
    % large intervals
    
    % note, we scale by dt later -- this is the key difference from
    % former code
    S = generate_integration_matrix(t,t);
    
    
    % n = 1
    Sstage{1} = generate_integration_matrix(t,param.corr.c*t(1));
    Lstage{1} = generate_interpolation_matrix(t,param.corr.c*t(1));

    for n = 2:opts.nquad
        Sstage{n} = generate_integration_matrix(t,t(n-1)+param.corr.c*(t(n)-t(n-1)));
        Lstage{n} = generate_interpolation_matrix(t,t(n-1)+param.corr.c*(t(n)-t(n-1)));
    end
end

        
for level = 2:opts.levels

    y{level} = zeros(m,opts.nquad);
    % correction loops -- can generalize to different tableaus later

    % compute f(t,y) at quadrature points
    f = zeros(m,opts.nquad);
    for n = 1:opts.nquad
        f(:,n) = ode(ti+param.dt*t(n),y{level-1}(:,n));
    end

    
    switch param.corr.type 
      case{0}
        % explicit
            
        K = zeros(m,param.corr.s);
        tk = ti;
        
        switch opts.grid
          case {1,4}
            % uniform nodes/gauss lobatto (includes left endpoint)
            y{level}(:,1) = yi;
        
          otherwise
            % again, handle first time step differently
            for s = 1:param.corr.s
                temp = yi + K(:,1:s-1)*param.corr.A(s,1:s-1)' + ...
                       param.dt * f * Sstage{1}(s,:)';
                K(:,s) = h0*ode(tk+param.corr.c(s)*h0,temp) - ...
                         h0*f*Lstage{1}(s,:)'; 
            end
            y{level}(:,1) = yi + K*param.corr.b' + param.dt*f*S(1,:)';
            
            tk = tk+h0;
        end
        for n = 2:opts.nquad
            K = zeros(m,param.corr.s);
            for s = 1:param.corr.s
                temp = y{level}(:,n-1) + K(:,1:s-1)*param.corr.A(s,1:s-1)' + ...
                   param.dt * f * (Sstage{n}(s,:)' - S(n-1,:)');

                K(:,s) = h(n-1)*ode(tk+param.corr.c(s)*h(n-1),temp) - ...
                         h(n-1)*f*Lstage{n}(s,:)';            

            end

            y{level}(:,n) = y{level}(:,n-1) + K*param.corr.b' + param.dt*f*(S(n,:)'-S(n-1,:)');
            tk = tk + h(n-1);
        end        
            
      case{1}
        % diagonally implicit

        K = zeros(m,param.corr.s);
        tk = ti;
        
        switch opts.grid
          case {1,4}
            % uniform nodes/gauss lobatto (includes left endpoint)
            y{level}(:,1) = yi;
        
          otherwise
            % again, handle first time step differently
            for s = 1:param.corr.s
                temp = yi + K(:,1:s-1)*param.corr.A(s,1:s-1)' + ...
                       param.dt * f * Sstage{1}(s,:)';

                func_hand = @(Ki) Ki - ...
                    h0*ode(tk+param.corr.c(s)*h0,...
                           temp+param.corr.A(s,s)*Ki) + ...
                    h0*f*Lstage{1}(s,:)'; 
                
                k_init = h0*ode(tk+param.corr.c(s)*h0,...
                                yi+param.corr.c(s)*param.dt*ode(tk,yi));

                [K(:,s),~,flag] = ...
                    fsolve(func_hand, k_init, optimset('Display','off',...
                                                       'TolFun',1e-9));
                if flag > 1
                    error('did not converge to root')
                end
                
            end
            y{level}(:,1) = yi + K*param.corr.b' + param.dt*f*S(1,:)';
            
            tk = tk+h0;
        end % switch
        
        for n = 2:opts.nquad
            K = zeros(m,param.corr.s);
            for s = 1:param.corr.s
                temp = y{level}(:,n-1) + K(:,1:s-1)*param.corr.A(s,1:s-1)' + ...
                   param.dt * f * (Sstage{n}(s,:)' - S(n-1,:)');


                func_hand = @(Ki) Ki - ...
                    h(n-1)*ode(tk+param.corr.c(s)*h(n-1),...
                               temp+param.corr.A(s,s)*Ki) + ...
                    h(n-1)*f*Lstage{n}(s,:)'; 
                
                k_init = h(n-1)*ode(tk+param.corr.c(s)*h(n-1),...
                                    y{level}(:,n-1)+param.corr.c(s)*param.dt*ode(tk,yi));

                [K(:,s),~,flag] = ...
                    fsolve(func_hand, k_init, optimset('Display','off',...
                                                       'TolFun',1e-9));                
                
                if flag > 1
                    error('did not converge to root')
                end
            end

            y{level}(:,n) = y{level}(:,n-1) + K*param.corr.b' + param.dt*f*(S(n,:)'-S(n-1,:)');
            tk = tk + h(n-1);
        end            
        
      otherwise
        % implicit
        error('fully implicit integrators not coded yet');
    end
end

if ~isfield(opts,'interp')
    % if opts.interp is not set, default to collocation solution
    opts.interp = 0;
end


        

switch opts.grid
  case{1,4} % uniform or gauss--lobatto
    sol = y{opts.levels}(:,end);

  otherwise % chebychev or gauss--legendre
    if opts.interp == 1
      % use interpolation!
      Lfinal = generate_interpolation_matrix(t,1); 
      sol = y{opts.levels}*Lfinal';
      
    else 
      % generate picard/collocation solution
      S_final = generate_integration_matrix(t,1);
    
      sol = yi;
      for n = 1:opts.nquad
          sol = sol +  param.dt*ode(ti + t(n)*param.dt,y{opts.levels}(:,n))*S_final(n);
      end
    end
    
end

