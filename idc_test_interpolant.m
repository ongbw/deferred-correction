% use idc to correct.
% note: very inefficient wrt memory

clear;

WRITETOFILE = 0;

tspan = [0 17.065];
y0 = [0.994; 
      0; 
      0;
      -2.0015851063790825];
ndim = length(y0);
ode = @(t,y) orbit(t,y);

% use ode45 to get adaptive time step
%[t,~] = ode45(ode,tspan,y0);
t = linspace(tspan(1),tspan(2),1e5);

% recompute provisional solution using
% forward euler predictor
dt = t(2:end) - t(1:end-1);
Nt = length(t);

y = zeros(Nt,ndim);
f = zeros(Nt,ndim);

y(1,:) = y0';
f(1,:) = ode(t(1),y(1,:));

for nt = 2:Nt
    y(nt,:) = y(nt-1,:) + dt(nt-1)*f(nt-1,:);
    f(nt,:) = ode(t(nt),y(nt,:));
end


figure(1);
plot(y(:,1),y(:,2))
title('predictor');
xlabel('y1')
ylabel('y2')
drawnow;


if WRITETOFILE
    % save predictor to file
    dlmwrite('adapt_predictor.dat',[y(1:1e3:end,1:2);y(end,1:2)]);
end

% number of correction levels
ncorrections = 2;


for ncorr = 1:ncorrections
    
    % construct spline for each dimension
    for n = 1:ndim
        isp{n} = splineint(t,f(:,n));
    end
    
    % initialize error vector
    %err  = zeros(Nt,ndim);
    ynew = zeros(Nt,ndim);
    fnew = zeros(Nt,ndim);
    ynew(1,:) = y0';
    fnew(1,:) = f(1,:);
    
    % forward euler correctors    
    for nt = 2:Nt
        %tk = t(nt-1);
        %dt = dt(nt-1);
        
        for n = 1:ndim
            %dt(nt-1)
            isp{n}.coefs(nt-1,:);
            %            pause
            %ppval(isp{n}.coefs(nt-1,:),t(nt)) 
            
            ynew(nt,n) = ynew(nt-1,n)  + ...
                dt(nt-1)*(fnew(nt-1,n) - f(nt-1,n)) + ...
                (polyval(isp{n}.coefs(nt-1,:),t(nt)) - polyval(isp{n}.coefs(nt-1,:),t(nt-1)));
            %err(nt,n) = err(nt-1,n) - (y(nt,n) - y(nt-1,n)) + ...
            %dt(nt-1)*(fnew(nt-1,n) - f(nt-1,n)) + ...
            %    (polyval(isp{n}.coefs(nt-1,:),t(nt)) - polyval(isp{n}.coefs(nt-1,:),t(nt-1)));
        end
        
        % temp = err(nt-1,:) + y(nt-1,:);

        % K = dt*ode(tk,temp);
        % for n = 1:ndim
        %     err(nt,n) = err(nt-1,n) + K(n) -  ...    
        %         dt*ppval(dsp{n},tk);
        % end
        %ynew(nt,:) = y(nt,:) + err(nt,:);
        fnew(nt,:) = ode(t(nt),ynew(nt,:));
    end
    
    %ynew = y + err;

    y = ynew;
    f = fnew;
    
    figure(ncorr+1);
    plot(y(:,1),y(:,2))
    title(sprintf('%d correction(s)',ncorr));
    xlabel('y1')
    ylabel('y2')
    drawnow;
    
    if WRITETOFILE
        % save correctors to file
        fid = sprintf('adapt_corr%d.dat',ncorr);
        dlmwrite(fid,[y(1:1e3:end,1:2);y(end,1:2)]);
    end

end
