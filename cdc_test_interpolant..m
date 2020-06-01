% use cdc to correct using an underlying interpolant.
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

% forward euler predictor
Nt = 1e5;
t = linspace(tspan(1),tspan(2),Nt);
dt = (tspan(2)-tspan(1))/(Nt-1);

y = zeros(Nt,ndim);
y(1,:) = y0';

for nt = 2:Nt
    y(nt,:) = y(nt-1,:) + dt*ode(t,y(nt-1,:));
end


figure(1);
plot(y(:,1),y(:,2))
title('predictor');
xlabel('y1')
ylabel('y2')
drawnow;

if WRITETOFILE
    % save predictor to file
    dlmwrite('predictor.dat',[y(1:1e3:end,1:2);y(end,1:2)]);
end

% number of correction levels
ncorrections = 1;


for ncorr = 1:ncorrections
    
    % construct spline for each dimension
    for n = 1:ndim
        dsp{n} = splinederiv(t,y(:,n));
    end
    
    % initialize error veector
    err = zeros(Nt,ndim);
    
    % forward euler correctors    
    for nt = 2:Nt
        tk = t(nt-1);
        dt = t(nt) - t(nt-1);
        temp = err(nt-1,:) + y(nt-1,:);

        K = dt*ode(tk,temp);
        for n = 1:ndim
            err(nt,n) = err(nt-1,n) + K(n) -  ...    
                dt*ppval(dsp{n},tk);
        end
    end
    
    ynew = y + err;

    
    y = ynew;
    
    figure(ncorr+1);
    plot(y(:,1),y(:,2))
    title(sprintf('%d correction(s)',ncorr));
    xlabel('y1')
    ylabel('y2')
    drawnow;
    
    if WRITETOFILE
        % save correctors to file
        fid = sprintf('corr%d.dat',ncorr);
        dlmwrite(fid,[y(1:1e3:end,1:2);y(end,1:2)]);
    end

end
