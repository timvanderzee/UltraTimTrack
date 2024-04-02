

figure(10)

geofeatures = handles.geofeatures;
norm_fun = @(c, x) c(1)*exp(-(x-c(2)).^2/(2*c(3)^2)) + c(4);

for i = 1:handles.NumFrames
   
    
%     bar(geofeatures(i).alpha.x, geofeatures(i).alpha.y);
%     hold on
% 
%     plot(geofeatures(i).alpha.x, norm_fun([geofeatures(i).alpha.A geofeatures(i).alpha.mu geofeatures(i).alpha.sigma geofeatures(i).alpha.b], geofeatures(i).alpha.x),'linewidth',2)
%     hold off
%     
%     ylim([0 500])
%     
%     
%     drawnow
%     pause
%     
%     
    mu(i) = geofeatures(i).alpha.mu;
    mu2(i) = geofeatures(i).alpha.median;
    sigma(i) = geofeatures(i).alpha.sigma;
    
%     sigma2(i) = std(geofeatures(i).alpha
    
end

figure(11)
subplot(121)
plot(mu); hold on
plot(mu2)

subplot(122)
plot(sigma)
%% refit
cost_fun = @(c, x, y) sum((y - norm_fun(c,x)).^2);
norm_fun = @(c, x) c(1)*exp(-(x-abs(c(2))).^2/(2*c(3)^2)) + c(4);

for i = 2:handles.NumFrames
   
    disp(i)
    a = geofeatures(i).alpha.x;
    w = geofeatures(i).alpha.y;

    C0 = [max(w)-min(w) 20 10 min(w)];
    C = fminsearch(@(p) cost_fun(p, a, w), C0);

    geofeatures(i).alpha.A = C(1);
    geofeatures(i).alpha.mu = C(2);
    geofeatures(i).alpha.sigma = abs(C(3));
    geofeatures(i).alpha.b = C(4);
    geofeatures(i).alpha.y = w;
    geofeatures(i).alpha.x = a;

    sigma(i) = geofeatures(i).alpha.sigma;
    mu(i) = geofeatures(i).alpha.mu;
end

%%
X = a;
COUNTS = w;

mS = sum(COUNTS.*X)/sum(COUNTS); % Identical to how you calculate the centre of mass.
stdS = sqrt(sum(COUNTS.*(X-mS).^2)/sum(COUNTS));




%%
x = mu(:);
x(:,2) = handles.Region.fas_pen(:) * 180/pi;

figure(10)
plot(diff(x(:,1))); hold on
plot(diff(x(:,2))); hold on
%%
fs = handles.FrameRate;
if ishandle(11),close(11);end
for i = 1
y = fft(x(:,i)-mean(x(:,i)));
f = (0:length(y)-1)*fs/length(y);

figure(11)
plot(f,abs(y)); hold on
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Magnitude')
end