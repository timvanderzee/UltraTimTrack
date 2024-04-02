

figure(10)

geofeatures = handles.geofeatures;
norm_fun = @(c, x) c(1)*exp(-(x-c(2)).^2/(2*c(3)^2));

for i = 1:handles.NumFrames
   
    
    bar(geofeatures(i).alpha.x, geofeatures(i).alpha.y);
    hold on

    plot(geofeatures(i).alpha.x, norm_fun([geofeatures(i).alpha.A geofeatures(i).alpha.mu geofeatures(i).alpha.sigma], geofeatures(i).alpha.x),'linewidth',2)
    hold off
    
    ylim([0 500])
    
    
    drawnow
    pause
    mu(i) = geofeatures(i).alpha.mu;
    sigma(i) = geofeatures(i).alpha.sigma;
    
%     sigma2(i) = std(geofeatures(i).alpha
    
end

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