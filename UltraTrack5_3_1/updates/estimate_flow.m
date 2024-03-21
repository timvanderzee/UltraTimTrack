sXY = 1; sS = 1;
af = affine_flow('image1', im1, 'image2', im2, ...
    'sigmaXY', sXY, 'sampleStep', sS);
af.sampleMethod = 'logpolar';
af.logWedges = 200;
af.sigmaRW = 0.1;
af.logRmin = 5;
af.logRmax = 100;
af.maxIter = 5;  % allow up to 5 iterations
af.velTol = 0.5; % but stop if centre flow speed drops below 1 pixel/frame
af = af.findFlow;
af.image2 = im2;    % second real image
af = af.findFlow;
ffound = af.flowStruct;
sigma = round(abs(mean([ffound.vx0;ffound.vy0])),0)