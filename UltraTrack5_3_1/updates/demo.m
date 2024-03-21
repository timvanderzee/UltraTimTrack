%% Affine Optic Flow Demonstration from David Young
%https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/27093/versions/6/previews/html/affine_flowdemo.html#8

%Demonstration of estimation of an affine (first-order) optic flow field 
%from a pair of images, using affine_flow.
%affine_flow does its own resampling, and you always give it the original 
%images to work on.

%Estimating the global flow field across the whole image is shown here, but 
%it can be also be used to make local estimates of the flow field. This is 
%done with the rectRegion or regionOfInterest options for conventional 
%sampling, or by setting logRmax and logCentre for log-polar sampling. 

%% Setting up images
%First read in two test images. These should be in the current directory.
%If features lie very roughly in a plane, the first-order flow provides an approximation to the true flow,
%despite not capturing its details or exact form.

% Convert to double for later processing. It's convenient to scale into the
% range 0 to 1 for imshow

% im1 = double(imread('maze1.png'))/256;
% im2 = double(imread('maze2.png'))/256;
load('C:\Users\brent\Desktop\DFG\Experiment\Mech\Pilot\04a\US\BH_A00W_MG_03.mat');
im1 = double(TVDdata.Im(:,:,195))/255;
im2 = double(TVDdata.Im(:,:,196))/255;
figure; imshow(im1);
figure; imshow(im2);

%% Measuring flow with conventional sampling
%Start by using the default rectilinear sampling.
%The motion between the test images is quite large, 
%so to measure it we need to smooth quite heavily. 
%We try a sigma of 25 pixels, and also sample every 25 pixels to cut down computation.

sXY = 15; sS = 15;
af = affine_flow('image1', im1, 'image2', im2, ...
    'sigmaXY', sXY, 'sampleStep', sS);
af = af.findFlow;
flow = af.flowStruct;

%% Display the estimated flow graphically.

%A set of flow vectors illustrating the estimated flow field is displayed 
%on the first image. This does not show the points at which the flow was 
%computed; these are just representative vectors at arbitrary points to 
%show the form of the first-order model.

%The flow vectors near the bottom of the image are larger than those higher 
%up. This is because the surface is closer to the camera in this region; 
%the flow vectors are like stereo disparities, larger (for parallel 
%cameras) for closer objects.

affine_flowdisplay(flow, im1, 50);

%Check whether the estimated flow registers the images.

%We can see whether the flow that has been found maps the first image onto 
%the second, by displaying edges. Here, the green edges are the original 
%first image, the blue edges are the second image, and the red edges are 
%the first image after warping by the flow.

%If the process has worked correctly, red and blue edges should be close 
%together. We do not expect exact overlap because the first-order flow can 
%only be an approximate, smooth, model of the true flow. Even the overall 
%form of the real flow is a perspective rather than an affine 
%transformation, and there is a lot of depth variation which adds 
%complexity. Nonetheless, the affine flow gives a respectable first 
%approximation.

affine_flowedgedisplay(flow, im1, im2);

%% Assessing accuracy using synthetic flow
%One way to assess the accuracy of the flow estimation is to synthesise the 
%second image by warping the first with a known flow field, then seeing if 
%we can recover the parameters of the field.

%First, we set some parameters, look at the flow field they generate, and 
%warp the first image according to this flow.

% The parameters for the test flow field
ftest.vx0 = 5;     % flow at centre of image (values changed later when
ftest.vy0 = 5;     % origin is moved)
ftest.d = 0.05;
ftest.r = -0.05;
ftest.s1 = 0.05;
ftest.s2 = -0.05;
% Shift origin to image origin
ftest = affine_flow.shift(ftest, -(size(im1,2)+1)/2, -(size(im1,1)+1)/2);

% Warp the first image
wtest = affine_flow.warp(ftest);    % convert to a warping matrix
ttest = maketform('affine', wtest); % make a transform structure
imtrans = imtransform_same(im1, ttest);

% Display the transformed image
imshow(imtrans);
% The flow is estimated as for the real image pair

%In general, sigmaXY (the smoothing constant) needs to be set either using 
%some prior knowledge of the likely maximum flow speed, or by experiment. 
%The spatial scale of the smoothing needs to be comparable to the maximum 
%optic flow speed.

%In this case it is easy to estimate the order of magnitude of the flow 
%from the known parameters and the image size (use the equation in 
%affine_flow's help) and a value of 10 for the smoothing constant seems 
%reasonable. It would be possible to refine this experimentally, given an 
%error measure, a set of images, and a set of flow fields.

%Since we know the images will be smoothed on a scale of 10 pixels, we can 
%cut down the computation by sampling the gradients every tenth row and 
%column before solving the least-squares problem. This is done with the 
%sampleStep parameter.

%We can see that the estimates are in the right region; whether they are 
%adequate for any task depends, of course, on the application.

af.sigmaXY = 10;
af.sampleStep = 10;
af.image2 = imtrans;
af = af.findFlow;
ffound = af.flowStruct;

% Inspect the values
disp('The test flow parameters:'); disp(ftest);
disp('The recovered flow parameters:'); disp(ffound);

% Show the input flow field in blue and the recovered flow field in red,
% with slightly different positioning so that both sets of vectors can be
% seen.
clf;
affine_flowdisplay(ftest, size(im1), 50, 'b');
affine_flowdisplay(ffound, size(im1), 45, 'r');

%% Finding flow using log-polar sampling
%If required, affine_flow resamples both images to log-polar coordinates, 
%and takes advantage of the change in resolution to estimate first-order 
%flow. To do this, it shifts the log-polar centre for the second image 
%using its initial estimate of the flow, so as to track the motion. This 
%produces low flow speeds at the high-resolution centre of the sampling 
%pattern, with the flow speed increasing (for pure first-order flow) 
%linearly with radius, matching the sample separation.

rmin = 3;                   % See help logsample for what these
rmax = min(size(im1))/2;    % parameters mean
xc = size(im1,2)/2;
yc = size(im1,1)/2;
nw = 300;
logim = logsample(im1, rmin, rmax, xc, yc, [], nw);
imshow(logim);

%As with the conventional approach, we need to set a smoothing constant. 
%We will also set some other sampling parameters - see the help information 
%for logsample.

%Smoothing is done on the resampled log-polar image, so the smoothing 
%constant can be smaller than in the conventional case - if the translation 
%is not too great, the expected norm of the first-order components divided 
%by the angle between wedges is a reasonable order of magnitude. (The angle 
%between wedges is 2*pi/logWedges.) Here, we adopt 200 wedges and a 
%smoothing constant of 4, but experimentation on a range of images and flow 
%fields could be done.

%We set logRmin to 5 - it isn't critical, but small rings suffer badly from 
%the mismatch between the original image resolution and the higher 
%resolution that the foveal centre of the log-polar image ought to have, 
%and there is little point in making logRmin too small.

%We set logRmax to 100 pixels. Setting it explicitly allows us to compare 
%the results without and with tracking.

%With no tracking (i.e. af.maxIter = 1), the function makes a single 
%estimate of the flow field, with both log-polar patterns centred at the 
%centres of the images (the default).

%We see that the estimate is reasonable, but not as accurate as the 
%conventional computation above. The lastSpeed field shows that the image 
%was estimated to be moving past the centre of the sampling pattern at 
%x pixels/frame - enough to destroy any contribution from the 
%high-resolution inner rings.

af.sampleMethod = 'logpolar';
af.logWedges = 200;
af.sigmaRW = 4;
af.logRmin = 5;
af.logRmax = 100;
af.maxIter = 5;  % allow up to 5 iterations
af.velTol = 0.5; % but stop if centre flow speed drops below 1 pixel/frame

af = af.findFlow;
ffound = af.flowStruct;

% Compare true flow and result
disp('Test flow parameters:'); disp(ftest);
disp('Recovered flow, log-polar sampling, tracking:'); disp(ffound);
disp('Tracking information:'); disp(af.trackInfo);

% and a visual check
affine_flowedgedisplay(ffound, im1, imtrans);

%We now perform the same operation but using the real image pair. Only a 
%qualitative comparison with the conventional approach is possible, but the 
%log-polar method appears to perform somewhat better in this case.
af.image2 = im2;    % second real image

af = af.findFlow;
ffound = af.flowStruct;

disp('Recovered flow for real images:'); disp(ffound);
disp('Tracking information:'); disp(af.trackInfo);

% and a visual check
affine_flowedgedisplay(ffound, im1, im2);