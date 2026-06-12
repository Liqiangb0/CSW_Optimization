%% Initialize trapezoidal corrugated steel web as initial configuration
% Plot original grayscale image and downsampled image
A = imread('BCSW1600_430x370x220x36mm.jpg');
figure();tiledlayout(2,1);
nexttile(2);
imshow(A);
nexttile(1);
B = 1 - double(imresize(rgb2gray(A),0.06))/255;
imagesc(1 - B);
colormap('gray');axis equal tight off;
%% Initialization
[nelx,nely] = size(B);nelz = 24;
Volfrac = mean(B(:));
penalK = 3;
rmin = 2;
filter = 2; filterBC = 'N';
eta = 0.5; beta = 2;
ocParam = [0.1,0.7,1.2];
itmax = 800;
Lx = 300;
penalG = 3;
nEig = 12;
pAgg = 50;
prSel = {['B','C','V'],[5,Volfrac]};
xPhysMat = repmat(B,1,1,nelz);
save('initialData.mat','xPhysMat');
CSWBuck(nelx,nely,nelz,...
    penalK,rmin,filter,filterBC,...
    eta,beta,ocParam,itmax,...
    Lx,penalG,nEig,pAgg,prSel,'initialData.mat');