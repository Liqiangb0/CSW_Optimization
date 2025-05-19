function loop = CSWBuck(nelx,nely,nelz,penalK,rmin,filter,filterBC,eta,beta,ocParam,itmax,Lx,penalG,nEig,pAgg,prSel,initialData)%Ver. 20241021
tic;updateRule = 'MMA';phiOutputInv = 50;isPlotBucklingMode = 0;isPlotCurrentDesign = 0;isPlotWaveform = 1;
%% PRE. 0) Create folder for saving results
meshtxt = ['Mesh' num2str(nelx) 'x' num2str(nely) 'x' num2str(nelz)];
dictionaryName = [meshtxt '_' char(datetime('now', 'Format', 'yyyyMMddHHmmss'))];
if ~exist(dictionaryName,'dir');mkdir(dictionaryName);end;addpath(dictionaryName);% Create a folder named dictinary
fileTitle = [dictionaryName '\' meshtxt '_' prSel{1}];
topHistoryFileID = fopen([fileTitle '_TopHistory.txt'],'w');% Create txt file to save top. history
fprintf(topHistoryFileID,'%10s\t%2s\t%2s\t%6s\t%6s\t%3s\t%4s\t%4s\t%6s\t%7s\t%7s\t%10s\t%10s\t%8s\t%9s\r\n',...
    'Iteration','g0','g1','penalK','penalG','eta','beta','pAgg','change','lambda1','lambda2','Compliance','Vol. Frac.','CPU time','Real time');
% create figure
[maxDU,threshold] = deal(Lx/5,0.5);                                        % set maximum display displacement & dispaly elements
f_History = figure('Name',[meshtxt '_' prSel{1} '_Obj&ConstValues']);
if isPlotCurrentDesign == 1,f_CurrentDesign = figure('Name',[meshtxt '_' prSel{1} '_Current Design']);end
if isPlotWaveform == 1,f_Waveform = figure('Name',[meshtxt '_' prSel{1} '_Waveform']); end
if any(prSel{1} == 'B')&&isPlotBucklingMode,f_BucklingMode = figure('Name',[meshtxt '_' prSel{1} '_Buckling Mode']);end
%% --------------------PRE. 1) Material and Continuation Parameters
[E0,Emin,nu] = deal(1,1e-6,0.3);                                           % Young's moduli & Poisson's ratio
penalCntK = {25,3,25,0.25};                                                % continuation schema on K-penal
penalCntG = {25,3,25,0.25};                                                % continuation schema on G-penal
betaCnt = {400,2,25,0.5};                                                  % continuation schema on beta
pAggCnt = {1e5,1e8,25,2e3};                                                % continuation schema on KS aggregation factor
cnt = @(v,vCn,l) v+ (l>=vCn{1}).*(v<vCn{2}).*(mod(l,vCn{3})==0).*vCn{4};   % function applying continuation
if prSel{1}(1) == 'V', volfrac = 1.0; else, volfrac = prSel{2}(end); end   % initialize volume fraction
save([fileTitle '_Input.mat'],'nelx','nely','nelz','penalK','rmin',...
    'filter','filterBC','eta','beta','ocParam','itmax','Lx','penalG',...
    'nEig','pAgg','prSel','updateRule','phiOutputInv','penalCntK',...
    'penalCntG','betaCnt','pAggCnt');               % save input parameters
%% ---------------------------------PRE. 2) Discretization Features
load(initialData);x = xPhysMat(:);                                         % Initialization of CSW
[Ly,Lz] = deal(nely/nelx*Lx,nelz/nelx*Lx);                                 % recovor Ly&Lz from aspect ratio
axislim = [-maxDU,Lx+maxDU;-maxDU,Ly+maxDU;-maxDU,Lz+maxDU;];              % display axis limits 
nEl = nelx * nely * nelz;                                                  % number of elements
nNode = (1+nely)*(1+nelx)*(1+nelz);                                        % number of nodes
nDof = nNode*3;                                                            % total number of DOFs
elNrs = reshape(1:nEl,nelx,nely,nelz);                                     % element numbering
nodeNrs = int32(reshape(1:nNode,nelx+1,nely+1,nelz+1)); clear nNode        % node numbering (int32)
dofNrs = int32(reshape(3*nodeNrs(1:nelx,1:nely,1:nelz)-2,nEl,1));          % first dof of element numbering
cMat = dofNrs + int32(repmat([0 1 2 3 4 5 3*(nelx+1)+[3 4 5 0 1 2]],1,2)+...
    [zeros(1,12) 3*(nelx+1)*(nely+1)*ones(1,12)]);  clear dofNrs           % connectivity matrix
% <<<<<<<<<<<<<<<<<<<<<<<<<<following code will be used for post-processing
elements = int32(reshape(nodeNrs(1:nelx,1:nely,1:nelz),nEl,1))+int32(...
    repmat([0 1 nelx+1+[1 0]],1,2)+[zeros(1,4),(nelx+1)*(nely+1)*ones(1,4)]);% element nodes numbering
[nodeX,nodeY,nodeZ] = ndgrid((0:nelx)*Lx/nelx,(0:nely)*Lx/nelx,(0:nelz)*Lx/nelx);
nodes = [nodeX(:), nodeY(:), nodeZ(:)];clear nodeX nodeY nodeZ;
[centerX,centerY,centerZ] = ndgrid((0.5:nelx)*Lx/nelx,(0.5:nely)*Lx/nelx,(0.5:nelz)*Lx/nelx);
efaces = [1, 2, 3, 4;... % face in the plane XOY
          5, 8, 7, 6;... % face in the plane XOY
          1, 5, 6, 2;... % face in the plane XOZ
          2, 6, 7, 3;... % face in the plane YOZ
          3, 7, 8, 4;... % face in the plane XOZ
          4, 8, 5, 1];   % face in the plane YOZ (Orientation points to the inside of the element)
% face numbering(https://abaqus-docs.mit.edu/2017/English/SIMACAEELMRefMap/simaelm-r-3delem.htm)
save([fileTitle '_MeshInfo'],'elements','nodes','efaces','centerX','centerY','centerZ');
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>Post-processing
% -------------------------------elemental stiffness matrix H8 element #3D#
Ke = 1/(1+nu)/(1-2*nu)/144*(...
    [32;  6;  6; -8; -6; -6;-10; -6; -3;  4;  6;  3;  4;  3;  6;-10; -3; -6; -8; -3; -3; -4;  3;  3;...
         32;  6;  6;  4;  3; -6;-10; -3; -6; -8; -6;  3;  4;  6;  3; -4;  3; -3; -8; -3; -3;-10; -6;...
             32;  6;  3;  4;  3;  3; -4;  3;  6;  4; -6; -6; -8; -6; -3;-10; -3; -3; -8; -3; -6;-10;...
                 32; -6; -6;  4; -6; -3;-10;  6;  3;-10;  3;  6;  4; -3; -6; -4; -3; -3; -8;  3;  3;...
                     32;  6;  6; -8; -6;  6;-10; -3; -3; -4;  3; -3;  4;  6;  3;-10; -6;  3; -8; -3;...
                         32; -3;  6;  4; -3;  3; -4;  6; -3;-10;  6; -6; -8;  3; -6;-10;  3; -3; -8;...
                             32;  6; -6; -8; -6;  6; -8; -3;  3; -4;  3; -3;  4;  3; -6;-10; -3;  6;...
                                 32; -6;  6;  4; -3; -3; -8;  3; -3;-10;  6;  3;  4; -6;  3; -4; -3;...
                                     32; -6; -3;  4;  3;  3; -8;  3;  6;-10;  6;  6; -8;  6;  3;-10;...
                                         32; -6;  6; -4; -3;  3; -8;  3; -3;-10;  3; -6;  4; -3;  6;...
                                             32; -6;  3;-10;  6;  3; -8;  3; -3; -4; -3; -3;  4; -6;...
                                                 32; -3;  6;-10; -3;  3; -8; -6;  3;-10; -6;  6; -8;...
                                                     32;  6; -6; -8; -6;  6;-10; -6;  3;  4;  6; -3;...
                                                         32; -6;  6;  4; -3; -6;-10;  3; -6; -8;  6;...
                                                             32; -6; -3;  4; -3; -3; -4; -3; -6;  4;...
                                                                 32; -6;  6;  4; -6;  3;-10;  6; -3;...
                                                                     32; -6;  6; -8;  6;  6;-10;  3;...
                                                                         32;  3; -6;  4;  3; -3; -4;...
                                                                             32;  6;  6; -8; -6; -6;...
                                                                                 32;  6;  6;  4;  3;...
                                                                                     32;  6;  3;  4;...
                                                                                         32; -6; -6;...
                                                                                             32;  6;...
                                                                                                 32]+nu*...
    [-48;  0;  0;  0; 24; 24; 12;  0; 12;  0;-24;  0;  0;  0;-24; 12; 12;  0; 12;  0;  0; 12;-12;-12;...
         -48;  0;-24;  0;  0;  0; 12; 12; 24;  0; 24;  0;  0;-24;-12; 12;-12;  0; 12;  0; 12; 12;  0;...
             -48;-24;  0;  0;-12;-12; 12;  0;-24;  0; 24; 24;  0;  0; 12; 12;  0;  0; 12; 12;  0; 12;...
                 -48;  0;  0;  0; 24;  0; 12;  0;-12; 12;-12;  0;  0;  0; 24; 12; 12; 12; 12;  0;  0;...
                     -48;  0;-24;  0; 24;  0; 12; 12; 12; 12;-12;  0;  0;-24;-12; 12;  0;  0; 12;  0;...
                         -48;  0;-24;  0; 12;-12; 12;  0; 12; 12;-24; 24;  0;-12;  0; 12;  0;  0; 12;...
                             -48;  0;  0;  0; 24;-24; 12;  0;  0; 12;-12; 12;  0;  0; 24; 12; 12;  0;...
                                 -48;  0;-24;  0;  0;  0; 12;  0; 12; 12;  0;  0;  0; 24;-12; 12; 12;...
                                     -48; 24;  0;  0;  0;  0; 12;-12;  0; 12;-24;-24;  0;  0;-12; 12;...
                                         -48;  0;  0; 12; 12;-12; 12;  0;  0; 12;-12;  0;  0;  0;-24;...
                                             -48;  0;-12; 12;  0;  0; 12;  0; 12; 12; 12;  0;  0; 24;...
                                                 -48; 12;  0; 12;  0;  0; 12;  0;-12; 12; 24;-24;  0;...
                                                     -48;  0;  0;  0; 24;-24; 12;  0;-12;  0;-24;  0;...
                                                         -48;  0;-24;  0;  0;  0; 12;-12; 24;  0;-24;...
                                                             -48; 24;  0;  0; 12; 12; 12;  0; 24;  0;...
                                                                 -48;  0;  0;  0; 24;  0; 12;  0; 12;...
                                                                     -48;  0;-24;  0;-24;  0; 12;-12;...
                                                                         -48;  0; 24;  0;-12; 12; 12;...
                                                                             -48;  0;  0;  0; 24; 24;...
                                                                                 -48;  0;-24;  0;  0;...
                                                                                     -48;-24;  0;  0;...
                                                                                         -48;  0;  0;...
                                                                                             -48;  0;...
                                                                                                 -48]);
% ---------------------------------element stiffness martix H8 element #3D#
Ke0(tril(ones(24))==1) = Ke';
Ke0 = reshape(Ke0,24,24);
Ke0 = Ke0 + Ke0' - diag(diag(Ke0));                                        % recover full element stiffness matrix Ke0
[sI,sII] = deal([]);
for j = 1:24                                  % loop on the DOFs of element
    sI = cat(2,sI,j:24);
    sII = cat(2,sII,repmat(j,1,24-j+1));
end
[iK,jK] = deal(cMat(:,sI)',cMat(:,sII)');clear sI sII
indexK = sort([iK(:),jK(:)],2,'descend');
if any(prSel{1}=='B') % >>>>>>>>>>>>>Perform Only If Buckling Is Active #B#
    Dmat0 = 1/((1+nu)*(1-2*nu))*...
        [1-nu     nu      nu             0           0           0;...
           nu	1-nu      nu             0           0           0;...
           nu     nu    1-nu             0           0        	 0;...
            0      0       0	(1-2*nu)/2           0       	 0;...
            0      0       0             0  (1-2*nu)/2	         0;...
            0      0       0             0           0	(1-2*nu)/2];       % non-dimensional elastic matrix
    xiG = sqrt(1/3)*[-1,1]; wxi = [1,1];
    [etaG,zetaG] = deal(xiG); [weta,wzeta] = deal(wxi);                    % Gauss nodes and weights
    [XIG,ETAG,ZETAG] = ndgrid(xiG,xiG,xiG);
    [WXI,WETA,WZETA] = ndgrid(wxi,wxi,wxi);
    Xie = [-1,-1,-1;...
            1,-1,-1;...
            1, 1,-1;...
           -1, 1,-1;...
           -1,-1, 1;...
            1,-1, 1;...
            1, 1, 1;...
           -1, 1, 1];                                                      % natural coordinates of the element e
    Xe = Xie .* Lx/nelx/2;                                                 % dimensions of the element (using for compute Jacobi matrix)
    lMat=zeros(6,9); lMat(1,1)=1; lMat(2,5)=1; lMat(3,9)=1;lMat(4,[2 4])=1;
    lMat(5,[6 8])=1; lMat(6,[3 7])=1;                                      % placement matrix (using for generating strain-displacement matrix)
    dN = @(xi,eta,zeta) [Xie(:,1)'       .* (1+Xie(:,2)'*eta) .* (1+Xie(:,3)'*zeta);...
                        (1+Xie(:,1)'*xi) .* Xie(:,2)'         .* (1+Xie(:,3)'*zeta);...
                        (1+Xie(:,1)'*xi) .* (1+Xie(:,2)'*eta) .* Xie(:,3)']/8;% shape fucntion logical derivatives (xi-, eta-, zeta-derivatives of shape function)
    B0 = @(gradN) lMat * kron(gradN,eye(3));                               % strain-displacement matrix
    %     [indexUniF,indexUniT] = deal([]);nNodeDOF = 3;nElNode = 8;nElDOF = nNodeDOF*nElNode;
    %     for i = 1:nNodeDOF:nElDOF
    %         indexUniT = cat(2,indexUniT,int32(((i-1)*(nElDOF-(i-2)/2)+1):3:i*(nElDOF-(i-1)/2)));
    %         indexUniF = cat(2,IndexUniF,int32((nElDOF*(i-1)+i):3:i*nElDOF));
    %     end
    indexUniF = [1,  4,  7, 10, 13, 16, 19, 22,...
                76, 79, 82, 85, 88, 91, 94,...
               151,154,157,160,163,166,...
               226,229,232,235,238,...
               301,304,307,310,...
               376,379,382,...
               451,454,...
               526];                                                       % Indexing of unique coefficients of the (full) element stress stiffness matrix
    indexUniT = [1,  4,  7, 10, 13, 16, 19, 22,...
                70, 73, 76, 79, 82, 85, 88,...
               130,133,136,139,142,145,...
               181,184,187,190,193,...
               223,226,229,232,...
               256,259,262,...
               280,283,...
               295];                                                       % Indexing of unique coefficients of the (lower tril) element stress stiffness matrix
    t2index = [2,  3,  4,  5,  6,  7,  8,...
              10, 11, 12, 13, 14, 15,...
              17, 18, 19, 20, 21,...
              23, 24, 25, 26,...
              28, 29, 30,...
              32, 33,...
              35];                                                         % times 2 columns
    [iG,jG] = deal(iK(indexUniT,:),jK(indexUniT,:));                       % indexing of unique G coefficients
    indexG = sort([iG(:),jG(:)],2,'descend');                              % indexing of G entries (lower half)
    [a1, a2] = deal(reshape(indexG(:,2),36,nEl)', reshape(indexG(:,1),36,nEl)'); % auxiliary set of indice (for compute array p)
    dZdu = zeros(36,24);                                                   % build U-derivative of matrix G
    for ii = 1:24 % --------------loop on the element displacment component
        tmp = 0;Uvec = zeros(24,1);Uvec(ii,1) = 1;                         % set a unit displacement component (Uvec:element node displacement)
        se = Dmat0*B0((dN(0,0,0)*Xe)\dN(0,0,0))*Uvec;                      % stresses at the element centre (vector)
        se = [se(1),se(4),se(6);...
              se(4),se(2),se(5);...
              se(6),se(5),se(3)];                                          % stress matrix
%         for i = 1:length(xiG)
%             for j = 1:length(etaG)
%                 for k = 1:length(zetaG)
%                     xi_ = xiG(i);eta_ = etaG(j); zeta_ = zetaG(k);         % current integration point
%                     w = wxi(i)*weta(j)*wzeta(k)*det(dN(xi_,eta_,zeta_)*Xe);% current integration weight
%                     gradN = (dN(xi_,eta_,zeta_)*Xe)\dN(xi_,eta_,zeta_);    % x-, y-, z-derivatives of shape function (shape function physical derivatives)
%                     B1 = [kron(gradN,[1,0,0]);kron(gradN,[0,1,0]);kron(gradN,[0,0,1])];% deformation gradient
%                     tmp = tmp+w*(B1'*kron(eye(3),se)*B1);                  % current contribution to dG0/du_i
%                 end
%             end
%         end
        for i = 1:length(XIG(:))
           xi_ = XIG(i); eta_ = ETAG(i);zeta_ = ZETAG(i);         % current intergration point
           w = WXI(i)*WETA(i)*WZETA(i)*det(dN(xi_,eta_,zeta_)*Xe);% current intergration weight
           gradN = (dN(xi_,eta_,zeta_)*Xe)\dN(xi_,eta_,zeta_);    % x-, y-, z-derivatives of shape function (shape function physical derivatives)
           B1 = [kron(gradN,[1,0,0]);kron(gradN,[0,1,0]);kron(gradN,[0,0,1])];% deformation gradient
           tmp = tmp+w*(B1'*kron(eye(3),se)*B1);                  % current contribution to dG0/du_i
        end
        dZdu(:,ii) = tmp(indexUniF)';                                      % extract indenpendent coefficients
    end
    dZdu(t2index,:) = 2*dZdu(t2index,:);                                   % times 2 for v-m-v product (ZTilde)
    fKS = @(p,v)max(v)+log(sum(exp(p*(v-max(v)))))/p;                      % K-S aggregation function
    dKS = @(p,v,dv)sum(exp(p.*(v-max(v)))'.*dv,2)./sum(exp(p*(v-max(v)))); % derivative of K-S function
end % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#B#
%% ---------------------PRE. 3) Loads, Supports And Passive Domains
%% Supports #CSWs#
fixedNodeSet1 = reshape(nodeNrs(:,1,:),[],1);                              % nodes in faces ABB'A'
fixedNodeSet2 = reshape(nodeNrs(:,end,:),[],1);                            % nodes in faces CDD'C' 
fixedNodeSet3 = reshape(nodeNrs(:,:,[1 end]),[],1);                        % nodes in faces ABCD & A'B'C'D'
fixedSet1 = [3*fixedNodeSet1;3*fixedNodeSet1-1;3*fixedNodeSet1-2;3*fixedNodeSet2-2];% restrain x-, y- and z-translations (ABB'A') and xtranslations (CDD'C')
fixedSet2 = [3*fixedNodeSet3-1;3*fixedNodeSet3-2];                         % restrain x- and y-translations ABCD & A'B'C'D' 
fixed = union(fixedSet1,fixedSet2);                                        % restrained DOFs 
free = setdiff(int32(1:nDof),fixed);                                       % set of free DOFs
clear fixed fixedNodeSet1 fixedNodeSet2 fixedNodeSet3 fixedSet1 fixedSet2 
%% Loads
waveform = xPhysMat(:,:,floor(nelz/2));
elLoadX = find(waveform(:,end)>0);elLoadZ = 1:nelz;
loadDof = squeeze(3*nodeNrs(elLoadX(1):elLoadX(end)+1,end,elLoadZ(1):elLoadZ(end)+1)); % Load DOFs !!!This line was modfied 2024/09/29!!!
nelF = size(loadDof)-1;                                                    % number of elements in the x,y,z direction of load area
magF = 1;                                                                  % total magnitude of load
modF = magF/(nelF(1)*nelF(2))*ones(size(loadDof));clear nelF               % modulus of the force density (force at single node)
[modF(1,:),modF(end,:)] = deal(modF(1,:)/2,modF(end,:)/2);
[modF(:,1),modF(:,end)] = deal(modF(:,1)/2,modF(:,end)/2);                 % consistent load on edge&corner nodes 
F = fsparse(loadDof(:),1,-modF(:),[nDof,1]);clear modF                     % Define Load Vector using fsparse
%% Active Domains
[pasS,pasV] = deal([],[]);                                                 % define passive domains S-solid x = 1, V-void x = 0
act = setdiff(int32((1:nEl)'),union(pasS,pasV));                           % set of active design variables
%% -----------------PRE. 4) Prepare Filter And Projection Operators
if filterBC == 'N', bcF = 'symmetric';else, bcF = 0;end                    % select filter BC (Neumann N;Dirichlect D)
tmpFit = -ceil(rmin)+1:ceil(rmin)-1;
[dx,dy,dz] = ndgrid(tmpFit,tmpFit,tmpFit);clear tmpFit                     % x-, y-, and z-distance between element e and i
h = max(0,rmin - sqrt(dx.^2+dy.^2+dz.^2));                                 % convolution kernel
Hs = imfilter(ones(nelx,nely,nelz),h,bcF);                                 % matrix of weights
dHs = Hs;                                                                  % initialize dHs
prj = @(v,eta,beta) (tanh(beta*eta)+tanh(beta*(v(:)-eta)))./...
    (tanh(beta*eta)+tanh(beta*(1-eta)));                                   % relaxed Heaviside projection
deta = @(v,eta,beta) -beta*csch(beta).*sech(beta*(v(:)-eta)).^2 .*...
    sinh(v(:)*beta).*sinh((1-v(:))*beta);                                  % projection eta-derivative
dprj = @(v,eta,beta) beta*(1-tanh(beta*(v-eta)).^2)./...
    (tanh(beta*eta)+tanh(beta*(1-eta)));                                   % projection x-derivative
%% ----------------PRE. 5) Allocate And Initialize Other Parameters
[x,dsK,dV] = deal(zeros(nEl,1));                                           % initialize vectors of size nEl*1
if any(prSel{1}=='B')
    [dsG,dmKS] = deal(zeros(nEl,1));                                       % initialize vectors of size nEl*1
    [phiDKphi,phiDGphi,adj] = deal(zeros(nEl,nEig));                       % initialize matrixes of size nEl*nEig
    phi=zeros(nDof,nEig); adjL = phi; adjV = phi;                          % initialize arrays of size nDof*nEig
    [plotL,plotR,muVec] = deal([],[],[]);
end
U=zeros(nDof,1);                                                           % initialize arrays of size nDof*1
dV(act,1) = 1/nEl;                                                         % derivative of volume fraction
[xpOld,loop,restartAs,ch] = deal(0,0,0,1);
if nargin > 16                                                             % additional paramter initialData
    load(initialData);x = xPhysMat(:);                                     % initialize design from saved data
else
    x(act)= volfrac*(nEl-length(pasV)-length(pasS))/length(act);           % volume fraction on "active" set
    x(pasS) = 1;                                                           % set x = 1 on "passive solid" set
end
xPhys = x; clear iK jK iG jG dx dy dz;                                     % initialize xPhys and free memory
fprintf('Initializaion took %0.6e s\n',toc);
% -------------------------------------------------Start Optimaization Loop
while loop < itmax && ch > 1e-6                                            % Two conditions for exiting the iteration
    T0 = cputime;tic;
    loop = loop + 1;                                                       % Update iteration counter
    %% -----------------------RL. 1) Compute Physical Density Field
    xTilde = imfilter(reshape(x,nelx,nely,nelz),h,bcF)./Hs;                % compute filtered field
    xPhys(act) = xTilde(act); clear xTilde                                 % modify active elements only
    if filter > 1                                                          % ft > 1, apply projection ()
        f = (mean(prj(xPhys,eta,beta))-volfrac)*(filter==3);               % function (volume of x-projected)
        while abs(f) > 1e-6 && prSel{1}(1) ~= 'V'                          % Newton loop for finding optimazied eta
            eta = eta - f/mean(deta(xPhys(:),eta,beta));
            f = mean(prj(xPhys,eta,beta)) - volfrac;
        end
        dHs = Hs./reshape(dprj(xPhys,eta,beta),nelx,nely,nelz);            % modification of the sensitivity
        xPhys = prj(xPhys,eta,beta);                                       % compute projected field
    end
    ch = max(abs(xPhys - xpOld));xpOld = xPhys;
    %% ----------------RL. 2) Setup And Solve Equilibrium Equations
    sK = (Emin+xPhys.^penalK*(E0-Emin));                                   % stiffness interpolation EK
    dsK(act) = penalK*(E0-Emin)*xPhys(act).^(penalK-1);                    % derivative of stiffness interpolation EK
    sK = reshape(Ke(:)*sK',length(Ke)*nEl,1);                              % lower parts of stiffness Ke after interpolation
    K = fsparse(indexK(:,1),indexK(:,2),sK,[nDof,nDof]);                   % assemble stiffness matrix fsparse
    K = K + K'- diag(diag(K));                                             % symmetrization of K
    dK = decomposition(K(free,free),'chol','lower');                       % decompose K and store factor
    U(free) = dK \ F(free);                                                % solve equilibrium system (a faster method)
    dc = -dsK.*sum((U(cMat)*Ke0).*U(cMat),2);                              % compute compliance sensitivity
    if any(prSel{1}=='B')% >>>>>>>>>>Perform Only If Buckling Is Active #B#
        %% --------------------RL. 3) Build Stress Stiffness Matrix
        sGP = (Dmat0*B0((dN(0,0,0)*Xe)\dN(0,0,0))*U(cMat)')';              % stresses at element centroids
        l = [1,2,3,4,5,6,7,8,2,3,4,5,6,7,8,3,4,5,6,7,8,4,5,6,7,8,5,6,7,8,6,7,8,7,8,8;...
             1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,5,5,5,5,6,6,6,7,7,8]';
        Z = zeros(nEl,36);                                                 % allocate array for compact storage of Ge coefficient
%         for i = 1:length(xiG)                 % loop over hex points
%             for j = 1:length(etaG)
%                 for k = 1:length(zetaG)
%                     xi_ = xiG(i); eta_ = etaG(j);zeta_ = zetaG(k);         % current intergration point
%                     w = wxi(i)*weta(j)*wzeta(k)*det(dN(xi_,eta_,zeta_)*Xe);% current intergration weight
%                     % --reduced represenation of strain-displacement matrix
%                     gradN = (dN(xi_,eta_,zeta_)*Xe)\dN(xi_,eta_,zeta_);    % shape function physical derivatives (x-,y- and z-derivatives)
%                     a_ = gradN(1,:); b_ = gradN(2,:);c_ = gradN(3,:);
%                     B = zeros(6,36);
%                     for jj = 1:36
%                         B(:,jj) = [a_(l(jj,1))*a_(l(jj,2));...
%                             b_(l(jj,1))*b_(l(jj,2));...
%                             c_(l(jj,1))*c_(l(jj,2));...
%                             b_(l(jj,2))*a_(l(jj,1))+ b_(l(jj,1))*a_(l(jj,2));...
%                             b_(l(jj,2))*c_(l(jj,1))+ b_(l(jj,1))*c_(l(jj,2));...
%                             c_(l(jj,2))*a_(l(jj,1))+ c_(l(jj,1))*a_(l(jj,2))];
%                     end
%                     Z = Z+sGP*B*w;                                         % current contribution to (unique ~= 0) elements of Ge
%                 end
%             end
%         end
        % Modified Version 20241013-15:46
        for i = 1:length(XIG(:))
            xi_ = XIG(i); eta_ = ETAG(i);zeta_ = ZETAG(i);         % current intergration point
            w = WXI(i)*WETA(i)*WZETA(i)*det(dN(xi_,eta_,zeta_)*Xe);% current intergration weight
            % --reduced represenation of strain-displacement matrix
            gradN = (dN(xi_,eta_,zeta_)*Xe)\dN(xi_,eta_,zeta_);    % shape function physical derivatives (x-,y- and z-derivatives)
            a_ = gradN(1,:); b_ = gradN(2,:);c_ = gradN(3,:);
            B  =   [a_(l(:,1)).*a_(l(:,2));...
                    b_(l(:,1)).*b_(l(:,2));...
                    c_(l(:,1)).*c_(l(:,2));...
                    b_(l(:,2)).*a_(l(:,1))+ b_(l(:,1)).*a_(l(:,2));...
                    b_(l(:,2)).*c_(l(:,1))+ b_(l(:,1)).*c_(l(:,2));...
                    c_(l(:,2)).*a_(l(:,1))+ c_(l(:,1)).*a_(l(:,2))];
            Z = Z+sGP*B*w;                                         % current contribution to (unique ~= 0) elements of Ge
        end
        sG0 = E0*xPhys.^penalG;                                            % stress interpolation EG
        dsG(act) = penalG * E0 * xPhys(act).^(penalG - 1);                 % derivative of stress interpolation EG
        sG = reshape((sG0.*Z)',36*nEl,1);                                  % unique coefficients of the stress stiffness Ge
        G = fsparse(indexG(:,1)  ,indexG(:,2)  ,sG,[nDof,nDof])+...
            fsparse(indexG(:,1)+1,indexG(:,2)+1,sG,[nDof,nDof])+...
            fsparse(indexG(:,1)+2,indexG(:,2)+2,sG,[nDof,nDof]);           % assemble global stress matrix G (fsparse)
        G = G + G' - diag(diag(G));                                        % symmetrization of G
        %% --------------- RL. 4) Sovle Buckling Eigenvalue Problem
        matFun = @(x) dK\(G(free,free)*x);                                 % matrix action function
        [eivecs,D] = eigs(matFun, length(free),nEig+4,'sa','Tolerance',1e-8);% compute eigenvalues
        [mu,ii] = sort(diag(-D),'descend');                                % sort of eigenvalues (mu = -D(i))
        eivSort = eivecs(:,ii(1:nEig));                                    % sort eigenvectors accordingly
        phi(free,:) = eivSort./sqrt(diag(eivSort'*K(free,free)*eivSort)'); % orthonomalize (phi'*K*phi = 1)
        %% -----------------------RL. 5) Sensitive Analysis Of BLFs
        dkeG = dsG .* Z;                                                   % x-derivative of Ge
        dkeG(:,t2index) = 2*dkeG(:,t2index);                               % *2 columns for v-m-v product
        for j = 1:nEig % loop on the eigenvalues included in the optimazation
            % 1) Term due to the elastic stiffness martirx (phi'*dK/dx*phi)
            t = phi(:,j);
            phiDKphi(:,j) = dsK.*sum((t(cMat)*Ke0).*t(cMat),2);
            % 2) Term due to the geometric stiffness martirx (phi'*dG/dx*phi)
            p = t(a1).*t(a2)+t(a1+1).*t(a2+1)+t(a1+2).*t(a2+2);
            phiDGphi(:,j) = sum(p.* dkeG,2);
            % 3)-------------------------------------Setup of adjoint loads
            tmp = zeros(nDof,1);
            for k = 1:24 %------contribution of each term dKg/du_i i = 1:nD
                tmp(cMat(:,k)) = tmp(cMat(:,k)) + (sG0.*p)*dZdu(:,k);
            end
            adjL(:,j) = tmp;
        end
        % --------Solve the adjoint problem and compute the term U'*dK/dx*V
        adjV(free,:) =dK\adjL(free,:);                                     % Method decomposition (faster)
        for j = 1:nEig
            vv = adjV(:,j);
            adj(:,j) = dsK .* sum((U(cMat)*Ke0).*vv(cMat),2);
        end
        % ---------Overall sensitivity expresstion for the "mu" eigenvalues
        dmu = -(phiDGphi + mu(1:nEig)' .* phiDKphi -adj );
    end%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<#B#
    %% ------------RL. 6) Select Objective Function And Constraints
    if loop==1, c0=F'*U; v0=mean(xPhys(:));end                             % initial compliance & volume fraction
    c = F'*U; v = mean(xPhys(:));                                          % current compliance & volume fraction
    switch prSel{1}
        case ['C','V'] %---------minimize compliance with volume constraint
            g0 = c/c0;
            dg0 = imfilter(reshape(dc/c0,nelx,nely,nelz)./dHs,h,bcF);      dg0 = repmat(sum(dg0,3),1,1,nelz);
            g1 = v/volfrac - 1;
            dg1 = imfilter(reshape(dV/volfrac,nelx,nely,nelz)./dHs,h,bcF); dg1 = repmat(sum(dg1,3),1,1,nelz);
            g1Vec = g1;                                                    % use for MMA
            dgiMat = dg1(:);                                               % use for MMA
        case ['V','C']% ---------minimize volume with compliance constraint
            g0 = mean(xPhys(:))./v0;
            dg0 = imfilter(reshape(dV/v0,nelx,nely,nelz)./dHs,h,bcF);      dg0 = repmat(sum(dg0,3),1,1,nelz);
            g1 = c/(prSel{2}*c0)-1;
            dg1 = imfilter(reshape(dc/(prSel{2}*c0),nelx,nely,nelz)./dHs,h,bcF);dg1 = repmat(sum(dg1,3),1,1,nelz);
            g1Vec = g1;                                                    % use for MMA
            dgiMat = dg1(:);                                               % use for MMA
        case ['B','C','V']% maximize BLF with compliance & volume constraints
            if loop==1, muKS0=fKS(pAgg,mu(1:nEig));g0=1;cMax=prSel{2}(1);
            else, g0=fKS(pAgg,mu(1:nEig))/muKS0;end % KS aggregation of mu (=1/lambda)
            dmKS = dKS(pAgg,mu(1:nEig),dmu);
            dg0 = imfilter(reshape(dmKS/muKS0,nelx,nely,nelz)./dHs,h,bcF); dg0 = repmat(sum(dg0,3),1,1,nelz);
            % -Constraint function: KS aggregation of compliance and volume
            g1Vec = [c;v]./[cMax*c0;volfrac] - 1;
            dg1c = imfilter(reshape(dc/(cMax*c0),nelx,nely,nelz)./dHs,h,bcF);dg1c = repmat(sum(dg1c,3),1,1,nelz);
            dg1V = imfilter(reshape(dV/volfrac,nelx,nely,nelz)./dHs,h,bcF);dg1V = repmat(sum(dg1V,3),1,1,nelz);
            dgiMat = [dg1c(:),dg1V(:)];                                    % use for KS or MMA
            g1 = fKS(pAgg,g1Vec);
            dg1 = dKS(pAgg,g1Vec,dgiMat);
            plotL(loop,:) = [1/g0/muKS0,1/mu(1)];strL='KS(-),\lambda_1(--)';
            plotR(loop,:) = [g1,g1Vec'];strR='g_1(-),gC(--),gV(.-)';
            muVec(loop,:) = mu';
        case ['V','C','B']% minimize volume with compliance & BLF constraints
            g0 = v/v0;
            dg0 = imfilter(reshape(dV/volfrac,nelx,nely,nelz)./dHs,h,bcF);
            % ---Constraint function: KS aggregation of BLFs and compliance
            muKS = fKS(pAgg,mu(1:nEig));
            dmKS = dKS(pAgg,mu(1:nEig),dmu);
            g1Vec = [prSel{2}(2)*muKS;c]./[1;prSel{2}(1)*c0] - 1;
            dg1c = imfilter(reshape(dc/(c0*prSel{2}(1)),nelx,nely,nelz)./dHs,h,bcF);
            dg1l = imfilter(reshape(dmKS*prSel{2}(2),nelx,nely,nelz)./dHs,h,bcF);
            dgiMat = [dg1l(:),dg1c(:)];
            g1 = fKS(pAgg,g1Vec);
            dg1 = dKS(pAgg,g1Vec,dgiMat);
            plotL(loop,:) = g0; strL = 'g_0';                              % store objective value
            plotR(loop,:) = [g1,g1Vec'];strR = 'g_1(-),gL(--),gc(.-)';
            muVec = cat(1,muVec,mu');                                      % store eigenvalues
    end
    %% ------------------------------RL. 7) Update Design Variables
    lm = ones(size(g1Vec))*nan;
    switch updateRule
        case 'ocUpdate' % <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<ocUpdate
            if loop==1, xOld = x(act); xOld1 = xOld; as = []; end          % initialize MMA history parameters
            [x0,as,lmid(1)] = ocUpdate(loop,x(act),dg0(act),g1,dg1(act),ocParam,xOld,xOld1,as,beta,restartAs);
            xOld1 = xOld;xOld = x(act);x(act)=x0; % >>>>>>>>>>>>>>>ocUpdate
        case 'MMA'% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<MMA
            if loop==1, xOld = x(act); xOld1 = xOld; low = []; upp = []; end          % initialize MMA history parameters
            [xMMA,~,~,lmid,~,~,~,~,~,low,upp] = mmasub_new(length(g1Vec),length(act),loop,...
                x(act),max(0,x(act)-ocParam(1)),min(1,x(act)+ocParam(1)),xOld,xOld1,...
                g0,dg0(act),g1Vec,dgiMat(act,:)',low,upp,1,zeros(size(g1Vec)),ones(size(g1Vec))*1e2,ones(size(g1Vec)),beta,restartAs);
            xOld1 = xOld;xOld = x(act);x(act)=xMMA;
    end
    if length(lmid) < 2, lm = [lmid;nan]; else,lm = lmid;end
    %% ---------RL. 8)Post-processing, Save, Plot and Print Results
    % Post-processing
    [Ux,Uy,Uz] = deal(U(1:3:end),U(2:3:end),U(3:3:end));                   % extract displacement component
    Utotal = sqrt(Ux.^2 + Uy.^2 + Uz.^2);
    scaleAuto = maxDU/max(abs(U),[],'all');
    nodes_deformed = nodes + scaleAuto*[Ux,Uy,Uz];
    vertices_var =  Utotal;
    eIndex = xPhys>threshold;
    faces = reshape(elements(eIndex,efaces')',4,[])';
    % save physical density
    xPhysMat = reshape(xPhys,nelx,nely,nelz);
    save([fileTitle '_xPhysMat_It' num2str(loop) '.mat'],'xPhysMat');
    % save displacment U
    save([fileTitle '_U_It' num2str(loop) '.mat'],'U');
    if isPlotCurrentDesign == 1% <<<<<<<<<<<<<<<<<<<<< figure CurrentDesign
        % calculated isometric surfaces and end caps.
        isovals = smooth3(xPhysMat,'box',1);
        [isof1,isov1] = isosurface(centerX,centerY,centerZ,isovals,0.5);
        [isof2,isov2] = isocaps(centerX,centerY,centerZ,isovals,0.5);
        figure(f_CurrentDesign);
        f_CurrentDesign.Name = [meshtxt '_' prSel{1} '_Current Design_It' num2str(loop)];
        tiledlayout(1,2,'TileSpacing','tight','Padding','tight');
        nexttile(1);cla;% ------------------------------------------subfig1
        patch('Faces',isof1,'Vertices',isov1,'FaceColor','b','EdgeColor','none');
        patch('Faces',isof2,'Vertices',isov2,'FaceColor','#D95319','EdgeColor','none');
        title(['It.' num2str(loop)]); view([135,25]);axis equal;grid on;
        xlabel('X');ylabel('Y');zlabel('Z');xlim(axislim(1,:));ylim(axislim(2,:));zlim(axislim(3,:));
        nexttile(2);cla;% ------------------------------------------subfig2
        patch('Faces',faces,'Vertices',nodes_deformed,'FaceVertexCData',vertices_var,'FaceColor','interp','Edgecolor','black');
        title({['Deformed It.' num2str(loop,"%3d")]; ['U scale factor = ' num2str(scaleAuto)]});
        grid on;axis equal;view([135 25]);
        xlabel('X');ylabel('Y');zlabel('Z');xlim(axislim(1,:));ylim(axislim(2,:));zlim(axislim(3,:));
    end % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> figure CurrentDesign
    if isPlotWaveform == 1% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<figure Waveform
    figure(f_Waveform);f_Waveform.Name = [meshtxt '_' prSel{1} '_Waveform_It' num2str(loop)];
    cla;imagesc(1-xPhysMat(:,:,floor(nelz/2)));colormap('gray');axis equal tight off; 
    end% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>figure Waveform 
    if any(prSel{1} == 'B') % Plot design, g0&g1 evolution, BLFs evolution
        % Save
        if loop == 1 || mod(loop,phiOutputInv) == 0,save([fileTitle '_phi_It' num2str(loop) '.mat'],'phi');end% save current eigenvector
        save([fileTitle '_mu.mat'],'muVec');                               % update ever loop
        save([fileTitle '_gi.mat'],'plotL','strL','plotR','strR');         % update ever loop
        % Plot % <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<figure History
        figure(f_History);
        f_History.Name = [meshtxt '_' prSel{1} '_Obj&ConstValues_It' num2str(loop)];
        tiledlayout(1,2,'TileSpacing','tight','Padding','tight');
        nexttile(1);cla;% ------------------------------------------subfig1
        yyaxis left; plot(1:loop,plotL);ylabel(strL);
        yyaxis right; plot(1:loop,plotR);ylabel(strR);title('Objective and constraint');
        nexttile(2);cla;% ------------------------------------------subfig2
        plot(1:loop,1./muVec(:,1:4));title('Lowest BLFs');
        % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>figure History
        if isPlotBucklingMode == 1% <<<<<<<<<<<<<<<<<<<<figure BucklingMode
            figure(f_BucklingMode);
            f_BucklingMode.Name = [meshtxt '_' prSel{1} '_Buckling Mode_It' num2str(loop)];
            tiledlayout(2,2,"TileSpacing","tight",'Padding','tight');
            for ppp = 1:4
                nexttile(ppp);cla;
                [phix,phiy,phiz] = deal(phi(1:3:end,ppp),phi(2:3:end,ppp),phi(3:3:end,ppp));
                vertices_var =  sqrt(phix.^2 + phiy.^2 + phiz.^2);
                scaleAuto = maxDU/max(abs(phi),[],'all');
                nodes_deformed = nodes + scaleAuto*[phix,phiy,phiz];
                patch('Faces',faces,'Vertices',nodes_deformed,'FaceVertexCData',vertices_var,'FaceColor','interp','Edgecolor','black');
                title(['Mode ' num2str(ppp)]);grid on;axis equal;view([90 0]);
                xlabel('X');ylabel('Y');zlabel('Z');xlim(axislim(1,:));ylim(axislim(2,:));zlim(axislim(3,:));
            end
        end % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>figure BucklingMode     
    end
    drawnow;
    % Print iteration infomation
    fprintf(topHistoryFileID,'%3i\t%6.4f\t%+0.2e\t%4.1f\t%4.1f\t%4.2f\t%4.1f\t%0.3e\t%0.3e\t%0.3e\t%0.3e\t%0.6e\t%6.4f\t%0.6e\t%0.3e\r\n',...
             loop,    g0,         g1,         penalK,      penalG,      eta,        beta,        pAgg,        ch,        lm(1),       lm(2),     c,        v,              cputime-T0,      toc);    
    fprintf('It.:%3i | g0:%6.4f | g1:%+0.2e | penK:%4.1f | penG:%4.1f | eta:%4.2f | beta:%4.1f | pAgg:%0.3e | ch:%0.3e | lm1:%0.3e | lm2:%0.3e | c:%0.2e | Volfrac:%6.4f | CPU time:%0.3e | Real time:%0.3e\n',...
             loop,    g0,         g1,         penalK,      penalG,      eta,        beta,        pAgg,        ch,        lm(1),       lm(2),       c,        v,              cputime-T0,      toc);% print result
    %% RL.END)apply continuation on penalization(s), beta & aggregation parameter(s)
    penalKold = penalK; penalGold = penalG; betaOld = beta;
    [penalK,penalG,beta,pAgg] = deal(cnt(penalK, penalCntK, loop),...
        cnt(penalG,penalCntG,loop),cnt(beta,betaCnt,loop),cnt(pAgg,pAggCnt,loop));
    if (beta-betaOld~= 0||penalK-penalKold~=0||penalG-penalGold~=0)
        restartAs = 1; else, restartAs = 0; end 
end
%%
fclose(topHistoryFileID);
saveas(f_History,[fileTitle '_It' num2str(loop) '_Obj&ConstValues.fig']);
if isPlotCurrentDesign == 1, saveas(f_CurrentDesign,[fileTitle '_It' num2str(loop) '.fig']);end
if isPlotWaveform == 1, saveas(f_Waveform,[fileTitle '_It' num2str(loop) '_Waveform.fig']);end
if any(prSel{1} == 'B')&&isPlotBucklingMode,saveas(f_BucklingMode,[fileTitle '_It' num2str(loop) '_Buckling Mode.fig']);end
end