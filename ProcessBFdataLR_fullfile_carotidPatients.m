%%% Script family members: CarotidFlow3TxL11_Beamform.m          (Beamformes and saves)
%%%                        CarotidFlow3TxL11_ProcessBFdata.m     (Loads bf data and does vector flow processing. Saves to file)
%%%                        CarotidFlow3TxL11_ExportData.m        (Exports data for particle visualization or USFig images)
%%%                        CarotidFlow3TxL11_Quantify.m          (PW Doppler and velocity traces)
% Ingvild K Ekroll, October 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
%clear variables
laptop = 0;
server = 1;

addpath(genpath('D:\MasterData\TilEmil\DopplerFunctions\'))
%addpath(genpath('C:\Users\ingvilek\utils\'))
%rmpath('C:\Users\ingvilek\vectorFlow\Code\ScannerCode\Tissue\')
%rmpath('C:\Users\ingvilek\usfig\')


%%{
recNo =1;
CarotidD = carotidPatientsL();
CarotidD = CarotidD{recNo};
datapath = CarotidD.data_proc_path %CarotidD.fastdrive_path  %CarotidD.data_proc_path;
metapath = CarotidD.data_main_path;
if isfield(CarotidD,'savepath')
    savepath = CarotidD.savepath;
else
    savepath = CarotidD.data_proc_path; 
end

if isfield(CarotidD,'powerthresh')
    powerthresh = CarotidD.powerthresh;
else
    powerthresh = 23;%NB! different powerthresh
end

application = 'carotid';
%}

filetype = CarotidD.fileType; %'raw'; % 'raw' or 'h5'
%g = gpuDevice;

load(fullfile(datapath,'bfPars'))
load(fullfile(metapath,'meta'))

tissueAdaptive = 0;
if tissueAdaptive
    inputVelocity = load(fullfile(datapath,'adaptiveFilterInputVelocities'));
     t_inputVelocity = linspace(0,length(inputVelocity.maxImageVels)/bfPars.PRF,length(inputVelocity.maxImageVels));
    inputVelocityF = inputVelocity.maxImageVels(inputVelocity.TissuePacket:90:end,:);
    tF = t_inputVelocity(inputVelocity.TissuePacket:90:end);
    figure()
    plot(inputVelocityF,'o')
end

useMemmap = 0;

%%   Load all data into one variable
nx = length(bfPars.x_axis);
nz = length(bfPars.z_axis);
nAngles = length(bfPars.rxVals);
bfPars.txRxSetName = CarotidD.txRxSetName;

% tStart = 0.5;
% tEnd = 1;

startFrame = 1 %CarotidD.frame_start;%floor(tStart*bfPars.PRF/CarotidD.packetsize) %CarotidD.frame_start;
endFrame = 3 %CarotidD.frame_end;% ceil(tEnd*bfPars.PRF/CarotidD.packetsize) %CarotidD.frame_end;
numSavedFrames = endFrame-startFrame+1;

LowResData = zeros( nz, nx, nAngles, CarotidD.packetsize*numSavedFrames );

if recNo == 201 %% Temporary hack
   
    load(fullfile(datapath,'LowResData')); % Loads previously file containing all timesteps
    
else
    ctr = 1;
    for ff = startFrame:endFrame
        fnameR = fullfile(datapath,sprintf([CarotidD.iq_file_suffix '_%i_r'],ff));
        fnameI = fullfile(datapath,sprintf([CarotidD.iq_file_suffix '_%i_i'],ff));
        tic
        if useMemmap, 
            mR = memmapfile( fnameR, 'Format', {CarotidD.bfFormat, [nz, nx, nAngles, CarotidD.packetsize], 'x'} );
            mI = memmapfile( fnameI, 'Format', {CarotidD.bfFormat, [nz, nx, nAngles, CarotidD.packetsize], 'x'} );

            LowResData(:,:,:,((ctr-1)*CarotidD.packetsize+1):ctr*CarotidD.packetsize) = single(mR.Data.x + 1i*mI.Data.x);
        else
            fid = fopen( fnameR, 'r');
            tempRdata = fread( fid, 'single');
            fclose(fid);
            fid = fopen( fnameI, 'r');
            tempIdata = fread( fid, 'single');
            fclose(fid);
            
            tempRdata2 = reshape( tempRdata, [nz, nx, nAngles, CarotidD.packetsize] );
            tempIdata2 = reshape( tempIdata, [nz, nx, nAngles, CarotidD.packetsize] );
            
            clear tempRdata;
            clear tempIdata
            for kk = 1:size( LowResData,3),
                LowResData(:,:,kk,((ctr-1)*CarotidD.packetsize+1):ctr*CarotidD.packetsize) = tempRdata2(:,:,kk,:) + 1i*tempIdata2(:,:,kk,:);
            end
        end
        toc
        ctr = ctr+1;
    end
end


%% Hack to look at data - find interesting frame
while 0
    figure(1);
    subplot(2,1,1); imagesc(20*log10(abs(squeeze(LowResData(:,:,5,100)))));
    [xin zin] = ginput(1);
    zin = floor(zin);
    xin = floor(xin);
    PWsig = squeeze(LowResData(zin,xin,5,:));

    skip=10;
    Nv=128;
    w=hamming(Nv);
    Nsmooth=3;

    Nt=length(PWsig);
    Nw=length(w);
    t=1:skip:Nt-Nw;

    iqM=zeros(Nw,length(t));
    for n=1:length(t),
        It=t(n):t(n)+Nw-1;
        iqM(:,n)=PWsig(It).*w;
    end;
    A=fft(iqM,Nv);
    P=A.*conj(A);
    P=filter2(ones(1,Nsmooth)/Nsmooth,P);
    PWspect = repmat( fftshift(P,1), 3, 1);

    vNyq = bfPars.PRF*bfPars.c/4/bfPars.f_demod;
    vaxis = linspace(-3*vNyq, 3*vNyq, size(PWspect,1)+1 ); vaxis(end) = [];

    figure(1); subplot(2,1,2); imagesc(t/bfPars.PRF,vaxis, 10*log10( abs( PWspect) ) );
    %caxis([-dyn 0]-gain);
    colormap( gray)
    xlabel('Time [s]');
    ylabel('Velocity [m/s]');
    title(sprintf('PW Doppler'));
end

%% Pick out frame for inspection
if 0 % all data are loaded
    tLook = 0.65
    fLook = tLook*bfPars.PRF
end

if 1
   tLook = 1.2;
   tStart = startFrame*CarotidD.packetsize/bfPars.PRF;
   t = tLook-tStart
   fLook = t*bfPars.PRF;
end


%% Sliding window processing parameters
processPacket = 270;
skip = 20;
nrLRFrames = size(LowResData,4);
maxNumProcessPackets = length(((processPacket/2):skip:(nrLRFrames-processPacket/2)));%length(1:skip:(nrLRFrames-processPacket))+1;

startPacket =1; %ceil(fLook/processPacket)
nrPackets = 1%maxNumProcessPackets; %1 

procParams.processPacket = processPacket;
procParams.skip = skip;
%procParams.vdtime = ((processPacket/2):skip:(nrLRFrames-processPacket/2))/bfPars.PRF;
procParams.vdtime = ((processPacket/2):skip:skip*(maxNumProcessPackets-1)+(processPacket/2))/bfPars.PRF;
%% Sliding window processing of data
Fctr = 1;
for kk = 1:1:nrPackets  %sliding window processing
   % packetData = LowResData(:,:,:,(startPacket-1)*processPacket +1: startPacket*processPacket +nrPackets*processPacket); %sliding window processing
    packetData = LowResData(:,:,:,(startPacket-1)*processPacket +1 +(kk-1)*skip : startPacket*processPacket +(kk-1)*skip); %sliding window processing

    
    %% Some settings
    % General
    showCompoundedFrames = 0; % Figure displaying all "B-mode" frames
    showHighPassFrames = 0; % Figure displaying all frames used for speckle tracking (after high pass filtering)
    useMedianFilter = 1;

    % Do you want the full aliasing correction procedure, or only simple least squares algorithm?
    doAliasingCorrection = CarotidD.doAliasingCorrection;

    if nrPackets<maxNumProcessPackets
        saveForMovie = 0;
    else
        saveForMovie = 1;
        disp('Saving data for later use');
    end
    
    % Beamforming parameters
    % bfPars;  NB: Comes with loaded data
    
    % Doppler ROI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    useMaxROI = 0; % If both options are 0, interactively chose ROI.
    useTrapezoidROI = 1;
   
    %¤¤¤¤¤¤¤ USED IF useMaxROI: ¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤ 
    switch application
        case 'carotid'
            MaxROIDepth = CarotidD.ROIdepthEnd; % [m]
            MinROIDepth = CarotidD.ROIdepthStart; % [m]
        case 'rat'
            MaxROIDepth = 0.025; % [m]
            MinROIDepth = 0.005; % [m]      
    end
    %¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤

    % Color Doppler/vector Doppler settings
    dropSamples = 1; %For R1 estimation, avoid clutter filter edge effects
    PowerThreshAfterHp = powerthresh; % For segmentation, in dB

    % Specle tracking settings
    useCorrThresh = 1;
    correlationThreshold = 0.1;

    % angle vectors
    txangles = bfPars.txVals;
    rxangles = bfPars.rxVals;

    % Filter specs
    hpPars.packetSize = processPacket;
    hpPars.typeFilter = 'FIR138_3'%CarotidD.typeFilter %'SVD3D' %CarotidD.typeFilter;
    hpPars.PRF = bfPars.PRF;
    hpPars.f_demod = bfPars.f_demod;
    hpPars.c = bfPars.c;
    if strcmpi(hpPars.typeFilter,'polyreg')
        hpPars.polyOrder = CarotidD.filterOrder; % If hpPars.typeFilter = 'polyreg'
    elseif startsWith(hpPars.typeFilter,'SVD')
        hpPars.maxClutDimFac = CarotidD.maxClutDimFac;
        hpPars.maxBloodDimFac = CarotidD.maxBloodDimFac;
        hpPars.CDdet_1D = 1; % Method for determining number of clutter space dimensions. Can be based on 1D (fast, less robust) or 2D infor (slower, more robust)     
        hpPars.estClutDim = 1;
        hpPars.estBloodDim = 0;
        hpPars.bloodFraction = 0.025;
    end
    
    %% Compound from low resolution images!
    preCData = packetData(:,:,(rxangles==0),:);
    recons_comp.data = sum( preCData, 3);
    
    remInds = CarotidD.reminds; 
    txangles(remInds) = [];
    rxangles(remInds) = [];
    packetData(:,:,remInds,:) = [];

    %% Clutter filter and R1 estimation    
    LowResR0hpMat = zeros(size(packetData,1),size(packetData,2),length(rxangles));
    if tissueAdaptive
        BWFac = 4;
        minStop = 20; %[mm/s]
        fixedStop = 172; %mm/s
    end
    
    clear b;
    for ii = 1:length(rxangles)
        if tissueAdaptive
            stopVel = max( max(inputVelocityF(LRframe,:) )*BWFac, minStop*1e-3);
            hpPars.stopband = stopVel; 
            %hpPars.stopband = fixedStop*1e-4;
        end
        
        %temp = permute( hp(permute(squeeze(packetData(:,:,ii,:)), [3 1 2]),hpPars), [2 3 1] );
       [temphp, ~, b(:,ii)] =  hp(permute(squeeze(packetData(:,:,ii,:)), [3 1 2]),hpPars);
       temp = permute(temphp,[2 3 1]);
        LowResR0hpMat(:,:,ii) = mean(conj(temp).*temp,3);
        %lowResPacketSize = size(temp,3);
        eval(sprintf('R%i = mean(conj(temp(:,:,1+dropSamples:end-2)).*temp(:,:,2+dropSamples:end-1),3);',ii))
    end
    clear temp;
    clear temphp;



    %% Clutter filter compounded data
    [iqhpComp,vThresh,b] = hp(permute(recons_comp.data(:,:,1,:), [4 1 2 3]),hpPars); 

    iqhpComp = permute(iqhpComp, [2 3 4 1]);

    maxIQHP = max(abs(iqhpComp(:)));
    maxIQ = max(abs(recons_comp.data(:)));


    if showHighPassFrames
        figure(3000);
        v = VideoWriter('hpFrames_FIR138_3_test.avi','Motion JPEG AVI');
        v.Quality = 95;
        open(v);
        for pp = 1:size(iqhpComp,4)
           titleString = sprintf('Basis for speckle tracking. Firing number %i',pp); 
           imagesc(bfPars.x_axis,bfPars.z_axis,squeeze(20*log10(abs(iqhpComp(:,:,1,pp))/maxIQHP))); caxis([-60 0])
           title(titleString);
           %pause()
           frame = getframe(gcf);
           writeVideo(v,frame);
        end
        close(v);
    end

     if showCompoundedFrames
        figure(2000);
        for pp = 1:size(recons_comp.data,4)
           titleString = sprintf('Compounded data original, Frame number %i',pp); 
           imagesc(bfPars.x_axis,bfPars.z_axis,squeeze(20*log10(abs(recons_comp.data(:,:,1,pp))/maxIQ))); caxis([-120 0])
           title(titleString);
           pause()
        end
    end
    
    R0Test = 10*log10(mean(abs(iqhpComp).^2,4));
    R0Mask = zeros(size(iqhpComp,1),size(iqhpComp,2));
    R0Mask(R0Test>PowerThreshAfterHp) = 1;

    % Compounded color flow processing --> Segmentation

    map1 = getColormapL();
    R1Test = mean(conj(iqhpComp(:,:,:,1+dropSamples:end-1)).*iqhpComp(:,:,:,2+dropSamples:end),4);
    Nx = 5;
    Nz = 5;
    myFilt = ones(Nx,Nz);
    Nthresh = 25;
    R0MaskConv =  conv2(R0Mask,myFilt,'same');

    R0MaskConv(R0MaskConv<Nthresh) = 0; 
    R0MaskConv(R0MaskConv~=0) = 1;

    R0Mask = R0MaskConv;
    R1Conv = conv2(R1Test,myFilt,'same');
    R1Masked = R1Conv.*R0Mask;

    [X,Z] = meshgrid( bfPars.x_axis, bfPars.z_axis);

    % Show image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(200); clf; 
    imagesc(bfPars.x_axis,bfPars.z_axis,angle(R1Masked)*bfPars.PRF/(2*pi)*bfPars.c/(2*bfPars.f_demod)); 
    colormap(map1);
    title('Compound CFI. "Axial velocity estimates"')
    hold on;

    if useMaxROI
        sprintf('Maximum ROI is used, based on hardcoded values')
        cfiZlimMin = MinROIDepth;
        cfiZlimMax = MaxROIDepth;
        xdrop = atan(max(txangles))*cfiZlimMax;
        cfiXlimMin = min(X(:))+xdrop; 
        cfiXlimMax = max(X(:)-xdrop);

        width = cfiXlimMax-cfiXlimMin;
        height = cfiZlimMax-cfiZlimMin;

        pos = [cfiXlimMin,cfiZlimMin,width,height];
        h = imrect(gca,pos);
        setColor(h,'y') 
        
        cfiMask = zeros(size(R0Mask));
        cfiMask(X<cfiXlimMax & X>cfiXlimMin & Z<cfiZlimMax & Z>cfiZlimMin) = 1;
        
    elseif useTrapezoidROI
        cfiMask = ones(size(R0Mask));

        for zz = 1:size(cfiMask,1)
           nzeros = ceil(zz*tan(max(abs(txangles))));
           cfiMask(zz,1:nzeros) = 0;
           cfiMask(zz,end-nzeros+1:end) = 0;
        end

        cfiMask(Z<MinROIDepth) = 0;
        
    else
        
        sprintf('Place rectangular ROI, double click when finished')
        h = imrect();
        setColor(h,'y') 
        wait(h);

        pos = getPosition(h);
        cfiZlimMin = pos(2);
        cfiZlimMax = cfiZlimMin + pos(4);
        xdrop = atan(max(rxangles))*cfiZlimMax;

        cfiXlimMinThresh = min(X(:))+xdrop; 
        cfiXlimMaxThresh = max(X(:)-xdrop);

        %Update ROI
        if abs(cfiXlimMinThresh)<abs(pos(1))
            setPosition(h,[cfiXlimMinThresh,cfiZlimMin,pos(3)-(abs(cfiXlimMinThresh-pos(1))),pos(4)])
        end

        pos = getPosition(h);

        if abs(cfiXlimMaxThresh)<abs(pos(1)+pos(3))
            setPosition(h,[pos(1),cfiZlimMin,cfiXlimMaxThresh-pos(1),pos(4)])
        end

        pos = getPosition(h);

        cfiXlimMin = pos(1);
        cfiXlimMax = pos(1)+pos(3);
        
        cfiMask = zeros(size(R0Mask));
        cfiMask(X<cfiXlimMax & X>cfiXlimMin & Z<cfiZlimMax & Z>cfiZlimMin) = 1;

        
        
        
    end


    R0Mask = R0Mask.*cfiMask;


    %% Make Matrix of Doppler data from given txrxPattern
    R1Multi = zeros(size(R1,1), size(R1,2),length(txangles)); 
    for rr = 1:length(txangles);
        eval(sprintf('R1Multi(:,:,rr) = R%i;',rr))% 
    end

    p.PRF = bfPars.PRF;

    myMask = R0Mask; % 

    clear R1Mat;
    for aa = 1:length(txangles)
        temp = R1Multi(:,:,aa);
        R1Mat(:,:,aa) = temp.*myMask;
    end

    Nx = 7;
    Nz = 7;

    myFilt = (1/(Nx*Nz))*ones(Nx,Nz);

    %Spatial averaging
    for aa = 1:length(txangles)
        R1Conv(:,:,aa) = conv2(R1Mat(:,:,aa),myFilt,'same');
    end

    R1ConvOrig = R1Conv;

    %% Vector Doppler estimation
    R1Conv = zeros( length( find( myMask == 1) ),length(txangles));
    for anglenr = 1:length(txangles),
        R1temp = R1ConvOrig(:,:,anglenr);
        R1Conv(:, anglenr) = R1temp(myMask == 1);
    end

    R1Mat = reshape(R1Conv,[1,size(R1Conv,1),size(R1Conv,2)]);
    p.txVals = txangles;
    p.rxVals = rxangles;
    p.application = application;
    p.doAliasingCorrection = doAliasingCorrection;

    tic
    [vectorEst, aliasPat, aliasMat,minRes] = vectorDoppler(R1Mat,p);
    vectorVel = vectorEst*bfPars.c/(2*bfPars.f_demod);
    t1 = toc;
    disp(sprintf('Finished with vector Doppler estimation, using %f s',t1));

    %% Resolving ambiguities

    % Some grids
    dx = bfPars.x_axis(2)-bfPars.x_axis(1);
    dz = bfPars.z_axis(2)-bfPars.z_axis(1);
    mymaskX = X(myMask == 1 );
    mymaskZ = Z(myMask == 1 );
    bfPars.mymaskX = mymaskX;
    bfPars.mymaskZ = mymaskZ;
    
    
   if doAliasingCorrection
        % Still ambiguity due to aliasing on all angles. Establish all candidate velocity vectors
        aMat = [-sin(txangles.')-sin(rxangles.') cos(txangles.')+cos(rxangles.')]./2; % Valid also when using multiple steering angles on tx and rx
        pseudoinvA = pinv(aMat);
        vecEstComp = aMat*vectorEst.';
        aliasfactTab = -2:2;
        vCand = zeros(2, size( vecEstComp, 2), length(aliasfactTab) );
        ictr = 1;
        for aliasfact = aliasfactTab,
            allCompMat = vecEstComp+ones(size(vecEstComp) )*p.PRF*aliasfact;
            vCand(:, :, ictr) = pseudoinvA*allCompMat;
            ictr = ictr+1;
        end


        %% Find correlation for all candidate velocity vectors using compounded data (GPU)   
        t2 = toc;
        ksize = 10;
        [finalEsts, corrValtot, chosenAliasPat] = blockMatching_noGPU( iqhpComp, vCand, aliasMat, aliasPat, bfPars,ksize); 
        t3 = toc;
        finalEstsX = finalEsts(1,:);
        finalEstsZ = finalEsts(2,:);
        finalEstsXNoThresh = finalEstsX;
        finalEstsZNoThresh = finalEstsZ;
        if useCorrThresh
            myMat = max( corrValtot, [], 1);
            finalEstsX(myMat < correlationThreshold ) = NaN;
            finalEstsZ(myMat < correlationThreshold ) = NaN;
        end

        disp(sprintf('Finished with ambiguity resolving, BM took %f s',t3-t2))
        
    else
        finalEstsX = vectorVel(:,1).';
        finalEstsZ = vectorVel(:,2).';
        finalEstsXNoThresh = finalEstsX;
        finalEstsZNoThresh = finalEstsZ;      
    end


    %% Interpolate to desired grids for saving and visualization
    FC = scatteredInterpolant(mymaskX,mymaskZ,double(finalEstsX.'),'linear','none');
    vxCorr = FC(X,Z);
    FC2 = scatteredInterpolant(mymaskX,mymaskZ,double(finalEstsZ.'),'linear','none');
    vzCorr = FC2(X,Z);

    FCNT = scatteredInterpolant(mymaskX,mymaskZ,double(finalEstsXNoThresh.'),'linear','none');
    vxCorrNT = FCNT(X,Z);
    FC2NT = scatteredInterpolant(mymaskX,mymaskZ,double(finalEstsZNoThresh.'),'linear','none');
    vzCorrNT = FC2NT(X,Z);
    
    
    if useMedianFilter
        vxCorrOrig = vxCorr;
        vzCorrOrig = vzCorr;
        if 0
            vxCorr = medfilt2(vxCorrOrig,[4,4],'symmetric');
            vzCorr = medfilt2(vzCorrOrig,[4,4],'symmetric');
        else
            vxCorr = mediannan(vxCorrOrig,5);
            vzCorr = mediannan(vzCorrOrig,5);
        end
        vxCorrNT = mediannan(vxCorrNT,5);
        vzCorrNT = mediannan(vzCorrNT,5);
    end
    

    if saveForMovie

            if kk == 1

                movieDataSliding.movieBmodeMat = zeros([size(R0Test,1) size(R0Test,2) maxNumProcessPackets]);
                movieDataSliding.movieCFIMatR1   = complex(zeros([size(R0Test,1) size(R0Test,2) maxNumProcessPackets]));
                movieDataSliding.movieCFIMatR0   = zeros([size(R0Test,1) size(R0Test,2) maxNumProcessPackets]);
                movieDataSliding.movieVxMat    = zeros([size(vxCorr,1) size(vxCorr,2) maxNumProcessPackets]);
                movieDataSliding.movieVzMat    = zeros([size(vxCorr,1) size(vzCorr,2) maxNumProcessPackets]);
                movieDataSliding.movieVxMatNT    = zeros([size(vxCorrNT,1) size(vxCorrNT,2) maxNumProcessPackets]);
                movieDataSliding.movieVzMatNT    = zeros([size(vxCorrNT,1) size(vzCorrNT,2) maxNumProcessPackets]);
                movieDataSliding.movieR0hpMat = zeros(size(LowResR0hpMat,1),size(LowResR0hpMat,2),length(rxangles),maxNumProcessPackets);
                movieDataSliding.movieR1Multi = complex(zeros(size(R1Multi,1),size(R1Multi,4),size(R1Multi,3),maxNumProcessPackets));
                movieDataSliding.filterMat = zeros(size(b,1),size(b,2),maxNumProcessPackets);
                movieDataSliding.ROImask = zeros(size(R0Mask,1),size(R0Mask,2),maxNumProcessPackets);

            end
            
          %  movieData = matfile(fullfile(savepath,'movieData.mat'),'Writable',true); % 

            movieDataSliding.movieBmodeMat(:,:,Fctr) = mean(abs(recons_comp.data(:,:,1,:)).^2,4);
            movieDataSliding.movieCFIMatR1(:,:,Fctr) = R1Test;
            movieDataSliding.movieCFIMatR0(:,:,Fctr) = mean(abs(iqhpComp).^2,4);
            movieDataSliding.movieVxMat   (:,:,Fctr) = vxCorr;
            movieDataSliding.movieVzMat   (:,:,Fctr) = vzCorr;
            movieDataSliding.movieVxMatNT (:,:,Fctr) = vxCorrNT;
            movieDataSliding.movieVzMatNT (:,:,Fctr) = vzCorrNT;
            movieDataSliding.movieR0hpMat (:,:,:,Fctr) = LowResR0hpMat;
            movieDataSliding.movieR1MultiMat(:,:,:,Fctr) = R1Multi;
            movieDataSliding.filterMat    (:,:,Fctr) = b;
            movieDataSliding.ROImask (:,:,Fctr) = R0Mask;

            movieDataSliding.movieGeoAxes =  [min(bfPars.x_axis) min(bfPars.z_axis),max(bfPars.x_axis) max(bfPars.z_axis)];
            movieDataSliding.movieVDxAxes = X;
            movieDataSliding.movieVDzAxes = Z;
            movieDataSliding.procParams = procParams;

%             if frame == endFrame
%                 save(fullfile(savepath,'MovieData_FromBFdata_FIR92_4'),'movie*','bfPars','hpPars','-v7.3')
%             end


    end

  %  reset(g); 
    
    disp(sprintf('Finished with %i of %i DopplerFrames',Fctr,nrPackets));

    Fctr = Fctr+1;
end

%% Visualization part
if ~saveForMovie
    % If only one frame.....:
    data = struct();
    R0B = mean(abs(recons_comp.data(:,:,1,:)).^2,4);
    R0BNorm = R0B/max(R0B(:));

    data.bmode = R0BNorm;
    data.R0 = mean(abs(iqhpComp).^2,4);
    data.R1 = R1Test;

    
    bmode = struct();

    bmodedata(:,:,1) = R0BNorm;
    bmodedata(:,:,2) = R0BNorm;

    bmode.data = bmodedata;
    bmode.times = [0 0.025];
    bmode.geo.axes = [min(bfPars.x_axis) min(bfPars.z_axis),max(bfPars.x_axis) max(bfPars.z_axis)];
    bmode.geo.tilt = 0;
    bmode.geo.gridSize = [size(R0BNorm) 1];
    bmode.params.preCompressed = 0;

    R1data(:,:,1) = R1Test;
    R1data(:,:,2) = R1Test;

    R0data(:,:,1) = mean(abs(iqhpComp).^2,4);
    R0data(:,:,2) = mean(abs(iqhpComp).^2,4);

    cfi = struct();


    cfi.data.R0 = R0data.*R0Mask;
    cfi.data.R1 = R1data.*R0Mask;
    cfi.times = [0 0.025];
    cfi.geo.axes = [min(bfPars.x_axis) min(bfPars.z_axis),max(bfPars.x_axis) max(bfPars.z_axis)];
    cfi.geo.tilt = 0;
    cfi.geo.gridSize = [size(R1Test) 1];
    cfi.params.VNyquist = bfPars.PRF*bfPars.c/4/bfPars.f_demod;
    cfi.params.PacketSize = processPacket;
    cfi.params.DynamicRange = 30;

    vxData(:,:,1) = vxCorr;
    vxData(:,:,2) = vxCorr;
    vzData(:,:,1) = vzCorr;
    vzData(:,:,2) = vzCorr;

    vvi.Vx = vxData;
    vvi.Vx(isnan(vvi.Vx))=0;
    vvi.Vz = vzData;
    vvi.Vz(isnan(vvi.Vz))=0;
    vvi.times = [0 0.025];
    vvi.geo.xAxis = unique(X(:));
    vvi.geo.zAxis = unique(Z(:));
    vvi.geo.gridSize = [size(vxCorr) 1];
    vvi.geo.tilt = 0;

    D = struct();
    D.bmode = bmode;
    D.cfi = cfi;
    D.vvi = vvi;

    % Figure setup
    p.scanType='linear';

    p.display.cfi =1;
    p.display.bmode = 1;
    p.display.vvi = 1;
    p.dispVVI = 0;
    p.showEcg = 0;

    p.bmNInterp = 2;
    p.cfNInterp = 2;
    p.cfSmoothR0 = 5;
    p.cfSmoothR1 = 5;

    p.cfColormap = struct('mapType','cfi1');
    p.cfVNyquist = D.cfi.params.VNyquist;

    p.arbToolbar = 0;
    p.bmGain = 0;
    p.bmDyn = 40;
    p.bgColor = [0.9400 0.9400 0.9400];

    % Add speckle tracking data
    % Spatial averaging 
    p.vvSSmoothSpace=[2 2]*1e-3;
    % Vector field visualization setup
    p.vvViz.decX =8; p.vvViz.decZ =4; 
    % Arrows
    p.vvViz.quivVis=1;
    p.vvViz.quivCol=[1 1 1]; % Color
    p.vvViz.quivLineWidth=1;
    p.vvViz.quivScaleArrows=1;
    p.vvViz.quivScaleFactor=1.5;  % Quiver 0.3
    % Streamline
    p.vvViz.streamVis = 0;
    p.vvViz.streamDensity=3;  
    p.vvViz.transStreamVis=0;
    p.vvViz.transStreamColor=[0 0 0];
    p.vvViz.transStreamAlpha=0.2;% Transparent streamline
    % Particle vis
    p.vvViz.particleVis = 1;                            % Turn on/off particle visualization
    p.vvViz.parDeltaT = (D.cfi.times(2)-D.cfi.times(1))/10;                         % Delta time between updates of particles (typically << frametime)
    p.vvViz.parNScat = 750;                             % Number of scatterers in initialized ROI or arbMap
    p.vvViz.parNTail = 2;                               % Number of previous positions in particle tail
    p.vvViz.parAlphaTail = 'exp';
    p.vvViz.parFadeOut = 1;                             % Fade out particles after a given time
    p.vvViz.parDuration = p.vvViz.parDeltaT*4;                         % Duration of particles before fade-out
    p.vvViz.parAlphaInvThres = 0.05;                    % Invalidate particles with alpha value below this threshold
    p.vvViz.parVisColor = [1,1,0];                      % 1 => use parColormap, [1,1,1] (color spec) => specific color
    vMax=1;
    p.vvViz.parVMax = vMax*100;                              % Maximum particle velocity for color scaling
    p.vvViz.parColor = [1,1,0];                         % Color of particles, -1=>use parColormap
    p.vvViz.parColormap = autumn(256);                     % Particle colormap
    p.vvViz.parVelThres = 0.0001*100;                        % Velocity threshold for particle removal, if parRemScat == 1
    p.vvViz.parArbDecRefill = 0.02;
    p.vvViz.parRefillType = 'arbmap';
    % p.vvViz.parLineVelScale = 0.8;% 
    p.vvColormap = struct('mapType','autumn'); 

    p.vvResample = 1;
    p.vvNInterpX=2;
    p.vvArb = 1; 

    p.bmTSmoothTime = 0;%0.01;
    p.cfTSmoothTime = 0;%(D.cfi.times(2)-D.cfi.times(1))*3;
    p.vvTSmoothTime = 0;%;(D.cfi.times(2)-D.cfi.times(1))*3;

    p.batchProcessing = 0;
    p.frameTime = 0;%D.cfi.times - D.cfi.times(1);


    p.vvDefaultAlpha = 1;
    p.vvViz.parAlphaMap=[];
    p.vvVelRange = [0 vMax];%%[-p.cfVNyquist, p.cfVNyquist];


    p.cfCFilterOrder = 2;
    p.cfCFilterType = 'gevu';
    p.leftright = 0;
    p.cfDefaultAlpha = 0.6;

    p.ylimVal = [0.005 0.032];
    p.xlimVal = [-0.018 0.018]; 

    p.cfGain =-3;

    p.cfOverrideArb = 0;          % => 1 to override arbitration map with supplied arbMap
    p.cfArbMap = R0Mask;                   % Arbitration map
    p.vvOverrideArb = 0;          % => 1 to override arbitration map with supplied arbMap
    p.vvArbMap = R0Mask;                   % Arbitration map
    p.vvArb = 1;                  % Arbitrate VVI based on CFI
    pvvTrim = 1;                 % Trim estimates near edges after arbitration
    p.vvSSmooth = 1;              % Spatial smoothing flag


    p.arbPars = struct('arbAlg','gevu','flowPri',0.5,'fPowCut',0.2,'rhoMin',0.01,'NSmoothArb',5,'NSmoothAlpha',5,'NSmoothR0',5,'NSmoothB',5,'arbErodeN',10);

     usFig = USFig(D,p);

else % Save data
    if ~exist(savepath),
         mkdir(savepath);
    end
    movieDataSliding.bfPars = bfPars;
    movieDataSliding.hpPars = hpPars;
    movieDataFileString = sprintf('movieDataSliding%s',hpPars.typeFilter);

    save(fullfile(savepath,movieDataFileString),'-struct','movieDataSliding','-v7.3'); 

  
end 

