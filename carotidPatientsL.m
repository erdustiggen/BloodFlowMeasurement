% Carotid data description file
function data = carotidPatientsL()


    %% 2nd may
    %%{
    bfDepth = 0.032
    ROIdepthStart = 0.005
    ROIdepthEnd = 0.032
    %}
    
    data = {};   
   
 
    %%
     recNo = 1;
    % Data info
    data{ recNo} = struct();
    data{ recNo}.date = '2018.04.04';
    data{ recNo}.description = 'Patient 1, Lars/Håvard, no stenoses, simple patient';
 %   data{ recNo}.data_main_path = ['G:' filesep 'data' filesep '2016.02.03 - Rat 3D - Verasonics' filesep];
    data{ recNo}.data_main_path = ['D:\MasterData\TilEmil\CarotidData\12-Apr-2018\13_34_57\'];
    data{ recNo}.data_path = ['D:\MasterData\TilEmil\CarotidData\12-Apr-2018\13_34_57\'];
    data{ recNo}.data_proc_path = [data{ recNo}.data_path 'bf' filesep];   
    data{ recNo}.chData_suffix = 'Carotid';
    data{ recNo}.iq_file_suffix = 'LowResBFdata_point1';
    data{ recNo}.bmode_file_suffix = 'BMODE';
    data{ recNo}.R0_file_suffix = 'R0';
    data{ recNo}.R1_file_suffix = 'R1';
    data{ recNo}.tdi_file_suffix = 'TDI';
    data{ recNo}.vd_file_suffix = 'VD';    
    data{ recNo}.meta_file = 'Carotid_meta';
    data{ recNo}.bfpars_file = 'bfPars';
    data{ recNo}.point_start = 1;
    data{ recNo}.point_end = 1;
    data{ recNo}.frame_start = 1;
    data{ recNo}.frame_end = 32;
    data{ recNo}.packetsize =270;
    data{ recNo}.ROIdepthStart = ROIdepthStart;
    data{ recNo}.ROIdepthEnd = ROIdepthEnd;
    data{ recNo}.fileType = 'raw';
    data{ recNo}.bfFormat = 'single';
    %data{ recNo}.rcvIx = 1;
    data{ recNo}.samplesPerWave = 4;
    data{ recNo}.glitchframes = [];
    data{ recNo}.rcvBufNum = 2;
    data{ recNo}.intBufNum = 3;
    
     % Processing setup
     data{ recNo}.powerthresh = 15; % Arbitration based on power in the compounded image after clutter filtering
     data{ recNo}.txVals = [-15 -15 -15 -15 -15 15 15 15 15]*pi/180; % NB number 5 and number 9 are only for B-mode images
     data{ recNo}.rxVals = [-15 -3.5 6 15 0 -6 3.5 15 0]*pi/180; % NB number 5 and number 9 are only for B-mode images

    txRx = []; txRx{1} = data{ recNo}.rxVals(1:5); txRx{2} = data{ recNo}.rxVals(6:9);   
    
    data{ recNo}.txinds = [1 3];
    data{ recNo}.reminds = [5 9];
     
    data{ recNo}.RxFnumDop = 1.5; % 
    data{ recNo}.RxFnumBmode = 0.9; % Want to increase resolution in B-mode over Doppler
    data{ recNo}.txRx = txRx; % Beamforming angles
    data{ recNo}.txRxSetName = 'Opt7Max15_2PW';
    data{ recNo}.doAliasingCorrection = 1;

    data{recNo}.demodulationFrequency =4.75e6; % 
    data{recNo}.bandPassVec = [2 2.2 7.3 7.4]*1e6; %
    data{recNo}.bandPassVecPW = [5 5.2 6.6 6.8]*1e6; % !!!! Yet to be decided!!!!
    data{recNo}.PWx = 0.0; % Where to beamform point for PWDoppler, x
    data{recNo}.PWz = 0.011;% Where to beamform point for PWDoppler, z
    
    data{ recNo}.bfDepth = bfDepth;
    
    data{recNo}.concFramesDop =1;
    data{recNo}.concFramesTissue = 1;
    
    data{recNo}.typeFilter = 'FIR138_3';
    data{recNo}.maxClutDimFac = 0.075;
    data{recNo}.maxBloodDimFac = 0.4;
    
   %export image setup
    data{ recNo}.xlimVal = [-2 2]*1e-2;                 
    data{ recNo}.ylimVal = [0.5 2.5]*1e-2;        
    data{ recNo}.parNScat = 500;                             % Number of scatterers in initialized ROI or arbMap
    data{ recNo}.parNTail = 6;                               % Number of previous positions in particle tail
    data{ recNo}.parAlphaTail = 'exp';
    data{ recNo}.parFadeOut = 2;                             % Fade out particles after a given time
    data{ recNo}.parDurationFac = 4;          % Duration of particles before fade-out
    data{ recNo}.vMax = 2.5;
    data{ recNo}.vvSSmoothSpace=[1.5 1.5]*1e-3;
    
    data{ recNo}.savepath = [data{ recNo}.data_proc_path data{recNo}.typeFilter];
    

    
    
    
    
    
    
    
    
    
    