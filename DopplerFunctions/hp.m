function [iqhp,vThresh,b] = hp(iq,v)

[packet,ranges,beams,angles] = size(iq);
v.vNyq = v.PRF*v.c/4/v.f_demod;
switch v.typeFilter
    case 'polyreg'
        % Fm = CoherentCompounding.Doppler.polyregMatrix(v.packetSize,v.polyOrder); % polynomial regression filter
        Fm = polyregMatrix(packet,v.polyOrder); % polynomial regression filter
        iq=reshape(iq,packet,ranges*beams*angles);%reshape to prepare for matrix muliplication
        iqhp=Fm*iq;
        iqhp=reshape(iqhp,packet,ranges,beams,angles);%reshape back to orig. size
        %iqhp = iqhp(2:end-1,:,:); %Avoid edge effects from filter
        [vThresh,b] = getpolyRegCutoff(v); %[m/s]
        
    case 'SVD' % Based on idea by Jerome Barangers et al, paper on adaptive SVD clutter filtering (review)
        M = reshape(iq,[size(iq,1), size(iq,2)*size(iq,3)]); 
        M = M.'; % Follow convention from paper (nx*nx, nt)
        clear iq;
        tic
        [U,S,V] = svd(M,'econ'); % Economy version where only the nececcary column vectors of U is returned (nt)
               
        Q = abs(U);
        meanQ = mean( Q, 1);
        stdQ = std( Q, [], 1);

        detrQ = bsxfun(@minus, Q, meanQ);
        sigmaNM = stdQ.'*stdQ;

        C = detrQ'*detrQ./sigmaNM/size( M,1); % This is the socalled similarity matrix
        
        CDdet_1D = 0; % Method for determining number of clutter space dimensions. Can be based on 1D (fast, less robust) or 2D infor (slower, more robust)     
       
        if CDdet_1D
            % Look at spread of energy along diagonal of matrix
            flipC = fliplr(C);
            moment = zeros( 1, size(M,2)/2);
            for CD = 1:size(M,2)/2 % All possible clutter dimensions
                mydiag = diag( flipC, 2*CD-1);
                ln = length( mydiag);
                halfln = (ln-1)/2;
                arm = (-halfln:halfln).';
                moment(CD) = sum( abs(arm).*mydiag )/sum( abs( arm) );
            end

            moment = fliplr( moment);

            % Find lokal minima. This will give the threshold for clutter space dimension
            % TODO: Make something robust.....
            invmoment = max(moment(:))-moment;
            MPP = 0.005;
            [mypeaks,locs] = findpeaks(invmoment,'MinPeakProminence',MPP);
            [mymaxpeak, locmax] = max(mypeaks);
            ix = locs(locmax); % This is the clutter space dimension threshold

            figure(100)
            plot( invmoment);
            hold on; plot(ix,mymaxpeak,'o')
        
        else
            
            lnC = size( C,1);
            CVec = C(:);
            CVec = CVec-mean(CVec);
            stdC = std(CVec);
            corrVal = zeros( lnC, lnC);
            for kk = 1:floor(lnC*v.maxClutDimFac)
                for mm = kk+1:floor(lnC*v.maxBloodDimFac)
                    compMat = zeros( lnC, lnC);
                    compMat(1:kk, 1:kk) = 1;
                    compMat(kk+1:mm, kk+1:mm) = 1;
                    compVec = compMat(:);
                    corrVal(kk, mm) = mean( CVec.*(compVec-mean(compVec) )/stdC/std(compVec) );
                end
                %clc
                %disp(kk);
            end

            [maxval, ind] = max( corrVal(:) );
            ix = mod( ind, lnC);
            disp([num2str( ix) ' eigenvectors!']);  
            
            if 1
               % Look at spatial singular vectors images
                for kk = 1:15
                    u = reshape(U(:,kk),ranges,beams);
                    figure(10); imagesc(10*log10(abs(u).^2)); title(sprintf('u %i',kk))
                    pause();
                end
 
            end
            
            
        end

        % Hard threshold for filtering out tissue signal
        sigma = diag(S); % Singular values in descending order
        sigma(1:ix) = 0;
      
        Fm = eye(packet).*sigma;
        
        Fiq = U*Fm*V';
      
        iqhp = reshape(Fiq.',packet,ranges,beams); 
        vThresh =0;
        b = ix;
        toc
    case 'eigHans' 
        % Author: Lasse L?vstakken
        % 2017.09.21  Optimalisering til ca 50% kj?retid  Hans Torp
        firstVector = v.firstVector;
        lastVector = v.lastVector;

        IQ = permute(iq,[2,3,1]);
   
        tic;
        ds = 1;
        if ds>1, IQ_eig = IQ(:,:,1:ds:end); else IQ_eig=IQ;end;
        ds_x = 1;
        ds_z = 1;
        if ds_x*ds_z>1,
            IQ_R = IQ_eig(1:ds_x:end, 1:ds_z:end, :);
        else
            IQ_R=IQ_eig;
        end;
        [NR,NB,NP] = size(IQ_R);
        IQ_R = reshape(IQ_R,NR*NB,NP);
        %IQ_R = permute(IQ_R,[2,1])';
        % Estimate correlation matrix.
        R = IQ_R'*IQ_R/(NB*NR);
        % Force Hermetian symmetric (numerical errors).
        R = 0.5*(R'+R);
        [V,D] = eig(R); % LAPACK ZHEEV is used for hermetian symmetric matrices
        % Flip to get largest first.
        V = fliplr(V);
        % Fm = zeros(NP);
        % % 1kHz
        % for n=firstVector:lastVector
        %     Fm = Fm + V(:,n)*V(:,n)'; % Outer product
        % end
        Vused=V(:,firstVector:lastVector);
        Fm=Vused*Vused';
        [NR,NB,NP] = size(IQ_eig);
        % IQ_eig = permute(reshape(IQ_eig,NR*NB,NP),[2,1]);
        % IQ_eig = Fm*IQ_eig;
        % IQF = permute(reshape(IQ_eig,NP,NR,NB),[2,3,1]);
        IQ_eig = reshape(IQ_eig,NR*NB,NP);
        IQ_eig = IQ_eig*Fm;
        IQF = reshape(IQ_eig,NR,NB,NP);
        toc
        if(0)
            D = rot90(D,2);
            lambda = abs(diag(D));
            lambdaLog = 10*log10(lambda);
            % Stem plot
            figure(70);
            clf;
            stem(lambdaLog);   
        end
        
        iqhp =permute(IQF,[3,1,2])
        vThresh = 0;
        b = 0;
        
        
    case  'SVD3D' % Based on idea by Jerome Barangers et al, paper on adaptive SVD clutter filtering (review)
        M = reshape(iq,[size(iq,1), size(iq,2)*size(iq,3)*size(iq,4)]); 
        M = M.'; % Follow convention from paper (nx*ny*nz, nt)
        clear iq;
        tic
        [U,S,V] = svd(M,'econ'); % Economy version where only the nececcary column vectors of U is returned (nt)
               
        Q = abs(U);
        meanQ = mean( Q, 1);
        stdQ = std( Q, [], 1);

        detrQ = bsxfun(@minus, Q, meanQ);
        sigmaNM = stdQ.'*stdQ;

        C = detrQ'*detrQ./sigmaNM/size( M,1); % This is the socalled similarity matrix
        
        CDdet_1D = v.CDdet_1D; % Method for determining number of clutter space dimensions. Can be based on 1D (fast, less robust) or 2D infor (slower, more robust)     
        estClutDim = v.estClutDim;
        estBloodDim = v.estBloodDim;
        bloodFraction = v.bloodFraction; % For testing purposes
         
        if CDdet_1D
              
            % Look at spread of energy along diagonal of matrix
            flipC = fliplr(C);
            moment = zeros( 1, size(M,2)/2);
            for CD = 1:size(M,2)/2 % All possible clutter dimensions
                mydiag = diag( flipC, 2*CD-1);
                ln = length( mydiag);
                halfln = (ln-1)/2;
                arm = (-halfln:halfln).';
                moment(CD) = sum( abs(arm).*mydiag )/sum( abs( arm) );
            end

            moment = fliplr( moment);

            
            if estClutDim
                % Find lokal minima. This will give the threshold for clutter space dimension
                invmoment = smooth(max(moment(:))-moment,7,'moving');
                MPP = 0.01;% 0.005;
                [mypeaks,locs] = findpeaks(invmoment,'MinPeakProminence',MPP);
                
                if isempty(mypeaks)
                     MPP = 0.005;% 0.005;
                     [mypeaks,locs] = findpeaks(invmoment,'MinPeakProminence',MPP);
                end
                
                if isempty(mypeaks)
                     MPP = 0.001;% 0.005;
                     [mypeaks,locs] = findpeaks(invmoment,'MinPeakProminence',MPP);
                end
                
             
                
                
                if 0
                    [mymaxpeak, locmax] = max(mypeaks);
                    ix = locs(locmax); % This is the clutter space dimension threshold
                    disp([num2str( ix) ' eigenvectors!']);  
                else
                    
                       if isempty(mypeaks)
                           ix = v.maxClutDimFac*packet;
                            sprintf('No peak found. Using max allowed clutter dim, at eig # %i',ix)

                       else
                            mymaxpeak = mypeaks(1);
                            ix = locs(1);
                            sprintf('Using first peak. Clutter peak at eig # %i',ix)
                       end
                end

                figure(100)
                plot( invmoment);
                hold on; plot(ix,mymaxpeak,'o'); hold off
                b = ix;
            elseif estBloodDim 
                MPP = 0.005;
                sMoment = smooth(moment,10,'moving');
                [mypeaks,locs] = findpeaks(sMoment,'MinPeakProminence',MPP);
                
                if 0
                    [mymaxpeak, locmax] = max(mypeaks);
                    bPeak = locs(locmax); % This should correspond to the eigenvector with maximum blood signal
                    disp('Using max peak')
                else
                    mymaxpeak = mypeaks(end);
                    bPeak = locs(end);
                    sprintf('Using last peak. Blood peak at eig # %i',bPeak)
                end
                figure(100)
                subplot(1,6,1);
                plot(moment,'r','LineWidth',2);
                hold on; plot(sMoment,'k','LineWidth',2);
                hold on; plot(bPeak,mymaxpeak,'og','LineWidth',2)
                hold off;
                
                bEigs = (bPeak-floor(bloodFraction*size(M,2)/2):bPeak+floor(bloodFraction*size(M,2)/2));
                ix = 1:size(M,2);
                ix(bEigs) = [];
                b = bPeak;
            end
        
        else
            
            lnC = size( C,1);
            Cd = C - eye(lnC);
            CVec = Cd(:);
            CVec = CVec-mean(CVec);
            stdC = std(CVec);
            corrVal = zeros( lnC, lnC);
            
             for kk = 1:floor(lnC*v.maxClutDimFac)
                for mm = kk+1:floor(lnC*v.maxBloodDimFac)
                    compMat = zeros( lnC, lnC);
                    compMat(1:kk, 1:kk) = 1;
                    compMat(kk+1:mm, kk+1:mm) = 1;
                    compMat(logical(eye(size(compMat)))) = 0;
                    compVec = compMat(:);
                    corrVal(kk, mm) = mean( CVec.*(compVec-mean(compVec) )/stdC/std(compVec) );
%                     corrVal(kk, mm) = mean( CVec(compVec == 1) );
                end
         
            end

            [maxval, ind] = max( corrVal(:) );
            ix = mod( ind, lnC);
            %ix = ix+2
            %disp([num2str( ix) ' eigenvectors! NB, have added 2 to threhold']);  
            disp([num2str( ix) ' eigenvectors!']);  

        end
        
        if 0          % Look at spatial singular vectors images
        for kk = 1:packet
            u = reshape(U(:,kk),ranges,beams,angles);
            figure(100); subplot(1,6,2); 
            imagesc(10*log10(interp2(sum(abs(u).^2,3),3))); title(sprintf('u %i',kk)); caxis([-50 -20]); 

            pause();
        end
        end

        % Hard threshold for filtering out tissue signal
        sigma = diag(S); % Singular values in descending order
        
        if estClutDim
            sigma(1:ix) = 0;
        elseif estBloodDim
            sigma(ix) = 0;
        end
        
        % Soft threshold, testing
        if isfield(v,'soft')
            if v.soft
                mywin = tukeywin(packet*2,0.6);
                sigma(ix+1:end) = sigma(ix+1:end).*mywin(packet+1:end-ix);
            end
        end
      
        Fm = eye(packet).*sigma;
        
        Fiq = U*Fm*V';
      
        iqhp = reshape(Fiq.',packet,ranges,beams,angles); 
        vThresh =0;
        toc
        
        
        if CDdet_1D && estBloodDim
           
            R0hp = iqhp.*conj(iqhp);
            figure(100);
            subplot(1,6,3);
            imagesc(10*log10(interp2(squeeze(mean(mean(R0hp,1),4)),3))); 
            title(sprintf('PD,eigs %i to %i',min(bEigs),max(bEigs)));
            colorbar;
            
        end
        
        
    case 'whitening'
                
        % 2017.10.07 Author: LHans Torp
        ds_x = 1;
        ds_z = 1;
        ds_y = 1;
        if ds_x*ds_z*ds_y>1,
            IQ_R = iq(1:ds_z:end,1:ds_x:end,1:ds_y:end, :);
        else
            IQ_R=iq;
        end;
        IQ_R=double(IQ_R);
        [NR,NB,NY,NP] = size(IQ_R);
        IQ_R = reshape(IQ_R,NR*NB*NY,NP);
        % Estimate correlation matrix.
        R = IQ_R'*IQ_R/(NB*NR*NY);
        % Force Hermetian symmetric (numerical errors).
        R = 0.5*(R'+R);    
        R=sqrtm(R);
        Fm=inv(R);  % whitening filter if we consider sqrtm
        [NR,NB,NY,NP] = size(iq);
        iq= reshape(iq,NR*NB*NY,NP);
        iq=iq*Fm;
        iqhp = reshape(iq,NR,NB,NY,NP);
        vThresh =0;
        b = [];
        
        
        
    case 'whitened_eig'
        
        firstVector = v.firstVector;
        lastVector = v.lastVector;
        
        ds_x = 1;
        ds_z = 1;
        ds_y = 1;
        if ds_x*ds_z*ds_y>1
            IQ_R = iq(1:ds_z:end,1:ds_x:end,1:ds_y:end, :);
        else
            IQ_R=iq;
        end
        IQ_R=double(IQ_R);
        [NR,NB,NY,NP] = size(IQ_R);
        IQ_R = reshape(IQ_R,NR*NB*NY,NP);
        
        [NR,NB,NY,NP] = size(iq);
        iq= reshape(iq,NR*NB*NY,NP);
        
        % Estimate correlation matrix.
        R = IQ_R'*IQ_R/(NB*NR*NY);
        % Force Hermetian symmetric (numerical errors).
        R = 0.5*(R'+R);    
        Rsqrt=sqrtm(R);
       
        Fm=inv(Rsqrt);  % whitening filter if we consider sqrtm
        iq_white=iq*Fm;
        
        % Eigenvector decomposition of correlation matrix
        [V,D] = eig(R);
        V = fliplr(V);
        Vused=V(:,firstVector:lastVector);
        Fm=Vused*Vused';
        iq = iq_white*Fm;
        
        iqhp = reshape(iq,NR,NB,NY,NP);
        vThresh =0;
        b = [];
        
        
     case  'SVD_white_test' % 
            
        clutterComps = v.clutterDims;

        M = reshape(iq,[size(iq,1), size(iq,2)*size(iq,3)*size(iq,4)]); 
        M = M.'; % Follow convention from paper (nx*ny*nz, nt)
        clear iq;
      
        [U,S,V] = svd(M,'econ'); % Economy version where only the nececcary column vectors of U is returned (nt)
        
        sigma = ones(length(diag(S)),1);
        sigma(1:clutterComps) = 0;
           
        Fm = eye(packet).*sigma;
        
        Fiq = U*Fm*V';
      
        iqhp = reshape(Fiq.',packet,ranges,beams,angles); 
        vThresh =0;
        b = [];   
        
    case 'eig'
        clutterComps = v.clutterDims;
        
        IQ_R = reshape(iq,packet,ranges*beams)'; % ! Something fishy here
        IQ_F = reshape(iq,packet,ranges*beams);
        R = IQ_R'*IQ_R/(beams*ranges);     
        % Force Hermetian symmetric (numerical errors).
        R = 0.5*(R'+R);            
        [V,D] = eig(R); 
        % Get largest eigenvectors first
        [Evals,Index] = sort(diag(D),'descend');
        V         = V(:,Index);       
%         D = fliplr(flipud(D));
%         V = fliplr(V);
        Fm = zeros(packet);   
    
        for n=clutterComps+1:packet,
            Fm = Fm + V(:,n)*V(:,n)'; % Outer product
        end        
    
        iqhp = Fm*IQ_F;
        iqhp = reshape(iqhp,packet,ranges,beams);
        vThresh = 0;
        b = 0;
        
        
    case 'hbfi'
        b = [-0.4909,0.6436,0.0773,-0.1226,-0.1077]; %regular bfi filter
        iqhp = filter(b,1,iq,[],1);
        iqhp = iqhp(length(b):end,:,:);
        vThresh = v.c*658.7/(2*v.f_demod); % -3dB cutoff [m/s]
        
    case 'hSTLow1'
        b =[-0.0860,-0.1868,-0.6281,0.6281,0.1868, 0.0860]; % Very low cut-off, zero-point at DC
        iqhp = filter(b,1,iq,[],1);
        iqhp = iqhp(length(b):end,:,:);
        vThresh = 0.04728; % -3dB cutoff [m/s]
        
    case 'hSTLow2'
        b = [-0.0995   -0.1080   -0.1144   -0.1184  0.8803   -0.1184   -0.1144   -0.1080   -0.0995]; % Cut-off 0.2 norm.,70dB damping
        iqhp = filter(b,1,iq,[],1);
        if angles == 1
             iqhp = iqhp(length(b):end,:,:);
        elseif angles>1
             iqhp = iqhp(length(b):end,:,:,:);
        end
        vThresh = 0.0539; % -3dB cutoff [m/s]
        
    case 'FIR92'

   b = [ 0.0118   -0.0029   -0.0028   -0.0028   -0.0028   -0.0029   -0.0030   -0.0030   -0.0029   -0.0027   -0.0024   -0.0019   -0.0013   -0.0005    0.0004    0.0015    0.0027    0.0039    0.0051...
         0.0063    0.0073    0.0082    0.0088    0.0092    0.0091    0.0086    0.0076    0.0061    0.0041    0.0016   -0.0015   -0.0051   -0.0091   -0.0135   -0.0182   -0.0231   -0.0281   -0.0331...
        -0.0380   -0.0427   -0.0470   -0.0509   -0.0542   -0.0568   -0.0588   -0.0600    0.9396   -0.0600   -0.0588   -0.0568   -0.0542   -0.0509   -0.0470   -0.0427   -0.0380   -0.0331   -0.0281...
        -0.0231   -0.0182   -0.0135   -0.0091   -0.0051   -0.0015    0.0016    0.0041    0.0061    0.0076    0.0086    0.0091    0.0092    0.0088    0.0082    0.0073    0.0063    0.0051    0.0039...
         0.0027    0.0015    0.0004   -0.0005   -0.0013   -0.0019   -0.0024   -0.0027   -0.0029   -0.0030   -0.0030   -0.0029   -0.0028   -0.0028   -0.0028   -0.0029    0.0118];
        
        iqhp = filter(b,1,iq,[],1);
        iqhp = iqhp(length(b):end,:,:);
        vThresh =  0.068*v.vNyq; % -3dB cutoff [m/s]
     
        
    case 'FIR92_2'
        
    b = [0.0054   -0.0056   -0.0038   -0.0028   -0.0022   -0.0018   -0.0015   -0.0012   -0.0008   -0.0003    0.0003    0.0010    0.0017    0.0025    0.0034    0.0042    0.0050    0.0058    0.0064...
         0.0069    0.0072    0.0073    0.0071    0.0066    0.0058    0.0046    0.0030    0.0010   -0.0013   -0.0040   -0.0070   -0.0103   -0.0139   -0.0176   -0.0215   -0.0255   -0.0295   -0.0334...
        -0.0371   -0.0407   -0.0439   -0.0467   -0.0492   -0.0511   -0.0525   -0.0534    0.9463   -0.0534   -0.0525   -0.0511   -0.0492   -0.0467   -0.0439   -0.0407   -0.0371   -0.0334   -0.0295...
        -0.0255   -0.0215   -0.0176   -0.0139   -0.0103   -0.0070   -0.0040   -0.0013    0.0010    0.0030    0.0046    0.0058    0.0066    0.0071    0.0073    0.0072    0.0069    0.0064    0.0058...
         0.0050    0.0042    0.0034    0.0025    0.0017    0.0010    0.0003   -0.0003   -0.0008   -0.0012   -0.0015   -0.0018   -0.0022   -0.0028   -0.0038   -0.0056    0.0054];

        iqhp = filter(b,1,iq,[],1);
        iqhp = iqhp(length(b):end,:,:);
        vThresh =  0.0613*v.vNyq; % -3dB cutoff [m/s]
        
        
    case 'FIR92AutoGen'
        b = createFIR92filter(v.Cutoff);
        iqhp = filter(b,1,iq,[],1);
        iqhp = iqhp(length(b):end,:,:);
        vThresh =  v.Cutoff*v.vNyq; % -3dB cutoff [m/s]
        
    case 'FIR92_5'
        
    b = [-0.01195   0.0077701   0.0056472   0.0040636   0.0028319   0.0017937 ...
        0.0008621 -3.4136e-05 -0.00092385  -0.0018065  -0.0026572  -0.0034437 ...
        -0.0041062  -0.0045899  -0.0048323  -0.0047829  -0.0043936  -0.0036403 ...
        -0.0025128  -0.0010364  0.00074624   0.0027567   0.0048978   0.0070417 ...
        0.0090532    0.010777    0.012061    0.012747    0.012697    0.011785 ...
        0.0099284   0.0070714   0.0032086  -0.0016335   -0.007376   -0.013897 ...
        -0.021018   -0.028539   -0.036227    -0.04384   -0.051115   -0.057787 ...
        -0.063615    -0.06841   -0.071961   -0.074135     0.92511   -0.074135 ...
        -0.071961    -0.06841   -0.063615   -0.057787   -0.051115    -0.04384 ...
        -0.036227   -0.028539   -0.021018   -0.013897   -0.007376  -0.0016335 ...
        0.0032086   0.0070714   0.0099284    0.011785    0.012697    0.012747 ...
        0.012061    0.010777   0.0090532   0.0070417   0.0048978   0.0027567  ...
        0.00074624  -0.0010364  -0.0025128  -0.0036403  -0.0043936  -0.0047829  ...
        -0.0048323  -0.0045899  -0.0041062  -0.0034437  -0.0026572  -0.0018065 ...
        -0.00092385 -3.4136e-05   0.0008621   0.0017937   0.0028319   0.0040636  ...
        0.0056472   0.0077701    -0.01195];
        
    iqhp = filter(b,1,iq,[],1);
    iqhp = iqhp(length(b):end,:,:);
    vThresh =  0.074*v.vNyq; % -6dB cutoff [m/s]
        
    case 'FIR96'
  b = [  -0.0068    0.0021    0.0020    0.0019    0.0020    0.0020    0.0019    0.0017    0.0014    0.0010    0.0003   -0.0004   -0.0012   -0.0021   -0.0030   -0.0038   -0.0044   -0.0048 ...
        -0.0049   -0.0046   -0.0038   -0.0027   -0.0012    0.0007    0.0029    0.0052    0.0075    0.0096    0.0114    0.0127    0.0133    0.0130    0.0117    0.0093    0.0058    0.0011 ...
        -0.0046   -0.0113   -0.0187   -0.0267   -0.0350   -0.0434   -0.0514   -0.0588   -0.0653   -0.0707   -0.0747   -0.0772    0.9220   -0.0772   -0.0747   -0.0707   -0.0653   -0.0588 ...
        -0.0514   -0.0434   -0.0350   -0.0267   -0.0187   -0.0113   -0.0046    0.0011    0.0058    0.0093    0.0117    0.0130    0.0133    0.0127    0.0114    0.0096    0.0075    0.0052 ...
         0.0029    0.0007   -0.0012   -0.0027   -0.0038   -0.0046   -0.0049   -0.0048   -0.0044   -0.0038   -0.0030   -0.0021   -0.0012   -0.0004    0.0003    0.0010    0.0014    0.0017 ...
        0.0019    0.0020    0.0020    0.0019    0.0020    0.0021   -0.0068];
        
        iqhp = filter(b,1,iq,[],1);
        iqhp = iqhp(length(b):end,:,:);
        vThresh =  0.0859*v.vNyq; % -3dB cutoff [m/s] (wstop = 0.03*vNyq,wpass = 0.1*vNyq, 100dB stopband)
        
    case 'FIR134'
    b = [0.0067   -0.0100   -0.0030    0.0007    0.0023    0.0027    0.0026    0.0021    0.0014    0.0007    0.0000   -0.0006   -0.0012   -0.0017   -0.0020   -0.0023   -0.0023   -0.0022 ...
        -0.0018   -0.0013   -0.0006    0.0002    0.0010    0.0019    0.0028    0.0034    0.0039    0.0041    0.0040    0.0035    0.0027    0.0016    0.0002   -0.0014   -0.0031   -0.0047 ...
       -0.0061   -0.0072   -0.0078   -0.0079   -0.0073   -0.0061   -0.0042   -0.0018    0.0011    0.0043    0.0076    0.0107    0.0134    0.0154    0.0165    0.0164    0.0150    0.0121 ...
        0.0077    0.0019   -0.0053   -0.0136   -0.0228   -0.0326   -0.0425   -0.0523   -0.0613   -0.0694   -0.0761   -0.0811   -0.0842    0.9148   -0.0842   -0.0811   -0.0761   -0.0694 ...
       -0.0613   -0.0523   -0.0425   -0.0326   -0.0228   -0.0136   -0.0053    0.0019    0.0077    0.0121    0.0150    0.0164    0.0165    0.0154    0.0134    0.0107    0.0076    0.0043 ...
        0.0011   -0.0018   -0.0042   -0.0061   -0.0073   -0.0079   -0.0078   -0.0072   -0.0061   -0.0047   -0.0031   -0.0014    0.0002    0.0016    0.0027    0.0035    0.0040    0.0041 ...
        0.0039    0.0034    0.0028    0.0019    0.0010    0.0002   -0.0006   -0.0013   -0.0018   -0.0022   -0.0023   -0.0023   -0.0020   -0.0017   -0.0012   -0.0006    0.0000    0.0007 ...
        0.0014    0.0021    0.0026    0.0027    0.0023    0.0007   -0.0030   -0.0100    0.0067 ];

        iqhp = filter(b,1,iq,[],1);
        iqhp = iqhp(length(b):end,:,:);
        vThresh =  0.0912*v.vNyq; % -3dB cutoff [m/s] (wstop = 0.05*vNyq,wpass = 0.1*vNyq)
        
    case 'FIR98'
        
    b = [0.0114   -0.0093   -0.0056   -0.0029   -0.0007    0.0010    0.0022    0.0031    0.0035    0.0032    0.0024    0.0011   -0.0006   -0.0023   -0.0039   -0.0049   -0.0051   -0.0044 ...
        -0.0027   -0.0002    0.0026    0.0054    0.0075    0.0084    0.0079    0.0057    0.0022   -0.0023   -0.0070   -0.0111   -0.0137   -0.0141   -0.0118   -0.0069    0.0002    0.0085 ...
         0.0167    0.0234    0.0269    0.0260    0.0197    0.0077   -0.0095   -0.0308   -0.0546   -0.0785   -0.1003   -0.1178   -0.1291    0.8670   -0.1291   -0.1178   -0.1003   -0.0785 ... 
        -0.0546   -0.0308   -0.0095    0.0077    0.0197    0.0260    0.0269    0.0234    0.0167    0.0085    0.0002   -0.0069   -0.0118   -0.0141   -0.0137   -0.0111   -0.0070   -0.0023 ...
         0.0022    0.0057    0.0079    0.0084    0.0075    0.0054    0.0026   -0.0002   -0.0027   -0.0044   -0.0051   -0.0049   -0.0039   -0.0023   -0.0006    0.0011    0.0024    0.0032 ...
         0.0035    0.0031    0.0022    0.0010   -0.0007   -0.0029   -0.0056   -0.0093    0.0114];
        
     
        iqhp = filter(b,1,iq,[],1);
        iqhp = iqhp(length(b):end,:,:);
        vThresh =  0.14*v.vNyq; % -3dB cutoff [m/s] (wstop = 0.1*vNyq,wpass = 0.15*vNyq, 70 dB stopband)
        
        
        
    case 'FIR88'
        b = [0.0024   -0.0104    0.0042    0.0052    0.0035    0.0016   -0.0002   -0.0016   -0.0027   -0.0034   -0.0036   -0.0031   -0.0020   -0.0004    0.0016    0.0036    0.0053    0.0063 ...
             0.0062    0.0049    0.0025   -0.0008   -0.0046   -0.0082   -0.0109   -0.0121   -0.0112   -0.0081   -0.0029    0.0039    0.0113    0.0182    0.0232    0.0251    0.0227    0.0155 ...
             0.0032   -0.0137   -0.0341   -0.0563   -0.0784   -0.0985   -0.1144   -0.1247    0.8718   -0.1247   -0.1144   -0.0985   -0.0784   -0.0563   -0.0341   -0.0137    0.0032    0.0155 ...
             0.0227    0.0251    0.0232    0.0182    0.0113    0.0039   -0.0029   -0.0081   -0.0112   -0.0121   -0.0109   -0.0082   -0.0046   -0.0008    0.0025    0.0049    0.0062    0.0063 ...
             0.0053    0.0036    0.0016   -0.0004   -0.0020   -0.0031   -0.0036   -0.0034   -0.0027   -0.0016   -0.0002    0.0016    0.0035    0.0052    0.0042   -0.0104    0.0024];

        iqhp = filter(b,1,iq,[],1);
        iqhp = iqhp(length(b):end,:,:);
        vThresh =  0.137*v.vNyq; % -3dB cutoff [m/s] (wstop = 0.08*vNyq,wpass = 0.15*vNyq)
        
        
    case 'FIR70'

   b = [ 0.0050   -0.0059   -0.0040   -0.0031   -0.0025   -0.0020   -0.0014   -0.0006    0.0004    0.0016    0.0030    0.0045    0.0060    0.0074    0.0086    0.0094    0.0098    0.0095    0.0085...
         0.0067    0.0040    0.0004   -0.0041   -0.0093   -0.0153   -0.0218   -0.0287   -0.0358   -0.0428   -0.0495   -0.0557   -0.0611   -0.0655   -0.0688   -0.0709    0.9285   -0.0709   -0.0688...
        -0.0655   -0.0611   -0.0557   -0.0495   -0.0428   -0.0358   -0.0287   -0.0218   -0.0153   -0.0093   -0.0041    0.0004    0.0040    0.0067    0.0085    0.0095    0.0098    0.0094    0.0086...
         0.0074    0.0060    0.0045    0.0030    0.0016    0.0004   -0.0006   -0.0014   -0.0020   -0.0025   -0.0031   -0.0040   -0.0059    0.0050];
        
        iqhp = filter(b,1,iq,[],1);
        iqhp = iqhp(length(b):end,:,:);
        vThresh =  0.0817*v.vNyq; % -3dB cutoff [m/s] (wstop = 0.02*vNyq,wpass = 0.1*vNyq)
        
    case 'FIR74'
        
    b = [0.0112   -0.0033   -0.0032   -0.0033   -0.0035   -0.0036   -0.0037   -0.0035   -0.0031   -0.0024   -0.0014   -0.0000    0.0016    0.0034    0.0053    0.0072    0.0089    0.0103    0.0112    0.0115    0.0109    0.0095    0.0070    0.0036...
        -0.0009   -0.0064   -0.0128   -0.0198   -0.0274   -0.0352   -0.0430   -0.0505   -0.0575   -0.0636   -0.0686   -0.0723   -0.0747    0.9246   -0.0747   -0.0723   -0.0686   -0.0636   -0.0575   -0.0505   -0.0430   -0.0352   -0.0274   -0.0198...
        -0.0128   -0.0064   -0.0009    0.0036    0.0070    0.0095    0.0109    0.0115    0.0112    0.0103    0.0089    0.0072    0.0053    0.0034    0.0016   -0.0000   -0.0014   -0.0024   -0.0031   -0.0035   -0.0037   -0.0036   -0.0035   -0.0033...
        -0.0032   -0.0033    0.0112];
    
        iqhp = filter(b,1,iq,[],1);
        iqhp = iqhp(length(b):end,:,:);
        
        vThresh = 0.0851*v.vNyq; %
        
        
        case 'FIR30'
     
    b = [-0.0118   -0.0031   -0.0000    0.0054    0.0121    0.0181    0.0209    0.0179    0.0071   -0.0119   -0.0382   -0.0690   -0.1001   -0.1269   -0.1450    0.8486   -0.1450   -0.1269   -0.1001   -0.0690   -0.0382   -0.0119    0.0071    0.0179...    
         0.0209    0.0181    0.0121    0.0054   -0.0000   -0.0031   -0.0118];
     
       
        iqhp = filter(b,1,iq,[],1);
        iqhp = iqhp(length(b):end,:,:);
        
        vThresh = 0.1726*v.vNyq; %( wstop = 0.05*vNyq,wpass = 0.21*vNyq)
        
        
    case 'FIR30_2'
    
    b = [-0.0054    0.0083    0.0091    0.0109    0.0117    0.0101    0.0051   -0.0040   -0.0171   -0.0338   -0.0528   -0.0723   -0.0905   -0.1052   -0.1149    0.8818   -0.1149   -0.1052   -0.0905   -0.0723   -0.0528   -0.0338   -0.0171   -0.0040...
        0.0051    0.0101    0.0117    0.0109    0.0091    0.0083   -0.0054];
          
        iqhp = filter(b,1,iq,[],1);
        iqhp = iqhp(length(b):end,:,:);
        
        vThresh = 0.14*v.vNyq; %( wstop = 0.02*vNyq,wpass = 0.18*vNyq)
     
        

    case 'FIR92_3'
    
    b =  [0.0118   -0.0029   -0.0028   -0.0028   -0.0028   -0.0029   -0.0030   -0.0030   -0.0029   -0.0027   -0.0024   -0.0019   -0.0013   -0.0005    0.0004    0.0015    0.0027 ...
    0.0039    0.0051    0.0063    0.0073    0.0082    0.0088    0.0092    0.0091    0.0086    0.0076    0.0061    0.0041    0.0016   -0.0015   -0.0051   -0.0091   -0.0135 ...
   -0.0182   -0.0231   -0.0281   -0.0331   -0.0380   -0.0427   -0.0470   -0.0509   -0.0542   -0.0568   -0.0588   -0.0600    0.9396   -0.0600   -0.0588   -0.0568   -0.0542 ...
   -0.0509   -0.0470   -0.0427   -0.0380   -0.0331   -0.0281   -0.0231   -0.0182   -0.0135   -0.0091   -0.0051   -0.0015    0.0016    0.0041    0.0061    0.0076    0.0086 ...
    0.0091    0.0092    0.0088    0.0082    0.0073    0.0063    0.0051    0.0039    0.0027    0.0015    0.0004   -0.0005   -0.0013   -0.0019   -0.0024   -0.0027   -0.0029 ...
   -0.0030   -0.0030   -0.0029   -0.0028   -0.0028   -0.0028   -0.0029    0.0118];


    iqhp = filter(b,1,iq,[],1);
    iqhp = iqhp(length(b):end,:,:);

    vThresh = 0.08*v.vNyq; %wstop = 0.02*vNyq; wpass = 0.08*vNyq
        
    
    case 'FIR92_4'
        
    b = [-0.0083    0.0039    0.0081    0.0004    0.0038   -0.0019    0.0007   -0.0037   -0.0017   -0.0046   -0.0025   -0.0037   -0.0010   -0.0008    0.0022    0.0030    0.0055 ...
        0.0057    0.0065    0.0051    0.0037    0.0005   -0.0024   -0.0062   -0.0089   -0.0112   -0.0115   -0.0105   -0.0072   -0.0025    0.0037    0.0102    0.0165    0.0213 ...
        0.0237    0.0227    0.0178    0.0087   -0.0044   -0.0209   -0.0396   -0.0594   -0.0786   -0.0956   -0.1089   -0.1174    0.8797   -0.1174   -0.1089   -0.0956   -0.0786 ...
        -0.0594   -0.0396   -0.0209   -0.0044    0.0087    0.0178    0.0227    0.0237    0.0213    0.0165    0.0102    0.0037   -0.0025   -0.0072   -0.0105   -0.0115   -0.0112 ...
        -0.0089   -0.0062   -0.0024    0.0005    0.0037    0.0051    0.0065    0.0057    0.0055    0.0030    0.0022   -0.0008   -0.0010   -0.0037   -0.0025   -0.0046   -0.0017 ...
        -0.0037    0.0007   -0.0019    0.0038    0.0004    0.0081    0.0039   -0.0083];
    
    
    iqhp = filter(b,1,iq,[],1);
    iqhp = iqhp(length(b):end,:,:);

    vThresh = 0.14*v.vNyq; %wstop = 0.08*vNyq; wpass = 0.14*vNyq
    
    case 'FIR92_6' %matches FIR138_3, wstop 0.045, wpass 0.08, Astop 69
        b = [-0.0142217288112245       0.00598581498191368       0.00504012241262655       0.00432719500222103       0.00373449194735445...
        0.00317160097953801       0.00256865363050219       0.00186370202528736       0.00103853189256352       8.5252173198542e-05...
     -0.000971002024746796      -0.00209097877421146      -0.00320834233235719      -0.00424217279819033      -0.00510303799105878...
       -0.0057014764247235      -0.00595468522558224      -0.00579182112962306      -0.00516281005188153       -0.0040474036619614...
      -0.00245850995408867      -0.00044708287530096       0.00190053731598859       0.00446023954852331        0.0070767622376189...
       0.00956913404482788        0.0117405794475341        0.0133877276561731        0.0143105532748599        0.0143260250120001...
         0.013281953734215        0.0110685685028143       0.00762842881726474       0.00295998634460994      -0.00287493712776024...
      -0.00975314549970265       -0.0174931810812164        -0.025861445566922       -0.0345849171098829       -0.0433466446993566...
       -0.0518196342425656       -0.0596632230767555       -0.0665803289144022       -0.0722877690466454        -0.076547281971846...
       -0.0791846779673601         0.919960316981144       -0.0791846779673601        -0.076547281971846       -0.0722877690466454...
       -0.0665803289144022       -0.0596632230767555       -0.0518196342425656       -0.0433466446993566       -0.0345849171098829...
        -0.025861445566922       -0.0174931810812164      -0.00975314549970265      -0.00287493712776024       0.00295998634460994...
       0.00762842881726474        0.0110685685028143         0.013281953734215        0.0143260250120001        0.0143105532748599...
        0.0133877276561731        0.0117405794475341       0.00956913404482788        0.0070767622376189       0.00446023954852331...
       0.00190053731598859      -0.00044708287530096      -0.00245850995408867       -0.0040474036619614      -0.00516281005188153...
      -0.00579182112962306      -0.00595468522558224       -0.0057014764247235      -0.00510303799105878      -0.00424217279819033...
      -0.00320834233235719      -0.00209097877421146     -0.000971002024746796       8.5252173198542e-05       0.00103853189256352...
       0.00186370202528736       0.00256865363050219       0.00317160097953801       0.00373449194735445       0.00432719500222103...
       0.00504012241262655       0.00598581498191368       -0.0142217288112245];
    
        iqhp = filter(b,1,iq,[],1);
        iqhp = iqhp(length(b):end,:,:);

        vThresh = 0.098*v.vNyq; %????????
        
    
    case 'FIR138' % Should be used with two transmit angles, similar observation time as when used with three angles and FIR92_4
        b = [-0.0163    0.0106    0.0073    0.0048    0.0028    0.0013    0.0001   -0.0008   -0.0016   -0.0022   -0.0026   -0.0028   -0.0029   -0.0029   -0.0027   -0.0023   -0.0017   -0.0010   -0.0002    0.0007    0.0016 ...
            0.0024    0.0032    0.0039    0.0043    0.0045    0.0044    0.0041    0.0034    0.0024    0.0012   -0.0002   -0.0017   -0.0033   -0.0048   -0.0061   -0.0071   -0.0078   -0.0080   -0.0077   -0.0069   -0.0055 ...
           -0.0036   -0.0013    0.0013    0.0042    0.0070    0.0098    0.0122    0.0141    0.0154    0.0157    0.0151    0.0134    0.0106    0.0065    0.0014   -0.0048   -0.0118   -0.0195   -0.0277   -0.0361   -0.0444 ...
           -0.0524   -0.0597   -0.0660   -0.0713   -0.0751   -0.0775    0.9217   -0.0775   -0.0751   -0.0713   -0.0660   -0.0597   -0.0524   -0.0444   -0.0361   -0.0277   -0.0195   -0.0118   -0.0048    0.0014    0.0065 ...
            0.0106    0.0134    0.0151    0.0157    0.0154    0.0141    0.0122    0.0098    0.0070    0.0042    0.0013   -0.0013   -0.0036   -0.0055   -0.0069   -0.0077   -0.0080   -0.0078   -0.0071   -0.0061   -0.0048 ...
           -0.0033   -0.0017   -0.0002    0.0012    0.0024    0.0034    0.0041    0.0044    0.0045    0.0043    0.0039    0.0032    0.0024    0.0016    0.0007   -0.0002   -0.0010   -0.0017   -0.0023   -0.0027   -0.0029 ...
           -0.0029   -0.0028   -0.0026   -0.0022   -0.0016   -0.0008    0.0001    0.0013    0.0028    0.0048    0.0073    0.0106   -0.0163 ];
    
       
        iqhp = filter(b,1,iq,[],1);
        iqhp = iqhp(length(b):end,:,:);

        vThresh = 0.14*v.vNyq; %????????
        
    case 'FIR138_2' % Wider stopband (wstop = 0.1 Nyquist ) than FIR138
        b = [-0.008230911208292   0.013677737052159   0.002315432349574  -0.002772127735242  -0.004190745808755  -0.003656865797247  -0.002258584801907  -0.000656289897754   0.000773873754742   0.001833779324514   0.002434430195933 ...
   0.002558105616276   0.002239935097602   0.001547735038285   0.000592373907412  -0.000485908218947  -0.001532267620955  -0.002382210801555  -0.002883343413527  -0.002932256415101  -0.002480137602476  -0.001559677092757 ...
  -0.000290114823546   0.001153092856976   0.002528473584716   0.003603438954240   0.004158015837248   0.004051381673564   0.003231644105433   0.001768436318758  -0.000153554698261  -0.002249136008113  -0.004172337079437 ...
  -0.005576330816383  -0.006164001691543  -0.005754062064439  -0.004309274472218  -0.001970426319432   0.000961438765796   0.004046870129072   0.006776388499930   0.008645605987166   0.009240968293367   0.008322441053120 ...
   0.005874213402211   0.002133668380914  -0.002421919924152  -0.007121985518996  -0.011193372773283  -0.013867321807791  -0.014516842924145  -0.012761976776499  -0.008561028547337  -0.002257436921919   0.005424241333902 ...
   0.013425379897192   0.020483477433190   0.025249802768168   0.026542762979240   0.023385426091555   0.015319047829725   0.002328728212808  -0.014965847926289  -0.035401819339027  -0.057421399900803  -0.079150710430929 ...
  -0.098627522345843  -0.114047328036113  -0.123942750987368   0.872648601668912  -0.123942750987368  -0.114047328036113  -0.098627522345843  -0.079150710430929  -0.057421399900803  -0.035401819339027  -0.014965847926289 ...
   0.002328728212808   0.015319047829725   0.023385426091555   0.026542762979240   0.025249802768168   0.020483477433190   0.013425379897192   0.005424241333902  -0.002257436921919  -0.008561028547337  -0.012761976776499 ...
  -0.014516842924145  -0.013867321807791  -0.011193372773283  -0.007121985518996  -0.002421919924152   0.002133668380914   0.005874213402211   0.008322441053120   0.009240968293367   0.008645605987166   0.006776388499930 ...
   0.004046870129072   0.000961438765796  -0.001970426319432  -0.004309274472218  -0.005754062064439  -0.006164001691543  -0.005576330816383  -0.004172337079437  -0.002249136008113  -0.000153554698261   0.001768436318758 ...
   0.003231644105433   0.004051381673564   0.004158015837248   0.003603438954240   0.002528473584716   0.001153092856976  -0.000290114823546  -0.001559677092757  -0.002480137602476  -0.002932256415101  -0.002883343413527 ...
  -0.002382210801555  -0.001532267620955  -0.000485908218947   0.000592373907412   0.001547735038285   0.002239935097602   0.002558105616276   0.002434430195933   0.001833779324514   0.000773873754742  -0.000656289897754 ...
  -0.002258584801907  -0.003656865797247  -0.004190745808755  -0.002772127735242   0.002315432349574   0.013677737052159  -0.008230911208292];
        
         iqhp = filter(b,1,iq,[],1);
        iqhp = iqhp(length(b):end,:,:,:);

        vThresh = 0.14*v.vNyq; %No, this is not correct
        
        
        
    case 'FIR138_3' % More narrow stopband: wstop = 0.045 vNYq, wpass = 0.08 vNyq, Astop = 69
        
       b = [0.013974507740562  -0.004610553690858  -0.003959583887660  -0.003421437717813  -0.002960559310507  -0.002546987603092  -0.002152407693477  -0.001752369109738  -0.001335666591534  -0.000886985569956 ...
      -0.000406397631234   0.000105239633353   0.000636228776527   0.001171627782228   0.001691442937575   0.002171473383531   0.002586613469418   0.002909053988966   0.003114143802294   0.003178699822964 ...
       0.003086226298184   0.002825822207217   0.002394929610520   0.001799244313407   0.001053101010350   0.000180580071057  -0.000786262822472  -0.001805922295758  -0.002831500727849  -0.003809890804760 ...
      -0.004687127561322  -0.005408537732333  -0.005923021230674  -0.006185023274673  -0.006155931408084  -0.005809908523186  -0.005131231888127  -0.004123857624046  -0.002803648953984  -0.001207644615326 ...
       0.000614581185567   0.002598863723559   0.004660609482286   0.006719040134487   0.008639305794054   0.010372033472662   0.011759174363705   0.012712332141398   0.013147237722652   0.012976492549904 ...
       0.012124612380978   0.010533980201414   0.008197019674078   0.005103604492666   0.001282897291549  -0.003220624315483  -0.008330358995150  -0.013945505108362  -0.019943186225561  -0.026180338307836 ...
      -0.032504133383748  -0.038751252020856  -0.044758360378768  -0.050362088368919  -0.055407601522240  -0.059751812084099  -0.063269556969409  -0.065858778396590  -0.067443137132184   0.932023291287906 ...
      -0.067443137132184  -0.065858778396590  -0.063269556969409  -0.059751812084099  -0.055407601522240  -0.050362088368919  -0.044758360378768  -0.038751252020856  -0.032504133383748  -0.026180338307836 ... 
      -0.019943186225561  -0.013945505108362  -0.008330358995150  -0.003220624315483   0.001282897291549   0.005103604492666   0.008197019674078   0.010533980201414   0.012124612380978   0.012976492549904 ... 
       0.013147237722652   0.012712332141398   0.011759174363705   0.010372033472662   0.008639305794054   0.006719040134487   0.004660609482286   0.002598863723559   0.000614581185567  -0.001207644615326 ...
      -0.002803648953984  -0.004123857624046  -0.005131231888127  -0.005809908523186  -0.006155931408084  -0.006185023274673  -0.005923021230674  -0.005408537732333  -0.004687127561322  -0.003809890804760 ...
      -0.002831500727849  -0.001805922295758  -0.000786262822472   0.000180580071057   0.001053101010350   0.001799244313407   0.002394929610520   0.002825822207217   0.003086226298184   0.003178699822964 ...
       0.003114143802294   0.002909053988966   0.002586613469418   0.002171473383531   0.001691442937575   0.001171627782228   0.000636228776527   0.000105239633353  -0.000406397631234  -0.000886985569956 ...
      -0.001335666591534  -0.001752369109738  -0.002152407693477  -0.002546987603092  -0.002960559310507  -0.003421437717813  -0.003959583887660  -0.004610553690858   0.013974507740562];
        
        
     iqhp = filter(b,1,iq,[],1);
     iqhp = iqhp(length(b):end,:,:,:);

     vThresh = []; %No
       
     
    case 'FIR138_4' %wstop 0.025, wpass 0.065
        
        b = [
   0.001509115042861  -0.000034525537241  -0.000068058043150  -0.000123441184113  -0.000199799236827  -0.000295772847031 ...
  -0.000409631826660  -0.000538837471052  -0.000680320688218  -0.000830280521160  -0.000984244713277  -0.001137115984294 ...
  -0.001283183529831  -0.001416312534708  -0.001529945280526  -0.001617424560718  -0.001672024421941  -0.001687164029112 ...
  -0.001657053829376  -0.001575488334601  -0.001438248371948  -0.001241557281363  -0.000982959153697  -0.000661777398953 ...
  -0.000278410320538   0.000164319152671   0.000661899818034   0.001207757782797   0.001793239408264   0.002407737460812 ...
   0.003038634878843   0.003671550350141   0.004290419378162   0.004877802000780   0.005415116787498   0.005882910701854 ...
   0.006261412573742   0.006530491058815   0.006670602795222   0.006662933083476   0.006488650693590   0.006137430238845 ...
   0.005590187740414   0.004838424012058   0.003874213151504   0.002692974003530   0.001293857245185  -0.000320330967458 ...
  -0.002142727120984  -0.004162520312261  -0.006364910204055  -0.008731233479658  -0.011239196381512  -0.013863068795619 ...
  -0.016574184652672  -0.019341150382749  -0.022130543077621  -0.024907398199958  -0.027634995587018  -0.030278486455168 ...
  -0.032800546367721  -0.035166399563847  -0.037342653999334  -0.039297982495329  -0.041004244412832  -0.042436247418235 ...
  -0.043572840812700  -0.044397087617836  -0.044896627480886   0.954936011298596  -0.044896627480886  -0.044397087617836 ...
  -0.043572840812700  -0.042436247418235  -0.041004244412832  -0.039297982495329  -0.037342653999334  -0.035166399563847 ...
  -0.032800546367721  -0.030278486455168  -0.027634995587018  -0.024907398199958  -0.022130543077621  -0.019341150382749 ...
  -0.016574184652672  -0.013863068795619  -0.011239196381512  -0.008731233479658  -0.006364910204055  -0.004162520312261 ...
  -0.002142727120984  -0.000320330967458   0.001293857245185   0.002692974003530   0.003874213151504   0.004838424012058 ...
   0.005590187740414   0.006137430238845   0.006488650693590   0.006662933083476   0.006670602795222   0.006530491058815 ...
   0.006261412573742   0.005882910701854   0.005415116787498   0.004877802000780   0.004290419378162   0.003671550350141 ...
   0.003038634878843   0.002407737460812   0.001793239408264   0.001207757782797   0.000661899818034   0.000164319152671 ...
  -0.000278410320538  -0.000661777398953  -0.000982959153697  -0.001241557281363  -0.001438248371948  -0.001575488334601 ...
  -0.001657053829376  -0.001687164029112  -0.001672024421941  -0.001617424560718  -0.001529945280526  -0.001416312534708 ...
  -0.001283183529831  -0.001137115984294  -0.000984244713277  -0.000830280521160  -0.000680320688218  -0.000538837471052 ...
  -0.000409631826660  -0.000295772847031  -0.000199799236827  -0.000123441184113  -0.000068058043150  -0.000034525537241 ...
   0.001509115042861];
     
    
     iqhp = filter(b,1,iq,[],1);
     iqhp = iqhp(length(b):end,:,:);

     vThresh = []; %No
     
     
    case 'FIR276' % Should be used with one transmit angle, similar observation time as when used with three angles and FIR92_4
        b = [ -0.0191    0.0062    0.0052    0.0044    0.0036    0.0029    0.0024    0.0018    0.0014    0.0010    0.0006    0.0003    0.0000   -0.0002   -0.0004   -0.0006   -0.0008   -0.0009   -0.0011   -0.0012   -0.0013 ...
           -0.0014   -0.0014   -0.0015   -0.0015   -0.0015   -0.0014   -0.0014   -0.0013   -0.0012   -0.0011   -0.0010   -0.0009   -0.0007   -0.0005   -0.0003   -0.0001    0.0001    0.0003    0.0006    0.0008    0.0010 ...
            0.0012    0.0014    0.0016    0.0018    0.0019    0.0021    0.0022    0.0022    0.0023    0.0023    0.0022    0.0022    0.0020    0.0019    0.0017    0.0015    0.0012    0.0009    0.0006    0.0003   -0.0001 ...
           -0.0005   -0.0009   -0.0013   -0.0017   -0.0020   -0.0024   -0.0027   -0.0030   -0.0033   -0.0035   -0.0037   -0.0039   -0.0040   -0.0040   -0.0040   -0.0038   -0.0037   -0.0034   -0.0031   -0.0027   -0.0023 ...
           -0.0018   -0.0013   -0.0007   -0.0000    0.0007    0.0014    0.0021    0.0028    0.0035    0.0042    0.0049    0.0055    0.0061    0.0066    0.0071    0.0074    0.0077    0.0078    0.0079    0.0078    0.0076 ...
            0.0072    0.0067    0.0061    0.0053    0.0043    0.0033    0.0020    0.0007   -0.0008   -0.0024   -0.0041   -0.0059   -0.0078   -0.0098   -0.0118   -0.0139   -0.0160   -0.0181   -0.0202   -0.0222   -0.0242 ...
           -0.0262   -0.0281   -0.0298   -0.0315   -0.0330   -0.0344   -0.0356   -0.0367   -0.0376   -0.0383   -0.0388   -0.0391    0.9608   -0.0391   -0.0388   -0.0383   -0.0376   -0.0367   -0.0356   -0.0344   -0.0330 ...
           -0.0315   -0.0298   -0.0281   -0.0262   -0.0242   -0.0222   -0.0202   -0.0181   -0.0160   -0.0139   -0.0118   -0.0098   -0.0078   -0.0059   -0.0041   -0.0024   -0.0008    0.0007    0.0020    0.0033    0.0043 ...
            0.0053    0.0061    0.0067    0.0072    0.0076    0.0078    0.0079    0.0078    0.0077    0.0074    0.0071    0.0066    0.0061    0.0055    0.0049    0.0042    0.0035    0.0028    0.0021    0.0014    0.0007 ...
           -0.0000   -0.0007   -0.0013   -0.0018   -0.0023   -0.0027   -0.0031   -0.0034   -0.0037   -0.0038   -0.0040   -0.0040   -0.0040   -0.0039   -0.0037   -0.0035   -0.0033   -0.0030   -0.0027   -0.0024   -0.0020 ...
           -0.0017   -0.0013   -0.0009   -0.0005   -0.0001    0.0003    0.0006    0.0009    0.0012    0.0015    0.0017    0.0019    0.0020    0.0022    0.0022    0.0023    0.0023    0.0022    0.0022    0.0021    0.0019 ...
            0.0018    0.0016    0.0014    0.0012    0.0010    0.0008    0.0006    0.0003    0.0001   -0.0001   -0.0003   -0.0005   -0.0007   -0.0009   -0.0010   -0.0011   -0.0012   -0.0013   -0.0014   -0.0014   -0.0015 ...
           -0.0015   -0.0015   -0.0014   -0.0014   -0.0013   -0.0012   -0.0011   -0.0009   -0.0008   -0.0006   -0.0004   -0.0002    0.0000    0.0003    0.0006    0.0010    0.0014    0.0018    0.0024    0.0029    0.0036 ...
            0.0044    0.0052    0.0062   -0.0191];

        iqhp = filter(b,1,iq,[],1);
        iqhp = iqhp(length(b):end,:,:);

        vThresh = 0.14*v.vNyq; %????????

        
          
    case 'fonFilt23'
    
    b =  [ -0.0026    0.0019    0.0129   -0.0067   -0.0206   -0.0156    0.0269    0.0613    0.0206   -0.1179   -0.2822    0.6440   -0.2822   -0.1179    0.0206    0.0613    0.0269 ...
            -0.0156   -0.0206   -0.0067    0.0129    0.0019 -0.0026];
    
    iqhp = filter(b,1,iq,[],1);
    iqhp = iqhp(length(b):end,:,:);

    vThresh = 0.39*v.vNyq; %fdatool: wstop = 0.15*vNyq; wpass = 0.45*vNyq
    
    case 'adaptiveFIR93'
     
        Fs=1;
        Hd=fdesign.highpass('Fst,Fp,Ast,Ap',v.stopband/2/v.vNyq,(v.stopband/v.vNyq+0.054)/2,70,1,Fs);
        d=design(Hd,'equiripple');
        b = d.Numerator;
    
        iqhp = filter(b,1,iq,[],1);
        iqhp = iqhp(length(b):end,:,:);
        vThresh = []; %No
        
    case 'lowFIR93' % Equal to lowest cutoff adaptive FIR for Ingvilds bifurcation
        
        b = [
          -0.011925193998808
          -0.000642038673895
          -0.000535245935194
          -0.000349455399694
          -0.000076555648469
           0.000275568123219
           0.000706473260132
           0.001201834495338
           0.001756690457315
           0.002350791050173
           0.002973741950271
           0.003599595425583
           0.004213072227120
           0.004783840362063
           0.005293445520966
           0.005709302310164
           0.006012105985762
           0.006169457521860
           0.006163664733480
           0.005964943024577
           0.005559642989763
           0.004922984988376
           0.004047809914178
           0.002916749775083
           0.001531565004264
          -0.000115425878858
          -0.002011487647337
          -0.004152732287728
          -0.006514388667572
          -0.009076745826816
          -0.011804025460354
          -0.014688879154183
          -0.017696821149016
          -0.020742493289271
          -0.023828260459422
          -0.026892680560760
          -0.029892630205109
          -0.032779936528326
          -0.035508677034626
          -0.038036149209270
          -0.040319242527647
          -0.042322364890432
          -0.044009209517448
          -0.045353845816081
          -0.046330076131927
          -0.046924025502099
           0.952877907483894
          -0.046924025502099
          -0.046330076131927
          -0.045353845816081
          -0.044009209517448
          -0.042322364890432
          -0.040319242527647
          -0.038036149209270
          -0.035508677034626
          -0.032779936528326
          -0.029892630205109
          -0.026892680560760
          -0.023828260459422
          -0.020742493289271
          -0.017696821149016
          -0.014688879154183
          -0.011804025460354
          -0.009076745826816
          -0.006514388667572
          -0.004152732287728
          -0.002011487647337
          -0.000115425878858
           0.001531565004264
           0.002916749775083
           0.004047809914178
           0.004922984988376
           0.005559642989763
           0.005964943024577
           0.006163664733480
           0.006169457521860
           0.006012105985762
           0.005709302310164
           0.005293445520966
           0.004783840362063
           0.004213072227120
           0.003599595425583
           0.002973741950271
           0.002350791050173
           0.001756690457315
           0.001201834495338
           0.000706473260132
           0.000275568123219
          -0.000076555648469
          -0.000349455399694
          -0.000535245935194
          -0.000642038673895
          -0.011925193998808];

        iqhp = filter(b,1,iq,[],1);
         iqhp = iqhp(length(b):end,:,:);
    
         vThresh = []; %No
        
         
    case 'highFIR93' % Equal to adaptive FIR frame 37 for Ingvild Bifurcation
        
        b = [-0.006096404013133
   0.014348589278224
   0.001694896914077
  -0.002920264795190
  -0.004254011143246
  -0.004280773880190
  -0.003745078007576
  -0.002913031847791
  -0.001842571784658
  -0.000580177489918
   0.000819190406033
   0.002252889134803
   0.003569046687366
   0.004610003059083
   0.005200901885699
   0.005196675644604
   0.004508994470569
   0.003113614506037
   0.001096288519572
  -0.001368542088921
  -0.004024451456639
  -0.006557054121968
  -0.008615174971930
  -0.009858467271977
  -0.009994039283383
  -0.008832068392328
  -0.006326819998224
  -0.002593618586867
   0.002101247059564
   0.007326449865009
   0.012505696202353
   0.016964350329035
   0.020039813786675
   0.021098163746773
   0.019594292876578
   0.015160632545542
   0.007690637395564
  -0.002688265741933
  -0.015611219674609
  -0.030419242259312
  -0.046243794539055
  -0.062153540777811
  -0.077028396444037
  -0.089887773692295
  -0.099806660537227
  -0.106052641568141
   0.891799573550450
  -0.106052641568141
  -0.099806660537227
  -0.089887773692295
  -0.077028396444037
  -0.062153540777811
  -0.046243794539055
  -0.030419242259312
  -0.015611219674609
  -0.002688265741933
   0.007690637395564
   0.015160632545542
   0.019594292876578
   0.021098163746773
   0.020039813786675
   0.016964350329035
   0.012505696202353
   0.007326449865009
   0.002101247059564
  -0.002593618586867
  -0.006326819998224
  -0.008832068392328
  -0.009994039283383
  -0.009858467271977
  -0.008615174971930
  -0.006557054121968
  -0.004024451456639
  -0.001368542088921
   0.001096288519572
   0.003113614506037
   0.004508994470569
   0.005196675644604
   0.005200901885699
   0.004610003059083
   0.003569046687366
   0.002252889134803
   0.000819190406033
  -0.000580177489918
  -0.001842571784658
  -0.002913031847791
  -0.003745078007576
  -0.004280773880190
  -0.004254011143246
  -0.002920264795190
   0.001694896914077
   0.014348589278224
  -0.006096404013133];
         

         iqhp = filter(b,1,iq,[],1);
         iqhp = iqhp(length(b):end,:,:);
    
         vThresh = []; %No

    case 'adaptiveFIR139'
        
        Fs=1;
        Hd=fdesign.highpass('Fst,Fp,Ast,Ap',v.stopband/2/v.vNyq,(v.stopband/v.vNyq+0.0355)/2,70,1,Fs);
        d=design(Hd,'equiripple');
        b = d.Numerator;
        
        iqhp = filter(b,1,iq,[],1);
        iqhp = iqhp(length(b):end,:,:);
        vThresh = []; %No
        
        
    case 'neo_1' % wstop 0.1, wpass 0.17, astop 80, length 80
       b = [0.005545706923913  -0.013085028665167   0.000534721349586   0.004744363751183   0.004918611388015   0.003489118898781   0.001480786871068  -0.000636923849776  -0.002584252303306  -0.004043188841420  -0.004708342260579 ...
          -0.004349734803425  -0.002898955627535  -0.000473141593234   0.002482647305280   0.005365158989120   0.007445695518378   0.008052892405559   0.006792451856590   0.003599916807490  -0.001080354860612  -0.006399883787872 ...
          -0.011179266071593  -0.014159623120739  -0.014273522489920  -0.010927082794540  -0.004213296772896   0.004983984346084   0.015034776219935   0.023784215761666   0.028868247234843   0.028174895067493   0.020234545093671 ...
           0.004602163492289  -0.017976779123005  -0.045558245942589  -0.075230066859058  -0.103506164459546  -0.126851520534738  -0.142230602743398   0.852402623508901  -0.142230602743398  -0.126851520534738  -0.103506164459546 ...
          -0.075230066859058  -0.045558245942589  -0.017976779123005   0.004602163492289   0.020234545093671   0.028174895067493   0.028868247234843   0.023784215761666   0.015034776219935   0.004983984346084  -0.004213296772896 ...
          -0.010927082794540  -0.014273522489920  -0.014159623120739  -0.011179266071593  -0.006399883787872  -0.001080354860612   0.003599916807490   0.006792451856590   0.008052892405559   0.007445695518378   0.005365158989120 ...
           0.002482647305280  -0.000473141593234  -0.002898955627535  -0.004349734803425  -0.004708342260579  -0.004043188841420  -0.002584252303306  -0.000636923849776   0.001480786871068   0.003489118898781   0.004918611388015 ...
           0.004744363751183   0.000534721349586  -0.013085028665167   0.005545706923913];

         iqhp = filter(b,1,iq,[],1);
         iqhp = iqhp(length(b):end,:,:);
    
         vThresh = []; %No
         
         
    case 'neo_2' % wstop 0.15, wpass 0.22, astop 80, length 80
        b = [ 0.005218525022673  -0.012870990866268   0.002835104661452   0.005955543736387   0.004005710879372   0.000703402495238  -0.002238085302377  -0.003977481706405  -0.004100848064498  -0.002567918187732   0.000151072730177 ...
           0.003178338280429   0.005303768834988   0.005526134976573   0.003438478800078  -0.000515032343730  -0.004963708975530  -0.008117423116830  -0.008376192561406  -0.005108328354003   0.000964516449089   0.007791996593694 ...
           0.012633175426034   0.013063844359155   0.008025194431668  -0.001443269514199  -0.012288974187625  -0.020258530585254  -0.021401753163502  -0.013672267010399   0.001863530062002   0.020775863608538   0.036172443924266 ...
           0.040596547839781   0.028406822257535  -0.002145827841500  -0.047815509010583  -0.100709556159676  -0.150087880160012  -0.185102428584608   0.802242194486113  -0.185102428584608  -0.150087880160012  -0.100709556159676 ...
          -0.047815509010583  -0.002145827841500   0.028406822257535   0.040596547839781   0.036172443924266   0.020775863608538   0.001863530062002  -0.013672267010399  -0.021401753163502  -0.020258530585254  -0.012288974187625 ...
          -0.001443269514199   0.008025194431668   0.013063844359155   0.012633175426034   0.007791996593694   0.000964516449089  -0.005108328354003  -0.008376192561406  -0.008117423116830  -0.004963708975530  -0.000515032343730 ...
           0.003438478800078   0.005526134976573   0.005303768834988   0.003178338280429   0.000151072730177  -0.002567918187732  -0.004100848064498  -0.003977481706405  -0.002238085302377   0.000703402495238   0.004005710879372 ...
           0.005955543736387   0.002835104661452  -0.012870990866268   0.005218525022673];

         iqhp = filter(b,1,iq,[],1);
         iqhp = iqhp(length(b):end,:,:);
    
         vThresh = []; %No

    
        
        
    case  'SVD3D_cleaner' % Based on idea by Jerome Barangers et al, paper on adaptive SVD clutter filtering (review)
        [packet,ranges,beams,angles] = size(iq);
        M = reshape(iq,[size(iq,1), size(iq,2)*size(iq,3)*size(iq,4)]); 
        M = M.'; % Follow convention from paper (nx*ny*nz, nt)
        clear iq;
        tic
        [U,S,V] = svd(M,'econ'); % Economy version where only the nececcary column vectors of U is returned (nt)
               
        Q = abs(U);
        meanQ = mean( Q, 1);
        stdQ = std( Q, [], 1);

        detrQ = bsxfun(@minus, Q, meanQ);
        sigmaNM = stdQ.'*stdQ;

        C = detrQ'*detrQ./sigmaNM/size( M,1); % This is the socalled similarity matrix
        
        CDdet_1D = v.VDdet_1D;  
        estClutDim = v.estClutDim;
        estBloodDim = v.estBloodDim;
        bloodFraction = v.bloodFraction; % For testing purposes
         
        if CDdet_1D % Ingvild og Jørgens supermetode ;)
              
            % Look at spread of energy along diagonal of matrix
            flipC = fliplr(C);
            moment = zeros( 1, size(M,2)/2);
            for CD = 1:size(M,2)/2 % All possible clutter dimensions
                mydiag = diag( flipC, 2*CD-1);
                ln = length( mydiag);
                halfln = (ln-1)/2;
                arm = (-halfln:halfln).';
                moment(CD) = sum( abs(arm).*mydiag )/sum( abs( arm) );
            end

            moment = fliplr( moment);

            if estClutDim % Estimates the end of the assumend clutter region (last eigenvector representing clutter)
                % Find lokal minima. This will give the threshold for clutter space dimension
                invmoment = smooth(max(moment(:))-moment,7,'moving');
                MPP = 0.01;
                [mypeaks,locs] = findpeaks(invmoment,'MinPeakProminence',MPP);
                
                if isempty(mypeaks)
                     MPP = 0.005;
                     [mypeaks,locs] = findpeaks(invmoment,'MinPeakProminence',MPP);
                end           
                if isempty(mypeaks)
                     MPP = 0.001;
                     [mypeaks,locs] = findpeaks(invmoment,'MinPeakProminence',MPP);
                end
                
                if isempty(mypeaks)
                   ix = v.maxClutDimFac*packet;
                    sprintf('No peak found. Using max allowed clutter dim, at eig # %i',ix)

                else
                    mymaxpeak = mypeaks(1);
                    ix = locs(1);
                    sprintf('Using first peak. Clutter peak at eig # %i',ix)
                end

                figure(100)
                plot( invmoment);
                hold on; plot(ix,mymaxpeak,'o'); hold off
                b = ix;
            elseif estBloodDim % Estimates the central blood eigenvector, and picks out a set of blood eigenvectors based on the bloodFraction parameter
                MPP = 0.005;
                sMoment = smooth(moment,10,'moving');
                [mypeaks,locs] = findpeaks(sMoment,'MinPeakProminence',MPP);
                
                mymaxpeak = mypeaks(end);
                bPeak = locs(end);
                sprintf('Using last peak. Blood peak at eig # %i',bPeak)
             
                figure(102)
                plot(moment,'r','LineWidth',2);
                hold on; plot(sMoment,'k','LineWidth',2);
                hold on; plot(bPeak,mymaxpeak,'og','LineWidth',2)
                hold off;
                
                bEigs = (bPeak-floor(bloodFraction*size(M,2)/2):bPeak+floor(bloodFraction*size(M,2)/2));
                ix = 1:size(M,2);
                ix(bEigs) = [];
                b = bPeak;
            end
        
        else % This is the original approach by Baranger, with fitting of two squares
            
            lnC = size( C,1);
            Cd = C - eye(lnC);
            CVec = Cd(:);
            CVec = CVec-mean(CVec);
            stdC = std(CVec);
            corrVal = zeros( lnC, lnC);
            
             for kk = 1:floor(lnC*v.maxClutDimFac)
                for mm = kk+1:floor(lnC*v.maxBloodDimFac)
                    compMat = zeros( lnC, lnC);
                    compMat(1:kk, 1:kk) = 1;
                    compMat(kk+1:mm, kk+1:mm) = 1;
                    compMat(logical(eye(size(compMat)))) = 0;
                    compVec = compMat(:);
                    corrVal(kk, mm) = mean( CVec.*(compVec-mean(compVec) )/stdC/std(compVec) );
                end
         
            end

            [maxval, ind] = max( corrVal(:) );
            ix = mod( ind, lnC);
 
            disp([num2str( ix) ' eigenvectors!']);  

        end
        
        if 0 % Use this to look at spatial singular vectors images
            for kk = 1:packet
                u = reshape(U(:,kk),ranges,beams,angles);
                figure(100); 
                imagesc(10*log10(interp2(sum(abs(u).^2,3),3))); title(sprintf('u %i',kk)); caxis([-50 -20]); 

                pause();
            end
        end

        % Hard threshold for filtering out tissue signal
        sigma = diag(S); % Singular values in descending order
        
        if estClutDim
            sigma(1:ix) = 0;
        elseif estBloodDim
            sigma(ix) = 0;
        end
       
        Fm = eye(packet).*sigma;
        Fiq = U*Fm*V';
      
        iqhp = reshape(Fiq.',packet,ranges,beams,angles); 
        vThresh =0;
        toc
        
        
        if CDdet_1D && estBloodDim % Shows a power doppler image of the selected set of eigenvectors
           
            R0hp = iqhp.*conj(iqhp);
            figure(101);
            imagesc(10*log10(interp2(squeeze(mean(mean(R0hp,1),4)),3))); 
            title(sprintf('PD,eigs %i to %i',min(bEigs),max(bEigs)));
            colorbar;
            
        end
        
    case 'no'
        iqhp = iq;
        vThresh = 0;
        
        
        
    otherwise
        disp('Unknown filter method')
        
end







end