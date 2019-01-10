function [spectrums, midptab] = VectorTrackingDoppler_func( IQData, params)
%
x = params.x; 
z = params.z;
trackveltab = params.trackveltab; %(-5:0.01:5);
halflen = params.halflen; %30; % vindulengden blir da 2*halflen +1
step = params.step; % bestemmer steglengde i glidende vindu
PRF = params.PRF; %meta.maxPRF_Rec/meta.nrAngRec;
Xin = params.x_axis; %bfPars.x_axis.';
Zin = params.z_axis; %bfPars.z_axis.';
f0 = params.f0;
anglefact = 1/cosd( params.steerangle);

if isfield(params,'xvel')
    fixed = 0;
    xvel = params.xvel; %Velocity in point (x,z)
    zvel = params.zvel; %Velocity in point (x,z)
else
    fixed = 1;
    fixedAngle = params.fixedAngle;
end

useGPU = params.useGPU;
packetSizeFrame = size( IQData,3);
origo = [x z];
windowfunc = hamming(halflen*2+1);
midptab = halflen+1:step:size( IQData,3)-halflen;

if fixed
    angletab = ones(1,params.nFrames)*fixedAngle;
else
    angletab = zeros( 1, size( xvel, 2)  );
    vtot = zeros( 1, size( xvel, 2)  );
    for kk = 1:size( xvel, 2),
        angletab( kk) = atand( xvel(kk)./zvel(kk) );
        vtot( kk) = sqrt( xvel(kk).^2 + zvel(kk).^2 );
    end
end

framenrIn = 1:length( angletab);
framenrOut = midptab/(packetSizeFrame)+0.5;

if length( angletab) > 2,
    currangletab = interp1( framenrIn, angletab, framenrOut, 'linear', 'extrap' );
else
    currangletab = angletab*ones( size( framenrOut) );
end

currangletab( isnan( currangletab) ) = 0;
currangletab = smooth(currangletab,0.1,'rlowess');


if 0,
    xlims = [ 0 max( midptab) - min( midptab ) ]/PRF;

    disp('PW Doppler');

    PWsig = complex( zeros( size( IQData,3), 1) );
    for kk = 1:size( IQData, 3),
        PWsig(kk) = interp2(Xin, Zin, ...
            IQData(:,:,kk), x, z, 'linear', 0 );
        if mod(kk, 500) == 0,
            clc
            kk
        end
    end


    % iq: coloumn vector of iq signal
    % skip: number of samples interval between each fft
    % w: window function
    % Nv number of velocity bins
    % Nsmooth: temporal averaging filter
    %  Hans Torp

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

    PWcorrstatic = 54;

    figure(10); imagesc( xlims, vaxis/cosd(PWcorrstatic), 10*log10( abs( PWspect) ) );
    caxis([-dyn 0]-gain);
    colormap( gray)
    xlabel('Time [s]');
    ylabel('Velocity [m/s]');
    title(sprintf('PW Doppler, corrected with %i',PWcorrstatic));
    set( gca, 'FontSize', 18);
    set( gcf, 'Position', [753 212 656 502] );

    PWspectAcorr = zeros( length( trackveltab), size( PWspect,2) ); 
    for kk = 1:size( PWspect, 2),
        vin = vaxis./cosd( currangletab(kk)-10 );
        vout = trackveltab;
        PWspectAcorr(:,kk) = interp1(vin, PWspect(:,kk), vout);
    end


    figure(11); imagesc( xlims, trackveltab, 10*log10( abs( PWspectAcorr) ) );
    caxis([-dyn 0]-gain);
    colormap( gray)
    xlabel('Time [s]');
    ylabel('Velocity [m/s]');
    title('PW Doppler')
    set( gca, 'FontSize', 18);
    set( gcf, 'Position', [753 212 656 502] );
    set(gca,'YDir','Normal')
end

tic

mm = 1;

gpuSamps = 500;

allangles = zeros( 1, length(midptab) );
if useGPU
    g = gpuDevice;
%     spectcomps = gpuArray( zeros( halflen*2+1, length( trackveltab), 1 ) );
%     spectrums = gpuArray( zeros( length( trackveltab), length(midptab), 1 ) );
    spectrums = zeros( length( trackveltab), length(midptab), 1 );


    Xin_GPU = gpuArray( Xin);
    Zin_GPU = gpuArray( Zin);

    anglefact_GPU = gpuArray( anglefact);
    f0_GPU = gpuArray( f0);
    IQData_GPU = gpuArray( IQData(:,:,1:min( gpuSamps+2*halflen, size( IQData,3) ) ) );
else
%     spectcomps = zeros( halflen*2+1, length( trackveltab), 1 );
    spectrums = zeros( length( trackveltab), length(midptab), 1 );
end

figure(3000); 
ctr = 1;
offset = 0;
for locmidp = midptab,

    if useGPU && locmidp > offset + halflen + gpuSamps,
        offset = offset + gpuSamps;
        IQData_GPU = gpuArray( IQData(:,:, offset+1:min( offset+gpuSamps+2*halflen, size( IQData,3) ) ) );
    end

    currangle = currangletab( ctr);
   % sprintf('Line 163: currangle is %g',currangle)
    ctr = ctr+ 1;
    
    trackcoords = zeros(halflen*2+1, length( trackveltab), length( currangle), 2 );
    for angleind = 1:length(currangle),
        track_angle = currangle(angleind);
        if fixed
             vel_vector = [sind(track_angle) cosd(track_angle)];
        else
             vel_vector = [xvel(angleind) zvel(angleind)]/sqrt( xvel(angleind).^2+zvel(angleind).^2)*sign(zvel(angleind));
        end
        vpos = (-halflen:halflen).'*trackveltab/PRF;
        trackcoords(:,:,angleind,:) = reshape( vpos(:)*vel_vector, ...
                    [halflen*2+1 length( trackveltab) 2] );
    end

 %       sprintf('Line 179: vel_vector is %g %g',vel_vector(1),vel_vector(2))

    loc_trackcoords = trackcoords + repmat( permute( origo, [1 3 4 2]), ...
        [halflen*2+1 length( trackveltab) length( currangle)] );

   
    trackinds = locmidp-halflen:locmidp+halflen;                
    if useGPU,
        loc_trackcoords_GPU = gpuArray( loc_trackcoords);
        spectcomps = interp3(Xin_GPU, Zin_GPU, 1:length( trackinds), ...
                IQData_GPU(:,:,trackinds-offset), ...
                loc_trackcoords_GPU(:,:,:,1), ...
                loc_trackcoords_GPU(:,:,:,2), ...
                repmat( (1:length( trackinds) ).', 1, length( trackveltab) ), 'linear', 0 ).* ...
                exp(2*pi*1i*f0_GPU*( loc_trackcoords_GPU(:,:,:,2)*2/1540*anglefact_GPU ) );
    else
        spectcomps = interp3(Xin, Zin, 1:length( trackinds), ...
                IQData(:,:,trackinds), ...
                loc_trackcoords(:,:,:,1), ...
                loc_trackcoords(:,:,:,2), ...
                repmat( (1:length( trackinds) ).', 1, length( trackveltab) ), 'spline', 0 ).* ...
                exp(2*pi*1i*f0*( loc_trackcoords(:,:,:,2)*2/1540*anglefact ) );
    end    
        
    spectrum = sum( spectcomps ...
        .*repmat(windowfunc, [1 size( spectcomps,2), size( spectcomps, 3)]), 1);
    spectrums(:,mm, :) = gather( spectrum );
    allangles(mm) = currangle;

    clc

    mm = mm + 1;
end

toc
