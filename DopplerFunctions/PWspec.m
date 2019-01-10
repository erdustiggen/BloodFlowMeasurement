function [PWspect, vaxis] = PWspec( PWsig, bfPars),
    disp('PW Doppler');
    
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
% 
%     PWcorrstatic = 54;
% 
%     figure(10); imagesc( xlims, vaxis/cosd(PWcorrstatic), 10*log10( abs( PWspect) ) );
%     caxis([-dyn 0]-gain);
%     colormap( gray)
%     xlabel('Time [s]');
%     ylabel('Velocity [m/s]');
%     title(sprintf('PW Doppler, corrected with %i',PWcorrstatic));
%     set( gca, 'FontSize', 18);
%     set( gcf, 'Position', [753 212 656 502] );
% 
%     PWspectAcorr = zeros( length( trackveltab), size( PWspect,2) ); 
%     for kk = 1:size( PWspect, 2),
%         vin = vaxis./cosd( currangletab(kk)-10 );
%         vout = trackveltab;
%         PWspectAcorr(:,kk) = interp1(vin, PWspect(:,kk), vout);
%     end
% 
% 
%     figure(11); imagesc( xlims, trackveltab, 10*log10( abs( PWspectAcorr) ) );
%     caxis([-dyn 0]-gain);
%     colormap( gray)
%     xlabel('Time [s]');
%     ylabel('Velocity [m/s]');
%     title('PW Doppler')
%     set( gca, 'FontSize', 18);
%     set( gcf, 'Position', [753 212 656 502] );
%     set(gca,'YDir','Normal')