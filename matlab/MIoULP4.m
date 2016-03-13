%% import data

clear

if 0    % change this value manually, depending on what you need
    [ans1, ans2] = uigetfile({'../data/*.mat'});  % choose file to load
    sett.data_file   = strcat(ans2,ans1);         % filepath of the data file

    load(sett.data_file)
else
    load('2016-03-07_TiS_data.mat')
    load('2016-03-07_TiS_spider_data.mat')
end

clear ans ans1 ans2

%% commenting section
sett.comment.graph      = 'set graphs visible or not';
sett.comment.pic_size   = 'standard resolution';
sett.comment.processed  = 'cropped resolution';
sett.comment.N_pic      = 'numbers of picture taken';
sett.comment.x          = 'list of x-index for pixel used - binning';
sett.comment.y          = 'list of y-index for pixel used - binning';
sett.comment.binx       = 'number of pixels binned in the x direction';
sett.comment.biny       = 'number of pixels binned in the y direction';
sett.comment.ND_len = '[nm] width of ND filters';
arm.comment.step_m  = '[nm] average step through the data';
ax.comment.c_const  = '[nm/fs] speed of light in nm/fs';
ax.comment.wl_hene  = '[nm] wavelenght of HeNe cw laser';
ax.comment.f_hene   = '[PHz] opt. frequency of HeNe cw laser (0.473933649 PHz)';
ax.comment.wl_800   = '[nm] wavelenght at 800 nm';
ax.comment.f_800    = '[PHz] opt. frequency of HeNe cw laser (0.374740572 PHz)';
ax.comment.wl_760   = '[nm] wavelenght at 760 nm';
ax.comment.f_760    = '[PHz] opt. frequency of HeNe cw laser (0.394736842 PHz)';
ax.comment.tau      = '[fs] vector of delays in fs, from arm.vect_c';
ax.comment.dtau     = '[fs] delay step in fs, from arm.step_m';
ax.comment.v        = 'vertical coordinate of dummy point';
ax.comment.h        = 'horizontal coordinate of dummy point';
ax.comment.cv       = 'vertical coordinate of the center point';
ax.comment.ch       = 'horizontal coordinate of the center point';
ax.comment.NF       = 'number of frequency points, next power of 2 from arm.N';
ax.comment.fr       = '[PHz] vector of frequencies, with length NF';
ax.comment.fr0      = '[PHz] vector of frequencies centered on 0 [ ax.fr0(NF/2+1) = 0 ]';
ax.comment.wlen     = '[nm] vector of wavelengths, from ax.fr0';
ax.comment.g        = '[a.u.] gaussian form, elevated to the 8 and centered at ax.f_760';
ax.comment.NDph     = 'phase as function of frequency, given by the filter considering schott/N-BK7 material from refractiveindex.com';
ax.comment.p        = 'linear fit in the bandwidth of the pulse';
ax.comment.NDph_2   = 'phase without the linear part';
in.comment.f        = '[a.u.] data fft on tau axis (3) and shifted to the center (fftshift)';
in.comment.fmod     = '[a.u.] data modulated by gaussian form';
in.comment.phase    = 'phase from the spider, interpolated on the ax.fr0 scale';
in.comment.int      = '[a.u.] intensity from the spider, interpolated on the ax.fr0 scale';
in.comment.fphased  = '[a.u.] data modulated by ND filter phase, spider phase and intensity (reference)';
in.comment.center   = '[a.u.] intensity at the center of the xy plane';
in.comment.center_n = 'intensity at the center normalised';
in.comment.avg      = '[a.u.] average intensity on the xy plane for each time/delay';
in.comment.avg_n    = 'average intensity on the xy plane normalised';
spi.ax.comment.fr   = '[PHz] vector of frequencies from the spider';

%% define script settings
% set graph visible or not
sett.graph          = true;
%sett.graph  = true;

%% set constants
ax.c_const  = 299.792458;                               % nm/fs
% wavelength and frequency for HeNe light (632.8 nm):
ax.wl_hene  = 632.8;                                    % nm
ax.f_hene   = ax.c_const/ax.wl_hene;                    % f = 0.473933649 PHz (or 10^15 Hz)
% wavelength and frequency relative to 800nm:
ax.wl_800   = 800.0;                                    % nm
ax.f_800    = ax.c_const/ax.wl_800;                     % f = 0.374740572 PHz (or 10^15 Hz)
% wavelength and frequency relative to 760nm:
ax.wl_760   = 760.0;                                    % nm
ax.f_760    = ax.c_const/ax.wl_760;                     % f = 0.394736842 PHz (or 10^15 Hz)

arm.step_m  = mean(diff(arm.vect_c));                   % nm

ax.tau      = arm.vect_c/ax.c_const;                    % nm/(nm/fs) = fs
ax.dtau     = arm.step_m/ax.c_const;                    % nm/(nm/fs) = fs

% define coordinates for a test point
ax.v        = ceil(sett.processed(1)/2); %center of the vertical axis
ax.h        = 120; %ceil(sett.processed(2)/2);
% define coordinates for the pixel at the center
ax.cv       = ceil(sett.processed(1)/2);
ax.ch       = ceil(sett.processed(2)/2);

if sett.graph
    %% plot graphs
    figure
    plot(ax.tau, squeeze(in.all(ax.v,ax.h, :)), '-*' )
    title('intensity'); xlabel('fs'); ylabel('a.u.');
end

%% FFT on data
% get vector of frequencies
ax.NF   = pow2( nextpow2(arm.N) ) ;                         % number of frequency points is next power of 2
ax.fr   = ( (1:1:ax.NF) ./ ax.dtau ) ./ ax.NF ;             % array of frequencies, with length NF
ax.fr0  = ( (-ax.NF/2:ax.NF/2-1) ./ ax.dtau ) ./ ax.NF ;    % array of frequencies centered on 0
ax.wlen = ax.c_const./ax.fr0;                               % array of wavelengths

% FFT on the data
in.f    = fftshift( fft(in.all, ax.NF, 3 ), 3);             % fft on data, with shifting

%clear in.all

if sett.graph
    %% plot graphs
    figure
    plot(ax.fr0, squeeze( abs( in.f(ax.v,ax.h,:) ) ) );
    hold on
    plot(ax.fr0, zeros(ax.NF,1));
    title('intensity spectrum (in Hz)'); xlabel('PHz'); ylabel('a.u.');
    
    figure
    plot(ax.fr0, unwrap( squeeze( angle( in.f(ax.v,ax.h,:) ) ) ) );
    title('phase spectrum (in Hz)'); xlabel('PHz'); ylabel('a.u.');
end

%% remove reference w/ gaussian

% generate gaussian mask
ax.g       = mygauss(-ax.NF/2:ax.NF/2-1,120,10)';% generate gaussian with FWHM 16
ax.g       = circshift(ax.g,floor(ax.f_760/max(ax.fr0)*ax.NF/2)); % shift the gaussian to the left

if sett.graph
    %% how does the gaussian act? Testing area
    figure
    plot(ax.fr0, squeeze(abs(in.f(ax.v,ax.h,:))) );
    title('spectrum and gaussian compared');
    xlabel('frequency [PHz]'); ylabel('intensity: a.u.');
    hold on;
    plot(ax.fr0, max(squeeze(abs(in.f(ax.v,ax.h,ax.NF/2+50:end)))).*ax.g,'r');

    % plot absolute value of fft modulated by gaussian
    figure
    plot(ax.fr0, abs(squeeze(in.f(ax.v,ax.h,:))).*ax.g );
    title('spectrum modulated by gaussian');
    xlabel('frequency [PHz]'); ylabel('intensity: a.u.');

    % plot phase value of fft modulated by gaussian
    figure
    plot(ax.fr0, unwrap(angle(squeeze(in.f(ax.v,ax.h,:)))) )
    title('phase NOT modulated by gaussian');
    xlabel('frequency [PHz]'); ylabel('intensity: a.u.');    
end

%% modulate the spectrum
ax.g        = reshape(ax.g, [1,1,ax.NF]);                   % reshape to use bsxfun
in.fmod     = bsxfun(@times,in.f, ax.g);                    % modulate spectrum using bsxfun

%clear in.f

if sett.graph
	%% spectrum in wavelength
    figure1 = figure;

    % Create axes
    axes1 = axes('Parent',figure1);
    title({'Spectrum of the ultrafast pulse'});
    %xlim(axes1,[1.5e-07 1.55e-06]);
    xlim(axes1,[150 1450]);
    %ylim(axes1,[-207.907338442333 12837.6767481348]);
    box(axes1,'on');
    hold(axes1,'on');

    % Create plot
    plot( ax.wlen, squeeze( abs( in.fmod(ax.v,ax.h,:) ) ) );

    % Create labels
    xlabel({'\lambda [nm]'});
    ylabel({'Intensity [a.u.]'}); 
    
    clear figure1 axes1
end  

%% add various phases
% interpole data to match ax.tau vector
spi.ax.fr   = spi.rec.of.omega'./(2*pi);
in.phase    = interp1(spi.ax.fr, spi.rec.of.phify(spi.rec.y0indo(1),:,1), ax.fr0, 'linear', 0);
in.int      = interp1(spi.ax.fr, spi.rec.of.Ify(spi.rec.y0indo(1),:,1),   ax.fr0, 'linear', 0);

sect = [(630:1024)]'; %[(318:415),(615:710)]';
% create vector of phase for ND filter
sett.ND_len = 5.41e6;  % nm
ax.NDph     = (2*pi/ax.c_const*sett.ND_len).*ax.fr0.*(real(n_BK7(ax.fr0)) -1 );
ax.p        = polyfit(ax.fr0(sect),ax.NDph(sect),1);
ax.NDph_2   = squeeze(ax.NDph) - (ax.fr0.*ax.p(1)+ax.p(2));
ax.NDph_3   = zeros(1,1024);
ax.NDph_3(sect)   = ax.NDph_2(sect);

if sett.graph
    figure
    plot(ax.fr0,unwrap(angle(squeeze(in.f(ax.cv,ax.ch,:)))))
    hold on
    plot(ax.fr0,unwrap(angle(squeeze(in.f(ax.cv,ax.ch,:)).*exp(1i.*ax.NDph_2'))))
    hold on
    plot(ax.fr0,unwrap(angle(squeeze(in.f(ax.cv,ax.ch,:)).*exp(1i.*ax.NDph'))))
    title('Phase');
    xlabel('Optical Frequency [PHz]')
end

% add phase to in.fmod_save
in.phase    = reshape(in.phase, [1,1,ax.NF]);
in.int      = reshape(in.int, [1,1,ax.NF]);
ax.NDph_3   = reshape(ax.NDph_3, [1,1,ax.NF]);
ax.NDph     = reshape(ax.NDph, [1,1,ax.NF]);

% which of the next two should I use?
%in.int(find(in.int==0)) = 1;
in.int      = in.int+eps;
%in.int(find(in.phase==0)) = inf;

in.fphased  = bsxfun(@times,in.fmod, exp(1i*in.phase).*exp(1i.*ax.NDph) ./sqrt(in.int));   

clear ax.g %in.fmod

%% inverse FFT
%in.t    = ifft( ifftshift(in.fmod, 3), arm.N, 3 );
in.t    = ifft( ifftshift(in.fphased, 3), arm.N, 3 );
%clear in.fphased

%% I(t)
in.center   = squeeze(abs(in.t(ax.cv,ax.ch,:)));
in.center_n = in.center./max(in.center);
in.avg      = squeeze(mean(mean(abs(in.t), 1), 2));
in.avg_n    = in.avg./max(in.avg);
%in.avg      = squeeze(mean(mean(abs(in.phased), 1), 2));

if 1
    %%
    figure
    plot(ax.tau, in.center_n)
    title('I(t)'); xlabel('Time [fs]'); ylabel('Normalized Intensity [a.u.]')
    grid on;
    hold on;
    plot(ax.tau, in.avg_n)
    % set(gcf,'color','w')
    % xlim( [400,1200]);
    legend('center','average');
end

%% unwrap phase % ! variables names have not been converted in this section !

%{ 
% % find max of intensity
% % define fittype
% ft = fittype( 'poly22' );
% % prepare arrays for fit
% s = squeeze((angle(in.t(:,:,1))));
% %[xData, yData, ~] = prepareSurfaceData( sett.x, sett.y, s );
% [xData, yData, ~] = prepareSurfaceData( (1:sett.processed(1))', (1:sett.processed(2)), s );
% 
% coeff = zeros(arm.N, 6);
% 
% %figure
% %for n=1:arm.N
% n=600;
%     
%     s = squeeze((angle(in.t(:,:,n))));
% 
%     s = unwrap(s,[],2);
%     s = unwrap(s,[],1);
%     
%     %% Fit model to data.
%     [ph(n).fit, ph(n).gof]   = fit( [xData, yData], s(:), ft );
%     coeff(n,1) = ph(n).fit.p00;
%     coeff(n,2) = ph(n).fit.p10;
%     coeff(n,3) = ph(n).fit.p01;
%     coeff(n,4) = ph(n).fit.p20;
%     coeff(n,5) = ph(n).fit.p11;
%     coeff(n,6) = ph(n).fit.p02;
%     
%     if sett.graph
%         %% Plot fit with data.
%         figure( 'Name', 'untitled fit 1' );
%         h = plot( fitresult, [xData, yData], zData );
%         legend( h, 'untitled fit 1', 's vs. x, y', 'Location', 'NorthEast' );
%         % Label axes
%         xlabel x
%         ylabel y
%         zlabel s
%         grid on
%     end
% 
%     n % print 'n' to undestand the progress of the script
% %end
% 
% clear xData yData
% clear xx yy ft s n
% 
% if sett.graph
%     %%
%     figure
%     plot(coeff(:,4));
%     hold on;
%     plot(coeff(:,6),'r');
%     
%     figure
%     plot(coeff(:,5));
%     
%     figure
%     plot(coeff(:,2));
%     hold on;
%     plot(coeff(:,3),'r');
% end
% 
% %% add spherical wave phase
% [xData, yData] = meshgrid(1:sett.processed(1),1:sett.processed(2));
% zData      = feval(ph(600).fit,xData,yData); %array of phase evaluated on every pixel
% zData      = reshape(zData, [sett.processed(1),sett.processed(2),1]);
% in.t1      = bsxfun(@times, in.t, exp(1i.*zData) );
% 
% clear xData yData zData in.t
% 
% %% add temporal phase
% % get spectral phase file
% [ans1, ans2] = uigetfile({'../data/*.mat'});  % choose file to load
% sett.data_file   = strcat(ans2,ans1);         % filepath of the data file
% 
% load(sett.data_file)
% 
% clear ans ans1 ans2
% %%
% % interpole data to match ax.tau vector
% t_ph        = rec.ot.t';
% % deprecated
% %angE        = angle( mean(rec.ot.Ety(rec.y0indo(1),:,:),3) ); % unwrap(angle( mean(rec.ot.Ey(:,:,:),3) ));
% %ph_t        = interp1( t_ph.*1e-15, angE , ax.tau, 'linear', 0);
% Etau        = interp1(t_ph.*1e-15,mean(rec.ot.Ety(rec.y0indo(1),:,:),3),ax.tau,'linear',0);
% ph_t        = squeeze(angle(Etau));
% %%
% % add phase to in.fmod_save
% ph_t        = reshape(ph_t, [1,1,ax.NF]);
% in.phased   = bsxfun(@times,in.t1, exp(1i*ph_t));   
% 
% clear in.t1
%}

%% isosurface plot % ! variables names have not been converted in this section !
%{
if 0
    %% example
    [xx, yy, zz] = meshgrid(-100:100, -100:100, -100:100);
    V = exp(-(xx.^2+yy.^2+zz.^2)/1000); % just a gauss shape in 3D
    p = patch(isosurface(xx, yy, zz, V, 0.5));
    isonormals(xx,yy,zz,V, p)
    p.FaceColor = 'red';
    p.EdgeColor = 'none';
    daspect([1 1 1])
    view(3) %or" camup([1 0 0 ]); campos([25 -55 5]) "
    camlight; lighting phong
end

% with real data
if 0
    I     = abs(in.phased(:,:,:));
    Inorm = I./max( max( max(I)));

    [xm,ym,zm] = meshgrid(sett.x,sett.y,ax.tau);
    sett.pp     = 2.2e-6; % pixel pitch in m
    is = patch(isosurface(xm.*pp, ym.*pp, zm.*pp, I, mean(mean(mean(I))) ) );
    is = patch(isosurface(xm.*pp, ym.*pp, zm.*pp, I,  max( max( max(I))) ) );
    is = patch(isosurface(xm.*pp, ym.*pp, zm.*pp, Inorm,  0.8 ) );
    % is = isosurface(xm,ym,zm,I(:,:,1:10), mean(mean(mean(I(:,:,1:10)))) );
    camlight; lighting phong
end

if 0
    figure
    for n = 1:1024
        pause(1);
        imagesc(I(:,:,n))
    end
end
%}
