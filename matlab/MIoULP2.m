%% dummy section
no

%% import data

[ans1, ans2] = uigetfile({'../data/*.mat'});  % choose file to load
sett.data_file   = strcat(ans2,ans1);         % filepath of the data file

load(sett.data_file)

clear ans ans1 ans2

%% define computing settings
% set number of parts to which divide the file
%sett.divide = 5;

% set graph visible or not
sett.graph  = false;
%sett.graph  = true;

%% set constants
c_const     = 299792458;                                % m/s
% wavelength and frequency for HeNe light (632.8 nm):
ax.wl_hene  = 632.8e-9;                                 % m
ax.f_hene = c_const/ax.wl_hene;                         % f = c/wlen = 0.473933649 × 10^15 Hz
% wavelength and frequency relative to 800nm:
ax.wl_800   = 800.0e-9;                                 % m
ax.f_800    = c_const/ax.wl_800;                        % f = c/wlen = 0.374740572 × 10^15 Hz
% wavelength and frequency relative to 760nm:
ax.wl_760   = 760.0e-9;                                 % m
ax.f_760    = c_const/ax.wl_760;                        % f = c/wlen = 0.394736842 × 10^14 Hz

arm.step_0  = mean(diff(arm.vect_0));

ax.tau      = arm.vect_0*0.001/c_const;                 % mm*(m/mm)/(m/s) = s
ax.dtau     = arm.step_0*0.001/c_const;                 % mm*(m/mm)/(m/s) = s

% define coordinates for a dummy point
ax.v        = ceil(sett.processed(1)/2); %center of the vertical axis
ax.h        = 120;%ceil(sett.processed(2)/2);
% define coordinates for the pixel at the center
ax.cv       = ceil(sett.processed(1)/2);
ax.ch       = ceil(sett.processed(2)/2);

if sett.graph
    %% plot graphs
    figure
    plot(ax.tau, squeeze(im_all(ax.v,ax.h, :)), '-*' )
end

%% FFT on data
% get vector of frequencies
ax.NF   = pow2( nextpow2(arm.N) ) ;                         % number of frequency points is next power of 2
ax.fr   = ( (1:1:ax.NF) ./ ax.dtau ) ./ ax.NF ;             % array of frequencies, with length NF
ax.fr0  = ( (-ax.NF/2:ax.NF/2-1) ./ ax.dtau ) ./ ax.NF ;    % array of frequencies centered on 0
ax.wlen = c_const./ax.fr0;                                  % array of wavelengths

% FFT on the data
im_f    = fftshift( fft(im_all, ax.NF, 3 ), 3);               % fft on data, with shifting

clear im_all

if sett.graph
    %% plot graphs
    figure
    plot(ax.fr0, squeeze( abs( im_f(ax.v,ax.h,:) ) ) );
    title('');
    figure
    plot(ax.fr0, unwrap( squeeze( angle( im_f(ax.v,ax.h,:) ) ) ) );
    title('');
end

%% remove reference w/ gaussian

% generate gaussian mask
%g       = gausswin(ax.NF, 16);                             % deprecated
g       = mygauss(1:ax.NF,16,8);% generate gaussian with FWHM 16
g       = circshift(g,floor(ax.f_760/max(ax.fr0)*ax.NF/2)); % shift the gaussian to the left

if sett.graph
    %% how does the gaussian act? Testing area
    figure
    plot(ax.fr0, squeeze(abs(im_f(ax.v,ax.h,:))) );
    title('spectrum and gaussian compared');
    xlabel('frequency [Hz]'); ylabel('intensity: a.u.');
    hold on;
    plot(ax.fr0, max(squeeze(abs(im_f(ax.v,ax.h,:)))).*g);

    % plot absolute value of fft modulated by gaussian
    figure
    plot(ax.fr0, abs(squeeze(im_f(ax.v,ax.h,:))).*g );
    title('spectrum modulated by gaussian');
    xlabel('frequency [Hz]'); ylabel('intensity: a.u.');

    % plot phase value of fft modulated by gaussian
    figure
    plot(ax.fr0, unwrap(angle(squeeze(im_f(ax.v,ax.h,:)))) )
    title('phase NOT modulated by gaussian');
    xlabel('frequency [Hz]'); ylabel('intensity: a.u.');
    
end

%% modulate the spectrum
g       = reshape(g, [1,1,ax.NF]);                          % reshape to use bsxfun
im_fmod = bsxfun(@times,im_f, g);                           % modulate spectrum using bsxfun

clear im_f g

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
    plot( ax.wlen.*1e9, squeeze( abs( im_fmod(ax.v,ax.h,:) ) ) );

    % Create labels
    xlabel({'\lambda [nm]'});
    ylabel({'Intensity [a.u.]'}); 
    
    clear figure1 axes1
end                                  % save data for later

%% inverse FFT
im_t    = ifft( ifftshift(im_fmod, 3), arm.N, 3 );
clear im_fmod   % not sure about this part,
                % don't we need this as it is A_1(x,y,w) ?
                % then we retrieve the x,y dependence of phi_r working on the phase
                % while we retrieve the w dependence of phi_r using the SEA-SPIDER

%% unwrap phase

% find max of intensity
% define fittype
ft = fittype( 'poly22' );
% prepare arrays for fit
s = squeeze((angle(im_t(:,:,1))));
%[xData, yData, ~] = prepareSurfaceData( sett.x, sett.y, s );
[xData, yData, zData] = prepareSurfaceData( (1:processed(1))', (1:processed(2)), s );

coeff = zeros(arm.N, 6);

%figure
for n=1:arm.N
    
    s = squeeze((angle(im_t(:,:,n))));

    s = unwrap(s,[],2);
    s = unwrap(s,[],1);
    
    %% Fit model to data.
    [ph(n).fit, ph(n).gof]   = fit( [xData, yData], s(:), ft );
    coeff(n,1) = ph(n).fit.p00;
    coeff(n,2) = ph(n).fit.p10;
    coeff(n,3) = ph(n).fit.p01;
    coeff(n,4) = ph(n).fit.p20;
    coeff(n,5) = ph(n).fit.p11;
    coeff(n,6) = ph(n).fit.p02;
    
    if sett.graph
        %% Plot fit with data.
        figure( 'Name', 'untitled fit 1' );
        h = plot( fitresult, [xData, yData], zData );
        legend( h, 'untitled fit 1', 's vs. x, y', 'Location', 'NorthEast' );
        % Label axes
        xlabel x
        ylabel y
        zlabel s
        grid on
    end

    n % print 'n' to undestand the progress
end

clear xData yData
clear xx yy ft s n

if sett.graph
    figure
    plot(coeff(:,4));
    hold on;
    plot(coeff(:,6));
end

%% add spherical wave phase
A = ?;%array of phase
A          = reshape(A, [sett.processed(1),sett.processed(2),1]);
im_t2      = bsxfun(@times, im_t, exp(1i.*A) );

%clear im_t

%% add temporal phase
% get spectral phase file
[ans1, ans2] = uigetfile({'../data/*.mat'});  % choose file to load
sett.data_file   = strcat(ans2,ans1);         % filepath of the data file

load(sett.data_file)

clear ans ans1 ans2

% interpole data to match ax.tau vector
t_ph        = rec.ot.t';
angE        = unwrap(2*pi*rand(2048,1)); % unwrap(angle( mean(rec.ot.Ey(:,:,:),3) ));
ph_t        = interp1( t_ph.*1e-15, angE , ax.tau, 'linear', 0);

% add phase to im_fmod_save
ph_t        = reshape(ph_t, [1,1,ax.NF]);
im_phased   = bsxfun(@times,im_t1, exp(1i*ph_t));   

%clear v im_t1

%% isosurface plot
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
    I     = abs(im_phased(:,:,:));
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