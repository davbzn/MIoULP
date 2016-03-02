%% dummy section
no

%% import data

[ans1, ans2] = uigetfile({'../data/*.mat'});  % choose file to load
sett.data_file   = strcat(ans2,ans1);         % filepath of the data file

load(sett.data_file)

clear ans ans1 ans2

%% set constants
c_const     = 299792458;                                % m/s
% wavelength and frequency for HeNe light (632.8 nm):
ax.wl_hene  = 632.8e-9;                                 % m
ax.f_hene = c_const/ax.wl_hene;                         % f = c/wlen = 0.473933649 × 10^15 Hz
% wavelength and frequency relative to 800nm:
ax.wl_800   = 800.0e-9;                                 % m
ax.f_800    = c_const/ax.wl_800;                        % f = c/wlen = 0.374740572 × 10^15 Hz

arm.step_0  = mean(diff(arm.vect_0));

ax.tau      = arm.vect_0*0.001/c_const;                 % mm*(m/mm)/(m/s) = s
ax.dtau     = arm.step_0*0.001/c_const;                 % mm*(m/mm)/(m/s) = s

% define coordinates for a dummy point
% in this case is at the center of the sensor
ax.v        = ceil(sett.processed(1)/2);
ax.h        = ceil(sett.processed(2)/2);
% define coordinates for the pixel at the center
ax.cv        = ceil(sett.processed(1)/2);
ax.ch        = ceil(sett.processed(2)/2);

if 0
    %% plot graphs
    figure
    plot(ax.tau, squeeze(im_all(ax.v,ax.h, :)), '-*' )
end

%% FFT on data
% get vector of frequencies
ax.NF   = pow2( nextpow2(arm.N) ) ;                         % number of frequency points is next power of 2
ax.fr   = ( (1:1:ax.NF) ./ ax.dtau ) ./ ax.NF ;             % array of frequencies, with length NF
ax.fr0  = ( (-ax.NF/2:ax.NF/2-1) ./ ax.dtau ) ./ ax.NF ;    % array of frequencies centered on 0
%wlen    = c_const./fr;                                     % array of wavelengths

% FFT on the data
im_f    = fftshift( fft(im_all, ax.NF, 3 ) );               % fft on data, with shifting

clear im_all

if 0
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
g       = gausswin(ax.NF, 16);                              % generate gaussian with FWHM 16
g       = circshift(g,floor(ax.f_800/max(ax.fr0)*ax.NF/2)); % shift the gaussian to the left

if 0
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

%% add spherical wave phase

%% inverse FFT
im_t    = ifft( fftshift(im_fmod), arm.N, 3 );
clear im_fmod   % not sure about this part,
                % don't we need this as it is A_1(x,y,w) ?
                % then we retrieve the x,y dependence of phi_r working on the phase
                % while we retrieve the w dependence of phi_r using the SEA-SPIDER

%% unwrap phase

% find max of intensity
% define fittype
ft = fittype( 'poly22' );
% prepare arrays for fit
%[xx, yy]    = meshgrid( (1:sett.processed(1))', 1:sett.processed(2));
[xx, yy]    = meshgrid( sett.x, sett.y);

%figure
for n=1:arm.N
    
    s = squeeze((angle(im_t(:,:,n))));

    s = unwrap(s,[],1);
    s = unwrap(s,[],2);
    
    [ph(n).fit, ph(n).gof]   = fit([xx(:), yy(:)], s(:), ft);
    
    n
end

clear xx yy ft s n