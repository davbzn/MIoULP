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
%sett.graph  = false;
sett.graph  = true;

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
g       = gausswin(ax.NF, 16);                              % generate gaussian with FWHM 16
g       = circshift(g,floor(ax.f_800/max(ax.fr0)*ax.NF/2)); % shift the gaussian to the left

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
end

%% add spherical wave phase

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
[xData, yData, zData] = prepareSurfaceData( sett.x, sett.y, s );
clear zData

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

    n
end

clear xData yData
clear xx yy ft s n

figure
plot(coeff(:,4));
hold on;
plot(coeff(:,6));