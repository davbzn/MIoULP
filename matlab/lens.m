%% set constants

c_const     = 299792458;                                % m/s
% frequency for HeNe light:
% c/freq = (3e8 (m / s)) / (0.633 um)   = 4.73933649 × 10^14 Hz
%                                       = 0.473933649 × 10^15 Hz

%% get data locations

% get location for data without lens
[ans1, ans2] = uigetfile({'../data/*.mat'});            % choose file to load
setting.data_file(1).str    = strcat(ans2,ans1);        % filepath of the data file

% get location for data with lens
[ans1, ans2] = uigetfile({'../data/*.mat'});            % choose file to load
setting.data_file(2).str    = strcat(ans2,ans1);        % filepath of the data file

clear ans ans1 ans2

%% load and elaborate data

for d=1:2
    %% load data
    
    load(setting.data_file(d).str)
    
    if d==1
        % allocate memory to get coefficients
        coeff = zeros(arm.N,6,2);
    end
    %% define some values
    
    % number of frequency points is next power of 2
    f_points    = pow2( nextpow2(arm.N) );
    arm.step_0  = mean(diff(arm.vect_0));
    dtau        = arm.step_0*0.001/c_const;             % mm*(m/mm)/(m/s) = s
    max_f       = 0.5/dtau-1/f_points/dtau;             % maximum frequency
    
    %% FFT on data

    im_f = fftshift( fft(im_all, f_points, 3 ) );       % fft, trying the whole 3D matrix

	clear im_all
    
    %% remove reference w/ gaussian

    % generate gaussian mask
    g   = gausswin(f_points, 16);                               % generate gaussian with FWHM 16
    g   = circshift(g,floor(0.47e15/max_f*f_points/2));         % shift the gaussian to the left
                                                                % at f= 4.7*1e14 Hz (633nm)
    % modulate the spectrum
    g   = reshape(g, [1,1,f_points]);                           % reshape to use bsxfun
    im_fmod = bsxfun(@times,im_f, g);
    clear im_f g
    
    %% inverse FFT
    im_t    = ifft( fftshift(im_fmod), arm.N, 3 );
    clear im_fmod
    
    %% unwrap and get coefficient
    
    % first index is VERTICAL axis and second is HORIZONTAL axis: s(V,H)
    lv  = length( im_t(:,1,1) );
    lh  = length( im_t(1,:,1) );
    
    for delay = 1:arm.N
        s   = squeeze((angle(im_t(:,:,delay))));
        
        % unwrap the center line along the horizontal axis
        sh  = unwrap( s( lv/2+1, :) );
        ph = polyfit(1:lh, sh, 2);
        
        sv  = unwrap( s( :, lh/2+1) );
        pv = polyfit(1:lv, sv', 2);
        
        % save coefficient in variable
        coeff(delay,1:3,d)  = ph;
        coeff(delay,4:6,d)  = pv;
        
        clear s sh sv ph pv
    end
    clear delay im_t %lv lh
    
end
clear d

%% evaluate radius
pp      = 2.8e-3; % 2.8um = 2.8e-3 mm % pixel pitch
conv    = 0.633e-3/(2*pi*1.5); % wavelength/2*pi*refractive_index % conversion factor
Dah     = -coeff(:,1,1)+coeff(:,1,2);
Dbh     = -coeff(:,2,1)+coeff(:,2,2);
Dav     = -coeff(:,4,1)+coeff(:,4,2);
Dbv     = -coeff(:,5,1)+coeff(:,5,2);
evpv    = lv*pp;% evaluation point
evph    = lh*pp;

radius = zeros(arm.N,2);
radius(:,1) = ( 1 + conv^2.*(2*Dah/pp^2*evph+Dbh/pp).^2).^1.5 ./(2*Dah/pp^2*conv);
%(abs(1+(coeff(:,2,1)-coeff(:,2,2)).^2)).^(1.5)./abs(coeff(:,1,1)-coeff(:,1,2));

radius(:,2) = ( 1 + conv^2.*(2*Dav/pp^2*evpv+Dbv/pp).^2).^1.5 ./(2*Dav/pp^2*conv);
%radius(:,2) = (abs(1+(coeff(:,5,1)-coeff(:,5,2)).^2)).^(1.5)./abs(coeff(:,4,1)-coeff(:,4,2));

figure
plot(radius(:,1))
title('Radius of curvature');
hold on
plot(radius(:,2))
xlabel('delay, normalised'); ylabel('mm')