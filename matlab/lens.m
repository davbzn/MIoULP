%% set constants

c_const     = 299792458;                                % m/s
% frequency for HeNe light:
% c/freq = (3e8 (m / s)) / (0.633 um)   = 4.73933649 × 10^14 Hz
%                                       = 0.473933649 × 10^15 Hz

%% get data locations

% get location for data without lens
[ans1, ans2] = uigetfile({'../data/*.mat'});                                % choose file to load
sett.data_file(1).str    = strcat(ans2,ans1);                            % filepath of the data file

% get location for data with lens
sett.nr_lens = 1;
for d=1:sett.nr_lens
    [ans1, ans2] = uigetfile({'../data/*.mat'});                            % choose file to load
    sett.data_file(1+sett.nr_lens).str    = strcat(ans2,ans1);        % filepath of the data file
end

clear ans ans1 ans2

%% load and elaborate data

for d=1:sett.nr_lens+1
    %% load data
    
    load(sett.data_file(d).str)
    
    if d==1
        % allocate memory to get coefficients
        coeff = zeros(arm.N,6,sett.nr_lens+1);
    end
    %% define or evaluate some values
    
    % number of frequency points is next power of 2
    ax.NF          = pow2( nextpow2(arm.N) );                  % number of points in the frequency domain
    arm.step_0  = mean(diff(arm.vect_0));
    axdtau        = arm.step_0*0.001/c_const;                 % mm*(m/mm)/(m/s) = s
    ax.max_f       = (0.5-1/ax.NF)/ax.dtau;                       % maximum frequency
    
    %% FFT on data

    im_f = fftshift( fft(im_all, ax.NF, 3 ) );                 % fft, trying the whole 3D matrix

	clear im_all
    
    %% remove reference w/ gaussian

    % generate gaussian mask
    g   = gausswin(ax.NF, 16);                               % generate gaussian with FWHM 16
    g   = circshift(g,floor(0.47e15/ax.max_f*ax.NF/2));         % shift the gaussian to the left
                                                                % at f= 4.7*1e14 Hz (633nm)
    % modulate the spectrum
    g   = reshape(g, [1,1,ax.NF]);                           % reshape to use bsxfun
    im_fmod = bsxfun(@times,im_f, g);
    clear im_f g
    
    %% inverse FFT
    im_t    = ifft( fftshift(im_fmod), arm.N, 3 );
    clear im_fmod
    
    %% unwrap and get coefficient
    
    % first index is VERTICAL axis and second is HORIZONTAL axis: s(V,H)
    lv  = length( im_t(:,1,1) );
    lh  = length( im_t(1,:,1) );
    
    % PROBLEM: should get a mean value of As and Bs
    for delay = 1:arm.N
        s   = squeeze((angle(im_t(:,:,delay))));
        
        % unwrap the center line along the horizontal axis
        sh  = unwrap( s( lv/2+1, :) );
        ph = polyfit(1:lh, sh, 2);
        
        sv  = unwrap( s( :, lh/2+1) );
        pv = polyfit(1:lv, sv', 2);
        
        % save coefficient in variable
        coeff(delay,1:3,d)  = ph; % (delay,(ah,bh,ch,av,bv,cv),file)
        coeff(delay,4:6,d)  = pv;
        
        clear s sh sv ph pv
    end
    %% test and clearing
    
    %{
    delay = 128;
    figure
    s   = squeeze((angle(im_t(:,:,delay))));
    plot(1:lh, unwrap(s( lv/2+1, :)), 'o');
    hold on
    plot(1:lh, coeff(delay,1,d).*(1:lh).^2 + coeff(delay,2,d).*(1:lh) + coeff(delay,3,d));
    figure
    plot(1:lv, unwrap(s( :, lh/2+1)), 'o');
    hold on
    plot(1:lv, coeff(delay,4,d).*(1:lv).^2 + coeff(delay,5,d).*(1:lv) + coeff(delay,6,d));
    %}
    
    clear delay im_t %lv lh
    
end
clear d

%% evaluate radius
nglass  = 1.5151; % from refractive.info
nair    = 1;
pp      = 2.8e-3; % 2.8um = 2.8e-3 mm % pixel pitch
ppv     = pp*1.2;
pph     = pp*1.6;
conv    = 0.633e-3/(2*pi); % wavelength/2*pi % conversion factor

% evaluation point
evpv    = (lv/2+1)*pp;
evph    = (lh/2+1)*pp;
%{
% create derivatives...
fdh     = conv*(2*evph/pph^2*(coeff(:,1,2)/nglass-coeff(:,1,1)/nair)+1/pph*(coeff(:,2,2)/nglass-coeff(:,2,1)/nair));
fdv     = conv*(2*evpv/ppv^2*(coeff(:,4,2)/nglass-coeff(:,4,1)/nair)+1/ppv*(coeff(:,5,2)/nglass-coeff(:,5,1)/nair));
sdh     = conv*(2/pph^2*(coeff(:,1,2)/nglass-coeff(:,1,1)/nair) );
sdv     = conv*(2/ppv^2*(coeff(:,4,2)/nglass-coeff(:,5,1)/nair) );

radius = zeros(arm.N, 2);
radius(:,1) = (abs( 1 + fdh.^2 )).^(1.5) ./ abs( sdh );
%radius(:,1) = ( abs(1 + conv^2.*(2*Dah/pp^2*evph+Dbh/pp).^2) ).^(3/2) ./abs(2*Dah/(pp^2)*conv);
%(abs(1+(coeff(:,2,1)-coeff(:,2,2)).^2)).^(1.5)./abs(coeff(:,1,1)-coeff(:,1,2));

radius(:,2) = (abs( 1 + fdv.^2 )).^(1.5) ./ abs( sdv );
%radius(:,2) = ( abs(1 + conv^2.*(2*Dav/pp^2*evpv+Dbv/pp).^2) ).^(3/2) ./abs(2*Dav/(pp^2)*conv);
%radius(:,2) =
%(abs(1+(coeff(:,5,1)-coeff(:,5,2)).^2)).^(1.5)./abs(coeff(:,4,1)-coeff(:,4,2));
%}

% create derivatives
fdh     = conv/(nglass-nair)*(2/pph^2 * evph * (coeff(:,1,2)-coeff(:,1,1)) + (coeff(:,2,2)-coeff(:,2,1))/pph );
fdv     = conv/(nglass-nair)*(2/ppv^2 * evpv * (coeff(:,4,2)-coeff(:,4,1)) + (coeff(:,5,2)-coeff(:,5,1))/ppv );
sdh     = conv/(nglass-nair)*(2/pph^2 * (coeff(:,1,2)-coeff(:,1,1)) );
sdv     = conv/(nglass-nair)*(2/ppv^2 * (coeff(:,4,2)-coeff(:,5,1)) );

radius = zeros(arm.N, 2);
% horizontal
radius(:,1) = (abs( 1 + fdh.^2 )).^(1.5) ./ abs( sdh );
% vertical
radius(:,2) = (abs( 1 + fdv.^2 )).^(1.5) ./ abs( sdv );

%% plot radii of curvature
figure
plot(radius(:,1)/58.3719); % horizontal
title('Radius of curvature'); % Create title
hold on
plot(radius(:,2)*58.3719); % vertical
xlabel('delay, normalised'); ylabel('mm') % Create labels
legend([],'horizontal','vertical'); % Create legend