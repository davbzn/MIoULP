%% change filepath
cd(fileparts(mfilename('fullpath')));

%% import data

[ans1, ans2] = uigetfile({'../data/*.mat'});  % choose file to load
setting.data_file   = strcat(ans2,ans1);         % filepath of the data file

load(setting.data_file)

clear ans ans1 ans2

%% set constants
c_const     = 299792458;                        % m/s
% frequency for HeNe light: c/freq = (3e8 (m / s)) / (0.633 um) = 4.73933649 × 10^14 hertz

chop = 20;
arm.step_0  = mean(diff(arm.vect_0));

tau         = arm.vect_0(chop:end)*0.001/c_const;      % mm*(m/mm)/(m/s) = s
dtau        = arm.step_0*0.001/c_const;       % mm*(m/mm)/(m/s) = s

figure(1)
plot(tau, squeeze(im_all(300,400, chop:end)), '-*' )
%data = squeeze(mean(mean(im_all(300,200:209,chop:end),1),2));

%% FFT on data
% get vector of frequencies
w_points = pow2( nextpow2(arm.N) ) ;            % number of frequency points is next power of 2
fr  = ( (1:1:w_points) ./ dtau ) ./ w_points ;  % create array of frequencies, with length w_points

% FT on the data
%im_w = fft(im_all(300,400,:), w_points );       % fft, trying with just a pixel for now
%im_w = fft(data, w_points );                    % fft, trying with just a pixel for now
im_w = fft(im_all, w_points, 3 );       % fft, trying the whole 3D matrix

lambdas = c_const./fr;                          % array of wavelengths

clear im_all

%% fiddle around with graphs
% graphs NOT circshifted
%figure(2)
%plot(fr, squeeze(real(im_w)))
%figure(3)
%plot(fr, squeeze(abs(im_w)))
%figure(4)
%plot(fr, squeeze(angle(im_w)))

% circshift & plot
fr0 = ( (-w_points/2:w_points/2-1) ./ dtau ) ./ w_points ;
im_w2 = fftshift(im_w);

%figure(5)
%plot(fr0, squeeze(real(im_w2)))
figure(6)
plot(fr0, squeeze(abs(im_w2(300,400,:))))
figure(7)
plot(fr0, squeeze(angle(im_w2(300,400,:))))

clear im_w

%% remove reference w/ gaussian

% generatre gaussian mask
w       = gausswin(w_points, 16);                               % generate gaussian with FWHM 16
w       = circshift(w,floor(0.47e15/max(fr0)*w_points/2));      % shift the gaussian to the left
                                                                % at f= 4.7*1e14 Hz (633nm)
w2 = reshape(w, [1,1,256]);

%{
% generate heaviside mask
h.width     = 35;
h.center    = 100;
heav        = heaviside((-w_points/2:w_points/2-1)'-h.center+h.width) - heaviside((-w_points/2:w_points/2-1)'-h.center-h.width);
%}

% how does the gaussian act?
figure(8)
plot(fr0, squeeze(abs(im_w2(300,400,:))) )
hold on
plot(fr0, max(squeeze(abs(im_w2(300,400,:)))).*w)
%plot(fr0, max(squeeze(abs(im_w2))).*heav)

% plot absolute value of fft modulated by gaussian
figure(9)
plot(fr0, abs(squeeze(im_w2(300,400,:))).*w )

% plot phase value of fft modulated by gaussian
figure(10)
plot(fr0, angle(squeeze(im_w2(300,400,:))) )

%% modulate the spectrum
%im_w3 = abs(squeeze(im_w2)) .* w    .* exp( 1i.*angle(squeeze(im_w2)).*w );
im_w4 = bsxfun(@times,im_w2, w2);% .* exp( 1i.*bsxfun(@times, angle(im_w2), w) );
clear im_w2
%% add spherical wave phase


%% fft back to time
% inverse fft, one pixel
%im_t    = ifft(im_w3, arm.N-chop+1 );

%figure(11)
%plot(tau, abs(im_t))
%plot(arm.vect_0(chop:end), abs(im_t))
%figure(12)
%plot(tau, unwrap(angle(im_t)))
%plot(arm.vect_0(chop:end), unwrap(angle(im_t)))

% ifft of the whole image
im_t2   = ifft(im_w4, arm.N-chop+1, 3 );
clear im_w4
im_t3   = fftshift(im_t2);
clear im_t2

chop2 = 1;
surface = unwrap(unwrap(-squeeze((angle(im_t3(chop2:end-chop2+1,chop2:end-chop2+1,50)))),[],1),[],2);
figure(13)
mesh(surface)
%plot(tau, abs(im_t2))
%plot(arm.vect_0(chop:end), abs(im_t2))
figure(14)
%plot(tau, unwrap(angle(im_t2)))
%plot(arm.vect_0(chop:end), unwrap(angle(im_t2)))

clear im_w4
