%% change filepath
cd(fileparts(mfilename('fullpath')));
no
%% import data

[ans1, ans2] = uigetfile({'../data/*.mat'});  % choose file to load
setting.data_file   = strcat(ans2,ans1);         % filepath of the data file

load(setting.data_file)

clear ans ans1 ans2

%% set constants
c_const     = 299792458;                        % m/s
% frequency for HeNe light: c/freq = (3e8 (m / s)) / (0.633 um) = 4.73933649 × 10^14 hertz

chop = 1;
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
im_w = fftshift( fft(im_all(:,:,chop:end), w_points, 3 ) );       % fft, trying the whole 3D matrix
%im_w = fftshift( fft(im_all(:,:,chop:end)-bsxfun(@times,mean(im_all,3),ones(1,1,arm.N-chop+1)), w_points, 3 ) );       % fft, trying the whole 3D matrix

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
%im_w2 = fftshift(im_w);

figure(5)
plot(fr0, squeeze(real(im_w(300,400,:))))
figure(6)
plot(fr0, squeeze(abs(im_w(300,400,:))))
figure(7)
plot(fr0, unwrap(squeeze(angle(im_w(300,400,:)))))

%% remove reference w/ gaussian

% generate gaussian mask
w       = gausswin(w_points, 18);                               % generate gaussian with FWHM 16
% shift the gaussian to the left
w       = circshift(w,floor(0.40e15/max(fr0)*w_points/2));      % at f= 4.0*1e14 Hz (750nm)
%w       = circshift(w,floor(0.47e15/max(fr0)*w_points/2));      % at f= 4.7*1e14 Hz (633nm)
w2 = reshape(w, [1,1,w_points]);

%{
% generate heaviside mask
h.width     = 35;
h.center    = 100;
heav        = heaviside((-w_points/2:w_points/2-1)'-h.center+h.width) - heaviside((-w_points/2:w_points/2-1)'-h.center-h.width);
%}

% how does the gaussian act?
figure(8)
plot(fr0, squeeze(abs(im_w(300,400,:))) )
hold on
plot(fr0, max(squeeze(abs(im_w(300,400,:)))).*w)
%plot(fr0, max(squeeze(abs(im_w2))).*heav)

% plot absolute value of fft modulated by gaussian
figure(9)
plot(fr0, abs(squeeze(im_w(300,400,:))).*w )

% plot phase value of fft modulated by gaussian
figure(10)
plot(fr0, unwrap(angle(squeeze(im_w(300,400,:)))) )

%% modulate the spectrum
%im_w3 = abs(squeeze(im_w2)) .* w    .* exp( 1i.*angle(squeeze(im_w2)).*w );
im_w4 = bsxfun(@times,im_w, w2);% .* exp( 1i.*bsxfun(@times, angle(im_w2), w) );
clear im_w
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


%% unwrap phase
%surface75 = squeeze((angle(im_t3(:,:,75))));
%s = surface75;

su = zeros(pic_size(1),pic_size(2),att.N-chop+1);

for n=1:arm.N-chop+1
    
    s = squeeze((angle(im_t3(:,:,n))));

    %s(301,401:end)  = unwrap(s(301,401:end));
    % unwrap the center line
    s(301,401:end)  = unwrap(s(301,401:end));
    s(301,1:401)    = flip( unwrap( flip(s(301,1:401)) ) );
    
    s(301:end,:) = unwrap(s(301:end,:),[],1);
    s(1:301,:) = flip( unwrap( flip(s(1:301,:),1) ,[],1), 1);
    
    s(301:end,401)  = unwrap(s(301:end,401));
    s(1:301,401)    = flip( unwrap( flip(s(1:301,401))));
    
    s(:,401:end) = unwrap(s(:,401:end),[],2);
    s(:,1:401) = flip( unwrap( flip(s(:,1:401),2) ,[],2), 2);
    su (:,:,n) = s;
end

%figure(17)
%mesh(s)
%plot(tau, abs(im_t2))
%plot(arm.vect_0(chop:end), abs(im_t2))
%figure(14)
%plot(tau, unwrap(angle(im_t2)))
%plot(arm.vect_0(chop:end), unwrap(angle(im_t2)))

%% 
% I=abs(im_t3(:,:,:));

% x = 1:1:600;
% y = 1:1:800;
% [xm,ym,zm] = meshgrid(x,y,tau);
% [xm,ym,zm] = meshgrid(y,x,tau(1:10));
% is = isosurface(xm,ym,zm,I, mean(mean(mean(I))) );
% is = isosurface(xm,ym,zm,I(:,:,1:10), mean(mean(mean(I(:,:,1:10)))) );

% for n = 1:256
% pause(1);
% imagesc(I(:,:,n))
% end