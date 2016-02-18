%% import data
%load('../data/hene_2016-02-17.mat')

%% set constants
c_const = 299792458;                        % m/s
% frequency for HeNe light: c/freq = (3e8 (m / s)) / (0.633 um) = 4.73933649 × 10^14 hertz
pos.N = 2*pos.n+1; % better use a power of 2. e.g. 256, 512, 1024, ...

tau         = positions*0.001/c_const;      % mm*(m/mm)/(m/s) = s
dtau        = pos.step*0.001/c_const;       % mm*(m/mm)/(m/s) = s

figure(1)
plot(tau, squeeze(im_all(300,400, :)) )

%% FFT on data
% get vector of frequencies
w_points = pow2( nextpow2(pos.N) ) ;            % number of frequency points is next power of 2
fr  = ( (1:1:w_points) ./ dtau ) ./ w_points ;  % create array of frequencies, with length w_points

% FT on the data
im_w = fft(im_all(300,400,:), w_points );       % fft, trying with just a pixel for now

lambdas = c_const./fr;                          % array of wavelengths

%% fiddle around with graphs
figure(2)
plot(fr, squeeze(real(im_w)))
figure(3)
plot(fr, squeeze(abs(im_w)))
figure(4)
plot(fr, squeeze(angle(im_w)))

% circshift & plot
fr0 = ( (-w_points/2:w_points/2-1) ./ dtau ) ./ w_points ;
im_w2 = fftshift(im_w);

figure(5)
plot(fr0, squeeze(real(im_w2)))
figure(6)
plot(fr0, squeeze(abs(im_w2)))
figure(7)
plot(fr0, squeeze(angle(im_w2)))

%% remove reference w/ gaussian

% generatre gaussian mask
alpha   = 12;                             % define gaussian parameter
w       = gausswin(w_points,alpha);     % generate gaussian
w       = circshift(w,155);             % shift the gaussian to the left    %2.5*10^(4)*circshift(w,90);

% generate heaviside mask
h.width     = 35;
h.center    = 155;
heav        = heaviside((-w_points/2:w_points/2-1)'-h.center+h.width) - heaviside((-w_points/2:w_points/2-1)'-h.center-h.width);

%N = 301;
%n = -(N-1)/2:(N-1)/2;
%stdev   = (w_points-1)/(2*alpha);
%y       = exp(-1/2*(n/stdev).^2); 
 
figure(8)                               % how does the gaussian act?
plot(fr, squeeze(abs(im_w2)) )
hold on
plot(fr, max(squeeze(abs(im_w2))).*w)
plot(fr, max(squeeze(abs(im_w2))).*heav)

%im_w_new=im_w2.*(w+1i*w);
%im_w_new = w.*(real(im_w2)+1i.*imag(im_w2));
% w=(2.5*10^(4))^(-1)*w;

figure(9)
plot(fr, abs(squeeze(im_w2)).*w )
%plot(fr, abs(im_w_new))

figure(10)
plot(fr, angle(squeeze(im_w2)).*w )
%plot(fr, angle(im_w_new))

im_w3 = abs(squeeze(im_w2)) .* w    .* exp( 1i.*angle(squeeze(im_w2)).*w );
im_w4 = abs(squeeze(im_w2)) .* heav .* exp( 1i.*angle(squeeze(im_w2)).*heav );

%% add spherical wave phase


%% fft back to time
im_t    = ifft(im_w3, pos.N );             % inverse fft
im_t2   = ifft(im_w4, pos.N );

figure(11)
plot(positions, abs(im_t))
figure(12)
plot(positions, abs(im_t2))