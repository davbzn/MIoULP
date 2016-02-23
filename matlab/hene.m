%%
z=Zaber.Controller;

%%

wc=CameraViewerStd('WinVideo.Controller','ImageViewerStd',figure,'Webcam Viewer');
set(wc.processor.flatten_colour,'Value',1);


%%
stepsize_mm = 630/5*1e-9*1e3;
z.move_absolute(2,20.013787+22*stepsize_mm, 'mm')

%%
stepsize_mm = 630/30*1e-9*1e3;
z.move_relative(2,stepsize_mm,'mm')

%% take image:
%im=wc.camera.take_snapshot;

%figureNF('image');clf
%imagesc(mean(im,3))

%%

% attuator parameters
att.center      = 20;                                   % mm
att.step_min    = 0.047625/1000;                        % um/1000 = mm
att.step        = att.step_min;                         % mm
att.N           = pow2(6); % step number as power of 2: 8 --> 256; 9 --> 512
att.min         = att.center -  att.N   *att.step/2;    % mm
att.max         = att.center + (att.N-1)*att.step/2;    % mm
att.vect        = [att.min:att.step:att.max];           % mm
att.vect_read   = [];

% arm parameters
arm.center      = 0;                                    % mm
arm.step_min    = 2*att.step_min;                       % mm
arm.step        = 2*att.step;                           % mm
arm.N           = att.N;    % step number
arm.min         = arm.center -  arm.N   *arm.step/2;    % mm from the center
arm.max         = arm.center + (arm.N-1)*arm.step/2;    % mm from the center
arm.vect        = [arm.min:arm.step:arm.max];           % mm
arm.vect_read   = [];

% camera
pic_size   = size(wc.camera.take_snapshot);
% resolution 800x600
% intensity 0-511

% prepare variables
im_store   = zeros(pic_size(1),pic_size(2),pic_size(3),10);
im_all     = zeros(pic_size(1),pic_size(2),att.N);
im         = zeros(pic_size(1),pic_size(2));

% move attuator before the start (solves the problem of inverting the direction)
z.move_absolute(2, att.vect(1)-0.2, 'mm');
z.move_absolute(2, att.vect(1)-0.1, 'mm');
%z.get_position(2,'mm')

% move attuator to the start
z.move_absolute(2, att.vect(1), 'mm');
pause(0.5);
att.vect_read(1) = z.get_position(2,'mm');
arm.vect_read(1) = 0;

for j=1:10
    im_store(:,:,:,j)=wc.camera.take_snapshot;
    pause(0.1);
end
    
% get average of ten photos and 3 pixels
im = mean(mean(im_store,4),3);
im_all(:,:,1) = im;

for n=2:att.N
    % move attuator to position
    %z.move_absolute(2, att.vect(n), 'mm');
    z.move_relative(2, att.step, 'mm');
    pause(0.1);
    att.vect_read(n) = z.get_position(2,'mm');
    arm.vect_read(n) = 2*(att.vect_read(n) - att.vect_read(n-1)) + arm.vect_read(n-1);

    % take ten photos to average
    for j=1:10
        im_store(:,:,:,j)=wc.camera.take_snapshot;
        %pause(0.1);
    end
    
    % get average of ten photos and 3 pixels
    im = mean(mean(im_store,4),3);
    im_all(:,:,n) = im;
    
end
clear im im_store n j ans

arm.vect_0 = arm.vect_read - mean(arm.vect_read);
