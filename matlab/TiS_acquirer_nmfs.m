%% start attuator (zaber) controller
z   = Zaber.Controller;

% commands are (usually ID=2)
% z.move_absolute(ID, position_in_mm, 'mm')
% z.move_relative(ID, position_in_mm, 'mm')

%% start camera controller
if 0
    wc=CameraViewerStd('WinVideo.Controller','ImageViewerStd',figure,'Webcam Viewer');
    set(wc.processor.flatten_colour,'Value',1);
else
    wc = CameraViewerStd('GENTL.Controller','ImageViewerStd',figure,'GENTL (AVT) Viewer');
    
    %% set binning
    wc.processor.binning.binx.Value = 4;
    wc.processor.binning.biny.Value = 4;
    set(wc.processor.binning.do,'Value',1)
end

%% take image for testing
if 0
    %%
    wc.camera.take_snapshot;
    im=wc.processor.image;
    x=wc.processor.x;
    y=wc.processor.y;
    figureNF('image');clf
    imagesc(x,y,im)
end

%% define settings, parameters and values

% general settings
sett.pic_size   = size(wc.camera.take_snapshot);
sett.processed  = size(wc.processor.image);
sett.N_pic      = 3;                                   % number of pictures
sett.x          = wc.processor.x;
sett.y          = wc.processor.y;
sett.binx       = wc.processor.binning.binx.Value;
sett.biny       = wc.processor.binning.biny.Value;

% actuator parameters
att.center      = 14.2e6;                               % nm
att.step_min    = 47.625;                               % nm
att.step        = att.step_min;                         % nm
att.N           = pow2(8); % step number as power of 2: 8 --> 256; 9 --> 512
att.min         = att.center -  att.N   *att.step/2;    % nm
att.max         = att.center + (att.N-1)*att.step/2;    % nm
att.vect        = att.min:att.step:att.max;             % nm
att.vect_read   = [];

% arm parameters
arm.center      = 0;                                    % nm
arm.step_min    = 2*att.step_min;                       % nm
arm.step        = 2*att.step;                           % nm
arm.N           = att.N;    % step number
arm.min         = arm.center -  arm.N   *arm.step/2;    % nm from the center
arm.max         = arm.center + (arm.N-1)*arm.step/2;    % nm from the center
arm.vect        = arm.min:arm.step:arm.max;             % nm
arm.vect_read   = [];

% prepare variables for processed images (grayscale)
in.store        = zeros(sett.processed(1), sett.processed(2), sett.N_pic);
in.all          = zeros(sett.processed(1), sett.processed(2), att.N);

%% start acquisition of data

% move attuator before the start (solves the problem of inverting the direction)
z.move_absolute(2, att.vect(1)-0.2, 'mm');
pause(0.5);
z.move_absolute(2, att.vect(1)-0.1, 'mm');
pause(0.5);

% move attuator to the start and read position
z.move_absolute(2, att.vect(1), 'mm');
pause(0.5);
att.vect_read(1) = z.get_position(2,'mm')*1e6;          % nm
arm.vect_read(1) = 0;                                   % nm

% take sett.N_pic photos to average
for j=1:sett.N_pic
    wc.camera.take_snapshot;
    in.store(:,:,j) = wc.processor.image;               % a.u.
end
    
% get average of sett.N_pic photos
in.all(:,:,1) = mean(in.store,3);                       % a.u.

for n=2:att.N
    % move attuator and read position
    z.move_relative(2, att.step, 'mm');
    att.vect_read(n) = z.get_position(2,'mm')*1e6;      % nm
    arm.vect_read(n) = 2*(att.vect_read(n) - att.vect_read(n-1)) + arm.vect_read(n-1); % nm

    % take sett.N_pic photos to average
    for j=1:sett.N_pic
        wc.camera.take_snapshot;
        in.store(:,:,j) = wc.processor.image;           % a.u.
    end
    
    % get average of sett.N_pic photos
    in.all(:,:,n) = mean(in.store,3);                   % a.u.
    
end
clear im_store n j ans

% shift arm.vect such that the n/2+1 position in equal to 0
arm.vect_c = arm.vect_read - arm.vect_read(att.N/2+1);  % nm

%% comment variables

sett.comment.pic_size       = 'size of standard resolution';
sett.comment.processed      = 'size of processed image';
sett.comment.N_pic          = 'number of pictures';
sett.comment.x              = 'x coordinates of pixels';
sett.comment.y              = 'y coordinates of pixels';
sett.comment.binx           = 'binning in the x direction';
sett.comment.biny           = 'binning in the y direction';

% actuator parameters
att.comment.center          = '[nm] attuator center position';
att.comment.step_min        = '[nm] minumum attuator step size';
att.comment.step            = '[nm] actual attuator step size, should be a multiple of att.step_min';
att.comment.N               = 'number of steps';
att.comment.min             = '[nm] min position of the attuator';
att.comment.max             = '[nm] max position of the attuator';
att.comment.vect            = '[nm] prefixed positons of the attuator';
att.comment.vect_read       = '[nm] actual positons of the attuator';

% arm parameters
arm.comment.center          = '[nm] arm center position';
arm.comment.step_min        = '[nm] minumum step size';
arm.comment.step            = '[nm] actual step size';
arm.comment.N               = 'number of steps';
arm.comment.min             = '[nm] min position of the arm';
arm.comment.max             = '[nm] max position of the arm';
arm.comment.vect            = '[nm] prefixed positons of the arm';
arm.comment.vect_read       = '[nm] actual positons of the arm, read by the attuator';
arm.comment.vect_c          = '[nm] actual positons of the arm, centered';

% prepare variables for processed images (grayscale)
in.comment.all              = '[a.u] intensity of interferogram as function of (x,y,delay)';