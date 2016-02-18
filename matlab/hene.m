%%
z=Zaber.Controller;

%%

wc=CameraViewerStd('WinVideo.Controller','ImageViewerStd',figure,'Webcam Viewer');
set(wc.processor.flatten_colour,'Value',1);


%%

z.move_absolute(2,22.208580, 'mm')


%% take image:
im=wc.camera.take_snapshot;

figureNF('image');clf
imagesc(mean(im,3))

%%
step_zaber = 0.047625/1000; % um/1000 = mm
pic_size   = size(wc.camera.take_snapshot);

pos.center = 23;
pos.step   = 2*step_zaber;
pos.n      = 150; % --> 2n+1
pos.min    = pos.center - pos.n*pos.step;
pos.max    = pos.center + pos.n*pos.step;

positions  = [pos.min:pos.step:pos.max];
N          = 2*pos.n + 1;
    
im_store   = zeros(pic_size(1),pic_size(2),pic_size(3),10);
im_all     = zeros(pic_size(1),pic_size(2),N);
im         = zeros(pic_size(1),pic_size(2));

% move attuator to the start
z.move_absolute(2,positions(1), 'mm');
pause(0.5);

for n=1:N
    pos_this = positions(n);
    z.move_absolute(2,pos_this, 'mm');

    for j=1:10
        im_store(:,:,:,j)=wc.camera.take_snapshot;
        pause(0.1);
    end
    im = mean(mean(im_store,4),3);
    im_all(:,:,n) = im;
    
end
