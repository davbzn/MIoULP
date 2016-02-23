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

%% delayscan
PDIscan=[];
PDIscan.start = 20.0;
PDIscan.stepsize_mm = 632.8/20 *1e-9*1e3; % lambda/10 devided by two -> delaystage
PDIscan.N = 60;
PDIscan.stop = PDIscan.start+PDIscan.stepsize_mm*PDIscan.N;
PDIscan.cycles = (PDIscan.stop-PDIscan.start)/(630/2*1e-9*1e3);
PDIscan.imperstep = 5;
PDIscan

%%
z.move_absolute(2, PDIscan.start-0.2, 'mm');
z.move_absolute(2, PDIscan.start-0.1, 'mm');
z.move_absolute(2, PDIscan.start, 'mm');
PDIscan.stagepos_mm          = zeros(PDIscan.N,1);
PDIscan.stagepos_mm_reported = zeros(PDIscan.N,1);
h = awaitbar(0,'Running scan, please wait...');
for n=1:PDIscan.N
    z.move_relative(2, PDIscan.stepsize_mm, 'mm');
    pause(0.1);
    PDIscan.stagepos_mm(n) = PDIscan.start+n*PDIscan.stepsize_mm;
    %PDIscan.stagepos_mm_reported(n) = z.get_position(2)*1e-15/2*SI('c')*1e3;
    PDIscan.stagepos_mm_reported(n) = z.get_position(2,'mm');
    for m=1:PDIscan.imperstep
        tmp(:,:,:,m)=wc.camera.take_snapshot;
    end
    im = mean(mean(tmp,3),4);
    if n==1
        PDIscan.I = zeros([size(im),PDIscan.N]);
        [PDIscan.x,PDIscan.y]=makexyvectors(PDIscan.I(:,:,1));
    end
    PDIscan.I(:,:,n) = im;
    abort = awaitbar(n/PDIscan.N,h,sprintf('step %i of %i',n,PDIscan.N));  % with update message
    if abort; break; end
end

%% save:
filename = [datestr(now,'yyyymmdd_HHMM'),'_PDIscan']
save(filename','PDIscan')

%%
pfit=polyfit(PDIscan.stagepos_mm,PDIscan.stagepos_mm_reported,1)
figureNF('stagepos');clf
plot(PDIscan.stagepos_mm, PDIscan.stagepos_mm_reported, '.',...
    PDIscan.stagepos_mm, polyval(pfit,PDIscan.stagepos_mm),'-')
title(sprintf('slope %.3f',pfit(1)))
xlabel('set position (mm)')
ylabel('reported position (mm)')





%%
mpos=maxpos2d_by_moment(PDIscan.x,PDIscan.y,mean(PDIscan.I,3))
x0=mpos(1);
y0=mpos(2);
x0=173
y0=314
figureNF('scan');clf
him=imagesc(PDIscan.x,PDIscan.y,PDIscan.I(:,:,3));
line(min_max(PDIscan.x),[1,1]*y0)
line([1,1]*x0,min_max(PDIscan.y))
htit=title('scan');
%%
for n=1:PDIscan.N
    set(him,'CData',PDIscan.I(:,:,n))
    set(htit,'String',sprintf('step %i of %i',n,PDIscan.N))
    pause(0.05)
    drawnow
end

%%
lineout = makecol(squeeze(PDIscan.I(find_index(PDIscan.y,y0),find_index(PDIscan.x,x0),:)));
lineout=lineout-min(lineout);
if 0
    delay_fs = makecol(2*(PDIscan.stagepos_mm-PDIscan.stagepos_mm(1))*1e-3./SI('c')*1e15)
    delay_fs = delay_fs;
else
    delay_fs = makecol(2*(PDIscan.stagepos_mm_reported-PDIscan.stagepos_mm_reported(1))*1e-3./SI('c')*1e15)
end
T0 = 2*pi/chng_rng(632.8)
cycles = delay_fs/T0
%wind = mygauss(delay_fs-mean(delay_fs),15,4);
figureNF('lineout');clf
%plot(delay_fs, lineout, '.-', delay_fs, lineout.*wind, 'g-')
plot(cycles, lineout, '.-')%, cycles, lineout.*wind, 'g-')
xlabel('cycles / T_0^{HeNe}')
%
omega_0 = chng_rng(632.8)
omega = dftConjAxis(delay_fs);
lineout_FFT = dftv(delay_fs,1,lineout,omega);

figureNF('lineoutFFT');clf
plot(omega./omega_0,abs(lineout_FFT),'.-')
xlim([0 5])
xlabel('\omega/omega_0^{HeNe}')

%%
scan_FFT1=PDIscan.I.*0;
for n=1:size(PDIscan.I,1)
    for m=1:size(PDIscan.I,2)
       scan_FFT1(n,m,:) = dftv(delay_fs,1,makecol(PDIscan.I(n,m,:)),omega);
    end
end
%%
figureNF('scan_FFT');clf
plot(omega./omega_0,abs(squeeze(scan_FFT1(find_index(PDIscan.y,y0),find_index(PDIscan.x,x0),:) )))


%% 2D stuff:
%PDIscan.I2 = reshape(PDIscan.I,size(PDIscan.I,3),size(PDIscan.I,1),size(PDIscan.I,2));
omega_0 = chng_rng(632.8)
omega = dftConjAxis(delay_fs);
scan_FFT = fftshift(fft(PDIscan.I,[],3));
%scan_FFT = fftshift(fft(PDIscan.I2,[],1));

filter = mygauss(omega-1,0.5,4);

%scan_FFT_filt = bsxfun(@times, scan_FFT, filter);
%scan_FFT_filt = btimes(scan_FFT, filter');


figureNF('scan_FFT');clf
plot((omega-mean(diff(omega)))./omega_0,abs(squeeze(scan_FFT(find_index(PDIscan.y,y0),find_index(PDIscan.x,x0),:))))
%plot((omega-mean(diff(omega)))./omega_0,abs(squeeze(scan_FFT(:,find_index(PDIscan.y,y0),find_index(PDIscan.x,x0)))))

%%
PDIscan.I2 = reshape(PDIscan.I,size(PDIscan.I,3),size(PDIscan.I,1),size(PDIscan.I,2));
scan_FFT2 = dftv(delay_fs,1,PDIscan.I2,omega,1);
figureNF('scan_FFT2');clf
plot(omega./omega_0,abs(squeeze(scan_FFT2(:,find_index(PDIscan.y,y0),find_index(PDIscan.x,x0)))))
