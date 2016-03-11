%% import data
% User must load manually the next list of files:
% 20160307_longscan2.mat
% 20160307_sett2.mat
% 20160307_att2.mat
% 20160307_arm2
% 20160307_130125_SEAF_neon2p4bar_805.0nm_10,50mradpfs.mat

%% change units of measure from [mm,s] to [nm,fs]

att.center = att.center*1e6;
att.step_min = att.step_min *1e6;
att.step = att.step *1e6;
att.min = att.min *1e6;
att.max = att.max *1e6;
att.vect = att.vect *1e6;
att.vect_read = att.vect_read *1e6;

arm.step_min = arm.step_min *1e6;
arm.step = arm.step*1e6;
arm.min = arm.min *1e6;
arm.max = arm.max *1e6;
arm.vect = arm.vect *1e6;
arm.vect_read = arm.vect_read *1e6;
arm.vect_0 = arm.vect_0 *1e6;

%% rename data
in.all = im_all;
spi.fileformat_version  = fileformat_version;
spi.proc_settings       = proc_settings;
spi.raw                 = raw;
spi.rec                 = rec;
spi.save_timestamp      = save_timestamp;
spi.valid               = valid;

arm.vect_c              = arm.vect_0;
arm = rmfield(arm, 'vect_0');

clear im_all fileformat_version proc_settings raw rec save_timestamp valid arm.vect_0

%% resize data
in.all = in.all(143:343,224:424,:);
sett.x = sett.x(143:343);
sett.y = sett.y(224:424);
sett.processed = [201,201];