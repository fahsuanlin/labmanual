close all; clear all;

subject='s033';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cmd=sprintf('!om_assemble -HeadMat %s_d10_wb_dec_072521_bem.geom %s_d10_wb_dec_072521_bem.cond %s_d10_wb_dec_072521.hm',subject, subject, subject);
eval(cmd);

cmd=sprintf('!om_assemble -DSM %s_d10_wb_dec_072521_bem.geom %s_d10_wb_dec_072521_bem.cond %s_d10_wb_dec_072521-lh.dip %s_d10_wb_dec_072521-lh.dsm brain',subject, subject, subject, subject);
eval(cmd);

cmd=sprintf('!om_assemble -DSM %s_d10_wb_dec_072521_bem.geom %s_d10_wb_dec_072521_bem.cond %s_d10_wb_dec_072521-rh.dip %s_d10_wb_dec_072521-rh.dsm brain',subject, subject, subject, subject);
eval(cmd);

cmd=sprintf('!om_assemble -h2ipm %s_d10_wb_dec_072521_bem.geom %s_d10_wb_dec_072521_bem.cond %s_d10_wb_dec_072521.seeg  %s_d10_wb_dec_072521.h2ipm',subject, subject, subject, subject);
eval(cmd);

cmd=sprintf('!om_assemble -ds2ipm %s_d10_wb_dec_072521_bem.geom %s_d10_wb_dec_072521_bem.cond %s_d10_wb_dec_072521-lh.dip %s_d10_wb_dec_072521.seeg seeg_wb_072521-lh.ds2ipm',subject, subject, subject, subject);
eval(cmd);

cmd=sprintf('!om_assemble -ds2ipm %s_d10_wb_dec_072521_bem.geom %s_d10_wb_dec_072521_bem.cond %s_d10_wb_dec_072521-rh.dip %s_d10_wb_dec_072521.seeg seeg_wb_072521-rh.ds2ipm',subject, subject, subject, subject);
eval(cmd);

cmd=sprintf('!om_minverser %s_d10_wb_dec_072521.hm %s_d10_wb_dec_072521.hm_inv',subject, subject);
eval(cmd);

cmd=sprintf('!om_gain -InternalPotential %s_d10_wb_dec_072521.hm_inv %s_d10_wb_dec_072521-lh.dsm %s_d10_wb_dec_072521.h2ipm seeg_wb_072521-lh.ds2ipm %s_d10_wb_dec_072521-lh.gain.mat',subject, subject, subject, subject);
eval(cmd);

cmd=sprintf('!om_gain -InternalPotential %s_d10_wb_dec_072521.hm_inv %s_d10_wb_dec_072521-rh.dsm %s_d10_wb_dec_072521.h2ipm seeg_wb_072521-rh.ds2ipm %s_d10_wb_dec_072521-rh.gain.mat',subject, subject, subject, subject);
eval(cmd);


