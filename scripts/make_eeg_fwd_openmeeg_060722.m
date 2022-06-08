close all; clear all;

fstem='180322_PYW_d10_wb_eeg_060722';


%#OPENMEEG routines
cmd=sprintf('!om_assemble -HeadMat %s_bem.geom %s_bem.cond %s.hm',fstem,fstem,fstem);
eval(cmd);

cmd=sprintf('!om_assemble -DSM %s_bem.geom %s_bem.cond %s-lh.dip %s-lh.dsm brain',fstem,fstem,fstem,fstem);
eval(cmd);

cmd=sprintf('!om_assemble -DSM %s_bem.geom %s_bem.cond %s-rh.dip %s-rh.dsm brain',fstem,fstem,fstem,fstem);
eval(cmd);

cmd=sprintf('!om_assemble -h2em %s_bem.geom %s_bem.cond %s.eegsensors %s.h2em',fstem,fstem,fstem,fstem);
eval(cmd);

cmd=sprintf('!om_minverser %s.hm %s.hm_inv', fstem,fstem);
eval(cmd);

cmd=sprintf('!om_gain -EEG %s.hm_inv %s-lh.dsm %s.h2em %s-lh.gain.mat',fstem,fstem,fstem,fstem);
eval(cmd);

cmd=sprintf('!om_gain -EEG %s.hm_inv %s-rh.dsm %s.h2em %s-rh.gain.mat',fstem,fstem,fstem,fstem);
eval(cmd);
                                       
