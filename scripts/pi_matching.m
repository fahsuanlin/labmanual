close all; clear all;

R=[5:100]; %ohm; loading resistance;
R0=50; %source impedance

f=300e6;% Hz; Larmor frequency

x=sqrt(R.*R0);
L=x./2./pi./f.*1e9; %nH
C=1./2./pi./f./x.*1e12; %pF


figure; hold on;

hh1=plotyy(R,L,R,C);
hh1c=get(hh1,'child');
etc_plotstyle

xlabel('loading resistance (Ohm)');
ylabel(hh1(1),'nH');
ylabel(hh1(2),'pF');

set(gcf,'pos',[1440         744         1000         600]);
hgexport(gcf,sprintf('pi_matching_7T.png'), hgexport('factorystyle'),'Format','png');



figure;
f=123.2e6;% Hz; Larmor frequency

x=sqrt(R.*R0);
L=x./2./pi./f.*1e9; %nH
C=1./2./pi./f./x.*1e12; %pF


hh2=plotyy(R,L,R,C);
hh2c=get(hh2,'child');
%set(hh2c{1},'color','b');
%set(hh2c{2},'color','r');
etc_plotstyle

xlabel('loading resistance (Ohm)');
ylabel(hh2(1),'nH');
ylabel(hh2(2),'pF');

set(gcf,'pos',[1440         744         1000         600]);
hgexport(gcf,sprintf('pi_matching_3T.png'), hgexport('factorystyle'),'Format','png');

return;

set(hh1(1),'ylim',[ 0 40])
set(hh1(2),'ylim',[ 10 190])
set(hh2(1),'ylim',[ 0 40])
set(hh2(2),'ylim',[ 10 190])
