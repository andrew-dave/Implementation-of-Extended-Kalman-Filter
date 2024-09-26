function [fxun,Ati,Uti] = assump_linear(Gi,Ri)

syms Px Py Pz r p y vx vy vz wx wy wz ax ay az ngx ngy ngz nax nay naz bgx bgy bgz bax bay baz nbgx nbgy nbgz nbax nbay nbaz;
pos = [Px Py Pz];
ori = [r p y];
lvel = [vx vy vz];
avel = [wx wy wz];
accl = [ax ay az];
ng = [ngx ngy ngz];
na = [nax nay naz];
bg = [bgx bgy bgz];
ba = [bax bay baz];
nbg = [nbgx nbgy nbgz];
nba = [nbax nbay nbaz];
%F=vertcat(vx,vy,vz,inv(G_sym)*[wx-ngx-nbgx;wy-ngy-nbgy;wz-ngz-nbgz],[0;0;-9.81]+R_sym*[ax-nax-nbax;ay-nay-nbay;az-naz-nbaz],bgx,bgy,bgz,bax,bay,baz);
fxun=vertcat(lvel.',inv(Gi)*((avel-ng-bg).'),[0;0;-9.81]+Ri*((accl-na-ba).'),nbg.',nba.');
% J=jacobian(F);
Ati=jacobian(fxun,[pos,ori,lvel,bg,ba]);
Uti=jacobian(fxun,[ng,na]);

