clear ;

m1=0.4 ; m2=1.0 ; m3=0.08 ; m4=0.12 ;
L1=0.781 ; L2=0.2 ; L3=0.4 ; L4=0.4 ;
I1=m1*L1*L1/12 ; I2=m2*L2*L2/12 ; I3=m3*L3*L3/12 ; I4=m4*L4*L4/12
L1G=L1/2 ; L2G=L2/2 ; L3G=L3/2 ; L4G=L4/2
L1R=L1*0.1 ;
g=9.81 ;

q10=0*%pi/180;
q20=5*%pi/180;
q30=45*%pi/180;
q40=45*%pi/180;

x10=(L1G-L1R)*sin(q10);
x20=(L1 -L1R)*sin(q10)-L2G*sin(q20);
x30=(L1 -L1R)*sin(q10)-L2 *sin(q20)-L3G*sin(q30);
x40=(L1 -L1R)*sin(q10)-L2 *sin(q20)-L3 *sin(q30)-L4G*sin(q40);

x0=(m1*x10+m2*x20+m3*x30+m4*x40)/(m1+m2+m3+m4);
disp("x0=",x0);

function xdot=four_link(t,x,para)   // 微分方程式を解くための関数
    q1=x(1);
    q2=x(2);
    q3=x(3);
    q4=x(4);
    q1d=x(5);
    q2d=x(6);
    q3d=x(7);
    q4d=x(8);
    
    tau2 = para(1)*t^3 + para(2)*t^2 + para(3)*t ;
    tau3 = para(4)*t^3 + para(5)*t^2 + para(6)*t ;
    tau4 = para(7)*t^3 + para(8)*t^2 + para(9)*t ;
        
    M411=(L1-L1R)^2+L1R^2+2*L1R*(L1-L1R)* cos(q1);
    M412=-(L1-L1R)*L2 *cos(q2-q1)-L1R*L2* cos(q2);
    M413=-(L1-L1R)*L3 *cos(q3-q1)-L1R*L3* cos(q3);
    M414=-(L1-L1R)*L4G*cos(q4-q1)-L1R*L4G*cos(q4);
    M422=L2^2 ;
    M423=L2*L3* cos(q3-q2);
    M424=L2*L4G*cos(q4-q2);
    M433=L3^2 ;
    M434=L3*L4G*cos(q4-q3);
    M444=L4G^2;
    M4=[M411 M412 M413 M414;M412 M422 M423 M424;M413 M423 M433 M434;M414 M424 M434 M444];
    
    M311=(L1-L1R)^2+L1R^2+2*L1R*(L1-L1R)*cos(q1);
    M312=-(L1-L1R)*L2 *cos(q2-q1)-L1R*L2 *cos(q2);
    M313=-(L1-L1R)*L3G*cos(q3-q1)-L1R*L3G*cos(q3);
    M322=L2^2 ;
    M323=L2*L3G*cos(q3-q2);
    M333=L3^2 ;
    M3=[M311 M312 M313 0;M312 M322 M323 0;M313 M323 M333 0;0 0 0 0]; 

    M211=(L1-L1R)^2+L1R^2+2*L1R*(L1-L1R)*cos(q1);
    M212=-(L1-L1R)*L2G*cos(q2-q1)-L1R*L2G*cos(q2);
    M222=L2G^2 ; ;
    M2=[M211 M212 0 0;M212 M222 0 0;0 0 0 0;0 0 0 0]; 

    M111=(L1G-L1R)^2+L1R^2+2*L1R*(L1G-L1R)*cos(q1);
    M1=[M111 0 0 0;0 0 0 0;0 0 0 0;0 0 0 0]; 

    I=[I1 0 0 0;0 I2 0 0;0 0 I3 0;0 0 0 I4];
    
    M=m1*M1+m2*M2+m3*M3+I ;

    h1=( m1*(L1R-L1G)+m2*(L1R-L1)+m3*(L1R-L1)+m4*(L1R-L1) )*L1R*sin(q1)*q1d^2 ;
    h1=h1+(m2*L2G+m3*L2+m4*L2)*( (L1-L1R)*sin(q2-q1)+L1R*sin(q2) )*q2d^2 ;
    h1=h1+(m3*L3G+m4*L3)*( (L1-L1R)*sin(q3-q1)+L1R*sin(q3) )*q3d^2 ;
    h1=h1+m4*L4G*( (L1-L1R)*sin(q4-q1)+L1R*sin(q4) )*q4d^2 ;
    
    h2=(m2*L2G+m3*L2+m4*L2)*(L1R-L1)*sin(q2-q1)*q1d^2 ;
    h2=h2-(m3*L3G+m4*L3)*L2*sin(q3-q2)*q3d^2 ;
    h2=h2-m4*L4G*L2*sin(q4-q2)*q4d^2 ;
    
    h3=(m3*L3G+m4*L3)*(L1R-L1)*sin(q3-q1)*q1d^2 ;
    h3=h3+(m3*L3G+m4*L3)*L2*sin(q3-q2)*q2d^2 ;
    h3=h3-m4*L4G*L3*sin(q4-q3)*q4d^2 ;
    
    h4=m4*(L1R-L1)*L4G*sin(q4-q1)*q1d^2;
    h4=h4+m4*L2*L4G*sin(q4-q2)*q2d^2 ;
    h4=h4+m4*L3*L4G*sin(q4-q3)*q3d^2 ;
    
    h=[h1;h2;h3;h4];

    G1=-( m1*(L1G-L1R)+m2*(L1-L1R)+m3*(L1-L1R)+m4*(L1-L1R) )*sin(q1)*g ;
    G2=(m2*L2G+m3*L2+m4*L2)*sin(q2)*g ;
    G3=(m3*L3G+m4*L3)*sin(q3)*g ;
    G4=m4*L4G*sin(q4)*g;
    G=[G1;G2;G3;G4];
   
    tau=[ -tau2 ; tau2-tau3 ; tau3-tau4 ; tau4 ];
    tmp=inv( M )*( tau - h - G );
    
    xdot(1)=q1d;
    xdot(2)=q2d;
    xdot(3)=q3d;
    xdot(4)=q4d;
    xdot(5)=tmp(1);
    xdot(6)=tmp(2);
    xdot(7)=tmp(3);
    xdot(8)=tmp(4);
endfunction

//x0=[q10;q20;q30;q40;0;0;0;0];
//t=0:0.001:0.6;              // 評価時間　1[ms]ごと
//t0=0;                       // 時刻の初期値
//x=ode(x0,t0,t,four_link);   // 微分方程式のソルバ

function f=cost_f(a) 
    x0=[q10;q20;q30;q40;0;0;0;0];
    t=0:0.001:0.6;              // 評価時間　1[ms]ごと
    t0=0;                       // 時刻の初期値
    x=ode(x0,t0,t,list(four_link,a));   // 微分方程式のソルバ
    
    f=0;
    idx=1:601;				// プロット 25[ms]=0.025[s]ごと
    xq1=x(1,idx)';
    xq2=x(2,idx)';
    xq3=x(3,idx)';
    xq4=x(4,idx)';
    tt=t(1,idx)'
    
    x0a=L1R*xq1;
    x00a=x0a-L1R*sin(xq1);
    x1a=L1R*xq1+(L1-L1R)*sin(xq1);
    x2a=L1R*xq1+(L1-L1R)*sin(xq1)-L2*sin(xq2);
    x3a=L1R*xq1+(L1-L1R)*sin(xq1)-L2*sin(xq2)-L3*sin(xq3);
    x4a=L1R*xq1+(L1-L1R)*sin(xq1)-L2*sin(xq2)-L3*sin(xq3)-L4*sin(xq4);

    z0a=L1R*ones(tt);
    z00a=z0a-L1R*cos(xq1);
    z1a=L1R+(L1-L1R)*cos(xq1);
    z2a=L1R+(L1-L1R)*cos(xq1)-L2*cos(xq2);
    z3a=L1R+(L1-L1R)*cos(xq1)-L2*cos(xq2)-L3*cos(xq3);
    z4a=L1R+(L1-L1R)*cos(xq1)-L2*cos(xq2)-L3*cos(xq3)-L4*cos(xq4);
    
    f1=0.001*(z4a-0.1)'*(z4a-0.1);
    f2=(x4a(601)-0.6)^2 + z4a(601)^2 ;
    f=f1+0.1*f2;
endfunction

function [f, g, ind]=cost_fg(x, ind)
    f = cost_f(x)
    g = numderivative(cost_f, x);
    mprintf("f = %f g = %f \n",f,g);
endfunction

A0=zeros(1,9);
[fout,Aopt] = optim(cost_fg, A0);

x0=[q10;q20;q30;q40;0;0;0;0];
t=0:0.001:0.6;              // 評価時間　1[ms]ごと
t0=0;                       // 時刻の初期値
x=ode(x0,t0,t,list(four_link,Aopt));   // 微分方程式のソルバ

idx=1:25:601;				// プロット 25[ms]=0.025[s]ごと
xq1=x(1,idx)';
xq2=x(2,idx)';
xq3=x(3,idx)';
xq4=x(4,idx)';
tt=t(1,idx)'

x0a=L1R*xq1;
x00a=x0a-L1R*sin(xq1);
x1a=L1R*xq1+(L1-L1R)*sin(xq1);
x2a=L1R*xq1+(L1-L1R)*sin(xq1)-L2*sin(xq2);
x3a=L1R*xq1+(L1-L1R)*sin(xq1)-L2*sin(xq2)-L3*sin(xq3);
x4a=L1R*xq1+(L1-L1R)*sin(xq1)-L2*sin(xq2)-L3*sin(xq3)-L4*sin(xq4);

z0a=L1R*ones(tt);
z00a=z0a-L1R*cos(xq1);
z1a=L1R+(L1-L1R)*cos(xq1);
z2a=L1R+(L1-L1R)*cos(xq1)-L2*cos(xq2);
z3a=L1R+(L1-L1R)*cos(xq1)-L2*cos(xq2)-L3*cos(xq3);
z4a=L1R+(L1-L1R)*cos(xq1)-L2*cos(xq2)-L3*cos(xq3)-L4*cos(xq4);

subplot(211)
sL=0.5;
for i=1:length(tt)
    xx=sL*(i-1)+[x00a(i) x0a(i) x1a(i) x2a(i) x3a(i) x4a(i)];
    zz=[z00a(i) z0a(i) z1a(i) z2a(i) z3a(i) z4a(i)]';
    plot2d(xx,zz);
    plot2d(sL*(i-1),0,-1);
    plot2d(sL*(i-1)+x0a(i),z0a(i),-2);
    plot2d(sL*(i-1)+x1a(i),z1a(i),-3);
    plot2d(sL*(i-1)+x2a(i),z2a(i),-4);
    plot2d(sL*(i-1)+x3a(i),z3a(i),-5);
    plot2d(sL*(i-1)+x4a(i),z4a(i),-6);
end
xgrid
ca=gca();
ca.data_bounds=[-1 -0.2 ; 13 0.9 ] ;

subplot(223)
plot2d(t',180/%pi*x(1,:)',1);
plot2d(t',180/%pi*x(2,:)',2);
plot2d(t',180/%pi*x(3,:)',3);
plot2d(t',180/%pi*x(4,:)',4);
legend(['q1';'q2';'q3';'q4']);
xgrid

subplot(224)
ttau2 = Aopt(1)*t.^3 + Aopt(2)*t.^2 + Aopt(3)*t ;
ttau3 = Aopt(4)*t.^3 + Aopt(5)*t.^2 + Aopt(6)*t ;
ttau4 = Aopt(7)*t.^3 + Aopt(8)*t.^2 + Aopt(9)*t ;
plot2d(t',ttau2',1);
plot2d(t',ttau3',2);
plot2d(t',ttau4',3);
legend(['tau2';'tau3';'tau4']);
xgrid
