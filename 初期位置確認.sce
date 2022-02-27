clear;//足裏円モデル
figure(51);
m3=0.1;
m2=1;
m1=m3*2;
M=m1+m2+m3;


Sth=75/180*%pi;//股傾斜角度
L3=0.8;
L1=L3*sin(Sth);
L2=0.2;
Lg=(L3*cos(Sth)+L2/tan(%pi/3))*tan(%pi/3);//初期位置での原点から足先端までの距離
l1g=0.4;
l2g=0.2;
l3g=0.4;
l1r=0.1;//足裏円形状の半径


//初期角度L2=0.4用
//th1=18.125/180*%pi;
//th2=33.664/180*%pi;
//th3=0/180*%pi;

//初期角度L2=0.2用
th1=6.037/180*%pi;
th2=33.747/180*%pi;
th3=0/180*%pi;

th12=th1+th2;
th123=th12+th3;
TH=[th1 th12 th123]'

thd1=0;
thd12=0;
thd123=0/180*%pi;

thdd1=0;
thdd12=0;
thdd123=0;
THdd=[thdd1 thdd12 thdd123]'
x0=[th1;th12;th123;0;0;0];


I1=m1*L1*L1/12;
I2=m2*L2*L2/12;
I3=m3*L3*L3/12;

g=9.80665;

xr=0;//足裏接点初期値


dt=0.001;

ke=1;
ku=0;
kud=0;


ztc=-3.557;//目標軌跡円の中心点（x=0とする）
x3=xr+(L1-l1r)*sin(th1)-L2*sin(th12)-L3*sin(th123);
z3=l1r+(L1-l1r)*cos(th1)-L2*cos(th12)-L3*cos(th123);
tr=( (z3-ztc)^2 + x3^2 )^0.5;//目標軌跡円の半径

//目標軌跡円を描く
hth=18/180*%pi;
dth=0.05/180*%pi;
for q=-(-0.25-dth):dth:hth
    xprt=[tr*cos(q-dth+%pi-%pi/2) , tr*cos(q+%pi-%pi/2)];
    zprt=[tr*sin(q-dth+%pi-%pi/2)+ztc , tr*sin(q+%pi-%pi/2)+ztc];
    plot(xprt,zprt,'--');
end

//地面
xplt=[-2,2];
zplt=[0,0];
plot(xplt,zplt,);



tu2=0;//トルクの初期値



//初期位置のリンク描写
 x1=xr+(L1-l1r)*sin(th1);
 x2=xr+(L1-l1r)*sin(th1)-L2*sin(th12);
 x3=xr+(L1-l1r)*sin(th1)-L2*sin(th12)-L3*sin(th123);
 z1=l1r+(L1-l1r)*cos(th1);
 z2=l1r+(L1-l1r)*cos(th1)-L2*cos(th12);
 z3=l1r+(L1-l1r)*cos(th1)-L2*cos(th12)-L3*cos(th123);
 
 //円弧とリンク１の交点を求める
 xonc=l1r*sin(th1-%pi)+xr;
 yonc=l1r*cos(th1-%pi)+l1r;

 xprt=[xonc,x1];
 zprt=[yonc,z1];
 plot(xprt,zprt,'b');
 xprt=[x1,x2];
 zprt=[z1,z2];
 plot(xprt,zprt,'b');
 xprt=[x2,x3];
 zprt=[z2,z3];
 plot(xprt,zprt,'b');
 
 
//リンク1を基準に円弧を描く
hth=30/180*%pi;
dth=1/60*%pi;
for q=-(hth-dth):dth:hth
    xprt=[l1r*cos(q-dth-th1-%pi/2)+xr , l1r*cos(q-th1-%pi/2)+xr];
    zprt=[l1r*sin(q-dth-th1-%pi/2)+l1r , l1r*sin(q-th1-%pi/2)+l1r];
    plot(xprt,zprt,'--');
end

 /*
 //初期位置の重心
 xg1=xr+(l1g-l1r)*sin(th1);
 xg2=xr+(L1-l1r)*sin(th1)-l2g*sin(th12);
 xg3=xr+(L1-l1r)*sin(th1)-L2*sin(th12)-l3g*sin(th123);
 zg1=l1r+(l1g-l1r)*cos(th1);
 zg2=l1r+(L1-l1r)*cos(th1)-l2g*cos(th12);
 zg3=l1r+(L1-l1r)*cos(th1)-L2*cos(th12)-l3g*cos(th123);
 
x0=(m1*xg1+m2*xg2+m3*xg3)/M;
z0=(m1*zg1+m2*zg2+m3*zg3)/M;
xplt=[x0,x0];
zplt=[1+z0,-1+z0];
plot(xplt,zplt,':');
xplt=[1+x0,-1+x0];
zplt=[z0,z0];
plot(xplt,zplt,':');
