clear;//足裏点モデル
figure(2);
m1=0.2;
m2=1;
m3=0.2;
M=m1+m2+m3;
I1=0.03;
I2=0.2;
I3=0.03;
Sth=15/180*%pi;//股傾斜角度
L3=0.8;
L1=0.8;//L3*cos(Sth);
L2=0.4;
Lg=(L3*sin(Sth)+L2/tan(%pi/3))*tan(%pi/3);
l1g=0.4;
l2g=0.2;
l3g=0.2;
g=9.80665;
k=8;//トルク制御ゲイン
//th1=10/180*%pi;
//th2=10.40/180*%pi;
//th3=40.074/180*%pi;

//th1=18.125/180*%pi;
//th2=33.664/180*%pi;
//th3=0/180*%pi;

th1=18.125/180*%pi;
th2=33.664/180*%pi;
th3=0/180*%pi;

th1=18.125/180*%pi;
th2=33.664/180*%pi;
th3=0/180*%pi;

th12=th1+th2;
th123=th12+th3;
TH=[th1 th12 th123]'


thd1=0;
thd12=0;
thd123=0/180*%pi;
THd=[thd1 thd12 thd123]'

thdd1=0;
thdd12=0;
thdd123=0;
THdd=[thdd1 thdd12 thdd123]'

mc2 = m2*L1*l2g + m3*L1*L2;
mc3 = m3*L2*l3g;
mc23 = m3*L1*l3g;

M0 = m1*[l1g^2 0 0;0 0 0;0 0 0] + m2*[L1^2 0 0;0 l2g^2 0;0 0 0] + m3*[L1^2 0 0;0 L2^2 0;0 0 l3g^2] + [I1 0 0;0 I2 0;0 0 I3];
Mc2=-mc2*[0 1 0;1 0 0;0 0 0];
Mc3=mc3*[0 0 0;0 0 1;0 1 0];
Mc23=-mc23*[0 0 1;0 0 0;1 0 0];


tu1=0;
tu2=-k*th123;
tu3=0;


THdd=inv(M0 + Mc2*cos(th2) + Mc3*cos(th3) + Mc23*cos(th2+th3)) * ([-tu2 tu2 0]' + mc2*sin(th2)*[-(thd12)^2 thd1^2 0]' + mc3*sin(th3)*[0 (thd123)^2 -(thd12)^2]' + mc23*sin(th2+th3)*[-(thd123)^2 0 thd1^2]' - [-m1*l1g-m2*L1-m3*L1 0 0;0 m2*l2g+m3*L2 0;0 0 m3*l3g]*[sin(th1) sin(th12) sin(th123)]'*g);
 
 dt=0.0001;
 Tth1=[th1];
 Tth12=[th12];
 Tth123=[th123];
 T=0:dt:5
 
 //初期位置のリンク描写
 x1=L1*sin(th1);
 z1=L1*cos(th1);
 x2=L1*sin(th1)-L2*sin(th12);
 z2=L1*cos(th1)-L2*cos(th12);
 x3=L1*sin(th1)-L2*sin(th12)-L3*sin(th123);
 z3=L1*cos(th1)-L2*cos(th12)-L3*cos(th123);
 
 xprt=[0,x1];
 zprt=[0,z1];
 plot(xprt,zprt,'--');
 xprt=[x1,x2];
 zprt=[z1,z2];
 plot(xprt,zprt,'--');
 xprt=[x2,x3];
 zprt=[z2,z3];
 plot(xprt,zprt,'--');
 
//地面
xplt=[-2,2];
zplt=[0,0];
plot(xplt,zplt,);
 
 //初期位置の重心
 xg1=l1g*sin(th1);
 zg1=l1g*cos(th1);
 xg2=L1*sin(th1)-l2g*sin(th12);
 zg2=L1*cos(th1)-l2g*cos(th12);
 xg3=L1*sin(th1)-L2*sin(th12)-l3g*sin(th123);
 zg3=L1*cos(th1)-L2*cos(th12)-l3g*cos(th123);
 
x0=(m1*xg1+m2*xg2+m3*xg3)/(m1+m2+m3);
z0=(m1*zg1+m2*zg2+m3*zg3)/(m1+m2+m3);
xplt=[x0,x0];
zplt=[1+z0,-1+z0];
plot(xplt,zplt,':');
xplt=[1+x0,-1+x0];
zplt=[z0,z0];
plot(xplt,zplt,':');
 
 count=0;
 for t=dt:dt:0.6
     THd_a=THd;
     THd=THd+THdd*dt;
     TH=TH+(THd_a+THd)/2*dt
     
     th1=TH(1,1);
     th12=TH(2,1);
     th123=TH(3,1);
     th2=th12-th1;
     th3=th123-th12;
     thd1=THd(1,1);
     thd12=THd(2,1);
     thd123=THd(3,1);
     
     Tth1=[Tth1 th1];
     Tth12=[Tth12 th12];
     Tth123=[Tth123 th123];
     
     tu2=-k*th123;
     THdd=inv(M0 + Mc2*cos(th2) + Mc3*cos(th3) + Mc23*cos(th2+th3)) * ([-tu2 tu2 0]' + mc2*sin(th2)*[-(thd12)^2 thd1^2 0]' + mc3*sin(th3)*[0 (thd123)^2 -(thd12)^2]' + mc23*sin(th2+th3)*[-(thd123)^2 0 thd1^2]' - [-m1*l1g-m2*L1-m3*L1 0 0;0 m2*l2g+m3*L2 0;0 0 m3*l3g]*[sin(th1) sin(th12) sin(th123)]'*g);
     
     if modulo(t,0.05)==0 then
          //リンク描写
          count=count+1;

            x1=L1*sin(th1);
            z1=L1*cos(th1);
            x2=L1*sin(th1)-L2*sin(th12);
            z2=L1*cos(th1)-L2*cos(th12);
            x3=L1*sin(th1)-L2*sin(th12)-L3*sin(th123);
            z3=L1*cos(th1)-L2*cos(th12)-L3*cos(th123);
 
            xprt=[0,x1];
            zprt=[0,z1];
            plot2d(xprt,zprt,count);
            xprt=[x1,x2];
            zprt=[z1,z2];
            plot2d(xprt,zprt,count);
            xprt=[x2,x3];
            zprt=[z2,z3];
            plot2d(xprt,zprt,count);
             xg1=l1g*sin(th1);
             
            //重心描写
            zg1=l1g*cos(th1);
            xg2=L1*sin(th1)-l2g*sin(th12);
            zg2=L1*cos(th1)-l2g*cos(th12);
            xg3=L1*sin(th1)-L2*sin(th12)-l3g*sin(th123);
            zg3=L1*cos(th1)-L2*cos(th12)-l3g*cos(th123);
 
            x0=(m1*xg1+m2*xg2+m3*xg3)/M;
            z0=(m1*zg1+m2*zg2+m3*zg3)/M;
            xplt=[x0,x0];
            zplt=[0.01+z0,-0.01+z0];
            plot2d(xplt,zplt,count);
            xplt=[0.01+x0,-0.01+x0];
            zplt=[z0,z0];
            plot2d(xplt,zplt,count);

          
     
    end
end
 
 //plot(T,Tth1);
 //plot(T,Tth12);
 //plot(T,Tth123);
 
 
