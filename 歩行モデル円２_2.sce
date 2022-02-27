clear;//足裏円モデル
figure(5);

m1=0.2;
m2=1;
m3=0.2;
M=m1+m2+m3;

Sth=15/180*%pi;//股傾斜角度
L3=0.8;
L1=0.781;//L3*cos(Sth);
L2=0.4;
Lg=(L3*sin(Sth)+L2/tan(%pi/3))*tan(%pi/3);
l1g=0.4;
l2g=0.2;
l3g=0.4;
l1r=0.1;//足裏円形状の半径

I1=m1*L1*L1/12;
I2=m2*L2*L2/12;
I3=m3*L3*L3/12;

g=9.80665;
k=0;//トルク制御ゲイン
k2=-20;

xr=0;//足裏接点初期値
xrd=0;//足裏接点初期速度

//ztc=-2.393;//目標軌跡円の中心点、（x=0とする）
ztc=-2.17
//tr=2.493;//目標軌跡円の半径
tr=2.3;


//目標軌跡円を描く
hth=30/180*%pi;
dth=5/180*%pi;
for q=-(hth-dth):dth:hth
    xprt=[tr*cos(q-dth+%pi-%pi/2) , tr*cos(q+%pi-%pi/2)];
    zprt=[tr*sin(q-dth+%pi-%pi/2)+ztc , tr*sin(q+%pi-%pi/2)+ztc];
    plot(xprt,zprt,'--');
end


//初期角度
//th1=10/180*%pi;
//th2=10.40/180*%pi;
//th3=40.074/180*%pi;

//th1=18.125/180*%pi;
//th2=33.664/180*%pi;
//th3=0/180*%pi;

//th1=18.125/180*%pi;
//th2=33.664/180*%pi;
//th3=0/180*%pi;

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
 plot(xprt,zprt,'--');
 xprt=[x1,x2];
 zprt=[z1,z2];
 plot(xprt,zprt,'--');
 xprt=[x2,x3];
 zprt=[z2,z3];
 plot(xprt,zprt,'--');
 
 
//リンク1を基準に円弧を描く
hth=30/180*%pi;
dth=1/180*%pi;
for q=-(hth-dth):dth:hth
    xprt=[l1r*cos(q-dth-th1-%pi/2)+xr , l1r*cos(q-th1-%pi/2)+xr];
    zprt=[l1r*sin(q-dth-th1-%pi/2)+l1r , l1r*sin(q-th1-%pi/2)+l1r];
    plot(xprt,zprt,'--');
end

 
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
 

 
//地面
xplt=[-2,2];
zplt=[0,0];
plot(xplt,zplt,);




//トルク制御(初期トルク)----------------------------------------------------------------------------------------------------------------------------------
tu1=0;
e=( ( (z3-ztc)^2 + x3^2 )^0.5 -tr)
ed=0;
tu2=k*e-k2*ed;
e_a=e;
tu3=0;

//thdd1,thdd12,thdd123を求める
M3=[L1^2 -L1*L2 -L1*l3g;-L1*L2 L2^2 L2*l3g;-L1*l3g L2*l3g l3g^2] + [1 0 0;0 0 0;0 0 0]*2*l1r*(L1-l1r)*(cos(th1)-1) + [0 1 0;1 0 0;0 0 0]*L2*(-L1*(cos(th2)-1)-l1r*(cos(th12)-cos(th2))) + [0 0 1;0 0 0;1 0 0]*l3g*(-L1*(cos(th2+th3)-1) - l1r*(cos(th123)-cos(th2+th3))) + [0 0 0;0 0 1;0 1 0]*L2*l3g*(cos(th3)-1);

M2=[L1^2 -L1*l2g 0;-L1*l2g l2g^2 0;0 0 0] + [1 0 0;0 0 0;0 0 0]*2*l1r*(L1-l1r)*(cos(th1)-1) + [0 1 0;1 0 0;0 0 0]*l2g*(-L1*(cos(th2)-1)-l1r*(cos(th12)-cos(th2)));

M1=[l1g^2 0 0;0 0 0;0 0 0] + [1 0 0;0 0 0;0 0 0]*2*l1r*(l1g-l1r)*(cos(th1)-1);

NL= [-(m1*l1g+m2*L1+m3*L1-M*l1r)*l1r*sin(th1);-(m2*l2g + m3*L2)*(L1-l1r)*sin(th2);-m3*l3g*(L1-l1r)*sin(th2+th3)]*thd1^2 + [(m2*l2g+m3*L2)*((L1-l1r)*sin(th2)+l1r*sin(th12));0;m3*L2*l3g*sin(th3)]*thd12^2 + [m3*l3g*((L1-l1r)*sin(th2+th3)+l1r*sin(th123));-m3*L2*l3g*sin(th3);0]*thd123^2;

THdd=inv(m1*M1+m2*M2+m3*M3+[I1 0 0;0 I2 0;0 0 I3])*([-tu2;tu2;0]-[-(m1*l1g+m2*L1+m3*L1-M*l1r) 0 0;0 m2*l2g+m3*L2 0;0 0 m3*l3g]*[sin(th1);sin(th12);sin(th123)]*g-NL);

 dt=0.0001;
 Tth1=[th1];
 Tth12=[th12];
 Tth123=[th123];
 T=0:dt:5

 count=0;
 for t=dt:dt:0.6
     //前の角加速度から現在の角速度，角度を求める
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
     
     xrd_a=xrd;//前ステップの足裏接点速度を保存
     xrd=l1r*thd1;//リンク１の角速度から足裏接点の速度を求める
     xr=xr+(xrd+xrd_a)/2*dt;//ステップ間の速度を現在と前ステップの速度の平均と考え足裏接点位置を求める
     
     
     //Tth1=[Tth1 th1];
     //Tth12=[Tth12 th12];
     //Tth123=[Tth123 th123];
     
     e=( ( (z3-ztc)^2 + x3^2 )^0.5 -tr)
     ed=(e-e_a)/dt;
     tu2=k*e-k2*ed;
     e_a=e;
     
     //現在の状態からトルクを修正---------------------------------------------------------------------------------------
     

     
     if modulo(t,0.05)==0 then
          //リンク描写
            count=count+1;//色設定

            //リンク描写
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
            plot2d(xprt,zprt,count);
            xprt=[x1,x2];
            zprt=[z1,z2];
            plot2d(xprt,zprt,count);
            xprt=[x2,x3];
            zprt=[z2,z3];
            plot2d(xprt,zprt,count);
            
            //リンク1を基準に円弧を描く
            hth=30/180*%pi;
            dth=5/180*%pi;
            for q=-(hth-dth):dth:hth
                xprt=[l1r*cos(q-dth-th1-%pi/2)+xr,l1r*cos(q-th1-%pi/2)+xr];
                zprt=[l1r*sin(q-dth-th1-%pi/2)+l1r,l1r*sin(q-th1-%pi/2)+l1r];
                plot2d(xprt,zprt,count);
            end
            
            //重心描写
            
            xg1=xr+(l1g-l1r)*sin(th1);
            xg2=xr+(L1-l1r)*sin(th1)-l2g*sin(th12);
            xg3=xr+(L1-l1r)*sin(th1)-L2*sin(th12)-l3g*sin(th123);
            zg1=l1r+(l1g-l1r)*cos(th1);
            zg2=l1r+(L1-l1r)*cos(th1)-l2g*cos(th12);
            zg3=l1r+(L1-l1r)*cos(th1)-L2*cos(th12)-l3g*cos(th123);
 
            x0=(m1*xg1+m2*xg2+m3*xg3)/M;
            z0=(m1*zg1+m2*zg2+m3*zg3)/M;
            xplt=[x0,x0];
            zplt=[0.01+z0,-0.01+z0];
            plot2d(xplt,zplt,count);
            xplt=[0.01+x0,-0.01+x0];
            zplt=[z0,z0];
            plot2d(xplt,zplt,count);
            
            

    end
                //現在の状態からthdd1,thdd12,thdd123を求める
    M3=[L1^2 -L1*L2 -L1*l3g;-L1*L2 L2^2 L2*l3g;-L1*l3g L2*l3g l3g^2] + [1 0 0;0 0 0;0 0 0]*2*l1r*(L1-l1r)*(cos(th1)-1) + [0 1 0;1 0 0;0 0 0]*L2*(-L1*(cos(th2)-1)-l1r*(cos(th12)-cos(th2))) + [0 0 1;0 0 0;1 0 0]*l3g*(-L1*(cos(th2+th3)-1) - l1r*(cos(th123)-cos(th2+th3))) + [0 0 0;0 0 1;0 1 0]*L2*l3g*(cos(th3)-1);

    M2=[L1^2 -L1*l2g 0;-L1*l2g l2g^2 0;0 0 0] + [1 0 0;0 0 0;0 0 0]*2*l1r*(L1-l1r)*(cos(th1)-1) + [0 1 0;1 0 0;0 0 0]*l2g*(-L1*(cos(th2)-1)-l1r*(cos(th12)-cos(th2)));

    M1=[l1g^2 0 0;0 0 0;0 0 0] + [1 0 0;0 0 0;0 0 0]*2*l1r*(l1g-l1r)*(cos(th1)-1);

    NL= [-(m1*l1g+m2*L1+m3*L1-M*l1r)*l1r*sin(th1);-(m2*l2g + m3*L2)*(L1-l1r)*sin(th2);-m3*l3g*(L1-l1r)*sin(th2+th3)]*thd1^2 + [(m2*l2g+m3*L2)*((L1-l1r)*sin(th2)+l1r*sin(th12));0;m3*L2*l3g*sin(th3)]*thd12^2 + [m3*l3g*((L1-l1r)*sin(th2+th3)+l1r*sin(th123));-m3*L2*l3g*sin(th3);0]*thd123^2;

    THdd=inv(m1*M1+m2*M2+m3*M3+[I1 0 0;0 I2 0;0 0 I3])*([-tu2;tu2;0]-[-(m1*l1g+m2*L1+m3*L1-M*l1r) 0 0;0 m2*l2g+m3*L2 0;0 0 m3*l3g]*[sin(th1);sin(th12);sin(th123)]*g-NL);
    

end
 //ca=gca();
 //ca.data_bounds=[-0.8,-0.2;0.8,1.4];
 //plot(T,Tth1);
 //plot(T,Tth12);
 //plot(T,Tth123);
//square(1,1,1,1);
