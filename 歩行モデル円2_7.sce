clear;//足裏円モデル
figure(6);
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

xr=0;//足裏接点初期値

dt=0.001;

ke=1;
ku=0;
kud=0;

//ztc=-2.17;//ztc=-2.393;//目標軌跡円の中心点、（x=0とする）
//tr=2.3;//tr=2.493;//目標軌跡円の半径

ztc=-2.37
tr=2.493;//目標軌跡円の半径


//目標軌跡円を描く
hth=30/180*%pi;
dth=5/180*%pi;
for q=-(hth-dth):dth:hth
    xprt=[tr*cos(q-dth+%pi-%pi/2) , tr*cos(q+%pi-%pi/2)];
    zprt=[tr*sin(q-dth+%pi-%pi/2)+ztc , tr*sin(q+%pi-%pi/2)+ztc];
    plot(xprt,zprt,'--');
end

//地面
xplt=[-2,2];
zplt=[0,0];
plot(xplt,zplt,);

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
//初期角度
th1=18.125/180*%pi;
th2=33.664/180*%pi;
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

tu2=0;//トルクの初期値

//コントローラー関数---------------------------------------------------------------------------------------------------------------------------------------
T=0.001;//予測する時間の長さ
dti=0.0001;//予測するステップの細かさ
//最適なトルク]を考える
function j=fj(x,u)
    th1i=x(1);
    th12i=x(2);
    th123i=x(3);
    th2i=th12i-th1i;
    th3i=th123i-th12i;
    thd1i=x(4);
    thd12i=x(5);
    thd123i=x(6);
    tu2i=u;

    THi=[th1i;th12i;th123i];
    THdi=[thd1i;thd12i;thd123i];
    THddi=THdd;
    
    je=0;//積分を行うための変数
    for t=dti:dti:T
        M3=[L1^2 -L1*L2 -L1*l3g;-L1*L2 L2^2 L2*l3g;-L1*l3g L2*l3g l3g^2] + [1 0 0;0 0 0;0 0 0]*2*l1r*(L1-l1r)*(cos(th1i)-1) + [0 1 0;1 0 0;0 0 0]*L2*(-L1*(cos(th2i)-1)-l1r*(cos(th12i)-cos(th2i))) + [0 0 1;0 0 0;1 0 0]*l3g*(-L1*(cos(th2i+th3i)-1) - l1r*(cos(th123i)-cos(th2i+th3i))) + [0 0 0;0 0 1;0 1 0]*L2*l3g*(cos(th3i)-1);

        M2=[L1^2 -L1*l2g 0;-L1*l2g l2g^2 0;0 0 0] + [1 0 0;0 0 0;0 0 0]*2*l1r*(L1-l1r)*(cos(th1i)-1) + [0 1 0;1 0 0;0 0 0]*l2g*(-L1*(cos(th2i)-1)-l1r*(cos(th12i)-cos(th2i)));

        M1=[l1g^2 0 0;0 0 0;0 0 0] + [1 0 0;0 0 0;0 0 0]*2*l1r*(l1g-l1r)*(cos(th1i)-1);

        NL= [-(m1*l1g+m2*L1+m3*L1-M*l1r)*l1r*sin(th1i);-(m2*l2g + m3*L2)*(L1-l1r)*sin(th2i);-m3*l3g*(L1-l1r)*sin(th2i+th3i)]*thd1i^2 + [(m2*l2g+m3*L2)*((L1-l1r)*sin(th2i)+l1r*sin(th12i));0;m3*L2*l3g*sin(th3i)]*thd12i^2 + [m3*l3g*((L1-l1r)*sin(th2i+th3i)+l1r*sin(th123i));-m3*L2*l3g*sin(th3i);0]*thd123i^2;

        THddi=inv(m1*M1+m2*M2+m3*M3+[I1 0 0;0 I2 0;0 0 I3])*([-tu2i;tu2i;0]-[-(m1*l1g+m2*L1+m3*L1-M*l1r) 0 0;0 m2*l2g+m3*L2 0;0 0 m3*l3g]*[sin(th1i);sin(th12i);sin(th123i)]*g-NL);

        //コントローラーの積分はode関数なしで
        //前の角加速度から現在の角速度，角度を求める
        THd_a=THdi;
        THdi=THdi+THddi*dti;
        THi=THi+(THd_a+THdi)/2*dti;

        th1i=THi(1,1);
        th12i=THi(2,1);
        th123i=THi(3,1);
        th2i=th12i-th1i;
        th3i=th123i-th12i;
        
        //xr=(th1i-x0(1,1))*l1r;//足裏接点は必要な際th1iから求める。
        //x3=xr+(L1-l1r)*sin(thi1)-L2*sin(th12i)-L3*sin(th123i);
        //z3=l1r+(L1-l1r)*cos(th1i)-L2*cos(th12i)-L3*cos(th123i);
        //e=( ( (z3-ztc)^2 + x3^2 )^0.5 -tr)
        //je=je+(e*dti)^2;//目標との差の2乗の積分
    end
    
        xr=(th1i-x0(1,1))*l1r;//足裏接点は必要な際th1iから求める。
        x3=xr+(L1-l1r)*sin(th1i)-L2*sin(th12i)-L3*sin(th123i);
        z3=l1r+(L1-l1r)*cos(th1i)-L2*cos(th12i)-L3*cos(th123i);
    //j=ke*je ;//+ ku*(u*T)^2 + kud*((u-tu2)/dt)^2;//評価関数
    j=(( ( (z3-ztc)^2 + x3^2 )^0.5 -tr))^2

endfunction



maxu=100;
minu=-100;
h=0.000000000001;
eta=1;//勾配法の学習率
function ans=gradent(left,right,x)
    f_left=fj(x,left);
    f_right=fj(x,right);
    minyU=left;
    if f_left>f_right then
        minyU=right;
    end
    

    a=(left+right)/2;
    f_t=0;
    t_h=0;
    while 1 
        //s=(fj(x,a+h/2)-fj(x,a-h/2))/h;
        f_t=fj(x,a);
        f_h=fj(x,a+h);
        //s=(fj(x,a1)-fj(x,a2))/h;
        sq=(f_h-f_t)/h;
        if sq<0.001 then        //傾きが0.0001より小さくなったら勾配法終了
            break;
        end
        a=a-eta*sq;//傾きと学習率に従って、測定点をずらす
        //disp(f_t);
    end
    
    if (left<a) && (a<right) && (minyU>f_t) then   //範囲外の解は入れない.かつ、いままで見つけた極小値で最も小さい
        minyU = a;
    end
    
    
    
    
        a=(right-left)/3+left;
    f_t=0;
    t_h=0;
    while 1 
        //s=(fj(x,a+h/2)-fj(x,a-h/2))/h;
        f_t=fj(x,a);
        f_h=fj(x,a+h);
        //s=(fj(x,a1)-fj(x,a2))/h;
        sq=(f_h-f_t)/h;
        if sq<0.0001 then        //傾きが0.0001より小さくなったら勾配法終了
            break;
        end
        a=a-eta*sq;//傾きと学習率に従って、測定点をずらす
        //disp(f_t);
    end
    
    if (left<a) && (a<right) && (minyU>f_t) then   //範囲外の解は入れない.かつ、いままで見つけた極小値で最も小さい
        minyU = a;
    end
    
    
    
            a=(right-left)*2/3+left;
    f_t=0;
    t_h=0;
    while 1 
        //s=(fj(x,a+h/2)-fj(x,a-h/2))/h;
        f_t=fj(x,a);
        f_h=fj(x,a+h);
        //s=(fj(x,a1)-fj(x,a2))/h;
        sq=(f_h-f_t)/h;
        if sq<0.0001 then        //傾きが0.0001より小さくなったら勾配法終了
            break;
        end
        a=a-eta*sq;//傾きと学習率に従って、測定点をずらす
        //disp(f_t);
    end
    
    if (left<a) && (a<right) && (minyU>f_t) then   //範囲外の解は入れない.かつ、いままで見つけた極小値で最も小さい
        minyU = a;
    end
    
    
    ans=minyU;
endfunction


//微分方程式
function xdot=flink(t,x)
    th1=x(1);
    th12=x(2);
    th123=x(3);
    thd1=x(4);
    thd12=x(5);
    thd123=x(6);
    th2=th12-th1;
    th3=th123-th12;
    
    //コントローラー
    tu2=gradent(minu,maxu,x);
    //tu2=ke*th12;

    M3=[L1^2 -L1*L2 -L1*l3g;-L1*L2 L2^2 L2*l3g;-L1*l3g L2*l3g l3g^2] + [1 0 0;0 0 0;0 0 0]*2*l1r*(L1-l1r)*(cos(th1)-1) + [0 1 0;1 0 0;0 0 0]*L2*(-L1*(cos(th2)-1)-l1r*(cos(th12)-cos(th2))) + [0 0 1;0 0 0;1 0 0]*l3g*(-L1*(cos(th2+th3)-1) - l1r*(cos(th123)-cos(th2+th3))) + [0 0 0;0 0 1;0 1 0]*L2*l3g*(cos(th3)-1);

    M2=[L1^2 -L1*l2g 0;-L1*l2g l2g^2 0;0 0 0] + [1 0 0;0 0 0;0 0 0]*2*l1r*(L1-l1r)*(cos(th1)-1) + [0 1 0;1 0 0;0 0 0]*l2g*(-L1*(cos(th2)-1)-l1r*(cos(th12)-cos(th2)));

    M1=[l1g^2 0 0;0 0 0;0 0 0] + [1 0 0;0 0 0;0 0 0]*2*l1r*(l1g-l1r)*(cos(th1)-1);

    NL= [-(m1*l1g+m2*L1+m3*L1-M*l1r)*l1r*sin(th1);-(m2*l2g + m3*L2)*(L1-l1r)*sin(th2);-m3*l3g*(L1-l1r)*sin(th2+th3)]*thd1^2 + [(m2*l2g+m3*L2)*((L1-l1r)*sin(th2)+l1r*sin(th12));0;m3*L2*l3g*sin(th3)]*thd12^2 + [m3*l3g*((L1-l1r)*sin(th2+th3)+l1r*sin(th123));-m3*L2*l3g*sin(th3);0]*thd123^2;

    THdd=inv(m1*M1+m2*M2+m3*M3+[I1 0 0;0 I2 0;0 0 I3])*([-tu2;tu2;0]-[-(m1*l1g+m2*L1+m3*L1-M*l1r) 0 0;0 m2*l2g+m3*L2 0;0 0 m3*l3g]*[sin(th1);sin(th12);sin(th123)]*g-NL);

    xdot(1)=thd1;
    xdot(2)=thd12;
    xdot(3)=thd123;
    xdot(4)=THdd(1);
    xdot(5)=THdd(2);
    xdot(6)=THdd(3);
    
endfunction

x0=[th1;th12;th123;0;0;0];
t=0:dt:0.3;              // 評価時間　1[ms]ごと
t0=0;                       // 時刻の初期値
x=ode(x0,t0,t,flink);       // 微分方程式のソルバ

idx=1:20:301;				// プロット 50[ms]=0.05[s]ごと
xq1=x(1,idx)';//行、列を指定して入れる
xq2=x(2,idx)';
xq3=x(3,idx)';

//リンク描写(17まで)
count=0;
for i=1:15
    
    th1=xq1(i);
    th12=xq2(i);
    th123=xq3(i);
    count=count+1;//色設定

    //リンク描写
            xr=(th1-x0(1,1))*l1r;
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
 
            xm=(m1*xg1+m2*xg2+m3*xg3)/M;
            zm=(m1*zg1+m2*zg2+m3*zg3)/M;
            xplt=[xm,xm];
            zplt=[0.01+zm,-0.01+zm];
            plot2d(xplt,zplt,count);
            xplt=[0.01+xm,-0.01+xm];
            zplt=[zm,zm];
            plot2d(xplt,zplt,count);
            
end
