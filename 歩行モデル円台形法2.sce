clear;//足裏円モデル
figure(1);
//subplot(311);
m3=0.1;
m2=0.4;
m1=m3*2;
M=m1+m2+m3;

Sth=75/180*%pi;//股傾斜角度
L3=0.8;
L1=L3*sin(Sth);
L2=0.4;
Lg=(L3*cos(Sth)+L2/tan(%pi/3))*tan(%pi/3);//初期位置での原点から足先端までの距離
l1g=0.4;
l2g=0.2;
l3g=0.4;
l1r=0.1;//足裏円形状の半径

//初期角度
th1=13.181/180*%pi;
th2=37.991/180*%pi;
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
x0=[th1;th12;th123;0;0;0];//初期状態

I1=m1*L1*L1/12;
I2=m2*L2*L2/12;
I3=m3*L3*L3/12;

g=9.80665;

xr=0;//足裏接点初期値
tu2=-2;//トルクの初期値（最初の勾配法の初期値になる）
tu2d=-1;

endT=0.3;//シミュレーション終了時刻
dt=0.001;//シミュレーションの周期

lookT=10;//描写するリンクの間隔


T=0.015;//予測する時間の長さ（現在から予測ホライゾンまでの時間）
dti=0.001;//予測時のステップの周期

maxu=200;//解を探す最大値(無効化)
minu=-200;//解を探す最小値(無効化)
eta=900000;//勾配法の学習率

E_sp=0.000000015; E_sp_q=E_sp^2;//勾配法の際、探索を終了する傾き、その2乗

h=0.0000000001;//傾きを計算する際の極小の値

//ゲイン一覧
ke=1;
ku=0;
kud=0;

ztc=-3.557;//目標軌跡円の中心点（x=0とする）

x3=xr+(L1-l1r)*sin(th1)-L2*sin(th12)-L3*sin(th123);
z3=l1r+(L1-l1r)*cos(th1)-L2*cos(th12)-L3*cos(th123);
tr=( (z3-ztc)^2 + x3^2 )^0.5;//目標軌跡円の半径（必ず目標軌跡状に足先端がくるようになる）

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
plot(xplt,zplt);



//コントローラー関数---------------------------------------------------------------------------------------------------------------------------------------
//予測の内部もode関数を使ったようがいいかもしれない。
//高速化をしないとやっていけないかもしれない。


//最適なトルクの時間微分を考える
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
    
    je=0;//積分を行うための変数
    for t=dti:dti:T
        M3=[L1^2 -L1*L2 -L1*l3g;-L1*L2 L2^2 L2*l3g;-L1*l3g L2*l3g l3g^2] + [1 0 0;0 0 0;0 0 0]*2*l1r*(L1-l1r)*(cos(th1i)-1) + [0 1 0;1 0 0;0 0 0]*L2*(-L1*(cos(th2i)-1)-l1r*(cos(th12i)-cos(th2i))) + [0 0 1;0 0 0;1 0 0]*l3g*(-L1*(cos(th2i+th3i)-1) - l1r*(cos(th123i)-cos(th2i+th3i))) + [0 0 0;0 0 1;0 1 0]*L2*l3g*(cos(th3i)-1);

        M2=[L1^2 -L1*l2g 0;-L1*l2g l2g^2 0;0 0 0] + [1 0 0;0 0 0;0 0 0]*2*l1r*(L1-l1r)*(cos(th1i)-1) + [0 1 0;1 0 0;0 0 0]*l2g*(-L1*(cos(th2i)-1)-l1r*(cos(th12i)-cos(th2i)));

        M1=[l1g^2 0 0;0 0 0;0 0 0] + [1 0 0;0 0 0;0 0 0]*2*l1r*(l1g-l1r)*(cos(th1i)-1);

        NL= [-(m1*l1g+m2*L1+m3*L1-M*l1r)*l1r*sin(th1i);-(m2*l2g + m3*L2)*(L1-l1r)*sin(th2i);-m3*l3g*(L1-l1r)*sin(th2i+th3i)]*thd1i^2 + [(m2*l2g+m3*L2)*((L1-l1r)*sin(th2i)+l1r*sin(th12i));0;m3*L2*l3g*sin(th3i)]*thd12i^2 + [m3*l3g*((L1-l1r)*sin(th2i+th3i)+l1r*sin(th123i));-m3*L2*l3g*sin(th3i);0]*thd123i^2;

        THddi=inv(m1*M1+m2*M2+m3*M3+[I1 0 0;0 I2 0;0 0 I3])*([-tu2i;tu2i;0]-[-(m1*l1g+m2*L1+m3*L1-M*l1r) 0 0;0 m2*l2g+m3*L2 0;0 0 m3*l3g]*[sin(th1i);sin(th12i);sin(th123i)]*g-NL);

        //前の角加速度から現在の角速度，角度を求める
        THd_a=THdi;
        THdi=THdi+THddi*dti;
        THi=THi+(THd_a+THdi)/2*dti;

        th1i=THi(1,1);
        th12i=THi(2,1);
        th123i=THi(3,1);
        th2i=th12i-th1i;
        th3i=th123i-th12i;
    end
    xr=(th1i-x0(1,1))*l1r;//足裏接点は必要な際th1iから求める。
    x3=xr+(L1-l1r)*sin(th1i)-L2*sin(th12i)-L3*sin(th123i);
    z3=l1r+(L1-l1r)*cos(th1i)-L2*cos(th12i)-L3*cos(th123i);
        
    j=ke*( ( (z3-ztc)^2 + x3^2 )^0.5 -tr)^2
endfunction








function ans=gradent(left,right,x)
    /*
    f_left=fj(x,left);
    f_right=fj(x,right);
    minyU=left;
    minyY=f_left
    if f_left>f_right then
        minyU=right;
        minyY=f_right;
    end
   */
    
    
    //(right+left)/2;
    a=tu2;//前回の解を初期値にして解を探す
    etai=eta;
    sq=0;
    while 1
        f_t=fj(x,a);
        sq_af=sq
        sq=(fj(x,a+h) -f_t)/h;
        printf("%lf ",sq);
        if (sq*sq)<E_sp_q then        //傾きが0.0000000001より小さくなったら勾配法終了
            break;
        end

        if sq*sq_af < 0 then//極値をこえる際に学習率を下げる
            etai=etai/2;
        end
        
        a=a-etai*sq;//傾きと学習率に従って、測定点をずらす
    end
    //disp(f_t);
    //if (left<a) && (a<right) && (minyY>f_t) then   //範囲外の解は入れない.かつ、いままで見つけた極小値で最も小さい      //上限をなくしました
        
        minyU = a;
        //minyY=f_t;
    //end
    
    /*
    //初めの答えは-1.74063235くらいになるはず
   XXU=-2.0969654702:0.00000000001:-2.0969654700;
   XXF=[];
    figure(1);
    for i=XXU
        XXF=[XXF fj(x,i)];
    end
    plot(XXU,XXF);
*/
    
    
    ans=minyU;
endfunction



TH=[th1;th12;th123];
THd=[0;0;0];

Xth=[TH];
Xtu2=[];

for t=dt:dt:endT

    tu2=gradent(minu,maxu,[TH;THd]);
    Xtu2=[Xtu2 tu2];
    printf("\nt=%lf tu2=%lf\n",t,tu2);
    //現在の状態からthdd1,thdd12,thdd123を求める
    M3=[L1^2 -L1*L2 -L1*l3g;-L1*L2 L2^2 L2*l3g;-L1*l3g L2*l3g l3g^2] + [1 0 0;0 0 0;0 0 0]*2*l1r*(L1-l1r)*(cos(th1)-1) + [0 1 0;1 0 0;0 0 0]*L2*(-L1*(cos(th2)-1)-l1r*(cos(th12)-cos(th2))) + [0 0 1;0 0 0;1 0 0]*l3g*(-L1*(cos(th2+th3)-1) - l1r*(cos(th123)-cos(th2+th3))) + [0 0 0;0 0 1;0 1 0]*L2*l3g*(cos(th3)-1);

    M2=[L1^2 -L1*l2g 0;-L1*l2g l2g^2 0;0 0 0] + [1 0 0;0 0 0;0 0 0]*2*l1r*(L1-l1r)*(cos(th1)-1) + [0 1 0;1 0 0;0 0 0]*l2g*(-L1*(cos(th2)-1)-l1r*(cos(th12)-cos(th2)));

    M1=[l1g^2 0 0;0 0 0;0 0 0] + [1 0 0;0 0 0;0 0 0]*2*l1r*(l1g-l1r)*(cos(th1)-1);

    NL= [-(m1*l1g+m2*L1+m3*L1-M*l1r)*l1r*sin(th1);-(m2*l2g + m3*L2)*(L1-l1r)*sin(th2);-m3*l3g*(L1-l1r)*sin(th2+th3)]*thd1^2 + [(m2*l2g+m3*L2)*((L1-l1r)*sin(th2)+l1r*sin(th12));0;m3*L2*l3g*sin(th3)]*thd12^2 + [m3*l3g*((L1-l1r)*sin(th2+th3)+l1r*sin(th123));-m3*L2*l3g*sin(th3);0]*thd123^2;

    THdd=inv(m1*M1+m2*M2+m3*M3+[I1 0 0;0 I2 0;0 0 I3])*([-tu2;tu2;0]-[-(m1*l1g+m2*L1+m3*L1-M*l1r) 0 0;0 m2*l2g+m3*L2 0;0 0 m3*l3g]*[sin(th1);sin(th12);sin(th123)]*g-NL);
    

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
     Xth=[Xth TH];

end


//リンク描写
idx=1:lookT:int(endT/dt)+1;				// プロット 20[ms]=0.02[s]ごと
xq1=Xth(1,idx)';//行、列を指定して入れる
xq2=Xth(2,idx)';
xq3=Xth(3,idx)';

count=0;
for i=1:int(endT/dt)/lookT+1
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

disp(Xtu2);
