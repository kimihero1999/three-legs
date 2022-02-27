figure(20);
clear ;

m1=0.4 ; m2=1.0 ; m3=0.1 ; m4=0.1 ;
L1=0.781 ; L2=0.2 ; L3=0.4 ; L4=0.4 ;
I1=m1*L1*L1/12 ; I2=m2*L2*L2/12 ; I3=m3*L3*L3/12 ; I4=m4*L4*L4/12
L1G=L1/2 ; L2G=L2/2 ; L3G=L3/2 ; L4G=L4/2
L1R=L1*0.05 ;
g=9.81 ;

q10=-5*%pi/180;
th20=5*%pi/180;
q20=q10+th20;
q30=45*%pi/180;
th40=0*%pi/180;
q40=q30+th40;

q3ref=-10*%pi/180;

x10=(L1G-L1R)*sin(q10);
x20=(L1 -L1R)*sin(q10)-L2G*sin(q20);
x30=(L1 -L1R)*sin(q10)-L2 *sin(q20)-L3G*sin(q30);
x40=(L1 -L1R)*sin(q10)-L2 *sin(q20)-L3 *sin(q30)-L4G*sin(q40);

x0=(m1*x10+m2*x20+m3*x30+m4*x40)/(m1+m2+m3+m4);
disp("x0=",x0);

A_th2  =120*%pi/180;
t_th2  =2.4 ;
om_th2 =2*%pi/t_th2;

A_th4  =120*%pi/180;
t_th4  =t_th2/2 ;
om_th4 =2*%pi/t_th4;

kp_q3=3;
kd_q3=0.1;

E_b=[1 0 0 0;1 0 1 0;0 1 0 0;0 1 0 1];



tu3=0;//トルクの初期値


T=0.02;//予測する時間の長さ（現在から予測ホライゾンまでの時間）
dti=0.001;//予測時のステップの周期

maxu=200;//解を探す最大値(無効化)
minu=-200;//解を探す最小値(無効化)
eta=150000;//勾配法の学習率

h=0.0000000001;//傾きを計算する際の極小の値

//ゲイン一覧
ke=1;



xr=0;//足裏接点初期値
ztc=-3.557;//目標軌跡円の中心点（x=0とする）
xx=xr+(L1-L1R)*sin(q10)-L2*sin(q20)-L3*sin(q30)-L4*sin(q40);
zz=L1R+(L1-L1R)*cos(q10)-L2*cos(q20)-L3*cos(q30)-L4*cos(q40);
tr=( (zz-ztc)^2 + xx^2 )^0.5;//目標軌跡円の半径（必ず目標軌跡状に足先端がくるようになる）


//トルクの評価
function j=fj(x,ft,u)//fuは現在の時間
    
    q1=x(1);
    q3=x(2);
    q1d=x(3);
    q3d=x(4);
    
    // Hip flexure joint に制御トルクを与える
    tau3 = u;
    
    Qi=[q1;q3];
    Qdi=[q1d;q2d];
    
    je=0;//積分を行うための変数
    
    for t=ft:dti:T+ft
        
        if t< t_th4/2 then
            th2   = -A_th2*sin(om_th2*t) + th20 ;
            th2d  = -om_th2*A_th2*cos(om_th2*t);
            th2dd =  om_th2^2*A_th2*sin(om_th2*t);
            
            th4   = A_th4*sin(om_th4*t) + th40 ;
            th4d  = om_th4*A_th4*cos(om_th4*t);
            th4dd = -om_th4^2*A_th4*sin(om_th4*t);
        else
            th2=0;th2d=0;th2dd=0;
            th4=0;th4d=0;th4dd=0;
        end

        q2=q1+th2;
        q2d=q1d+th2d;
        q4=q3+th4;
        q4d=q3d+th4d;    
    
        M411=(L1-L1R)^2+L1R^2+2*L1R*(L1-L1R)*cos(q1);
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
    
        M_b=E_b'*M*E_b ;
        h_b=E_b'*h ;
        G_b=E_b'*G ;



        tau_b=[ -tau3 ; tau3 ];
        tmp=inv( M_b(1:2,1:2) )*( tau_b - M_b(1:2,3:4)*[th2dd;th4dd] - h_b(1:2) - G_b(1:2) );

        //前の角加速度から現在の角速度，角度を求める
        Qd_a=Qdi;
        Qdi=Qdi+tmp*dti;
        Qi=Qi+(Qd_a+Qdi)/2*dti;
        
        
        q1=x(1,:)';
        th2=-max(0,A_th2*sin(om_th2*t)) + th20;
        q2=q1+th2;
        q3=x(2,:)';
        th4= max(0,A_th4*sin(om_th4*t)) + th40;
        q4=q3+th4;
        
        
        xq1=Qi(1,1);
        
        th2=-max(0,A_th2*sin(om_th2*t)) + th20;
        xq2=q1+th2;
        
        xq3=Qi(2,1);
        
        th4= max(0,A_th4*sin(om_th4*t)) + th40;
        xq4=q3+th4;
        
        xr=(xq1-q10)*L1R;//足裏接点は必要な際th1iから求める。
        x4a=xr+(L1-L1R)*sin(xq1)-L2*sin(xq2)-L3*sin(xq3)-L4*sin(xq4);
        z4a=L1R+(L1-L1R)*cos(xq1)-L2*cos(xq2)-L3*cos(xq3)-L4*cos(xq4);
        e_af=e;//初期のeを定義しておく必要がある。
        e=( ( (z4a-ztc)^2 + x4a^2 )^0.5 -tr)^2
        je=je+(+e_af)/2;
    end
    
    j=ke*je
endfunction






function ans=gradent(left,right,t,x)
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
    a=tu3;//前回の解を初期値にして解を探す
    etai=eta;
    sq=0;
    while 1
        f_t=fj(x,t,a);
        sq_af=sq
        sq=(fj(x,t,a+h) -f_t)/h;
        printf("%lf ",sq);
        if (sq*sq)^0.5<0.0000000001 then        //傾きが0.0000000001より小さくなったら勾配法終了
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
























function xdot=two_link(t,x)   // 微分方程式を解くための関数
    q1=x(1);
    q3=x(2);
    q1d=x(3);
    q3d=x(4);

    if t< t_th4/2 then
        th2   = -A_th2*sin(om_th2*t) + th20 ;
        th2d  = -om_th2*A_th2*cos(om_th2*t);
        th2dd =  om_th2^2*A_th2*sin(om_th2*t);

        th4   = A_th4*sin(om_th4*t) + th40 ;
        th4d  = om_th4*A_th4*cos(om_th4*t);
        th4dd = -om_th4^2*A_th4*sin(om_th4*t);
    else
        th2=0;th2d=0;th2dd=0;
        th4=0;th4d=0;th4dd=0;
    end
    
    q2=q1+th2;
    q2d=q1d+th2d;
    q4=q3+th4;
    q4d=q3d+th4d;    
    
    M411=(L1-L1R)^2+L1R^2+2*L1R*(L1-L1R)*cos(q1);
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
    
    M_b=E_b'*M*E_b ;
    h_b=E_b'*h ;
    G_b=E_b'*G ;
    
    xr=(q1-q10)*L1R;//足裏接点は必要な際th1iから求める。
    x4a=xr+(L1-L1R)*sin(q1)-L2*sin(q2)-L3*sin(q3)-L4*sin(q4);
    z4a=L1R+(L1-L1R)*cos(q1)-L2*cos(q2)-L3*cos(q3)-L4*cos(q4);
    e=( ( (z4a-ztc)^2 + x4a^2 )^0.5 -tr)^2
    
// Hip flexure joint に制御トルクを与える
    tau3 = gradent(-100,100,t,x) ;

    tau_b=[ -tau3 ; tau3 ];
    tmp=inv( M_b(1:2,1:2) )*( tau_b - M_b(1:2,3:4)*[th2dd;th4dd] - h_b(1:2) - G_b(1:2) );
    
    xdot(1)=q1d;
    xdot(2)=q3d;
    xdot(3)=tmp(1);
    xdot(4)=tmp(2);
endfunction




x0=[q10;q30;0;0];
t=0:0.001:1.0;              // 評価時間　1[ms]ごと
t0=0;                       // 時刻の初期値
x=ode(x0,t0,t,two_link);       // 微分方程式のソルバ
q1=x(1,:)';
th2=-max(0,A_th2*sin(om_th2*t')) + th20;
q2=q1+th2;
q3=x(2,:)';
th4= max(0,A_th4*sin(om_th4*t')) + th40;
q4=q3+th4;

idx=1:25:1001;				// プロット 25[ms]=0.025[s]ごと
xq1=q1(idx);
xq2=q2(idx);
xq3=q3(idx);
xq4=q4(idx);
tt=t(1,idx)';

x0a=L1R*xq1;
x00a=x0a-L1R*sin(xq1);
x1a=L1R*xq1+(L1-L1R)*sin(xq1);///////////////////////////////////////////////////
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
ca.data_bounds=[-1 -0.2 ; 22 0.9 ] ;

subplot(212)
plot2d(t',180/%pi*q1,1);
plot2d(t',180/%pi*q2,2);
plot2d(t',180/%pi*q3,3);
plot2d(t',180/%pi*q4,4);
plot2d(t',180/%pi*(q3-q2),5);
legend(['q1';'q2';'q3';'q4';'th3']);
xgrid

