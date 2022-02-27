clear;//特論課題
//プロット画面のクリア
figure(1);clf(1,"reset");
M = 0.79 ;              // 台車の質量[kg]
m = 0.04 ;              // 振子の質量[kg]
b  = 0.027 ;            // 車輪の半径[m]
L = 0.2 ;               // 振子の長さ = 2L[m]
G  = 11.6 ;             // ギア比
R  = 27 ;               // モータ内部抵抗[Ω]
Kt = 2.1703e-2 ;        // トルク定数[N*m]
Ke = 2.1703e-2 ;        // 逆起電力定数[V]
J1 = 5.7622e-6 ;        // 車輪の重心周りの慣性M[kg*m^2] 
I = 1/3*(m*L^2) ;      // 振子の重心周りの慣性M[kg*m^2]
g  = 9.81 ;             // 重力加速度[kg/s^2]

MM=[M+m+J1*(G/b)^2 m*L ; m*L I+m*L^2] ;
KK=[0 0;0 m*g*L] ;
CC=[-Kt*Ke/R*(G/b)^2 0;0 0] ;
NN=[G*Kt/b/R;0];

A=[zeros(2,2) eye(2,2) ; inv(MM)*KK inv(MM)*CC ] ;
B=[zeros(2,1) ; inv(MM)*NN ] ;
C=[1 0 0 0;0 1 0 0];

qq=0.1
rr=10
Q=[qq 0 0 0;0 qq 0 0;0 0 qq 0;0 0 0 qq];
R=rr;
P=ricc(A,B*inv(R)*B',Q,'cont');
Kfb=inv(R)*B'*P;
disp("Kfb2=",Kfb);

AA=A'*P+P*A+Q-P*B*inv(R)*B'*P;//AA＝0をみたすPを探すのがリカッチ方程式

disp("Kfb=",Kfb);
//Kfb=ppol(A,B,[-20 -20 -20 -20]);
//オブザーバーゲインを決定
HH=0.5
Lhat =  ppol( A' ,C',[-5000*HH,-160*HH,-4*HH,-7*HH]);

Lobs=Lhat';
disp("Lobs=",Lobs);

//実際の制御対象
M = 0.79*1  ;              // 台車の質量[kg]
m = 0.04 ;              // 振子の質量[kg]
L = 0.2*1 ;               // 振子の長さ = 2L[m]
I = 1/3*(m*L^2) ;      // 振子の重心周りの慣性M[kg*m^2]
//20%くらいの誤差に耐えるゲイン(台車の質量と振り子の長さに誤差が出る想定)

function xdot=pendulum(t,x)
    xhat = x(5:8);
    //状態推測値を用いたフィードバック
    v = -Kfb*xhat
    
    //操作量の飽和
    //if v>12.0 then v=12.0 ;
    //elseif v<-12.0 then v=-12.0;
    //end
    
        r = x(1); phi = x(2); rd = x(3); phid = x(4);
        Mtmp =[M+m+J1*(G/b)^2 m*L*cos(phi) ; m*L*cos(phi) I+m*L^2] ;
        tmp1=[m*L*phid^2*sin(phi) - Kt*Ke/R*(G/b)^2*rd + G*Kt/b/R*v ; m*g*L*sin(phi)] ;
        tmp2=inv(Mtmp)*tmp1 ;
        xdot(1)=x(3) ; xdot(2)=x(4) ; xdot(3)=tmp2(1) ; xdot(4)=tmp2(2);
        
        //オブザーバーの計算式
        u=v;
        y=[r;phi];//測定信号
        xhat=x(5:8);
        xhat_dot = A*xhat + B*v + Lobs*(y-C*xhat);
        xdot(5:8)=xhat_dot;
endfunction

x0=[0.05;5*%pi/180;0;0];   // 初期条件
xhat0=[0;0;0;0];
t0=0;                   // 初期時間
t=0:0.00001:10;             // 解を計算する時間

// 常微分方程式ソルバ
// ode は,dy/dt=f(t,y) , y(t0)=y0で定義された明示的なODEシステムを解くための標準関数です. 
x=ode([x0;xhat0],t0,t,pendulum);

figure(1);clf(1);
subplot(321) ;plot2d(t',x([1,5],:)'); xgrid(1);
ylabel("$r,\hat{r}$");
ca = gca();ca.font_size = 2 ;ca.y_label.font_size = 4 ;
subplot(322) ;plot2d(t',x([2,6],:)'); xgrid(1);
ylabel("$\phi,\hat{\phi}$");
ca = gca();ca.font_size = 2 ;ca.y_label.font_size = 4 ;
subplot(323) ;plot2d(t',x([3,7],:)'); xgrid(1);
ylabel("$\dot{r},\hat{\dot{r}$");
ca = gca();ca.font_size = 2 ;ca.y_label.font_size = 4 ;
subplot(324) ;plot2d(t',x([4,8],:)'); xgrid(1);
ylabel("$\dot{\phi},\hat{\dot{\phi}}$");
ca = gca();ca.font_size = 2 ;ca.y_label.font_size = 4 ;
//操作量の飽和
vv=-Kfb*x(5:8,:);
//vv=min(vv,12);
//vv=max(vv,-12);


subplot(325) ;plot2d(t,vv); xgrid(1);
ylabel("$v$");
ca = gca();ca.font_size = 2 ;ca.y_label.font_size = 4 ;
//タイトル
subplot(326);
xtitle("台車の質量+20%，振り子の長さ-20%の時");
ca=gca();ca.title.font_size=3;ca.title.position=[0.1,0.5];
