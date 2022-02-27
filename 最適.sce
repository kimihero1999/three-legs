clear ;//最適
a=100;
b=2900;


A=[-a 0;1 0];
B=[b;0];
poles=[-100+100*%i,-100-100*%i]//Kに使われているからいらないやつ
K=ppol(A,B,poles);//Kは使ってないぞー



Q=[1 0;0 1];
R=1;
P=ricc(A,B*inv(R)*B',Q,'cont');
Kfb2=inv(R)*B'*P;
disp("Kfb2=",Kfb2);

AA=A'*P+P*A+Q-P*B*inv(R)*B'*P;//AA＝0をみたすPを探すのがリカッチ方程式
//0になっているか確認するのもあり

k1=Kfb2(1,1);
k2=Kfb2(1,2);

function xdot=motor(t,x)   // 微分方程式を解くための関数
    omega=x(1);                 // 変位
    theta=x(2);              // 速度
    u=-k1*omega-k2*theta;
    omega_dot=-a*omega+b*u;
    
    xdot(1)=omega_dot;           // 関数の返り値
    xdot(2)=omega;
endfunction

x0=[0;1];                   // 状態変数の初期値
t=0:0.001:0.1;                 // 評価時間
t0=0;                       // 時刻の初期値
x=ode(x0,t0,t,motor);      // 微分方程式のソルバ

subplot(211)
plot(t,x(2,:));             // 第一変数の時刻歴をプロット
subplot(212)
u=-k1*x(1,:)-k2*x(2,:);
plot(t,u);


disp("AA",AA);
