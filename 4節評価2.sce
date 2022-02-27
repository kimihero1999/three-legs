clear ;
figure(2);
subplot(221);
m4=0.1; m3=0.1; m2=2; m1=(m3+m4)*2; msum=m1+m2+m3+m4;//各関節質量と全体重量
g=9.80665 ;

Sth=75/180*%pi;//股傾斜角度
L4=0.4; L3=0.4; L1=(L3+L4)*sin(Sth); L2=0.2;//各関節長（スタンスレッグの股傾斜）
Lg=(L3*cos(Sth)+L2/tan(%pi/3))*tan(%pi/3);//初期位置での原点から足先端までの距離
L1G=L1/2 ; L2G=L2/2 ; L3G=L3/2 ; L4G=L4/2
L1R=L1*0.05 ;//スタンスレッグの足裏半径
I1=m1*L1*L1/12 ; I2=m2*L2*L2/12 ; I3=m3*L3*L3/12 ;I4=m4*L4*L4/12;//各関節の慣性モーメント
//初期角度
th10=6.037/180*%pi; th20=33.747/180*%pi; th30=0/180*%pi; th40=0/180*%pi;
//初期絶対角度角度
q10=th10; q20=q10+th20; q30=q20+th30; q40=q30+th40;
//初期角速度は0のため定義していない

x10=(L1G-L1R)*sin(q10);
x20=(L1 -L1R)*sin(q10)-L2G*sin(q20);
x30=(L1 -L1R)*sin(q10)-L2 *sin(q20)-L3G*sin(q30);
x40=(L1 -L1R)*sin(q10)-L2 *sin(q20)-L3 *sin(q30)-L4G*sin(q40);

xr=0;//足裏接点初期値
tu2=0; tu3=0; tu4=0;//トルクの初期値（最初の勾配法の初期値になる）

endT=0.36;//シミュレーション終了時刻
dt=0.001;//シミュレーションの周期
lookT=20;//描写するリンクの間隔  0.03秒間隔
T=0.01;//予測する時間の長さ（現在から予測ホライゾンまでの時間）
dti=0.001;//予測時のステップの周期

maxu=200;//解を探す最大値(無効化)
minu=-200;//解を探す最小値(無効化)
eta20=200; eta30=200; eta40=200;//勾配法の学習率(初期値)

E_sp=0.000000015; E_sp_q=E_sp^2;//勾配法の際、探索を終了する傾き、その2乗

h=0.0000000001;//傾きを計算する際の極小の値

//ゲイン一覧
ke=4; ktu=0.01; kud=0;

ztc=-3.557;//目標軌跡円の中心点（x=0とする）

x40=(L1 -L1R)*sin(q10)-L2 *sin(q20)-L3 *sin(q30)-L4*sin(q40);
z40=L1R+(L1-L1R)*cos(q10)-L2*cos(q20)-L3*cos(q30)-L4*cos(q40);
tr=( (z40-ztc)^2 + x40^2 )^0.5;//目標軌跡円の半径（必ず目標軌跡状に足先端がくるようになる）

//目標軌跡円を描く
hth=10/180*%pi;//円の幅
dth=1/180*%pi;//線の細かさ（角度）
for q=-(hth-dth):dth:hth
    xprt=[tr*cos(q-dth+%pi-%pi/2) , tr*cos(q+%pi-%pi/2)];
    zprt=[tr*sin(q-dth+%pi-%pi/2)+ztc , tr*sin(q+%pi-%pi/2)+ztc];
    plot(xprt,zprt,'--');
end

//地面
xplt=[-0.8,0.8];
zplt=[0,0];
plot(xplt,zplt);





function j=fj(x,u)   //予測と、その評価（台形法）
    q1=x(1);
    q2=x(2);
    q3=x(3);
    q4=x(4);
    q1d=x(5);
    q2d=x(6);
    q3d=x(7);
    q4d=x(8);
    
    Q=[q1;q2;q3;q4];
    Qd=[q1d;q2d;q3d;q4d];
    
    tu2=u(1)
    tu3=u(2)
    tu4=u(3)
    je=0;
    for t=dti:dti:T 
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
    
        H=[h1;h2;h3;h4];

        G1=-( m1*(L1G-L1R)+m2*(L1-L1R)+m3*(L1-L1R)+m4*(L1-L1R) )*sin(q1)*g ;
        G2=(m2*L2G+m3*L2+m4*L2)*sin(q2)*g ;
        G3=(m3*L3G+m4*L3)*sin(q3)*g ;
        G4=m4*L4G*sin(q4)*g;
        G=[G1;G2;G3;G4];
        aaa=inv(M);
        //printf("aaa=%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n\n",aaa(1,1),aaa(1,2),aaa(1,3),aaa(1,4),aaa(2,1),aaa(2,2),aaa(2,3),aaa(2,4),aaa(3,1),aaa(3,2),aaa(3,3),aaa(3,4),aaa(4,1),aaa(4,2),aaa(4,3),aaa(4,4));
        Qdd=inv(M)*([tu2;tu2-tu3;tu3-tu4;tu4]-H-G);
    
        Qd_a=Qd;
        Qd=Qd+Qdd*dti;
        Q=Q+(Qd_a+Qd)/2*dti;
    
        q1=Q(1,1);
        q2=Q(2,1);
        q3=Q(3,1);
        q4=Q(4,1);
        q1d=Qd(1,1);
        q2d=Qd(2,1);
        q3d=Qd(3,1);
        q4d=Qd(4,1);
        
        xr=(q1-q10)*L1R;//足裏接点は必要な際th1iから求める。
        x4a_af=x4a;
        x4a=xr+(L1-L1R)*sin(q1)-L2*sin(q2)-L3*sin(q3)-L4*sin(q4);
        z4a_af=x4a
        z4a=L1R+(L1-L1R)*cos(q1)-L2*cos(q2)-L3*cos(q3)-L4*cos(q4);
        e_af=e;//初期のeを定義しておく必要がある。
        e=( ( (z4a-ztc)^2 + x4a^2 )^0.5 -tr)^2//スイングレッグ先端点と目標起動の最短距離の2乗
        je=je+(e+e_af)/2;//eを台形法で積分
    end
    
    j=ke*je+ktu*(tu2^2+tu3^2+tu4^2);//
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
    tu2i=tu2; tu3i=tu3; tu4i=tu4; //前回の解を初期値にして解を探す

    sq2=0; sq3=0; sq4=0;//初期のsq2_afに0を入れるためのもの

    eta2=eta20; eta3=eta30; eta4=eta40;//勾配法の学習率

    for i=1:200 //最大200回の探索で終了する
        //各変数で評価関数を偏微分（ひとつ前の値を取って置く）
        f=fj(x,[tu2i tu3i tu4i]);
        
        sq2_af = sq2;
        sq2=(fj(x,[tu2i+h tu3i tu4i])-f)/h;

        sq3_af = sq3;
        sq3=(fj(x,[tu2i tu3i+h tu4i])-f)/h;

        sq4_af = sq4;
        sq4=(fj(x,[tu2i tu3i tu4i+h])-f)/h;
  
        printf("%lf %lf %lf   ",sq2,sq3,sq4);
        if (sq2*sq2)<E_sp_q && (sq3*sq3)<E_sp_q && (sq4*sq4)<E_sp_q then        //傾きがE_spより小さくなったら勾配法終了
            break;
        end
        
        if sq2*sq2_af < 0 then//極値をこえる際に学習率を下げる
            eta2=eta2/2;
        end

        if sq3*sq3_af < 0 then//極値をこえる際に学習率を下げる
            eta3=eta3/2;
        end

        if sq4*sq4_af < 0 then//極値をこえる際に学習率を下げる
            eta4=eta4/2;
        end
        tu2i=tu2i-eta2*sq2;//傾きと学習率に従って、測定点をずらす
        tu3i=tu3i-eta3*sq3;
        tu4i=tu4i-eta4*sq4;
    end
    //disp(f_t);
    //if (left<a) && (a<right) && (minyY>f_t) then   //範囲外の解は入れない.かつ、いままで見つけた極小値で最も小さい      //上限をなくしました
        

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
    
    
    ans(1)=tu2i;
    ans(2)=tu3i;
    ans(3)=tu4i;
endfunction

x0=[q10;q30;0;0];
t=0:0.001:1.0;              // 評価時間　1[ms]ごと
t0=0;                       // 時刻の初期値

q1=q10;
q2=q20;
q3=q30;
q4=q40;
q1d=0;
q2d=0;
q3d=0;
q4d=0;

Q=[q1;q2;q3;q4];
Qd=[q1d;q2d;q3d;q4d];

XQ=[Q];//角度の履歴を保存する行列（初期状態を1列目に代入）
XQd=[Qd];//角速度の履歴を保存する行列（初期状態を1列目に代入）
Xtu2=[];Xtu3=[];Xtu4=[];//tu1,tu2,tu3の履歴を保存する行列

tt=dt:dt:endT//初めは1ステップ目から計算し、終了時間まで計算をするのでdt秒からdt秒おきにendTまで計算する
//シミュレーション
for t=tt
    printf("t=%lf   ",t);
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

    H=[h1;h2;h3;h4];

    G1=-( m1*(L1G-L1R)+m2*(L1-L1R)+m3*(L1-L1R)+m4*(L1-L1R) )*sin(q1)*g ;
    G2=(m2*L2G+m3*L2+m4*L2)*sin(q2)*g ;
    G3=(m3*L3G+m4*L3)*sin(q3)*g ;
    G4=m4*L4G*sin(q4)*g;
    G=[G1;G2;G3;G4];





    //評価関数内の積分の初期値を求めておく
    xr=(q1-q10)*L1R;//足裏接点は必要な際th1iから求める。
    x4a=xr+(L1-L1R)*sin(q1)-L2*sin(q2)-L3*sin(q3)-L4*sin(q4);
    z4a=L1R+(L1-L1R)*cos(q1)-L2*cos(q2)-L3*cos(q3)-L4*cos(q4);
    e=( ( (z4a-ztc)^2 + x4a^2 )^0.5 -tr)^2
    

    A=gradent(minu,maxu,[q1;q2;q3;q4;q1d;q2d;q3d;q4d]);
    tu2=A(1);
    tu3=A(2);
    tu4=A(3);
    Xtu2=[Xtu2 tu2];
    Xtu3=[Xtu3 tu3];
    Xtu4=[Xtu4 tu4];
    

    Qdd=inv(M)*([tu2;tu2-tu3;tu3-tu4;tu4]-H-G);

    Qd_a=Qd;
    Qd=Qd+Qdd*dti;
    Q=Q+(Qd_a+Qd)/2*dti;

    q1=Q(1,1);
    q2=Q(2,1);
    q3=Q(3,1);
    q4=Q(4,1);
    q1d=Qd(1,1);
    q2d=Qd(2,1);
    q3d=Qd(3,1);
    q4d=Qd(4,1);

    XQ=[XQ Q];//角度の履歴を保存
    XQd=[XQd Qd];//角速度の履歴を保存
    printf("\n");
end


//リンク描写
idx=1:lookT:int(endT/dt)+1;				// プロット 20[ms]=0.02[s]ごと
xq1=XQ(1,idx)';//行、列を指定して入れる
xq2=XQ(2,idx)';
xq3=XQ(3,idx)';
xq4=XQ(4,idx)';
count=1;//色変更のための変数
for i=1:int(endT/dt)/lookT+1
    q1=xq1(i);
    q2=xq2(i);
    q3=xq3(i);
    q4=xq4(i);

    //リンク描写
            
            xr=(q1-q10)*L1R;
            x1=xr+(L1-L1R)*sin(q1);
            x2=xr+(L1-L1R)*sin(q1)-L2*sin(q2);
            x3=xr+(L1-L1R)*sin(q1)-L2*sin(q2)-L3*sin(q3);
            x4=xr+(L1-L1R)*sin(q1)-L2*sin(q2)-L3*sin(q3)-L4*sin(q4);

            z1=L1R+(L1-L1R)*cos(q1);
            z2=L1R+(L1-L1R)*cos(q1)-L2*cos(q2);
            z3=L1R+(L1-L1R)*cos(q1)-L2*cos(q2)-L3*cos(q3);
            z4=L1R+(L1-L1R)*cos(q1)-L2*cos(q2)-L3*cos(q3)-L4*cos(q4);

            //円弧とリンク１の交点を求める
            xonc=L1R*sin(q1-%pi)+xr;
            yonc=L1R*cos(q1-%pi)+L1R;
            
            xprt=[xonc,x1];
            zprt=[yonc,z1];
            plot2d(xprt,zprt,count);
            xprt=[x1,x2];
            zprt=[z1,z2];
            plot2d(xprt,zprt,count);
            xprt=[x2,x3];
            zprt=[z2,z3];
            plot2d(xprt,zprt,count);
            xprt=[x3,x4];
            zprt=[z3,z4];
            plot2d(xprt,zprt,count);
            
            //リンク1を基準に円弧を描く
            hth=30/180*%pi;
            dth=5/180*%pi;
            for q=-(hth-dth):dth:hth
                xprt=[L1R*cos(q-dth-q1-%pi/2)+xr,L1R*cos(q-q1-%pi/2)+xr];
                zprt=[L1R*sin(q-dth-q1-%pi/2)+L1R,L1R*sin(q-q1-%pi/2)+L1R];
                plot2d(xprt,zprt,count);
            end
            
            //重心描写
            xg1=xr+(L1G-L1R)*sin(q1);
            xg2=xr+(L1-L1R)*sin(q1)-L2G*sin(q2);
            xg3=xr+(L1-L1R)*sin(q1)-L2*sin(q2)-L3G*sin(q3);
            xg4=xr+(L1-L1R)*sin(q1)-L2*sin(q2)-L3*sin(q3)-L4G*sin(q4);
            zg1=L1R+(L1G-L1R)*cos(q1);
            zg2=L1R+(L1-L1R)*cos(q1)-L2G*cos(q2);
            zg3=L1R+(L1-L1R)*cos(q1)-L2*cos(q2)-L3G*cos(q3);
            zg4=L1R+(L1-L1R)*cos(q1)-L2*cos(q2)-L3*cos(q3)-L4G*cos(q4);
 
            xm=(m1*xg1+m2*xg2+m3*xg3+m4*xg4)/msum;
            zm=(m1*zg1+m2*zg2+m3*zg3+m4*zg4)/msum;
            xplt=[xm,xm];
            zplt=[0.005+zm,-0.005+zm];
            plot2d(xplt,zplt,count);
            xplt=[0.005+xm,-0.005+xm];
            zplt=[zm,zm];
            plot2d(xplt,zplt,count);

            count=count+1;//色変更
end

Xtu2=[Xtu2 0];
Xtu3=[Xtu3 0];
Xtu4=[Xtu4 0];

tt=[0 tt];//グラフは時刻0から描くので0を挿入
index=1:endT/dt+1;
subplot(222);   title('トルク　tu2（黄）,tu3（緑）,tu4（青）'); xlabel('t[s]'); ylabel('τ[Nm]');
                plot2d(tt,Xtu2,7);
                plot2d(tt,Xtu3,3);
                plot2d(tt,Xtu4,2);
                xgrid();
                

subplot(223);   title('角度　q1（赤）,q2（黄）,q3（緑）,q4（青）'); xlabel('t[s]'); ylabel('θ[rad]');
                plot2d(tt,XQ(1,index),5); 
                plot2d(tt,XQ(2,index),7);
                plot2d(tt,XQ(3,index),3);
                plot2d(tt,XQ(4,index),2);
                xgrid();
subplot(224);   title('角速度　qd1（赤）,qd2（黄）,qd3（緑）,qd4（青）'); xlabel('t[s]'); ylabel('θdot[rad/s]');
                plot2d(tt,XQd(1,index),5); 
                plot2d(tt,XQd(2,index),7);
                plot2d(tt,XQd(2,index),3);
                plot2d(tt,XQd(3,index),2);
                xgrid();
