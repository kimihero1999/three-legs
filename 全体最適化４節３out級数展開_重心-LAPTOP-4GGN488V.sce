//3足歩行ロボット側面モデル 4節３出力 級数展開

//すべての時間で最適化する

clear ;

id=4;//figureの番号
figure(id);

Sth=75/180*%pi;//股傾斜角度
L4=0.4; L3=0.4; L1=(L3+L4)*sin(Sth); L2=0.2;//各関節長（スタンスレッグの股傾斜）
Lg=(L3*cos(Sth)+L2/tan(%pi/3))*tan(%pi/3);//初期位置での原点から足先端までの距離
L1G=L1/2 ; L2G=L2/2 ; L3G=L3/2 ; L4G=L4/2;//各関節の重心位置までの距離
L1R=L1*0.05;//スタンスレッグの足裏半径

g=9.80665 ;
m4=0.1; m3=0.1; m2=2; m1=(m3+m4)*2; msum=m1+m2+m3+m4;//各関節質量と全体重量
I1=m1*L1*L1/12 ; I2=m2*L2*L2/12 ; I3=m3*L3*L3/12 ;I4=m4*L4*L4/12;//各関節の慣性モーメント
xr0=0;//足裏接点初期値
th10=6.037/180*%pi; th20=33.747/180*%pi; th30=0/180*%pi; th40=0/180*%pi;//初期角度
q10=th10; q20=q10+th20; q30=q20+th30; q40=q30+th40;//初期絶対角度角度

//学習率
eta000=0.000001;
eta10=eta000;
eta20=eta000;
eta30=eta000;
eta00=[eta10 eta10 eta10 eta10 eta10 eta10 eta10;eta20 eta20 eta20 eta20 eta20 eta20 eta20;eta30 eta30 eta30 eta30 eta30 eta30 eta30];

E_sq=0.0000001; E_sq_q=E_sq^2;//勾配法の際、探索を終了する傾き、その2乗

//A0=[0 0 0 0 0 0 0;0 0 0 0 0 0 0;0 0 0 0 0 0 0];//初期値


A0=[00000.2144 00000.1688 00000.4671 00000.2450 00000.1039 00000.0889 -0000.2974;00000.5630 00000.1104 00000.0009 -0000.0016 00000.0053 00000.1579 00000.1149;00000.0048 00000.1423 00000.0369 00000.0466 00000.0312 00000.0287 00000.0140];//初期値+-型


h=0.0000000001;//傾きを計算する際の極小の値

//ゲイン一覧
ke=1;   kfp=2000; kbk=10000; kmt=10;
//ken=0.00000000000001;
ken=0;
right=200;//トルクの上限[n*m](現在無効)
left=200;//トルクの下限[n*m](現在無効)

ztc=-3.557;//目標軌跡円の中心点（x=0とする）


x40=(L1 -L1R)*sin(q10)-L2 *sin(q20)-L3 *sin(q30)-L4*sin(q40);
z40=L1R+(L1-L1R)*cos(q10)-L2*cos(q20)-L3*cos(q30)-L4*cos(q40);
tr=( (z40-ztc)^2 + x40^2 )^0.5;//目標軌跡円の半径（必ず目標軌跡状に足先端がくるようになる）

//目標軌跡円を描写
hth=10/180*%pi;//円の幅
dth=1/180*%pi;//線の細かさ（角度）
for q=-(hth-dth):dth:hth
    xprt=[tr*cos(q-dth+%pi-%pi/2) , tr*cos(q+%pi-%pi/2)];
    zprt=[tr*sin(q-dth+%pi-%pi/2)+ztc , tr*sin(q+%pi-%pi/2)+ztc];
    plot(xprt,zprt,'--');
end

xr0=0;
//重心
xg1=xr0+(L1G-L1R)*sin(q10);
xg2=xr0+(L1-L1R)*sin(q10)-L2G*sin(q20);
xg3=xr0+(L1-L1R)*sin(q10)-L2*sin(q20)-L3G*sin(q30);
xg4=xr0+(L1-L1R)*sin(q10)-L2*sin(q20)-L3*sin(q30)-L4G*sin(q40);
zg1=L1R+(L1G-L1R)*cos(q10);
zg2=L1R+(L1-L1R)*cos(q10)-L2G*cos(q20);
zg3=L1R+(L1-L1R)*cos(q10)-L2*cos(q20)-L3G*cos(q30);
zg4=L1R+(L1-L1R)*cos(q10)-L2*cos(q20)-L3*cos(q30)-L4G*cos(q40);

xm0=(m1*xg1+m2*xg2+m3*xg3+m4*xg4)/msum;
zm0=(m1*zg1+m2*zg2+m3*zg3+m4*zg4)/msum;

mztc=zm0;//重心目標軌道の中心点（x=0とする）
mtr=xm0;//重心目標軌道の半径

//目標軌跡円を描写
hth=100/180*%pi;//円の幅
dth=2/180*%pi;//線の細かさ（角度）
for q=-(hth-dth):dth:hth
    xprt=[mtr*cos(q-dth+%pi-%pi/2) , mtr*cos(q+%pi-%pi/2)];
    zprt=[mtr*sin(q-dth+%pi-%pi/2)+mztc , mtr*sin(q+%pi-%pi/2)+mztc];
    plot(xprt,zprt,'--');
end


//地面を描写
xplt=[-0.8,0.8];
zplt=[0,0];
plot(xplt,zplt);

lookT=20;

FT=1;
dt=0.001;
TT=dt:dt:FT;//初めは1ステップ目から計算し、終了時間まで計算をするのでdt秒からdt秒おきにendTまで計算する
tu_i=1:FT/dt; //tuのindexとして使うため、TTと同じ数のindexを用意

//関数内では記録をつけずに評価に必要な値のみ計算していく

//周期1sからの正弦波を作っておく。
TTp=TT*%pi*2;
TTp2=TTp*2;
TTp4=TTp2*2;
TTp_s=sin(TTp);
TTp2_s=sin(TTp2);
TTp4_s=sin(TTp4);
TTp_c=cos(TTp);
TTp2_c=cos(TTp2);
TTp4_c=cos(TTp4);

function j=fj(A)   //3つのトルクの5次曲線の係数の行列Aを受け取り、評価を出力する（評価関数）
    tu2T =A(1,1) + A(1,2)*TTp_s + A(1,3)*TTp_c + A(1,4)*TTp2_s + A(1,5)*TTp2_c + A(1,6)*TTp4_s + A(1,7)*TTp4_c;
    tu3T =A(2,1) + A(2,2)*TTp_s + A(2,3)*TTp_c + A(2,4)*TTp2_s + A(2,5)*TTp2_c + A(2,6)*TTp4_s + A(2,7)*TTp4_c;
    tu4T =A(3,1) + A(3,2)*TTp_s + A(3,3)*TTp_c + A(3,4)*TTp2_s + A(3,5)*TTp2_c + A(3,6)*TTp4_s + A(3,7)*TTp4_c;
    
    q1=q10;q2=q20;q3=q30;q4=q40;
    q1d=0;q2d=0;q3d=0;q4d=0;
    Q=[q10;q20;q30;q40];
    Qd=[0;0;0;0];
    
    //評価関数内の積分の初期値
    x4a=x40;
    e=0;
    mt=0;
    
    
    endT=FT;//歩行が終了した時間を格納する変数
    je=0;
    en=0;
    bakc=0;
    jmt=0;
    
    for i=tu_i;
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
        Qdd=inv(M)*([tu2T(i);tu2T(i)-tu3T(i);tu3T(i)-tu4T(i);tu4T(i)]-H-G);
    
        Qd_a=Qd;
        Qd=Qd+Qdd*dt;
        Q=Q+(Qd_a+Qd)/2*dt;
    
        q1=Q(1,1);
        q2=Q(2,1);
        q3=Q(3,1);
        q4=Q(4,1);
        q1d=Qd(1,1);
        q2d=Qd(2,1);
        q3d=Qd(3,1);
        q4d=Qd(4,1);
        
        
        
        
        
        xr=(q1-q10)*L1R;//足裏接点は必要な際th1iから求める。
        
        //重心
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
        
        mt_af=mt;//初期のmtを定義しておく必要がある。
        mt=( ( (xm-mztc)^2 + zm^2 )^0.5 -mtr)^2
        jmt=jmt+(mt+mt_af)/2;
        
        x4a_af=x4a;
        x4a=xr+(L1-L1R)*sin(q1)-L2*sin(q2)-L3*sin(q3)-L4*sin(q4);
        z4a=L1R+(L1-L1R)*cos(q1)-L2*cos(q2)-L3*cos(q3)-L4*cos(q4);
        
        if x4a<x4a_af then
            bakc=bakc+(x4a-x4a_af)^2;
        end

        e_af=e;//初期のeを定義しておく必要がある。
        e=( ( (z4a-ztc)^2 + x4a^2 )^0.5 -tr)^2
        je=je+(e+e_af)/2;
        
        en=en+abs(tu2T(i)*q2d)+abs(tu3T(i)*q3d)+abs(tu3T(i)*q3d);
    end
    Je=ke*je;
    printf(" Je=%f ",Je);
    
    Jene=ken*(en/dt-2000000)^2;
    printf(" Jene=%f ",Jene);
    
    Jbk=kbk*bakc;
    printf(" Jbk=%f ",Jbk);

    Jfp=kfp*((x4a+x40)^2+z4a^2);
    printf(" Jfp=%f ",Jfp);
    
    Jmt=kmt*jmt;
    printf(" Jmt=%f ",Jmt);
    
    j=Je+Jene+Jfp+Jbk+Jmt;
    printf(" x4a=%f||",x4a);
    
    
    
    printf(" E=%f||",en/dt);
endfunction



function ans=gradent(left,right)//勾配法で評価関数から最適解を求める
    eta=eta00;//勾配法の学習率
    /*
    //解の上限と下限を決めるプログラム
    f_left=fj(x,left);
    f_right=fj(x,right);
    minyU=left;
    minyY=f_left
    if f_left>f_right then
        minyU=right;
        minyY=f_right;
    end
   */

    sq=[0 0 0 0 0 0 0;0 0 0 0 0 0 0;0 0 0 0 0 0 0];//初期のsq_afに0を入れるためのもの
    
    
    AA=A0;
    for ii=1:50 //最大200回の探索で終了する
        
        //各変数で評価関数を偏微分（ひとつ前の値を取って置く）
        f=fj(AA);
        sq_af=sq;
        for i=1:3
            for j=1:7
                AAh=AA;
                AAh(i,j)=AAh(i,j)+h;
                sq(i,j)=(fj(AAh)-f)/h;
            end
        end
        printf("\n%4d ",ii);
        for i=1:3
            for j=1:7
                printf("%+6.5f ",sq(i,j));
            end
            printf(";");
        end
        
        /*
        if (sq(1,1)^2)<E_sq_q && (sq(1,2)^2)<E_sq_q && (sq(1,3)^2)<E_sq_q && (sq(1,4)^2)<E_sq_q && (sq(1,5)^2)<E_sq_q && (sq(1,1)^2)<E_sq_q && (sq(2,2)^2)<E_sq_q && (sq(2,3)^2)<E_sq_q && (sq(2,4)^2)<E_sq_q && (sq(2,5)^2)<E_sq_q && (sq(3,1)^2)<E_sq_q && (sq(3,2)^2)<E_sq_q && (sq(3,3)^2)<E_sq_q && (sq(3,4)^2)<E_sq_q && (sq(3,5)^2)<E_sq_q then
            break;//傾きがE_sqより小さくなったら勾配法終了
        end
        */

        
        //極値をこえる際に学習率を下げる
        for i=1:3
            for j=1:7
                if sq(i,j)*sq_af(i,j) < 0 then
                    eta(i,j)=eta(i,j)*2/3;
                else
                    eta(i,j)=eta00(i,j);
                end
            end
        end
        
        
        AA=AA-eta.*sq/(sum(sq.^2))^0.5*f;//傾きと学習率に従って、測定点をずらす
    end

    printf("\nsq=");disp(sq);

    printf("\nAA=");disp(AA);
    
    printf("\n");
    for i=1:3
        for j=1:7
            printf("%010.4f ",AA(i,j));
        end
        printf(";");
    end
    
    ans=AA;
endfunction


AE=gradent(left,right);

tu2T =AE(1,1) + AE(1,2)*TTp_s + AE(1,3)*TTp_c + AE(1,4)*TTp2_s + AE(1,5)*TTp2_c + AE(1,6)*TTp4_s + AE(1,7)*TTp4_c;
tu3T =AE(2,1) + AE(2,2)*TTp_s + AE(2,3)*TTp_c + AE(2,4)*TTp2_s + AE(2,5)*TTp2_c + AE(2,6)*TTp4_s + AE(2,7)*TTp4_c;
tu4T =AE(3,1) + AE(3,2)*TTp_s + AE(3,3)*TTp_c + AE(3,4)*TTp2_s + AE(3,5)*TTp2_c + AE(3,6)*TTp4_s + AE(3,7)*TTp4_c;


q1=q10;q2=q20;q3=q30;q4=q40;
q1d=0;q2d=0;q3d=0;q4d=0;
Q=[q10;q20;q30;q40];
Qd=[0;0;0;0];

XQ=[Q];//角度の履歴を保存する行列（初期状態を1列目に代入）
XQd=[Qd];//角速度の履歴を保存する行列（初期状態を1列目に代入）

for i=tu_i;
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
    
    Qdd=inv(M)*([tu2T(i);tu2T(i)-tu3T(i);tu3T(i)-tu4T(i);tu4T(i)]-H-G);
    
    Qd_a=Qd;
    Qd=Qd+Qdd*dt;
    Q=Q+(Qd_a+Qd)/2*dt;

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
end


//リンク描写
idx=1:lookT:int(FT/dt)+1;				// プロット 20[ms]=0.02[s]ごと
xq1=XQ(1,idx)';//行、列を指定して入れる
xq2=XQ(2,idx)';
xq3=XQ(3,idx)';
xq4=XQ(4,idx)';
count=1;//色変更のための変数
for i=1:int(FT/dt)/lookT+1
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
    zplt=[0.01+zm,-0.01+zm];
    plot2d(xplt,zplt,count);
    xplt=[0.01+xm,-0.01+xm];
    zplt=[zm,zm];
    plot2d(xplt,zplt,count);

    count=count+1;//色変更
end
figure(id+100);

tu2T=[tu2T 0];
tu3T=[tu3T 0];
tu4T=[tu4T 0];

tt=[0 TT];//グラフは時刻0から描くので0を挿入
index=1:FT/dt+1;

//リンク先端の速度
Xd=[];
Zd=[];
XZd=[];
for i=index
    q1=XQ(1,i);
    q2=XQ(2,i);
    q3=XQ(3,i);
    q4=XQ(4,i);
    q1d=XQd(1,i);
    q2d=XQd(2,i);
    q3d=XQd(3,i);
    q4d=XQd(4,i);
    
    Xdi=[L1R+(L1-L1R)*cos(q1) -L2*cos(q2) -L3*cos(q3) -L4*cos(q4)]*[q1d;q2d;q3d;q4d];
    Zdi= [-(L1-L1R)*sin(q1) L2*sin(q2) L3*sin(q3) L4*sin(q4)]*[q1d;q2d;q3d;q4d];
    XZdi=(Xdi^2+Zdi^2)^0.5;
    Xd=[Xd Xdi];
    Zd=[Zd Zdi];
    XZd=[XZd XZdi];
end

//重心の速度
xmd=[0];
zmd=[0];
xr=0;
xg1=xr+(L1G-L1R)*sin(q10);
xg2=xr+(L1-L1R)*sin(q10)-L2G*sin(q20);
xg3=xr+(L1-L1R)*sin(q10)-L2*sin(q20)-L3G*sin(q30);
xg4=xr+(L1-L1R)*sin(q10)-L2*sin(q20)-L3*sin(q30)-L4G*sin(q40);
zg1=L1R+(L1G-L1R)*cos(q10);
zg2=L1R+(L1-L1R)*cos(q10)-L2G*cos(q20);
zg3=L1R+(L1-L1R)*cos(q10)-L2*cos(q20)-L3G*cos(q30);
zg4=L1R+(L1-L1R)*cos(q10)-L2*cos(q20)-L3*cos(q30)-L4G*cos(q40);
xm_af=(m1*xg1+m2*xg2+m3*xg3+m4*xg4)/msum;
zm_af=(m1*zg1+m2*zg2+m3*zg3+m4*zg4)/msum;
for i=2:FT/dt+1
    q1=XQ(1,i);
    q2=XQ(2,i);
    q3=XQ(3,i);
    q4=XQ(4,i);
    
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
    xmd=[xmd (xm-xm_af)/dt];
    zmd=[zmd (zm-zm_af)/dt];
    xm_af=xm;
    zm_af=zm;
end
xzmd=(xmd.^2+zmd.^2).^0.5;

subplot(211);title('先端速度 Xd(赤),Zd(青),XZd(黒)'); xlabel('t[s]'); ylabel('v[m/s]');
plot2d(tt,Xd,5);//赤
plot2d(tt,Zd,2);//青
plot2d(tt,XZd,1);//黒
/*
for i=lookT*dt:lookT*dt:FT//リンクを描写した時間に縦線を引く
    xplt=[i,i];
    zplt=[0.3,-0.3];
    plot(xplt,zplt,':b');
end
*/
//xgrid();

subplot(212);title('重心速度 XMd(ピンク),ZMd(水色),XZMd(オレンジ)'); xlabel('t[s]'); ylabel('v[m/s]');
plot2d(tt,xmd,6);//ピンク
plot2d(tt,zmd,4);//水色
plot2d(tt,xzmd,32);//オレンジ
/*
for i=lookT*dt:lookT*dt:FT//リンクを描写した時間に縦線を引く
    xplt=[i,i];
    zplt=[0.3,-0.3];
    plot(xplt,zplt,':b');
end
*/
//xgrid();
                
                
                
figure(id+200);
subplot(221);   title('トルク tu2（黒）,tu3（緑）,tu4（青）'); xlabel('t[s]'); ylabel('[Nm]');
plot2d(tt,tu2T,1);
plot2d(tt,tu3T,3);
plot2d(tt,tu4T,2);
xgrid();
                
subplot(222);   title('角度 q1（黒）,q2（赤）,q3（緑）,q4（青）'); xlabel('t[s]'); ylabel('θ[rad]');
plot2d(tt,XQ(1,index),1); 
plot2d(tt,XQ(2,index),5);
plot2d(tt,XQ(3,index),3);
plot2d(tt,XQ(4,index),2);
xgrid();
                
subplot(223);   title('角速度 qd1（黒）,qd2（赤）,qd3（緑）,qd4（青）'); xlabel('t[s]'); ylabel('θdot[rad/s]');
plot2d(tt,XQd(1,index),1); 
plot2d(tt,XQd(2,index),5);
plot2d(tt,XQd(3,index),3);
plot2d(tt,XQd(4,index),2);
xgrid();
