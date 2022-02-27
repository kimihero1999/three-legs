
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
                
