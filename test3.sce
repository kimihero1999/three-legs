

etai=0.1;
sp=0.4;
sp_af=-0.4;
for i=1:100000
    etai=etai+sp*sp_af-((sp^2)^0.5)*((sp_af^2)^0.5);
end





