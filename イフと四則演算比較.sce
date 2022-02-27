etai=0.1;
sp=0.4;
sp_af=-0.4;
/*
for i=1:1000000
    etai=etai+sp*sp_af-((sp^2)^0.5)*((sp_af^2)^0.5);
end
*/

for i=1:1000000
    if sp*sp_af<0 then
        etai=etai/2
    end
end



