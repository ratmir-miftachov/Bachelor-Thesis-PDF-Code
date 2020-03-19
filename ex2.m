vec = rand(1,1000);

for i=1:1000
   if vec(i)>= 0.5
       vec(i)=1;
   elseif 0.25<vec(i)<0.5
       vec(i)=99;
   else vec(i)=0;
   end
       
end