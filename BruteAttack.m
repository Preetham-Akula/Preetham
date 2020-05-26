function out1 = BruteAttack(z)
a=length(z);
input('Do you want to continue? PRESS ENTER');
for k=1:256
    e=randi([0 1],1,a);
    c=xor(z,e)
end

out1=z;
        
        
