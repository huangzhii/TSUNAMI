function F = sample_ecdf(f,x)

F=zeros(1,200001);

t = 0:0.01:2000;

for i = 1 : length(x)-1
    t1 = ceil(x(i)*100);
    t2 = floor(x(i+1)*100);
    if t2<t1
        t2=t1;
    end
    if(t2 < 200001)
        F(t1:t2)=f(i);
    else
        F(t1:200000) = f(i);
        F(200001)=1;
    end
end
F(t2:end)=1;