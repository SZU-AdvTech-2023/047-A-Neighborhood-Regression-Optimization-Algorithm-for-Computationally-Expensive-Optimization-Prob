function f=expensive_benchmark_func(x,problem)
%% Ellipsoid
if problem==1    
   D=length(x);
   A=1:D;           
   f=sum(A.*x.*x);
%% Rosenbrock function
elseif problem==2
        D=length(x);
        s=0;
        for i=1:D-1
            s=s+100*(x(i).^2-x(i+1))^2+(1-x(i))^2;
        end
        f=s;
%% Ackley function
elseif problem==3
       D=length(x);
       s1=x*x';
       s2=sum(cos(2*pi*x));      
       f=(20-20*exp(-0.2*sqrt(1/D*s1)))+(exp(1)-exp(1/D*s2));
%% Griewank function
elseif problem==4
         D=length(x);
         s1=0;s2=1;
         for i=1:D
             s1=s1+x(i).^2;
             s2=s2*cos(x(i)/sqrt(i));
         end
         f=1/4000*s1-s2+1;
         
%% Rastrigin function         
elseif problem==5
        D=length(x);
        s=0;
        for i=1:D
            s=s+(x(i).^2-10*cos(2*pi*x(i))+10);
        end
        f=s;
else    
     error('there are only six funtions, problem value should be set to 1 to 5')
end
             
            
   

