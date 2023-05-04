g=2;
tol=0.001;
max1=50;
a=1;
b=1;
h1=1/4; %5x5
h2=1/8; %9x9
U1=dirich(@f1,@f2,@f3,@f4,a,b,h1,tol,max1);
U2=dirich(@f1,@f2,@f3,@f4,a,b,h2,tol,max1);
[Y1,X1]=meshgrid(0:h1:1, 0:h1:1);
subplot(1,2,1)
surf(X1,Y1,U1)
hold on
[Y2,X2]=meshgrid(0:h2:1, 0:h2:1);
subplot(1,2,2)
surf(X2,Y2,U2)

function F1=f1(x)
F1 = x.^2;
end
function F2=f2(x)
F2 = (x-1).^2;
end
function F3=f3(y)
F3 = y.^2;
end
function F4=f4(y)
F4 = (y-1).^2;
end

function U=dirich(f1,f2,f3,f4,a,b,h,tol,max1)
n=fix(a/h)+1;
m=fix(b/h)+1;
ave=(a*(feval(f1,0)+feval(f2,0))+b*(feval(f3,0)+feval(f4,0)))/(2*a+2*b);
U=ave*ones(n,m);
U(1,1:m)=feval(f3,0:h:(m-1)*h)';
U(n,1:m)=feval(f4,0:h:(m-1)*h)';
U(1:n,1)=feval(f1,0:h:(n-1)*h);
U(1:n,m)=feval(f2,0:h:(n-1)*h);
U(1,1)=(U(1,2)+U(2,1))/2;
U(1,m)=(U(1,m-1)+U(2,m))/2;
U(n,1)=(U(n-1,1)+U(n,2))/2;
U(n,m)=(U(n-1,m)+U(n,m-1))/2;
w=4/(2+sqrt(4-(cos(pi/(n-1))+cos(pi/(m-1)))^2));
err=1;
cnt=0;
while((err>tol)&(cnt<=max1))
    err=0;
    for j=2:m-1
        for i=2:n-1
            relx=w*(U(i+1,j)+U(i-1,j)+U(i,j+1)+U(i,j-1)-4*U(i,j)-h^2*2)/4;
            U(i,j)=U(i,j)+relx;
            if (err<=abs(relx))
                err=abs(relx);
            end
        end
    end
    cnt=cnt+1;
end
%U=flipud(U');
end
