function assg = bifurcation(K_b,tau_b)
BNN=zeros(20,20);
BNL=zeros(20,20);
X_0=zeros(20,1);
X_0(1,1) = 0.2; %inital condition
u_0=0.5; %mean flow velocity 
r=1.4; %gamma
c_0=399.6;
ZZ = [];
for qq = 1:length(tau_b)
    for zz = 1:length(K_b)
        K = K_b(zz);
        xf=0.29; %flame location
        tau= tau_b(zz); %time delay during each iteration
        w=@(J) J*pi;
        Beta=@(J) sqrt(3)/(r*u_0/c_0)*K*J*pi*sin(J*pi*xf);
        A1=zeros(20,20); A2=zeros(20,20); A3=zeros(20,20);
        Bmat=zeros(20,1); Umat=zeros(20,1); Pmat=zeros(20,1); %initializing the matrices

        J=1;
        for i=1:2:20
            Umat(i,1)=cos(J*pi*xf);
            J=J+1;
        end

        J=1;
        for i=2:2:20
            Pmat(i,1)=cos(J*pi*xf);
            J=J+1;
        end

        J=1;
        for i=2:2:20
            Bmat(i,1)=-Beta(J);
            J=J+1;
        end

        j=2;
        for i=1:2:20
            A1(i,j)=-1;
            j=j+2;
        end
        j=1;
        J=1;
        for i=2:2:20
            A1(i,j)=(w(J)^2);
            J=J+1;
            j=j+2;
        end
        c1=0.23;
        c2=0.06;
        zeta=@(J)(c1*(w(J)/pi)+c2*sqrt(pi/w(J)))/(2*pi);
        j=2;
        J=1;
        for i=2:2:20
            A1(i,j)=2*J*pi*zeta(J);
            J=J+1;
            j=j+2;
        end

        A2=Bmat*Umat';
        A3=Bmat*Pmat';

        BNN=A1-A2+tau.*A3;

        u=[];
        u(1)=0;
        i=2;
        Xi=X_0;
        T=[];
        T(1)=0;
        for t=linspace(0.01,40,250) %the time interval
            T(i)=t;
            X=rk43(tau,BNN,T(i),T(i-1),0.001,Xi,A2,A3,Umat,Pmat);     
            Xi=X;
            n1(i)=X(1,1);
            n2(i)=X(3,1);
            i=i+1;
        end

        R = n1(150:200); %considering the last 50 data points to calucalte the |U_1|
        ZZ(qq,zz) = max(R) - abs(min(R)); %difference between the max and min value of nita_1
    end
end


%% Hystheris - control paramter is varied backwards

BNN=zeros(20,20);
BNL=zeros(20,20);
X_0=zeros(20,1);
X_0(1,1) = 2.01;
ZZ_h = [];
K_b = fliplr(K_b); %reversing the control paramter
tau_b = fliplr(tau_b);%reversing the time delay 
for qq = 1:length(tau_b)
    for zz = 1:length(K_b)
        K = K_b(zz);

        xf=0.29;
        tau= tau_b(zz);
        w=@(J) J*pi;
        Beta=@(J) sqrt(3)/(r*u_0/c_0)*K*J*pi*sin(J*pi*xf);
        A1=zeros(20,20); A2=zeros(20,20); A3=zeros(20,20);Bmat=zeros(20,1);
        Umat=zeros(20,1);Pmat=zeros(20,1);


        J=1;
        for i=1:2:20
            Umat(i,1)=cos(J*pi*xf);
            J=J+1;
        end

        J=1;
        for i=2:2:20
            Pmat(i,1)=cos(J*pi*xf);
            J=J+1;
        end

        J=1;
        for i=2:2:20
            Bmat(i,1)=-Beta(J);
            J=J+1;
        end

        j=2;
        for i=1:2:20
            A1(i,j)=-1;
            j=j+2;
        end
        j=1;
        J=1;
        for i=2:2:20
            A1(i,j)=(w(J)^2);
            J=J+1;
            j=j+2;
        end
        c1=0.23;
        c2=0.06;
        zeta=@(J)(c1*(w(J)/pi)+c2*sqrt(pi/w(J)))/(2*pi); %damping term is defined here
        j=2;
        J=1;
        for i=2:2:20
            A1(i,j)=2*J*pi*zeta(J);
            J=J+1;
            j=j+2;
        end


        A2=Bmat*Umat';
        A3=Bmat*Pmat';
        BNN=A1-A2+tau.*A3;
        u=[];
        u(1)=0;
        i=2;
        Xi=X_0;
        T=[];
        T(1)=0;
        for t=linspace(0.01,40,200)
            T(i)=t;
            X=rk43(tau,BNN,T(i),T(i-1),0.001,Xi,A2,A3,Umat,Pmat);
            Xi=X;
            n1(i)=X(1,1);
            n2(i)=X(3,1);
            i=i+1;
        end

        R = n1(150:200);%considering the last 50 data points to calucalte the |U_1|
        ZZ_h(qq,zz) = max(R) - abs(min(R)); %difference between the max and min U_1 during each iteration
    end
end
%% Plots
figure(1)
plot(fliplr(K_b),ZZ(1,:),'o')
hold on
plot(K_b,ZZ_h(1,:),'ro')

% Bonus plot - using surf function and meshgrid
figure(2)
[X,Y] = meshgrid(fliplr(tau_b),fliplr(K_b));
surf(X,Y,ZZ)
[X,Y] = meshgrid((tau_b),(K_b));
CO(:,:,1) = zeros(20); % red
CO(:,:,2) = ones(20).*linspace(0.5,0.6,20); % green
CO(:,:,3) = ones(20).*linspace(0,1,20); % blue
hold on
surf(X,Y,ZZ_h,CO)

%% RK4 Function
    function Xs=rk43(tau,Bnn,T,to,h,xo,A2,A3,Umat,Pmat)
        if((T-tau)<0)
            Udtua= @(X) Umat'*X-tau.*Pmat'*X;
            f=@(t,X) -Bnn*X;
        end
        if((T-tau)>=0)
            Udtua=@(X) Umat'*X-tau.*Pmat'*X;
            Bnl=@(X) 3/4.*Udtua(xo).*A2-3/4*tau.*Udtua(xo).*A3;
            f=@(t,X) -Bnn*X-Bnl(X)*X;
        end

        n=(T-to)/h;
        for i=1:n
            k1=f(to,xo);
            k2=f(to+h/2,xo+0.5.*h.*k1);
            k3=f(to+h/2,xo+0.5.*h.*k2);
            k4=f(to+h,xo+h.*k3);
            k=(h/6).*(k1+2.*k2+2.*k3+k4);
            Xs=xo+k;
            xo=Xs;
            U=Udtua(Xs);
            to=to+h;
        end
    end

end