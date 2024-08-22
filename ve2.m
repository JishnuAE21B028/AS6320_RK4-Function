% AS6320 Assignment 3 - Rijke tube model eq analysis
% Aswin Balaji T (AE20B016)
% Submitted on 13-05-2023

%% Parameters
K_b = linspace(0.0001,0.1,20); % range of control parameter defined
tau_b = linspace(0.45,0.1,20); % range of time delay defined

%% Velocity and Energy plots
counter = 3;
for ct = 1:2
    tt = 0.1; % initial condition
    K = K_b(ct);
    Bnn=zeros(20,20);
    bnl=zeros(20,20);
    Xo=zeros(20,1);
    Xo(1,1) = tt; % initial condition
    uo=0.5; % mean flow velocity
    r=1.4; % gamma
    co=399.6; 
    xf=0.29; % flame location
    w=@(J) J*pi;
    Beta=@(J) sqrt(3)/(r*uo/co)*K*J*pi*sin(J*pi*xf);
    A1=zeros(20,20);
    A2=zeros(20,20);
    A3=zeros(20,20);
    Bmat=zeros(20,1);
    Umat=zeros(20,1);
    Pmat=zeros(20,1);

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
    zeta=@(J)(c1*(w(J)/pi)+c2*sqrt(pi/w(J)))/(2*pi); % damping term is defined here
    j=2;
    J=1;
    for i=2:2:20
        A1(i,j)=2*J*pi*zeta(J);
        J=J+1;
        j=j+2;
    end

    A2=Bmat*Umat';

    A3=Bmat*Pmat';

    Bnn=A1-A2+0.45.*A3;
    u=[];
    u(1)=0;
    i=2;
    Xi=Xo;
    T=[];
    T(1)=0;
    for t=linspace(0.01,40,200)
        T(i)=t;
        Xs=Xi;
        if((t-0.45)<0)
            Udtua= @(X) Umat'*X-0.45.*Pmat'*X;
            f=@(t,X) -Bnn*X;
        else
            Udtua=@(X) Umat'*X-0.45.*Pmat'*X;
            Bnl=@(X) 3/4.*Udtua(Xs).*A2-3/4*0.45.*Udtua(Xs).*A3;
            f=@(t,X) -Bnn*X-Bnl(Xs)*X;
        end

        n=(t-0)/0.001;
        for k=1:n
            k1=f(t,Xs);                                                              
            k2=f(t+0.5*0.001,Xs+0.5*0.001*k1);
            k3=f(t+0.5*0.001,Xs+0.5*0.001*k2);
            k4=f(t+0.001,Xs+0.001*k3);
            k=(0.001/6).*(k1+2.*k2+2.*k3+k4);
            Xs=Xs+k;
            t=t+0.001;
        end
        u(i)=Xs(1);  
        i=i+1;
    end
    figure(counter);
    plot(T,u,'LineWidth',1);
    xlabel('Time (t)');
    ylabel("u'");
    ylim([-0.3, 0.3]);
    counter = counter + 1;
end
