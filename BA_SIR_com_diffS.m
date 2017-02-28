function [t,S,I] = BA_SIR_com_diffS(N,n,tau,gamma,MaxTime,Type)

%ʵ��SIS��BA�ޱ�������ϴ����������ļ�FreeScale������adjacent_matrix�����丳ֵ��G���ɡ�
%�����ļ����ƺ�ע�͵���SIR����ģ�ͣ���������о�һ�¡�
%
% Program_7_7( N, n, tau, gamma, MaxTime, Type)
%      This is the MATLAB version of program 7.7 from page 280 of
% "Modeling Infectious Disease in humans and animals"
% by Keeling & Rohani.
%
% It is an SIR disease spread through a network. Allowed
% network types are 'Random','Spatial','Lattice' and 'SmallWorld'
%
% We assume N individuals, each with an averge of n contacts.

% In this model we define an individual by their status flag:
% Status=1   =>  Susceptible
% Status=2   =>  Infectious
% Status=0   =>  Recovered


% Sets up default parameters if necessary.
if nargin == 0
    N=50;
    n=4;
    tau=1;
    gamma=0.1;
    MaxTime=30;
    Type='ba';
    com_num=10;
end

% Checks all the parameters are valid
CheckGreater(N,0,'Number of individuals N');
CheckGreater(n,0,'Number of neighbours n');
CheckGreater(tau,0,'tau');
CheckGreater(gamma,0,'gamma');
CheckGreater(MaxTime,0,'MaxTime');

%Initialise the Network
% (X,Y) is location, G is the network graph matrix
% this means we use S and I for the number of susceptibles and infecteds

[X,Y,G,N]=Create_Network(N,n,Type);
%��kmeans����10������
[idx,c] = kmeans(G,com_num);
%X  Y�ֱ�Ϊ�������������� 
Status=1+0*X; Status(1)=2;
Rate=0*X; Rate(1)=gamma; Rate(find(G(:,1)))=tau; 
%RateΪ������ Rate(1)Ϊ��һ��Ԫ�أ�����Ϊgamma��
%G(:,1)��ʾȡ��һ�е�����Ԫ�أ���ôfind(G(:,1))�����ҵ���һ���в�Ϊ�������ֵ
t=0; i=1; S=N-1; I=1;R=0;
%���������ͱ�����S״̬��������Rate��Ϊ10
x=0;
while(x<N/10)
tempt=randint(1,1,[1,N]);
Rate(tempt)=10;
x=x+1;
end

% The main iteration
subplot(2,1,1);
%find(G) ����G�еķ���Ԫ�أ��С��С�ֵ
[j,k,s]=find(G);
plot([X(j) X(k)]',[Y(j) Y(k)]','-k');
hold on
Col=[0.7 0.7 0.7; 0 1 0; 1 0 0];
for k=1:N
    H(k)=plot(X(k),Y(k),'ok','MarkerFaceColor','g');
end
set(H(1),'MarkerFaceColor','r');
hold off;
drawnow;


while (t<MaxTime & I(end)>0)
    [step,Rate,Status,e]=Iterate(Rate,Status,G,tau,gamma,idx);
    %Status�Ǹ�����length(find(Status==1))Ϊ������ֵΪ1������
    i=i+1;
    t(i)=t(i-1)+step;
    S(i)=length(find(Status==1));
    I(i)=length(find(Status==2));
    R(i)=length(find(Status==0));
    set(H(e),'MarkerFaceColor',Col(Status(e)+1,:));

%     subplot(2,1,1);
%     [j,k,s]=find(G);
%     plot([X(j) X(k)]',[Y(j) Y(k)]','-k');
%     hold on
%     s=find(Status==0); plot(X(s),Y(s),'ok','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',8);
%     s=find(Status==1); plot(X(s),Y(s),'ok','MarkerFaceColor','g','MarkerSize',8);
%     s=find(Status==2); plot(X(s),Y(s),'ok','MarkerFaceColor','r','MarkerSize',10);
%     hold off;
if(i==60)
    X=[length(find(Status==1)),length(find(Status==2)),length(find(Status==0))];
    colormap(cool);
    pie(X,{'S','I','R'});
end
    subplot(4,1,3);
    plot(t,S,'-g',t,R,'-k');
    ylabel 'Susceptibles and Recovereds'
    subplot(4,1,4);
    plot(t,I,'-r');
    ylabel 'Infecteds'
    xlabel 'Time'
    drawnow;
    [];
end


% Create the Network
function [X,Y,G,N]=Create_Network(N,n,Type);

if n > (N-1)
    error('Impossible to have an average of %d contacts with a population size of %d',n,N);
end

G=sparse(1,1,0,N,N);
%X Y �ֱ�ΪN��Ԫ�ص���������������
X=rand(N,1); Y=rand(N,1);
switch Type
    
%random��������˼�룺�������i��j��ʹ��G(i,j)���Ϊ1���Խڵ�Ķ�*�ڵ���������ѭ����
    case {'Random','random'}
        contacts=0;
        while(contacts<n*N)
            i=ceil(rand(1,1)*N); j=ceil(rand(1,1)*N);
            %rand(1,1)�������1*1��0~1֮�������ceil(rand(1,1)*N)����1~100֮�������
            if i~=j & G(i,j)==0
                contacts=contacts+2;
                G(i,j)=1; G(j,1)=1;
            end
        end
              
    case {'Spatial','spatial'}
        D=(X*ones(1,N) - ones(N,1)*X').^2 + (Y*ones(1,N) - ones(N,1)*Y').^2;
        Prob=exp(-D*5)-rand(N,N); Prob=triu(Prob,1)-1e100*tril(ones(N,N),0);
        [y i]=sort(reshape(Prob,N*N,1));
        p=i(end+[(1-n*N/2):1:0]);
        i=1+mod(p-1,N); j=1+floor((p-1)/N);
        G=sparse([i; j],[j; i],1,N,N);

    case {'Lattice','lattice'}
        if mod(sqrt(N),1)~=0
            warning('N=%d is not a square number, rounding to %d',N,round(sqrt(N)).^2);
            N=round(sqrt(N)).^2;
        end
        [X,Y]=meshgrid([0:(sqrt(N)-1)]/(sqrt(N)-1));
        X=reshape(X,N,1); Y=reshape(Y,N,1);
        D=(X*ones(1,N) - ones(N,1)*X').^2 + (Y*ones(1,N) - ones(N,1)*Y').^2;
        Prob=triu(D,1)+1e100*tril(ones(N,N),0);
        [y i]=sort(reshape(Prob,N*N,1));
        p=i(1:(n*N/2 - 2*sqrt(N)));
        i=1+mod(p-1,N); j=1+floor((p-1)/N);
        G=sparse([i; j],[j; i],1,N,N);

    case {'SmallWorld','Smallworld','smallworld'}
        if mod(sqrt(N),1)~=0
            warning('N=%d is not a square number, rounding to %d',N,round(sqrt(N)).^2);
            N=round(sqrt(N)).^2;
        end
        [X,Y]=meshgrid([0:(sqrt(N)-1)]/(sqrt(N)-1));
        X=reshape(X,N,1); Y=reshape(Y,N,1);
        D=(X*ones(1,N) - ones(N,1)*X').^2 + (Y*ones(1,N) - ones(N,1)*Y').^2;
        Prob=triu(D,1)+1e100*tril(ones(N,N),0);
        [y i]=sort(reshape(Prob,N*N,1));
        p=i(1:(n*N/2 - 2*sqrt(N)));
        i=1+mod(p-1,N); j=1+floor((p-1)/N);
        G=sparse([i; j],[j; i],1,N,N);
        % Already made a lattice, now add R random connections
        R=10;
        for k=1:R
            i=1; j=1;
            while (i==j | G(i,j)==1)
                i=ceil(rand(1,1)*N); j=ceil(rand(1,1)*N);
            end
            G(i,j)=1; G(j,1)=1;
        end
%���±�������BA�ޱ������Ĵ��룬�������ģ�    
    case {'BA','ScaleFree','Scalefree','scalefree','BAscalefree','ba'}
        m0= 3; m= 3;%��ʼ����������
        adjacent_matrix = sparse( m0, m0);% ��ʼ���ڽӾ���
        for i = 1: m0            
            for j = 1:m0
                if j ~= i      %ȥ��ÿ���������γɵĻ�
                    adjacent_matrix(i,j) = 1;%������ʼ�ڽӾ���3��ͬ��ͬ�����ĵ�����
                end
            end
        end
        adjacent_matrix =sparse(adjacent_matrix);%�ڽӾ���ϡ�軯
        node_degree = zeros(1,m0+1);                %��ʼ����Ķ�
        node_degree(2: m0+1) = sum(adjacent_matrix);%�Զ�ά��������չ
        for iter= 4:N
            iter;                                %�ӵ�
            total_degree = 2*m*(iter- 4)+6;%���������д˵�Ķ�֮��
            cum_degree = cumsum(node_degree);%��������е�ĶȾ���
            choose= zeros(1,m);%��ʼ��ѡ�����
    % ѡ����һ�����µ������ӵĶ���
            r1= rand(1)*total_degree;%�����ɵ������ĸ���
            for i= 1:iter-1
                if (r1>=cum_degree(i))&( r1<cum_degree(i+1))%ѡȡ�ȴ�ĵ�
                    choose(1) = i;
                    break
                end
            end
            % ѡ���ڶ������µ������ӵĶ���
            r2= rand(1)*total_degree;
            for i= 1:iter-1
                if (r2>=cum_degree(i))&(r2<cum_degree(i+1))
                    choose(2) = i;
                    break
                end
            end
            while choose(2) == choose(1)%��һ����͵ڶ�������ͬ�Ļ�����������
                r2= rand(1)*total_degree;
                for i= 1:iter-1
                    if (r2>=cum_degree(i))&(r2<cum_degree(i+1))
                        choose(2) = i;
                        break
                    end
                end
            end
            % ѡ�����������µ������ӵĶ���
            r3= rand(1)*total_degree;
            for i= 1:iter-1
                if  (r3>=cum_degree(i))&(r3<cum_degree(i+1))
                    choose(3) = i;
                    break
                end
            end
            while (choose(3)==choose(1))|(choose(3)==choose(2))
                r3= rand(1)*total_degree;
                for i=1:iter-1
                    if (r3>=cum_degree(i))&(r3<cum_degree(i+1))
                        choose(3) = i;
                        break
                    end
                end
            end
            %�µ���������, ���ڽӾ�����и���
            for k = 1:m
                adjacent_matrix(iter,choose(k)) = 1;
                adjacent_matrix(choose(k),iter) = 1;
            end
            node_degree=zeros(1,iter+1);
            node_degree(2:iter+1) = sum(adjacent_matrix);
        end
        G = adjacent_matrix;        
        
    otherwise
        error('Do not recognise network type %s',Type);
end



%Do the Up-Dating.
function [step,Rate,Status,Event]=Iterate(Rate,Status,G,tau,gamma,idx);

%ʹ��Kmeans��������������c~c9�ֱ�����������ͬ�����Ľ�㣬km�������������ı�ǩ
km=idx;
c1=find(km==1);
c2=find(km==2);
c3=find(km==3);
c4=find(km==4);
c5=find(km==5);
c6=find(km==6);
c7=find(km==7);
c8=find(km==8);
c9=find(km==9);
c10=find(km==10);

Sum=sum(Rate); Cum=cumsum(Rate);

Event=min(find(Cum>rand(1,1)*Sum));
%��Ȼ������ģ�����Event�����ʵĵ���ͻȻ����Ԫ���������ҹؼ���Event��֤��ΪRate��0Ԫ�ص������������Ϊ����ҵ�һ����ڵ������Ľڵ�����
Status(Event)=mod(Status(Event)+1,3);
%���Status(Event)����ΪS��=1��ʱ����ô��ֵ2����תcase2���������ΪI��=2������ֵ0����תcase 0���ֱ���ֵ1������ΪS
contacts=find(G(:,Event) & Status==1);
%���I���࣬�Ӵ��������
switch Status(Event)
    case 1
        
    case 2
        switch km(Event)
            case 1
                f1=find(Status(c1)==1);
                f2=find(Status(c1)==2);
                if(f1)             
                    Rate(c1(f1))=Rate(c1(f1))+1;
                end
                if(f2)
                    Rate(c1(f2))=Rate(c1(f2))-0.01;
                end
            case 2
                f1=find(Status(c2)==1);
                f2=find(Status(c2)==2);
                if(f1)             
                    Rate(c2(f1))=Rate(c2(f1))+1;
                end
                if(f2)
                    Rate(c2(f2))=Rate(c2(f2))-0.01;
                end
            case 3
                f1=find(Status(c3)==1);
                f2=find(Status(c3)==2);
                if(f1)             
                    Rate(c3(f1))=Rate(c3(f1))+1;
                end
                if(f2)
                    Rate(c3(f2))=Rate(c3(f2))-0.01;
                end
            case 4
                f1=find(Status(c4)==1);
                f2=find(Status(c4)==2);
                if(f1)             
                    Rate(c4(f1))=Rate(c4(f1))+1;
                end
                if(f2)
                    Rate(c4(f2))=Rate(c4(f2))-0.01;
                end
            case 5
                f1=find(Status(c5)==1);
                f2=find(Status(c5)==2);
                if(f1)             
                    Rate(c5(f1))=Rate(c5(f1))+1;
                end
                if(f2)
                    Rate(c5(f2))=Rate(c5(f2))-0.01;
                end
            case 6
                f1=find(Status(c6)==1);
                f2=find(Status(c6)==2);
                if(f1)             
                    Rate(c6(f1))=Rate(c6(f1))+1;
                end
                if(f2)
                    Rate(c6(f2))=Rate(c6(f2))-0.01;
                end
            case 7
                f1=find(Status(c7)==1);
                f2=find(Status(c7)==2);
                if(f1)             
                    Rate(c7(f1))=Rate(c7(f1))+1;
                end
                if(f2)
                    Rate(c7(f2))=Rate(c7(f2))-0.01;
                end
            case 8
                f1=find(Status(c8)==1);
                f2=find(Status(c8)==2);
                if(f1)             
                    Rate(c8(f1))=Rate(c8(f1))+1;
                end
                if(f2)
                    Rate(c8(f2))=Rate(c8(f2))-0.01;
                end
            case 9
                f1=find(Status(c9)==1);
                f2=find(Status(c9)==2);
                if(f1)             
                    Rate(c9(f1))=Rate(c9(f1))+1;
                end
                if(f2)
                    Rate(c9(f2))=Rate(c9(f2))-0.01;
                end
            case 10
                f1=find(Status(c10)==1);
                f2=find(Status(c10)==2);
                if(f1)             
                    Rate(c10(f1))=Rate(c10(f1))+1;
                end
                if(f2)
                    Rate(c10(f2))=Rate(c10(f2))-0.01;
                end
        end
 %  ����switch�Ĵ��a�������^ c1���и�Ⱦ�ĽY�c�����Nԓ��^�еĽY�c�ĸ�Ⱦ���ʶ�����1     
        Rate(Event)=gamma;
        Rate(contacts)=Rate(contacts)+tau;
        %�����д���ʵ���趨I�Ļָ��ʺ�S�ĸ�Ⱦ��
        
        %��ʱ�Ȳ�Ҫ �����ûʲôӰ��G(Event,:)=0;
    case 0
        switch km(Event)
            case 1
                f1=find(Status(c1)==1);
                f2=find(Status(c1)==2);
                if(f1)             
                    Rate(c1(f1))=Rate(c1(f1))-1;
                end
                if(f2)
                    Rate(c1(f2))=Rate(c1(f2))+0.01;
                end
            case 2
                f1=find(Status(c2)==1);
                f2=find(Status(c2)==2);
                if(f1)             
                    Rate(c2(f1))=Rate(c2(f1))-1;
                end
                if(f2)
                    Rate(c2(f2))=Rate(c2(f2))+0.01;
                end
            case 3
                f1=find(Status(c3)==1);
                f2=find(Status(c3)==2);
                if(f1)             
                    Rate(c3(f1))=Rate(c3(f1))-1;
                end
                if(f2)
                    Rate(c3(f2))=Rate(c3(f2))+0.01;
                end
            case 4
                f1=find(Status(c4)==1);
                f2=find(Status(c4)==2);
                if(f1)             
                    Rate(c4(f1))=Rate(c4(f1))-1;
                end
                if(f2)
                    Rate(c4(f2))=Rate(c4(f2))+0.01;
                end
            case 5
                f1=find(Status(c5)==1);
                f2=find(Status(c5)==2);
                if(f1)             
                    Rate(c5(f1))=Rate(c5(f1))-1;
                end
                if(f2)
                    Rate(c5(f2))=Rate(c5(f2))+0.01;
                end
            case 6
                f1=find(Status(c6)==1);
                f2=find(Status(c6)==2);
                if(f1)             
                    Rate(c6(f1))=Rate(c6(f1))-1;
                end
                if(f2)
                    Rate(c6(f2))=Rate(c6(f2))+0.01;
                end
            case 7
                f1=find(Status(c7)==1);
                f2=find(Status(c7)==2);
                if(f1)             
                    Rate(c7(f1))=Rate(c7(f1))-1;
                end
                if(f2)
                    Rate(c7(f2))=Rate(c7(f2))+0.01;
                end
            case 8
                f1=find(Status(c8)==1);
                f2=find(Status(c8)==2);
                if(f1)             
                    Rate(c8(f1))=Rate(c8(f1))-1;
                end
                if(f2)
                    Rate(c8(f2))=Rate(c8(f2))+0.01;
                end
            case 9
                f1=find(Status(c9)==1);
                f2=find(Status(c9)==2);
                if(f1)             
                    Rate(c9(f1))=Rate(c9(f1))-1;
                end
                if(f2)
                    Rate(c9(f2))=Rate(c9(f2))+0.01;
                end
            case 10
                f1=find(Status(c10)==1);
                f2=find(Status(c10)==2);
                if(f1)             
                    Rate(c10(f1))=Rate(c10(f1))-1;
                end
                if(f2)
                    Rate(c10(f2))=Rate(c10(f2))+0.01;
                end
        end
        % For SIR type dynamics we require the following 2 lines
        Rate(Event)=0;
        Rate(contacts)=Rate(contacts)-tau;

        % For SIS type dynamics we require the following 3 lines
% 	     Status(Event)=1;
%          Rate(Event)=tau*length(find(G(:,Event) & Status==2));
%          Rate(contacts)=Rate(contacts)-tau;

end
Rate=Rate.*sign(Status);
%sign(Status)Ԫ��ֵΪ��������1��Ϊ��������-1��Ϊ0.����0
%a.*b Ϊa��b��ӦԪ�����
%log()��eΪ�� ����step�����Ǹ�ƨ�ð��� step����ѭ��ʱ�䣬maxtime���
step=-log(rand(1,1))/Sum;

% Does a simple check on the value
function []=CheckGreaterOrEqual(Parameter, Value, str)

m=find(Parameter<Value);
if length(m)>0
    error('Parameter %s(%g) (=%g) is less than %g',str,m(1),Parameter(m(1)),Value);
end

function []=CheckGreater(Parameter, Value, str)

m=find(Parameter<=Value);
if length(m)>0
    error('Parameter %s(%g) (=%g) is less than %g',str,m(1),Parameter(m(1)),Value);
end

