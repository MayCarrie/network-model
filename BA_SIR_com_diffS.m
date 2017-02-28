function [t,S,I] = BA_SIR_com_diffS(N,n,tau,gamma,MaxTime,Type)

%实现SIS在BA无标度网络上传播，利用文件FreeScale中生成adjacent_matrix，把其赋值给G即可。
%！！文件中似乎注释掉了SIR传播模型，这个可以研究一下。
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
%用kmeans划分10个社区
[idx,c] = kmeans(G,com_num);
%X  Y分别为列向量和行向量 
Status=1+0*X; Status(1)=2;
Rate=0*X; Rate(1)=gamma; Rate(find(G(:,1)))=tau; 
%Rate为列向量 Rate(1)为第一个元素，令其为gamma，
%G(:,1)表示取第一列的所有元素，那么find(G(:,1))就是找到第一列中不为零的索引值
t=0; i=1; S=N-1; I=1;R=0;
%区别主动和被动的S状态，主动的Rate设为10
x=0;
while(x<N/10)
tempt=randint(1,1,[1,N]);
Rate(tempt)=10;
x=x+1;
end

% The main iteration
subplot(2,1,1);
%find(G) 返回G中的非零元素，行、列、值
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
    %Status是个矩阵，length(find(Status==1))为矩阵中值为1的数量
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
%X Y 分别为N个元素的列向量和行向量
X=rand(N,1); Y=rand(N,1);
switch Type
    
%random网络整体思想：随机产生i和j，使得G(i,j)随机为1，以节点的度*节点数来控制循环。
    case {'Random','random'}
        contacts=0;
        while(contacts<n*N)
            i=ceil(rand(1,1)*N); j=ceil(rand(1,1)*N);
            %rand(1,1)随机产生1*1阶0~1之间的数，ceil(rand(1,1)*N)就是1~100之间的整数
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
%以下便是生成BA无标度网络的代码，本大侠的！    
    case {'BA','ScaleFree','Scalefree','scalefree','BAscalefree','ba'}
        m0= 3; m= 3;%初始化网络数据
        adjacent_matrix = sparse( m0, m0);% 初始化邻接矩阵
        for i = 1: m0            
            for j = 1:m0
                if j ~= i      %去除每个点自身形成的环
                    adjacent_matrix(i,j) = 1;%建立初始邻接矩阵，3点同均同其他的点相连
                end
            end
        end
        adjacent_matrix =sparse(adjacent_matrix);%邻接矩阵稀疏化
        node_degree = zeros(1,m0+1);                %初始化点的度
        node_degree(2: m0+1) = sum(adjacent_matrix);%对度维数进行扩展
        for iter= 4:N
            iter;                                %加点
            total_degree = 2*m*(iter- 4)+6;%计算网络中此点的度之和
            cum_degree = cumsum(node_degree);%求出网络中点的度矩阵
            choose= zeros(1,m);%初始化选择矩阵
    % 选出第一个和新点相连接的顶点
            r1= rand(1)*total_degree;%算出与旧点相连的概率
            for i= 1:iter-1
                if (r1>=cum_degree(i))&( r1<cum_degree(i+1))%选取度大的点
                    choose(1) = i;
                    break
                end
            end
            % 选出第二个和新点相连接的顶点
            r2= rand(1)*total_degree;
            for i= 1:iter-1
                if (r2>=cum_degree(i))&(r2<cum_degree(i+1))
                    choose(2) = i;
                    break
                end
            end
            while choose(2) == choose(1)%第一个点和第二个点相同的话，重新择优
                r2= rand(1)*total_degree;
                for i= 1:iter-1
                    if (r2>=cum_degree(i))&(r2<cum_degree(i+1))
                        choose(2) = i;
                        break
                    end
                end
            end
            % 选出第三个和新点相连接的顶点
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
            %新点加入网络后, 对邻接矩阵进行更新
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

%使用Kmeans方法划分社区，c~c9分别存放了所属不同社区的结点，km则是所属社区的标签
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
%虽然是随机的，但是Event会大概率的等于突然变大的元素索引，且关键是Event保证不为Rate中0元素的索引。该语句为随机找到一个与节点相连的节点索引
Status(Event)=mod(Status(Event)+1,3);
%如果Status(Event)本身为S（=1）时，那么赋值2，跳转case2；如果本身为I（=2），则赋值0，跳转case 0，又被赋值1，即变为S
contacts=find(G(:,Event) & Status==1);
%如果I增多，接触率则会变大
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
 %  以上switch的代碼：如果社區 c1中有感染的結點，那麼該社區中的結點的感染概率都增加1     
        Rate(Event)=gamma;
        Rate(contacts)=Rate(contacts)+tau;
        %上两行代码实现设定I的恢复率和S的感染率
        
        %暂时先不要 ，大概没什么影响G(Event,:)=0;
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
%sign(Status)元素值为正：返回1；为负，返回-1；为0.返回0
%a.*b 为a和b对应元素相乘
%log()以e为底 所以step到底是干屁用啊？ step控制循环时间，maxtime相关
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

