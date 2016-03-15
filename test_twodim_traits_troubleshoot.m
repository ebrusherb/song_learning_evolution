alpha=.5;
x=.25*ones(4,1);
y=[.2 .3];
z=[.3 .2];
%% recursion equations

F=zeros(4,4);
prefs=ones(4,4);
prefs([1 3 5 8 10 11 14 16])=(1+alpha);
F=prefs.*repmat(reshape(x,1,[]),4,1);
F=F./repmat(sum(F,2),1,4);
F=repmat([col(y); col(z)],1,4).*F;

newx=col(sum(F,1));
newy=[.5*(sum(F(:,1))+sum(F(:,2))); .5*(sum(F(:,3))+sum(F(:,3)))];
newz=[.5*(sum(F(:,1))+sum(F(:,3))); .5*(sum(F(:,2))+sum(F(:,4)))];

%%
alpha=2;
% x=[.49;.01;.01;.49];
% y=[.5 0];
% z=[0 .5];

x=[.49;.01;.01;.49];
y=[.49 .01];
z=[.01 .49];

T=1000;
popdyns=zeros(8,T+1);
popdyns(1:4,1)=col(x);
popdyns(5:6,1)=col(y);
popdyns(7:8,1)=col(z);

for t=1:T
    x=popdyns(1:4,t);
    y=popdyns(5:6,t);
    z=popdyns(7:8,t);
    F=zeros(4,4);
    prefs=ones(4,4);
    prefs([1 3 5 8 10 11 14 16])=(1+alpha);
    F=prefs.*repmat(reshape(x,1,[]),4,1);
    F=F./repmat(sum(F,2),1,4);
    F=repmat([col(y); col(z)],1,4).*F;

    newx=col(sum(F,1));
    newy=[.5*(sum(F(:,1))+sum(F(:,2))); .5*(sum(F(:,3))+sum(F(:,4)))];
    newz=[.5*(sum(F(:,1))+sum(F(:,3))); .5*(sum(F(:,2))+sum(F(:,4)))];
    popdyns(1:4,t+1)=col(newx);
    popdyns(5:6,t+1)=col(newy);
    popdyns(7:8,t+1)=col(newz);
end

subplot(2,1,1)
plot(1:(T+1),popdyns(1:4,:))
set(gca,'ylim',[0 1])
legend({'x11','x12','x21','x22'})
subplot(2,1,2)
plot(1:(T+1),popdyns(5:8,:))
legend({'y1','y2','z1','z2'})


%%
alpha1=2.01;
alpha2=2;
x=rand(8,1);
x=x/sum(x);
x=[.5 0 0 0 0 0 0 .5];
r=.5;
y=r*[.9 .1];
z=(1-r)*[0.1 .9];

T=100;
popdyns=zeros(12,T+1);
popdyns(1:8,1)=col(x);
popdyns(9:10,1)=col(y);
popdyns(11:12,1)=col(z);


prefs=ones(4,4);
prefs([1 5 10 14])=1+alpha1;
prefs([3 8 11 16])=1+alpha2;
prefs=repmat(prefs,1,2);

T1offspring=ones(4,8);
T1offspring(3:4,1:4)=1/2;
T1offspring(1:2,5:8)=1/2;
T1offspring(3:4,5:8)=0;
T2offspring=ones(4,8);
T2offspring(1:2,1:4)=0;
T2offspring(3:4,1:4)=1/2;
T2offspring(1:2,5:8)=1/2;

check=zeros(T+1,1);
for t=1:T
    x=popdyns(1:8,t);
    x=round(x*1E8)/1E8;
    x=x/sum(x);
%     sum(x)
    y=popdyns(9:10,t);
    z=popdyns(11:12,t);
    pickmales=prefs.*repmat(reshape(x,1,[]),4,1);
    F1=pickmales./repmat(sum(pickmales,2),1,8);
    F2=pickmales./repmat([1+alpha1*sum(x([1 2 5 6])); 1+alpha1*sum(x([3 4 7 8])); 1+alpha2*sum(x([1 3 5 7])); 1+alpha2*sum(x([2 4 6 8]))],1,8);
    F=repmat([col(y); col(z)],1,8).*F1;
%     sum(sum(F))
    check(t)=sum(sum(F1))-sum(sum(F2));
    
    T1now=T1offspring.*F;
    T2now=T2offspring.*F;
    newx=zeros(1,8);
    newx(1:4)=sum(T1now(:,1:4),1)+sum(T1now(:,5:8),1);
    newx(5:8)=sum(T2now(:,1:4),1)+sum(T2now(:,5:8),1);
    
    newy=[sum(sum(T1now(:,[1 2 5 6]))) sum(sum(T1now(:,[3 4 7 8])))];
    newz=[sum(sum(T2now(:,[1 3 5 7]))) sum(sum(T2now(:,[2 4 6 8])))];

    popdyns(1:8,t+1)=col(newx);
    popdyns(9:10,t+1)=col(newy);
    popdyns(11:12,t+1)=col(newz);
end

eightcols=cbrewer('qual','Set1',8);
subplot(2,1,1)
cla
hold on
l=zeros(1,8);
for i=1:8
    l(i)=plot(1:(T+1),popdyns(i,:),'Color',eightcols(i,:));
end
set(gca,'ylim',[0 1])
legend(l,{'T1s11','T1s12','T1s21','T1s22','T2s11','T2s12','T2s21','T2s22'})
hold off
subplot(2,1,2)
plot(1:(T+1),popdyns(9:12,:))
legend({'y1','y2','z1','z2'})

%%
alpha1=2.01;
alpha2=2;
x=rand(8,1);
x=x/sum(x);
x=[.5 0 0 0 0 0 0 .5];
r=.5;
y=r*[.9 .1];
z=(1-r)*[0.1 .9];

T=100;
popdyns=zeros(12,T+1);
popdyns(1:8,1)=col(x);
popdyns(9:10,1)=col(y);
popdyns(11:12,1)=col(z);


prefs=ones(4,4);
prefs([1 5 10 14])=1+alpha1;
prefs([3 8 11 16])=1+alpha2;
prefs=repmat(prefs,1,2);

T1offspring=ones(4,8);
T1offspring(3:4,1:4)=1/2;
T1offspring(1:2,5:8)=1/2;
T1offspring(3:4,5:8)=0;
T2offspring=ones(4,8);
T2offspring(1:2,1:4)=0;
T2offspring(3:4,1:4)=1/2;
T2offspring(1:2,5:8)=1/2;

check=zeros(T+1,1);
for t=1:T
    x=popdyns(1:8,t);
    x=round(x*1E8)/1E8;
    x=x/sum(x);
%     sum(x)
    y=popdyns(9:10,t);
    z=popdyns(11:12,t);
    pickmales=prefs.*repmat(reshape(x,1,[]),4,1);
    F1=pickmales./repmat(sum(pickmales,2),1,8);
    F2=pickmales./repmat([1+alpha1*sum(x([1 2 5 6])); 1+alpha1*sum(x([3 4 7 8])); 1+alpha2*sum(x([1 3 5 7])); 1+alpha2*sum(x([2 4 6 8]))],1,8);
    F=repmat([col(y); col(z)],1,8).*F1;
%     sum(sum(F))
    check(t)=sum(sum(F1))-sum(sum(F2));
    
    T1now=T1offspring.*F;
    T2now=T2offspring.*F;
    newx=zeros(1,8);
    newx(1:4)=sum(T1now(:,1:4),1)+sum(T1now(:,5:8),1);
    newx(5:8)=sum(T2now(:,1:4),1)+sum(T2now(:,5:8),1);
    
    newy=[sum(sum(T1now(:,[1 2 5 6]))) sum(sum(T1now(:,[3 4 7 8])))];
    newz=[sum(sum(T2now(:,[1 3 5 7]))) sum(sum(T2now(:,[2 4 6 8])))];

    popdyns(1:8,t+1)=col(newx);
    popdyns(9:10,t+1)=col(newy);
    popdyns(11:12,t+1)=col(newz);
end

eightcols=cbrewer('qual','Set1',8);
subplot(2,1,1)
cla
hold on
l=zeros(1,8);
for i=1:8
    l(i)=plot(1:(T+1),popdyns(i,:),'Color',eightcols(i,:));
end
set(gca,'ylim',[0 1])
legend(l,{'T1s11','T1s12','T1s21','T1s22','T2s11','T2s12','T2s21','T2s22'})
hold off
subplot(2,1,2)
plot(1:(T+1),popdyns(9:12,:))
legend({'y1','y2','z1','z2'})



