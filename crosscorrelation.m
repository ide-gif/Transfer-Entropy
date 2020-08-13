%% ���ݑ��։��(���n������̂܂܉�́j
XC=zeros(numcell,numcell,2);
maxlag=4;   %���ւ��Ƃ�͈́i�t���[���j
for i=1:numcell
    for j=1:numcell
        XC(i,j,1)=max(xcorr(Int(:,i),Int(:,j),maxlag,'coeff'));
        XC(i,j,2)=finddelay(Int(:,i),Int(:,j),maxlag);    %�l�����̂Ƃ�i�̂ق����������΂��Ă���
    end
end
for i=1:numcell
    XC(i,i,1)=NaN;
end

%fisher��Z�ϊ�
zXC=zeros(numcell,numcell);
for i=1:numcell
    for j=1:numcell
        zXC(i,j)=(1/2)*(log((1+XC(i,j,1))/(1-XC(i,j,1))));
    end
end

%zXC�̕���+1.96SD���傫�����̂������Ƃ݂Ȃ�
FCXC=zXC>(nanmean(zXC,'all')+1.96*nanstd(zXC,0,'all'));
[XCout,XCin]=find(FCXC==1);

numXCin=zeros(1,numcell);    %i�Ԗڂ̍זE�̓��͐�
numXCout=zeros(1,numcell);   %i�Ԗڂ̍זE�̏o�͐�
for i=1:numcell
    numXCin(i)=nnz(XCin(:)==i);
    numXcout(i)=nnz(XCout(:)==i);
end

