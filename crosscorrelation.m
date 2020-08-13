%% 相互相関解析(時系列をそのまま解析）
XC=zeros(numcell,numcell,2);
maxlag=4;   %相関をとる範囲（フレーム）
for i=1:numcell
    for j=1:numcell
        XC(i,j,1)=max(xcorr(Int(:,i),Int(:,j),maxlag,'coeff'));
        XC(i,j,2)=finddelay(Int(:,i),Int(:,j),maxlag);    %値が正のときiのほうが早く発火している
    end
end
for i=1:numcell
    XC(i,i,1)=NaN;
end

%fisherのZ変換
zXC=zeros(numcell,numcell);
for i=1:numcell
    for j=1:numcell
        zXC(i,j)=(1/2)*(log((1+XC(i,j,1))/(1-XC(i,j,1))));
    end
end

%zXCの平均+1.96SDより大きいものを結合とみなす
FCXC=zXC>(nanmean(zXC,'all')+1.96*nanstd(zXC,0,'all'));
[XCout,XCin]=find(FCXC==1);

numXCin=zeros(1,numcell);    %i番目の細胞の入力数
numXCout=zeros(1,numcell);   %i番目の細胞の出力数
for i=1:numcell
    numXCin(i)=nnz(XCin(:)==i);
    numXcout(i)=nnz(XCout(:)==i);
end

