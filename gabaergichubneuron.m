%% GABAergic Hub Neuronより
FCGABA=zeros(numcell,numcell);  %機能的結合があるとき値を入れる
for i=1:numcell
    for j=i+1:numcell
        distribution=[];    %細胞iに対する細胞jの発火分布
        for k=5:numframe-5
            if ras(k,i)==1
                for l=-4:4
                    if (ras(k+l,j))==1
                        lag=(l);
                        distribution=[distribution lag];
                    end
                end
            end
        end
        if (ttest(distribution)==1)&&(kstest(((distribution-mean(distribution))/std(distribution))==1))
            if mean(distribution)>0 
                FCGABA(i,j)=1;  %i→jで発火
            else
                FCGABA(j,i)=1;  %j→iで発火
            end
        end        
    end
end
[GABAout,GABAin]=find(FCGABA==1);

numGABAin=zeros(1,numcell);    %i番目の細胞の入力数
numGABAout=zeros(1,numcell);   %i番目の細胞の出力数
for i=1:numcell
    numGABAin(i)=nnz(GABAin(:)==i);
    numGABAout(i)=nnz(GABAout(:)==i);
end
