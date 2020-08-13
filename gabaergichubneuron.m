%% GABAergic Hub Neuron���
FCGABA=zeros(numcell,numcell);  %�@�\�I����������Ƃ��l������
for i=1:numcell
    for j=i+1:numcell
        distribution=[];    %�זEi�ɑ΂���זEj�̔��Ε��z
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
                FCGABA(i,j)=1;  %i��j�Ŕ���
            else
                FCGABA(j,i)=1;  %j��i�Ŕ���
            end
        end        
    end
end
[GABAout,GABAin]=find(FCGABA==1);

numGABAin=zeros(1,numcell);    %i�Ԗڂ̍זE�̓��͐�
numGABAout=zeros(1,numcell);   %i�Ԗڂ̍זE�̏o�͐�
for i=1:numcell
    numGABAin(i)=nnz(GABAin(:)==i);
    numGABAout(i)=nnz(GABAout(:)==i);
end
