i=3;
j=20;

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
ttest(distribution)
kstest(((distribution-mean(distribution))/std(distribution)))
figure
histogram(distribution);