%% 時系列グラフを表示するプログラム
i=20;
plotT=Time;
plotT(numframe)=[];
figure
for n=i:i+10
    plot(plotT,Int(:,i));
end
