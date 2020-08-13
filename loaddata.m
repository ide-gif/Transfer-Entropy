%% データの読み込み
dt=0.2;                         %sec
dt=dt/60;  %フレーム当たりの時間（分）
[numframe,back] = size(data);
numcell = back-2;
Time = data(:,1)*dt;
Intback=data(:,back);

%1列目、最終列を削除
data(:,back)=[];
data(:,1)=[];

%蛍光強度の正規化
Int=data-Intback;
Int=Int./Int(1,:);
AllInt=mean(Int,2); %全細胞の平均の蛍光輝度

%平均蛍光を算出
G=mean(Int,2);
%% ラスターの作製

%閾値の設定（微分値が平均＋3sdを越えたとき）
D=diff(Int);    %微分値をDに格納
Dshre=mean(D)+(3*std(D));   %閾値設定
Int(numframe,:)=[]; %最後のフレームを削除

%ラスターの作製
ras = D>Dshre;  %閾値を越えたとき値を付与（理論値）
ras = double(ras);  %数値に変更
NA = ras;   %NA (neural activity)は微分値が越えている個所をすべて対象とする

%1が続く箇所を削除
for i=1:numframe-2
    for j=1:numcell
        if (ras(numframe-i,j)==1)&&(ras(numframe-(i+1),j)==1)
            ras(numframe-i,j)=0;
        end
    end
end

%発火数の解析
fire=sum(ras);  %発火数
firingrate=fire/(dt*numframe);  %1分間当たりの発火頻度
[firerank,firerankcell]=sort(fire,'descend');   %発火数を降順にソート
fireave=sum(fire)/(numcell*Time(numframe));     %1細胞の1分当たりの発火数

%% ネットワークの活動回数を調べる
bin=round(1/(60*dt));          %1 secとするのに必要なビン数
rasframe=size(ras);     %numframe-1の値
rasframe=rasframe(1);
newras=zeros(round(rasframe/bin),numcell);
for i=1:(round(rasframe/bin))
    for j=1:bin
        newras(i,:)=newras(i,:)+ras(5*i-(j-1),:);
    end
end
newras(newras>1)=1;
newras=sum(newras.');
numevent=nnz(newras/numcell>0.7)/(dt*numframe);

%% Generalized transfer entropy GTE(J→I)をGTE(i,j)に入れる
% 各ベクトルをつくる 
k=2;    %何フレーム前までさかのぼるか
S=1; %S=1は同フレームを考慮する、S=0は考慮しない

cell1=NA+1;    %itを格納
cell1back=zeros(numframe-1,k,numcell);    %it-1(k)を格納
cell2back=zeros(numframe-1,k,numcell);    %jt-1+S(k)を格納
for t=1+k:numframe-1
    for i=1:k
        cell1back(t,i,:)=NA(t-i,:);    %tからiフレーム前のrasのデータを入れる
        cell2back(t,i,:)=NA(t-i+S,:);
    end
end

%cell1back(cell2back)がとりうる値は(0,0),(0,1),(1,0),(1,1)。
%2進法で考えて、それぞれを0,1,2,3という値に変換してrecell1backという配列に格納
recell1back=zeros(numframe-1,numcell);
recell2back=zeros(numframe-1,numcell);

recell1back(:,:)=cell1back(:,1,:);
recell2back(:,:)=cell2back(:,1,:);
for i=1:numcell
    for j=1:k-1
        recell1back(:,i)=recell1back(:,i)+(cell1back(:,j+1,i)*(2^j));
        recell2back(:,i)=recell2back(:,i)+(cell2back(:,j+1,i)*(2^j));
    end
end

recell1back=recell1back+1;  %計算しやすいように1を加算する
recell2back=recell2back+1;


% 確率を求める
dim=2^k;
P1 = zeros(2,dim,dim);  %it,it-1(k),jt-1+S(k)
P2 = zeros(dim,dim);    %it-1(k),jt-1+S(k)
P3 = zeros(2,dim);      %it(k),it-1(k)
P4 = zeros(1,dim);        %it-1(k)

GTE=zeros(numcell,numcell); 
for c1=1:numcell
    for c2=1:numcell
        for i=1:2
            for j=1:dim
                for k=1:dim
                    P1(i,j,k)=nnz(cell1(:,c1)==i&recell1back(:,c1)==j&recell2back(:,c2)==k);
                    P2(j,k)=nnz(recell1back(:,c1)==j&recell2back(:,c2)==k);
                    P3(i,j)=nnz(cell1(:,c1)==i&recell1back(:,c1)==j);
                    P4(j)=nnz(recell1back(:,c1)==j);
                    if P1(i,j,k)~=0  %0*log(0)及び0*log(1/0)=0なので、P1=0を除外する
                        %c2→c1のgeneralized transfer entropy
                        GTE(c1,c2)=GTE(c1,c2)+(P1(i,j,k)/numframe)*log((P1(i,j,k)*P4(j))/(P3(i,j)*P2(j,k)));
                    end
                end
            end
        end
    end
end

for i=1:numcell
    GTE(i,i)=NaN;
end

GTEsort=reshape(GTE,1,numcell*numcell);     %GTEを1次元配列に（ソートするため）
GTEsort=sort(GTEsort,'descend','MissingPlacement','last');      %降順にソート

shre=round((numcell*numcell-numcell)*0.1);      %top10%以上の値をもつとき結合とみなす

FCGTE=GTE>=GTEsort(shre);       %top10%以上の値をFCGTEに格納（logical配列）
[GTEin,GTEout] = find(FCGTE==1);    %結合を持つ細胞の抽出

numGTEin=zeros(1,numcell);    %i番目の細胞の入力数
numGTEout=zeros(1,numcell);   %i番目の細胞の出力数
for i=1:numcell
    numGTEin(i)=nnz(GTEin(:)==i);
    numGTEout(i)=nnz(GTEout(:)==i);
end

meanGTE=mean(numGTEin); %1細胞が持つ結合数の平均(top10%なので細胞数*0.1と同じような数になる）
norin=numGTEin()/meanGTE;    %平均値で正規化
norout=numGTEout()/meanGTE;

%平均以上の結合数を持つ細胞の抽出(high connectivity cell)
highinROI=find(norin()>=1);  %ROI番号
highoutROI=find(norout()>1);
highinscore=norin(highinROI);
highoutscore=norout(highoutROI);

%高結合性細胞の抽出 (highconnectivity neuron)
highROI=[highinROI highoutROI];
highROI=sort(highROI);
for i=2:length(highROI)
    if highROI(i)==highROI(i-1)
        highROI(i)=nan;
    end
end
highROI=rmmissing(highROI);

%その他ニューロンの抽出 (other neuron)
otherROI=(1:numcell);
for i=1:length(highROI)
    if nnz(otherROI()==highROI(i))~=0
        otherROI(otherROI()==highROI(i))=nan;
    end
end
otherROI=rmmissing(otherROI);
otherin=norin(otherROI);
otherout=norout(otherROI);

%HCニューロンをさらに2つにわける。1と、norin,noroutの最大値との平均値を閾値とする。
topinROI=find(norin>mean([max(norin),1]));        %Top input neuronのROI番号
topoutROI=find(norout>mean([max(norout),1]));     %Top output neuronのROI番号
topROI=[topinROI,topoutROI];
topROI=sort(topROI);
for i=2:length(topROI)
    if topROI(i)==topROI(i-1)
        topROI(i)=nan;
    end
end
topROI=rmmissing(topROI);
%% 情報を格納
connection=[GTEin GTEout];  %1行目out→2行目in
cennection={'in','out';GTEin,GTEout};
norGTEscore={'in','out';norin.',norout.'};
numGTE={'in','out';numGTEin.',numGTEout.'};
GTEdata={'結合の有無(i,j)でj→iへの結合',...
    '結合を持つROI',...
    '結合数',...
    '正規化した結合数',...
    'HC neuron',...
    'Highinneuron',...
    'Highoutneuron',...
    'Other neuron',...
    '1細胞の1分当たりの発火数',...
    'ネットワークの活動回数（/min)',...
    'TopinputneuronのROI',...
    'TopoutputneuronのROI',...
    'TopneuronのROI';...
    FCGTE,connection,numGTE,norGTEscore,...
    highROI.',highinROI.',highoutROI.',otherROI.',...
    fireave,numevent,topinROI.',topoutROI.',topROI.'};
filename='data.mat';
savefile=strcat(savedir,filename);
save(savefile,'GTEdata');

%% グラフを出力
showfigure

%% エクセルへの出力
dataname='results.xlsx';
savefile=strcat(savedir,dataname);

type=repmat("HC",length(highROI),1);
type=[type;repmat("Other",length(otherROI),1)];
ROI=[highROI.';otherROI.'];
ROItype=[ROI,type];
subject={'ROI','type'};
writecell(subject,savefile,'Range','A1');
writematrix(ROItype,savefile,'Range','A2');
writecell(inf,savefile,'Range','C1');

HCratio={'HCの数','Otherの数';...
    length(highROI),length(otherROI)};
writecell(HCratio,savefile,'Range','F1');
writematrix(GTEdata{1,9},savefile,'Range','H1');
writematrix(GTEdata{1,10},savefile,'Range','I1');
writematrix(GTEdata{2,9},savefile,'Range','H2');
writematrix(GTEdata{2,10},savefile,'Range','I2');
writecell(inf,savefile,'Range','J1');
