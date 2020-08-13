%% �f�[�^�̓ǂݍ���
dt=0.2;                         %sec
dt=dt/60;  %�t���[��������̎��ԁi���j
[numframe,back] = size(data);
numcell = back-2;
Time = data(:,1)*dt;
Intback=data(:,back);

%1��ځA�ŏI����폜
data(:,back)=[];
data(:,1)=[];

%�u�����x�̐��K��
Int=data-Intback;
Int=Int./Int(1,:);
AllInt=mean(Int,2); %�S�זE�̕��ς̌u���P�x

%���όu�����Z�o
G=mean(Int,2);
%% ���X�^�[�̍쐻

%臒l�̐ݒ�i�����l�����ρ{3sd���z�����Ƃ��j
D=diff(Int);    %�����l��D�Ɋi�[
Dshre=mean(D)+(3*std(D));   %臒l�ݒ�
Int(numframe,:)=[]; %�Ō�̃t���[�����폜

%���X�^�[�̍쐻
ras = D>Dshre;  %臒l���z�����Ƃ��l��t�^�i���_�l�j
ras = double(ras);  %���l�ɕύX
NA = ras;   %NA (neural activity)�͔����l���z���Ă���������ׂđΏۂƂ���

%1�������ӏ����폜
for i=1:numframe-2
    for j=1:numcell
        if (ras(numframe-i,j)==1)&&(ras(numframe-(i+1),j)==1)
            ras(numframe-i,j)=0;
        end
    end
end

%���ΐ��̉��
fire=sum(ras);  %���ΐ�
firingrate=fire/(dt*numframe);  %1���ԓ�����̔��Εp�x
[firerank,firerankcell]=sort(fire,'descend');   %���ΐ����~���Ƀ\�[�g
fireave=sum(fire)/(numcell*Time(numframe));     %1�זE��1��������̔��ΐ�

%% �l�b�g���[�N�̊����񐔂𒲂ׂ�
bin=round(1/(60*dt));          %1 sec�Ƃ���̂ɕK�v�ȃr����
rasframe=size(ras);     %numframe-1�̒l
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

%% Generalized transfer entropy GTE(J��I)��GTE(i,j)�ɓ����
% �e�x�N�g�������� 
k=2;    %���t���[���O�܂ł����̂ڂ邩
S=1; %S=1�͓��t���[�����l������AS=0�͍l�����Ȃ�

cell1=NA+1;    %it���i�[
cell1back=zeros(numframe-1,k,numcell);    %it-1(k)���i�[
cell2back=zeros(numframe-1,k,numcell);    %jt-1+S(k)���i�[
for t=1+k:numframe-1
    for i=1:k
        cell1back(t,i,:)=NA(t-i,:);    %t����i�t���[���O��ras�̃f�[�^������
        cell2back(t,i,:)=NA(t-i+S,:);
    end
end

%cell1back(cell2back)���Ƃ肤��l��(0,0),(0,1),(1,0),(1,1)�B
%2�i�@�ōl���āA���ꂼ���0,1,2,3�Ƃ����l�ɕϊ�����recell1back�Ƃ����z��Ɋi�[
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

recell1back=recell1back+1;  %�v�Z���₷���悤��1�����Z����
recell2back=recell2back+1;


% �m�������߂�
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
                    if P1(i,j,k)~=0  %0*log(0)�y��0*log(1/0)=0�Ȃ̂ŁAP1=0�����O����
                        %c2��c1��generalized transfer entropy
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

GTEsort=reshape(GTE,1,numcell*numcell);     %GTE��1�����z��Ɂi�\�[�g���邽�߁j
GTEsort=sort(GTEsort,'descend','MissingPlacement','last');      %�~���Ƀ\�[�g

shre=round((numcell*numcell-numcell)*0.1);      %top10%�ȏ�̒l�����Ƃ������Ƃ݂Ȃ�

FCGTE=GTE>=GTEsort(shre);       %top10%�ȏ�̒l��FCGTE�Ɋi�[�ilogical�z��j
[GTEin,GTEout] = find(FCGTE==1);    %���������זE�̒��o

numGTEin=zeros(1,numcell);    %i�Ԗڂ̍זE�̓��͐�
numGTEout=zeros(1,numcell);   %i�Ԗڂ̍זE�̏o�͐�
for i=1:numcell
    numGTEin(i)=nnz(GTEin(:)==i);
    numGTEout(i)=nnz(GTEout(:)==i);
end

meanGTE=mean(numGTEin); %1�זE�����������̕���(top10%�Ȃ̂ōזE��*0.1�Ɠ����悤�Ȑ��ɂȂ�j
norin=numGTEin()/meanGTE;    %���ϒl�Ő��K��
norout=numGTEout()/meanGTE;

%���ψȏ�̌����������זE�̒��o(high connectivity cell)
highinROI=find(norin()>=1);  %ROI�ԍ�
highoutROI=find(norout()>1);
highinscore=norin(highinROI);
highoutscore=norout(highoutROI);

%���������זE�̒��o (highconnectivity neuron)
highROI=[highinROI highoutROI];
highROI=sort(highROI);
for i=2:length(highROI)
    if highROI(i)==highROI(i-1)
        highROI(i)=nan;
    end
end
highROI=rmmissing(highROI);

%���̑��j���[�����̒��o (other neuron)
otherROI=(1:numcell);
for i=1:length(highROI)
    if nnz(otherROI()==highROI(i))~=0
        otherROI(otherROI()==highROI(i))=nan;
    end
end
otherROI=rmmissing(otherROI);
otherin=norin(otherROI);
otherout=norout(otherROI);

%HC�j���[�����������2�ɂ킯��B1�ƁAnorin,norout�̍ő�l�Ƃ̕��ϒl��臒l�Ƃ���B
topinROI=find(norin>mean([max(norin),1]));        %Top input neuron��ROI�ԍ�
topoutROI=find(norout>mean([max(norout),1]));     %Top output neuron��ROI�ԍ�
topROI=[topinROI,topoutROI];
topROI=sort(topROI);
for i=2:length(topROI)
    if topROI(i)==topROI(i-1)
        topROI(i)=nan;
    end
end
topROI=rmmissing(topROI);
%% �����i�[
connection=[GTEin GTEout];  %1�s��out��2�s��in
cennection={'in','out';GTEin,GTEout};
norGTEscore={'in','out';norin.',norout.'};
numGTE={'in','out';numGTEin.',numGTEout.'};
GTEdata={'�����̗L��(i,j)��j��i�ւ̌���',...
    '����������ROI',...
    '������',...
    '���K������������',...
    'HC neuron',...
    'Highinneuron',...
    'Highoutneuron',...
    'Other neuron',...
    '1�זE��1��������̔��ΐ�',...
    '�l�b�g���[�N�̊����񐔁i/min)',...
    'Topinputneuron��ROI',...
    'Topoutputneuron��ROI',...
    'Topneuron��ROI';...
    FCGTE,connection,numGTE,norGTEscore,...
    highROI.',highinROI.',highoutROI.',otherROI.',...
    fireave,numevent,topinROI.',topoutROI.',topROI.'};
filename='data.mat';
savefile=strcat(savedir,filename);
save(savefile,'GTEdata');

%% �O���t���o��
showfigure

%% �G�N�Z���ւ̏o��
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

HCratio={'HC�̐�','Other�̐�';...
    length(highROI),length(otherROI)};
writecell(HCratio,savefile,'Range','F1');
writematrix(GTEdata{1,9},savefile,'Range','H1');
writematrix(GTEdata{1,10},savefile,'Range','I1');
writematrix(GTEdata{2,9},savefile,'Range','H2');
writematrix(GTEdata{2,10},savefile,'Range','I2');
writecell(inf,savefile,'Range','J1');
