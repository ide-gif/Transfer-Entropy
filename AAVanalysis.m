%% ����
%highindata�ɂ�sourceROI�ō����͐��̍זE��target�łǂ̂悤�Ȑ��������o��
%highinnROI1�ɂ�source��ROI�ԍ����Atranshighin�ɂ�target�ƑΉ�����ROI�ԍ��������Ă���
%8,9,10,21,22�s�ڂ�ύX����
%loaddata��1�s�ڂ�target��DIV�ɕύX����
clearvars
%%
source='DIV11';
target='DIV9';
datainf='191228-3-area1';

filename='C:\analysisresults\AAV\';   
dir=strcat(filename,datainf,'\');   %�f�[�^���i�[���Ă���f�B���N�g��
filename=strcat(filename,datainf,'\',datainf,'.xlsx');  %ROI�ԍ������Ă�t�@�C��
datamat='\data.mat';
sourcedir=strcat(dir,source,datamat);
targetdir=strcat(dir,target,datamat);
%ROI�̐��ڃf�[�^��ǂݍ���
ROItrans=readtable(filename,'PreserveVariableName',true);
%source��ROI�ԍ���target��ROI�ԍ��ɕϊ�
sourceROI=ROItrans.DIV11;
targetROI=ROItrans.DIV9;

%source��GTE�f�[�^��g�ݍ���
sourcefile=load(sourcedir);
sourcedata=sourcefile.GTEdata;

%loaddata�̃f�B���N�g����ύX����
targetfile=load(targetdir);
targetdata=targetfile.GTEdata;

%% 1�߂ł�highconnectivity��2�߂ł�highconnectivity�����ׂ�
transhighin=zeros(1,length(sourcedata{2,6}));    %�ϊ����source�̍��������j���[����
transhighout=zeros(1,length(sourcedata{2,7}));
transother=zeros(1,length(sourcedata{2,8}));
for i=1:length(sourceROI)
    transhighin(sourcedata{2,6}==sourceROI(i))=targetROI(i);
    transhighout(sourcedata{2,7}==sourceROI(i))=targetROI(i);
    transother(sourcedata{2,8}==sourceROI(i))=targetROI(i);
end

%% �����ō���������ROI�ԍ��𒊏o����
%source�cDIV13,target�cDIV11
%�g������otherROI1,otherROI,numcell,length(FCGTE),targetROI,sourceROI

j=1;
bothhighROI=[;];
for i=1:length(targetdata{2,5})
    if nnz(sourcedata{2,5}()==sourceROI(targetROI==targetdata{2,5}(i)))~=0
        bothhighROI(1,j)=sourceROI(targetROI==targetdata{2,5}(i));
        bothhighROI(2,j)=targetdata{2,5}(i);
        j=j+1;
    end
end
bothhighROI={strcat(source,'highROI'),strcat(target,'highROI'),'�ǂ����ł�����������ROI�ԍ�';...
    bothhighROI(1,:).',bothhighROI(2,:).',0};


%% �v�����]�[��
%�~�O���t���쐬���邽�߂�table
valueset = 1:3;
type={'Undetected' 'Other' 'Highconnectivity'};

%source�ɂ�����highin��target�ł̐���
sourcehighin=zeros(1,length(transhighin));
for i=1:length(transhighin)
    if nnz(targetdata{2,6}(:)==transhighin(i))==1
        sourcehighin(i)=3;  %��������
    elseif nnz(targetdata{2,7}(:)==transhighin(i))==1
        sourcehighin(i)=3;
    elseif nnz(targetdata{2,8}(:)==transhighin(i))==1
        sourcehighin(i)=2;  %���̑�
    else
        sourcehighin(i)=1;  %�����o
    end
end
highindata=categorical(sourcehighin(:),valueset,type,'Ordinal',true);
figure
pie(highindata);
title('Highin')
%source��highin�̂����A������target��high,other,undetected���Z�o����
highinratio=[(nnz(sourcehighin()==3)),...
    (nnz(sourcehighin()==2)),...
    (nnz(sourcehighin()==1))];

%source�ɂ�����highout��target�ł̐���
sourcehighout=zeros(1,length(transhighout));
for i=1:length(transhighout)
    if nnz(targetdata{2,6}(:)==transhighout(i))==1
        sourcehighout(i)=3;
    elseif nnz(targetdata{2,7}(:)==transhighout(i))==1
        sourcehighout(i)=3;
    elseif nnz(targetdata{2,8}(:)==transhighout(i))==1
        sourcehighout(i)=2;
    else
        sourcehighout(i)=1;
    end
end
highoutdata=categorical(sourcehighout(:),valueset,type,'Ordinal',true);
figure
pie(highoutdata);
title('Highout')

%source��highout�̂����A������target��high,other,undetected���Z�o����
highoutratio=[(nnz(sourcehighout()==3)),...
    (nnz(sourcehighout()==2)),...
    (nnz(sourcehighout()==1))];


%source�ɂ�����Other��target�ł̐���
sourceother=zeros(1,length(transother));
for i=1:length(transother)
    if nnz(targetdata{2,6}(:)==transother(i))==1
        sourceother(i)=3;
    elseif nnz(targetdata{2,7}(:)==transother(i))==1
        sourceother(i)=3;
    elseif nnz(targetdata{2,8}(:)==transother(i))==1
        sourceother(i)=2;
    else
        sourceother(i)=1;
    end
end
otherdata=categorical(sourceother(:),valueset,type,'Ordinal',true);
figure
pie(otherdata);
title('Other')

%source��other�̂����A������target��high,other,undetected���Z�o����
otherratio=[(nnz(sourceother()==3)),...
    (nnz(sourceother()==2)),...
    (nnz(sourceother()==1))];

cellprop={'highin','Highout','Other',strcat(source,'�̍זE��',target,'�łǂ̂悤�Ȑ�����1�����o2���̑�3��������'),0,0,0;...
    sourcehighin.',sourcehighout.',sourceother.',0,0,0,0;...
    highinratio,highoutratio,otherratio,'�����̎Z�o','��������','���̑�','�����o�̏���'};

A=repmat("highin",length(sourcehighin),1);
A=[A;repmat("highout",length(sourcehighout),1)];
A=[A;repmat("other",length(sourceother),1)];
cat=[sourcehighin.';sourcehighout.';sourceother.'];
valueset=[3,2,1];
type={'HC','Other','Undetected'};
cat=categorical(cat,valueset,type,'Ordinal',true);
cat={A,cat};

%% �������ۑ�����Ă��邩���ׂ�
len=length(sourcedata{2,1});
transGTE=zeros(length(targetdata{2,1}));    %source�̋@�\�I�������i�[�������́iROI�ԍ���target)

for i=1:len
    for j=1:len
        if sourcedata{2,1}(i,j)==1
            a=targetROI(sourceROI==i);
            b=targetROI(sourceROI==j);
            if (isnan(a)==0)&&(isnan(b)==0)
                transGTE(a,b)=1;
            end
        end
    end
end

k=1;
for i=1:length(targetdata{2,1})
    for j=1:length(targetdata{2,1})
        if (targetdata{2,1}(i,j)==1)&&(transGTE(i,j)==1)
            savedoutROI(k)=j;   %ROI�ԍ���target�̂��̂ł��邱�Ƃɒ���
            savedinROI(k)=i;
            k=k+1;
        end
    end
end

%% �ۑ����ꂽ�����Ɋւ��זE�̐����𒲂ׂ�
len=length(savedinROI);
sourcesavedin=zeros(1,len);
sourcesavedout=zeros(1,len);
for i=1:len
    sourcesavedin(i)=sourceROI(targetROI==savedinROI(i));
    sourcesavedout(i)=sourceROI(targetROI==savedoutROI(i));
end
savedcon={strcat(source,'in'),strcat(source,'out'),...
    strcat(target,'in'),strcat(target,'out'),'�ۑ�����Ă��������Ƃ��̐���';...
    sourcesavedin.',sourcesavedout.',savedinROI.',savedoutROI.',0};


%�g�����́cbothdhighROI,savedcon
sourceintype=zeros(1,len);
sourceouttype=zeros(1,len);
targetintype=zeros(1,len);
targetouttype=zeros(1,len);

for i=1:len
    if nnz(bothhighROI{2,1}()==savedcon{2,1}(i))~=0
        sourceintype(i)=1;  %1��saved high connectivity
    elseif nnz(sourcedata{2,8}==savedcon{2,1}(i))==0
        sourceintype(i)=2;  %2�͕Е��݂̂�high connectivity
    else
        sourceintype(i)=3;  %3�͗����Ƃ���other neuron
    end
end
for i=1:len
    if nnz(bothhighROI{2,1}()==savedcon{2,2}(i))~=0
        sourceouttype(i)=1;
    elseif nnz(sourcedata{2,8}==savedcon{2,2}(i))==0
        sourceouttype(i)=2;
    else
        sourceouttype(i)=3;
    end
end
for i=1:len
    if nnz(bothhighROI{2,2}()==savedcon{2,3}(i))~=0
        targetintype(i)=1;
    elseif nnz(targetdata{2,5}==savedcon{2,3}(i))~=0
        targetintype(i)=2;
    else
        targetintype(i)=3;
    end
end
for i=1:len
    if nnz(bothhighROI{2,2}()==savedcon{2,4}(i))~=0
        targetouttype(i)=1;
    elseif nnz(targetdata{2,5}==savedcon{2,4}(i))~=0
        targetouttype(i)=2;
    else
        targetouttype(i)=3;
    end
end

%�e�����̃^�C�v�̊������Z�o����
sourcecontype=zeros(1,6);
targetcontype=zeros(1,6);
SHCSHC=zeros();
for i=1:len
    if (sourceintype(i)==1)&&(sourceouttype(i)==1)
        sourcecontype(1)=sourcecontype(1)+1;    %SHC-SHC
        SHCSHC=[SHCSHC,i];                      %SHC-SHC����������savedcon�̍s�ԍ�
    elseif (sourceintype(i)==1)&&(sourceouttype(i)==2)
        sourcecontype(2)=sourcecontype(2)+1;    %SHC-HC
    elseif (sourceintype(i)==1)&&(sourceouttype(i)==3)
        sourcecontype(3)=sourcecontype(3)+1;    %SHC-Other
    elseif (sourceintype(i)==2)&&(sourceouttype(i)==1)
        sourcecontype(2)=sourcecontype(2)+1;    %SHC-HC
    elseif (sourceintype(i)==2)&&(sourceouttype(i)==2)
        sourcecontype(4)=sourcecontype(4)+1;    %HC-HC
    elseif (sourceintype(i)==2)&&(sourceouttype(i)==3)
        sourcecontype(5)=sourcecontype(5)+1;    %HC-Other
    elseif (sourceintype(i)==3)&&(sourceouttype(i)==1)
        sourcecontype(3)=sourcecontype(3)+1;    %SHC-Other
    elseif (sourceintype(i)==3)&&(sourceouttype(i)==2)
        sourcecontype(5)=sourcecontype(5)+1;    %HC-Other
    elseif (sourceintype(i)==3)&&(sourceouttype(i)==3)
        sourcecontype(6)=sourcecontype(6)+1;    %Other-Other
    end
end
SHCSHC(1)=[];

for i=1:len
    if (targetintype(i)==1)&&(targetouttype(i)==1)
        targetcontype(1)=targetcontype(1)+1;
    elseif (targetintype(i)==1)&&(targetouttype(i)==2)
        targetcontype(2)=targetcontype(2)+1;
    elseif (targetintype(i)==1)&&(targetouttype(i)==3)
        targetcontype(3)=targetcontype(3)+1;
    elseif (targetintype(i)==2)&&(targetouttype(i)==1)
        targetcontype(2)=targetcontype(2)+1;
    elseif (targetintype(i)==2)&&(targetouttype(i)==2)
        targetcontype(4)=targetcontype(4)+1;
    elseif (targetintype(i)==2)&&(targetouttype(i)==3)
        targetcontype(5)=targetcontype(5)+1;
    elseif (targetintype(i)==3)&&(targetouttype(i)==1)
        targetcontype(3)=targetcontype(3)+1;
    elseif (targetintype(i)==3)&&(targetouttype(i)==2)
        targetcontype(5)=targetcontype(5)+1;
    elseif (targetintype(i)==3)&&(targetouttype(i)==3)
        targetcontype(6)=targetcontype(6)+1;
    end
end

type={'SHC','HC','Other'}; %saved high connectivity(SHC)
valueset=[1,2,3];
sourceintype=categorical(sourceintype.',valueset,type,'Ordinal',true);
sourceouttype=categorical(sourceouttype.',valueset,type,'Ordinal',true);
targetintype=categorical(targetintype.',valueset,type,'Ordinal',true);
targetouttype=categorical(targetouttype.',valueset,type,'Ordinal',true);

savedcon{3,1}=sourceintype;
savedcon{3,2}=sourceouttype;
savedcon{3,3}=targetintype;
savedcon{3,4}=targetouttype;
%% ����p��dataframe������
savedcontype={'SHC-SHC','SHC-HC','SHC-Other',...
    'HC-HC','HC-Other','Other-Other'};
%source�̑S�����^�C�v
valueset=[8,6,5,4,3,2];
sourceallcontype=zeros(length(sourcedata{2,2}),2);
for i=1:length(sourcedata{2,2})
    for j=1:2
        if nnz(sourcedata{2,2}(i,j)==bothhighROI{2,1}())~=0
            sourceallcontype(i,j)=4;
        elseif nnz(sourcedata{2,2}(i,j)==sourcedata{2,5}())~=0
            sourceallcontype(i,j)=2;
        else
            sourceallcontype(i,j)=1;
        end
    end
end
sourceallcontype=sum(sourceallcontype,2);
allcontyperatio=[nnz(sourceallcontype()==8),...
    nnz(sourceallcontype()==6),...
    nnz(sourceallcontype()==5),...
    nnz(sourceallcontype()==4),...
    nnz(sourceallcontype()==3),...
    nnz(sourceallcontype()==2)];

sourceallcontype=categorical(sourceallcontype,valueset,savedcontype,'Ordinal',true);
A=repmat("�S����",length(sourceallcontype),1);

%�ۑ����ꂽ�����^�C�v(source)
valueset=[1,2,3,4,5,6];
sourcesavedcontype=horzcat(ones(1,sourcecontype(1)),...
    repmat(2,1,sourcecontype(2)),...
    repmat(3,1,sourcecontype(3)),...
    repmat(4,1,sourcecontype(4)),...
    repmat(5,1,sourcecontype(5)),...
    repmat(6,1,sourcecontype(6)));
sourcesavedcontype=categorical(sourcesavedcontype.',valueset,savedcontype,'Ordinal',true);
subject=repmat("�ۑ�����",length(sourcesavedcontype),1);
contype1=[sourceallcontype;sourcesavedcontype];
contype2=[A;subject];
contype={contype1,contype2};

%% �f�[�^�̕ۑ�
DIV=strcat(source,'to',target);
inf={'DIV','inf';...
    DIV,datainf};
celltype={'��������','���̑�','�����o'};
filename=strcat(source,'to',target,'.xlsx');
savefile=strcat(dir,filename);
writematrix(bothhighROI{1,1},savefile,'Sheet',1,'Range','A1');
writematrix(bothhighROI{1,2},savefile,'Sheet',1,'Range','B1');
writematrix(bothhighROI{1,3},savefile,'Sheet',1,'Range','C1');
writematrix(bothhighROI{2,1},savefile,'Sheet',1,'Range','A2');
writematrix(bothhighROI{2,2},savefile,'Sheet',1,'Range','B2');

%source��target�ł̐���
subject={'sourcetype','targettype'};
writecell(subject,savefile,'Sheet',2,'Range','A1');
writematrix(cat{1,1},savefile,'Sheet',2,'Range','A2');
writematrix(cat{1,2},savefile,'Sheet',2,'Range','B2');
writecell(inf,savefile,'Sheet',2,'Range','C1');
writematrix(cellprop{1,1},savefile,'Sheet',2,'Range','F1');
writematrix(cellprop{1,2},savefile,'Sheet',2,'Range','I1');
writematrix(cellprop{1,3},savefile,'Sheet',2,'Range','L1');
writecell(celltype,savefile,'Sheet',2,'Range','F2');
writecell(celltype,savefile,'Sheet',2,'Range','I2');
writecell(celltype,savefile,'Sheet',2,'Range','L2');
writematrix(cellprop{3,1},savefile,'Sheet',2,'Range','F3');
writematrix(cellprop{3,2},savefile,'Sheet',2,'Range','I3');
writematrix(cellprop{3,3},savefile,'Sheet',2,'Range','L3');
writecell(inf,savefile,'Sheet',2,'Range','o2');

%�����^�C�v
subject={'category';'�S����';'�ۑ�����'};
writecell(subject,savefile,'Sheet',3,'Range','F1');
writecell(savedcontype,savefile,'Sheet',3,'Range','G1');
writematrix(allcontyperatio,savefile,'Sheet',3,'Range','G2');
writematrix(sourcecontype,savefile,'Sheet',3,'Range','G3');
writecell(inf,savefile,'Sheet',3,'Range','M1');
%����p�̃f�[�^
subject=["Type","category"];
writematrix(subject,savefile,'Sheet',3,'Range','A1');
writematrix(contype{1,1},savefile,'Sheet',3,'Range','A2');
writematrix(contype{1,2},savefile,'Sheet',3,'Range','B2');
writecell(inf,savefile,'Sheet',3,'Range','C1');
