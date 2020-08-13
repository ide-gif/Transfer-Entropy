%% �f�[�^�t�@�C���̓ǂݍ���
clearvars
dir1='C:\analysisresults\immunostaining\����\ERK';
filename='\191228DIV12-area2-1_ERK';
IHCname='pERK/ERK';   %�p�����R�̃O���t�By���̖��O�ɂȂ�

filename1=strcat(dir1,filename,'.xlsx');
datamat='\data.mat';
dir=strcat(dir1,datamat);
%GTE�f�[�^��g�ݍ���
file=load(dir);
data=file.GTEdata;
topROI=data{2,13}.';
highROI=data{2,5};
otherROI=data{2,8};
for i=1:length(highROI)
    if nnz(highROI(i)==topROI)==1
        highROI(i)=nan;
    end
end
highROI=rmmissing(highROI);
%% �G�N�Z���t�@�C����top,high,other��ROI�ԍ�����������
Top=array2table(topROI,'VariableNames',{'TC'});
High=array2table(highROI,'VariableNames',{'HC'});
other=array2table(otherROI,'VariableNames',{'Other'});
writetable(Top,filename1,'Range','A1');
writetable(High,filename1,'Range','B1');
writetable(other,filename1,'Range','C1');

%% ���ʂ�ǂݍ����beeswarmplot�ŏo��

IHCdata=readtable(filename1,...
    'Sheet',2,...
    'Range','D:F');

plotdata={IHCdata.TC,IHCdata.HC,IHCdata.Other};
figure
IHC=plotSpread(plotdata,...
    'xnames',{'TC','HC','Other'},...
    'ylabel',IHCname,...
    'distributioncolors','black');
set(findall(1,'type','line'),...
    'markerSize',16);
%set(findall(1,'type','text'),...
%    'fontSize',16);
filename2=strcat(dir1,filename,'.fig');

savefig(filename2);
