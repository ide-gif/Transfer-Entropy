%% AAV�̂Ƃ��͂������I
%�ǂݍ��ރf�[�^
clearvars
data = readmatrix('C:\Ide_AQUACOSMOS\191228AAV-3\area1\DIV11\Results.csv');
DIV='DIV11';                     %DIV
datainf='191228-3-area1';       %�����肵���̂�
dir='C:\analysisresults\AAV\';  %AAV��immunostaining
inf={'DIV','inf';DIV,datainf};
savedir=strcat(dir,datainf,'\',DIV,'\'); %�ۑ�������f�B���N�g��


%% immunostaining�̂Ƃ��͂������I�I
clearvars
data = readmatrix('C:\Ide_AQUACOSMOS\191016\HippocampusDIV9-1\Results.csv');
DIV='ERK';                     %CREB,ERK,4E-BP1�ȂǗp�����^���p�N��
datainf='�O��';       %�O���A�����A���
dir='C:\analysisresults\immunostaining\';  %AAV��immunostaining
inf={'DIV','inf';DIV,datainf};
savedir=strcat(dir,datainf,'\',DIV,'\'); %�ۑ�������f�B���N�g��
