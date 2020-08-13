%% AAVのときはこっち！
%読み込むデータ
clearvars
data = readmatrix('C:\Ide_AQUACOSMOS\191228AAV-3\area1\DIV11\Results.csv');
DIV='DIV11';                     %DIV
datainf='191228-3-area1';       %いつ測定したのか
dir='C:\analysisresults\AAV\';  %AAVかimmunostaining
inf={'DIV','inf';DIV,datainf};
savedir=strcat(dir,datainf,'\',DIV,'\'); %保存をするディレクトリ


%% immunostainingのときはこっち！！
clearvars
data = readmatrix('C:\Ide_AQUACOSMOS\191016\HippocampusDIV9-1\Results.csv');
DIV='ERK';                     %CREB,ERK,4E-BP1など用いたタンパク質
datainf='前期';       %前期、中期、後期
dir='C:\analysisresults\immunostaining\';  %AAVかimmunostaining
inf={'DIV','inf';DIV,datainf};
savedir=strcat(dir,datainf,'\',DIV,'\'); %保存をするディレクトリ
