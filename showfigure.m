%% �e��ݒ�
fontsize=12;    %���̃t�H���g�T�C�Y
axisfontsize=16;    %�����x���̃t�H���g�T�C�Y
linewid=2;  %臒l�̐��̑���
%% high connectivity neuron�̊���
type={'Other','HC'};
valueset=[1,2];
a=[length(otherROI),length(highROI)];
a=zeros(1,length(otherROI));
a=a()+1;
b=zeros(1,length(highROI));
b=b()+2;
c=[a,b];
plotdata=categorical(c,valueset,type,'Ordinal',true);
pie(plotdata);
filename='HCratio.fig';
savefile=strcat(savedir,filename);
savefig(savefile);
%% GTE�̃q�X�g�O����
figure
nbin=ceil(1+log2(numcell*numcell-numcell));
histogram(GTE,nbin,'Normalization','probability');
ax=gca;
ax.FontSize=fontsize;
xlabel('TE score','FontSize',axisfontsize);
ylabel('Probability','FontSize',axisfontsize);
hold on
xline(GTEsort(shre),'--r',...
    'LineWidth',linewid);
hold off
filename='GTEhistogram.fig';
savefile=strcat(savedir,filename);
savefig(savefile);

%% GTE�̃q�[�g�}�b�v
figure
for i=1:numcell
    GTE(i,i)=0;
 end
imagesc(GTE);
colormap(flipud(hot));
ax=gca;
ax.FontSize=fontsize;
xlabel('Cell Number','FontSize',axisfontsize);
ylabel('Cell Number','FontSize',axisfontsize);
colorbar;
filename='GTEheatmap.fig';
savefile=strcat(savedir,filename);
savefig(savefile);

%% FCGTE�̃q�[�g�}�b�v
figure
figFCGTE=double(FCGTE);
imagesc(FCGTE);
colormap(flipud(hot));
ax=gca;
ax.FontSize=fontsize;
xlabel('Cell Number','FontSize',axisfontsize);
ylabel('Cell Number','FontSize',axisfontsize);
filename='Connectivitytable.fig';
savefile=strcat(savedir,filename);
savefig(savefile);

%% in�y��out�̃q�X�g�O����
figure
nbin=ceil(1+log2(numcell));
tiledlayout(1,2);
nexttile
histogram(norin,nbin,'Normalization','probability');
ax=gca;
ax.FontSize=fontsize;
xlabel({'Normarized','in-degree'},'FontSize',axisfontsize);
ylabel('Probability','FontSize',axisfontsize);
hold on
xline(1,'--r',...
    'LineWidth',linewid);
hold off
nexttile
histogram(norout,nbin,'Normalization','probability');
ax=gca;
ax.FontSize=fontsize;
xlabel({'Normarized','out-degree'},'FontSize',axisfontsize);
ylabel('Probability','FontSize',axisfontsize);
hold on
xline(1,'--r',...
    'LineWidth',linewid);
hold off
filename='inout-degree.fig';
savefile=strcat(savedir,filename);
savefig(savefile);

%% ���X�^�[�v���b�g
figure
imagesc(NA.');
colormap(flipud(hot));
xticks([1,numframe-1]);
yticks([1,numcell]);
xticklabels({'0',round(Time(numframe))});
yticklabels({'1',numcell});
ax=gca;
ax.FontSize=fontsize;
xlabel('Time (min)','FontSize',axisfontsize);
ylabel('Cell Number','FontSize',axisfontsize);
filename='rasterplot.fig';
savefile=strcat(savedir,filename);
savefig(savefile);

%% �l�b�g���[�N�̊���
figure
plot(Time,AllInt,'r');
ax=gca;
ax.FontSize=fontsize;
xlabel('Time (min)','FontSize',axisfontsize);
ylabel('F_{Ave}','FontSize',axisfontsize);
filename='Networkintensity.fig';
savefile=strcat(savedir,filename);
savefig(savefile);

figure
b=bar(newras*100/numcell);
xticklabels({});
ylim([0,100]);
b.FaceColor='r';
ax=gca;
ax.FontSize=fontsize;
ylabel('Active cells (%)','FontSize',axisfontsize);
filename='ActiveCell.fig';
savefile=strcat(savedir,filename);
savefig(savefile);
