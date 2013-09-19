values=zeros(2,1);
 objectivenums = [1 2 3 4];
 
 k=0;
 labels = cell({});
 options = struct;

 for i = 1:2
     for j=1:2
         options.model_id = i;
         options.objective=j;
         tempres=runsim(options);
         k=k+1;
         values(k,1)=tempres.Fmin_mindistance;
         values(k,2)=tempres.Fmin_maxdistance;
         labels{k} = [tempres.modelname,' ',tempres.objectivename];

     end
 end
 

range = values(:,2)-values(:,1);
values_diff = [values(:,1) range];


figure('name','objective comparison')
bh = bar(values_diff,'stacked');
%set(gca,'ylim',[1 1100])
set(gca,'xtick',[1 2 3 4],'xticklabel',labels)
set(bh(1),'FaceColor','none','EdgeColor','none')

%This works!