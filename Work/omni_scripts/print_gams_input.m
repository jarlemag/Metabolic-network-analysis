function print_gams_input(file,Aeb,Aec,Agb,Agc,Alb,Alc,be,bg,bl,lbb,lbc,ubb,ubc,cb,cc,intvarnames,contvarnames);
%PRINT_GAMS_INPUT
%
% print_gams_input(file,Aeb,Aec,Agb,Agc,Alb,Alc,be,bg,bl,lbb,lbc,ubb,ubc,cb,cc,intvarnames,contvarnames);
%
% Markus Herrgard 5/1/03

[ne,nb] = size(Aeb);
[ne,nc] = size(Aec);
[ng,nb] = size(Agb);
[nl,nb] = size(Alb);

fid = fopen(file,'w');

fprintf(fid,'Sets\n\n');

fprintf(fid,'ie Rows in Ae /ce1*ce%d/\n',ne);
fprintf(fid,'ig Rows in Ag /cg1*cg%d/\n',ng);
fprintf(fid,'il Rows in Al /cl1*cl%d/\n',nl);

% Binary variable names
fprintf(fid,'jb Bin vars /');
if (length(intvarnames) == 0)
    fprintf(fid,'%s\n',['y1*y' num2str(nb) '/']);
    for i = 1:nb
        bvarname{i} = ['y' num2str(i)];
    end
else
for i = 1:nb
    if (i <= length(intvarnames))
        if (i == 1) 
            fprintf(fid,'%s',intvarnames{i});
        elseif (mod(i,100) == 0)
            fprintf(fid,'\n%s',intvarnames{i});
        else
            fprintf(fid,',%s',intvarnames{i});
        end
        bvarname{i} = intvarnames{i};
    else
        if (i == 1)
            fprintf(fid,'y%d',i);
        elseif (mod(i,100) == 0)
            fprintf(fid,'\ny%d',i);
        else
            fprintf(fid,',y%d',i);
        end
        bvarname{i} = ['y' num2str(i)];
    end
end
fprintf(fid,'/\n');
end
% Continous variable names
fprintf(fid,'jc Cont vars /');
if (length(contvarnames) == 0)
    fprintf(fid,'%s\n',['x1*x' num2str(nc) '/']);
    for i = 1:nc
        cvarname{i} = ['x' num2str(i)];
    end
else
for i = 1:nc
    if (i <= length(contvarnames))
        if (i == 1)
            fprintf(fid,'%s',contvarnames{i});
        elseif (mod(i,100) == 0)
            fprintf(fid,'\n%s',contvarnames{i});
        else
            fprintf(fid,',%s',contvarnames{i});
        end
        cvarname{i} = contvarnames{i};
    else
        if (i == 1)
            fprintf(fid,'x%d',i);
        elseif (mod(i,100) == 0)
            fprintf(fid,'\nx%d',i);
        else
            fprintf(fid,',x%d',i);
        end
        cvarname{i} = ['x' num2str(i)];
    end
end
fprintf(fid,'/;\n\n');
end

fprintf(fid,'Parameters\n\n');
         
% A eq bin
fprintf(fid,'Aeb(ie,jb)\n/ ');
[ii,jj,val] = find(Aeb);
for i = 1:length(ii)
    fprintf(fid,'ce%d.%s %f\n',ii(i),bvarname{jj(i)},val(i));
end
fprintf(fid,'/\n');

% A eq cont
fprintf(fid,'Aec(ie,jc)\n/ ');
[ii,jj,val] = find(Aec);
for i = 1:length(ii)
    fprintf(fid,'ce%d.%s %f\n',ii(i),cvarname{jj(i)},val(i));
end
fprintf(fid,'/\n');

% A ge bin
fprintf(fid,'Agb(ig,jb)\n/ ');
[ii,jj,val] = find(Agb);
for i = 1:length(ii)
    fprintf(fid,'cg%d.%s %f\n',ii(i),bvarname{jj(i)},val(i));
end
fprintf(fid,'/\n');

% A ge cont
fprintf(fid,'Agc(ig,jc)\n/ ');
[ii,jj,val] = find(Agc);
for i = 1:length(ii)
    fprintf(fid,'cg%d.%s %f\n',ii(i),cvarname{jj(i)},val(i));
end
fprintf(fid,'/\n');

% A le bin
fprintf(fid,'Alb(il,jb)\n/ ');
[ii,jj,val] = find(Alb);
for i = 1:length(ii)
    fprintf(fid,'cl%d.%s %f\n',ii(i),bvarname{jj(i)},val(i));
end
fprintf(fid,'/\n');

% A le cont
fprintf(fid,'Alc(il,jc)\n/ ');
[ii,jj,val] = find(Alc);
for i = 1:length(ii)
    fprintf(fid,'cl%d.%s %f\n',ii(i),cvarname{jj(i)},val(i));
end
fprintf(fid,'/\n');

% b eq
fprintf(fid,'be(ie)\n/ ');
for i = 1:length(be)
    if (be(i) ~= 0)
        fprintf(fid,'ce%d %f\n',i,be(i));
    end
end
fprintf(fid,'/\n');

% b ge
fprintf(fid,'bg(ig)\n/ ');
for i = 1:length(bg)
    if (bg(i) ~= 0)
        fprintf(fid,'cg%d %f\n',i,bg(i));
    end
end
fprintf(fid,'/\n');

% b le
fprintf(fid,'bl(il)\n/ ');
for i = 1:length(bl)
    if (bl(i) ~= 0)
        fprintf(fid,'cl%d %f\n',i,bl(i));
    end
end
fprintf(fid,'/\n');

% lb bin
fprintf(fid,'lbb(jb)\n/ ');
for i = 1:length(lbb)
    if (lbb(i) ~= 0)
        fprintf(fid,'%s %f\n',bvarname{i},lbb(i));
    end
end
fprintf(fid,'/\n');

% lb cont
fprintf(fid,'lbc(jc)\n/ ');
for i = 1:length(lbc)
    if (lbc(i) ~= 0)
        fprintf(fid,'%s %f\n',cvarname{i},lbc(i));
    end
end
fprintf(fid,'/\n');

% ub bin
fprintf(fid,'ubb(jb)\n/ ');
for i = 1:length(ubb)
    if (ubb(i) ~= 0)
        fprintf(fid,'%s %f\n',bvarname{i},ubb(i));
    end
end
fprintf(fid,'/\n');

% ub cont
fprintf(fid,'ubc(jc)\n/ ');
for i = 1:length(ubc)
    if (ubc(i) ~= 0)
        fprintf(fid,'%s %f\n',cvarname{i},ubc(i));
    end
end
fprintf(fid,'/\n');

% c bin
fprintf(fid,'cb(jb)\n/ ');
for i = 1:length(cb)
    if (cb(i) ~= 0)
        fprintf(fid,'%s %f\n',bvarname{i},cb(i));
    end
end
fprintf(fid,'/\n');

% c cont
fprintf(fid,'cc(jc)\n/ ');
for i = 1:length(cc)
    if (cc(i) ~= 0)
        fprintf(fid,'%s %f\n',cvarname{i},cc(i));
    end
end
fprintf(fid,'/\n');

fclose(fid);