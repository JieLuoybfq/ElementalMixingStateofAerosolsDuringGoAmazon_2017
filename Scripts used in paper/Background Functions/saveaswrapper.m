W = whos;
structname = 'wrapper'; %input('What do you want to name the struct?','s');
% numexcep = input('How many variables would you like to exclude?');
% if numexcep > 0
%     exceptions = cell(1,numexcep);
%     for j = 1:numexcep
%     exceptions{j} = input('Which variables do you want leave out?','s');
%     end
% end
for i = 1:length(W)

        eval([structname,'.(W(i).name) = eval(W(i).name);']);
        clear(W(i).name)

end
varlist = {'i','W','structname'};
clear(varlist{:});
clear varlist