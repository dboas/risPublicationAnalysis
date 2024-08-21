%% 
% Load Records
% Note that I sometimes just export the current year and previous year. 
% I indicate this with a _(YYYY-1).ris at the end of the file.
% Or if it is just hat year then I put YYYY.ris.
% I can reconstruct all prior years easily if needed either by cutting and
% pasting from prior exports or just exporting new from readcube

rec = loadRIS('export_2024-8-20_2017.ris', 1);


%%
% Adjaceny Matrix for Desired Authors
yearTarget = 2024;

% What do I do about Faculty who have left NPC?
% I guess I can leave them at least for a few years since they're
% interactions take awhile to publish.
% But I do ever remove them, then it messes up things if I go back in
% time... this would suggest I need a year range for membership for each author... hmm...
namesAuthors = [{'Awad, L'},{'Bifano, T'},{'Bigio, I'},{'Boas, D'},{'Chen, A'},... %5
    {'Chen, J'},{'Cheng, J'},{'Chung, D'},{'Cronin-Golomb, A'},{'Cruz-Martin, A'},... %10
    {'Davison, I'},{'Demas, J'},{'Dennis, A'},{'Devor, A'},{'Economo, M'},{'Emili, A'},...%16 {'Fawcett, H'},{'Ferre, C'},...
    {'Gabel, C'},{'Gill, S'},... % 18   {'Greer, D'},...
    {'Goldstein, L'},{'Han, X'},{'Hasselmo, M'},{'Howe, M'},{'Kilic, K'},{'Kiran, S'},... %  24  {'Kumar, D'},...
    {'Lewis, L'},{'Ling, S'},{'Mertz, J'},{'Ndao, A'},{'Nia, H'},... % 29   {'O''Shea, T'},...
    {'Ramachandran, S'},{'Ramirez, S'},{'Roblyer, D'},{'Sander, M'},... % 33
    {'Scott, B'},{'Sen, K'},{'Somers, D'},{'Stern, C'},{'Stern, R'},{'Tager-Flusberg, H'},... % 39
    {'Tian, L'},{'White, J'},{'Wolozin, B'},... % 42
    {'Yucel, M'},... %43
    {'Cheng, X'},{'Ferre, C'},{'DePasquale, B'},{'Gavornik, J'},{'Kumar, D'},{'O''Shea, T'},... %49
    {'Russek, S'},{'Thunemann, M'},{'Wallace, M'},{'Wang, T'},{'Yang, C'},{'Younger, M'},... %55
    {'Stangl, M','Ellis, T','Greer, D','Luebke, J','Rosene, D','Stephen, E','TCW, J'},... %62
    {'Yücel, M'},{'Cruz-Martín, A'},{'Kılıç, K'},{'Chen, I'}];

% second column of duplicates must go in descending order. And duplicates must all be at end
duplicateNames = [5 66;... % Chen
                  23 65;... % Kilic
                  10 64;... % Cruz-Martin
                  43 63];   % Yucel

% make sure I manually delete abstracts and other non-sense inreadcube
% before exporting the ris file
excludeJournals = {'bioRxiv','Res. Sq.'};

[AuthorAdj,AuthorNames,AuthorNumPubs,AuthorNumCollabPubs, rec1] = loadRIS_Adjacency_NPC(rec,namesAuthors,yearTarget,excludeJournals);

% remove duplicate names
lstCollabPubs = find(rec1.isCollaborative==1);
for ii=1:size(duplicateNames,1)
    AuthorAdj(duplicateNames(ii,1),:) = AuthorAdj(duplicateNames(ii,1),:) + AuthorAdj(duplicateNames(ii,2),:);
    AuthorAdj(:,duplicateNames(ii,1)) = AuthorAdj(:,duplicateNames(ii,1)) + AuthorAdj(:,duplicateNames(ii,2));
    AuthorNumCollabPubs(duplicateNames(ii,1)) = AuthorNumCollabPubs(duplicateNames(ii,1)) + AuthorNumCollabPubs(duplicateNames(ii,2));
end
for ii=1:size(duplicateNames,1)
    AuthorAdj(duplicateNames(ii,2),:) = [];
    AuthorAdj(:,duplicateNames(ii,2)) = [];
    AuthorNumCollabPubs(duplicateNames(ii,2)) = [];
    
    for jj=1:length(lstCollabPubs)
        [lia,locb] = ismember(duplicateNames(ii,2),rec1.CollabAuthors{lstCollabPubs(jj)});
        if lia==1
            rec1.CollabAuthors{lstCollabPubs(jj)}(locb) = duplicateNames(ii,1);
        end
    end
end

nAuthors = size(AuthorAdj,1);

%%
% histogram number of collaborations
figure(1)
hist(sum(full(AuthorAdj)>0,2),[0:1:20])
xlim([-1 18])
set(gcf,'color',[1 1 1])
set(gca,'fontsize',20)
ylabel('# Faculty')
xlabel('# Collaborations')




%%
% List Authors and who they have published with

% number of collaborative papers each year
disp(sprintf('\n\n'));
lstCollabPubs = find(rec1.isCollaborative==1);
nPubsPerYear = [];
for iY = 2017:2024
    nPubsPerYear(end+1) = length(find(rec1.Year(lstCollabPubs)==iY));
    disp(sprintf('Year %4d had %2d collaborative publications',iY,nPubsPerYear(end)))
end
disp(sprintf('\n\n'));

% report names based on number of collaborative papers
% numCPub = sum(full(AuthorAdj),2); % FIXME BUT THIS DOUBLE COUNTS PAPERS WITH 3 OR MORE COLLABORATIVE AUTHORS
numCPub = AuthorNumCollabPubs;
[foo,isort] = sort(numCPub,'descend');

for ii=1:nAuthors
    lst = find(AuthorAdj(isort(ii),:)>0);
    foos = sprintf('%s(%d,%d) \t- ',AuthorNames{isort(ii)}, AuthorNumCollabPubs(isort(ii)), length(lst) );
    
    [foo,isort2] = sort( full(AuthorAdj(isort(ii),lst)), 'descend' );
    
    for jj=1:length(lst)
        foos = sprintf( '%s%s(%d); ',foos,AuthorNames{lst(isort2(jj))},full(AuthorAdj(isort(ii),lst(isort2(jj)))) );
    end
    disp(foos)
end

%%
% list collaborative publications

% number of collaborative papers each year
disp(sprintf('\n\n'));
lstCollabPubs = find(rec1.isCollaborative==1);
nPubsPerYear = [];
for iY = 2017:2024
    nPubsPerYear(end+1) = length(find(rec1.Year(lstCollabPubs)==iY));
    disp(sprintf('Year %4d had %2d collaborative publications',iY,nPubsPerYear(end)))
end


disp(sprintf('\n\n'))
disp('Collaborative Publications')
disp('==========================');
for ii=1:length(lstCollabPubs)
    
    strAuthors = '';
    for jj=1:length(rec1.CollabAuthors{lstCollabPubs(ii)})
        kk = rec1.CollabAuthors{lstCollabPubs(ii)}(jj);
        strAuthors = sprintf('%s, %s', strAuthors, replace(AuthorNames{kk},',',''));
    end
    
    % strAuthors(3:end)
%    disp( sprintf('%d) %s, "%s", %s, %d.', ii, rec1.AuthorList{lstCollabPubs(ii)}(3:end), rec1.Title{lstCollabPubs(ii)}, rec1.Journal{lstCollabPubs(ii)}, rec1.Year(lstCollabPubs(ii)) ));
    disp( sprintf('%s, "%s", %s, %d.', rec1.AuthorList{lstCollabPubs(ii)}(3:end), rec1.Title{lstCollabPubs(ii)}, rec1.Journal{lstCollabPubs(ii)}, rec1.Year(lstCollabPubs(ii)) ));
    
end

%%
% list collaborative publications sorted by key authors

disp(sprintf('\n\n'))
disp('Collaborative Publications')
disp('==========================');

lstCollabPubsTmp = lstCollabPubs;

for ii=20 %43 %41

    lstCollabPubsTmp = lstCollabPubs;
    
    lstCAuthors = setdiff(find(AuthorAdj(ii,:)>0),ii); 
    
    for jj=1:length(lstCAuthors)
        
        % list papers with 2 collaborative authors
        ll = 0;
        lstRemove = [];
        for kk=1:length(lstCollabPubsTmp)
            foo = rec1.CollabAuthors{lstCollabPubsTmp(kk)};
            if length(foo)==2 && ismember(ii,foo) 
                if ismember(lstCAuthors(jj),foo)
                    ll = ll + 1;
                    if ll==1
                        disp( sprintf('\n%s, %s', replace(AuthorNames{ii},',',''), replace(AuthorNames{lstCAuthors(jj)},',','') ) )
                        disp( '-----------------' );
                    end
                    disp( sprintf('%d) %s, "%s", %s, %d.', ll, rec1.AuthorList{lstCollabPubsTmp(kk)}(3:end), rec1.Title{lstCollabPubsTmp(kk)}, rec1.Journal{lstCollabPubsTmp(kk)}, rec1.Year(lstCollabPubsTmp(kk)) ));
                    lstRemove = [lstRemove kk];
                end
            end        
        end
        lstCollabPubsTmp( lstRemove ) = [];
        
    end
        
        % list papers with 3 or collaborative authors
        ll = 0;
        lstRemove = [];
        for kk=1:length(lstCollabPubsTmp)
            foo = rec1.CollabAuthors{lstCollabPubsTmp(kk)};
            if length(foo)>2 && ismember(ii,foo)>0 
                ll = ll + 1;
                if ll==1
                    disp( sprintf('\n%s et al', replace(AuthorNames{ii},',','') ) )
                    disp( '-----------------' );
                end
                disp( sprintf('%d) %s, "%s", %s, %d.', ll, rec1.AuthorList{lstCollabPubsTmp(kk)}(3:end), rec1.Title{lstCollabPubsTmp(kk)}, rec1.Journal{lstCollabPubsTmp(kk)}, rec1.Year(lstCollabPubsTmp(kk)) ));
                lstRemove = [lstRemove kk];
            end        
        end
        lstCollabPubsTmp( lstRemove ) = [];
        
       
end

%%
% Load adjacency for each year

for ii=1:2
    [AuthorAdjAll{ii},AuthorNames,boo,boo1,rec1] = loadRIS_Adjacency_NPC(rec,namesAuthors,2022+ii);
    
    % remove duplicate names
    lstCollabPubs = find(rec1.isCollaborative==1);
    for jj=1:size(duplicateNames,1)
        AuthorAdjAll{ii}(duplicateNames(jj,1),:) = AuthorAdjAll{ii}(duplicateNames(jj,1),:) + AuthorAdjAll{ii}(duplicateNames(jj,2),:);
        AuthorAdjAll{ii}(:,duplicateNames(jj,1)) = AuthorAdjAll{ii}(:,duplicateNames(jj,1)) + AuthorAdjAll{ii}(:,duplicateNames(jj,2));
    end
    for jj=1:size(duplicateNames,1)
        AuthorAdjAll{ii}(duplicateNames(jj,2),:) = [];
        AuthorAdjAll{ii}(:,duplicateNames(jj,2)) = [];
    end
    
    foo = sum(AuthorAdjAll{ii},2);
    lstAll{ii} = find(foo>0);
end


%%
% create a graph from adjacency matrix
iYear = 2;

lstA = lstAll{iYear};
AuthorAdj = AuthorAdjAll{iYear};

nAuthors = size(AuthorAdjAll{iYear},1);
nAuthorsA = length(lstA);

g=graph(AuthorAdjAll{iYear}(lstA,lstA));

groups = zeros(nAuthors,1);
groups(lstA) = conncomp(g)'; % get the connected components in the graph
% figure(1)
% plot(g)

% get edges and weights
foo = g.Edges;
boo = foo.Variables;
edges = lstA(boo(:,1:2));
weights = boo(:,3);

% all nodes have 1/r^2 repulsion and 1/r attraction
% attraction with other nodes is weighted by their number of pubs (their
%     adjacency weight)
% repulsions are all unit weight
% Ultimately I want node circle diameters to be given by number of pubs for
%     that node
% I need a 1/r attraction at the origin to center everything

%%
% start with positions for each connected group around the unit circle

nGroups = length(unique(groups));

pos = 0.2*randn(nAuthors,2) + 3*[cos(2*3.14159*groups/nGroups) sin(2*3.14159*groups/nGroups)];

figure(1)
clf

%%
% evolve the node positions

if 0 % randomize positions a little
    pos = pos + rand(size(pos)) * 1;
end


for ii = 1:100
    nIterMax = 10;
    pos = graphForces( pos, AuthorAdjAll{iYear}, lstA, nIterMax );
    
    figure(2)
    plot( pos(lstA,1), pos(lstA,2), '.' )
    set(gca,'fontsize',16)
    set(gcf,'color',[1 1 1])
    if 0
        set(gca,'ytick',[])
        set(gca,'xtick',[])
    end
    hold on
    for iE = 1:size(edges,1)
        he = plot( [pos(edges(iE,1),1) pos(edges(iE,2),1)], [pos(edges(iE,1),2) pos(edges(iE,2),2)] );
        switch weights(iE)
            case 1
                set(he,'color',[1 1 1]*0.5)
                set(he,'linewidth',1)
                set(he,'linestyle','--')
            case 2
                set(he,'color',[1 1 1]*0)
                set(he,'linewidth',1)
                set(he,'linestyle','-')
            case 3
                set(he,'color',[0 1 0])
                set(he,'linewidth',1)
            case {4,5}
                set(he,'color',[0 1 1])
                set(he,'linewidth',2)
            case {6,7}
                set(he,'color',[0 0 1])
                set(he,'linewidth',3)
            case {8,9}
                set(he,'color',[1 0 1])
                set(he,'linewidth',4)
            otherwise
                set(he,'color',[1 0 0])
                set(he,'linewidth',5)
        end
%        set(he,'linewidth',log(weights(iE))+1)
    end
    
    for jj=1:nAuthors
        if ismember(jj,lstA)
            if 1
                he = text( pos(jj,1), pos(jj,2), sprintf('%s',AuthorNames{jj}) );
            else
                he = text( pos(jj,1), pos(jj,2), sprintf('%s',AuthorNames{jj}) );
            end
            switch length(find(AuthorAdjAll{iYear}(jj,:)>0))
                case {1,2}
                    set(he,'fontsize',8)
                case {2,3}
                    set(he,'fontsize',10)
                case {4,5}
                    set(he,'fontsize',12)
                case {6,7}
                    
                    set(he,'fontsize',14)
                otherwise
                    set(he,'fontsize',16)
            end
        end
    end
    
    hold off
    title( sprintf('Year %d',2016+iYear) )
    xlabel( sprintf('nIter = %d',ii) )
    
    pause(0.01)
end

%%
% Increment the Adjacency Matrix
iYear = iYear + 1;

% new nodes appearing
lstNew = setdiff(lstAll{iYear},lstAll{iYear-1});
%AuthorAdj = AuthorAdj1996;
lstA = lstAll{iYear};

% set positions of new nodes based on mean of connected nodes
pos(lstNew,:) = 0;
for ii=1:length(lstNew)
    lstJJ = find(AuthorAdjAll{iYear}(lstNew(ii),:)>0);
    lstJJ = setdiff(lstJJ,lstNew);
    if length(lstJJ)>1
        pos(lstNew(ii),:) = mean(pos(lstJJ,:),1) + 0.1*randn(1,2);
    elseif length(lstJJ)==1
        pos(lstNew(ii),:) = mean(pos(lstJJ,:),1) + 0.1*randn(1,2);
    else
        pos(lstNew(ii),:) = [6 3] + 0.1*randn(1,2);
    end
end

% update edges and weights
g=graph(AuthorAdjAll{iYear}(lstA,lstA));
foo = g.Edges;
boo = foo.Variables;
edges = lstA(boo(:,1:2));
weights = boo(:,3);

 