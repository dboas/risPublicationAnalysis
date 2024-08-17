% This calculates how many papers an author has published with another
% author. 
%
%
%
%  TO DO
% Can I do it by year
%rec0.uni
%
%
function [AuthorAdj,AuthorNames,AuthorNumPubs,AuthorNumCollabPubs,rec] = loadRIS_Adjacency_NPC( rec0, namesAuthors, yearTarget, excludeJournals )


nRec = 0;

    

% and sort in descending order
nAuthors = length(namesAuthors);

yearLst = unique(rec0.Year);
yearMin = min(yearLst);
yearMax = max(yearLst);
nYears = yearMax - yearMin + 1;
if ~exist('yearTarget')
    yearTarget = yearMax;
end

if ~exist('excludeJournals')
    excludeJournals = {};
end

idxOrder = [1:length(namesAuthors)];

% old legacy code... must have re-ordered things
uniqueAuthors = cell(1,nAuthors);
for ii=1:nAuthors
    uniqueAuthors{ii} = namesAuthors{ii}; %rec0.uniqueAuthors{lstAuthors(idxOrder(ii))};
end
AuthorNames = uniqueAuthors;

AuthorNumPubs = zeros(nAuthors,1);
AuthorNumCollabPubs = zeros(nAuthors,1);

AuthorAdj = sparse(nAuthors,nAuthors);
uniqueAuthorListRec = [];

% Load the file to create adjacency matrix
fid = fopen( rec0.filenm, 'r');
tline = fgetl( fid );
yearLast = [];
nYearNotDefined = 0;
while ischar(tline)
    if length(tline)>5
        % Start of a new record
        if strcmpi(tline(1:5),'TY  -')
            if nRec>0
                if length(rec.Journal)==nRec
                    nRec = nRec + 1;
                end
            else
                nRec = nRec + 1;
            end
            rec.email{nRec} = '';
            rec.Title{nRec} = '';
            rec.nAuthors(nRec) = 0;
            rec.AuthorList{nRec} = '';
            rec.isCollaborative(nRec) = 0;
            rec.CollabAuthors{nRec} = [];
            uniqueAuthorFlag = 0;
            uniqueAuthorListRec = [];  % track the unique authors for a given publication
        end
        
        % Year
        if strcmpi(tline(1:5),'Y1  -') | strcmpi(tline(1:5),'PY  -')
            if length(tline)<10
                foo = [];
            else
                foo = str2num(tline(7:10));
            end
            if ~isempty(foo)
                rec.Year(nRec) = foo;
                yearLast = foo;
            else
                rec.Year(nRec) = yearLast;
                nYearNotDefined = nYearNotDefined + 1;
            end
        end   
        
        % Journal
        if strcmpi(tline(1:5),'JA  -')
            rec.Journal{nRec} = tline(7:end);
        elseif strcmpi(tline(1:5),'T2  -')
            rec.Journal{nRec} = tline(7:end);
        elseif strcmpi(tline(1:5),'J2  -')
            rec.Journal{nRec} = tline(7:end);
        elseif strcmpi(tline(1:5),'JO  -')
            rec.Journal{nRec} = tline(7:end);
        end   
        
        % Title
        if strcmpi(tline(1:5),'T1  -') | strcmpi(tline(1:5),'TI  -')
            rec.Title{nRec} = tline(7:end);
        end   
        
        % Address
        % Parse it into a known country
        % get email address
%         if strcmpi(tline(1:5),'AD  -')
%             rec.Country(nRec) = -1;
%             for ii=1:nCountry
%                 for jj=1:length(lstCountry{ii})
%                     if ~isempty(regexpi(tline(6:end),lstCountry{ii}{jj}))
%                         rec.Country(nRec) = ii;
%                         break
%                     end
%                 end
%             end
%             if rec.Country(nRec)==-1
%                 disp(tline)
%                 rec.Country(nRec)=0;
%             end
%             
%             ii=regexpi(tline,'@');
%             lst=regexpi(tline,' ');
%             if ~isempty(ii) & ~isempty(lst)
%                 jj=find(lst<ii(1),1,'last');
%                 if jj==length(lst)
%                     rec.email{nRec} = tline((lst(jj)+1):end);
%                 else
%                     rec.email{nRec} = tline((lst(jj)+1):lst(jj+1)-1);
%                 end
%             end
%         end
        
        % Keywords from Abstract and Title
%         if strcmpi(tline(1:5),'N2  -') | strcmpi(tline(1:5),'AB  -')
%             nKey = 0;
%             rec.Keywords(nRec,1) = 0;
%             for ii=1:nKeywords
%                 for jj=1:length(lstKeywords{ii})
%                     foos = [rec.Title{nRec} ' ' tline(6:end)];
%                     if ~isempty(regexpi(foos,lstKeywords{ii}{jj}))
%                         nKey = nKey + 1;
%                         rec.Keywords(nRec,nKey) = ii;
%                         break
%                     end
%                 end
%             end
%         end
            
        % Author
        if (strcmpi(tline(1:5),'A1  -') | strcmpi(tline(1:5),'AU  -'))
            rec.nAuthors(nRec) = rec.nAuthors(nRec) + 1;
            ii = regexpi(tline(7:end),',');
            foos = tline(6+[1:ii+2]);
            
            rec.AuthorList{nRec} = sprintf('%s, %s', rec.AuthorList{nRec}, replace(foos,',','') );
            
            % if the publication author matches an NPC faculty, then flag it
            flag = 0;
            for ii=1:length(uniqueAuthors)
                if strcmpi(foos,uniqueAuthors{ii})  
                    flag = 1;
                    break
                end
            end
            % if the publication journal is on the exclude list, then flag
            flagInclude = 1;
            for jj=1:length(excludeJournals)
                if strcmpi(rec0.Journal{nRec},excludeJournals{jj})
                    flagInclude = 0;
                end
            end
            if yearTarget > 0 % then look for collaborative publications up to and including target year
                if flag==1 & rec0.Year(nRec)<=yearTarget & flagInclude  % if an NPC faculty and published <= targetYear then check if it is a collaborative paper
                    AuthorNumPubs(ii) = AuthorNumPubs(ii) + 1; % one more pub for this author

                    uniqueAuthorListRec(end+1) = ii; % keep track of unique NPC authors for the given publication

                    rec.CollabAuthors{nRec}(end+1) = ii;

                    if length(uniqueAuthorListRec)>1 % if we now have 2 or more unique NPC authors on a given pub, then add to the Adj matrix
                        AuthorNumCollabPubs(ii) = AuthorNumCollabPubs(ii) + 1; % one more collaborative pub for this author
                        if length(uniqueAuthorListRec)==2 % also one more collab pub for first author
                            AuthorNumCollabPubs(uniqueAuthorListRec(1)) = AuthorNumCollabPubs(uniqueAuthorListRec(1)) + 1; % 
                        end
                        for jj=1:(length(uniqueAuthorListRec)-1)
                            AuthorAdj(ii,uniqueAuthorListRec(jj)) = AuthorAdj(ii,uniqueAuthorListRec(jj)) + 1;
                            AuthorAdj(uniqueAuthorListRec(jj),ii) = AuthorAdj(uniqueAuthorListRec(jj),ii) + 1;
                        end

                        rec.isCollaborative(nRec) = 1; % mark the publication as collaborative
                    end                       
                end
            else % look for collaborative publications only in the target year
                % same code as above but ==yearTarget
                if flag==1 & rec0.Year(nRec)==-yearTarget & flagInclude  % if an NPC faculty and published <= targetYear then check if it is a collaborative paper
                    AuthorNumPubs(ii) = AuthorNumPubs(ii) + 1; % one more pub for this author

                    uniqueAuthorListRec(end+1) = ii; % keep track of unique NPC authors for the given publication

                    rec.CollabAuthors{nRec}(end+1) = ii;

                    if length(uniqueAuthorListRec)>1 % if we now have 2 or more unique NPC authors on a given pub, then add to the Adj matrix
                        AuthorNumCollabPubs(ii) = AuthorNumCollabPubs(ii) + 1; % one more collaborative pub for this author
                        if length(uniqueAuthorListRec)==2 % also one more collab pub for first author
                            AuthorNumCollabPubs(uniqueAuthorListRec(1)) = AuthorNumCollabPubs(uniqueAuthorListRec(1)) + 1; % 
                        end
                        for jj=1:(length(uniqueAuthorListRec)-1)
                            AuthorAdj(ii,uniqueAuthorListRec(jj)) = AuthorAdj(ii,uniqueAuthorListRec(jj)) + 1;
                            AuthorAdj(uniqueAuthorListRec(jj),ii) = AuthorAdj(uniqueAuthorListRec(jj),ii) + 1;
                        end

                        rec.isCollaborative(nRec) = 1; % mark the publication as collaborative
                    end                       
                end
            end

        end
            
    end
    tline = fgetl( fid );
end
fclose(fid);

disp( sprintf(' Number of records with non-defined year = %d',nYearNotDefined) )



    