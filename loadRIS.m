%
% ( near infrared | ( NIRS | fnirs ) ) ( brain dp:2010 )
% try modifying this to reject (NOT) photo acoustic, fluorescence
%
% autoradiography and "blood flow"
% autoradiography and glucose
% fmri and brain
% resting state fmri
% ((NIRS or "near infrared") and brain) NOT (fluorescence OR "photoacoustic")
%
%
%
function rec = loadRIS( filenm, flagAuthor )


nRec = 0;
if ~exist('flagAuthor','var')
    flagAuthor = 0;
end
    
%%%%%%%%%%%
% COUNTRIES
lstCountry = {{'Japan','Tokyo','Tokai','Jichi','Hitachi','Osaka','Sendai'},...
    {'China','Beijing','Wuhan'},...
    {'Singapore'},...
    {'Taiwan'},...
    {'Korea'},...
    {'U.K.','United Kingdom','UK','England','Ireland'},...
    {'Germany','Berlin','Munich','wuerzburg','Mannheim'},...
    {'Italy'},...
    {'Switzerland'},...
    {'France'},...
    {'Spain'},...
    {'Netherlands'},...
    {'Austria'},...
    {'Finland'},...
    {'Poland'},...
    {'Sweden','Stockholm'},...
    {'Norway'},...
    {'Belgium'},...
    {'Slovenia'},...
    {'Greece'},...
    {'Hungary'},...
    {'Turkey','Istanbul','Bogazici'},...
    {'Denmark'},...
    {'USA','United States','New York','Illinois','Pennsylvania','Virginia','Tufts','Miami'},...
    {'Canada','British Columbia'},...
    {'Brazil'},...
    {'Israel'},...
    {'New Zealand'},...
    {'Australia'},...
    {'South Africa'},...
    };
nCountry = length(lstCountry);
lstContinent = {'Asia','Europe','America','Other'};
lstContinentIdx = {[1 2 3 4 5],[6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23],[24 25 26],[27 28 29 30]};


%%%%%%%%%%%
% KEYWORDS
lstKeywords = {...
    {'Motion Artifact','Motion correction'},...
    {'Functional Connectivity','Resting State'},...
    {'Mobility','gait','wireless','walk','wearable','portable'},...
    {'Prefrontal Cortex','PFC'},...
    {'Infants','development','baby'},...
    {'Motor Sensory Cortex','motor','somatosensory','finger','median nerve'},...
    {'Visual Cortex','visual stimulation','occipital'},...
    {'Psychiatry','Psychiatric','depression','depressive','abuse','schizophrenia','addiction'},...
    {'Verbal Fluency'},...
    {'Pain'},...
    {'Epilepsy','seizure','ictal'},...
    {'TMS/DCS','TMS','transcranial magnetic stimulation','direct current stimulation'},...
    {'BCI','brain computer interface','brain machine interface'},...
    {'NVC','neuro vascular coupling','neuro-vascular coupling','neurovascular coupling'},...
    {'Stroke','Ischemia'},...
    {'Cross-talk','cross talk'},...
    {'fMRI'},...
    {'Physiologic Interference','Contamination'},...
    {'Oscillations','reactivity','autoregulation'},...
    {'Cognitive'},...
    {'Alzheimer','Alzheimer''s'},...
    {'EEG'},...
    {'Schizophrenia'},...
    {'Depression','depressive','abuse','schizophrenia','addiction'},...
    {'Hyperscanning','hyper-scanning'},...
    {'Anatomical guidance','atlas'}...
    {'short separation'}...
    };
nKeywords = length(lstKeywords);

rec.uniqueAuthors = {};
rec.uniqueAuthorsNum = [];
rec.uniqueAuthorsByYear = zeros(2025,1);
rec.filenm = filenm;

fid = fopen( filenm, 'r');
tline = fgetl( fid );
yearLast = [];
nYearNotDefined = 0;
while ischar(tline)
    if length(tline)>5
        % Start of a new record
        if strcmpi(tline(1:5),'TY  -')
            if nRec>0
                if length(rec.Journal)==nRec % if the prior record did not have a 
                    % journal then this results in the new record over-writing the 
                    % previous one that did not have a journal
                    nRec = nRec + 1;
                end
            else
                nRec = nRec + 1;
            end
            rec.email{nRec} = '';
            rec.Title{nRec} = '';
            rec.nAuthors(nRec) = 0;
            uniqueAuthorFlag = 0;
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
            if uniqueAuthorFlag>0
                rec.uniqueAuthorsByYear(rec.Year(nRec)) = rec.uniqueAuthorsByYear(rec.Year(nRec)) + uniqueAuthorFlag;
            end
        end   
        
        % Journal
        if strcmpi(tline(1:5),'JA  -')
            rec.Journal{nRec} = tline(7:end);
        elseif strcmpi(tline(1:5),'J2  -')
            rec.Journal{nRec} = tline(7:end);
        elseif strcmpi(tline(1:5),'JO  -')
            rec.Journal{nRec} = tline(7:end);
        elseif strcmpi(tline(1:5),'T2  -')
            rec.Journal{nRec} = tline(7:end);
        end   
        
        % Title
        if strcmpi(tline(1:5),'T1  -') | strcmpi(tline(1:5),'TI  -')
            rec.Title{nRec} = tline(7:end);
        end   
        
        % Address
        % Parse it into a known country
        % get email address
        if strcmpi(tline(1:5),'AD  -')
            rec.Country(nRec) = -1;
            for ii=1:nCountry
                for jj=1:length(lstCountry{ii})
                    if ~isempty(regexpi(tline(6:end),lstCountry{ii}{jj}))
                        rec.Country(nRec) = ii;
                        break
                    end
                end
            end
            if rec.Country(nRec)==-1
                disp(tline)
                rec.Country(nRec)=0;
            end
            
            ii=regexpi(tline,'@');
            lst=regexpi(tline,' ');
            if ~isempty(ii) & ~isempty(lst)
                jj=find(lst<ii(1),1,'last');
                if jj==length(lst)
                    rec.email{nRec} = tline((lst(jj)+1):end);
                else
                    rec.email{nRec} = tline((lst(jj)+1):lst(jj+1)-1);
                end
            end
        end
        
        % Keywords from Abstract and Title
        if strcmpi(tline(1:5),'N2  -') | strcmpi(tline(1:5),'AB  -')
            nKey = 0;
            rec.Keywords(nRec,1) = 0;
            for ii=1:nKeywords
                for jj=1:length(lstKeywords{ii})
                    foos = [rec.Title{nRec} ' ' tline(6:end)];
                    if ~isempty(regexpi(foos,lstKeywords{ii}{jj}))
                        nKey = nKey + 1;
                        rec.Keywords(nRec,nKey) = ii;
                        break
                    end
                end
            end
        end
            
        % Author
        if (strcmpi(tline(1:5),'A1  -') | strcmpi(tline(1:5),'AU  -')) & flagAuthor==1
            rec.nAuthors(nRec) = rec.nAuthors(nRec) + 1;
            ii = regexpi(tline(7:end),',');
            foos = tline(6+[1:ii+2]);
            flag = 0;
            for ii=1:length(rec.uniqueAuthors)
                if strcmpi(foos,rec.uniqueAuthors{ii})
                    flag = 1;
                    break
                end
            end
            if flag==1
                rec.uniqueAuthorsNum(ii) = rec.uniqueAuthorsNum(ii) + 1;
            else
                rec.uniqueAuthors{end+1} = foos;
                rec.uniqueAuthorsNum(end+1) = 1;
                uniqueAuthorFlag = uniqueAuthorFlag + 1;
            end
        end
            
    end
    tline = fgetl( fid );
end
fclose(fid);

disp( sprintf(' Number of records with non-defined year = %d',nYearNotDefined) )

%%%%%%%%%%%%%%%%
% Some Figures
figure(1);
x=[1993:2017];
hist(rec.Year,x)
y=histc(rec.Year,x);
x1 = x(1:end-1);
y1 = y(1:end-1);
p=polyfit(x1,log(y1),1);
figure(4)
plot(x1,y1,'.',x1,exp(p(1)*x1+p(2)),'-')
p=polyfit(x(1:10),log(y(1:10)),1);
figure(6)
plot(x(1:10),y(1:10),'.',x(1:10),exp(p(1)*x(1:10)+p(2)),'-')

if 0
    figure(2)
    countryLabel{1} = 'Unk';
    for ii=1:nCountry
        countryLabel{ii+1} = lstCountry{ii}{1}(1:3);
    end
    hist(rec.Country,[0:1:length(countryLabel)-1])
    set(gca,'xtick',[0:length(countryLabel)-1])
    set(gca,'xticklabel',countryLabel)
    
    figure(3)
    lst = find(rec.Year<=2002);
    y=[];
    y(:,1) = histc(rec.Country(lst),[0:length(countryLabel)-1]);
    lst = find(rec.Year>2002);
    y(:,2) = histc(rec.Country(lst),[0:length(countryLabel)-1]);
    bar(y)
    set(gca,'xtick',[1:length(countryLabel)])
    set(gca,'xticklabel',countryLabel)
    legend('<=2002','>2002')
    
    figure(5)
    y=histc(rec.Country,[0:1:length(countryLabel)-1]);
    yy=[];
    for ii=1:length(lstContinent)
        yy(ii) = 0;
        for jj=1:length(lstContinentIdx{ii})
            yy(ii) = yy(ii) + y(lstContinentIdx{ii}(jj));
        end
    end
    yy(end) = yy(end) + y(1);
    bar(yy)
    set(gca,'xtick',[1:length(lstContinent)])
    set(gca,'xticklabel',lstContinent)
    
end


%%%%%%%%%%%
% Journals
uniqueJournals = unique(lower(rec.Journal));
for ii=1:length(uniqueJournals)
    nJournalPubs(ii) = length(find(strcmpi(rec.Journal,uniqueJournals{ii})==1));
end
[foo,idxOrder] = sort(nJournalPubs,2,'descend');
disp(' ')
disp(' ')
disp(' ')
disp('======================================')
disp('Number of Publications in Each Journal')
disp(sprintf('   %d unique journals',length(uniqueJournals)))
disp('======================================')
for ii=1:length(uniqueJournals)
    disp(sprintf('%2d) %3d  -  %s',ii,nJournalPubs(idxOrder(ii)),uniqueJournals{idxOrder(ii)}))
end


%%%%%%%%%%%
% Keywords
for ii=1:nKeywords
    [ir,ic] = find(rec.Keywords==ii);
    nKeywordArticles(ii) = length(ir);
    jj = 0;
    dy = 1;
    for iYear = 1993:dy:2016
        jj = jj + 1;
        lst = find(rec.Year(ir)>=iYear & rec.Year(ir)<(iYear+dy));
        nKeywordPerYearGroup(ii,jj) = length(lst);
    end
end
[foo,idxOrder] = sort(nKeywordArticles,2,'descend');
disp(' ')
disp(' ')
disp(' ')
disp('========================================')
disp('Number of Publications with Each Keyword')
disp('========================================')
for ii=1:nKeywords
    foos = '';
%    for jj=1:size(nKeywordPerYearGroup,2)
    for jj=(size(nKeywordPerYearGroup,2)-4):size(nKeywordPerYearGroup,2)
        foos = sprintf('%s,%d',foos,nKeywordPerYearGroup(idxOrder(ii),jj));
    end
    disp(sprintf('%2d) %3d  -  %s    (%s)',ii,nKeywordArticles(idxOrder(ii)),lstKeywords{idxOrder(ii)}{1},foos(2:end)))
end
disp(sprintf('    %3d  -  None',length(find(rec.Keywords(:,1)==0))))



%%%%%%%%%%%
% Authors
if flagAuthor==1
    disp(' ')
    disp(' ')
    disp(' ')
    disp('======================================')
    disp('Number of Publications for Each Author')
    disp(sprintf('    %d unique authors',length(rec.uniqueAuthors)))
    disp('======================================')
    [foo,idxOrder] = sort(rec.uniqueAuthorsNum,2,'descend');
    for ii=1:200
        num = rec.uniqueAuthorsNum(idxOrder(ii));
        foos = rec.uniqueAuthors{idxOrder(ii)};
        disp(sprintf('%2d) %3d  -  %s',ii,num,foos))
    end
    
    for ii=1:200
        topNames{ii} = rec.uniqueAuthors{idxOrder(ii)};
    end
    
end
    