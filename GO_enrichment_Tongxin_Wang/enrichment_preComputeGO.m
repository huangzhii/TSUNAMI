clc
clear

% precompute the getrelatives function for all the genes in the platform
%load('GOannotation.mat');
folder = 'data';
geoID = 'GSE31684'; %GSE33382 %GSE39058 %GSE18842
GPLFile = 'GPL570-55999.txt'; %GPL10295.txt %GPL14951-11332.txt %GPL570-55999.txt
[probeID, geneSymbol, geneBankAcc] = GPLRead_570(fullfile(folder,geoID,GPLFile));

%%%%%%%%%%
uniGene = unique(geneSymbol);

GO = geneont('live',true);
totalGO = GO.Terms(end).id;
HGann = goannotread('goa_human.gaf','Fields',{'DB_Object_Symbol','GOid'});
HGmap = containers.Map();

% seperate 3 parts of GO
mf_mask = strcmp(get(GO.terms,'ontology'),'molecular function');
bp_mask = strcmp(get(GO.terms,'ontology'),'biological process');
cc_mask = strcmp(get(GO.terms,'ontology'),'cellular component');

mf_GO = GO(GO.terms(mf_mask));
bp_GO = GO(GO.terms(bp_mask));
cc_GO = GO(GO.terms(cc_mask));

for i = 1:numel(HGann)
    key = HGann(i).DB_Object_Symbol;
%     if ( isempty(mf_GO(HGann(i).GOid).Term) && isempty(bp_GO(HGann(i).GOid).Term) ...
%         && isempty(cc_GO(HGann(i).GOid).Term))
%         continue;
%     end
    if (isKey(HGmap,key))
        HGmap(key) = [HGmap(key) HGann(i).GOid];
    else
        HGmap(key) = HGann(i).GOid;
    end
end

% take 1 level up and 0 level down when calculating relatives
% use map to represent relatives
mf_GOrelatives = containers.Map();
bp_GOrelatives = containers.Map();
cc_GOrelatives = containers.Map();
fprintf('Calculating GO relatives...\n');
for i = 1:length(uniGene)
    fprintf('%d/%d\n',i,length(uniGene));
    if isKey(HGmap,uniGene{i})
        mf_GOrelatives(uniGene{i}) = getrelatives(mf_GO,unique(HGmap(uniGene{i})),'height',1,'depth',0);
        bp_GOrelatives(uniGene{i}) = getrelatives(bp_GO,unique(HGmap(uniGene{i})),'height',1,'depth',0);
        cc_GOrelatives(uniGene{i}) = getrelatives(cc_GO,unique(HGmap(uniGene{i})),'height',1,'depth',0);
    end
end

%precalculate total gene count
fprintf('Calculating totalGenesCount...\n');
mf_totalGenesCount = zeros(totalGO,1); % a vector of GO term counts for the entire chip.
bp_totalGenesCount = zeros(totalGO,1);
cc_totalGenesCount = zeros(totalGO,1);

for i = 1:length(uniGene)
    fprintf('%d/%d\n',i,length(uniGene));
    if isKey(mf_GOrelatives,uniGene{i})
        mf_goid = mf_GOrelatives(uniGene{i});
        mf_totalGenesCount(mf_goid) = mf_totalGenesCount(mf_goid) + 1;
    end
    if isKey(bp_GOrelatives,uniGene{i})
        bp_goid = bp_GOrelatives(uniGene{i});
        bp_totalGenesCount(bp_goid) = bp_totalGenesCount(bp_goid) + 1;
    end
    if isKey(cc_GOrelatives,uniGene{i})
        cc_goid = cc_GOrelatives(uniGene{i});
        cc_totalGenesCount(cc_goid) = cc_totalGenesCount(cc_goid) + 1;
    end
end

save(fullfile(folder,geoID,[GPLFile(1:end-4),'_GOrelatives.mat']), ...
    'mf_GOrelatives','bp_GOrelatives','cc_GOrelatives',...
    'mf_totalGenesCount', 'bp_totalGenesCount', 'cc_totalGenesCount',...
    '-v7.3');

