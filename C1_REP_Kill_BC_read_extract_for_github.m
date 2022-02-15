% This code extracts Barcode of length 11 bases from R1 (forward reads) and R2(reverse reads) 
% which is between V1 and V2, then filtered through InfusionBarcode inventory, 
% also above cutoff (10). 
%% Infusion Barcode
[~, ~, raw] = xlsread('Infusion_Barcode_20K.xlsx','Sheet1');
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
stringVectors = string(raw(:,1));
stringVectors(ismissing(stringVectors)) = '';
m1 = size(stringVectors);
for k=1:m1
    Sequ1= stringVectors(k);
    barcode(k) = extractBetween(Sequ1,"AAACGTCTTGCTCGAG","GTGGCGGCCGCTCTAGAAC");%R1
end
InfusionBarcode = barcode';
%%
clearvars -except InfusionBarcode
cutoff = 10; 
warning('off','all')
warning( 'off', 'MATLAB:xlswrite:AddSheet' ) ;
 
for i = 1:78 
    clearvars -except InfusionBarcode i cutoff
    tic
for j=1:2 % j = 1 for R1 and j = 2 for R2
    if j ==1
filename_fastq = '1_S1_L001_R1_001.fastq';
    elseif j==2
        filename_fastq = '1_S1_L001_R2_001.fastq';
    end
read_sequ= fastqread(filename_fastq);
temp = struct2cell(read_sequ.').'; 
[m1,~] = size(temp);
Sequ = temp(1:m1,2);
[tot_seq, ~]  = size(temp);
BC_extracted = strings([m1,1]);
BC_Backward = strings([m1,1]);
for k=1:m1
    Sequ1= string(Sequ(k));
    if j==1
        BC = extractBetween(Sequ1,"AAACGTCTTGCTCGAG","GTGGCGGCCGCTCTAGAAC");%R1
        TT = ~isempty(BC);
        BC_L = length(char(BC));
        if (TT == 0||BC_L~=11)
            BC = "";
        end
        BC_extracted(k) = BC;
    else
        BC = extractBetween(Sequ1,"GTTCTAGAGCGGCCGCCAC","CTCGAGCAAGACGTTT");%R2 
        TT = ~isempty(BC);
        BC_L = length(char(BC));
        if (TT == 0||BC_L~=11)
            BC = "";
        end
        BC_extracted(k) = BC;
    end
end
BC_extracted = BC_extracted(~cellfun(@isempty, BC_extracted));
[BC_unq,ia_BC,ic_BC] = unique(BC_extracted);
BC_counts = accumarray(ic_BC,1);
Barcode = [BC_unq, BC_counts];
if j==1 %of R1
    R1_BC = Barcode; 
    R1_unq_c = double(R1_BC(:,2));
    R1_total = sum(R1_unq_c);
    for k=1:length(R1_BC)
        R1 = R1_BC(k,1);
        R1 = convertStringsToChars(R1);
        R1_com1 = seqrcomplement(R1);
        R1_com1 = strcat(R1_com1);
        R1_RC(k,:)= cellstr(R1_com1);
    end
    R1_RC = [R1_RC,R1_BC(:,2)];
else %of R2
    R2_BC = Barcode;
    R2_unq_c = double(R2_BC(:,2));
    R2_total = sum(R2_unq_c);
    for k=1:length(R2_BC)
        R2 = R2_BC(k,1);
        R2 = convertStringsToChars(R2);
        R2_com2 = seqrcomplement(R2);
        R2_com2 = strcat(R2_com2);
        R2_RC(k,:)= cellstr(R2_com2);
    end
    R2_RC = [R2_RC,R2_BC(:,2)];
end
end
TF_R1BC_R2RC = ismember((R1_BC(:,1)),(R2_RC(:,1)));
R1_BC = R1_BC(TF_R1BC_R2RC,:);
TF_InfusionBarcode = ismember((R1_BC(:,1)),InfusionBarcode);
R1_BC = R1_BC(TF_InfusionBarcode,:);
idx = find(double(R1_BC(:,2))>= cutoff);
R1_BC = R1_BC(idx,:);
R1_BC = table(R1_BC(:,1),double(R1_BC(:,2)));
end

