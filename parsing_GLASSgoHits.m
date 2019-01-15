
close all;
clear;clc
addpath('D:\matlab_myfunction');

cd D:\Borellia_sRNA\GLASSgo\novel_miRNA_NCI60

%%
fname='novel_miRNA_NCI60.fasta';
fa=readtext(fname,'\t');
miRNA=[fa(1:2:end),fa(2:2:end)];
miRNA(:,1)=regexprep(miRNA(:,1),'>','');
len_miRNA=strlength(miRNA(:,2));

n=length(miRNA);
fieldnames={'ID','length_miRNA','strain','chr','strand','start-end','length_hits','p.c.VAL','taxID'};
hits=fieldnames;
for i=1:n
    if mod(i,10)==0
        disp([i,n]);
    end
    
    fname=strcat('novel_miRNA_NCI60_',num2str(i));
    text = fileread(fname);
    TextAsCells = regexp(text, '\n', 'split');
    
    if isempty(TextAsCells{1})
        disp(strcat('No hits:',num2str(i)));
    end
    
    if ~isempty(TextAsCells{1})
        TextAsCells(1:2)=[];
        ind=strfindidx('>',TextAsCells);        
        
        hits_i=cell(length(ind),length(fieldnames));
        for j=1:length(ind)
            str=TextAsCells{ind(j)};
            
            strain=regexp(str,'(\s[A-Za-z]+\s[-_:;+,>/\''\\\[\]\(\)\.\w\s]*p.c.VAL:)','match');
            strain=regexprep(strain,'-p.c.VAL:','');
            
            chr=regexp(str, '([-\.\w]+:[c\d]*-\d+)','match');
            chr=cellsplit(chr,':','s');
            
            strand={'+'};
            if strfind(chr{2},'c')
                strand={'-'};
                chr(2)=regexprep(chr(2),'c','');
            end
            
            se=cellsplit(chr(2),'-','n');
            len_hit=num2cell(abs(se(1)-se(2))+1);
            
            pident=regexp(str, '(p.c.VAL:\d+.\d+%)','match');
            pident=cellsplit(pident,':','s');
            pident=regexprep(pident(2),'%','');
            
            taxID=regexp(str,'(taxID:\d+)','match');
            taxID=cellsplit(taxID,':','s');
            taxID=taxID(2);
            
            hits_i(j,:)=[miRNA(i,1),len_miRNA(i),strain,chr(1),strand,chr(2),len_hit,pident,taxID];
        end

        hits=[hits;hits_i];
    end
end

cell2txt('GLASSgo_novel_miRNA_NCI60.txt',hits,'\t');


















