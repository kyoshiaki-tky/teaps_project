function [tableFinal] = importStr(filename, split)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%ref
%http://questionbox.jp.msn.com/qa3181359.html
%filename='CLs.csv';
%split = ',';

fid=fopen(filename);
table0=textscan(fid,'%s', 'delimiter', '', 'MultipleDelimsAsOne',0);
fclose(fid);
table0=table0{:};
[sizeTableRow sizeTableColumn] = size(table0);
%{
  table1 = char(table0);
  table1 = [repmat(split, sizeTableRow,1) table1 repmat(split, sizeTableRow,1)];

  table2 = mat2cell(table2, ones(1,size(table2,1)));
  index = strfind(table2, split);
%}  
  
for k = 1: sizeTableRow
  table1 = char(table0(k));
  table1 = [split table1 split];
  index = strfind(table1, split); %search split
  %index = findstr(table1, split); %search split
  
  splitline = cell((length(index) - 1), 1);
  for count = 1:(length(index) - 1)
    data = table1((index(count) + 1):(index(count + 1) - 1));

    if isempty(data)
      splitline(count)=cellstr('0');
    else
      splitline(count) = {data};
    end

  end
  tableFinal(k,:)=splitline';
end
  

end

