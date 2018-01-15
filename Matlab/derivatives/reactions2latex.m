run models/sbml/lacy_lacz2.m


fp = fopen('table.tex','w');

fprintf(fp,'\\begin{center} \n');
fprintf(fp,'\\begin{tabular}{ c | c | c } \n');
fprintf(fp,'\\hline \n');

M=length(reaction);
for i=1:M
    
    r = strrep(reaction{i},'->','$\longrightarrow$');
    
%     r = [ '\ce' reaction{i} '}' ];
    
    fprintf(fp,' $R_{%d}$  &  %s  &  %f  %s  \\hline \n',i,r,rate(i),'\\');
    
    
end


fprintf(fp,'\\end{tabular} \n');
fprintf(fp,'\\end{center} \n');

  

fprintf(fp,'\n \n \n \n');




fprintf(fp,'\\begin{center} \n');
fprintf(fp,'\\begin{tabular}{ c | c } \n');
fprintf(fp,'\\hline \n');

N=length(species);
for i=1:M
    
    fprintf(fp,'  %s  &  %f  %s  \\hline \n',species{i},rate(i),'\\');
    
end


fprintf(fp,'\\end{tabular} \n');
fprintf(fp,'\\end{center} \n');








fclose(fp);
