function GenerateNonLinearSystemFunction(function_filename,system,variable_string)

% This function generates the actual code of the function m-file that is 
% required in order to numerically solve the given system of nonlinear
% equations.

filename = strcat([function_filename '.m']);
heading_linestring = strcat(['function F = ' function_filename '(' variable_string ')\n\n']);
fid = fopen(filename,'w');
fprintf(fid,heading_linestring);
eqs_num = length(system);
linestring = strcat(['F = zeros(' num2str(eqs_num) ');\n']);
fprintf(fid,linestring);
for equation_index = 1:1:eqs_num
    equation = system(equation_index);
    linestring = char(equation);
    for variables_num = 1:1:eqs_num
        substr = strcat([variable_string num2str(variables_num)]);
        newstr = strcat([variable_string '(' num2str(variables_num) ')']);
        linestring = strrep(linestring,substr,newstr);
    end;
    linestring = strcat(['F(' num2str(equation_index) ') = ' linestring ';\n']);
    linestring = strrep(linestring,'==','-');
    fprintf(fid,linestring);
end;    
fclose(fid);

end

