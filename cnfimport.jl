using SparseArrays

#clause_length = 3


colInfo = Int64[]
rowInfo = Int64[]
dataInfo = Int64[]
index= 0



file = open("cnfProblems//uf20-91//uf20-01.cnf")
for line in eachline(file)
    if startswith(line,"c") || isempty(line)
        continue
    elseif startswith(line,"p")
        global problem_info = line
    else
        global index +=1
        junk = split(line)
        for i in junk
            junk2 = tryparse(Int,i)
            if junk2 != 0 && !isnothing(junk2)
                push!(rowInfo, index)
                push!(colInfo,abs(junk2))
                push!(dataInfo,-sign(junk2))
            end
        end
    end
    #colInfo = zeros(parse(Int,problemType[4])*clause_length)
end
close(file)

A = transpose(sparse(colInfo,rowInfo,dataInfo))
