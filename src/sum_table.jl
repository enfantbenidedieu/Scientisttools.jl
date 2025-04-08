# Load packages
using DataFrames

# Create sum table
function sum_table(X,Y,col)
    stats = combine(groupby(hcat(X,select(Y,col)),col),names(X) .=> sum,renamecols=false)
    #rename first columns
    rename!(stats,col => :Modalities)
    #sortby Modalities
    sort!(stats,:Modalities)
    return stats
end;
