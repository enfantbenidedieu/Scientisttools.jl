# Set directory
basename = joinpath(@__DIR__,"..","data")

"""
    load_decathlon(; choice="all")


"""
function load_decathlon(; choice="all")::DataFrame
    # Check if choice in vector
    if !(choice in ["actif","ind_sup","quanti_sup","quali_sup","all"]) throw(ArgumentError("'choice' should be one of 'actif', 'ind_sup', 'quanti_sup', 'quali_sup', 'all'.")) end;
    
    # Import actif elements
    if choice in ["actif","ind_sup","quanti_sup","quali_sup"]
        if choice === "actif"
            data,columns = XLSX.readtable(joinpath(string(basename,"/decathlon.xlsx")),"Feuil1");
        elseif choice === "ind_sup"
            data,columns = XLSX.readtable(joinpath(string(basename,"/decathlon.xlsx")),"Feuil2");
        elseif choice === "quanti_sup"
            data,columns = XLSX.readtable(joinpath(string(basename,"/decathlon.xlsx")),"Feuil3");
        else
            data,columns = XLSX.readtable(joinpath(string(basename,"/decathlon.xlsx")),"Feuil4");
        end;
        donnee = DataFrame(data,columns)
    else
        data, columns = XLSX.readtable(joinpath(string(basename,"/decathlon.xlsx")),"Feuil1")
        data2, columns2 = XLSX.readtable(joinpath(string(basename,"/decathlon.xlsx")),"Feuil2")
        data3, columns3 = XLSX.readtable(joinpath(string(basename,"/decathlon.xlsx")),"Feuil3")
        data4, columns4 = XLSX.readtable(joinpath(string(basename,"/decathlon.xlsx")),"Feuil4")
        donnee, donnee2 = DataFrame(data,columns), DataFrame(data2,columns2)
        donnee3, donnee4 = DataFrame(data3,columns3), DataFrame(data4,columns4)
        donnee = leftjoin(leftjoin(vcat(donnee,donnee2),donnee3,on=:Athlètes),donnee4,on=:Athlètes)
    end;
    return donnee
end;

"""
    load_children(; choice="all")


...
# Description


# Examples


# Author(s)

...
"""
function load_children(; choice="all")
    # Check if choice in vector
    if !(choice in ["actif","row_sup","col_sup","quali_sup","all"]) throw(ArgumentError("'choice' should be one of 'actif', 'row_sup', 'col_sup', 'quali_sup', 'all'.")) end;
    
    # Import actif elements
    if choice in ["actif","row_sup","col_sup","quali_sup"]
        if choice === "actif"
            data, columns = XLSX.readtable(joinpath(string(basename,"/children.xlsx")),"Feuil1");
        elseif choice === "row_sup"
            data, columns = XLSX.readtable(joinpath(string(basename,"/children.xlsx")),"Feuil2");
        elseif choice === "col_sup"
            data, columns = XLSX.readtable(joinpath(string(basename,"/children.xlsx")),"Feuil3");
        else
            data, columns = XLSX.readtable(joinpath(string(basename,"/children.xlsx")),"Feuil4");
        end;
        donnee = DataFrame(data,columns)
    else
        data,columns = XLSX.readtable(joinpath(string(basename,"/children.xlsx")),"Feuil1")
        data2,columns2 = XLSX.readtable(joinpath(string(basename,"/children.xlsx")),"Feuil2")
        data3,columns3 = XLSX.readtable(joinpath(string(basename,"/children.xlsx")),"Feuil3")
        data4,columns4 = XLSX.readtable(joinpath(string(basename,"/children.xlsx")),"Feuil4")
        donnee,donnee2 = DataFrame(data,columns),DataFrame(data2,columns2)
        donnee3,donnee4 = DataFrame(data3,columns3),DataFrame(data4,columns4)
        donnee = leftjoin(leftjoin(vcat(donnee,donnee2),donnee3,on=:rownames),donnee4,on=:rownames)
    end;
end;

"""
    load_races_canines(; choice="all")


...

...
"""
function load_races_canines(; choice="all")
    # Check if choice in vector
    if !(choice in ["actif","ind_sup","quanti_sup","quali_sup","all"]) throw(ArgumentError("'choice' should be one of 'actif', 'ind_sup', 'quali_sup', 'quanti_sup', 'all'.")) end;
    
    if choice in ["actif","ind_sup","quanti_sup","quali_sup"]
        if choice === "actif"
            data,columns = XLSX.readtable(joinpath(string(basename,"/races_canines.xlsx")),"Feuil1");
        elseif choice === "ind_sup"
            data,columns = XLSX.readtable(joinpath(string(basename,"/races_canines.xlsx")),"Feuil2");
        elseif choice === "quali_sup"
            data,columns = XLSX.readtable(joinpath(string(basename,"/races_canines.xlsx")),"Feuil3");
        else
            data,columns = XLSX.readtable(joinpath(string(basename,"/races_canines.xlsx")),"Feuil4");
        end;
        donnee = DataFrame(data,columns)
    else
        data, columns = XLSX.readtable(joinpath(string(basename,"/races_canines.xlsx")),"Feuil1")
        data2, columns2 = XLSX.readtable(joinpath(string(basename,"/races_canines.xlsx")),"Feuil2")
        data3, columns3 = XLSX.readtable(joinpath(string(basename,"/races_canines.xlsx")),"Feuil3")
        data4, columns4 = XLSX.readtable(joinpath(string(basename,"/races_canines.xlsx")),"Feuil4")
        donnee, donnee2 = DataFrame(data,columns), DataFrame(data2,columns2)
        donnee3, donnee4 = DataFrame(data3,columns3), DataFrame(data4,columns4)
        donnee = leftjoin(leftjoin(vcat(donnee,donnee2),donnee3,on=:Chien),donnee4,on=:Chien)
    end;
    return donnee
end;
"""



"""
function load_autos2005(;choice="all")
    # Check if choice in vector
    if !(choice in ["actif","ind_sup","quanti_sup","quali_sup","all"]) throw(ArgumentError("'choice' should be one of 'actif', 'ind_sup', 'quanti_sup', 'quali_sup', 'all'.")) end;
    
    # Import actif elements
    if choice in ["actif","ind_sup","quanti_sup","quali_sup"]
        if choice === "actif"
            data,columns = XLSX.readtable(joinpath(string(basename,"/autos2005.xlsx")),"Feuil1");
        elseif choice === "ind_sup"
            data,columns = XLSX.readtable(joinpath(string(basename,"/autos2005.xlsx")),"Feuil2");
        elseif choice === "quanti_sup"
            data,columns = XLSX.readtable(joinpath(string(basename,"/autos2005.xlsx")),"Feuil3");
        else
            data,columns = XLSX.readtable(joinpath(string(basename,"/autos2005.xlsx")),"Feuil4");
        end;
        donnee = DataFrame(data,columns)
    else
        data, columns = XLSX.readtable(joinpath(string(basename,"/autos2005.xlsx")),"Feuil1")
        data2, columns2 = XLSX.readtable(joinpath(string(basename,"/autos2005.xlsx")),"Feuil2")
        data3, columns3 = XLSX.readtable(joinpath(string(basename,"/autos2005.xlsx")),"Feuil3")
        data4, columns4 = XLSX.readtable(joinpath(string(basename,"/autos2005.xlsx")),"Feuil4")
        donnee, donnee2 = DataFrame(data,columns), DataFrame(data2,columns2)
        donnee3, donnee4 = DataFrame(data3,columns3), DataFrame(data4,columns4)
        donnee = leftjoin(leftjoin(vcat(donnee,donnee2),donnee3,on=:Modele),donnee4,on=:Modele)
    end;
    return donnee
end;
"""
    get_dataset(dataset_name;...)



"""
function get_dataset(dataset_name::AbstractString;choice::AbstractString="all")
    if dataset_name == "decathlon"
        return load_decathlon(choice=choice)
    elseif dataset_name == "children"
        return load_children(choice=choice)
    elseif dataset_name == "races_canines"
        return load_races_canines(choice=choice)
    elseif dataset_name == "autos2005"
        return load_autos2005(choice=choice)
    elseif dataset_name == "wine"
        data, columns = XLSX.readtable(joinpath(string(basename,"/wine.xlsx")),"Feuil1")
        return DataFrame(data,columns)
    end;
end;