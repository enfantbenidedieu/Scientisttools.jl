"""
    geom_circle(gg,x,y;...)

Draw circle

...
# Description

# Author(s)
Duvérier DJIFACK ZEBAZE djifacklab@gmail.com

...
"""
function geom_circle(gg,x,y;center=(0,0),diameter=2,npoints=100,color="black")
    r = diameter / 2
    tt = LinRange(0,2*pi,npoints)
    xx, yy = center[1] .+ r * cos.(tt), center[2] .+ r * sin.(tt)
    dat = DataFrame(Dict(x => xx, y => yy))
    return gg.geom_path(data=dat,gg.aes(x=x,y=y),color=color)
end;
"""
    geom_segment(gg,data;...)

Draw segment

...
# Description 


# Author(s)
Duvérier DJIFACK ZEBAZE djifacklab@gmail.com

...
"""
function geom_segment(gg,data;color="black")
    function add_zeros(x,y;x_lab="Dim.1",y_lab="Dim.2")
        return DataFrame(Dict(x_lab => [0,x], y_lab => [0,y]))    
    end;
    x_dim,y_dim = names(data)
    dat = reduce(vcat,[add_zeros(data[i,1],data[i,2],x_lab=x_dim,y_lab=y_dim) for i in 1:nrow(data)]);
    return gg.geom_line(data=dat,gg.aes(x = x_dim, y= y_dim),color=color)
end;