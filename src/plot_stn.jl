function plot_stn(stn;filename="stn.pdf")
	nr_vertices = nv(stn)
	x = zeros(Int32,nr_vertices)
	y = zeros(Int32,nr_vertices)
		
	for i in 1:nr_vertices
           x[i] = get_prop(stn,i,:pos)[1]
           y[i] = get_prop(stn,i,:pos)[2]
	end

	w = [get_prop(stn,edge,:weight) for edge in collect(edges(stn))]
	
	gp = gplot(stn,x,y,edgelinewidth = w .* 5, EDGELINEWIDTH=1)
	draw(PDF(filename),gp)

end
