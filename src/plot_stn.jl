function plot_stn(stn;filename="stn.pdf",nodesize=1,nodefillc="orange",linetype="straight",max_edgelinewidth=1)
	nr_vertices = nv(stn)
	x = []
	y = []
		
	for i in 1:nr_vertices
		push!(x,stn[i][:x])
		push!(y,stn[i][:y])
	end
	x = Int32.(x)
	y = Int32.(y)
	
	w = weight_matrix(stn)
	
	gp = gplot(stn,x,-y,
		edgelinewidth = w.nzval .^ 1/4,
		nodesize=nodesize,
		nodefillc=nodefillc,
		linetype=linetype,
		EDGELINEWIDTH=max_edgelinewidth)
	
	draw(PDF(filename),gp)

end
