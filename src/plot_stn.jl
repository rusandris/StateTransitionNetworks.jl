function plot_stn(stn;filename="stn.pdf",nodesize=1,nodefillc="orange",linetype="straight",max_edgelinewidth=1,nodelabels=false)
	nr_vertices = nv(stn)
	x = []
	y = []
		
	for i in 1:nr_vertices
		ilabel = label_for(stn,i)
		push!(x,stn[ilabel][:x])
		push!(y,stn[ilabel][:y])
	end
	x = Int32.(x)
	y = Int32.(y)
	
	w = weight_matrix(stn)
	
	gp = gplot(stn,x,-y,
		edgelinewidth = reduce(vcat, [w[i,:].nzval for i in 1:nv(stn)]) .^ 0.33,
		nodelabel = if nodelabels 1:nv(stn) else nothing end,
		nodesize=nodesize,
		nodefillc=nodefillc,
		linetype=linetype,
		EDGELINEWIDTH=max_edgelinewidth)
	
	draw(PDF(filename),gp)

end
