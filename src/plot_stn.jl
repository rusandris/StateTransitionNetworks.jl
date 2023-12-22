export plot_stn
function plot_stn(stn;filename="stn.pdf",nodesize=1,nodefillc="orange",linetype="straight",max_edgelinewidth=1,weight_exponent=1,nodelabels=false,kwargs...)
	
	prob_states,pos_states = get_state_distribution(stn)
	#vertex alphas could be a function of their probability
	
	x = pos_states[:,1]
	y = pos_states[:,2]
	
	w = get_weight_matrix(stn)
	
	gp = gplot(stn,x,-y;kwargs...,
		edgelinewidth = reduce(vcat, [w[i,:].nzval for i in 1:nv(stn)]) .^ weight_exponent,
		nodelabel = if nodelabels 1:nv(stn) else nothing end,
		nodesize=nodesize,
		nodefillc=nodefillc,
		linetype=linetype,
		EDGELINEWIDTH=max_edgelinewidth)
	
	draw(PDF(filename),gp)

end
