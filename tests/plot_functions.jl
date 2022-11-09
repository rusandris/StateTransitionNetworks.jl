function plot_phase_portrait(timeseries::Matrix,i,j;filename)
	dim = size(timeseries)[2]
	pl = plot(timeseries[:,i],timeseries[:,j],
	lc=:gray10,
	lw=2,xlabel="$i",
	ylabel="$j",
	framestyle=:origin,
	xguidefontsize=18,
	yguidefontsize=18,
	tickfontsize=10,
	legend=false,
	dpi=300)
	
	savefig(pl,filename)
end

function plot_timeseries(timeseries::Matrix;filename,dims)
	dim = size(timeseries)[2]
	pl = plot(timeseries[:,dims],
	legend=false,
	lw=1,
	xguidefontsize=18,
	yguidefontsize=18,
	tickfontsize=10,
	dpi=300)
	
	savefig(pl,filename)
end

function plot_psection(timeseries::Matrix;plane,idxs,filename)
	dim = size(timeseries)[2]
	psection = poincaresos(Dataset(timeseries), plane;idxs=idxs); 
	i = idxs[1]
	j = idxs[2]
	
	
	pl = plot(psection[:,1],psection[:,2],
	st=:scatter,
	ms=1,
	markerstrokewidth=0.001,
	markeralpha=0.5,
	mc=:gray10,
	xlabel="$i",
	ylabel="$j",
	framestyle=:origin,
	xguidefontsize=18,
	yguidefontsize=18,
	tickfontsize=10,
	legend=false,
	dpi=300)
	
	savefig(pl,filename)
end

function plot_degree_distribution(stn;bins=nothing,xaxis,yaxis,mc,filename)
	if bins != nothing
		h = fit(Histogram,outdegree(stn),nbins=bins)
	else
		h = fit(Histogram,outdegrees(stn))
	end
	
	h = normalize(h,mode=:probability)
	bin_centers = h.edges[1][1:end-1] + diff(h.edges[1])/2
    dd = plot(bin_centers,h.weights,
    st=:scatter,
    ms=3,
    mc=:gray10,
    xaxis=xaxis,
    yaxis=yaxis,
    legend=false)

	savefig(dd,filename)
end

function plot_degree_distribution!(stn;plot,bins=nothing,xaxis,yaxis,mc,label)
	if bins != nothing
		h = fit(Histogram,outdegree(stn),nbins=bins)
	else
		h = fit(Histogram,outdegrees(stn))
	end
	
	h = normalize(h,mode=:probability)
	bin_centers = h.edges[1][1:end-1] + diff(h.edges[1])/2
	
    dd = plot!(plot,bin_centers,h.weights,
    st=:scatter,
    mc=mc,
    ms=4,
    xaxis=xaxis,
    markerstrokealpha=0.001,
    yaxis=yaxis,
    label = label,
    xlabel=L"$k$",
	ylabel=L"$P(k)$",
	xguidefontsize=18,
	yguidefontsize=18,
	tickfontsize=10,
	dpi=300)

end

function histogram_degree_distribution(stn;bins=:auto,xaxis,yaxis,filename)
	h = histogram(outdegree(stn),
	normalize=:probability,
	bins = bins,
	xaxis = xaxis,
	yaxis = yaxis,
	xlabel=L"$k$",
	ylabel=L"$P(k)$",
	xguidefontsize=18,
	yguidefontsize=18,
	tickfontsize=10,
	legend=false)
	
	savefig(h,filename)
end


