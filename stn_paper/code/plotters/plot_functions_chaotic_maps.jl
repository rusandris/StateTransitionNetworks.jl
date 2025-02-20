function plot_measure(ps,ms,special_ps;λs=[],labels,vertical_lw,
    red_params=[],inset_param_range,inset_data=[],inset_box,inset_ticks,inset_tickfontsize,inset_ylims,alphas,marker_colors,orders,kwargs...)


    pl = plot(;kwargs...)


    ylims = kwargs[:ylims]


    #plot measure curves
    for i in 1:length(orders)
        plot!(pl,ps,ms[:,i],lw=2,la=alphas[i],lc=:gray50,label = labels[i])
    end


    #plot lyap_exps
    if !isempty(λs)
        plot!(pl,ps,λs,lw=1,lc=:gray10,label=L"\lambda")
    end

    #plot vertical lines on S plot
    for (i,p) in enumerate(special_ps)
        plot!(pl,[p,p],[ylims[1],ylims[end]],ls=:dash,lw=vertical_lw,lc=marker_colors[i],label=false)
    end

    #plot inset
    inset_xticks,inset_yticks = inset_ticks
    plot!(pl, inset=inset_box, subplot=2, framestyle=:box,xticks=inset_xticks,yticks=inset_yticks,ylims=inset_ylims,tickfontsize=inset_tickfontsize)

    

    fi = findfirst(p -> p >= inset_param_range[1],ps)
    li = findlast(p -> p <= inset_param_range[end],ps)
    inset_indices = fi:li
    
    if isempty(inset_data)
        plot!(pl[2], ps[inset_indices], ms[inset_indices,2*length(orders)], label=nothing ,st=:scatter, markerstrokewidth=0, mc=:gray60, ms=4, ma=1) #plot highest order numerical
        plot!(pl[2], ps[inset_indices], ms[inset_indices,length(orders)], label=nothing, lw=1, color=:gray10) #plot highest order analytical
    else
        ps_inset = inset_data[:,1]
        ms_inset = inset_data[:,2]
        rw_vals_inset = inset_data[:,3]
     
        plot!(pl[2], ps_inset, rw_vals_inset, label=nothing ,st=:scatter, markerstrokewidth=0, mc=:gray60, ms=4, ma=1) #plot highest order numerical
        plot!(pl[2], ps_inset, ms_inset, label=nothing, lw=1, color=:gray10) #plot highest order analytical
    end



    if !isempty(red_params)
        fi = findfirst(p -> p >= red_params[1],ps)
        li = findlast(p -> p <= red_params[end],ps)

        plot!(pl[2], ps[fi:li], ms[fi:li,length(orders)], label=nothing, lw=1, color=:red)
    end

    return pl
end

function plot_orbit_diagram(od,ps,special_ps;ms_od,ma_od,marker_shapes,marker_size,marker_colors,marker_offset=0.1,kwargs...)
    #--------------------------plot orbit diagram------------------------
    od_plot = plot(;kwargs...)
    n = length(od[1])
    od_ylims = plot_params_od[:ylims]

    for (i, p) in enumerate(ps)

        plot!(od_plot,fill(p,n),od[i],
        st=:scatter,
        ms = ms_od,
        ma = ma_od,
        mc =:gray10,
        markerstrokewidth=0.0,
        legend=false)       

    end

    #plot markers on od
    plot!(od_plot,special_ps,fill(od_ylims[2]-marker_offset,length(special_ps)),
    markershape=marker_shapes,
    mc=marker_colors,
    ms=marker_size,
    st=:scatter,
    markerstrokewidth=0.5,
    markerstrokecolor=:gray10,
    xformatter=:none)

    #plot vertical lines on od
    for (i,p) in enumerate(special_ps)
        plot!(od_plot,[p,p],[od_ylims[1],od_ylims[2]-marker_offset],ls=:dash,lw=2,lc=marker_colors[i])
    end

    return od_plot

end



function plot_renyi_spectrums(qs,renyi_spectrums,Ss,ps,ref_vals;marker_shapes,marker_colors,alphas=ones(length(ps)),marker_size=10,annotation_fontsize=22,annotation_factor=0.9,labels,plot_ref_values=true,kwargs...)

    pl = plot(;kwargs...) 

    #max_renyi = maximum(renyi_spectrums)

    spec_ylims = kwargs[:ylims]

    #vertical lines
    plot!(pl,[0.0,0.0],[0.0,spec_ylims[end]],ls=:dash,lc=:gray10,label="")
    plot!(pl,[1.0,1.0],[0.0,spec_ylims[end]],ls=:dash,lc=:gray10,label="")

    for (i,p) in enumerate(ps)

        plot!(qs,renyi_spectrums[:,i],lw=3,lc=marker_colors[i],alpha = alphas[i],label="")

        S = Ss[i]
        
        plot!([1.0],[S],
        ms=marker_size,
        markershape = marker_shapes[i],
        alpha = alphas[i],
        st=:scatter,
        mc = marker_colors[i],
        markerstrokewidth=0.5,
        markerstrokecolor=:gray10,
        label=labels[i])
        
    end

    if plot_ref_values == true
        h_top,λ = ref_vals

        #point for reference topological entropy
        plot!([0.0],[h_top],st=:scatter,mc=:gray10,ms=4,markerstrokewidth=0.0,label="") #L"h_t=0.46469")
        annotate!(0.05, annotation_factor*h_top , text(L"h_t", :gray10, :left, annotation_fontsize))
        
        #point for lyap exp 
        plot!([1.0],[λ],st=:scatter,mc=:gray10,ms=4,markerstrokewidth=0.0,label="")
        annotate!(1.05, annotation_factor*λ, text(L"\lambda", :gray10, :left, annotation_fontsize))

        
        #other annotations, arrows
        #annotate!(0.0, spec_ylims[end], text(L"\bar{K}_0", :gray10, :left,annotation_fontsize))
        #annotate!(1.0, spec_ylims[end], text(L"\bar{K}_1", :gray10, :left, annotation_fontsize))
        #quiver!([1.3], [λ - 0.1], quiver=([-0.3], [0.11]),color=:gray10)
        #quiver!([0.4], [h_top - 0.1], quiver=([-0.3], [0.06]),color=:gray10)
    end

    return pl

end


function plot_renyi_spectrums!(pl,qs,renyi_spectrums,Ss,ps,ref_vals;marker_shapes,marker_colors,alphas,labels,plot_ref_values=true,kwargs...)
    
    max_renyi = maximum(renyi_spectrums)

    for (i,p) in enumerate(ps)
        #plot renyi spectrum curve
        plot!(pl,qs,renyi_spectrums[:,i],lw=2,lc=marker_colors[i],alpha = alphas[i],label="")

        #plot measure at q=1
        S = Ss[i]
        
        plot!(pl,[1.0],[S],
        ms=6,
        markershape = marker_shapes[i],
        alpha = alphas[i],
        st=:scatter,
        mc = marker_colors[i],
        markerstrokewidth=0.0,
        label=labels[i])
        
    end

    if plot_ref_values == true
        h_top,λ = ref_vals

        plot!(pl,[0.0],[h_top],st=:scatter,mc=:gray10,ms=4,markerstrokewidth=0.0,label="") #L"h_t=0.46469")

        annotate!(pl,0.0, max_renyi, text(L"\bar{K}_0", :gray10, :left, 18))
        annotate!(pl,1.0, max_renyi, text(L"\bar{K}_1", :gray10, :left, 18))
        annotate!(pl,0.1, h_top , text(L"h_t", :gray10, :left, 18))
        #quiver!([0.4], [h_top - 0.1], quiver=([-0.3], [0.06]),color=:gray10)

        plot!(pl,[1.0],[λ],st=:scatter,mc=:gray10,ms=4,markerstrokewidth=0.0,label="")#L"S_{KS}=$(S_KS_approx)")
        annotate!(pl,1.1, λ, text(L"\lambda", :gray10, :left, 18))
        #quiver!([1.3], [λ - 0.1], quiver=([-0.3], [0.11]),color=:gray10)
    end

end


