using DynamicalSystems
using Plots,LaTeXStrings
pyplot()


guidefontsize=20
tickfontsize=15
legendfontsize=15

function plane_surface(x,y;plane)
	a = plane[1:3]
	b = plane[4]
	z = b/a[3] - a[1]*x/a[3] - a[2]*y/a[3]
end


plane = (2,0.0) #y = 0
plot_plane = (0,1,0.01,0.0) #close plane hack


ds = PredefinedDynamicalSystems.roessler()
set_parameter!(ds,2,0.42)

traj, = trajectory(ds,100;Ttr = 500)
psection_plus = DynamicalSystemsBase.poincaresos(traj, plane; save_idxs=[1,2,3],warning=true,direction=1)
psection_minus = DynamicalSystemsBase.poincaresos(traj, plane; save_idxs=[1,2,3],warning=true,direction=-1)


pl = plot(traj[:,1],traj[:,2],traj[:,3],
	xlabel=L"x",
	ylabel=L"y",
	zlabel=L"z",
	guidefontsize=guidefontsize,
	legendfontsize=legendfontsize,
	tickfontsize=tickfontsize,
	lc=:gray10,
	zticks=[0,6,12],
	label="",
	lw=1,
	la=0.9,
	xticks=[-5,0,5],
	yticks=[-5,0,5],
	cam=(-33,29),dpi=300)

plot!(pl,psection_plus[:,1],psection_plus[:,2],psection_plus[:,3],
	st=:scatter,
	mc=:red,
	ms=5,
	ma=1,
	guidefontsize=guidefontsize,
	legendfontsize=legendfontsize,
	tickfontsize=tickfontsize,
	markerstrokewidth=0.00,
	label="",
	cam=(-33,29))
	
plot!(pl,psection_minus[:,1],psection_minus[:,2],psection_minus[:,3],
	st=:scatter,
	mc=:blue,
	ms=5,
	guidefontsize=guidefontsize,
	legendfontsize=legendfontsize,
	tickfontsize=tickfontsize,
	markerstrokewidth=0.00,
	label="",
	cam=(-33,29))

y = [-1,0.001]
x = [-10,12]

surface!(pl,x,y,(x,y)->plane_surface(x,y;plane=plot_plane),c=:greens,alpha=0.4,zlims=(0,12),colorbar=false,clim=(0,0))

savefig(pl,"roessler3D.pdf")

