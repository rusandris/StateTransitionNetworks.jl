using Documenter,StateTransitionNetworks


makedocs(
         sitename = "StateTransitionNetworks.jl",
         modules  = [StateTransitionNetworks],
         pages=[
                "Home" => "index.md",
                "Creating STNs" => "create_stn.md",
                "Network measures" => "network_measures.md",
                "API" => "api.md",
               ])
