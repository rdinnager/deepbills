library(JuliaCall)
julia_setup()
julia_library("PANDA")

julia_command('test = load("data/clads_Birds.jld2", "output")')
julia_command('save_ClaDS_in_R(test, "data/clad_Birds.rData")')