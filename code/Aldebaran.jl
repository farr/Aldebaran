module Aldebaran

using Kepler
using CARMAKepler

function load_timeseries(before_song=true, song=true)
    ts = []
    ys = []
    dys = []

    if before_song 
        for i in 3:9
            data = readdlm("/Users/farr/Documents/Research/Aldebaran/data/$(i)/table$(i).dat")
            push!(ts, data[:,1])
            push!(ys, data[:,2])
            push!(dys, data[:,3])
        end
    end

    if song
        data = readdlm("/Users/farr/Documents/Research/Aldebaran/data/song/tablesong.dat")
        push!(ts, data[:,1])
        push!(ys, data[:,2])
        push!(dys, data[:,3])
    end

    (ts, ys, dys)
end

const P_min = 300.0
const P_max = 1200.0

const K_min = 70.0
const K_max = 280.0

const mu_Hz = 1e-6*3600.0*24.0 # per day

# We think nu_max is ~2 mu_Hz
const f_min = 1.0*mu_Hz
const f_max = 5.0*mu_Hz

end
