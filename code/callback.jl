function callback(pts, lnprobs, n, thin)
    EnsembleSampler.basic_callback(pts, lnprobs, n, thin)
    h5open("state.save.tmp", "w") do f
        f["pts", "shuffle", (), "compress", 3] = pts
        f["lnprobs", "shuffle", (), "compress", 3] = lnprobs
        attrs(f)["nsteps"] = n
        attrs(f)["thin"] = thin
    end
    mv("state.save.tmp", "state.save.hdf5", remove_destination=true)
end
