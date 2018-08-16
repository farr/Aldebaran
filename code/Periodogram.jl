module Periodogram

function generate_mean_vectors(ts)
    nv = size(ts, 1)

    ntot = sum([size(t,1) for t in ts])

    mvs = zeros(ntot, nv)

    i = 1
    for j in eachindex(ts)
        t = ts[j]

        k = i + size(t, 1) - 1
        mvs[i:k, j] = 1.0

        i = k+1
    end

    mvs
end

function response_matrix(ts, ys, dys, f)
    mvs = generate_mean_vectors(ts)

    all_ts = vcat(ts...)
    all_ys = vcat(ys...)
    all_dys = vcat(dys...)

    cv = cos.(all_ts.*(2.0*pi*f))
    sv = sin.(all_ts.*(2.0*pi*f))

    all_ts, all_ys, all_dys, hcat(mvs, cv, sv)
end

function response_matrix(ts, ys, dys, fs::Array{Float64, 1})
    mvs = generate_mean_vectors(ts)

    all_ts = vcat(ts...)
    all_ys = vcat(ys...)
    all_dys = vcat(dys...)

    csvs = zeros(size(all_ts, 1), 2*size(fs, 1))

    for i in eachindex(all_ts)
        for j in eachindex(fs)
            csvs[i,2*j-1] = cos(2.0*pi*fs[j]*all_ts[i])
            csvs[i,2*j] = sin(2.0*pi*fs[j]*all_ts[i])
        end
    end

    all_ts, all_ys, all_dys, hcat(mvs, csvs)
end

function pgram(ts, ys, dys, f)
    all_ts, all_ys, all_dys, A = response_matrix(ts, ys, dys, f)

    for i in eachindex(all_ys)
        all_ys[i] *= 1.0./all_dys[i]
        A[i,1:end] *= 1.0./all_dys[i]
    end

    A\all_ys
end

function residual(ts, ys, dys, f, pgram)
    all_ts, all_ys, all_dys, A = response_matrix(ts, ys, dys, f)

    all_rs = all_ys - A*pgram

    rs = Array{Float64, 1}[]
    i = 1
    for j in eachindex(ts)
        k = i + size(ys[j],1) - 1
        push!(rs, all_rs[i:k])
        i = k+1
    end

    rs
end

function basis_pursuit_pgram(ts, ys, dys, fs_init, df_window, nbasis)
    fspeak = Float64[]

    while length(fspeak) < nbasis
        As = Float64[]
        Bs = Float64[]

        sel = trues(length(fs_init))
        for f in fspeak
            sel = sel .& (abs.(fs_init - f) .> df_window)
        end

        for f in fs_init[sel]
            ff = vcat(fspeak, [f])
            c = pgram(ts, ys, dys, ff)
            push!(As, c[end-1])
            push!(Bs, c[end])
        end
        power = sqrt.(As.*As + Bs.*Bs)

        push!(fspeak, fs_init[sel][indmax(power)])
    end

    fspeak
end

end
