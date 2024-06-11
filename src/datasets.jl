_file(x) = joinpath(CACHE[], string(x, ".arrow"))

clear_scratchspaces!() = Scratch.clear_scratchspaces!(@__MODULE__)

const datasets = 
    CSV.read(
        IOBuffer(
"""
dsname,filename,version,sha2
box,tkxnh,1,4ac2038735c5286ca0d6a706b4feef5b34bd93560bc4cabc7223addf2366e4c0
elstongrizzle,5vrbw,1,f33e08ad5a91ab5dd9889953dbb97cc52a299fdac9db8b37ec4f87ad2dacadbd
oxboys,cz6g3,1,079f46e404e43c1848a65d3a4e5b6a14b2cb17565c77819d7ee3effe72f5ebd0
sizespeed,kazgm,2,d3512795afc101445fdd78e183919378ef7d993d66b4fd55307a0bf57fba1e0b
ELP_ldt_item,c6gxd,1,f851910e1435659ca662ad49cfb8deb6b7cf287e4ce4969103dba11b32ab2e6c
ELP_ldt_subj,rqenu,2,d9c88915681b64fc9db975f9bb2d6f402058fee5cb35887f9de7d07776efdd56
ELP_ldt_trial,3evhy,2,57a83679f8f458b1d9bb56e05099a556ffdaf15af67464e9ee608c844fc4fa9c
movies,,1,0bdd5e9565328f07bbcfd83c634951ab2933fc160cc07e9bad8200b4c84d90ee
ratings,,1,a2187528c1c6b2a87b1cc09845c0a5591a2fb01088ed7173c187be4a9bee83e2
fggk21,vwecy,1,0fa959f095f8b92135496b6f8c8a8b5a3e896e8875f0ba6928bd074559d8a796
fggk21_Child,c2fmn,1,61c91e00336e6f804e9f6b86986ebb4a14561cc4908b3a21cb27c113d2b51a5c
fggk21_Score,7fqx3,1,99d73ee705aaf5f4ee696eadbba992d0113ba6f467ce337a62a63853e4617400
kkl15,p8cea,2,90d7bb137c8613d7a15c8597c461aee7c7cb0f0989a07c80fc93e1fbe2e5c156
kwdyz11,4cv52,3,2fa23aa8aa25e1adb10183c8d29646ae0d19d6baef9d711c9906f7fa1b225571
"""
        ),
        Table;
        downcast=true,
        pool=false,
    )

if @isdefined(_cacheddatasets)
    empty!(_cacheddatasets)    # start from an empty cache in case datasets has changed
else
    const _cacheddatasets = Dict{Symbol, Arrow.Table}()
end

"""
    dataset(name::Union(Symbol, AbstractString))

Return as an `Arrow.Table` the dataset named `name`.

Available dataset names, their versions, the filenames on the osf.io site and an SHA2 checksum of their contents
are in the table `datasets`.

The files are cached in the scratchspace for this package.  The name of this directory is the value of `CACHE[]`.
"""
function dataset(nm::AbstractString)
    return get!(_cacheddatasets, Symbol(nm)) do  # retrieve from cache if available, otherwise
        # check for nm in datasets table first so MMDS can be overridden
        rows = filter(==(nm) ∘ getproperty(:dsname), datasets)
        if isempty(rows)
            nm in MMDS || error("Dataset '$nm' is not available")
            MixedModelsDatasets.dataset(nm)
        else
            row = only(rows)       # check that there is only one matching row and extract it
            fnm = _file(nm)
            if !isfile(fnm) || row.sha2 ≠ bytes2hex(open(sha2_256, fnm))
                if ismissing(row.filename)
                    load_quiver()  # special-case `ratings` and `movies`
                else
                    Downloads.download(
                        string("https://osf.io/", row.filename, "/download?version=", row.version),
                        fnm,
                    )
                end
            end
            if row.sha2 ≠ bytes2hex(open(sha2_256, fnm))
                throw(error("Invalid checksum on downloaded $nm dataset, version $(row.version)"))
            end
            Arrow.Table(fnm)
        end
    end
end
dataset(nm::Symbol) = dataset(string(nm))
