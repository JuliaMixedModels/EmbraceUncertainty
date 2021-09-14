function stripcodechunks(
    fn::AbstractString;                                # source filename
    chunkre::Tuple{Regex,Regex} = ( r"^<<"m, r"^@$"m), # begin/end chunk delimiters
    fn_ext::AbstractString = ".Rnw",                   # extension for source
    chnk_ext::AbstractString = ".txt",                 # extension for file of chunks
    body_ext::AbstractString = ".tex",
)
    spltnm = splitext(fn)
    if last(spltnm) ≠ fn_ext
        throw(ArgumentError("fn = $fn does not end in $fn_ext"))
    end
    basenm = first(spltnm)
    iochnk = open(string(basenm, chnk_ext), "w")
    iobody = open(string(basenm, body_ext), "w")
    str = read(fn, String)
    offend = 1
    chnkcount = 0
    mstart = match(first(chunkre), str, offend)
    while !isnothing(mstart)  # found an R chunk
        offstart = mstart.offset
        chnkcount += 1
        write(iobody, SubString(str, offend, offstart - 1))
        write(iobody, "\n % chunk $chnkcount here\n")
        write(iochnk, "\n%%%% Chunk $chnkcount\n")
        mend = match(last(chunkre), str, offstart)
        if isnothing(mend)
            throw(error("Code chunk starting at $offset does not appear to end"))
        end
        offend = mend.offset + 1
        write(iochnk, SubString(str, offstart, offend))
        mstart = match(first(chunkre), str, offend)
    end
    write(iobody, SubString(str, offend, length(str)))
    close(iobody)
    close(iochnk)
end

