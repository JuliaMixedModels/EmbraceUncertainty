const ML_25M_URL = "https://files.grouplens.org/datasets/movielens/ml-25m.zip"

const metadata = Dict{String,String}("url" => ML_25M_URL)

function create_arrow(fname, df)
    arrowfile = _file(splitext(basename(fname))[1])
    Arrow.write(arrowfile, df; compress=:lz4, metadata)
    return arrowfile
end

function extract_csv(zipfile, fname; kwargs...)
    file = only(filter(f -> endswith(f.name, fname), zipfile.files))
    return CSV.read(file, DataFrame; delim=',', header=1, kwargs...)
end

const GENRES = ["Action", "Adventure", "Animation",
                "Children", "Comedy", "Crime",
                "Documentary", "Drama",
                "Fantasy", "Film-Noir",
                "Horror",
                "IMAX",
                "Musical", "Mystery",
                "Romance",
                "Sci-Fi",
                "Thriller",
                "War", "Western"]

function load_quiver()
    @info "Downloading data"
    quiver = String[]
    open(Downloads.download(ML_25M_URL), "r") do io
        zipfile = ZipFile.Reader(io)
        @info "Extracting and saving ratings"
        ratings = DataFrame(
            extract_csv(
                zipfile,
                "ratings.csv";
                types=[Int32, Int32, Float32, Int32],
                pool=[true, true, true, false],
            )
        )
        ratings.movieId = tagpad(ratings.movieId, "M")
        ratings.userId = tagpad(ratings.userId, "U")
        push!(quiver, create_arrow("ratings.csv", ratings))
        @info "Extracting movies that are in the ratings table"
        movies = extract_csv(zipfile, "movies.csv"; types=[Int32,String,String], pool=false)
        movies.movieId = tagpad(movies.movieId, "M")
        links = extract_csv(zipfile, "links.csv"; types=[Int32,Int32,Int32])
        links.movieId = tagpad(links.movieId, "M")
        movies = leftjoin!(
            leftjoin!(
                sort!(combine(groupby(ratings, :movieId), nrow => :nrtngs), :nrtngs; rev=true),
                movies,
                on=:movieId,
            ),
            links;
            on=:movieId,
        )
        disallowmissing!(movies; error=false)
        movies.nrtngs = Int32.(movies.nrtngs)
        for g in GENRES
            setproperty!(movies, replace(g, "-" => ""), contains.(movies.genres, g))
        end
        select!(movies, Not("genres"))  # now drop the original genres column
        push!(quiver, create_arrow("movies.csv", movies))
        @info "Extracting and saving README"
        readme = only(filter(f -> endswith(f.name, "README.txt"), zipfile.files))
        open(joinpath(CACHE[], "README.txt"), "w") do io
            write(io, read(readme))
        end

        return nothing
    end

    return quiver
end

function readme()
    return Markdown.parse_file(joinpath(CACHE[], "README.txt"))
end
