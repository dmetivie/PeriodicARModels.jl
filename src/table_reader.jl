#Packages
using CSV, DataFrames, Dates, DataFramesMeta

"""
    truncate_MV(df)

Truncates the dataframe of temperatures series `df` to remove missing values. If there is a missing value in the last 1000 values, it removes all the values after. If there is a missing value before, the function removes all the values before.
"""
function truncate_MV(df_full)
    df, f = copy(df_full), true
    while f
        if any(diff(df.DATE[end-1000:end]) .> Day(1)) #If the last MV is close to the end the series...
            df = @chain df begin
                @transform(:diff = [diff(:DATE); Day(1)])
                @aside beg = _.DATE[findlast(_.diff .> Day(1))]
                @subset(:DATE .<= beg) #in a first time we keep what is before it
                @select Not(:diff)
            end
        elseif any(diff(df.DATE) .> Day(1)) #Else we keep what is after.
            df = @chain df begin
                @transform(:diff = [diff(:DATE); Day(1)])
                @aside beg = _.DATE[findlast(_.diff .> Day(1))]
                @subset(:DATE .> beg)
                @select Not(:diff)
            end
            f = false #to stop the loop
        else
            f = false
        end
    end
    if Date(year(df.DATE[1]), month(df.DATE[1])) + Month(1) - df.DATE[1] < Day(20) #If there are not enough days in the first month I remove it.
        df = df[df.DATE.>=(Date(year(df.DATE[1]), month(df.DATE[1]))+Month(1)), :]
    end
    if df.DATE[end] - Date(year(df.DATE[end]), month(df.DATE[end])) < Day(19) #The same thing for the last month.
        df = df[df.DATE.<Date(year(df.DATE[end]), month(df.DATE[end])), :]
    end
    return df
end

"""
    extract_series(file::String; year=nothing)

Extracts the temperatures time series from `file` and return a dataframe. You can precise the year with the argument `year`. If `truncate_results = true`, truncates the series with the function [`truncate_MV`](@ref).
If the file name starts with "TN", "TG" or "TX", the function considers that it is an ECA&D file and returns a dataframe with a date column and one temperature column of the type precised in the name of the file.
Otherwise it considers that it is a file from the DRIAS database so it returns a dataframe with one column for the date and three columns for the three types of daily temperatures.
"""
function extract_series(file::String; truncate_results=true, year=nothing)
    file_name = splitdir(file)[end]
    type_data = file_name[1:2]
    if type_data ∈ ["TN", "TG", "TX"]#startswith(file_name, "TN") || startswith(file_name, "TG") || startswith(file_name, "TX") #If it's an ECA&D file
        table = CSV.read(file, DataFrame, normalizenames=true, skipto=22, header=21, ignoreemptyrows=true, dateformat="yyyymmdd", types=Dict(:DATE => Date))
        df = table[table[!, end].==0, [end - 2, end - 1]]
        df[!, 2] /= 10
    elseif type_data == "er" #If it's a era5 file
        df = CSV.read(file, DataFrame, normalizenames=true, ignoreemptyrows=true, dateformat="yyyy-mm-dd", types=Dict(:DATE => Date))
    elseif type_data == "RC" #If it's a rcp file
        df = CSV.read(file, DataFrame, normalizenames=true, comment="#", ignoreemptyrows=true, dateformat="yyyymmdd", types=Dict(:DATE => Date))
    else #If it's a DRIAS file
        table = CSV.read(file, DataFrame, normalizenames=true, skipto=49, header=48, ignoreemptyrows=true, dateformat="yyyymmdd", types=Dict(:DATE => Date))
        df = @chain table begin
            @rsubset :TN != -999.99
            @rsubset :TG != -999.99
            @rsubset :TX != -999.99
        end
    end
    # df.DATE = Date.(string.(df.DATE), dateformat"yyyymmdd")
    df = isnothing(year) ? df : df[Date(year).<=df.DATE.<Date(year + 1), :]
    return truncate_results ? truncate_MV(df) : df
end


"""
    Common_indexes(series_vec::AbstractVector{DataFrame})
    Common_indexes(files...)

From different ECA&D-type dataframes (with 2 colums) in `series_vec` or extracted in `files` return a tuple with the series of their mutual date and their associated temperatures in a matrix.
For example if we consider a dataframe `df_TN` of daily minimal temperatures from 1985 to 2010 and another one `df_TX` of daily maximal temperatures from 1970 to 2005, `Common_indexes(df_TN,df_TX)` will return a tuple `(date_vec, x)` where :
- `date_vec` is the dates vector (each day as `Date` object) for the period 1985-2005.
- `x` is a matrix where the first column is the series of daily minimal temperatures and the second the series of daily maximal temperatures, both corresponding to the dates in date_vec.
"""
function Common_indexes(series_vec::AbstractVector{DataFrame})
    Date_vecs = [series.DATE for series in series_vec]
    if all(y -> y == Date_vecs[1], Date_vecs)
        Temps_vecs = [series[:, 2] for series in series_vec]
        return Date_vecs[1], stack(Temps_vecs)
    else #If the timelines are differents, we take the common timeline of the two series.
        date_vec = maximum(Date_vec[1] for Date_vec in Date_vecs):minimum(Date_vec[end] for Date_vec in Date_vecs)
        x = stack(series[:, 2][findfirst(series.DATE .== date_vec[1]):findfirst(series.DATE .== date_vec[end])] for series in series_vec)
        return date_vec, x
    end
end
Common_indexes(files...) = Common_indexes(extract_series.(collect(files))) #println(files == (joinpath(StationsPath,"TN_Montpellier.txt"), joinpath(StationsPath,"TX_Montpellier.txt")))
Common_indexes(files::AbstractVector{String}) = Common_indexes(extract_series.(files))



