module fitsIO
    using FITSIO, TextParse, AstroLib, DataFrames, Dates, WeakRefStrings, Printf
    using ProgressMeter, PrettyTables, Dierckx, Gnuplot, QuadGK

    import Base.write

    export adaptive_bin,
        create_dataframe,
        find_fits,
        get_current_date,
        get_hdu_col_names,
        get_hdu_keys,
        get_header_values,
        get_key_values,
        image_to_table,
        is_known,
        mag_lim,
        new_spectra_df!,
        prepare_stack,
        rebin,
        remove_missing!,
        stack_spectra!,
        sub_in_df!,
        write

    """
        adstring_temp(ra, dec)

    Piccolo hack per far funzionare correttamente adstring() e far si che non impazzisca con
    coordinate del tipo -0.qualcosa.
    """
    function adstring_temp(ra, dec)
        coord = adstring(ra, dec, precision = 1)
        if dec > 0
            return coord
        else
            ras, decs = split(coord, "  ")
            decs = '-'*strip(decs, ['-', '+'])
            return string(ras, "  ", decs)
        end
    end

    """
        get_hdu_keys(hdu::TableHDU)

    Ritorna tutte le chiavi per l'HDU.
    """
    function get_hdu_keys(hdu::TableHDU)
        return keys(read_header(hdu))
    end

    # ------------------------------ ** ------------------------------ #

    """
        get_key_values(hdu::TableHDU, keys::Tuple)

    Ritorna una coppia di (value[key], comment[value]) per ogni key.
    Se serve una sola chiave, passa (key, ) invece che (key) per evitare casotto. Usa più
    che puoi le keywords invece degli indici che son difficili da maneggiare.
    """
    function get_key_values(hdu::TableHDU, keys::Tuple)
        return [read_key(hdu, key) for key in keys]
    end

    # ------------------------------ ** ------------------------------ #

    """
        get_hdu_col_names(hdu::TableHDU)

    Ritorna i nomi delle colonne del'hdu.
    """
    function get_hdu_col_names(hdu::TableHDU)
        return FITSIO.colnames(hdu)
    end

    # ------------------------------ ** ------------------------------ #
    # Si occupa della conversione di tipi nella lettura dei dataframe

    """
        date_to_string(d)

    Converte un tipo Dates in una string, più semplice da gestire.
    """
    function date_to_string(d::TimeType)
        yy = string(year(d))
        mm = lpad(month(d), 2, "0")
        dd = lpad(day(d), 2, "0")
        h = lpad(hour(d), 2, "0")
        m = lpad(minute(d), 2, "0")
        s = lpad(minute(d), 2, "0")
        return yy*'-'*mm*'-'*dd*'T'*h*':'*m*':'*s
    end

    # ------------------------------ ** ------------------------------ #

    """
        get_current_date()

    Ritorna la data attuale come stringa.
    """
    function get_current_date()
        date_to_string(now())
    end
"""
        create_dataframe_internal(hdu::TableHDU, internal_id = "myid")

    Crea e ritorna un dataframe dato un `TableHDU` da un file Fits.
    """
    function create_dataframe_internal(hdu::TableHDU, internal_id = "myid")
        col_data = []
        col_names = get_hdu_col_names(hdu)
        push!(col_data, collect(1:length(read(hdu, col_names[1]))))
        # Creo la colonna con id progressivi per i join, in quanto mi risulta impossbile
        #  preservare il numero di dati in caso di spettri multipli
        for col_name in col_names
            push!(col_data, read(hdu, col_name))
        end
        col_names = vcat(internal_id, col_names)
        df = DataFrame(col_data, Symbol.(col_names))
        return df
    end

    # ------------------------------ ** ------------------------------ #

    """
        prepare_text(file_in, separator::Char = '|', file_out::AbstractString = "temp.gitignore"), cchar::Char = '#')

    Riscrive il file di input in modo che sia ancora un csv, ma senza spazi tra
    una entry e l'altra, e con `|` (pipe) come delimitatore consistente ovunque.
    """
    function prepare_text(file_in, separator::Char = '|', file_out::AbstractString = "temp.gitignore",
                            cchar::Char = '#')
        format_separator = separator
        if occursin(separator, "^.[\$()|*+?")
            format_separator = "\\$separator"
        end
        temp_file = open(file_out, "w")

        regex_entry = Regex("(?<=$format_separator)[^$format_separator]*")
        regex_line_start = Regex("^($format_separator)")

        for line in eachline(file_in)
            # Sistema la linea in modo da poter usare sempre la stessa regex
            if !occursin(regex_line_start, line)
                line = "$separator$line"
            end

            reg_match = eachmatch(regex_entry, line)
            for (i, match) in enumerate(reg_match)
                if i < length(collect(reg_match))
                    write(temp_file, strip(match.match)*'|')
                else
                    write(temp_file, strip(match.match))
                end
            end
            write(temp_file, '\n')
        end
        close(temp_file)
    end

    # ------------------------------ ** ------------------------------ #

    # Nota bene: assume che la prima riga siano gli header del csv
    """
        create_dataframe_internal(file_in; separator::Char = '|', cchar::Char = '#', hexists::Bool = true, dict = nothing, internal_id = "myid")

    Crea il dataframe, dato il file pulito preparato da `prepare_text`.
    """
    function create_dataframe_internal(file_in; separator::Char = '|', cchar::Char = '#',
                hexists::Bool = true, dict = nothing, internal_id = "myid")
        (data, col) = csvread(file_in, separator, commentchar = cchar, header_exists = hexists)
        data_ = []
        # Cambia da WeakStringRef in Stringhe normali
        for item in data
            if isa(item, StringArray)
                push!(data_, convert(Array{String}, item))
            else
                push!(data_, item)
            end
        end

        col_names = vcat(Symbol(internal_id), Symbol.(col))
        test = collect(1:length(data[1]))
        ndata = pushfirst!(data_, test)
        data_frame = DataFrame(ndata, col_names)
        if dict !== nothing
            for col in col_names
                if haskey(dict, col)
                    data_frame[!, col] = dict[col].(data_frame[:, col])
                end
            end
        end
        return data_frame
    end


    # ------------------------------ ** ------------------------------ #

    """
        find_fits(parent_folder::AbstractString, file_out::AbstractString = "/tmp/fits_files.txt")

    Ritorna tutti i file fits che ci sono all'interno della cartella indicata, e ne scrive i
    nomi su file.
    """
    function find_fits(parent_folder::AbstractString)
        run(pipeline(`find $parent_folder -iname "*.fits"`, stdout = "/tmp/fits_files.txt"))
        files = read_keyword("/tmp/fits_files.txt")
        return files
    end

    # ------------------------------ ** ------------------------------ #


    """
        create_dataframe(file_in::AbstractString, file_type::Char; hdu_number::Int = 2, separator::Char = '|', cchar::Char = '#', hexists::Bool = true, file_out::AbstractString = "temp", dict = nothing)

    filetype: `t` per testo, `f` per fits. Si occupa di creare il dataframe una volta
    che viene fornito file e tipo di file in input. Nel caso sia un file di testo,
    il separatore default è la `|` pipe, un commento è segnalato da `cchar` e la
    presenza/assenza di header (nomi colonna) è segnalata da `hexists`.
    """
    function create_dataframe(file_in::AbstractString, file_type::Char; hdu_number::Int = 2,
        separator::Char = '|', cchar::Char = '#', hexists::Bool = true,
        file_out::AbstractString = "temp", dict = nothing)
        # separator::Char = '│' <- Se servisse l'altro carattere
        if file_type == 't'
            prepare_text(file_in, separator, file_out, cchar)
            new_frame = create_dataframe_internal(file_out, cchar = cchar, hexists = hexists, dict = dict)
            run(`rm -f temp`)
            return new_frame
        elseif file_type == 'f'
            f = FITS(file_in)
            return_ = create_dataframe_internal(f[hdu_number])
            close(f)
            return return_
        end
    end

    # ------------------------------ ** ------------------------------ #

    """
        write(df::DataFrame, file_name::AbstractString, internal_id = "myid")

    Scrive il file fits dal dataframe.
    """
    function write(df::DataFrame, file_name::AbstractString, internal_id = "myid")
        df = deepcopy(df)

        if internal_id in names(df)
            select!(df, Not(internal_id))
        end

        data = []
        field_names = names(df)
        for symbol in names(df)
            push!(data, df[:, symbol])
        end
        f = FITS(file_name, "w")
        write(f, field_names, data)
        close(f)
    end

    # ------------------------------ ** ------------------------------ #

    """
        is_close(coord_one::Tuple, coord_two::Tuple, max_distance::Int64 = 1)

    Determina se due sorgenti sono vicine tra loro. Se lo sono ritorna `true`, altrimenti
    `false`. Il grosso del lavoro viene fatto da `add_new_survey`.
    """
    function is_close(coord_one::Tuple, coord_two::Tuple, max_distance::Int64 = 1)
        distance = gcirc(2, coord_one, coord_two)
        if distance < max_distance
            return true
        end
        return false
    end

    # ------------------------------ ** ------------------------------ #

    """
        is_known(QSO_coord::Tuple, QSO_survey::DataFrame)

    Passando le coordinate di un oggetto (`QSO_coord`) e un dataframe con tutte le
    coordinate degli oggetti noti ritorna `(true, match_index)` se l'oggetto è noto,
    `(false, 0)` altrimenti.
    """
    function is_known(qso_coord::Tuple, qso_db::DataFrame)
        for (i, row) in enumerate(eachrow(qso_db))
            if is_close(qso_coord, (row[:RAd], row[:DECd]))
                return (true, i)
            end
        end
        return (false, 0)
    end

    # ------------------------------ ** ------------------------------ #

    """
    Estende collect in modo da permettere di ottenere i tipi di unione. Se hai per es.
    `Union{Float, Int}` ritorna un array `[Float, Int]`.
    """
    Base.collect(t::Union{Type, DataType, Union{}}) = _collect(t, [])
    _collect(t::Type, list) = t<:Union{} ? push!(list, t) : _collect(t.b, push!(list, t.a))
    _collect(t::Union{DataType,Core.TypeofBottom}, list) = push!(list, t)

    # ------------------------------ ** ------------------------------ #

    """
        remove_missing!(df::DataFrame)

    Rimuove i `missing` sostituendoli con una entry apposita. Per i float passa un `NaN`, per
    una stringa passa ` `, mentre per gli int passa `0`.
    """
    function remove_missing!(df::DataFrame)

        function parse_missing(x::Missing, type)
            if type == String
                return " "
            elseif type == Float64
                return NaN
            elseif type == Int64
                return 0
            end
        end
        parse_missing(x::Number, type) = x
        parse_missing(x::String, type) = x

        col_names = names(df)
        for name in col_names
            col_type = eltype(df[!, name])
            if col_type == Missing
                select!(df, Not(name))
            elseif isa(col_type, Union)
                non_missing = collect(col_type)[2]
                df[!, name] = parse_missing.(df[!, name], non_missing)
            end
        end
    end

    # ------------------------------ ** ------------------------------ #

    """
    write(table_::DataFrame, file_out::String, mod::String = t, internal_id::String = "myid")

    Stampa a file un DataFrame.
    """
    function write(table_::DataFrame, file_out::String, mod::Char,
        internal_id::String = "myid")
        
        table = deepcopy(table_)
        if internal_id in names(table)
            select!(table, Not(internal_id))
        end

        open(file_out, "w") do f
            pretty_table(f, table; tf = markdown)
        end
    end

    # ------------------------------ ** ------------------------------ #

    """
        shift_normalize_spectra(spec::DataFrame, z::Union{Float64, Int64};
             normalize_at::Union{Float64, Int64} = 1450)

    Shifta lo spettro al restframe e lo normalizza alla lunghezza d'onda indicata.
    Lo spettro che viene passato deve essere già opportunamente rebinnato.
    """
    function shift_normalize_spectra(spec_::DataFrame, z::Union{Float64, Int64};
        normalize_at::Union{Float64, Int64} = 1450)

        step = minimum(spec_[!, :wpix])
        spec = copy(spec_)

        spec[!, :wave] = spec[!, :wave]./(1 + z)
        # Cerco di evitare che eventuali crolli a 1450 mi uccidando lo stacking
        new_wave = spec[!, :wave][1]:35:spec[!, :wave][end]
        smoothed_spec = rebin(new_wave, spec)
        spline = Spline1D(smoothed_spec[!, :wave], smoothed_spec[!, :flux])
        spec[!, :flux] ./= spline(normalize_at)
        #spec[!, :wpix] ./= spline(normalize_at)
        #spec[!, :flux] ./= normalizing_flux(spec, normalize_at, step)
        #spec[!, :wpix] ./= normalizing_flux(spec, normalize_at, step)

        return spec
    end

    # ------------------------------ ** ------------------------------ #

    #TODO: Testa con intervalli irregolari
    """
        rebin(new_wave::AbstractVector, spec::DataFrame)

    Rebinna basandosi su un nuovo intervallo di lunghezze d'onda fornito.
    `new_wave` è il nuovo intervallo di lunghezze d'onda, da passare per esempio come
    start:step:end, spec lo spettro iniziale.
    Le unità sono arbitrarie. Ritorna un nuovo spettro come dataframe.

    È necessario testare se funziona con intervalli non regolari.
    """
    function rebin(new_wave::AbstractVector, spec::DataFrame)
        wave = spec[:, :wave]
        flux = spec[:, :flux]
        dw = (wave[2] - wave[1])/2
        ndw = (new_wave[2] - new_wave[1])/2

        transform_matrix = zeros(Float64, (length(new_wave), length(wave)))
        sum_check = ndw/dw

        for (j, lj) in enumerate(new_wave)
            nw_lower = lj - ndw
            nw_upper = lj + ndw
            for (i, li) in enumerate(wave)
                w_lower = li - dw
                w_upper = li + dw
                if isapprox(w_lower, nw_lower) && isapprox(w_upper, nw_upper)
                    transform_matrix[j, i] = 1
                elseif w_lower > nw_lower && w_upper < nw_upper
                     transform_matrix[j, i] = 1
                elseif w_lower < nw_upper && w_upper > nw_upper
                    transform_matrix[j, i] = (nw_upper - w_lower)/dw
                elseif w_lower < nw_lower && w_upper > nw_upper
                    transform_matrix[j, i] = (w_upper - nw_lower)/dw
                else
                    continue
                end

                if isapprox(sum(transform_matrix[j, :]), sum_check)
                    break
                end
            end
        end

        new_flux = transform_matrix * flux
        for j in 1:length(new_wave)
            new_flux[j] = new_flux[j]/sum(transform_matrix[j, :])
        end

        return DataFrame([new_wave, new_flux], [:wave, :flux])
    end

    # ------------------------------ ** ------------------------------ #

    """
        prepare_stack(spec::DataFrame, z::Union{Float64, Int64}, rebin_method; grid_low_ = 900, grid_high_ = 2000, grid_step_ = 2, normalize_at_ = 1450)

    Prepara uno spettro per lo stacking, shiftandolo al restframe, normalizzandolo e
    rebinnandolo. spec è lo spettro in ingresso, z il redshift osservato, `rebin_method`
    permette di scegliere se rebinnare tramite interpolazione oppure altro metodo.
    `grid_low`, `grid_step` e `grid_high` indicano gli estremi e lo step per il rebinning.
    `normalize_at` infine indica a quel lunghezza d'onda normalizzare lo spettro.
    Ritorna un nuovo oggetto.
    """
    function prepare_stack(spec::DataFrame, z::Union{Float64, Int64}, rebin_method = rebin;
        grid_low_::Union{Float64, Int64} = 900, grid_high_::Union{Float64, Int64} = 2000,
        grid_step_::Union{Float64, Int64} = 2, normalize_at_::Union{Float64, Int64} = 1450)

        spec_ = shift_normalize_spectra(spec, z; normalize_at = normalize_at_)

        new_wave = grid_low_:grid_step_:grid_high_
        rebinned_spec = rebin_method(new_wave, spec_)
        return rebinned_spec
    end

    # ------------------------------ ** ------------------------------ #

    """
        stack_spectra!(stack::DataFrame, count::DataFrame, spec_list::Array{DataFrame})

    Crea il vero stacking. `stack` deve essere un dataframe con `wave` e `flux` (estremi a
    piacere), ma con la colonna `flux` azzerata. `spec_list` è la lista di spettri già
    rebinnati, `count` è un dataframe `wave`:`count` che tiene conto di quanti spettri
    contribuiscono ad una certa frequenza. `limit` indica se escludere uno spettro nel caso
    in cui sia stato ridotto male.
    """
    function stack_spectra!(stack::DataFrame, count::DataFrame, spec_list::Array{DataFrame})
        if length(spec_list) > 0
            spec = pop!(spec_list)

            for (rstack, rc, rspec) in zip(eachrow(stack), eachrow(count), eachrow(spec))
                if isnan(rspec[:flux])
                    rstack[:flux] += 0. # Lo scrivo esplicitamente per chiarezza.
                else
                    rstack[:flux] += rspec[:flux]
                    rc[:count] += 1
                end
            end

            stack_spectra!(stack, count, spec_list)
        else
            stack[!, :flux] ./= count[!, :count]
            return stack
        end
    end

    # ------------------------------ ** ------------------------------ #

    """
        sub_in_df(df::DataFrame, f, list)

    Data una lista con tanti elementi quanti le colonne, sostituiee per
    ogni colonna il corrispondente elemento della lista quando viene
    soddisfatta la condizione indicata dalla funzione f.
    """
    function sub_in_df!(df::DataFrame, f, list)
        for (i, name) in enumerate(names(df))
            for (n, val) in enumerate(df[!, name])
                if f(val)
                    df[!, name][n] = list[i]
                end
            end
        end
    end

    # ------------------------------ ** ------------------------------ #

    """
        mag_lim(l; target = .34, y = 1, bs_ = 0.05, s = 100000, tol = 1e-7, k = 3, user_max_mag = nothing, offset = 0)
    
    Determina la magnitudine limite. l è lista da plottare e di cui determinare la
    magnitudine limite, target è il valore a cui deve convergere la bisezione.
    y determina se mostrare il plot di controllo o meno, bs è il bin size dell'istogramma,
    s è il parametro che regola la smoothness dello spline, tol è la tolleranza per
    la bisezione, k l'ordine dello spline, user_max_mag è necessario nel caso si debba
    settare a mano la magnitudine massima, offset serve per dare un offset se necessario.
    """
    function mag_lim(l; target = .34, y = 1, bs_ = 0.05, s = 0, tol = 1e-7, k = 3, user_max_mag = nothing, offset = 0, mult = 2.355, user_sigma = nothing)
        si = 0.01 # sampling interval

        if iszero(s)
            s = length(l)
        end

        spline_is_zero = 0

        check_spline(x) = x >= 0 ? x : 0

        h = hist(l, bs = bs_)

        spl = Spline1D(h.bins, h.counts; s = s, k = k)
        point_list = h.bins[1]:si:h.bins[end]
        spline_points = check_spline.(spl(point_list))

        max_id = findall(spline_points .== maximum(spline_points))[1]
        max_mag = point_list[max_id]

        if !isnothing(user_max_mag)
            max_mag = user_max_mag
        end

        int_x = findall(point_list .> max_mag)
        int_approx = sum(spline_points[int_x])*si

        int = 0
        sigma = 0

        i = 0
        inc = 0.1

        while (abs((int/2)/int_approx - target) > tol && inc > 1e-7)
            int, t = quadgk(x -> spl(x), max_mag, max_mag + i)
            if (int/2)/int_approx > target
                i -= inc
                inc = inc/10
            end
            i += inc
        end
        sigma = i

        new_counts = copy(spline_points)
        for i in 1:length(point_list)
            if point_list[i] <= max_mag
                new_counts[i] = 0
            end
        end

        gauss(x, m, s) = 1/(sqrt(2*3.1415)*s) * exp(-((x - m)/s)^2/2) * (2*sum(new_counts))*si
        gauss_x = 10:si:30

        if !isnothing(user_sigma)
            sigma = user_sigma
        end

        gauss_points = gauss.(gauss_x, max_mag, sigma)
        gauss_points_reverse = maximum(gauss_points) .- gauss.(gauss_x, max_mag + mult*sigma, sigma)

        var = max_mag + mult*sigma
        y_var = 0.01

        if y == 1
        @gp "set style fill solid 0.4"
        @gp :- "set key left box"
        @gp :- h.bins h.counts "w boxes lc 'grey60' t 'Distribuzione per skymapper g'"
        @gp :- max_mag spline_points[max_id] "pt 7 ps 0.7 lc 'red' lw 0.1 t 'Massimo della distribuzione'"
        @gp :- var y_var "pt 7 ps 0.7 lc 'blue' lw 0.1 t 'Magnitudine limite'"
        @gp :- point_list spline_points "w l lc 'red' dashtype 2 t 'Spline polinomiale'"
        @gp :- gauss_x gauss_points "w l lc 'black' notitle"
        @gp :- gauss_x gauss_points_reverse "w l lc 'blue' notitle"
        save(term="pngcairo size 720, 420 fontscale 0.8", output = "output.png")
    end

        return max_mag + offset, sigma, max_mag + mult*sigma + offset
        end

    # ------------------------------ ** ------------------------------ #

    """
        adaptive_bin(id_list, z_list; min_number::Int64 = 50, min_width::Float64 = 0.2)

    Costruisce le classi per la prf tramite bin adattivo. Ritorna una tupla, il cui elemento
    1 è il dizionario delle classi, e il secondo elemento è l'intervallo di redshift per ogni
    classe.
    Ovviamente `id_list` e `z_list` sono la lista di id e redshift, `min_number` e `min_width` sono
    il minimo numero di elementi da mettere in una classe e la minima larghezza della classe
    stessa.
    """
    function adaptive_bin(id_list, z_list; min_number::Int64 = 50, min_width::Float64 = 0.2)
        d = Dict{Int64, Array{Int64}}()
        bin_width_dict = Dict{Int64, String}()

        i = 1
        z_fill = Array{Float64, 1}()
        id_fill = Array{Int64, 1}()
        last_id = length(id_list)

        for (en_id, (id, z)) in enumerate(zip(id_list, z_list))
            cond = length(id_fill) < min_number || abs(minimum(z_fill) - maximum(z_fill)) <= min_width
            if (en_id != last_id) && cond
                push!(z_fill, z)
                push!(id_fill, id)
            else
                push!(z_fill, z)
                push!(id_fill, id)
                d[i] = id_fill
                bin_width_dict[i] = string(round(minimum(z_fill), digits = 2))* " - " * string(round(maximum(z_fill), digits = 2))
                i += 1
                z_fill = Array{Float64, 1}()
                id_fill = Array{Int64, 1}()
            end
        end

        return d, bin_width_dict
    end

    # ------------------------------ ** ------------------------------ #

    """
        build_class_dataframe(class_dict)

    Costruisce il dataframe con le classi (ordinate in ordine sensato, cioè a basso class_id corrisponde
    basso z) e ritorna il detto dataframe. Durante il processo le classi vengono riordinate in modo che
    indice di classe e redshift abbiano andamento coerente.
    """
    function build_class_dataframe(class_dict)
        max_class_id = maximum(keys(class_dict)) + 1

        id_list = Array{Int64, 1}()
        class_list = Array{Int64, 1}()

        for key in keys(class_dict)
            for id in class_dict[key]
                push!(id_list, id)
                push!(class_list, max_class_id - key)
            end
        end
        return sort!(DataFrame(id = id_list, class = class_list), :id)
    end

    # ------------------------------ ** ------------------------------ #

    """
        function adaptive_class(df::DataFrame; mn = 50, mw = 0.2)

    Dato un dataframe con (:myid, :z), anche con z non ordinati ritorna un dataframe (:myid, :class).
    mn è il minimo numero di oggetti necessari per bin, mentre mw è la minima ampiezza del bin.
    Ritorna il dataframe con le classi e un dizionario che contiene l'ampiezza in z per ogni intervallo.
    """
    function adaptive_bin(df::DataFrame; mn = 50, mw = 0.2)
        sorted_df = sort(df, :z, rev = true)
        class_dict, z_class_bin = adaptive_bin(sorted_df.myid, sorted_df.z, min_number = mn, min_width = mw)
        class_df = build_class_dataframe(class_dict)

        sorted_z_class_bin = Dict{Int64, String}()
        max_class = maximum(keys(z_class_bin)) + 1
        for key in keys(z_class_bin)
            sorted_z_class_bin[max_class - key] = z_class_bin[key]
        end

        return class_df, sorted_z_class_bin
    end

    # ------------------------------ ** ------------------------------ #
end