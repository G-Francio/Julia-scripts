module fitsIO
    using FITSIO, TextParse, AstroLib, DataFrames, Dates, WeakRefStrings
    using PrettyTables, Dierckx, Gnuplot, QuadGK, MyAstroUtils, StatsBase

    import Base.write, Base.isnan

    export adaptive_bin,
        add_iid!,
        change_missing!,
        conversion_dict,
        create_dataframe,
        find_fits,
        generate_dataset!,
        generate_QSG,
        get_current_date,
        get_hdu_col_names,
        get_hdu_keys,
        get_header_values,
        get_key_values,
        image_to_table,
        is_known,
        isnan,
        login_qub,
        mag_lim,
        new_spectra_df!,
        parse_galex_results!,
        prepare_stack,
        prepMatchFiles!,
        rebin,
        remove_missing!,
        sedplot,
        sedplotMagLim,
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
        get_hdu_keys(hdu::Union{TableHDU, ASCIITableHDU})

    Ritorna tutte le chiavi per l'HDU.
    """
    function get_hdu_keys(hdu::Union{TableHDU, ASCIITableHDU})
        return keys(read_header(hdu))
    end

    # ------------------------------ ** ------------------------------ #

    """
        get_key_values(hdu::Union{TableHDU, ASCIITableHDU}, keys::Tuple)

    Ritorna una coppia di (value[key], comment[value]) per ogni key.
    Se serve una sola chiave, passa (key, ) invece che (key) per evitare casotto. Usa più
    che puoi le keywords invece degli indici che son difficili da maneggiare.
    """
    function get_key_values(hdu::Union{TableHDU, ASCIITableHDU}, keys::Tuple)
        return [read_key(hdu, key) for key in keys]
    end

    # ------------------------------ ** ------------------------------ #

    """
        get_hdu_col_names(hdu::Union{TableHDU, ASCIITableHDU})

    Ritorna i nomi delle colonne del'hdu.
    """
    function get_hdu_col_names(hdu::Union{TableHDU, ASCIITableHDU})
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
        create_dataframe_internal(hdu::Union{TableHDU, ASCIITableHDU})

    Crea e ritorna un dataframe dato un `TableHDU` da un file Fits.
    """
    function create_dataframe_internal(hdu::Union{TableHDU, ASCIITableHDU})
        col_data = []
        col_names = get_hdu_col_names(hdu)

        for col_name in col_names
            push!(col_data, read(hdu, col_name))
        end

        df = DataFrame(col_data, Symbol.(col_names))
        return df
    end

    # ------------------------------ ** ------------------------------ #

    """
        prepare_text(file_in, separator::Char = '|', file_out::AbstractString = "temp.gitignore"), cchar::Char = '#')

    Riscrive il file di input in modo che sia ancora un csv, ma senza spazi tra
    una entry e l'altra, e con `|` (pipe) come delimitatore consistente ovunque.
    """
    function prepare_text(file_in, separator::Char = '|', file_out::AbstractString = "temp.gitignore", cchar::Char = '#')

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
        create_dataframe_internal(file_in; separator::Char = '|', cchar::Char = '#', hexists::Bool = true, dict = nothing)

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

        col_names = Symbol.(col)
        data_frame = DataFrame(data_, col_names)

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
    function write(df::DataFrame, file_name::AbstractString)
        df = deepcopy(df)

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
            pretty_table(f, table; tf = tf_markdown)
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
    #function sub_in_df!(df::DataFrame, f, list)
    #    for (i, name) in enumerate(names(df))
    #        for (n, val) in enumerate(df[!, name])
    #            if f(val)
    #                df[!, name][n] = list[i]
    #            end
    #        end
    #    end
    #end

    function sub_in_df!(df::DataFrame, f, list)
        subList = copy(list)

        function fixNaNReplaceList(r, l)
            for (n, (er, el)) in enumerate(zip(r, l))
                t = eltype(er)
                if t <: Integer && isnan(el)
                    l[n] = 0
                end
            end
            return l
        end

        for row in eachrow(df)
            name  = string.(findall(f, row))
            index = indexin(name, names(df))
            subList = fixNaNReplaceList(row[index], list[index])
            row[index] = subList
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
    function mag_lim(l; target = .34, y = 1, bs_ = 0.05, s = 0, tol = 1e-7, k = 3, user_offset = 0, user_max_mag = nothing, user_sigma = nothing)
        spline_is_zero = 0
        si             = 0.01
        # ^ sampling interval

        function check_spline(x)
            if x >= 0
                return x
            else
                return 0
            end
        end

        if iszero(s)
            s = length(l)
        end

        h = hist(l, bs = bs_)

        spl           = Spline1D(h.bins, h.counts; s = s, k = k)
        point_list    = h.bins[1]:si:h.bins[end]
        spline_points = check_spline.(spl(point_list))

        max_id  = argmax(spline_points)
        max_mag = point_list[max_id]

        if !isnothing(user_max_mag)
            max_mag = user_max_mag
        end

        int_x      = findall(point_list .> max_mag)
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
        gauss_x = 10:si:25

        if !isnothing(user_sigma)
            sigma = user_sigma
        end

        gauss_points         = gauss.(gauss_x, max_mag, sigma)

        if y == 1
            @gp "set style fill solid 0.4"
            @gp :- "set key left box"
            @gp :- h.bins h.counts "w boxes lc 'grey60' t 'Distribuzione per banda'"
            @gp :- max_mag spline_points[max_id] "pt 7 ps 0.7 lc 'red' lw 0.1 t 'Massimo della distribuzione'"
            @gp :- point_list spline_points "w l lc 'red' dashtype 2 t 'Spline polinomiale'"
            @gp :- gauss_x gauss_points "w l lc 'black' notitle"
            save(term = "pngcairo size 720, 420 fontscale 0.8", output = "output.png")
        end

        return max_mag + user_offset, sigma
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
        adaptive_bin(df::DataFrame; mn = 50, mw = 0.2)

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

    """
        login_qub()
    
    Fornisce il login al DB.
    """
    function login_qub()
        cred = Array{String, 1}()
        
        open("/home/francio/.cred.auth") do file
            for ln in eachline(file)
                push!(cred, ln)
            end
        end

        user = string(cred[1])
        pass = string(cred[2])
        MyAstroUtils.DBconnect("127.0.0.1", user = user, passwd = pass, dbname = "Qubrics")
    end

    # ------------------------------ ** ------------------------------ #

    function parse_entry(x::Nothing, type)
        if type == String
            return "Nothing"
        elseif type <: Number
            return convert(type, -999)
        else 
            return x
        end
    end

    function parse_entry(x::Missing, type)
        if type == String
            return " "
        elseif type <: AbstractFloat
            return convert(type, NaN)
        elseif type <: Integer
            return convert(type, 0)
        end
    end

    function parse_entry(x::Number, type)
        x
    end

    function parse_entry(x::String, type)
        x
    end

    # ------------------------------ ** ------------------------------ #

    """
        remove_missing!(df::DataFrame)

    Rimuove i `missing` ed i `nothing` sostituendoli con una entry apposita. 
    I `nothing` diventano -999 (vario tipo, stringe incluse).
    I `missing` diventano il tipo neutro corrispondente (0, " ", ...).
    Le colonne di soli missing sono rimosse.
    """
    function remove_missing!(df::DataFrame)
        col_names = names(df)
        for name in col_names
            col_type = eltype(df[!, name])
            if col_type == Missing
                select!(df, Not(name))
            elseif col_type == Nothing
                df[!, name] = ["-999" for i in 1:length(df[!, name])]
            elseif isa(col_type, Union)
                non_missing = collect(col_type)[2]
                df[!, name] = parse_entry.(df[!, name], non_missing)
            end
        end
    end

    # ------------------------------ ** ------------------------------ #

    """
        change_missing!(df::DataFrame)

    Rimuove i `missing` ed i `nothing` sostituendoli con una entry apposita. 
    I `nothing` diventano -999 (vario tipo, stringe incluse).
    I `missing` diventano il tipo neutro corrispondente (0, " ", ...).
    Le colonne di soli `missing` o `nothing` sono trasformate in " " oppure "-999".
    """
    function change_missing!(df::DataFrame)
        col_names = names(df)
        for name in col_names
            col_type = eltype(df[!, name])
            if col_type == Missing
                df[!, name] = [" " for i in 1:length(df[!, name])]
            elseif col_type == Nothing
                df[!, name] = ["-999" for i in 1:length(df[!, name])]
            elseif isa(col_type, Union)
                non_missing = collect(col_type)[2]
                df[!, name] = parse_entry.(df[!, name], non_missing)
            end
        end
    end

    # ------------------------------ ** ------------------------------ #

    isnan(x::String) = false

    # ------------------------------ ** ------------------------------ #

    """
        generate_QSG(ms::DataFrame, star_selection::Union{Array{Int, 1}, Nothing} = nothing; n_star::Int = 5000, return_star_id::Bool = false)

    Produce i tre fits con QSO, GAL e STAR (in quest'ordine). Nel caso non sia passato una lista di indici le stelle sono scelte casualmente
    in numero  n_star (default a 5000). Il dataframe di partenza deve avere le colonne `otype`, `raj2000` e `dej2000` contemporaneamente,
    """
    function generate_QSG(ms::DataFrame, star_selection::Union{Array{Int, 1}, Nothing} = nothing; n_star::Int = 5000, return_star_id::Bool = false)
        @assert "otype" in names(ms)
        @assert "raj2000" in names(ms)
        @assert "dej2000" in names(ms)

        qso_indexes = findall(x -> x == "QSO",  ms.otype)
        gal_indexes = findall(x -> x == "GAL",  ms.otype)
        str_indexes = findall(x -> x == "STAR", ms.otype)

        if isnothing(star_selection)
            star_selection = sample(1:length(str_indexes), n_star, replace = false)
        end

        qso = ms[qso_indexes, :]
        gal = ms[gal_indexes, :]
        str = ms[str_indexes, :][star_selection, :]

        if return_star_id
            return qso, gal, str, star_selection
        else
            return qso, gal, str
        end
    end

    # ------------------------------ ** ------------------------------ #

    """
        prepMatchFiles!(qso::DataFrame, gal::DataFrame, str::DataFrame, folder::String = "/home/francio/gen")

    Produce i file pronti per essere caricati su MAST per il crossmatch. Sono da passare le coordinate (in ordine) di qso, gal
    e star. Produce tre file di testo puliti e formattati come necessario.
    """
    function prepMatchFiles!(qso::DataFrame, gal::DataFrame, str::DataFrame, folder::String = "/home/francio/gen")
        write(qso, "$folder/QSO_$(get_current_date()).fits")
        write(gal, "$folder/GAL_$(get_current_date()).fits")
        write(str, "$folder/STAR_$(get_current_date()).fits")

        qso = hcat(DataFrame(iid = collect(1:length(qso[!, 1]))), qso)
        gal = hcat(DataFrame(iid = collect(1:length(gal[!, 1]))), gal)
        str = hcat(DataFrame(iid = collect(1:length(str[!, 1]))), str)

        select!(qso, [:iid, :raj2000, :dej2000])
        select!(gal, [:iid, :raj2000, :dej2000])
        select!(str, [:iid, :raj2000, :dej2000])

        write(qso, "$folder/coord/QSO_$(get_current_date()).txt",  't')
        write(gal, "$folder/coord/GAL_$(get_current_date()).txt",  't')
        write(str, "$folder/coord/STAR_$(get_current_date()).txt", 't')

        run(`/home/francio/exe/script/prepare_xmatch.sh $folder/coord`)
    end

    # ------------------------------ ** ------------------------------ #

    """
        parse_galex_results(path::String)

    Prende i file con il risultato del crossmatch e sistema rimuovendo i doppioni. Ritorna un dataframe.
    """
    function parse_galex_results!(df::DataFrame)
        @assert "iid" in names(df)

        # specifico per Galex!
        sub_in_df!(df, x -> x == -999., fill(NaN, length(eachcol(df))))

        sort!(df, :dstArcSec)
        unique!(df, :iid)
        sort!(df, :iid)

        return df
    end

    # ------------------------------ ** ------------------------------ #

    """
        galex_mag_lim(df::DataFrame, 
        m_nuv_ais::Float64, m_fuv_ais::Float64,
        m_nuv_mis::Float64, m_fuv_mis::Float64,
        e_m_nuv_ais::Float64, e_m_fuv_ais::Float64)
        e_m_nuv_mis::Float64, e_m_fuv_mis::Float64)

    Sostituisce in Galex la magnitdini limite determinate con gli errori corrispondenti.
    MIS ha la priorità (se non viene visto nel MIS non viene visto nemmeno in AIS). Le magnitudini sono determinate
    con il solito metodo.

    Siccome un -999.0 identifica anche una misura a detector spento, dopo aver sostituo le magnidutini rimetto a NaN
    le misure che non sono state effettivamente prese.
    """
    function galex_mag_lim!(df::DataFrame, 
        m_nuv_ais::Float64, m_fuv_ais::Float64,
        m_nuv_mis::Float64, m_fuv_mis::Float64,
        e_m_nuv_ais::Float64, e_m_fuv_ais::Float64,
        e_m_nuv_mis::Float64, e_m_fuv_mis::Float64)

        for row in eachrow(df)
            if isnan(row.nuv_mag) && row.mis_nuv == 1
                row.nuv_mag = m_nuv_mis
                row.nuv_magerr = e_m_nuv_mis
            elseif isnan(row.nuv_mag) && row.mis_nuv == 0 && row.ais_nuv == 1
                row.nuv_mag = m_nuv_ais
                row.nuv_magerr = e_m_nuv_ais
            end
            
            if isnan(row.fuv_mag) && row.mis_fuv == 1
                row.fuv_mag = m_fuv_mis
                row.fuv_magerr = e_m_fuv_mis
            elseif isnan(row.fuv_mag) && row.mis_fuv == 0 && row.ais_fuv == 1
                row.fuv_mag = m_fuv_ais
                row.fuv_magerr = e_m_fuv_ais
            end

            if iszero(row.nuv_exptime)
                row.nuv_mag = NaN
                row.nuv_magerr = NaN
            end
            
            if iszero(row.fuv_exptime)
                row.fuv_mag = NaN
                row.fuv_magerr = NaN
            end
        end
    end

    # ------------------------------ ** ------------------------------ #

    """
        clean_dataset!(df::DataFrame)

    Pulisce un dataframe dalle colonne non utili per la PRF.
    """
    function clean_dataset!(df::DataFrame)
        to_remove = Array{String, 1}()
        for name in names(df)
            if name in ("RAs", "DECs", "obs", "otype", "z", "ebmv_sfd", "ext_iz", "matched")
                push!(to_remove, name)
            end
        end
        select!(df, Not(to_remove))
        return df
    end

    # ------------------------------ ** ------------------------------ #

    """
        generate_dataset!(df::DataFrame, g_data::DataFrame, tiles_ais::DataFrame, tiles_mis::DataFrame, l_galex::Array{Float64, 1}, l_mag::Array{Float64, 1}, l_err::Array{Float64, 1})

    Dato un DataFrame con coordinate, magnitudini ed errori aggiunge Galex e produce il file testuale per la PRF.
    `df` è contiene le informazioni generali, `tiles_xxx` le informazioni da Galex (non processate), `g_data` le
    informazini da Galex, `l_mag`, `l_err` le lista di magnitudini limite ed errori. 
    Le magnitudini limite vanno date secondo l'ordine di cui a seguire:

    E' necessario avere :iid nel file originale per il join con i risultati da Galex.
    NOTA BENISSIMO: l'input in questo caso è un file con oggetti tutti dello stesso tipo (es: tutti QSO o GAL).
    """
    function generate_dataset!(df::DataFrame, g_data::DataFrame, tiles_ais::DataFrame, tiles_mis::DataFrame,
        l_galex::Array{Float64, 1}, l_mag::Array{Float64, 1}, l_err::Array{Float64, 1})
        @assert "iid" in names(df)
        @assert "iid" in names(g_data)

        clean_dataset!(df)
        parse_galex_results!(g_data)

        g_footprint = DataFrame(
            ais_nuv = zeros(Int16, length(df[!, 1])),
            ais_fuv = zeros(Int16, length(df[!, 1])),
            mis_nuv = zeros(Int16, length(df[!, 1])),
            mis_fuv = zeros(Int16, length(df[!, 1]))
        )
        
        match_ais_nuv = unique!(xmatch(df.raj2000, df.dej2000, tiles_ais.ra_cent, tiles_ais.dec_cent, 2304.)[1])
        match_ais_fuv = unique!(xmatch(df.raj2000, df.dej2000, tiles_ais.ra_cent, tiles_ais.dec_cent, 2232.)[1])
        match_mis_nuv = unique!(xmatch(df.raj2000, df.dej2000, tiles_mis.ra_cent, tiles_mis.dec_cent, 2304.)[1])
        match_mis_fuv = unique!(xmatch(df.raj2000, df.dej2000, tiles_mis.ra_cent, tiles_mis.dec_cent, 2232.)[1])

        g_footprint.ais_nuv[match_ais_nuv] .= 1
        g_footprint.ais_fuv[match_ais_fuv] .= 1
        g_footprint.mis_nuv[match_mis_nuv] .= 1
        g_footprint.mis_fuv[match_mis_fuv] .= 1

        df = hcat(df, g_footprint)
        df_galex = leftjoin(df, g_data, on = :iid)
        remove_missing!(df_galex)

        galex_mag_lim!(df_galex, l_galex...)

        select!(df_galex, Not([:iid, :dstArcSec, :survey, :nuv_exptime, :fuv_exptime, :ais_nuv, :ais_fuv, :mis_nuv, :mis_fuv]))
        
        mag = select(df_galex, [:u_psf, :v_psf, :g_psf, :r_psf, :i_psf, :z_psf,
                                :phot_g_mean_mag, :phot_bp_mean_mag, :phot_rp_mean_mag, 
                                :h_m, :j_m, :k_m,
                                :w1mpro, :w2mpro, :w3mpro, :w4mpro,
                                :nuv_mag, :fuv_mag])
        err = select(df_galex, [:e_u_psf, :e_v_psf, :e_g_psf, :e_r_psf, :e_i_psf, :e_z_psf,
                                :phot_g_mean_mag_error, :phot_bp_mean_mag_error, :phot_rp_mean_mag_error, 
                                :h_cmsig, :j_cmsig, :k_cmsig,
                                :w1sigmpro, :w2sigmpro, :w3sigmpro, :w4sigmpro,
                                :nuv_magerr, :fuv_magerr])
                                
        sub_in_df!(mag, isnan, l_mag)
        sub_in_df!(err, isnan, l_err)

        return mag, err
    end

    # ------------------------------ ** ------------------------------ #

    """
        add_iid!(df::DataFrame)
    Genera la colonna :iid in un DataFrame.
    """
    function add_iid!(df::DataFrame)
        df[!, :iid] = collect(1:length(df[!, 1]))
        
    end

    # ------------------------------ ** ------------------------------ #

    function sedplot(skym1::DataFrameRow; z_est=nothing, z_spec=nothing, path::String="SED")
        function wise_vega2AB()
            [2.699, 3.339, 5.174, 6.620]
        end

        function addflux(skym::DataFrame)
            flux = DataFrame(object_id=skym.object_id)

            # http://wise2.ipac.caltech.edu/docs/release/allsky/expsup/sec4_4h.html#conv2ab
            wise_f0 = [309.540, 171.787, 31.674, 8.363]
            # Restore Vega magnitudes
            wise_dm = wise_vega2AB()
            w1mpro = skym.w1mpro .- wise_dm[1]
            w2mpro = skym.w2mpro .- wise_dm[2]
            w3mpro = skym.w3mpro .- wise_dm[3]
            w4mpro = skym.w4mpro .- wise_dm[4]

            flux[!,:w1]     = wise_f0[1] .* 10 .^( -w1mpro ./ 2.5)
            flux[!,:w2]     = wise_f0[2] .* 10 .^( -w2mpro ./ 2.5)
            flux[!,:w3]     = wise_f0[3] .* 10 .^( -w3mpro ./ 2.5)
            flux[!,:w4]     = wise_f0[4] .* 10 .^( -w4mpro ./ 2.5)
            flux[!,:err_w1] = wise_f0[1] .* 10 .^(-(w1mpro.-skym.w1sigmpro) ./ 2.5) .- flux[:,:w1]
            flux[!,:err_w2] = wise_f0[2] .* 10 .^(-(w2mpro.-skym.w2sigmpro) ./ 2.5) .- flux[:,:w2]
            flux[!,:err_w3] = wise_f0[3] .* 10 .^(-(w3mpro.-skym.w3sigmpro) ./ 2.5) .- flux[:,:w3]
            flux[!,:err_w4] = wise_f0[4] .* 10 .^(-(w4mpro.-skym.w4sigmpro) ./ 2.5) .- flux[:,:w4]

            # https://old.ipac.caltech.edu/2mass/releases/allsky/doc/sec6_4a.html
            twomass_f0 = [1594., 1024., 666.7]
            flux[!,:J]     = twomass_f0[1] .* 10 .^( -skym.j_m ./ 2.5)
            flux[!,:H]     = twomass_f0[2] .* 10 .^( -skym.h_m ./ 2.5)
            flux[!,:K]     = twomass_f0[3] .* 10 .^( -skym.k_m ./ 2.5)
            flux[!,:err_J] = twomass_f0[1] .* 10 .^(-(skym.j_m.-skym.j_cmsig) ./ 2.5) .- flux[:,:J]
            flux[!,:err_H] = twomass_f0[2] .* 10 .^(-(skym.h_m.-skym.h_cmsig) ./ 2.5) .- flux[:,:H]
            flux[!,:err_K] = twomass_f0[3] .* 10 .^(-(skym.k_m.-skym.k_cmsig) ./ 2.5) .- flux[:,:K]

            # http://skymapper.anu.edu.au/data-release/dr1/#Filters
            AB_f0 = 3631.
            flux[!,:u]     = AB_f0[1] .* 10 .^( -skym.u_psf ./ 2.5)
            flux[!,:v]     = AB_f0[1] .* 10 .^( -skym.v_psf ./ 2.5)
            flux[!,:g]     = AB_f0[1] .* 10 .^( -skym.g_psf ./ 2.5)
            flux[!,:r]     = AB_f0[1] .* 10 .^( -skym.r_psf ./ 2.5)
            flux[!,:i]     = AB_f0[1] .* 10 .^( -skym.i_psf ./ 2.5)
            flux[!,:z]     = AB_f0[1] .* 10 .^( -skym.z_psf ./ 2.5)
            flux[!,:err_u] = AB_f0[1] .* 10 .^(-(skym.u_psf.-skym.e_u_psf) ./ 2.5) .- flux[:,:u]
            flux[!,:err_v] = AB_f0[1] .* 10 .^(-(skym.v_psf.-skym.e_v_psf) ./ 2.5) .- flux[:,:v]
            flux[!,:err_g] = AB_f0[1] .* 10 .^(-(skym.g_psf.-skym.e_g_psf) ./ 2.5) .- flux[:,:g]
            flux[!,:err_r] = AB_f0[1] .* 10 .^(-(skym.r_psf.-skym.e_r_psf) ./ 2.5) .- flux[:,:r]
            flux[!,:err_i] = AB_f0[1] .* 10 .^(-(skym.i_psf.-skym.e_i_psf) ./ 2.5) .- flux[:,:i]
            flux[!,:err_z] = AB_f0[1] .* 10 .^(-(skym.z_psf.-skym.e_z_psf) ./ 2.5) .- flux[:,:z]
            
            return flux[1,:]
        end

        flux = addflux(DataFrame(skym1))
        file = "$(path)/" * string(skym1[:object_id]) * ".png"
        #isfile(file)  &&  return nothing
        w = [  3500,   3800,   5200,   6200,   7600,   9000, 12_500, 16_500, 21_500,  34_350,  46_000, 115_600, 220_800]
        m = [    :u,     :v,     :g,     :r,     :i,     :z,     :J,     :H,     :K,     :w1,     :w2,     :w3,     :w4]
        e = [:err_u, :err_v, :err_g, :err_r, :err_i, :err_z, :err_J, :err_H, :err_K, :err_w1, :err_w2, :err_w3, :err_w4]
        l = [    "u",    "v",    "g",    "r",    "i",    "z",    "J",    "H",    "K",    "W1",    "W2",    "W3",    "W4"]
        f = 3e18 ./ w
        mm = collect(flux[m]);  mm .*= 1e-10 .* f;
        ee = collect(flux[e]);  ee .*= 1e-10 .* f;
        tmp = [mm.-ee mm.+ee]
        tmp = tmp[findall(.!isnan.(tmp))]
        yr = [extrema(tmp)...] .* [0.5, 2]
        #println(flux[[:object_id, m...]])
        @gp    :sed xlab="Obs .freq. [Hz]" ylab="{/Symbol n} F_{/Symbol n} [10^{-13} erg s^{-1} cm^{-2}]" xlog=true ylog=true xr=[1e13, 3e15] yr=yr "set clip two" "set bars 0" :-
        tit = "id=" * string(flux[:object_id]) * ", i=" * string(round(skym1[:i_psf]*100)/100)
        isnothing(z_est)   ||  (tit *= ", z_{est}="  *string(round(z_est *100)/100))
        isnothing(z_spec)  ||  (tit *= ", z_{spec}=" *string(round(z_spec*100)/100))
        @gp :- :sed title=tit
        for i in 1:length(l)
            @gp :- :sed "set label '" * l[i] * "' at first " * string(f[i]) * ", first  " * string(mm[i]) * " offset character 0.5, -0.5" :-
        end
        for z in [0, 1, 2, 3, 4, 5]
            n = 3.e18 / 1216 / (1+z)
            @gp :- :sed n .* [1,1] [0.1, 1000]  "w l notit dt 3 lc rgb 'grey'"
            @gp :- :sed "set label 'Ly{/Symbol a} at z=" * string(z) * "' textcolor 'grey' at first " * string(n) * ", screen 0.14 rotate by 90 offset character -1,0" :-
            n = 3.e18 / 1e4 / (1+z)
            @gp :- :sed n .* [1,1] [0.1, 1000]  "w l notit dt 3 lc rgb 'grey'"
            @gp :- :sed "set label '1{/Symbol m}m at z=" * string(z) * "' textcolor 'grey' at first " * string(n) * ", screen 0.14 rotate by 90 offset character -1,0" :-
        end
        @gp :- :sed f mm ee "w yerr notit pt 7 lc rgb 'red'"
        save(:sed, term="pngcairo enhanced size 750,500", output=file)
        save(:sed, file*".gp")
    end

    # ------------------------------ ** ------------------------------ #

    function sedplotMagLim(_skym1::DataFrameRow; z_est=nothing, z_spec=nothing, path::String="SED")
        function wise_vega2AB()
            [2.699, 3.339, 5.174, 6.620]
        end
    
        function addflux(skym::DataFrame)
            _flux = DataFrame(object_id=skym.object_id)
    
            # http://wise2.ipac.caltech.edu/docs/release/allsky/expsup/sec4_4h.html#conv2ab
            wise_f0 = [309.540, 171.787, 31.674, 8.363]
            # Restore Vega magnitudes
            wise_dm = wise_vega2AB()
            w1mpro = skym.w1mpro .- wise_dm[1]
            w2mpro = skym.w2mpro .- wise_dm[2]
            w3mpro = skym.w3mpro .- wise_dm[3]
            w4mpro = skym.w4mpro .- wise_dm[4]
    
            _flux[!,:w1]     = wise_f0[1] .* 10 .^( -w1mpro ./ 2.5)
            _flux[!,:w2]     = wise_f0[2] .* 10 .^( -w2mpro ./ 2.5)
            _flux[!,:w3]     = wise_f0[3] .* 10 .^( -w3mpro ./ 2.5)
            _flux[!,:w4]     = wise_f0[4] .* 10 .^( -w4mpro ./ 2.5)
            _flux[!,:err_w1] = wise_f0[1] .* 10 .^(-(w1mpro.-skym.w1sigmpro) ./ 2.5) .- _flux[:,:w1]
            _flux[!,:err_w2] = wise_f0[2] .* 10 .^(-(w2mpro.-skym.w2sigmpro) ./ 2.5) .- _flux[:,:w2]
            _flux[!,:err_w3] = wise_f0[3] .* 10 .^(-(w3mpro.-skym.w3sigmpro) ./ 2.5) .- _flux[:,:w3]
            _flux[!,:err_w4] = wise_f0[4] .* 10 .^(-(w4mpro.-skym.w4sigmpro) ./ 2.5) .- _flux[:,:w4]
    
            # https://old.ipac.caltech.edu/2mass/releases/allsky/doc/sec6_4a.html
            twomass_f0 = [1594., 1024., 666.7]
            _flux[!,:J]     = twomass_f0[1] .* 10 .^( -skym.j_m ./ 2.5)
            _flux[!,:H]     = twomass_f0[2] .* 10 .^( -skym.h_m ./ 2.5)
            _flux[!,:K]     = twomass_f0[3] .* 10 .^( -skym.k_m ./ 2.5)
            _flux[!,:err_J] = twomass_f0[1] .* 10 .^(-(skym.j_m.-skym.j_cmsig) ./ 2.5) .- _flux[:,:J]
            _flux[!,:err_H] = twomass_f0[2] .* 10 .^(-(skym.h_m.-skym.h_cmsig) ./ 2.5) .- _flux[:,:H]
            _flux[!,:err_K] = twomass_f0[3] .* 10 .^(-(skym.k_m.-skym.k_cmsig) ./ 2.5) .- _flux[:,:K]
    
            # http://skymapper.anu.edu.au/data-release/dr1/#Filters
            AB_f0 = 3631.
            _flux[!,:u]     = AB_f0[1] .* 10 .^( -skym.u_psf ./ 2.5)
            _flux[!,:v]     = AB_f0[1] .* 10 .^( -skym.v_psf ./ 2.5)
            _flux[!,:g]     = AB_f0[1] .* 10 .^( -skym.g_psf ./ 2.5)
            _flux[!,:r]     = AB_f0[1] .* 10 .^( -skym.r_psf ./ 2.5)
            _flux[!,:i]     = AB_f0[1] .* 10 .^( -skym.i_psf ./ 2.5)
            _flux[!,:z]     = AB_f0[1] .* 10 .^( -skym.z_psf ./ 2.5)
            _flux[!,:err_u] = AB_f0[1] .* 10 .^(-(skym.u_psf.-skym.e_u_psf) ./ 2.5) .- _flux[:,:u]
            _flux[!,:err_v] = AB_f0[1] .* 10 .^(-(skym.v_psf.-skym.e_v_psf) ./ 2.5) .- _flux[:,:v]
            _flux[!,:err_g] = AB_f0[1] .* 10 .^(-(skym.g_psf.-skym.e_g_psf) ./ 2.5) .- _flux[:,:g]
            _flux[!,:err_r] = AB_f0[1] .* 10 .^(-(skym.r_psf.-skym.e_r_psf) ./ 2.5) .- _flux[:,:r]
            _flux[!,:err_i] = AB_f0[1] .* 10 .^(-(skym.i_psf.-skym.e_i_psf) ./ 2.5) .- _flux[:,:i]
            _flux[!,:err_z] = AB_f0[1] .* 10 .^(-(skym.z_psf.-skym.e_z_psf) ./ 2.5) .- _flux[:,:z]
            
            return _flux[1,:]
        end
    
        m = [:u_psf,   :v_psf,   :g_psf,   :r_psf,   :i_psf,   :z_psf,   :h_m,     :j_m,     :k_m,     :w1mpro,    :w2mpro,    :w3mpro,    :w4mpro]
        e = [:e_u_psf, :e_v_psf, :e_g_psf, :e_r_psf, :e_i_psf, :e_z_psf, :h_cmsig, :j_cmsig, :k_cmsig, :w1sigmpro, :w2sigmpro, :w3sigmpro, :w4sigmpro]
    
        skym1  = deepcopy(_skym1)
        skym1_ = deepcopy(_skym1)
    
        for (_m, _e) in zip(m, e)
            if _skym1[_e] > 5
                skym1[_e]  = NaN
                skym1_[_e] = NaN
            else
                skym1_[_m] = NaN
                skym1_[_e] = NaN
            end
        end
    
        flux       = addflux(DataFrame(skym1))
        fluxLimMag = addflux(DataFrame(skym1_))
    
        file = "$(path)/" * string(skym1[:object_id]) * ".png"
        #isfile(file)  &&  return nothing
        w = [  3500,   3800,   5200,   6200,   7600,   9000, 12_500, 16_500, 21_500,  34_350,  46_000, 115_600, 220_800]
        m = [    :u,     :v,     :g,     :r,     :i,     :z,     :J,     :H,     :K,     :w1,     :w2,     :w3,     :w4]
        e = [:err_u, :err_v, :err_g, :err_r, :err_i, :err_z, :err_J, :err_H, :err_K, :err_w1, :err_w2, :err_w3, :err_w4]
        l = [    "u",    "v",    "g",    "r",    "i",    "z",    "J",    "H",    "K",    "W1",    "W2",    "W3",    "W4"]
        f = 3e18 ./ w
    
        mm = collect(flux[m]);  mm .*= 1e-10 .* f;
        ee = collect(flux[e]);  ee .*= 1e-10 .* f;
    
        mm_ = collect(fluxLimMag[m]);  mm_ .*= 1e-10 .* f;
        ee_ = collect(fluxLimMag[e]);  ee_ .*= 1e-10 .* f;
    
        tmp = [mm.-ee mm.+ee]
        tmp = tmp[findall(.!isnan.(tmp))]
        yr = [extrema(tmp)...] .* [0.5, 2]
        #println(flux[[:object_id, m...]])
        @gp    :sed xlab="Obs .freq. [Hz]" ylab="{/Symbol n} F_{/Symbol n} [10^{-13} erg s^{-1} cm^{-2}]" xlog=true ylog=true xr=[1e13, 3e15] yr=yr "set clip two" "set bars 0" :-
        tit = "id=" * string(flux[:object_id]) * ", i=" * string(round(skym1[:i_psf]*100)/100)
        isnothing(z_est)   ||  (tit *= ", z_{est}="  *string(round(z_est *100)/100))
        isnothing(z_spec)  ||  (tit *= ", z_{spec}=" *string(round(z_spec*100)/100))
        @gp :- :sed title=tit

        for i in 1:length(l)
            @gp :- :sed "set label '" * l[i] * "' at first " * string(f[i]) * ", first  " * string(mm[i]) * " offset character 0.5, -0.5" :-
        end
        for z in [0, 1, 2, 3, 4, 5]
            n = 3.e18 / 1216 / (1+z)
            @gp :- :sed n .* [1,1] [0.1, 1000]  "w l notit dt 3 lc rgb 'grey'"
            @gp :- :sed "set label 'Ly{/Symbol a} at z=" * string(z) * "' textcolor 'grey' at first " * string(n) * ", screen 0.14 rotate by 90 offset character -1,0" :-
            n = 3.e18 / 1e4 / (1+z)
            @gp :- :sed n .* [1,1] [0.1, 1000]  "w l notit dt 3 lc rgb 'grey'"
            @gp :- :sed "set label '1{/Symbol m}m at z=" * string(z) * "' textcolor 'grey' at first " * string(n) * ", screen 0.14 rotate by 90 offset character -1,0" :-
        end
        @gp :- :sed f mm  ee  "w yerr notit pt 7 lc rgb 'red'"
        @gp :- :sed f mm_ ee_ "w yerr notit pt 7 lc rgb 'blue'"
        save(:sed, term="pngcairo enhanced size 750,500", output=file)
        save(:sed, file*".gp")
    end
end
