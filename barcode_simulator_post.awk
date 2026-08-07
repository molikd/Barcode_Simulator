#!/usr/bin/awk -f

# Barcode Simulator Post-processing
#
# Dependency-free replacement for Barcode_Simulator_Post.R and
# Barcode_Simulator_Post_Single.R. Numeric results are emitted as tidy TSV
# instead of being rendered directly as R plots.

BEGIN {
    stderr = "/dev/stderr"
    set_defaults()
    parse_arguments()

    if (show_help) {
        usage()
        exit 0
    }

    validate_options()
    prepare_output()
    if (metadata_file != "")
        load_metadata(metadata_file)

    if (mode == "distances")
        run_distances(0)
    else if (mode == "all-distances")
        run_distances(1)
    else if (mode == "mantel")
        run_mantel()
    else
        fail("unknown mode: " mode)

    close_output()
    exit 0
}

function set_defaults() {
    mode = ""
    output_file = ""
    metadata_file = "meta.data.file.txt"
    jaccard_mode = "quantitative"
    mash_values = "distance"
    show_help = 0
}

function parse_arguments(    i, argument, option, value, equals_at) {
    for (i = 1; i < ARGC; i++) {
        argument = ARGV[i]
        if (argument == "--")
            continue

        if (argument == "-h" || argument == "--help") {
            show_help = 1
            continue
        }

        if (substr(argument, 1, 2) == "--") {
            equals_at = index(argument, "=")
            if (equals_at) {
                option = substr(argument, 1, equals_at - 1)
                value = substr(argument, equals_at + 1)
            } else {
                option = argument
                if (i + 1 >= ARGC)
                    fail("missing value for " option)
                value = ARGV[++i]
            }
            set_option(option, value)
            continue
        }

        if (mode == "") {
            mode = normalize_mode(argument)
            continue
        }

        input_file[++input_count] = argument
    }
    ARGC = 1
}

function set_option(option, value) {
    if (option == "--output")
        output_file = value
    else if (option == "--meta")
        metadata_file = (value == "none" ? "" : value)
    else if (option == "--jaccard")
        jaccard_mode = tolower(value)
    else if (option == "--mash-values")
        mash_values = tolower(value)
    else
        fail("unknown option: " option)
}

function normalize_mode(value) {
    value = tolower(value)
    gsub(/_/, "-", value)
    if (value == "all-distances" || value == "distances" || value == "mantel")
        return value
    return value
}

function validate_options() {
    if (mode == "")
        fail("an analysis mode is required")
    if (jaccard_mode != "quantitative" && jaccard_mode != "binary")
        fail("--jaccard must be 'quantitative' or 'binary'")
    if (mash_values != "distance" && mash_values != "similarity")
        fail("--mash-values must be 'distance' or 'similarity'")

    if (output_file == "") {
        if (mode == "mantel")
            output_file = "mantel.tsv"
        else if (mode == "distances")
            output_file = "distances.tsv"
        else
            output_file = "all_distances.tsv"
    }
}

function prepare_output() {
    if (output_file != "-") {
        printf "%s", "" > output_file
        close(output_file)
    }
}

function run_distances(all_points,    i, otu_file, mash_file, dataset, type_name, samples_count, mean_otu, mean_mash, correlation, a, b, key) {
    if (input_count == 0)
        discover_files("otu")
    if (input_count == 0)
        fail("no OTU tables found; pass one or more *_otu_table* files")

    sort_inputs()
    if (all_points)
        emit("dataset\ttype\tsample_a\tsample_b\totu_jaccard\tmash_distance")
    else
        emit("dataset\ttype\tsamples\tpairs\tmean_otu_jaccard\tmean_mash_distance\tmantel_r")

    for (i = 1; i <= input_count; i++) {
        otu_file = input_file[i]
        if (index(otu_file, "_otu_table") == 0)
            fail("distances modes expect an OTU table, received: " otu_file)
        mash_file = paired_mash_path(otu_file)

        ensure_matrix(otu_file)
        ensure_matrix(mash_file)
        compare_matrices(otu_file, mash_file)

        dataset = dataset_name(otu_file)
        type_name = lookup_type(otu_file)
        samples_count = common_sample_count(otu_file, mash_file)

        if (!all_points) {
            mean_otu = comparison_n ? comparison_sum_x / comparison_n : "NA"
            mean_mash = comparison_n ? comparison_sum_y / comparison_n : "NA"
            correlation = correlation_value()
            emit(dataset "\t" type_name "\t" samples_count "\t" comparison_n \
                 "\t" format_number(mean_otu) "\t" format_number(mean_mash) \
                 "\t" format_number(correlation))
        } else {
            clear_array(common_samples)
            collect_common_samples(otu_file, mash_file)
            sort_named_values(common_samples, sorted_sample)
            for (a = 1; a <= sorted_count; a++) {
                for (b = a + 1; b <= sorted_count; b++) {
                    key = otu_file SUBSEP sorted_sample[a] SUBSEP sorted_sample[b]
                    if (!(key in matrix_value))
                        continue
                    key = mash_file SUBSEP sorted_sample[a] SUBSEP sorted_sample[b]
                    if (!(key in matrix_value))
                        continue
                    emit(dataset "\t" type_name "\t" sorted_sample[a] "\t" sorted_sample[b] \
                         "\t" format_number(matrix_value[otu_file SUBSEP sorted_sample[a] SUBSEP sorted_sample[b]]) \
                         "\t" format_number(matrix_value[mash_file SUBSEP sorted_sample[a] SUBSEP sorted_sample[b]]))
                }
            }
            clear_array(sorted_sample)
        }
    }
}

function run_mantel(    i, j, file_a, file_b, samples_count) {
    if (input_count == 0)
        discover_files("matrices")
    if (input_count == 0)
        fail("no distance inputs found; pass *_otu_table* and/or *_mash_dists* files")

    sort_inputs()
    for (i = 1; i <= input_count; i++)
        ensure_matrix(input_file[i])

    emit("file_a\tfile_b\ttype_a\ttype_b\tsamples\tpairs\tpearson_r")
    for (i = 1; i <= input_count; i++) {
        for (j = i; j <= input_count; j++) {
            file_a = input_file[i]
            file_b = input_file[j]
            compare_matrices(file_a, file_b)
            samples_count = common_sample_count(file_a, file_b)
            emit(file_a "\t" file_b "\t" matrix_type[file_a] "\t" matrix_type[file_b] \
                 "\t" samples_count "\t" comparison_n "\t" format_number(correlation_value()))
        }
    }
}

function ensure_matrix(path) {
    if (path in matrix_loaded)
        return
    if (index(path, "_otu_table") > 0)
        load_otu_matrix(path)
    else if (index(path, "_mash_dists") > 0)
        load_mash_matrix(path)
    else
        fail("cannot infer input type from filename: " path)
    matrix_loaded[path] = 1
}

function load_otu_matrix(path,    status, line, header_seen, fields, field_count, sample_count, row, i, j, value_i, value_j, pair, numerator, denominator, sample_a, sample_b, temporary_sample, otu_value, seen_sample) {
    header_seen = 0
    sample_count = 0
    row = 0
    clear_array(otu_sample)
    clear_array(pair_numerator)
    clear_array(pair_denominator)

    while ((status = getline line < path) > 0) {
        sub(/\r$/, "", line)
        if (line ~ /^[[:space:]]*$/)
            continue

        if (!header_seen) {
            if (line ~ /^#[[:space:]]*Constructed/)
                continue
            if (line ~ /^#/ && line !~ /^#[[:space:]]*OTU[[:space:]]*ID/)
                continue

            field_count = split_fields(line, fields)
            if (field_count < 3)
                fail("OTU table needs at least two sample columns: " path)
            for (i = 2; i <= field_count; i++) {
                sample_a = normalize_sample(fields[i])
                if (sample_a == "")
                    fail("empty sample name in OTU table: " path)
                if (sample_a in seen_sample)
                    fail("duplicate normalized sample name '" sample_a "' in " path)
                seen_sample[sample_a] = 1
                otu_sample[++sample_count] = sample_a
                matrix_sample[path SUBSEP sample_a] = 1
            }
            header_seen = 1
            continue
        }

        if (line ~ /^#/)
            continue
        field_count = split_fields(line, fields)
        if (field_count < sample_count + 1)
            fail("short OTU-table row in " path ": " line)
        row++

        for (i = 1; i <= sample_count; i++) {
            if (!is_number(fields[i + 1]) || fields[i + 1] + 0 < 0)
                fail("OTU counts must be non-negative numbers in " path ": " fields[i + 1])
            otu_value[i] = fields[i + 1] + 0
        }

        for (i = 1; i <= sample_count; i++) {
            value_i = otu_value[i]
            for (j = i + 1; j <= sample_count; j++) {
                value_j = otu_value[j]
                pair = i SUBSEP j
                if (jaccard_mode == "binary") {
                    if (value_i > 0 && value_j > 0)
                        pair_numerator[pair]++
                    if (value_i > 0 || value_j > 0)
                        pair_denominator[pair]++
                } else {
                    pair_numerator[pair] += (value_i < value_j ? value_i : value_j)
                    pair_denominator[pair] += (value_i > value_j ? value_i : value_j)
                }
            }
        }
    }
    close(path)

    if (status < 0)
        fail("could not read OTU table: " path)
    if (!header_seen || row == 0)
        fail("OTU table has no data rows: " path)

    for (i = 1; i <= sample_count; i++) {
        for (j = i + 1; j <= sample_count; j++) {
            pair = i SUBSEP j
            denominator = pair_denominator[pair]
            numerator = pair_numerator[pair]
            sample_a = otu_sample[i]
            sample_b = otu_sample[j]
            if (sample_a > sample_b) {
                temporary_sample = sample_a
                sample_a = sample_b
                sample_b = temporary_sample
            }
            matrix_value[path SUBSEP sample_a SUBSEP sample_b] = \
                denominator == 0 ? 0 : 1 - numerator / denominator
        }
    }
    matrix_type[path] = "otu"
    matrix_sample_count[path] = sample_count

    clear_array(otu_sample)
    clear_array(pair_numerator)
    clear_array(pair_denominator)
}

function load_mash_matrix(path,    status, line, fields, field_count, sample_a, sample_b, value, temporary_sample, row) {
    row = 0
    while ((status = getline line < path) > 0) {
        sub(/\r$/, "", line)
        if (line ~ /^[[:space:]]*$/ || line ~ /^#/)
            continue
        field_count = split_fields(line, fields)
        if (field_count < 3)
            fail("Mash row needs at least three fields in " path ": " line)
        if (!is_number(fields[3]))
            fail("non-numeric Mash value in " path ": " fields[3])

        sample_a = normalize_sample(fields[1])
        sample_b = normalize_sample(fields[2])
        value = fields[3] + 0
        if (mash_values == "similarity")
            value = 1 - value
        if (value < 0 || value > 1)
            fail("Mash distance outside [0,1] in " path ": " value)

        matrix_sample[path SUBSEP sample_a] = 1
        matrix_sample[path SUBSEP sample_b] = 1
        if (sample_a == sample_b)
            continue
        if (sample_a > sample_b) {
            temporary_sample = sample_a
            sample_a = sample_b
            sample_b = temporary_sample
        }
        matrix_value[path SUBSEP sample_a SUBSEP sample_b] = value
        row++
    }
    close(path)

    if (status < 0)
        fail("could not read Mash distances: " path)
    if (row == 0)
        fail("Mash file has no pairwise distances: " path)

    matrix_type[path] = "mash"
    matrix_sample_count[path] = count_matrix_samples(path)
}

function compare_matrices(file_a, file_b,    a, b, key_a, key_b, x, y) {
    comparison_n = 0
    comparison_sum_x = 0
    comparison_sum_y = 0
    comparison_sum_xx = 0
    comparison_sum_yy = 0
    comparison_sum_xy = 0

    clear_array(common_samples)
    collect_common_samples(file_a, file_b)
    sort_named_values(common_samples, sorted_sample)
    for (a = 1; a <= sorted_count; a++) {
        for (b = a + 1; b <= sorted_count; b++) {
            key_a = file_a SUBSEP sorted_sample[a] SUBSEP sorted_sample[b]
            key_b = file_b SUBSEP sorted_sample[a] SUBSEP sorted_sample[b]
            if (!(key_a in matrix_value) || !(key_b in matrix_value))
                continue
            x = matrix_value[key_a]
            y = matrix_value[key_b]
            comparison_n++
            comparison_sum_x += x
            comparison_sum_y += y
            comparison_sum_xx += x * x
            comparison_sum_yy += y * y
            comparison_sum_xy += x * y
        }
    }
    clear_array(sorted_sample)
}

function correlation_value(    denominator_x, denominator_y) {
    if (comparison_n < 2)
        return "NA"
    denominator_x = comparison_n * comparison_sum_xx - comparison_sum_x * comparison_sum_x
    denominator_y = comparison_n * comparison_sum_yy - comparison_sum_y * comparison_sum_y
    if (denominator_x <= 0 || denominator_y <= 0)
        return "NA"
    return (comparison_n * comparison_sum_xy - comparison_sum_x * comparison_sum_y) / \
           sqrt(denominator_x * denominator_y)
}

function collect_common_samples(file_a, file_b,    key, pieces, sample) {
    clear_array(common_samples)
    for (key in matrix_sample) {
        split(key, pieces, SUBSEP)
        if (pieces[1] != file_a)
            continue
        sample = pieces[2]
        if ((file_b SUBSEP sample) in matrix_sample)
            common_samples[sample] = 1
    }
}

function common_sample_count(file_a, file_b,    sample, count) {
    collect_common_samples(file_a, file_b)
    count = 0
    for (sample in common_samples)
        count++
    return count
}

function count_matrix_samples(path,    key, pieces, count) {
    count = 0
    for (key in matrix_sample) {
        split(key, pieces, SUBSEP)
        if (pieces[1] == path)
            count++
    }
    return count
}

function load_metadata(path,    status, line, fields, field_count, header_seen, file_column, type_column, i, file_name, type_name) {
    header_seen = 0
    while ((status = getline line < path) > 0) {
        sub(/\r$/, "", line)
        if (line ~ /^[[:space:]]*$/)
            continue
        field_count = split_fields(line, fields)
        if (!header_seen) {
            for (i = 1; i <= field_count; i++) {
                if (tolower(fields[i]) == "file")
                    file_column = i
                if (tolower(fields[i]) == "type" || tolower(fields[i]) == "treatment")
                    type_column = i
            }
            if (!file_column || !type_column)
                fail("metadata needs File and Type/Treatment columns: " path)
            header_seen = 1
            continue
        }
        file_name = fields[file_column]
        type_name = fields[type_column]
        metadata_type[file_name] = type_name
        metadata_type[base_name(file_name)] = type_name
    }
    close(path)
    if (status < 0) {
        # The historical R scripts assumed this file, but post-processing is
        # still useful without it. Missing default metadata is non-fatal.
        if (path != "meta.data.file.txt")
            fail("could not read metadata: " path)
        metadata_file = ""
    }
}

function lookup_type(path,    name, dataset) {
    if (path in metadata_type)
        return metadata_type[path]
    name = base_name(path)
    if (name in metadata_type)
        return metadata_type[name]
    dataset = dataset_name(path)
    if (dataset in metadata_type)
        return metadata_type[dataset]
    return "unknown"
}

function paired_mash_path(otu_path,    mash_path) {
    mash_path = otu_path
    sub(/_otu_table[^\/]*$/, "_mash_dists.txt", mash_path)
    return mash_path
}

function dataset_name(path,    name) {
    name = base_name(path)
    sub(/_otu_table.*$/, "", name)
    sub(/_mash_dists.*$/, "", name)
    return name
}

function normalize_sample(value) {
    gsub(/^[[:space:]]+|[[:space:]]+$/, "", value)
    value = base_name(value)
    sub(/\.msh$/, "", value)
    sub(/\.fasta$/, "", value)
    sub(/\.fa$/, "", value)
    sub(/\.fastq$/, "", value)
    sub(/\.fq$/, "", value)
    return value
}

function base_name(path) {
    sub(/^.*\//, "", path)
    return path
}

function split_fields(line, fields,    count) {
    clear_array(fields)
    if (index(line, "\t") > 0)
        count = split(line, fields, "\t")
    else
        count = split(line, fields, /[[:space:]]+/)
    return count
}

function is_number(value) {
    return value ~ /^[-+]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$/
}

function discover_files(kind,    command, path, status) {
    if (kind == "otu")
        command = "find . ! -name . -prune -type f -name '*_otu_table*' -print"
    else
        command = "find . ! -name . -prune -type f \\( -name '*_otu_table*' -o -name '*_mash_dists*' \\) -print"

    while ((status = command | getline path) > 0) {
        sub(/^\.\//, "", path)
        input_file[++input_count] = path
    }
    close(command)
    if (status < 0)
        fail("could not discover input files")
}

function sort_inputs(    i, j, value) {
    for (i = 2; i <= input_count; i++) {
        value = input_file[i]
        j = i - 1
        while (j >= 1 && input_file[j] > value) {
            input_file[j + 1] = input_file[j]
            j--
        }
        input_file[j + 1] = value
    }
}

function sort_named_values(source, destination,    key, i, j, value) {
    clear_array(destination)
    sorted_count = 0
    for (key in source)
        destination[++sorted_count] = key
    for (i = 2; i <= sorted_count; i++) {
        value = destination[i]
        j = i - 1
        while (j >= 1 && destination[j] > value) {
            destination[j + 1] = destination[j]
            j--
        }
        destination[j + 1] = value
    }
}

function format_number(value) {
    if (value == "NA")
        return "NA"
    if (value > -0.0000000000005 && value < 0.0000000000005)
        value = 0
    return sprintf("%.12g", value)
}

function emit(line) {
    if (output_file == "-")
        print line
    else
        print line >> output_file
}

function close_output() {
    if (output_file != "-")
        close(output_file)
}

function clear_array(array,    key) {
    for (key in array)
        delete array[key]
}

function usage() {
    print "Usage: Barcode_Simulator_Post MODE [options] [files...]"
    print ""
    print "Modes:"
    print "  distances       summarize paired OTU and Mash distance files"
    print "  all-distances   emit the former scatterplot's paired point data"
    print "  mantel          correlate every unique pair of supplied matrices"
    print ""
    print "Options:"
    print "  --output FILE             output TSV ('-' for stdout)"
    print "  --meta FILE               File/Type metadata table ('none' to disable)"
    print "  --jaccard MODE            quantitative (R-compatible) or binary"
    print "  --mash-values MODE        distance (standard Mash) or similarity"
    print "  -h, --help                show this help"
    print ""
    print "If files are omitted, *_otu_table* and *_mash_dists* files are discovered"
    print "in the current directory. Distances modes derive each Mash filename by"
    print "replacing _otu_table... with _mash_dists.txt."
}

function fail(message) {
    print "Barcode_Simulator_Post: " message > stderr
    print "Try 'Barcode_Simulator_Post --help' for usage." > stderr
    exit 2
}
