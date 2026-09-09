#!/usr/bin/awk -f

# Barcode Distance
#
# Calculate distance matrices directly from simulated barcode FASTA files.
# This is intentionally self-contained: no OTU picker, Mash sketch, R package,
# plotting library, or other bioinformatics program is required.

BEGIN {
    stderr = "/dev/stderr"
    initialize_metrics()
    parse_arguments()

    if (show_help) {
        usage()
        exit 0
    }

    set_output_defaults()
    validate_options()

    if (manifest_file != "")
        load_manifest(manifest_file)
    else
        load_positional_dataset()

    validate_design()
    prepare_outputs()
    analyze_design()
    close(output_file)
    close(summary_file)
    if (write_graph)
        write_svg_graph(graph_file)
    report()
    if (run_metadata != "")
        write_run_metadata()
    exit 0
}

function initialize_metrics() {
    metric_count = 5
    metric_name[1] = "jaccard"
    metric_name[2] = "ruzicka"
    metric_name[3] = "bray-curtis"
    metric_name[4] = "cosine"
    metric_name[5] = "hellinger"

    write_matrices = 1
    write_graph = 1
    manifest_file = ""
    output_file = ""
    summary_file = ""
    matrix_prefix = ""
    graph_file = ""
    dataset_label = "dataset"
    show_help = 0
    loglevel = 2
    run_metadata = ""
}

function parse_arguments(    argument_index, argument, option, value, equals_at) {
    for (argument_index = 1; argument_index < ARGC; argument_index++) {
        argument = ARGV[argument_index]
        if (argument == "--")
            continue

        if (argument == "-h" || argument == "--help") {
            show_help = 1
            continue
        }
        if (argument == "--no-matrices") {
            write_matrices = 0
            continue
        }
        if (argument == "--no-graph") {
            write_graph = 0
            continue
        }

        if (substr(argument, 1, 2) == "--") {
            equals_at = index(argument, "=")
            if (equals_at) {
                option = substr(argument, 1, equals_at - 1)
                value = substr(argument, equals_at + 1)
            } else {
                option = argument
                if (argument_index + 1 >= ARGC)
                    fail("missing value for " option)
                value = ARGV[++argument_index]
            }
            set_option(option, value)
            continue
        }

        positional_file[++positional_count] = argument
    }
    ARGC = 1
}

function set_option(option, value) {
    if (option == "--manifest")
        manifest_file = value
    else if (option == "--output")
        output_file = value
    else if (option == "--summary")
        summary_file = value
    else if (option == "--matrix-prefix")
        matrix_prefix = value
    else if (option == "--graph")
        graph_file = value
    else if (option == "--dataset")
        dataset_label = value
    else if (option == "--loglevel")
        parse_loglevel(value)
    else if (option == "--run-metadata")
        run_metadata = value
    else
        fail("unknown option: " option)
}

function set_output_defaults(    name) {
    if (manifest_file != "") {
        name = base_name(manifest_file)
        sub(/[.]gz$/, "", name)
        sub(/-manifest[.]tsv$/, "", name)
        sub(/[.]tsv$/, "", name)
    } else
        name = safe_name(dataset_label)

    if (name == "")
        name = "barcode"
    analysis_prefix = name

    if (output_file == "")
        output_file = analysis_prefix "-distances.tsv"
    if (summary_file == "")
        summary_file = analysis_prefix "-distance-summary.tsv"
    if (matrix_prefix == "")
        matrix_prefix = analysis_prefix
    if (graph_file == "")
        graph_file = analysis_prefix "-distance-summary.svg"
}

function validate_options() {
    if (manifest_file != "" && positional_count)
        fail("--manifest cannot be combined with positional FASTA files")
    if (manifest_file == "" && positional_count == 0)
        fail("pass --manifest FILE or at least two FASTA files")
    if (dataset_label == "")
        fail("--dataset cannot be empty")
    if (index(dataset_label, SUBSEP))
        fail("--dataset contains an unsupported control character")
    if (matrix_prefix == "" && write_matrices)
        fail("--matrix-prefix cannot be empty")
    if (path_key(output_file) == path_key(summary_file))
        fail("pairwise and summary outputs must be different files")
    if (write_graph && (path_key(graph_file) == path_key(output_file) || \
        path_key(graph_file) == path_key(summary_file)))
        fail("the SVG graph must not overwrite a TSV output")
    if (manifest_file != "" && is_output_path(manifest_file))
        fail("an output path would overwrite the input manifest")
}

function prepare_outputs() {
    truncate_file(output_file)
    print "dataset\ttreatment\treplicate\tsample_a\tsample_b\treads_a\treads_b\tshared_sequences\tunion_sequences\tjaccard\truzicka\tbray_curtis\tcosine\thellinger" >> output_file
    close(output_file)

    truncate_file(summary_file)
    print "dataset\ttreatment\treplicate\tsamples\tpairs\tmetric\tmean_distance\tsd_distance" >> summary_file
    close(summary_file)
}

function parse_loglevel(value) {
    value = tolower(value)
    if (value == "error")
        loglevel = 0
    else if (value == "warning" || value == "warn")
        loglevel = 1
    else if (value == "info")
        loglevel = 2
    else if (value == "debug")
        loglevel = 3
    else
        fail("--loglevel must be error, warning, info, or debug (received '" value "')")
}

function log_message(level, text) {
    if (level == "debug" && loglevel < 3)
        return
    if (level == "info" && loglevel < 2)
        return
    if (level == "warning" && loglevel < 1)
        return
    print text > stderr
}

function shell_quote(text) {
    gsub(/'/, "'\\''", text)
    return "'" text "'"
}

function json_escape(text) {
    gsub(/\\/, "\\\\", text)
    gsub(/"/, "\\\"", text)
    gsub(/\n/, "\\n", text)
    gsub(/\r/, "\\r", text)
    gsub(/\t/, "\\t", text)
    return text
}

function run_stamp(    command, status, line) {
    command = "date -u +%Y-%m-%dT%H:%M:%SZ"
    if ((status = (command | getline line)) > 0) {
        close(command)
        return line
    }
    close(command)
    return "unknown"
}

function run_platform(    command, status, line) {
    command = "uname -srm"
    if ((status = (command | getline line)) > 0) {
        close(command)
        return line
    }
    close(command)
    return "unknown"
}

# Transparent input reader: plain files use getline, *.gz files stream
# through gzip -dc so compressed inputs just work.
function open_read(path) {
    read_is_pipe = (path ~ /[.]gz$/)
    if (read_is_pipe)
        read_source = "gzip -dc -- " shell_quote(path)
    else
        read_source = path
}

function read_next(    status) {
    if (read_is_pipe)
        status = (read_source | getline read_line)
    else
        status = (getline read_line < read_source)
    return status
}

function read_close() {
    close(read_source)
    read_is_pipe = 0
}

function write_run_metadata(    json) {
    json = "{\"tool\": \"Barcode_Distance\", "
    json = json "\"timestamp\": \"" json_escape(run_stamp()) "\", "
    json = json "\"platform\": \"" json_escape(run_platform()) "\", "
    if (manifest_file != "")
        json = json "\"manifest\": \"" json_escape(manifest_file) "\", "
    else
        json = json "\"positional_inputs\": " positional_count ", "
    json = json "\"datasets\": " analyzed_datasets ", "
    json = json "\"samples\": " analyzed_samples ", "
    json = json "\"pairwise_output\": \"" json_escape(output_file) "\", "
    json = json "\"summary_output\": \"" json_escape(summary_file) "\""
    if (write_matrices)
        json = json ", \"matrix_prefix\": \"" json_escape(matrix_prefix) "\""
    if (write_graph)
        json = json ", \"graph\": \"" json_escape(graph_file) "\""
    json = json "}"
    truncate_file(run_metadata)
    print json >> run_metadata
    close(run_metadata)
    log_message("info", "  run metadata: " run_metadata)
}

function load_manifest(path,    status, line, fields, field_count, header_seen, file_column, treatment_column, replicate_column, sample_column, field_index, file_path, treatment, replicate, sample, dataset, sample_key) {
    manifest_directory = directory_name(path)
    input_path_seen[path_key(path)] = 1
    open_read(path)
    while ((status = read_next()) > 0) {
        line = read_line
        sub(/\r$/, "", line)
        if (line ~ /^[[:space:]]*$/)
            continue

        field_count = split(line, fields, "\t")
        if (!header_seen) {
            for (field_index = 1; field_index <= field_count; field_index++) {
                if (tolower(fields[field_index]) == "file")
                    file_column = field_index
                else if (tolower(fields[field_index]) == "treatment")
                    treatment_column = field_index
                else if (tolower(fields[field_index]) == "replicate")
                    replicate_column = field_index
                else if (tolower(fields[field_index]) == "sample")
                    sample_column = field_index
            }
            if (!file_column || !treatment_column || !replicate_column || !sample_column)
                fail("manifest needs file, treatment, replicate, and sample columns: " path)
            header_seen = 1
            continue
        }

        file_path = resolve_manifest_path(fields[file_column])
        treatment = fields[treatment_column]
        replicate = fields[replicate_column]
        sample = fields[sample_column]
        if (file_path == "" || treatment == "" || replicate == "" || sample == "")
            fail("manifest contains an empty required value: " line)
        if (index(sample, SUBSEP) || index(treatment, SUBSEP) || index(replicate, SUBSEP))
            fail("manifest value contains an unsupported control character: " line)
        if (is_output_path(file_path))
            fail("an output path would overwrite an input FASTA: " file_path)
        input_path_seen[path_key(file_path)] = 1

        dataset = treatment "-r" format_replicate(replicate)
        register_dataset(dataset, treatment, replicate)
        sample_key = dataset SUBSEP sample
        if (sample_key in design_sample_seen)
            fail("duplicate sample '" sample "' in dataset " dataset)
        design_sample_seen[sample_key] = 1
        design_sample_count[dataset]++
        dataset_sample[dataset SUBSEP design_sample_count[dataset]] = sample
        dataset_file[dataset SUBSEP design_sample_count[dataset]] = file_path
    }
    read_close()

    if (status < 0)
        fail("could not read manifest: " path)
    if (!header_seen)
        fail("manifest is empty: " path)
}

function load_positional_dataset(    file_index, file_path, sample, sample_key) {
    register_dataset(dataset_label, dataset_label, 1)
    for (file_index = 1; file_index <= positional_count; file_index++) {
        file_path = positional_file[file_index]
        if (is_output_path(file_path))
            fail("an output path would overwrite an input FASTA: " file_path)
        input_path_seen[path_key(file_path)] = 1
        sample = sample_name_from_path(file_path)
        sample_key = dataset_label SUBSEP sample
        if (sample_key in design_sample_seen)
            fail("duplicate normalized sample name '" sample "'")
        design_sample_seen[sample_key] = 1
        design_sample_count[dataset_label]++
        dataset_sample[dataset_label SUBSEP design_sample_count[dataset_label]] = sample
        dataset_file[dataset_label SUBSEP design_sample_count[dataset_label]] = file_path
    }
}

function register_dataset(dataset, treatment, replicate,    treatment_key) {
    if (!(dataset in dataset_seen)) {
        dataset_seen[dataset] = 1
        dataset_name[++dataset_count] = dataset
        dataset_treatment[dataset] = treatment
        dataset_replicate[dataset] = replicate
    } else if (dataset_treatment[dataset] != treatment || dataset_replicate[dataset] != replicate)
        fail("dataset name collision in design: " dataset)

    treatment_key = treatment
    if (!(treatment_key in treatment_seen)) {
        treatment_seen[treatment_key] = 1
        treatment_name[++treatment_count] = treatment
    }
}

function validate_design(    dataset_index, dataset) {
    if (dataset_count == 0)
        fail("the design contains no datasets")
    for (dataset_index = 1; dataset_index <= dataset_count; dataset_index++) {
        dataset = dataset_name[dataset_index]
        if (design_sample_count[dataset] < 2)
            fail("dataset " dataset " needs at least two samples")
    }
    if (path_key(output_file) in input_path_seen)
        fail("pairwise output would overwrite an input FASTA: " output_file)
    if (path_key(summary_file) in input_path_seen)
        fail("summary output would overwrite an input FASTA: " summary_file)
    if (write_graph && (path_key(graph_file) in input_path_seen))
        fail("SVG output would overwrite an input FASTA: " graph_file)
}

function analyze_design(    dataset_index, dataset, sample_index, sample, path) {
    for (dataset_index = 1; dataset_index <= dataset_count; dataset_index++) {
        dataset = dataset_name[dataset_index]
        current_dataset = dataset
        current_sample_count = design_sample_count[dataset]
        log_message("debug", "Barcode_Distance: analyzing dataset " dataset)

        for (sample_index = 1; sample_index <= current_sample_count; sample_index++) {
            sample = dataset_sample[dataset SUBSEP sample_index]
            path = dataset_file[dataset SUBSEP sample_index]
            current_sample[sample_index] = sample
            load_fasta(path, sample)
        }
        log_message("debug", "Barcode_Distance: dataset " dataset ": " \
                    current_sample_count " samples, " feature_count " distinct features")

        calculate_dataset_distances(dataset)
        if (write_matrices)
            write_dataset_matrices(dataset)
        write_dataset_summary(dataset)
        analyzed_datasets++
        analyzed_samples += current_sample_count
        clear_dataset_state()
    }
}

function load_fasta(path, sample,    status, line, sequence, header_seen, records, quality) {
    open_read(path)

    # Find the first record; blank and ";" comment lines are skipped.
    line = ""
    while ((status = read_next()) > 0) {
        line = read_line
        sub(/\r$/, "", line)
        if (line ~ /^[[:space:]]*$/ || line ~ /^;/)
            continue
        break
    }
    if (status < 0)
        fail("could not read FASTA: " path)
    if (status == 0)
        fail("FASTA contains no records: " path)

    if (line ~ /^@/) {
        load_fastq_records(path, sample, line)
        read_close()
        return
    }

    sequence = ""
    header_seen = 0
    records = 0
    while (1) {
        if (line ~ /^>/) {
            if (header_seen) {
                add_sequence(sample, sequence, path)
                records++
            }
            header_seen = 1
            sequence = ""
        } else if (line ~ /^;/ || line ~ /^[[:space:]]*$/) {
            # skip comments and blanks
        } else {
            if (!header_seen)
                fail("FASTA sequence appears before its first header: " path)
            gsub(/[[:space:]]/, "", line)
            sequence = sequence toupper(line)
        }
        if ((status = read_next()) <= 0)
            break
        line = read_line
        sub(/\r$/, "", line)
    }
    read_close()

    if (status < 0)
        fail("could not read FASTA: " path)
    if (header_seen) {
        add_sequence(sample, sequence, path)
        records++
    }
    if (records == 0)
        fail("FASTA contains no records: " path)
}

function load_fastq_records(path, sample,    status, line, sequence, quality, records) {
    # FASTQ inputs use positional four-line records; the header line was
    # already consumed by the caller. Quality strings must agree with their
    # sequence (wrapped FASTQ is not supported); qualities themselves are
    # ignored because distances use complete sequences as features.
    records = 0
    while (1) {
        if ((status = read_next()) <= 0)
            fail("truncated FASTQ record in " path)
        sequence = read_line
        sub(/\r$/, "", sequence)
        gsub(/[[:space:]]/, "", sequence)
        if ((status = read_next()) <= 0)
            fail("truncated FASTQ record in " path)
        if (read_line !~ /^\+/)
            fail("expected '+' line in FASTQ record: " path)
        if ((status = read_next()) <= 0)
            fail("truncated FASTQ record in " path)
        quality = read_line
        sub(/\r$/, "", quality)
        if (length(quality) != length(sequence))
            fail("FASTQ quality length differs from sequence length (wrapped FASTQ is not supported): " path)
        add_sequence(sample, toupper(sequence), path)
        records++

        if ((status = read_next()) <= 0)
            break
        line = read_line
        sub(/\r$/, "", line)
        if (line ~ /^[[:space:]]*$/)
            fail("blank line inside FASTQ data: " path)
        if (line !~ /^@/)
            fail("expected FASTQ header starting with '@': " path)
    }

    if (status < 0)
        fail("could not read FASTQ: " path)
    if (records == 0)
        fail("FASTQ contains no records: " path)
}

function add_sequence(sample, sequence, path,    feature_key) {
    if (sequence == "")
        fail("empty FASTA sequence in " path)
    if (sequence !~ /^[ACGTRYSWKMBDHVN.-]+$/)
        fail("invalid DNA/IUPAC sequence in " path ": " sequence)

    if (!(sequence in feature_seen)) {
        feature_seen[sequence] = 1
        feature_name[++feature_count] = sequence
    }
    feature_key = sample SUBSEP sequence
    feature_abundance[feature_key]++
    sample_total[sample]++
}

function calculate_dataset_distances(dataset,    a, b, sample_a, sample_b, feature_index, feature, count_a, count_b, shared, union_count, minimum_sum, maximum_sum, absolute_sum, total_sum, dot_product, square_a, square_b, hellinger_square, proportion_a, proportion_b, value, denominator, metric_index, metric) {
    for (a = 1; a <= current_sample_count; a++) {
        sample_a = current_sample[a]
        for (b = a + 1; b <= current_sample_count; b++) {
            sample_b = current_sample[b]
            shared = 0
            union_count = 0
            minimum_sum = 0
            maximum_sum = 0
            absolute_sum = 0
            total_sum = 0
            dot_product = 0
            square_a = 0
            square_b = 0
            hellinger_square = 0

            for (feature_index = 1; feature_index <= feature_count; feature_index++) {
                feature = feature_name[feature_index]
                count_a = feature_abundance[sample_a SUBSEP feature] + 0
                count_b = feature_abundance[sample_b SUBSEP feature] + 0
                if (count_a > 0 && count_b > 0)
                    shared++
                if (count_a > 0 || count_b > 0)
                    union_count++
                minimum_sum += count_a < count_b ? count_a : count_b
                maximum_sum += count_a > count_b ? count_a : count_b
                absolute_sum += count_a > count_b ? count_a - count_b : count_b - count_a
                total_sum += count_a + count_b
                dot_product += count_a * count_b
                square_a += count_a * count_a
                square_b += count_b * count_b
                proportion_a = count_a / sample_total[sample_a]
                proportion_b = count_b / sample_total[sample_b]
                hellinger_square += (sqrt(proportion_a) - sqrt(proportion_b)) ^ 2
            }

            set_distance("jaccard", a, b, union_count ? 1 - shared / union_count : 0)
            set_distance("ruzicka", a, b, maximum_sum ? 1 - minimum_sum / maximum_sum : 0)
            set_distance("bray-curtis", a, b, total_sum ? absolute_sum / total_sum : 0)
            denominator = sqrt(square_a) * sqrt(square_b)
            set_distance("cosine", a, b, denominator ? 1 - dot_product / denominator : 0)
            set_distance("hellinger", a, b, sqrt(hellinger_square) / sqrt(2))

            print dataset "\t" dataset_treatment[dataset] "\t" dataset_replicate[dataset] \
                  "\t" sample_a "\t" sample_b "\t" sample_total[sample_a] \
                  "\t" sample_total[sample_b] "\t" shared "\t" union_count \
                  "\t" format_number(pair_distance["jaccard" SUBSEP a SUBSEP b]) \
                  "\t" format_number(pair_distance["ruzicka" SUBSEP a SUBSEP b]) \
                  "\t" format_number(pair_distance["bray-curtis" SUBSEP a SUBSEP b]) \
                  "\t" format_number(pair_distance["cosine" SUBSEP a SUBSEP b]) \
                  "\t" format_number(pair_distance["hellinger" SUBSEP a SUBSEP b]) >> output_file
            for (metric_index = 1; metric_index <= metric_count; metric_index++) {
                metric = metric_name[metric_index]
                value = pair_distance[metric SUBSEP a SUBSEP b]
                metric_sum[metric] += value
                metric_sum_square[metric] += value * value
                metric_pairs[metric]++
            }
        }
    }
}

function set_distance(metric, a, b, value) {
    value = clamp_unit(value)
    pair_distance[metric SUBSEP a SUBSEP b] = value
    pair_distance[metric SUBSEP b SUBSEP a] = value
}

function write_dataset_matrices(dataset,    metric_index, metric, path, row, column, matrix_key) {
    for (metric_index = 1; metric_index <= metric_count; metric_index++) {
        metric = metric_name[metric_index]
        path = matrix_prefix "-" safe_name(dataset) "-" metric "-matrix.tsv"
        matrix_key = path_key(path)
        if (matrix_key in matrix_path_seen)
            fail("matrix filename collision: " path)
        if (matrix_key in input_path_seen)
            fail("matrix output would overwrite an input FASTA: " path)
        if (is_output_path(path))
            fail("matrix output would overwrite another analysis output: " path)
        matrix_path_seen[matrix_key] = 1
        truncate_file(path)

        printf "%s", "sample" >> path
        for (column = 1; column <= current_sample_count; column++)
            printf "\t%s", current_sample[column] >> path
        printf "\n" >> path

        for (row = 1; row <= current_sample_count; row++) {
            printf "%s", current_sample[row] >> path
            for (column = 1; column <= current_sample_count; column++) {
                if (row == column)
                    printf "\t0" >> path
                else
                    printf "\t%s", format_number(pair_distance[metric SUBSEP row SUBSEP column]) >> path
            }
            printf "\n" >> path
        }
        close(path)
        written_matrices++
    }
}

function write_dataset_summary(dataset,    metric_index, metric, count, mean, variance, deviation, treatment) {
    treatment = dataset_treatment[dataset]
    for (metric_index = 1; metric_index <= metric_count; metric_index++) {
        metric = metric_name[metric_index]
        count = metric_pairs[metric]
        mean = count ? metric_sum[metric] / count : 0
        variance = count ? metric_sum_square[metric] / count - mean * mean : 0
        if (variance < 0 && variance > -0.000000000001)
            variance = 0
        deviation = variance > 0 ? sqrt(variance) : 0

        print dataset "\t" treatment "\t" dataset_replicate[dataset] \
              "\t" current_sample_count "\t" count "\t" metric \
              "\t" format_number(mean) "\t" format_number(deviation) >> summary_file
        graph_sum[treatment SUBSEP metric] += mean
        graph_count[treatment SUBSEP metric]++
    }
}

function write_svg_graph(path,    cell_width, cell_height, left_margin, top_margin, right_margin, bottom_margin, width, height, metric_index, treatment_index, metric, treatment, mean, x, y, color, legend_index, legend_x, legend_y, label_x, label_y) {
    cell_width = 30
    cell_height = 42
    left_margin = 145
    top_margin = 75
    right_margin = 155
    bottom_margin = 145
    width = left_margin + treatment_count * cell_width + right_margin
    if (width < 760)
        width = 760
    height = top_margin + metric_count * cell_height + bottom_margin

    truncate_file(path)
    print "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" >> path
    print "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" width "\" height=\"" height "\" viewBox=\"0 0 " width " " height "\" role=\"img\" aria-labelledby=\"title desc\">" >> path
    print "<title id=\"title\">Mean barcode distance by treatment</title>" >> path
    print "<desc id=\"desc\">Generated entirely by POSIX AWK from complete simulated barcode sequence counts.</desc>" >> path
    print "<rect width=\"100%\" height=\"100%\" fill=\"#ffffff\"/>" >> path
    print "<text x=\"" left_margin "\" y=\"32\" font-family=\"sans-serif\" font-size=\"20\" font-weight=\"bold\" fill=\"#222222\">Mean barcode distance by treatment</text>" >> path
    print "<text x=\"" left_margin "\" y=\"54\" font-family=\"sans-serif\" font-size=\"12\" fill=\"#555555\">Direct full-sequence distances; no Mash, OTU picker, R, or plotting library</text>" >> path

    for (metric_index = 1; metric_index <= metric_count; metric_index++) {
        metric = metric_name[metric_index]
        y = top_margin + (metric_index - 1) * cell_height
        print "<text x=\"" left_margin - 10 "\" y=\"" y + 25 "\" text-anchor=\"end\" font-family=\"sans-serif\" font-size=\"13\" fill=\"#222222\">" xml_escape(metric_label(metric)) "</text>" >> path

        for (treatment_index = 1; treatment_index <= treatment_count; treatment_index++) {
            treatment = treatment_name[treatment_index]
            mean = graph_sum[treatment SUBSEP metric] / graph_count[treatment SUBSEP metric]
            x = left_margin + (treatment_index - 1) * cell_width
            color = distance_color(mean)
            print "<rect x=\"" x "\" y=\"" y "\" width=\"" cell_width - 1 "\" height=\"" cell_height - 1 "\" fill=\"" color "\" stroke=\"#ffffff\" stroke-width=\"1\"><title>" xml_escape(treatment) " / " xml_escape(metric_label(metric)) ": " format_number(mean) "</title></rect>" >> path
        }
    }

    label_y = top_margin + metric_count * cell_height + 8
    for (treatment_index = 1; treatment_index <= treatment_count; treatment_index++) {
        treatment = treatment_name[treatment_index]
        label_x = left_margin + (treatment_index - 1) * cell_width + 18
        print "<text transform=\"translate(" label_x "," label_y ") rotate(55)\" text-anchor=\"start\" font-family=\"sans-serif\" font-size=\"11\" fill=\"#333333\">" xml_escape(treatment) "</text>" >> path
    }

    legend_x = left_margin + treatment_count * cell_width + 35
    legend_y = top_margin
    print "<text x=\"" legend_x "\" y=\"" legend_y - 12 "\" font-family=\"sans-serif\" font-size=\"12\" fill=\"#333333\">distance</text>" >> path
    for (legend_index = 0; legend_index <= 10; legend_index++) {
        mean = legend_index / 10
        y = legend_y + (10 - legend_index) * 16
        print "<rect x=\"" legend_x "\" y=\"" y "\" width=\"18\" height=\"16\" fill=\"" distance_color(mean) "\"/>" >> path
        if (legend_index % 5 == 0)
            print "<text x=\"" legend_x + 25 "\" y=\"" y + 12 "\" font-family=\"sans-serif\" font-size=\"11\" fill=\"#333333\">" sprintf("%.1f", mean) "</text>" >> path
    }
    print "</svg>" >> path
    close(path)
}

function clear_dataset_state() {
    clear_array(current_sample)
    clear_array(feature_seen)
    clear_array(feature_name)
    clear_array(feature_abundance)
    clear_array(sample_total)
    clear_array(pair_distance)
    clear_array(metric_sum)
    clear_array(metric_sum_square)
    clear_array(metric_pairs)
    feature_count = 0
    current_sample_count = 0
}

function is_output_path(path) {
    path = path_key(path)
    return path == path_key(output_file) || path == path_key(summary_file) || \
           (write_graph && path == path_key(graph_file))
}

function path_key(path) {
    while (substr(path, 1, 2) == "./")
        path = substr(path, 3)
    gsub(/\/+/, "/", path)
    return path
}

function resolve_manifest_path(path) {
    if (substr(path, 1, 1) == "/" || index(path, "/") || manifest_directory == ".")
        return path
    return manifest_directory "/" path
}

function directory_name(path) {
    if (!index(path, "/"))
        return "."
    sub(/\/[^\/]*$/, "", path)
    return path == "" ? "/" : path
}

function base_name(path) {
    sub(/^.*\//, "", path)
    return path
}

function sample_name_from_path(path,    name) {
    name = base_name(path)
    sub(/[.]gz$/, "", name)
    sub(/[.](fasta|fa|fna|fas|fastq|fq)$/, "", name)
    if (name == "")
        fail("could not derive a sample name from: " path)
    return name
}

function format_replicate(value) {
    if (value ~ /^[0-9]+$/)
        return sprintf("%03d", value + 0)
    return safe_name(value)
}

function safe_name(value) {
    gsub(/[^A-Za-z0-9_.-]/, "_", value)
    return value
}

function clamp_unit(value) {
    if (value < 0)
        return 0
    if (value > 1)
        return 1
    return value
}

function format_number(value) {
    if (value > -0.0000000000005 && value < 0.0000000000005)
        value = 0
    return sprintf("%.12g", value)
}

function metric_label(metric) {
    if (metric == "jaccard") return "Jaccard"
    if (metric == "ruzicka") return "Ruzicka"
    if (metric == "bray-curtis") return "Bray-Curtis"
    if (metric == "cosine") return "Cosine"
    if (metric == "hellinger") return "Hellinger"
    return metric
}

function distance_color(value,    red, green, blue) {
    value = clamp_unit(value)
    red = int(245 - 126 * value)
    green = int(248 - 204 * value)
    blue = int(255 - 96 * value)
    return sprintf("#%02x%02x%02x", red, green, blue)
}

function xml_escape(value,    result, character_index, character) {
    result = ""
    for (character_index = 1; character_index <= length(value); character_index++) {
        character = substr(value, character_index, 1)
        if (character == "&") result = result "&amp;"
        else if (character == "<") result = result "&lt;"
        else if (character == ">") result = result "&gt;"
        else if (character == "\"") result = result "&quot;"
        else if (character == "'") result = result "&apos;"
        else result = result character
    }
    return result
}

function truncate_file(path) {
    printf "%s", "" > path
    close(path)
}

function clear_array(array,    key) {
    for (key in array)
        delete array[key]
}

function report() {
    log_message("info", "Barcode Distance complete")
    log_message("info", "  datasets: " analyzed_datasets)
    log_message("info", "  samples: " analyzed_samples)
    log_message("info", "  pairwise distances: " output_file)
    log_message("info", "  summary: " summary_file)
    if (write_matrices)
        log_message("info", "  matrices: " written_matrices " (prefix " matrix_prefix ")")
    if (write_graph)
        log_message("info", "  graph: " graph_file)
}

function usage() {
    print "Usage: Barcode_Distance --manifest FILE [options]"
    print "       Barcode_Distance [options] SAMPLE1.fasta SAMPLE2.fasta [...]"
    print ""
    print "Calculate distances directly from complete barcode sequences and render"
    print "an SVG heatmap using POSIX AWK only. No Mash, R, OTU picker, or plotting"
    print "library is required."
    print ""
    print "Inputs:"
    print "      --manifest FILE       Barcode_Experiment manifest (multiple datasets)"
    print "      --dataset NAME        label for positional FASTA inputs"
    print ""
    print "Outputs:"
    print "      --output FILE         pairwise wide TSV"
    print "      --summary FILE        per-dataset metric summary TSV"
    print "      --matrix-prefix TEXT  prefix for square metric matrices"
    print "      --graph FILE          AWK-generated SVG distance heatmap"
    print "      --no-matrices         suppress square matrix files"
    print "      --no-graph            suppress the SVG heatmap"
    print "      --loglevel LEVEL      error, warning, info (default), or debug"
    print "      --run-metadata FILE   write a JSON run sidecar"
    print "  -h, --help                show this help"
    print ""
    print "Inputs may be plain or gzip-compressed (.gz) FASTA or FASTQ. FASTQ"
    print "records use positional four-line records; qualities must agree with"
    print "their sequence (wrapped FASTQ is not supported) and are otherwise"
    print "ignored."
    print ""
    print "Metrics: binary Jaccard, quantitative Jaccard/Ruzicka, Bray-Curtis,"
    print "cosine distance, and Hellinger distance. All use full barcode sequences"
    print "as features; abundance-aware metrics use per-sample sequence counts."
}

function fail(message) {
    print "Barcode_Distance: " message > stderr
    print "Try 'Barcode_Distance --help' for usage." > stderr
    exit 2
}
