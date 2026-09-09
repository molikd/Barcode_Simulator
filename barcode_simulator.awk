#!/usr/bin/awk -f

# Barcode Simulator
#
# A dependency-free, POSIX-AWK rewrite of the original Bash simulator. The
# historical command-line options are retained; clearer aliases are provided
# for new code.

BEGIN {
    stderr = "/dev/stderr"

    set_defaults()
    parse_arguments()

    if (show_help) {
        usage()
        exit 0
    }

    validate_options()
    if (dry_run) {
        dry_run_report()
        exit 0
    }
    if (!no_disk_check)
        disk_preflight(estimate_worst_bytes())
    if (seed == "") {
        srand()
        seed_random(int(rand() * 2147483646) + 1)
        reported_seed = "automatic"
    } else {
        seed_random(seed)
        reported_seed = seed
    }

    if (reuse_file != "")
        load_reuse_file(reuse_file)
    else
        generate_master_list()

    if (save_file != "")
        save_master_list(save_file)

    if (!only_master)
        generate_fastas()

    report()
    if (run_metadata != "")
        write_run_metadata()
    exit 0
}

function set_defaults() {
    num_fasta = 1
    total_clusters = 1
    sequences_per_cluster = 1
    min_per_file = 1
    max_per_file = 1
    min_length = 100
    max_length = 100
    min_differences = 0
    max_differences = 0
    project_name = "default"
    reuse_file = ""
    reference_file = ""
    save_file = ""
    header_suffix = "-run"
    only_master = 0
    show_help = 0
    seed = ""
    loglevel = 2
    dry_run = 0
    no_disk_check = 0
    output_format = "fasta"
    use_gzip = 0
    max_bytes_per_fasta = 0
    run_metadata = ""
}

function parse_arguments(    i, arg, option, value, equals_at) {
    for (i = 1; i < ARGC; i++) {
        arg = ARGV[i]

        if (arg == "--")
            continue

        if (substr(arg, 1, 2) == "--") {
            equals_at = index(arg, "=")
            if (equals_at) {
                option = substr(arg, 1, equals_at - 1)
                value = substr(arg, equals_at + 1)
            } else {
                option = arg
                value = ""
            }

            if (option == "--help") {
                show_help = 1
                continue
            }
            if (option == "--only") {
                only_master = 1
                continue
            }
            if (option == "--dry-run") {
                dry_run = 1
                continue
            }
            if (option == "--no-disk-check") {
                no_disk_check = 1
                continue
            }
            if (option == "--gzip") {
                use_gzip = 1
                continue
            }

            if (!equals_at) {
                if (i + 1 >= ARGC)
                    fail("missing value for " option)
                value = ARGV[++i]
            }
            set_option(option, value)
            continue
        }

        if (substr(arg, 1, 1) == "-") {
            option = substr(arg, 1, 2)
            value = substr(arg, 3)

            if (option == "-h" && value == "") {
                show_help = 1
                continue
            }
            if (option == "-O" && value == "") {
                only_master = 1
                continue
            }
            if (length(option) != 2)
                fail("unknown option: " arg)
            if (value == "") {
                if (i + 1 >= ARGC)
                    fail("missing value for " option)
                value = ARGV[++i]
            }
            set_option(option, value)
            continue
        }

        fail("unexpected positional argument: " arg)
    }

    # Prevent AWK from treating command-line options as input files.
    ARGC = 1
}

function set_option(option, value) {
    if (option == "-f" || option == "--num-fasta")
        num_fasta = value
    else if (option == "-o" || option == "--total-otus" || option == "--total-clusters")
        total_clusters = value
    else if (option == "-s" || option == "--number-of-seq-per-otu" || option == "--sequences-per-cluster")
        sequences_per_cluster = value
    else if (option == "-W" || option == "--max-number-otus-per-file" || option == "--max-sequences-per-file")
        max_per_file = value
    else if (option == "-w" || option == "--min-number-otus-per-file" || option == "--min-sequences-per-file")
        min_per_file = value
    else if (option == "-L" || option == "--max-length-of-sequences" || option == "--max-length")
        max_length = value
    else if (option == "-l" || option == "--min-length-of-sequences" || option == "--min-length")
        min_length = value
    else if (option == "-D" || option == "--max-injected-differences-in-otus" || option == "--max-differences")
        max_differences = value
    else if (option == "-d" || option == "--min-injected-differences-in-otus" || option == "--min-differences")
        min_differences = value
    else if (option == "-p" || option == "--project-name")
        project_name = value
    else if (option == "-r" || option == "--reuse")
        reuse_file = normalize_optional_file(value)
    else if (option == "-R" || option == "--gen-ref")
        reference_file = normalize_optional_file(value)
    else if (option == "-S" || option == "--save")
        save_file = normalize_optional_file(value)
    else if (option == "-a" || option == "--add")
        header_suffix = value
    else if (option == "--seed")
        seed = value
    else if (option == "--loglevel")
        parse_loglevel(value)
    else if (option == "--format")
        set_output_format(value)
    else if (option == "--max-bytes-per-fasta")
        max_bytes_per_fasta = value
    else if (option == "--run-metadata")
        run_metadata = value
    else
        fail("unknown option: " option)
}

function normalize_optional_file(value) {
    if (value == "false")
        return ""
    return value
}

function validate_options() {
    require_positive_integer("number of FASTA files", num_fasta)
    require_positive_integer("total clusters", total_clusters)
    require_positive_integer("sequences per cluster", sequences_per_cluster)
    require_positive_integer("minimum sequences per file", min_per_file)
    require_positive_integer("maximum sequences per file", max_per_file)
    require_positive_integer("minimum sequence length", min_length)
    require_positive_integer("maximum sequence length", max_length)
    require_nonnegative_integer("minimum differences", min_differences)
    require_nonnegative_integer("maximum differences", max_differences)
    require_nonnegative_integer("maximum bytes per FASTA", max_bytes_per_fasta)
    if (seed != "")
        require_nonnegative_integer("seed", seed)

    if (min_per_file > max_per_file)
        fail("minimum sequences per file cannot exceed maximum")
    if (min_length > max_length)
        fail("minimum sequence length cannot exceed maximum")
    if (min_differences > max_differences)
        fail("minimum differences cannot exceed maximum")
    if (reuse_file == "" && max_differences > min_length)
        fail("maximum differences cannot exceed the shortest generated sequence")
    if (project_name == "")
        fail("project name cannot be empty")
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

function set_output_format(value) {
    value = tolower(value)
    if (value != "fasta" && value != "fastq")
        fail("--format must be 'fasta' or 'fastq' (received '" value "')")
    output_format = value
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

function format_bytes(count,    unit) {
    if (count < 1024)
        return sprintf("%d B", count)
    count /= 1024
    if (count < 1024)
        return sprintf("%.1f KB", count)
    count /= 1024
    if (count < 1024)
        return sprintf("%.1f MB", count)
    count /= 1024
    if (count < 1024)
        return sprintf("%.1f GB", count)
    count /= 1024
    return sprintf("%.1f TB", count)
}

function output_directory(path,    directory) {
    if (!index(path, "/"))
        return "."
    directory = path
    sub(/\/[^\/]*$/, "", directory)
    return directory == "" ? "/" : directory
}

function disk_free_bytes(directory,    command, status, line, seen_header, fields, available) {
    command = "df -P " shell_quote(directory)
    seen_header = 0
    available = -1
    while ((status = (command | getline line)) > 0) {
        if (!seen_header) {
            seen_header = 1
            continue
        }
        if (split(line, fields) >= 4 && fields[4] ~ /^[0-9]+$/)
            available = fields[4] * 1024
        break
    }
    close(command)
    if (status < 0 || !seen_header)
        return -1
    return available
}

function disk_preflight(needed,    free, margin) {
    margin = int(needed * 1.1) + 1
    free = disk_free_bytes(".")
    if (free < 0) {
        log_message("warning", "Barcode_Simulator: could not determine free disk space; continuing")
        return
    }
    if (margin > free)
        fail("insufficient disk space: need ~" format_bytes(margin) ", have " \
             format_bytes(free) " (use --no-disk-check to override)")
    log_message("info", "  disk check: need ~" format_bytes(margin) ", free " format_bytes(free))
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
# through gzip -dc so compressed reuse files just work.
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

# Output stream: plain files truncate-then-append, --gzip streams records
# through gzip so compressed outputs need no post-processing pass.
function out_begin(path) {
    out_path = path
    if (use_gzip) {
        out_cmd = "gzip -9 -- > " shell_quote(path)
    } else {
        truncate_file(path)
        out_cmd = ""
    }
    out_bytes = 0
    out_files[++out_file_count] = path
}

function out_line(text) {
    if (out_cmd == "")
        print text >> out_path
    else
        print text | out_cmd
    out_bytes += length(text) + 1
}

function out_end() {
    if (out_cmd == "")
        close(out_path)
    else
        close(out_cmd)
    out_cmd = ""
}

function require_positive_integer(name, value) {
    if (value !~ /^[0-9]+$/ || value + 0 < 1)
        fail(name " must be a positive integer (received '" value "')")
}

function require_nonnegative_integer(name, value) {
    if (value !~ /^[0-9]+$/)
        fail(name " must be a non-negative integer (received '" value "')")
}

function generate_master_list(    cluster, copy, sequence_length, base, differences) {
    if (reference_file != "")
        truncate_file(reference_file)

    for (cluster = 1; cluster <= total_clusters; cluster++) {
        sequence_length = random_integer(min_length, max_length)
        base = random_dna(sequence_length)

        if (reference_file != "") {
            print ">REF-" cluster >> reference_file
            print base >> reference_file
        }

        for (copy = 1; copy <= sequences_per_cluster; copy++) {
            differences = random_integer(min_differences, max_differences)
            master[++master_count] = mutate(base, differences)
        }
    }

    if (reference_file != "")
        close(reference_file)
}

function random_dna(sequence_length,    i, sequence, bases) {
    bases = "ACGT"
    sequence = ""
    for (i = 1; i <= sequence_length; i++)
        sequence = sequence substr(bases, random_integer(1, 4), 1)
    return sequence
}

function mutate(sequence, changes,    changed, position, old_base, new_base, bases, i) {
    for (i in mutation_positions)
        delete mutation_positions[i]

    bases = "ACGT"
    changed = 0
    while (changed < changes) {
        position = random_integer(1, length(sequence))
        if (position in mutation_positions)
            continue

        old_base = substr(sequence, position, 1)
        do {
            new_base = substr(bases, random_integer(1, 4), 1)
        } while (new_base == old_base)

        sequence = substr(sequence, 1, position - 1) new_base substr(sequence, position + 1)
        mutation_positions[position] = 1
        changed++
    }
    return sequence
}

function load_reuse_file(path,    status, line, mode, sequence) {
    mode = 0
    sequence = ""

    open_read(path)
    while ((status = read_next()) > 0) {
        line = read_line
        sub(/\r$/, "", line)
        if (line ~ /^[[:space:]]*$/)
            continue

        if (mode == 0)
            mode = (line ~ /^>/ ? 2 : 1)

        if (mode == 1) {
            if (line ~ /^>/) {
                read_close()
                fail("reuse file mixes plain sequences and FASTA records: " path)
            }
            gsub(/[[:space:]]/, "", line)
            add_reused_sequence(line, path)
        } else if (line ~ /^>/) {
            if (sequence != "")
                add_reused_sequence(sequence, path)
            sequence = ""
        } else {
            gsub(/[[:space:]]/, "", line)
            sequence = sequence line
        }
    }
    read_close()

    if (status < 0)
        fail("could not read reuse file: " path)
    if (mode == 2 && sequence != "")
        add_reused_sequence(sequence, path)

    if (master_count == 0)
        fail("reuse file contains no sequences: " path)
}

# Read-only pass over a reuse file for --dry-run estimates. Sets the
# reuse_count, reuse_total_length, and reuse_max_length globals.
function count_reuse_sequences(path,    status, line, mode, sequence) {
    mode = 0
    sequence = ""
    reuse_count = 0
    reuse_total_length = 0
    reuse_max_length = 0

    open_read(path)
    while ((status = read_next()) > 0) {
        line = read_line
        sub(/\r$/, "", line)
        if (line ~ /^[[:space:]]*$/)
            continue

        if (mode == 0)
            mode = (line ~ /^>/ ? 2 : 1)

        if (mode == 1) {
            if (line ~ /^>/) {
                read_close()
                fail("reuse file mixes plain sequences and FASTA records: " path)
            }
            gsub(/[[:space:]]/, "", line)
            reuse_count++
            reuse_total_length += length(line)
            if (length(line) > reuse_max_length)
                reuse_max_length = length(line)
        } else if (line ~ /^>/) {
            if (sequence != "") {
                reuse_count++
                reuse_total_length += length(sequence)
                if (length(sequence) > reuse_max_length)
                    reuse_max_length = length(sequence)
            }
            sequence = ""
        } else {
            gsub(/[[:space:]]/, "", line)
            sequence = sequence line
        }
    }
    read_close()

    if (status < 0)
        fail("could not read reuse file: " path)
    if (mode == 2 && sequence != "") {
        reuse_count++
        reuse_total_length += length(sequence)
        if (length(sequence) > reuse_max_length)
            reuse_max_length = length(sequence)
    }
    if (reuse_count == 0)
        fail("reuse file contains no sequences: " path)
}

function add_reused_sequence(sequence, path) {
    sequence = toupper(sequence)
    if (sequence !~ /^[ACGTRYSWKMBDHVN.-]+$/)
        fail("reuse file contains a non-DNA sequence: " path)
    master[++master_count] = sequence
}

function save_master_list(path,    i) {
    truncate_file(path)
    for (i = 1; i <= master_count; i++)
        print master[i] >> path
    close(path)
}

function generate_fastas(    fasta_number, requested, part, output_file, i, j, temporary, header, sequence) {
    for (fasta_number = num_fasta; fasta_number >= 1; fasta_number--) {
        requested = random_integer(min_per_file, max_per_file)
        if (requested > master_count)
            requested = master_count

        for (i = 1; i <= master_count; i++)
            sample_order[i] = i
        for (i = master_count; i > 1; i--) {
            j = random_integer(1, i)
            temporary = sample_order[i]
            sample_order[i] = sample_order[j]
            sample_order[j] = temporary
        }

        part = 0
        output_file = split_path(fasta_number, part)
        out_begin(output_file)
        log_message("debug", "Barcode_Simulator: writing " output_file)
        for (i = 1; i <= requested; i++) {
            header = "FASTA-" fasta_number "_" i header_suffix
            sequence = master[sample_order[i]]
            if (max_bytes_per_fasta > 0 && out_bytes > 0 && \
                out_bytes + record_bytes(header, sequence) > max_bytes_per_fasta) {
                out_end()
                part++
                output_file = split_path(fasta_number, part)
                out_begin(output_file)
                log_message("debug", "Barcode_Simulator: writing " output_file)
            }
            write_record(header, sequence)
        }
        out_end()

        for (i = 1; i <= master_count; i++)
            delete sample_order[i]
    }
}

function split_path(fasta_number, part,    path) {
    path = project_name "-" fasta_number
    if (max_bytes_per_fasta > 0)
        path = path ".part" sprintf("%03d", part)
    return path sequence_extension()
}

function sequence_extension(    extension) {
    extension = (output_format == "fastq" ? ".fastq" : ".fasta")
    if (use_gzip)
        extension = extension ".gz"
    return extension
}

function record_bytes(header, sequence) {
    if (output_format == "fastq")
        return length(header) + 2 * length(sequence) + 6
    return length(header) + length(sequence) + 3
}

function write_record(header, sequence) {
    if (output_format == "fastq") {
        out_line("@" header)
        out_line(sequence)
        out_line("+")
        out_line(quality_string(length(sequence)))
    } else {
        out_line(">" header)
        out_line(sequence)
    }
}

function quality_string(sequence_length,    block) {
    # Uniform synthetic Q40 qualities. Barcode_Simulator models sequences,
    # not instrument error, so every base receives the same high score.
    if (sequence_length > length(quality_block)) {
        block = sprintf("%*s", sequence_length, "")
        gsub(/ /, "I", block)
        quality_block = block
    }
    return substr(quality_block, 1, sequence_length)
}

function estimate_master() {
    # Sets the est_master, est_avg_len, and est_max_len globals.
    if (reuse_file != "") {
        count_reuse_sequences(reuse_file)
        est_master = reuse_count
        est_avg_len = int(reuse_total_length / reuse_count)
        est_max_len = reuse_max_length
    } else {
        est_master = total_clusters * sequences_per_cluster
        est_avg_len = int((min_length + max_length) / 2)
        est_max_len = max_length
    }
}

function record_overhead(    header) {
    header = "FASTA-1_1" header_suffix
    if (output_format == "fastq")
        return length(header) + 6
    return length(header) + 3
}

function estimate_expected_bytes() {
    estimate_master()
    return int(num_fasta * ((min_per_file + max_per_file) / 2) * \
               (est_avg_len + record_overhead())) + \
           est_master * (est_avg_len + 1) + \
           total_clusters * (est_max_len + 16)
}

function estimate_worst_bytes() {
    estimate_master()
    return int(num_fasta * max_per_file * (est_max_len + record_overhead())) + \
           est_master * (est_max_len + 1) + \
           total_clusters * (est_max_len + 16)
}

function dry_run_report(    expected, worst, free) {
    expected = estimate_expected_bytes()
    worst = estimate_worst_bytes()
    print "Barcode_Simulator dry run (no files written)" > stderr
    print "  format: " output_format (use_gzip ? " (gzipped)" : "") > stderr
    print "  master sequences: ~" est_master > stderr
    if (!only_master) {
        print "  sequence files: " num_fasta \
              (max_bytes_per_fasta > 0 ? \
               " (split past " format_bytes(max_bytes_per_fasta) ")" : "") > stderr
        print "  sequences per file: ~" int((min_per_file + max_per_file) / 2) \
              " expected, " max_per_file " worst case" > stderr
    }
    print "  expected output: ~" format_bytes(expected) > stderr
    print "  worst-case output: ~" format_bytes(worst) > stderr
    free = disk_free_bytes(".")
    if (free < 0)
        print "  free disk space: unknown" > stderr
    else
        print "  free disk space: " format_bytes(free) > stderr
}

function write_run_metadata(    json, i, outputs) {
    outputs = ""
    for (i = 1; i <= out_file_count; i++)
        outputs = outputs (i > 1 ? ", " : "") "\"" json_escape(out_files[i]) "\""
    json = "{\"tool\": \"Barcode_Simulator\", "
    json = json "\"timestamp\": \"" json_escape(run_stamp()) "\", "
    json = json "\"platform\": \"" json_escape(run_platform()) "\", "
    json = json "\"seed\": \"" json_escape(reported_seed) "\", "
    json = json "\"parameters\": {" \
        "\"num_fasta\": " num_fasta ", " \
        "\"total_clusters\": " total_clusters ", " \
        "\"sequences_per_cluster\": " sequences_per_cluster ", " \
        "\"min_sequences_per_file\": " min_per_file ", " \
        "\"max_sequences_per_file\": " max_per_file ", " \
        "\"min_length\": " min_length ", " \
        "\"max_length\": " max_length ", " \
        "\"min_differences\": " min_differences ", " \
        "\"max_differences\": " max_differences ", " \
        "\"project_name\": \"" json_escape(project_name) "\", " \
        "\"format\": \"" output_format "\", " \
        "\"gzip\": " (use_gzip ? "true" : "false") ", " \
        "\"max_bytes_per_fasta\": " max_bytes_per_fasta ", " \
        "\"reuse_file\": \"" json_escape(reuse_file) "\", " \
        "\"reference_file\": \"" json_escape(reference_file) "\", " \
        "\"save_file\": \"" json_escape(save_file) "\", " \
        "\"only_master\": " (only_master ? "true" : "false") "}, "
    json = json "\"master_sequences\": " master_count ", "
    json = json "\"sequence_files\": [" outputs "]"
    if (reference_file != "")
        json = json ", \"reference\": \"" json_escape(reference_file) "\""
    if (save_file != "")
        json = json ", \"master_file\": \"" json_escape(save_file) "\""
    json = json "}"
    truncate_file(run_metadata)
    print json >> run_metadata
    close(run_metadata)
    log_message("info", "  run metadata: " run_metadata)
}

function truncate_file(path) {
    printf "%s", "" > path
    close(path)
}

function random_integer(minimum, maximum) {
    return minimum + int(random_unit() * (maximum - minimum + 1))
}

function seed_random(value) {
    random_state = value % 2147483647
    if (random_state <= 0)
        random_state += 2147483646
}

function random_unit() {
    random_state = (random_state * 48271) % 2147483647
    return random_state / 2147483647
}

function report() {
    log_message("info", "Barcode Simulator complete")
    log_message("info", "  seed: " reported_seed)
    log_message("info", "  master sequences: " master_count)
    if (!only_master)
        log_message("info", "  sequence files: " out_file_count)
}

function usage() {
    print "Usage: Barcode_Simulator [options]"
    print ""
    print "Generate related DNA barcode/gene-copy sequences and sample them into FASTA files."
    print ""
    print "Core options:"
    print "  -f, --num-fasta N                    FASTA files to generate (default: 1)"
    print "  -o, --total-otus N                   clusters/OTUs to generate (default: 1)"
    print "      --total-clusters N               clearer alias for --total-otus"
    print "  -s, --number-of-seq-per-otu N        variants per cluster (default: 1)"
    print "  -w, --min-number-otus-per-file N     minimum sequences per output (default: 1)"
    print "  -W, --max-number-otus-per-file N     maximum sequences per output (default: 1)"
    print "  -l, --min-length-of-sequences N      minimum generated length (default: 100)"
    print "  -L, --max-length-of-sequences N      maximum generated length (default: 100)"
    print "  -d, --min-injected-differences-in-otus N"
    print "                                           minimum substitutions (default: 0)"
    print "  -D, --max-injected-differences-in-otus N"
    print "                                           maximum substitutions (default: 0)"
    print "  -p, --project-name PREFIX            output prefix (default: default)"
    print "  -a, --add TEXT                       append TEXT to FASTA headers (default: -run)"
    print "      --seed N                         random seed for reproducible runs"
    print "      --loglevel LEVEL                 error, warning, info (default), or debug"
    print "      --format FMT                     fasta (default) or fastq (synthetic Q40)"
    print "      --gzip                           compress sequence outputs with gzip"
    print "      --max-bytes-per-fasta N          split outputs past N bytes (0 = unlimited)"
    print "      --run-metadata FILE              write a JSON run sidecar"
    print "      --dry-run                        estimate outputs and disk needs; write nothing"
    print "      --no-disk-check                  skip the disk-space preflight check"
    print ""
    print "Master/reference options:"
    print "  -r, --reuse FILE                     reuse FASTA or one-sequence-per-line input"
    print "  -R, --gen-ref FILE                   write one reference per generated cluster"
    print "  -S, --save FILE                      save master sequences, one per line"
    print "  -O, --only                           build/reuse master sequences only"
    print "  -h, --help                           show this help"
}

function fail(message) {
    print "Barcode_Simulator: " message > stderr
    print "Try 'Barcode_Simulator --help' for usage." > stderr
    exit 2
}
