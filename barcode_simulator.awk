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
    if (seed == "") {
        srand()
        reported_seed = "automatic"
    } else {
        srand(seed)
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

    while ((status = getline line < path) > 0) {
        sub(/\r$/, "", line)
        if (line ~ /^[[:space:]]*$/)
            continue

        if (mode == 0)
            mode = (line ~ /^>/ ? 2 : 1)

        if (mode == 1) {
            if (line ~ /^>/)
                fail("reuse file mixes plain sequences and FASTA records: " path)
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

    if (status < 0)
        fail("could not read reuse file: " path)
    if (mode == 2 && sequence != "")
        add_reused_sequence(sequence, path)
    close(path)

    if (master_count == 0)
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

function generate_fastas(    fasta_number, requested, output_file, i, j, temporary) {
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

        output_file = project_name "-" fasta_number ".fasta"
        truncate_file(output_file)
        for (i = 1; i <= requested; i++) {
            print ">FASTA-" fasta_number "_" i header_suffix >> output_file
            print master[sample_order[i]] >> output_file
        }
        close(output_file)

        for (i = 1; i <= master_count; i++)
            delete sample_order[i]
    }
}

function truncate_file(path) {
    printf "%s", "" > path
    close(path)
}

function random_integer(minimum, maximum) {
    return minimum + int(rand() * (maximum - minimum + 1))
}

function report() {
    print "Barcode Simulator complete" > stderr
    print "  seed: " reported_seed > stderr
    print "  master sequences: " master_count > stderr
    if (!only_master)
        print "  FASTA files: " num_fasta > stderr
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
