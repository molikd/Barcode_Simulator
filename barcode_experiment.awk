#!/usr/bin/awk -f

# Barcode Experiment
#
# Generate controlled metabarcode experiments across the five effects used in
# Molik, Pfrender, and Emrich (2020): abundance distribution, a shared conserved
# region, sequence errors, variable sequence length, and variable sample depth.

BEGIN {
    stderr = "/dev/stderr"
    set_defaults()
    detect_profiles()
    parse_arguments()

    if (show_help) {
        usage()
        exit 0
    }

    validate_options()
    prepare_outputs()
    run_design()
    report()
    exit 0
}

function set_defaults() {
    project_name = "experiment"
    samples = 68
    replicates = 1
    clusters = 68
    variants_per_cluster = 10

    fixed_depth = 136
    min_depth = 14
    max_depth = 1360

    baseline_length = 500
    min_variable_length = 350
    max_variable_length = 500

    min_errors = 1
    max_errors = 10
    conserved_sequence = "CGTCACACTTCATGATGGAATTGA"

    seed = 1
    factorial = 0
    selected_code = 0
    selected_effects = 0
    write_truth = 1
    manifest_file = ""
    truth_file = ""
    header_modifier = "sequence"
    show_help = 0
}

function detect_profiles(    i) {
    for (i = 1; i < ARGC; i++) {
        if (ARGV[i] == "--paper") {
            paper_profile = 1
            factorial = 1
            samples = 68
            replicates = 10
            clusters = 68
            variants_per_cluster = 10
            fixed_depth = 1360
            min_depth = 140
            max_depth = 13600
            baseline_length = 500
            min_variable_length = 350
            max_variable_length = 500
            min_errors = 1
            max_errors = 10
        }
    }
}

function parse_arguments(    i, argument, option, value, equals_at) {
    for (i = 1; i < ARGC; i++) {
        argument = ARGV[i]
        if (argument == "--")
            continue

        if (substr(argument, 1, 2) == "--") {
            equals_at = index(argument, "=")
            if (equals_at) {
                option = substr(argument, 1, equals_at - 1)
                value = substr(argument, equals_at + 1)
            } else {
                option = argument
                value = ""
            }

            if (option == "--help") {
                show_help = 1
                continue
            }
            if (option == "--factorial") {
                factorial = 1
                continue
            }
            if (option == "--paper")
                continue
            if (option == "--no-truth") {
                write_truth = 0
                continue
            }
            if (is_long_effect_flag(option)) {
                add_effect_flag(option)
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

        if (substr(argument, 1, 1) == "-") {
            option = substr(argument, 1, 2)
            value = substr(argument, 3)

            if (option == "-h" && value == "") {
                show_help = 1
                continue
            }
            if (is_short_effect_flag(option) && value == "") {
                add_effect_flag(option)
                continue
            }
            if (value == "") {
                if (i + 1 >= ARGC)
                    fail("missing value for " option)
                value = ARGV[++i]
            }
            set_option(option, value)
            continue
        }

        fail("unexpected positional argument: " argument)
    }
    ARGC = 1
}

function is_long_effect_flag(option) {
    return option == "--abundance" || option == "--abundances" || \
           option == "--abundences" || option == "--conserved" || \
           option == "--errors" || option == "--lengths" || \
           option == "--variable-depth" || option == "--picks"
}

function is_short_effect_flag(option) {
    return option == "-A" || option == "-C" || option == "-E" || \
           option == "-L" || option == "-N"
}

function add_effect_flag(option) {
    if (option == "-A" || option == "--abundance" || \
        option == "--abundances" || option == "--abundences")
        add_effect_bit(1)
    else if (option == "-C" || option == "--conserved")
        add_effect_bit(2)
    else if (option == "-E" || option == "--errors")
        add_effect_bit(4)
    else if (option == "-L" || option == "--lengths")
        add_effect_bit(8)
    else if (option == "-N" || option == "--variable-depth" || option == "--picks")
        add_effect_bit(16)
}

function add_effect_bit(bit) {
    if (int(selected_code / bit) % 2 == 0)
        selected_code += bit
    selected_effects = 1
}

function set_option(option, value) {
    if (option == "-P" || option == "--project" || option == "--project-name")
        project_name = value
    else if (option == "-M" || option == "--modifier" || option == "--seq-mod")
        header_modifier = value
    else if (option == "--effects")
        parse_effect_list(value)
    else if (option == "--samples" || option == "--num-fasta")
        samples = value
    else if (option == "--replicates")
        replicates = value
    else if (option == "--clusters" || option == "--total-otus")
        clusters = value
    else if (option == "--variants-per-cluster" || option == "--number-of-seq-per-otu")
        variants_per_cluster = value
    else if (option == "--fixed-depth")
        fixed_depth = value
    else if (option == "--min-depth")
        min_depth = value
    else if (option == "--max-depth")
        max_depth = value
    else if (option == "--baseline-length")
        baseline_length = value
    else if (option == "--min-variable-length")
        min_variable_length = value
    else if (option == "--max-variable-length")
        max_variable_length = value
    else if (option == "--min-errors")
        min_errors = value
    else if (option == "--max-errors")
        max_errors = value
    else if (option == "--conserved-sequence")
        conserved_sequence = toupper(value)
    else if (option == "--seed")
        seed = value
    else if (option == "--manifest")
        manifest_file = value
    else if (option == "--truth") {
        truth_file = value
        write_truth = 1
    } else
        fail("unknown option: " option)
}

function parse_effect_list(list,    items, count, i, effect) {
    gsub(/[[:space:]]/, "", list)
    list = tolower(list)
    if (list == "" || list == "none" || list == "baseline") {
        selected_effects = 1
        return
    }
    if (list == "all") {
        selected_code = 31
        selected_effects = 1
        return
    }

    count = split(list, items, ",")
    for (i = 1; i <= count; i++) {
        effect = items[i]
        if (effect == "abundance" || effect == "abundances" || effect == "a")
            add_effect_bit(1)
        else if (effect == "conserved" || effect == "c")
            add_effect_bit(2)
        else if (effect == "errors" || effect == "error" || effect == "e")
            add_effect_bit(4)
        else if (effect == "lengths" || effect == "length" || effect == "l")
            add_effect_bit(8)
        else if (effect == "depth" || effect == "variable-depth" || \
                 effect == "picks" || effect == "n")
            add_effect_bit(16)
        else
            fail("unknown effect in --effects: " effect)
    }
    selected_effects = 1
}

function validate_options() {
    require_positive_integer("samples", samples)
    require_positive_integer("replicates", replicates)
    require_positive_integer("clusters", clusters)
    require_positive_integer("variants per cluster", variants_per_cluster)
    require_positive_integer("fixed depth", fixed_depth)
    require_positive_integer("minimum depth", min_depth)
    require_positive_integer("maximum depth", max_depth)
    require_positive_integer("baseline length", baseline_length)
    require_positive_integer("minimum variable length", min_variable_length)
    require_positive_integer("maximum variable length", max_variable_length)
    require_positive_integer("minimum errors", min_errors)
    require_positive_integer("maximum errors", max_errors)
    require_nonnegative_integer("seed", seed)

    if (factorial && selected_effects)
        fail("--factorial cannot be combined with individual effect selections")
    if (min_depth > max_depth)
        fail("minimum depth cannot exceed maximum depth")
    if (min_variable_length > max_variable_length)
        fail("minimum variable length cannot exceed maximum")
    if (min_errors > max_errors)
        fail("minimum errors cannot exceed maximum")
    if (max_errors > baseline_length)
        fail("maximum errors cannot exceed the baseline sequence length")
    if ((factorial || (has_effect(selected_code, 4) && has_effect(selected_code, 8))) && \
        max_errors > min_variable_length)
        fail("maximum errors cannot exceed the shortest variable-length sequence")
    if (clusters < 3 && (factorial || has_effect(selected_code, 1)))
        fail("the abundance effect requires at least three clusters")
    if (project_name == "")
        fail("project name cannot be empty")
    if (conserved_sequence !~ /^[ACGTRYSWKMBDHVN.-]+$/)
        fail("conserved sequence must contain DNA/IUPAC characters")
}

function require_positive_integer(name, value) {
    if (value !~ /^[0-9]+$/ || value + 0 < 1)
        fail(name " must be a positive integer (received '" value "')")
}

function require_nonnegative_integer(name, value) {
    if (value !~ /^[0-9]+$/)
        fail(name " must be a non-negative integer (received '" value "')")
}

function prepare_outputs() {
    if (manifest_file == "")
        manifest_file = project_name "-manifest.tsv"
    if (truth_file == "")
        truth_file = project_name "-truth.tsv"

    truncate_file(manifest_file)
    print "file\ttreatment\treplicate\tsample\treads\tabundance\tconserved\terrors\tlengths\tvariable_depth\tclusters\tvariants_per_cluster\tseed\tfixed_depth\tmin_depth\tmax_depth\tbaseline_length\tmin_variable_length\tmax_variable_length\tmin_errors\tmax_errors\tconserved_sequence" >> manifest_file
    close(manifest_file)

    if (write_truth) {
        truncate_file(truth_file)
        print "file\ttreatment\treplicate\tsample\tread\tcluster\tvariant\tabundance_tier\tsequence_length\tsubstitutions" >> truth_file
        close(truth_file)
    }
}

function run_design(    replicate, code, first_code, last_code) {
    if (factorial) {
        first_code = 0
        last_code = 31
    } else {
        first_code = selected_code
        last_code = selected_code
    }

    for (replicate = 1; replicate <= replicates; replicate++) {
        generate_templates(replicate)
        for (code = first_code; code <= last_code; code++)
            generate_treatment(replicate, code)
        clear_array(template_sequence)
        clear_array(variable_length)
    }
}

function generate_templates(replicate,    cluster, template_length) {
    template_length = baseline_length
    if (max_variable_length > template_length)
        template_length = max_variable_length

    seed_random(derive_seed(replicate, 0, 17))
    for (cluster = 1; cluster <= clusters; cluster++)
        template_sequence[cluster] = random_dna(template_length)

    seed_random(derive_seed(replicate, 0, 31))
    for (cluster = 1; cluster <= clusters; cluster++)
        variable_length[cluster] = random_integer(min_variable_length, max_variable_length)
}

function generate_treatment(replicate, code,    treatment, reference_file, sample) {
    current_abundance = has_effect(code, 1)
    current_conserved = has_effect(code, 2)
    current_errors = has_effect(code, 4)
    current_lengths = has_effect(code, 8)
    current_depth = has_effect(code, 16)
    treatment = treatment_code(code)

    build_variant_pool(replicate, code)
    reference_file = project_name "-" treatment "-r" pad3(replicate) "-reference.fasta"
    write_reference_file(reference_file, treatment)

    seed_random(derive_seed(replicate, code, 53))
    for (sample = 1; sample <= samples; sample++)
        write_sample(replicate, sample, treatment)

    clear_array(variant_sequence)
    clear_array(variant_substitutions)
    clear_array(cluster_reference)
    clear_array(cluster_tier)
}

function build_variant_pool(replicate, code,    cluster, variant, variable, length_value, differences, key, pool_code) {
    pool_code = (current_errors ? 1 : 0) + (current_lengths ? 2 : 0)
    seed_random(derive_seed(replicate, pool_code, 71))
    set_abundance_tiers()

    for (cluster = 1; cluster <= clusters; cluster++) {
        length_value = current_lengths ? variable_length[cluster] : baseline_length
        variable = substr(template_sequence[cluster], 1, length_value)
        cluster_reference[cluster] = (current_conserved ? conserved_sequence : "") variable

        for (variant = 1; variant <= variants_per_cluster; variant++) {
            differences = current_errors ? random_integer(min_errors, max_errors) : 0
            key = cluster SUBSEP variant
            variant_sequence[key] = (current_conserved ? conserved_sequence : "") mutate(variable, differences)
            variant_substitutions[key] = differences
        }
    }
}

function set_abundance_tiers(    cluster) {
    high_clusters = int(clusters / 6)
    if (high_clusters < 1)
        high_clusters = 1
    medium_clusters = int(clusters * 2 / 6)
    if (medium_clusters < 1)
        medium_clusters = 1
    low_clusters = clusters - high_clusters - medium_clusters
    if (low_clusters < 1)
        fail("could not allocate abundance tiers; use at least three clusters")

    for (cluster = 1; cluster <= clusters; cluster++) {
        if (cluster <= high_clusters)
            cluster_tier[cluster] = "high"
        else if (cluster <= high_clusters + medium_clusters)
            cluster_tier[cluster] = "middling"
        else
            cluster_tier[cluster] = "low"
    }
}

function write_reference_file(path, treatment,    cluster) {
    truncate_file(path)
    for (cluster = 1; cluster <= clusters; cluster++) {
        print ">cluster-" cluster "|tier=" cluster_tier[cluster] "|effects=" treatment >> path
        print cluster_reference[cluster] >> path
    }
    close(path)
}

function write_sample(replicate, sample, treatment,    reads, path, read_number, cluster, variant, key, tier, sequence) {
    reads = current_depth ? random_integer(min_depth, max_depth) : fixed_depth
    path = project_name "-" treatment "-r" pad3(replicate) "-s" pad3(sample) ".fasta"
    truncate_file(path)

    for (read_number = 1; read_number <= reads; read_number++) {
        cluster = choose_cluster()
        variant = random_integer(1, variants_per_cluster)
        key = cluster SUBSEP variant
        tier = cluster_tier[cluster]
        sequence = variant_sequence[key]

        print ">" header_modifier "-" pad3(sample) "_read-" pad6(read_number) \
              "|cluster=" cluster "|variant=" variant "|tier=" tier \
              "|effects=" treatment "|replicate=" replicate >> path
        print sequence >> path

        if (write_truth)
            print path "\t" treatment "\t" replicate "\t" sample "\t" read_number \
                  "\t" cluster "\t" variant "\t" tier "\t" length(sequence) \
                  "\t" variant_substitutions[key] >> truth_file
    }
    close(path)
    if (write_truth)
        close(truth_file)

    print path "\t" treatment "\t" replicate "\t" sample "\t" reads \
          "\t" current_abundance "\t" current_conserved "\t" current_errors \
          "\t" current_lengths "\t" current_depth "\t" clusters \
          "\t" variants_per_cluster "\t" seed "\t" fixed_depth "\t" min_depth \
          "\t" max_depth "\t" baseline_length "\t" min_variable_length \
          "\t" max_variable_length "\t" min_errors "\t" max_errors \
          "\t" conserved_sequence >> manifest_file
    close(manifest_file)
    generated_samples++
    generated_reads += reads
}

function choose_cluster(    tier, start, count) {
    if (!current_abundance)
        return random_integer(1, clusters)

    # Each abundance tier receives equal expected read mass. Because the tiers
    # contain clusters in an approximate 1:2:3 ratio, individual clusters have
    # high, middling, and low expected abundance.
    tier = random_integer(1, 3)
    if (tier == 1) {
        start = 1
        count = high_clusters
    } else if (tier == 2) {
        start = high_clusters + 1
        count = medium_clusters
    } else {
        start = high_clusters + medium_clusters + 1
        count = low_clusters
    }
    return start + random_integer(0, count - 1)
}

function treatment_code(code,    result) {
    result = ""
    if (has_effect(code, 1)) result = result "A"
    if (has_effect(code, 2)) result = result "C"
    if (has_effect(code, 4)) result = result "E"
    if (has_effect(code, 8)) result = result "L"
    if (has_effect(code, 16)) result = result "N"
    return result == "" ? "O" : result
}

function has_effect(code, bit) {
    return int(code / bit) % 2
}

function random_dna(sequence_length,    sequence, bases, i) {
    sequence = ""
    bases = "ACGT"
    for (i = 1; i <= sequence_length; i++)
        sequence = sequence substr(bases, random_integer(1, 4), 1)
    return sequence
}

function mutate(sequence, changes,    changed, position, old_base, new_base, bases, i) {
    clear_array(mutation_positions)
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

function derive_seed(replicate, code, stream,    value) {
    value = seed + replicate * 104729 + code * 1009 + stream * 7919
    return value % 2147483647
}

function pad3(value) {
    return sprintf("%03d", value)
}

function pad6(value) {
    return sprintf("%06d", value)
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
    print "Barcode Experiment complete" > stderr
    print "  design: " (paper_profile ? "Molik et al. (2020) 2^5 factorial" : \
          (factorial ? "full 2^5 factorial" : treatment_code(selected_code))) > stderr
    print "  samples: " generated_samples > stderr
    print "  reads: " generated_reads > stderr
    print "  manifest: " manifest_file > stderr
    if (write_truth)
        print "  truth: " truth_file > stderr
}

function usage() {
    print "Usage: Barcode_Experiment [effects] [options]"
    print ""
    print "Generate metabarcode datasets with the five effects from Molik et al. (2020)."
    print ""
    print "Effect design:"
    print "  -A, --abundance          high/middling/low cluster abundance tiers"
    print "  -C, --conserved          prepend a shared conserved sequence"
    print "  -E, --errors             inject 1-10 substitutions per variant"
    print "  -L, --lengths            vary cluster length from 350-500 bp"
    print "  -N, --variable-depth     vary sample depth from 14-1360 reads"
    print "      --effects LIST       comma-separated effect names, 'all', or 'none'"
    print "      --factorial          generate baseline plus all 31 effect combinations"
    print "      --paper              Section 2.3 profile: 32 conditions x 10 replicates"
    print ""
    print "Design size and naming:"
    print "  -P, --project PREFIX     output prefix (default: experiment)"
    print "  -M, --modifier TEXT      FASTA header prefix (default: sequence)"
    print "      --samples N          samples per treatment (default: 68)"
    print "      --replicates N       replicates per treatment (default: 1)"
    print "      --clusters N         species/sequence clusters (default: 68)"
    print "      --variants-per-cluster N  variants per cluster (default: 10)"
    print "      --seed N             design seed (default: 1)"
    print ""
    print "Effect levels:"
    print "      --fixed-depth N      baseline reads per sample (default: 136)"
    print "      --min-depth N        variable-depth minimum (default: 14)"
    print "      --max-depth N        variable-depth maximum (default: 1360)"
    print "      --baseline-length N  baseline sequence length (default: 500)"
    print "      --min-variable-length N  variable-length minimum (default: 350)"
    print "      --max-variable-length N  variable-length maximum (default: 500)"
    print "      --min-errors N       error-effect minimum (default: 1)"
    print "      --max-errors N       error-effect maximum (default: 10)"
    print "      --conserved-sequence DNA  shared prefix"
    print ""
    print "The --paper profile uses 68 samples, 68 clusters, 10 variants per cluster,"
    print "1360 fixed reads, 140-13600 variable reads, and the published sequence"
    print "length/error levels. Explicit size and effect-level options may override it."
    print ""
    print "Provenance outputs:"
    print "      --manifest FILE      sample-level design table"
    print "      --truth FILE         read-level ground-truth table"
    print "      --no-truth           do not write read-level truth"
    print "  -h, --help               show this help"
}

function fail(message) {
    print "Barcode_Experiment: " message > stderr
    print "Try 'Barcode_Experiment --help' for usage." > stderr
    exit 2
}
