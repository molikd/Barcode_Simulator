#!/usr/bin/awk -f

# Build the three paper-style analysis inputs from Barcode_Experiment FASTAs:
# cluster counts (OTU proxy), exact-sequence counts (ASV proxy), and exact
# k-mer Jaccard distances (deterministic stand-in for a Mash sketch).

BEGIN {
    stderr = "/dev/stderr"
    if (project == "")
        project = "paper"
    if (k == "")
        k = 5
    if (k !~ /^[0-9]+$/ || k + 0 < 1)
        fail("k must be a positive integer")
}

FNR == 1 {
    start_file(FILENAME)
}

/^>/ {
    current_cluster = header_value($0, "cluster")
    if (current_cluster == "")
        fail("missing cluster annotation in " FILENAME ": " $0)
    next
}

/^[[:space:]]*$/ {
    next
}

{
    if (current_cluster == "")
        fail("sequence encountered before a FASTA header in " FILENAME)

    sequence = toupper($0)
    otu_count[current_dataset SUBSEP current_cluster SUBSEP current_sample]++
    asv_count[current_dataset SUBSEP sequence SUBSEP current_sample]++

    otu_taxon_key = current_dataset SUBSEP current_cluster
    if (!(otu_taxon_key in otu_taxon_seen)) {
        otu_taxon_seen[otu_taxon_key] = 1
        otu_taxon_count[current_dataset]++
        otu_taxon_index = otu_taxon_count[current_dataset]
        otu_taxon[current_dataset SUBSEP otu_taxon_index] = current_cluster
    }

    asv_taxon_key = current_dataset SUBSEP sequence
    if (!(asv_taxon_key in asv_taxon_seen)) {
        asv_taxon_seen[asv_taxon_key] = 1
        asv_taxon_count[current_dataset]++
        asv_taxon_index = asv_taxon_count[current_dataset]
        asv_taxon[current_dataset SUBSEP asv_taxon_index] = sequence
    }

    if (length(sequence) < k)
        fail("sequence is shorter than k in " FILENAME)
    for (position = 1; position <= length(sequence) - k + 1; position++) {
        word = substr(sequence, position, k)
        kmer_key = current_dataset SUBSEP current_sample SUBSEP word
        if (!(kmer_key in kmer_seen)) {
            kmer_seen[kmer_key] = 1
            kmer_total_key = current_dataset SUBSEP current_sample
            kmer_count[kmer_total_key]++
            kmer_index_value = kmer_count[kmer_total_key]
            kmer_name[kmer_total_key SUBSEP kmer_index_value] = word
        }
    }
}

END {
    if (failed)
        exit 2
    if (dataset_count == 0)
        fail("no FASTA inputs were read")

    metadata_file = project "-analysis-meta.tsv"
    truncate_file(metadata_file)
    print "File\tTreatment" >> metadata_file
    close(metadata_file)

    for (dataset_index = 1; dataset_index <= dataset_count; dataset_index++) {
        dataset = dataset_name[dataset_index]
        sort_scoped(sample_name, sample_count[dataset], dataset)
        sort_scoped(otu_taxon, otu_taxon_count[dataset], dataset)
        sort_scoped(asv_taxon, asv_taxon_count[dataset], dataset)

        otu_file = dataset "_otu_table.txt"
        asv_file = dataset "_asv_otu_table.txt"
        mash_file = dataset "_mash_dists.txt"
        asv_mash_file = dataset "_asv_mash_dists.txt"

        write_count_table(otu_file, dataset, "otu")
        write_count_table(asv_file, dataset, "asv")
        write_kmer_distances(mash_file, dataset)
        copy_kmer_distances(mash_file, asv_mash_file)

        treatment = dataset
        sub("^" project "-", "", treatment)
        sub(/-r[0-9][0-9][0-9]$/, "", treatment)
        print otu_file "\t" treatment "-OTU" >> metadata_file
        print asv_file "\t" treatment "-ASV" >> metadata_file
        close(metadata_file)
    }
}

function start_file(path,    name, dataset_key, sample_key) {
    name = path
    sub(/^.*\//, "", name)
    if (name !~ /-r[0-9][0-9][0-9]-s[0-9][0-9][0-9][.]fasta$/)
        fail("unexpected experiment FASTA name: " path)

    current_sample = name
    sub(/[.]fasta$/, "", current_sample)
    current_dataset = current_sample
    sub(/-s[0-9][0-9][0-9]$/, "", current_dataset)
    current_cluster = ""

    dataset_key = current_dataset
    if (!(dataset_key in dataset_seen)) {
        dataset_seen[dataset_key] = 1
        dataset_name[++dataset_count] = current_dataset
    }

    sample_key = current_dataset SUBSEP current_sample
    if (!(sample_key in sample_seen)) {
        sample_seen[sample_key] = 1
        sample_count[current_dataset]++
        sample_index_value = sample_count[current_dataset]
        sample_name[current_dataset SUBSEP sample_index_value] = current_sample
    }
}

function header_value(header, wanted,    fields, count, field_index, prefix) {
    sub(/^>/, "", header)
    count = split(header, fields, /[|]/)
    prefix = wanted "="
    for (field_index = 1; field_index <= count; field_index++) {
        if (substr(fields[field_index], 1, length(prefix)) == prefix)
            return substr(fields[field_index], length(prefix) + 1)
    }
    return ""
}

function write_count_table(path, dataset, kind,    sample_index, taxon_index, sample, taxon, value) {
    truncate_file(path)
    printf "%s", "#OTU ID" >> path
    for (sample_index = 1; sample_index <= sample_count[dataset]; sample_index++)
        printf "\t%s", sample_name[dataset SUBSEP sample_index] >> path
    printf "\n" >> path

    if (kind == "otu") {
        for (taxon_index = 1; taxon_index <= otu_taxon_count[dataset]; taxon_index++) {
            taxon = otu_taxon[dataset SUBSEP taxon_index]
            printf "cluster-%s", taxon >> path
            for (sample_index = 1; sample_index <= sample_count[dataset]; sample_index++) {
                sample = sample_name[dataset SUBSEP sample_index]
                value = otu_count[dataset SUBSEP taxon SUBSEP sample] + 0
                printf "\t%d", value >> path
            }
            printf "\n" >> path
        }
    } else {
        for (taxon_index = 1; taxon_index <= asv_taxon_count[dataset]; taxon_index++) {
            taxon = asv_taxon[dataset SUBSEP taxon_index]
            printf "asv-%03d", taxon_index >> path
            for (sample_index = 1; sample_index <= sample_count[dataset]; sample_index++) {
                sample = sample_name[dataset SUBSEP sample_index]
                value = asv_count[dataset SUBSEP taxon SUBSEP sample] + 0
                printf "\t%d", value >> path
            }
            printf "\n" >> path
        }
    }
    close(path)
}

function write_kmer_distances(path, dataset,    a, b, sample_a, sample_b, intersection, union, distance, kmer_index, word, total_key_a, total_key_b) {
    truncate_file(path)
    for (a = 1; a <= sample_count[dataset]; a++) {
        sample_a = sample_name[dataset SUBSEP a]
        for (b = a + 1; b <= sample_count[dataset]; b++) {
            sample_b = sample_name[dataset SUBSEP b]
            intersection = 0
            total_key_a = dataset SUBSEP sample_a
            total_key_b = dataset SUBSEP sample_b
            for (kmer_index = 1; kmer_index <= kmer_count[total_key_a]; kmer_index++) {
                word = kmer_name[total_key_a SUBSEP kmer_index]
                if ((dataset SUBSEP sample_b SUBSEP word) in kmer_seen)
                    intersection++
            }
            union = kmer_count[total_key_a] + kmer_count[total_key_b] - intersection
            distance = union == 0 ? 0 : 1 - intersection / union
            print sample_a ".fasta.msh\t" sample_b ".fasta.msh\t" \
                  sprintf("%.12g", distance) "\t0\t1/1" >> path
        }
    }
    close(path)
}

function copy_kmer_distances(source, destination,    status, line) {
    truncate_file(destination)
    while ((status = getline line < source) > 0)
        print line >> destination
    close(source)
    close(destination)
    if (status < 0)
        fail("could not copy k-mer distances from " source)
}

function sort_scoped(values, count, scope,    i, j, value) {
    for (i = 2; i <= count; i++) {
        value = values[scope SUBSEP i]
        j = i - 1
        while (j >= 1 && values[scope SUBSEP j] > value) {
            values[scope SUBSEP (j + 1)] = values[scope SUBSEP j]
            j--
        }
        values[scope SUBSEP (j + 1)] = value
    }
}

function truncate_file(path) {
    printf "%s", "" > path
    close(path)
}

function fail(message) {
    print "build_paper_matrices.awk: " message > stderr
    failed = 1
    exit 2
}
