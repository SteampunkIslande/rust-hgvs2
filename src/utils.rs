use std::io::BufRead;

#[derive(Default)]
pub struct GeneDescription {
    chrom: String,
    start: i64,
    end: i64,
    id: String,
    strand: String,
    cds_start: i64,
    cds_end: i64,
    gene_name: String,
    exons: Vec<(i64, i64)>,
}

pub fn read_refgene(infile: impl std::io::Read) -> impl Iterator<Item = GeneDescription> {
    read_genepred(infile, true)
}

/// GenePred extension format:
///     http://genome.ucsc.edu/FAQ/FAQformat.html#GenePredExt
///     Column definitions:
///     0. string name;                 "Name of gene (usually transcript_id from GTF)"
///     1. string chrom;                "Chromosome name"
///     2. char[1] strand;              "+ or - for strand"
///     3. uint txStart;                "Transcription start position"
///     4. uint txEnd;                  "Transcription end position"
///     5. uint cdsStart;               "Coding region start"
///     6. uint cdsEnd;                 "Coding region end"
///     7. uint exonCount;              "Number of exons"
///     8. uint[exonCount] exonStarts;  "Exon start positions"
///     9. uint[exonCount] exonEnds;    "Exon end positions"
///     10. uint id;                    "Unique identifier"
///     11. string name2;               "Alternate name (e.g. gene_id from GTF)"
pub fn read_genepred(
    infile: impl std::io::Read,
    skip_first_column: bool,
) -> impl Iterator<Item = GeneDescription> {
    let file = std::io::BufReader::new(infile);
    file.lines().map_while(Result::ok).filter_map(move |line| {
        if line.starts_with("#") {
            return None;
        } else {
            let row: Vec<&str> = line
                .split("\t")
                .skip(if skip_first_column { 1 } else { 0 })
                .collect();

            // Originally skips trailing comma, but using filter_map this is not needed
            let exon_starts: Vec<i64> = row[8]
                .split(",")
                .filter_map(|s| s.parse::<i64>().ok())
                .collect();
            let exon_ends: Vec<i64> = row[8]
                .split(",")
                .filter_map(|s| s.parse::<i64>().ok())
                .collect();
            let exons: Vec<(i64, i64)> = exon_starts.into_iter().zip(exon_ends).collect();

            return Some(GeneDescription {
                chrom: row[1].into(),
                start: row[3].parse().ok()?,
                end: row[4].parse().ok()?,
                id: row[0].into(),
                strand: row[2].into(),
                cds_start: row[5].parse().ok()?,
                cds_end: row[6].parse().ok()?,
                gene_name: row[11].into(),
                exons: exons,
            });
        }
    })
}
