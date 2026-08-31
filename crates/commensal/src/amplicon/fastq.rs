/// Read sequences out of a FASTQ stream, one per record.
///
/// Records are four lines. Nothing here validates the quality line beyond
/// skipping it, since dereplication reads the sequence alone.
pub fn read_fastq<R: std::io::BufRead>(reader: R) -> Result<Vec<Vec<u8>>, crate::Error> {
    let mut out = Vec::new();
    let mut lines = reader.lines();
    while let Some(header) = lines.next() {
        let header = header?;
        if header.is_empty() {
            continue;
        }
        let Some(sequence) = lines.next() else { break };
        let sequence = sequence?;
        let _plus = lines.next();
        let _quality = lines.next();
        out.push(sequence.into_bytes());
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fastq_parsing_takes_the_sequence_line() {
        let text = "@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+\nIIII\n";
        let got = read_fastq(std::io::BufReader::new(text.as_bytes())).expect("parses");
        assert_eq!(got, vec![b"ACGT".to_vec(), b"TTTT".to_vec()]);
    }
}
