export default function FullLengthSequences({
  nucleotideVDJ,
  aminoAcidVDJ,
}: {
  nucleotideVDJ: string;
  aminoAcidVDJ: string;
}) {
  return (
    <>
      <h1>Full Length Sequences</h1>
      <p>Fasta Formatted</p>
      <h2>Nucleotide VDJ</h2>
      <p>{nucleotideVDJ}</p>
      <h2>Amino-Acid VDJ</h2>
      <p>{aminoAcidVDJ}</p>
    </>
  );
}
