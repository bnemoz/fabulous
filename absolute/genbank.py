from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from datetime import datetime

def create_genbank_record(
    sequence: str,
    record_id: str = "example_id",
    name: str = "Example",
    description: str = "Synthetic GenBank Record",
    organism: str = "synthetic construct",
    features: list = None
) -> SeqRecord:
    """
    Create a GenBank-style Biopython SeqRecord.

    Parameters:
    - sequence (str): The DNA or protein sequence.
    - record_id (str): Accession or record ID.
    - name (str): Short name (max 10 chars for GenBank).
    - description (str): Description of the record.
    - organism (str): Source organism.
    - features (list): Optional list of Bio.SeqFeature objects.

    Returns:
    - SeqRecord object that can be exported as GenBank.
    """
    seq = Seq(str(sequence))
    record = SeqRecord(
        seq,
        id=record_id,
        name=name[:10],  # GenBank name limit
        description=description,
        annotations={
            "molecule_type": "DNA",
            "topology": "linear",
            "data_file_division": "SYN",
            "date": datetime.today().strftime("%d-%b-%Y"),
            "organism": organism,
            "taxonomy": organism.split()  # crude taxonomy
        }
    )

    if features:
        record.features.extend(features)

    return record


def create_gb_from_ab(ab, to_file:str=None):

    features = []
    features.append(SeqFeature(location=FeatureLocation(ab['v_sequence_start'], ab['v_sequence_end']),
                                type='gene',
                                id='V gene',
                                qualifiers={"standard_name":"V gene"},
                                sub_features=None,)
                   )
    
    if ab['locus'] == 'IGH':
        d_start = ab['sequence_oriented'].find(ab['d_sequence'])
        d_end = d_start + len(ab['d_sequence'])
        features.append(SeqFeature(location=FeatureLocation(d_start, d_end),
                                type='gene',
                                id='D gene',
                                qualifiers={"standard_name":"D gene"},
                                sub_features=None,)
                       )
        
    features.append(SeqFeature(location=FeatureLocation(ab['j_sequence_start'], ab['j_sequence_end']),
                                type='gene',
                                id='J gene',
                                qualifiers={"standard_name":"J gene"},
                                sub_features=None,)
                   )

    for domain in ['fwr1', 'cdr1', 'fwr2', 'cdr2', 'fwr3', 'cdr3', 'fwr4']:
        domain_name =  domain.upper()[:3] + ('H' if ab['locus'] == 'IGH' else 'L') + domain[-1]
        domain_start = ab['sequence_oriented'].find(ab[domain])
        domain_end = domain_start + + len(ab[domain])
        features.append(SeqFeature(location=FeatureLocation(domain_start, domain_end),
                                type='domain' if domain.startswith('f') else 'region',
                                id=domain_name,
                                qualifiers={"standard_name":domain_name},
                                sub_features=None,)
                   )

    junction_start = ab['sequence_oriented'].find(ab['junction'])
    junction_end = junction_start + len(ab['junction'])
    features.append(SeqFeature(location=FeatureLocation(junction_start, junction_end),
                                type='misc',
                                id='Junction',
                                qualifiers={"standard_name":"Junction", "translation":ab['junction_aa']},
                                sub_features=None,)
                   )
    
    features.append(SeqFeature(location=FeatureLocation(ab['v_sequence_start'], ab['j_sequence_end']),
                                type='CDS',
                                id='VDJ',
                                qualifiers={"standard_name":"VDJ", "gene":"VDJ", "codon_start":1, "translation":ab['sequence_aa']},
                                sub_features=None,)
                   )

    if ab['leader']:
        sp_start = ab['sequence_oriented'].find(ab['signal_peptide'])
        sp_end = sp_start + len(ab['signal_peptide'])
        features.append(SeqFeature(location=FeatureLocation(sp_start, sp_end),
                                type='sig_peptide',
                                id='Signal Peptide',
                                qualifiers={"standard_name":"Signal Peptide", "translation":ab['signal_peptide_aa']},
                                sub_features=None,)
                   )

    if ab['trailer']:
        if ab['locus'] == 'IGH':
            for domain in ['CH1', 'CH2', 'CH3', 'CH4']:
                domain_start = ab['sequence_oriented'].find(ab[domain])
                domain_end = domain_start + + len(ab[domain])
                features.append(SeqFeature(location=FeatureLocation(domain_start, domain_end),
                                        type='domain',
                                        id=domain,
                                        qualifiers={"standard_name":domain},
                                        sub_features=None,)
                           )
        

    record = create_genbank_record(
        sequence = ab['sequence_oriented'],
        record_id = ab.id,
        name = ab.id,
        description = f"Antibody {ab.id} ({'heavy' if ab['locus'] == 'IGH' else 'light'} chain)",
        organism = ab['species'],
        features = features,
    )

    if to_file:
        with open(to_file, "w") as output_handle:
            SeqIO.write(record, output_handle, "genbank")
        return to_file
        
    else:
        return record