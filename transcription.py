from Bio.Seq import MutableSeq, Seq


def transcribe(dna):
    """Transcribe a DNA sequence into RNA.
    If given a string, returns a new string object.
    Given a Seq or MutableSeq, returns a new Seq object
    with an RNA alphabet.
    Trying to transcribe a protein or RNA sequence
    raises an exception.
    e.g.
     >>> transcribe("ACTGN")
    'ACUGN'
    """


    if isinstance(dna, Seq):
        return dna.transcribe()
    elif isinstance(dna, MutableSeq):
        return dna.toseq().transcribe()
    else:
        return dna.replace("T", "U").replace("t", "u")
