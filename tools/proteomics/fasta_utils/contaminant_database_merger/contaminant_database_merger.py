"""
Contaminant Database Merger
===========================
Append common contaminant proteins (cRAP) to a target FASTA database with a configurable
prefix, and remove duplicates.

Usage
-----
    python contaminant_database_merger.py --input target.fasta --add-crap --prefix CONT_ --output merged.fasta
    python contaminant_database_merger.py --input target.fasta --contaminants custom.fasta --output merged.fasta
"""

from typing import List

import click
import pyopenms as oms

# Built-in common contaminant sequences (subset of cRAP database)
BUILTIN_CONTAMINANTS = [
    ("P00761", "IVGGYTCGANTVPYQVSLNSGYHFCGGSLINSQWVVSAAHCYKSGIQVRLGEDNINVVEG"
               "NEQFISASKSIVHPSYNSNTLNNDIMLIKLKSAASLNSRVASISLPTSCASAGTQCLISGWGNTK"),
    ("P00766", "IVGGYTCGANTVPYQVSLNSGYHFCGGSLINSQWVVSAAHCYKSGIQVRLGEDNINVVEG"
               "NEQFISASKSIVHPSYNSNTLNNDIMLIKL"),
    ("P02769", "MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPF"
               "DEHVKLVNELTEFAKTCVADESHAGCEK"),
    ("P00432", "MDSASKLLSALFLGALLGASCIAAPPGAQKGESVTLNCKPVLDDFEPATQRFNGNIFHYP"
               "NTISATGRWNKEENAISEMFQNHFTKSNP"),
    ("P02533", "MSRSSFRRGSGSRSGSRSSSYSLGSRSGGFSSSGGFGGSRSLYGLGASRSSGSSYGLGGG"
               "SSSGGSTGGIRATSGFASRSSGGGYSSSGGFSG"),
    ("P35527", "SFSSRSASCISGGYRGSGGRSYSCGSCGISGGYRGSGGRSYSSGSCGISGGYRGSGGRSYS"
               "CSISCGIASGGYRGSGGRSYSCGSCGISGG"),
    ("P04264", "SRQFSSRSGYRSGGSYGGGSSGGGSISGSSYGSRSGSYRSGGSSGGSYGSRSGSYRSGGS"
               "GGSYGGSRSGSYRSGGSSGGSYGSRSGS"),
    ("P13645", "MTTSQYGRRSSQYGSYSQSTSYRGSGSRSSGYRSGGSSGYSSSGGYRSGGSSGGSYGSRS"
               "GSYRSGGSSGGSYGSRSGSYRSGGSSGGSYG"),
    ("P02768", "DAHKSEVAHRF KDLGEENFKALVLIAFAQYLQQCPFEDHVKLVNEVTEFAKTCVADESAA"
               "NCDKSLHTLFGDKLCTVATLRETYGEMADCCAK"),
    ("P01966", "VLSAADKGNVKAAWGKVGGHAAEYGAEALERMFLSFPTTKTYFPHFDLSH"),
    ("P02662", "RPKHPIKHQGLPQEVLNENLLRFFVAPFPEVFGKEKVNELSKDIGSESTEDQAMEDIKQ"),
    ("P02663", "KNTMEHVSSSEESIISQETYKQEKNMAINPSKENLCSTFCKEVVRNAN"),
    ("P02666", "RELEELNVPGEIVESLSSSEESITRINKKIEKFQSEEQQQTEDELQDKIHPFAQTQSLV"),
    ("P02754", "MKCLLLALALTCGAQALIVTQTMKGLDIQKVAGTWYSLAMAASDISLLDAQSAPLRVYVE"
               "ELKPTPEGDLEILLQKWENGECAQKKIIAE"),
    ("P00711", "KQFTKCELSQLLKDIDGYGGIALPELICTMFHTSGYDTQAIVENNESTEYGLFQISNKLWC"
               "KSSQVPQSRNICDISCDKFL"),
    ("P02787", "VPDKTVRWCAVSEHEATKCQSFRDHMKSVIPSDGPSVACVKKASYLDCIRAIAANEADAV"
               "TLDAGLVYDAYLAPNNLKPVVAEF"),
    ("P01012", "GSIGAASMEFCFDVFKELKVHHANENIFYCPIAIMSALAMVYLGAKDSTRTQINKVVRFD"
               "KLPGFGDSIEAQCGTSVNVHSSLRDILNQI"),
    ("P06702", "MTCKMSQLERNIETIINTFHQYSVKLGHPDTLNQGEFKELVRKDLQNFLKKENKNEKVI"
               "EHIMEDLDTNADKQLSFEEFIMLMARLTWASH"),
    ("P62894", "MAAAKKAVDKIKKLFLKFPEVKNEDLGAQTMFNLFDKPQSAGLCGAGGRPVLAG"),
    ("Q32MB2", "MKTFFIFTLTLAISATSAQQNNPFIFNEKYTMVSVLSKDPNCNKVVIGTDTQQYYSNAC"
               "GILLNCTGIDLFKDKPV"),
]


def get_builtin_contaminants(prefix: str = "CONT_") -> List[oms.FASTAEntry]:
    """Return a list of built-in contaminant FASTA entries with the given prefix."""
    entries = []
    for acc, seq in BUILTIN_CONTAMINANTS:
        e = oms.FASTAEntry()
        e.identifier = f"{prefix}{acc}"
        e.sequence = seq.replace(" ", "")
        e.description = f"Contaminant {acc}"
        entries.append(e)
    return entries


def load_fasta(input_path: str) -> List[oms.FASTAEntry]:
    """Load entries from a FASTA file."""
    entries = []
    fasta_file = oms.FASTAFile()
    fasta_file.load(input_path, entries)
    return entries


def save_fasta(entries: List[oms.FASTAEntry], output_path: str) -> None:
    """Save entries to a FASTA file."""
    fasta_file = oms.FASTAFile()
    fasta_file.store(output_path, entries)


def deduplicate(entries: List[oms.FASTAEntry]) -> List[oms.FASTAEntry]:
    """Remove duplicate entries based on identifier."""
    seen = set()
    result = []
    for entry in entries:
        key = entry.identifier.split()[0]
        if key not in seen:
            seen.add(key)
            result.append(entry)
    return result


def merge_contaminants(
    input_path: str,
    output_path: str,
    add_crap: bool = True,
    contaminants_path: str = None,
    prefix: str = "CONT_",
) -> dict:
    """Merge contaminant sequences into a target FASTA file.

    Returns statistics about the merge.
    """
    target_entries = load_fasta(input_path)
    contaminant_entries = []

    if add_crap:
        contaminant_entries.extend(get_builtin_contaminants(prefix))

    if contaminants_path:
        custom = load_fasta(contaminants_path)
        for entry in custom:
            if not entry.identifier.startswith(prefix):
                entry.identifier = f"{prefix}{entry.identifier}"
            contaminant_entries.append(entry)

    merged = target_entries + contaminant_entries
    deduped = deduplicate(merged)
    save_fasta(deduped, output_path)

    return {
        "target_count": len(target_entries),
        "contaminant_count": len(contaminant_entries),
        "merged_count": len(merged),
        "deduplicated_count": len(deduped),
    }


@click.command(help="Append contaminant proteins to a FASTA database.")
@click.option("--input", "input", required=True, help="Input target FASTA file")
@click.option("--add-crap", is_flag=True, help="Add built-in cRAP contaminants")
@click.option("--contaminants", default=None, help="Custom contaminant FASTA file")
@click.option("--prefix", default="CONT_", help="Prefix for contaminant accessions (default: CONT_)")
@click.option("--output", required=True, help="Output merged FASTA file")
def main(input, add_crap, contaminants, prefix, output) -> None:
    if not add_crap and not contaminants:
        raise click.UsageError("At least one of --add-crap or --contaminants is required.")

    stats = merge_contaminants(input, output, add_crap, contaminants, prefix)
    print(f"Target: {stats['target_count']}, Contaminants: {stats['contaminant_count']}, "
          f"Merged (dedup): {stats['deduplicated_count']}")


if __name__ == "__main__":
    main()
