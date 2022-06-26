from promoter_extraction.check_distance_between_genes import get_genes_with_potential_promoter

from promoter_extraction.sequence_io import parse_gb_file, read_genome_file


def prom_extract_multi(in_gb: str, prom_len: int, file_out: str, full_genome: str, distance_threshold: int = 300):
    """
    Function to extract sequences upstream of the gene (promoters)

    :param in_gb: file path of the genome genbank file
    :param prom_len: desired promoter length
    :param file_out: name of the output file
    :param full_genome: complete genome file
    :return str:
    """

    prom_out = ""

    gb_record = parse_gb_file(in_gb)

    # get the genes with potential promoters
    genes_with_potential_promoter = get_genes_with_potential_promoter(gb_record, distance_threshold=distance_threshold)

    for feature in gb_record.features:
        if feature.type == "gene":
            gene_start = feature.location.start.position  # identifies the start position of the gene
            gene_end = feature.location.end.position  # identifies the end position of the gene
            promoter_start = gene_start - prom_len #identifies the start of the promoter
            promoter_end = gene_end + prom_len #identifies the end of the promoter

            if feature.qualifiers.get('locus_tag')[0] in genes_with_potential_promoter:

                if feature.strand == -1:
                    feat_loc = str(feature.location)
                    potential_promoter = full_genome[gene_end:promoter_end]
                    current_promoter = potential_promoter.reverse_complement()
                    if 'db_xref' not in feature.qualifiers.keys():
                        prom_out += "Promoter reverse complement" + ";" + \
                                    str(feature.qualifiers.get('locus_tag'))[2:-2] + ";" + feat_loc + ";"
                    else:
                        prom_out += "Promoter reverse complement" + ";" + \
                                    str(feature.qualifiers.get('locus_tag'))[2:-2] + "___" + \
                                    feature.qualifiers['db_xref'][0] + ";" + feat_loc + ";"
                    prom_out += str(current_promoter.seq) + "\n"

                elif feature.strand == 1:
                    feat_loc = str(feature.location)
                    current_promoter = full_genome[promoter_start:gene_start]
                    if 'db_xref' not in feature.qualifiers.keys():
                        prom_out += "Promoter" + ";" + str(feature.qualifiers.get('locus_tag'))[2:-2] + ";" + feat_loc + ";"
                    else:
                        prom_out += "Promoter" + ";" + str(feature.qualifiers.get('locus_tag'))[2:-2] + "___" + \
                                    feature.qualifiers['db_xref'][0] + ";" + feat_loc + ";"
                    prom_out += str(current_promoter.seq) + "\n"

    file = open(file_out, 'w')
    file.write(prom_out)
    file.close()


if __name__ == "__main__":
    fasta_full_genome = read_genome_file("s_thermophilus_complete_genome.fasta")
    prom_extract_multi(in_gb = "s_thermophilus_genome_annotation.gb",
                       prom_len = 300,
                       file_out = "filtered_promoters_300_bp.fasta",
                       full_genome = fasta_full_genome)
