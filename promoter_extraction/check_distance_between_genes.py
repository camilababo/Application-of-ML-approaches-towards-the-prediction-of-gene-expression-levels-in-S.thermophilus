from typing import Tuple, List

from Bio.SeqRecord import SeqRecord

from promoter_extraction.sequence_io import parse_gb_file


def sort_tuple_list(tuple_list: List[Tuple[SeqRecord, int, int]]) -> List[Tuple[SeqRecord, int, int]]:
    """
    Function to sort a list of tuples by the first element of the tuple
    :param tuple_list: list of tuples
    :return list:-
    """
    return sorted(tuple_list, key=lambda tup: tup[1])


def create_tuple_list_of_gene_coordinates(gb_record: SeqRecord) -> List[Tuple[SeqRecord, int, int]]:
    """
    Function to create a list of tuples containing the gene object, the gene start position and end position

    :param gb_record: Genbank file
    :return list:
    """

    tuple_list = []

    for feature in gb_record.features:
        if feature.type == "gene":  # if the feature is a gene
            sequence_start = feature.location.start.position
            sequence_end = feature.location.end.position  # identifies the start position of the gene

            tuple_list.append((feature, sequence_start, sequence_end))

    return tuple_list


def filter_potencial_promoters(tuple_list: List[Tuple[SeqRecord, int, int]],
                               distance_threshold: int) -> List[str]:
    """
    Function to filter the list of tuples by the gene start position.

    :param tuple_list: list of tuples
    :param distance_threshold: distance threshold
    :return list:
    """

    filtered_list = []

    for i in range(len(tuple_list)):
        if i + 1 < len(tuple_list):
            start_next_gene = tuple_list[i + 1][1]
            end_current_gene = tuple_list[i][2]
            distance = start_next_gene - end_current_gene
            if distance > distance_threshold:
                gene_potencial_with_promoter = tuple_list[i + 1][0]
                filtered_list.append(gene_potencial_with_promoter.qualifiers.get('locus_tag')[0])

    return filtered_list


def get_genes_with_potential_promoter(gb_record: SeqRecord, distance_threshold: int) -> List[str]:
    """
    Function to get the genes with potential promoters

    :return list:
    """

    tuple_list = create_tuple_list_of_gene_coordinates(gb_record)
    tuple_list = sort_tuple_list(tuple_list)
    filtered_list = filter_potencial_promoters(tuple_list, distance_threshold)

    return filtered_list


if __name__ == "__main__":
    seq_record = parse_gb_file("../sequence_st.gb")
    tuple_list = create_tuple_list_of_gene_coordinates(seq_record)
    sorted_tuple_list = sort_tuple_list(tuple_list)
    for i in [100, 200, 300, 400]:
        filtered_list = filter_potencial_promoters(sorted_tuple_list, distance_threshold=i)

        print(f"{len(filtered_list)} genes with a distance of {i} bp")
