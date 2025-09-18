__author__ = 'jayna'

from typing import Dict, List, Any

import Bio.SeqRecord
from Bio import SeqIO, Seq
from datetime import datetime, date, time
import os

month_to_num = {'Jan': 1, 'Feb': 2, 'Mar': 3, 'Apr': 4, 'May': 5, 'Jun': 6, 'Jul': 7, 'Aug': 8, 'Sep': 9, 'Oct': 10,
                'Nov': 11, 'Dec': 12}

# gb_file = "/Users/jayna/Downloads/sequence-5.gb"
# test_fastoutfile = "/Users/jayna/Downloads/sequence-3.out.fasta"
#
#
# with open('hantavirus_segment_dict.txt', 'r') as file:
#     data = file.read()
#
# hantavirus_segments = eval(data)
#
# hantavirus_min_lengths = {"S": 845, "M": 2200, "L": 4200}

### TAXONOMY LIST

filename = "/Users/jayna/Evolve.Zoo Dropbox/Jayna Raghwani/PycharmProjects/processGenbank/probePanel/taxonomy_result.txt"


def gb_to_fasta__henipavirus(gb_file, out_fastafile):
    ### GENBANK FILE ###
    # create a variable to read your genbank file (i.e. gb_file). Note "rU", which indicates the file is readable.
    input_handle = open(gb_file, "rU")
    count = 0

    ### OUTPUT FASTA FILE ###
    # create a variable to write data to your output file (i.e. out_fastafile). Note "w", which indicates the file is writable.

    G_output_handle = open(out_fastafile[0], "w")
    F_output_handle = open(out_fastafile[1], "w")

    ### METADATA FILE ###
    # out_tablefile is a string variable that corresponds to the filename of the metadata file. Here I am using str_replace so that the filename is similar to the output fasta filena,e
    out_tablefile = out_fastafile[0].replace(".fasta", "_metadata.csv")

    # Similar to output_handle, create a variable so that you can write to the metadata file.
    out_tablefile_handle = open(out_tablefile, "w")

    # Write out the column names to the metadata file.
    out_tablefile_handle.write(
        "Sequence Name,Accession_no,Strain,Organism,Host,Country,Collection_date,Year,Sequence_length,Pubmed_id\n")

    ### PROCESS GENBANK FILE ###

    for index, record in enumerate(SeqIO.parse(input_handle, "genbank")):

        # extract genus of virus sequence
        taxonomy = record.annotations["taxonomy"]
        genus = taxonomy[len(taxonomy) - 1]
        print(genus)

        glycoprotein_product = []
        fusion_product = []

        # extract sequence corresponding to coding region
        for feature in record.features:

            if feature.type == "CDS":  # this selects only the coding region of the genbank file

                seq = feature.location.extract(
                    record).seq  # this extracts nucleotide sequence (should be in frame, i.e. no stop codons)

                # get the length of nucleotide sequence and store this value in seq_length variable
                seq_length = len(seq)

                # get the name of the product encoded by the CDS (will be polyprotein for west nile virus, but for other viruses could have multiple products)
                gene_product = feature.qualifiers["product"][0]

                # get the protein accession number and store it in the protein_id variable.
                # This corresponds to the translated amino acid sequence of the nucleotide sequence
                protein_id = feature.qualifiers["protein_id"][0]

                # prints out the four variables to screen
                # print(seq_length, gene_product, protein_id)
                if (gene_product in ["attachment glycoprotein", "G protein", "cell attachment protein"]):
                    glycoprotein_product = seq

                if (gene_product in ["fusion protein", "F protein"]):
                    fusion_product = seq

            # print(feature.qualifiers["protein_id"])
            # seq = feature.location.extract(record).seq

        # For each record in the genbank file, extracts the different subfields (or qualifiers) in the FEATURES field
        features_qualifiers = record.features[0].qualifiers

        # strain_name = record.annotations['accessions'][0]

        # Metadata / fields we are interested in extracting per sequence
        d = ""
        name = ""
        year = ""
        host = "NA"
        isolate = "NA"
        country = "NA"
        t_date = "NA"
        accession_no = record.annotations['accessions'][0]
        organism = features_qualifiers['organism'][0]

        # Get the Pubmed id of the publication associated with the record/sequence
        reference = record.annotations['references'][0]
        pubmed_id = reference.pubmed_id

        # Extracting which host the sequence was obtained from. This could be stored in /host= field or /isolation source
        if 'host' in features_qualifiers.keys():

            host = features_qualifiers['host'][0]
            host = str(host).replace(' ', ':')
            # [optional] print out host to screen to see how the above line changes the output

        elif 'isolation_source' in features_qualifiers.keys():
            host = features_qualifiers['isolation_source'][0]
            host = str(host).replace(' ', ':')

        # Extract the country information of where the sequence was sampled from
        if 'geo_loc_name' in features_qualifiers.keys():
            country = features_qualifiers['geo_loc_name'][0]
            country_parts = country.split(':')
            country = country_parts[0]
            country = country.replace(" ", "")

        # Extract the date of collection of the sequence
        if 'collection_date' in features_qualifiers.keys():
            d = features_qualifiers['collection_date'][0]

        # Extract isolate name
        if 'strain' in features_qualifiers.keys():
            isolate = features_qualifiers['strain'][0]

        # The Date field may need further processing.
        # If the date string is longer than 4 characters than in addition to year, likely contains information about month and day of sampling
        if (len(d) > 4):

            # splits date string by "-"
            parts = d.split('-')

            # I would suggest printing out the parts variable above to better understand what the above line does

            # This section undertakes further processing to identify which components correspond to month and year parts of the date,
            # and use this information to convert the string date e.g. 2002-01-03 or Dec-2001 into decimalised date e.g. 2002.013
            if (len(parts) == 2):

                # print time.strptime(d, "%b %y")
                # print record.features[2].qualifiers['isolate'][0]

                if parts[0].isdigit() == False:

                    # parts[0] must contain a string month.
                    # parts[1] must contain a integer year.
                    # no day information so assume date is in middle of the month.

                    year = int(parts[1])
                    d_date = datetime(int(parts[1]), month_to_num[parts[0]], 15)
                    days_in_year = datetime.strftime(d_date, '%j')
                    total_datetime = datetime(int(parts[1]), 12, 31)
                    total_days_in_year = date.strftime(total_datetime, '%j')
                    decimal_date = float(parts[1]) + float(days_in_year) / float(total_days_in_year)

                elif parts[0].isdigit():
                    # parts[0] must contain a string year.
                    year = int(parts[0])
                    d_date = datetime(int(parts[0]), int(parts[1]), 15)
                    days_in_year = datetime.strftime(d_date, '%j')
                    total_datetime = datetime(int(parts[0]), 12, 31)
                    total_days_in_year = date.strftime(total_datetime, '%j')
                    decimal_date = float(parts[0]) + float(days_in_year) / float(total_days_in_year)

                else:
                    # parts[1] must contain a string year.
                    year = int(parts[1])
                    d_date = datetime(int(parts[1]), month_to_num[parts[0]], 15)
                    days_in_year = datetime.strftime(d_date, '%j')
                    total_datetime = datetime(int(parts[1]), 12, 31)
                    total_days_in_year = date.strftime(total_datetime, '%j')
                    decimal_date = float(parts[1]) + float(days_in_year) / float(total_days_in_year)

                t_date = str('%.4f' % decimal_date)
            elif (len(parts) == 1):
                year = parts[0]
            else:

                if parts[1].isdigit():
                    # parts[0] must contain a integer year.
                    # parts[1] must contain a integer month.
                    # parts[2] must contain a integer day.

                    year = int(parts[0])

                    d_date = datetime(int(parts[0]), int(parts[1]), int(parts[2]))
                    days_in_year = datetime.strftime(d_date, '%j')
                    total_datetime = datetime(int(parts[0]), 12, 31)
                    total_days_in_year = date.strftime(total_datetime, '%j')
                    decimal_date = float(parts[0]) + float(days_in_year) / float(total_days_in_year)

                    t_date = str('%.4f' % decimal_date)

                else:
                    # parts[0] must contain a integer day.
                    # parts[1] must contain a string month.
                    # parts[2] must contain a integer year.

                    year = int(parts[2])
                    d_date = datetime(int(parts[2]), month_to_num[parts[1]], int(parts[0]))
                    days_in_year = datetime.strftime(d_date, '%j')
                    total_datetime = datetime(int(parts[2]), 12, 31)
                    total_days_in_year = date.strftime(total_datetime, '%j')
                    decimal_date = float(parts[2]) + float(days_in_year) / float(total_days_in_year)

                    t_date = str('%.4f' % decimal_date)

                    # print t_date
                    # print ">"+name
                    # print record.seq

        # if the date is already a number then it must correspond to year
        elif (is_number(d) == True):
            t_date = float(d)
            year = t_date

        # Filtering sequences by sequence length. Here we have chosen 9000 bp and above
        # We only write out these sequences to fasta file
        # if (gene_product):
        #

        name = accession_no + "_" + country + "_" + host + "_" + str(t_date) + "|" + str.replace(genus, " ", "..")
        G_output_handle.write(">" + name + "\n")
        F_output_handle.write(">" + name + "\n")

        G_output_handle.write(str(glycoprotein_product) + "\n")
        F_output_handle.write(str(fusion_product) + "\n")

        G_output_handle.flush()
        F_output_handle.flush()

        count += 1

        # This line writes out all the metadata collected for each record to the metadata output file (out_tablefile)
        out_tablefile_handle.write(name + "," +
                                   accession_no + "," + isolate + "," + organism + "," + host + "," + country + "," + str(
            t_date) + "," + str(year) + "," + str(
            seq_length) + "," + str(pubmed_id) + "\n")

    # close all the files once finished processing the genbank file.
    input_handle.close()
    out_tablefile_handle.close()
    G_output_handle.close()
    F_output_handle.close()


def process_gb(gb_file, working_dir, virus_family, min_length, host_taxafile):
    ### GENBANK FILE ###
    # create a variable to read your genbank file (i.e. gb_file). Note "rU", which indicates the file is readable.
    input_handle = open(gb_file, "rU")
    count = 0

    ### METADATA FILE ###
    # out_tablefile is a string variable that corresponds to the filename of the metadata file. Here I am using str_replace so that the filename is similar to the output fasta filena,e

    out_tablefile = working_dir + "/" + virus_family + "_metadata.csv"

    try:
        os.makedirs(working_dir + "/" + virus_family)

    except FileExistsError:
        # directory already exists
        pass

    # Similar to output_handle, create a variable so that you can write to the metadata file.
    out_tablefile_handle = open(out_tablefile, "w")

    # Write out the column names to the metadata file.
    out_tablefile_handle.write(
        "Sequence Name,Accession_no,Strain,Organism,Species/Genus,Host,Country,Collection_date,Year,Sequence_length,Pubmed_id\n")

    with open(host_taxafile) as f:
        rodentia_taxa = f.read().splitlines()

    ### PROCESS GENBANK FILE ###

    species_dict = {}

    for index, record in enumerate(SeqIO.parse(input_handle, "genbank")):

        # extract genus of virus sequence
        taxonomy = record.annotations["taxonomy"]

        if (len(taxonomy) >= 9):
            species = taxonomy[8]
        else:
            species = taxonomy[len(taxonomy) - 1]

        # For each record in the genbank file, extracts the different subfields (or qualifiers) in the FEATURES field
        features_qualifiers = record.features[0].qualifiers

        # strain_name = record.annotations['accessions'][0]

        # Metadata / fields we are interested in extracting per sequence
        d = ""
        name = ""
        year = ""
        host = "NA"
        isolate = "NA"
        country = "NA"
        t_date = "NA"
        accession_no = record.annotations['accessions'][0]
        organism = features_qualifiers['organism'][0]

        # Get the Pubmed id of the publication associated with the record/sequence
        reference = record.annotations['references'][0]
        pubmed_id = reference.pubmed_id

        # Extracting which host the sequence was obtained from. This could be stored in /host= field or /isolation source
        if 'host' in features_qualifiers.keys():

            host = features_qualifiers['host'][0]

            if (host in rodentia_taxa):
                host = str(host).replace(' ', ':')
            else:
                continue

            # [optional] print out host to screen to see how the above line changes the output

        elif 'isolation_source' in features_qualifiers.keys():
            host = features_qualifiers['isolation_source'][0]

            if (host in rodentia_taxa):
                host = str(host).replace(' ', ':')
            else:
                continue

        if 'Homo' in host:
            continue

        if "strain" in host or "NA" in host:
            continue

        # Extract the country information of where the sequence was sampled from
        if 'geo_loc_name' in features_qualifiers.keys():
            country = features_qualifiers['geo_loc_name'][0]
            country_parts = country.split(':')
            country = country_parts[0]
            country = country.replace(" ", ":")

        # Extract the date of collection of the sequence
        if 'collection_date' in features_qualifiers.keys():
            d = features_qualifiers['collection_date'][0]

        # Extract isolate name
        if 'strain' in features_qualifiers.keys():
            isolate = features_qualifiers['strain'][0]

        if 'note' in features_qualifiers.keys():
            notes = features_qualifiers['note'][0]

            if 'experiment' in notes or 'vaccine' in notes:
                continue

        if 'lab_host' in features_qualifiers.keys():
            continue

        # The Date field may need further processing.
        t_date, year = process_date(d)

        name = accession_no + "_" + country + "_" + host + "|" + str(t_date) + "|" + str.replace(species, " ",
                                                                                                 "..") + "|" + str.replace(
            organism, " ", ":")

        count += 1

        tmp_min_length = min_length

        for f in record.features:

            if (f.type == "CDS"):

                seq = f.location.extract(record).seq

                if (len(str(seq)) > tmp_min_length):
                    tmp_min_length = len(str(seq))
                    record.seq = seq

        N_percentage = str(record.seq).count('N') / len(str(record.seq))

        if (N_percentage >= 0.1):
            continue

        new_record = Bio.SeqRecord.SeqRecord(name=name, seq=record.seq)

        if (species not in species_dict.keys()):
            species_dict[species] = []
            species_dict[species].append(new_record)

            # This line writes out all the metadata collected for each record to the metadata output file (out_tablefile)
            out_tablefile_handle.write(
                name + "," +
                accession_no + "," +
                isolate + "," +
                organism + "," +
                species + "," +
                host + "," +
                country.replace(":", " ") + "," +
                str(t_date) + "," +
                str(year) + "," +
                str(len(record.seq)) + "," +
                str(pubmed_id) + "\n")

            out_tablefile_handle.flush()

        else:
            species_dict[species].append(new_record)

            # This line writes out all the metadata collected for each record to the metadata output file (out_tablefile)
            out_tablefile_handle.write(
                name + "," +
                accession_no + "," +
                isolate + "," +
                organism + "," +
                species + "," +
                host.replace(":", " ") + "," +
                country.replace(":", " ") + "," +
                str(t_date) + "," +
                str(year) + "," +
                str(len(record.seq)) + "," +
                str(pubmed_id) + "\n")

            out_tablefile_handle.flush()

    # Write out sequences to file
    try:
        os.makedirs(working_dir + "/" + virus_family)

    except FileExistsError:
        # directory already exists
        pass

    for s in species_dict.keys():

        species_name = (s.replace(" ", "_"))

        ### OUTPUT FASTA FILE ###
        # create a variable to write data to your output file (i.e. out_fastafile). Note "w", which indicates the file is writable.

        fastoutfile = working_dir + "/" + virus_family + "/" + str(species_name) + ".fasta"

        out_fastoutfile_handle = open(fastoutfile, "w")

        seq_records = species_dict[s]

        for record in seq_records:
            out_fastoutfile_handle.write(">" + str(record.name) + "\n")
            out_fastoutfile_handle.write(str(record.seq) + "\n")

        out_fastoutfile_handle.close()

    # close all the files once finished processing the genbank file.

    input_handle.close()
    out_tablefile_handle.close()


def process_gb_seg(gb_file, working_dir, virus_family, segment_length, segments, host_taxafile):
    product_gene_list = []
    ### GENBANK FILE ###
    # create a variable to read your genbank file (i.e. gb_file). Note "rU", which indicates the file is readable.
    input_handle = open(gb_file, "rU")
    count = 0

    ### METADATA FILE ###
    # out_tablefile is a string variable that corresponds to the filename of the metadata file. Here I am using str_replace so that the filename is similar to the output fasta filena,e
    out_tablefile = working_dir + "/" + virus_family + "_metadata.csv"

    try:
        os.makedirs(working_dir + "/" + virus_family)

    except FileExistsError:
        # directory already exists
        pass

    # Similar to output_handle, create a variable so that you can write to the metadata file.
    out_tablefile_handle = open(out_tablefile, "w")

    # Write out the column names to the metadata file.
    out_tablefile_handle.write(
        "Sequence Name,Accession_no,Segment,Strain,Organism,Species/Genus,Host,Country,Collection_date,Year,Sequence_length,Pubmed_id\n")

    with open(host_taxafile) as f:
        rodentia_taxa = f.read().splitlines()

    ### PROCESS GENBANK FILE ###

    species_dict_by_seg = {}

    for seg in segments.keys():
        species_dict_by_seg[seg] = {}

    for index, record in enumerate(SeqIO.parse(input_handle, "genbank")):

        # extract genus of virus sequence
        taxonomy = record.annotations["taxonomy"]

        # print(len(taxonomy), taxonomy)

        species_or_genus = "NA"
        taxonomy_length = len(taxonomy)
        if taxonomy_length == 8:
            species_or_genus = "unclassified " + taxonomy[7]

        elif taxonomy_length == 9:
            species_or_genus = taxonomy[8]

        elif taxonomy_length == 10:

            if "Orthohantavirus" in record.features[0].qualifiers['organism'][0]:
                species_or_genus = record.features[0].qualifiers['organism'][0]
            else:
                species_or_genus = taxonomy[9]

        elif taxonomy_length == 11:
            species_or_genus = taxonomy[10]

        # For each record in the genbank file, extracts the different subfields (or qualifiers) in the FEATURES field
        features_qualifiers = record.features[0].qualifiers

        # strain_name = record.annotations['accessions'][0]

        # Metadata / fields we are interested in extracting per sequence
        d = ""
        name = ""
        year = ""
        host = "NA"
        isolate = "NA"
        country = "NA"
        t_date = "NA"
        accession_no = record.annotations['accessions'][0]
        organism = features_qualifiers['organism'][0]

        # Get the Pubmed id of the publication associated with the record/sequence
        reference = record.annotations['references'][0]
        pubmed_id = reference.pubmed_id

        # Extracting which host the sequence was obtained from. This could be stored in /host= field or /isolation source
        if 'host' in features_qualifiers.keys():

            host = features_qualifiers['host'][0]

            if (host in rodentia_taxa):
                host = str(host).replace(' ', ':')
            else:
                continue

            # [optional] print out host to screen to see how the above line changes the output

        elif 'isolation_source' in features_qualifiers.keys():
            host = features_qualifiers['isolation_source'][0]

            if (host in rodentia_taxa):
                host = str(host).replace(' ', ':')
            else:
                continue

        if 'Homo' in host:
            continue

        if "strain" in host or "NA" in host:
            continue

        # Extract the country information of where the sequence was sampled from
        if 'geo_loc_name' in features_qualifiers.keys():
            country = features_qualifiers['geo_loc_name'][0]
            country_parts = country.split(':')
            country = country_parts[0]
            country = country.replace(" ", ":")

        # Extract the date of collection of the sequence
        if 'collection_date' in features_qualifiers.keys():
            d = features_qualifiers['collection_date'][0]

        # Extract isolate name
        if 'strain' in features_qualifiers.keys():
            isolate = features_qualifiers['strain'][0]

        if 'note' in features_qualifiers.keys():
            notes = features_qualifiers['note'][0]

            if 'experiment' in notes or 'vaccine' in notes:
                continue

        if 'lab_host' in features_qualifiers.keys():
            continue

        # The Date field may need further processing.
        t_date, year = process_date(d)

        # species_or_genus = get_virus_taxonomic_rank(record, taxonomic_level)

        name = accession_no + "_" + country + "_" + host + "|" + str(t_date) + "|" + str.replace(species_or_genus, " ",
                                                                                                 "..") + "|" + str.replace(
            organism, " ", ":")

        count += 1

        segment_records = {}
        for seg in segments.keys():
            segment_records[seg] = None

        for f in record.features:
            if f.type == "CDS":
                # print(index, f.qualifiers.keys())

                product_or_gene = "NULL"
                if ("product" in f.qualifiers.keys()):

                    product_or_gene = f.qualifiers["product"][0]
                    product_gene_list.append(f.qualifiers["product"][0])

                elif ("gene" in f.qualifiers.keys()):
                    product_or_gene = f.qualifiers["gene"][0]
                    product_gene_list.append(f.qualifiers["gene"][0])

                if (product_or_gene in ["hypothetical protein", "polyprotein"]):

                    if ("segment" in record.features[0].qualifiers.keys()):
                        product_or_gene = record.features[0].qualifiers["segment"][0] + " protein"

                for seg in segments.keys():

                    tmp_min_length = segment_length[seg]

                    if product_or_gene in segments[seg]:

                        seq = f.location.extract(record).seq

                        if (len(str(seq)) > tmp_min_length):
                            tmp_min_length = len(str(seq))

                            N_percentage = str(seq).count('N') / len(str(seq))

                            if (N_percentage >= 0.1):
                                continue

                            else:

                                segment_records[seg] = Bio.SeqRecord.SeqRecord(name=name, seq=seq)

                        else:
                            continue

        #     new_record = Bio.SeqRecord.SeqRecord(name=name, seq=record.seq)
        #

        for seg in segments.keys():

            seg_record = segment_records[seg]

            if (seg_record is not None):

                if (species_or_genus not in species_dict_by_seg[seg].keys()):

                    species_dict_by_seg[seg][species_or_genus] = []
                    species_dict_by_seg[seg][species_or_genus].append(seg_record)

                    # This line writes out all the metadata collected for each record to the metadata output file (out_tablefile)
                    out_tablefile_handle.write(
                        name + "," +
                        accession_no + "," +
                        seg + "," +
                        isolate + "," +
                        organism + "," +
                        species_or_genus + "," +
                        host + "," +
                        country.replace(":", " ") + "," +
                        str(t_date) + "," +
                        str(year) + "," +
                        str(len(seg_record.seq)) + "," +
                        str(pubmed_id) + "\n")

                    out_tablefile_handle.flush()

                else:
                    species_dict_by_seg[seg][species_or_genus].append(seg_record)

                    # This line writes out all the metadata collected for each record to the metadata output file (out_tablefile)
                    out_tablefile_handle.write(
                        name + "," +
                        accession_no + "," +
                        seg + "," +
                        isolate + "," +
                        organism + "," +
                        species_or_genus + "," +
                        host.replace(":", " ") + "," +
                        country.replace(":", " ") + "," +
                        str(t_date) + "," +
                        str(year) + "," +
                        str(len(seg_record.seq)) + "," +
                        str(pubmed_id) + "\n")

                    out_tablefile_handle.flush()


    for seg in segments.keys():
        try:
            os.makedirs(working_dir + "/" + virus_family + "/" + seg)

        except FileExistsError:
            # directory already exists
            pass

        for s in species_dict_by_seg[seg].keys():

            species_name = (s.replace(" ", "_"))
            fastoutfile = working_dir + "/" + virus_family + "/" + seg + "/" + str(species_name) + ".fasta"

            out_fastoutfile_handle = open(fastoutfile, "w")

            seq_records = species_dict_by_seg[seg][s]

            for record in seq_records:
                out_fastoutfile_handle.write(">" + str(record.name) + "\n")
                out_fastoutfile_handle.write(str(record.seq) + "\n")

            out_fastoutfile_handle.close()

    # close all the files once finished processing the genbank file.

    input_handle.close()
    out_tablefile_handle.close()


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def process_date(d):
    t_date = "NA"
    year = "NA"

    # If the date string is longer than 4 characters than in addition to year, likely contains information about month and day of sampling
    if (len(d) > 4):

        # splits date string by "-"
        parts = d.split('-')

        # I would suggest printing out the parts variable above to better understand what the above line does

        # This section undertakes further processing to identify which components correspond to month and year parts of the date,
        # and use this information to convert the string date e.g. 2002-01-03 or Dec-2001 into decimalised date e.g. 2002.013
        if (len(parts) == 2):

            # print time.strptime(d, "%b %y")
            # print record.features[2].qualifiers['isolate'][0]

            if parts[0].isdigit() == False:

                # parts[0] must contain a string month.
                # parts[1] must contain a integer year.
                # no day information so assume date is in middle of the month.

                year = int(parts[1])
                d_date = datetime(int(parts[1]), month_to_num[parts[0]], 15)
                days_in_year = datetime.strftime(d_date, '%j')
                total_datetime = datetime(int(parts[1]), 12, 31)
                total_days_in_year = date.strftime(total_datetime, '%j')
                decimal_date = float(parts[1]) + float(days_in_year) / float(total_days_in_year)

            elif parts[0].isdigit():
                # parts[0] must contain a string year.
                year = int(parts[0])
                d_date = datetime(int(parts[0]), int(parts[1]), 15)
                days_in_year = datetime.strftime(d_date, '%j')
                total_datetime = datetime(int(parts[0]), 12, 31)
                total_days_in_year = date.strftime(total_datetime, '%j')
                decimal_date = float(parts[0]) + float(days_in_year) / float(total_days_in_year)

            else:
                # parts[1] must contain a string year.
                year = int(parts[1])
                d_date = datetime(int(parts[1]), month_to_num[parts[0]], 15)
                days_in_year = datetime.strftime(d_date, '%j')
                total_datetime = datetime(int(parts[1]), 12, 31)
                total_days_in_year = date.strftime(total_datetime, '%j')
                decimal_date = float(parts[1]) + float(days_in_year) / float(total_days_in_year)

            t_date = str('%.4f' % decimal_date)
        elif (len(parts) == 1):
            year = parts[0]
        else:
            # check if month is in a number or string format
            if parts[1].isdigit():
                # parts[0] must contain a integer year.
                # parts[1] must contain a integer month.
                # parts[2] must contain a integer day.

                year = int(parts[0])

                d_date = datetime(int(parts[0]), int(parts[1]), int(parts[2]))
                days_in_year = datetime.strftime(d_date, '%j')
                total_datetime = datetime(int(parts[0]), 12, 31)
                total_days_in_year = date.strftime(total_datetime, '%j')
                decimal_date = float(parts[0]) + float(days_in_year) / float(total_days_in_year)

                t_date = str('%.4f' % decimal_date)

            else:
                # parts[0] must contain a integer day.
                # parts[1] must contain a string month.
                # parts[2] must contain a integer year.

                if (parts[1] in month_to_num and is_number(parts[2]) == True):
                    year = int(parts[2])
                    d_date = datetime(int(parts[2]), month_to_num[parts[1]], int(parts[0]))
                    days_in_year = datetime.strftime(d_date, '%j')
                    total_datetime = datetime(int(parts[2]), 12, 31)
                    total_days_in_year = date.strftime(total_datetime, '%j')
                    decimal_date = float(parts[2]) + float(days_in_year) / float(total_days_in_year)

                    t_date = str('%.4f' % decimal_date)

                # print t_date
                # print ">"+name
                # print record.seq

    # if the date is already a number then it must correspond to year
    elif (is_number(d) == True):
        t_date = float(d)
        year = t_date

    return (t_date, year)
