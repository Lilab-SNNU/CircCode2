from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import yaml
import argparse
import time
import os

def make_junction_seq(raw_circrnas, circrna_gtf, raw_reads, tmp_file_path, final_clean_fastq_reads):

    print(get_time(), 'make junction seq...')


    max_reads_lenth = 50


    # get the biggest length of reads, assemble the junction file 
    global junc_site
    junc_site = max_reads_lenth
    print(junc_site)

    read_name = str(raw_reads).split('/')[-1].split('.')[0]
    circ_name = read_name + '_' + str(raw_circrnas).split('/')[-1].split('.')[0]

    han_raw_circ = open(raw_circrnas)

    global true_junc_file
    true_junc_file = tmp_file_path + '/' + circ_name + '_true_junc.fa'
    data_junc = []

    # make the true junction file
    for seq_record in SeqIO.parse(han_raw_circ, 'fasta'):
        seq_true_junc = seq_record.seq[-max_reads_lenth:]+seq_record.seq[:max_reads_lenth]
        rec_junc_true = SeqRecord(seq_true_junc,
                                  id = seq_record.id,
                                  description = '')
        data_junc.append(rec_junc_true)
    han_junc = open(true_junc_file, 'w+')
    SeqIO.write(data_junc, han_junc, 'fasta')

    han_junc.close()
    han_raw_circ.close()

def make_junc_index_align(raw_reads, raw_circrnas, thread, tmp_file_path, soft_log_folder, final_clean_fastq_reads):

    read_name = str(raw_reads).split('/')[-1].split('.')[0]
    circ_name = read_name + '_' + str(raw_circrnas).split('/')[-1].split('.')[0]

    # make junction index
    print(get_time(), 'make junction index...')
    junc_index_log = soft_log_folder + '/' + read_name +'_junc_index_log'
    junc_index = tmp_file_path+'/'+circ_name+'_true_junc_index/' + circ_name+'_true_junc_index'
    junc_index_folder = tmp_file_path+'/'+circ_name+'_true_junc_index'
    if not os.path.exists(junc_index_folder):
         os.makedirs(junc_index_folder)
    os.system('./requiredSoft/STAR --runMode genomeGenerate --runThreadN {} --genomeDir {} --genomeFastaFiles {} --limitGenomeGenerateRAM 98419370026 --genomeSAindexNbases 10 > {}'.format(thread, tmp_file_path+'/'+circ_name+'_true_junc_index/', true_junc_file, junc_index_log))
    

    # map the junction to filter reads
    print(get_time(), 'map the filter reads to the junction of circRNA...')
    global out_name
    out_name = tmp_file_path +'/'+ circ_name
    junc_map_reads_log = soft_log_folder + '/' + read_name +'_junc_map_to_reads_log'
    os.system('./requiredSoft/STAR --runThreadN {} --outSAMtype BAM SortedByCoordinate --genomeDir {} --readFilesIn {} --outFileNamePrefix {} > {}'.format(thread, tmp_file_path+'/'+circ_name+'_true_junc_index/', final_clean_fastq_reads, out_name, junc_map_reads_log))

def bamtobed_fq():

    print(get_time(), 'transform the bamflie to bedfile and fastq')
    bamfile = out_name + 'Aligned.sortedByCoord.out.bam'
    global bedfile
    bedfile = bamfile + '.bamtobedresult.bed'
    os.system('bedtools bamtobed -cigar -i {} > {}'.format(bamfile, bedfile))
    global align_fq_file
    align_fq_file = bamfile + '.align_bamtofq.fastq'
    os.system('bedtools bamtofastq -i {} -fq {}'.format(bamfile, align_fq_file))

def filter_coding_circ_and_reads_index(raw_circrnas, raw_reads, tmp_file_path, reads_junc_num):

    reads_lenth_dic = {}

    with open(align_fq_file) as han_fastq_reads:
         final_clean_reads_record = SeqIO.parse(han_fastq_reads, 'fastq')
         for seq_record in final_clean_reads_record:
             reads_lenth_dic[seq_record.id] = len(seq_record.seq)

    print(get_time(), 'filter coding ablity by ribo-seq data...')
    coding_circ_reads_site = {}
    han_bed = open(bedfile)
    #han_not = open('align.txt','w+')    

    # save the circ which reads cross the junction, save the reads site
    for line in han_bed:
        list_line = line.split('\t')
        circ_id = list_line[0]
        read_id = list_line[3]
        reads_strand = list_line[5]
        cigar_value = list_line[6]
        # refer genome true site, junc_site = the longest read lenth
        start_site = eval(list_line[1]) + 1
        stop_site = eval(list_line[2])
        # reads true site
        for i in cigar_value:
            if i == 'S':
                mov_num = -int(cigar_value[:cigar_value.index(i)])
                break
            elif i == 'M':
                mov_num = 0
                break
        reads_start = start_site + mov_num
        reads_stop = reads_start + reads_lenth_dic[read_id] - 1
        if start_site <= junc_site and stop_site > junc_site and reads_strand == '+' and 'N' not in cigar_value and 1 >= (stop_site-start_site+1)/(reads_lenth_dic[read_id]) > 0.9:
             
            if (stop_site-start_site)/(reads_lenth_dic[read_id])>=1:
                han_not.write(circ_id+'\t'+read_id+'\n')
            if circ_id in coding_circ_reads_site.keys():
                # true site / junc_site = left(-1,-2,-3...), right(1,2,3...)
                coding_circ_reads_site[circ_id] += [read_id, reads_start-junc_site-1, reads_stop-junc_site]
            if circ_id not in coding_circ_reads_site.keys():
                coding_circ_reads_site[circ_id] = [read_id, reads_start-junc_site-1, reads_stop-junc_site]

    # filter by the num of reads cross junction
    coding_circ_list = []
    # write the reads site to new file
    read_name = str(raw_reads).split('/')[-1].split('.')[0]
    reads_cover_jun_file = tmp_file_path + '/' + read_name + '_read_cover_jun.txt'
    han_reads_cover = open(reads_cover_jun_file, 'w+')
    for key in coding_circ_reads_site.keys():
        if int(len(coding_circ_reads_site[key])/3) >= reads_junc_num:
            coding_circ_list.append(key)
        han_reads_cover.write(key + '\t' +str(int(len(coding_circ_reads_site[key])/3)))
        for item in coding_circ_reads_site[key]:
            han_reads_cover.write('\t'+str(item))
        han_reads_cover.write('\n')

    # write the circRNA with ribo_seq data to new file
    han_raw_circ = open(raw_circrnas)
    read_name = str(raw_reads).split('/')[-1].split('.')[0]
    coding_circ_file = tmp_file_path + '/' + read_name + '_true.fa'
    coding_circ_seq_list = []
    for seq_record in SeqIO.parse(han_raw_circ, 'fasta'):
        if seq_record.id in coding_circ_list:
            junc_reads_count = str(int(len(coding_circ_reads_site[seq_record.id])/3))
            rec_coding = SeqRecord(seq_record.seq,
                            id=seq_record.id+'$'+junc_reads_count,
                            description='')
            coding_circ_seq_list.append(rec_coding)
    han_coding_circ = open(coding_circ_file,'w+')
    SeqIO.write(coding_circ_seq_list, han_coding_circ, 'fasta')
    han_bed.close()
    han_reads_cover.close()
    han_coding_circ.close()
    han_not.close()

def get_time():
    nowtime = time.strftime('[%Y-%m-%d %H:%M:%S]', time.localtime())
    return nowtime

def main():
    parse = argparse.ArgumentParser(description='This script helps to make junction index for each circ!')
    parse.add_argument('-y', '--yaml', required=True, help='please input the config.yaml file')
    args = parse.parse_args()

    yamlfile = args.yaml
    con_file = open(yamlfile)
    fileload = yaml.full_load(con_file)

    raw_reads = fileload['raw_reads']
    transcript_fasta = fileload['transcript_fasta']
    ribosome_fasta = fileload['ribosome_fasta']
    tmp_file_path = fileload['tmp_file_location']
    thread = fileload['thread']
    ribotype = fileload['ribotype']
    raw_circrnas = fileload['circrnas']
    circrna_gtf = fileload['circrna_gtf']

    reads_junc_num = 5

    if not os.path.exists(tmp_file_path):
         os.makedirs(tmp_file_path)     

    soft_log_folder = tmp_file_path + '/' + 'log'

    if not os.path.exists(soft_log_folder):
         os.makedirs(soft_log_folder)

    for item in raw_reads:
        raw_reads = item
        read_name = str(raw_reads).split('/')[-1].split('.')[0]

        final_clean_fastq_reads = tmp_file_path + '/' + read_name + '_final_clean.fastq'

        make_junction_seq(raw_circrnas, circrna_gtf, raw_reads, tmp_file_path, final_clean_fastq_reads)
        make_junc_index_align(raw_reads, raw_circrnas, thread, tmp_file_path, soft_log_folder, final_clean_fastq_reads)
        bamtobed_fq()
        filter_coding_circ_and_reads_index(raw_circrnas, raw_reads, tmp_file_path, reads_junc_num)
    print(get_time(), 'filter_coding_reads_circ finished!!!')

if __name__ == '__main__':
    main()

