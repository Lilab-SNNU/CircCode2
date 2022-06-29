from Bio import SeqIO
#from Bio.Seq import Seq
#from Bio.SeqRecord import SeqRecord
import yaml
import argparse
import time
import os


def predict_adapter(raw_reads, ribotype, tmp_file_path):

    read_name = str(raw_reads).split('/')[-1].split('.')[0]
    fastq_reads = raw_reads

    if ribotype == 'sra' or raw_reads.split('.')[-1] == 'sra':
        print(get_time(), '{} sra to fastq...'.format(read_name))
        os.system('fastq-dump {} -O {}'.format(raw_reads, tmp_file_path))
        fastq_reads = tmp_file_path + '/' + read_name + '.fastq'

    # obtain the sample of refer reads from raw reads
    print(get_time(), 'obtain refer reads from {}...'.format(read_name))
    
    raw_read_seq_num = 4000000
    refer_fastq_reads = tmp_file_path + '/' + read_name + '_refer.fastq'
    os.system('head -n {} {} > {}'.format(raw_read_seq_num, fastq_reads, refer_fastq_reads))

    # the num of these reads to be scaned
    refer_reads_num = 1000
    
    refer_reads_list = []
    seq_index = 0
    with open(refer_fastq_reads) as han_fastq:
        raw_read_seq_record = SeqIO.parse(han_fastq, 'fastq')
        for line in raw_read_seq_record:
            seq_index += 1
            if seq_index <= refer_reads_num:
                refer_reads_list.append(line.seq)
            elif seq_index > refer_reads_num:
                break

    # window scanning and predict adapters
    print(get_time(), 'analysis refer reads...')
    predict_adapter_window = 10
    adapter_start_rate_list = [0.8, 0.7, 0.6]
    
    adapter_start_pinlv_dic = {}
    global adapter_full_dic
    adapter_full_dic = {}

    for adapter_start_rate in adapter_start_rate_list:
        for str_read_seq in refer_reads_list:
            for site in range(0,len(str_read_seq)-predict_adapter_window):
                tmp_seq = str(str_read_seq[site:site+predict_adapter_window])
                tmp_count_file = os.popen('grep -c {} {}'.format(tmp_seq, refer_fastq_reads))
                tmp_seq_count = eval(tmp_count_file.read().replace('\n', ''))
                tmp_count_file.close()
                pinlv = tmp_seq_count/raw_read_seq_num*4
                if 1 >= pinlv >= adapter_start_rate:
                    adapter_start_pinlv_dic[tmp_seq] = pinlv
                    if tmp_seq not in adapter_full_dic.keys():
                        adapter_full_dic[tmp_seq] = str_read_seq[site:]
                    elif tmp_seq in adapter_full_dic.keys():
                        # find the full adapter for each adapter start
                        if len(str_read_seq[site:]) > len(adapter_full_dic[tmp_seq]):
                            adapter_full_dic[tmp_seq] = str_read_seq[site:]
                    break
        if len(adapter_full_dic) != 0:
            break

    if len(adapter_full_dic) == 0:
        xxx = 0
        #print(get_time(), 'no adapter in these reads')
    elif len(adapter_full_dic) != 0:
        # find the longest full adapter
        adapter_full_list = []
        not_pipei_list = []
        pipei_lenth = 10
        for value_full in adapter_full_dic.values():
            adapter_full_list.append(value_full)
        moban = max(adapter_full_list)
        for item in adapter_full_list:
            for i in range(len(item)-pipei_lenth):
                if moban.find(item[i:i+pipei_lenth]) == 0 and item[0:predict_adapter_window] not in moban:
                    moban = item[:i] + moban
                elif moban.find(item[i:i+pipei_lenth]) == len(moban)-pipei_lenth:
                    moban += item[i+pipei_lenth:]
                else:
                    not_pipei_list.append(item)
        global adapter_full_seq_list
        adapter_full_seq_list = []
        adapter_full_seq_list.append(moban)
        #print(get_time(), 'adapter_start in reads rate: {}'.format(adapter_start_pinlv_dic))
        #print(get_time(), 'adapter_full is: {}'.format(adapter_full_dic))
        #print(get_time(), 'adapter full seq is: {}'.format(moban))

def make_mrna_rrna_index(transcript_fasta, ribosome_fasta, thread, tmp_file_path, raw_reads, soft_log_folder):

    read_name = str(raw_reads).split('/')[-1].split('.')[0]
    transcript_name = str(transcript_fasta).split('/')[-1].split('.')[0]
    ribosome_name = str(ribosome_fasta).split('/')[-1].split('.')[0]

    print(get_time(), 'make rRNA index...')
    rRNA_index_log = soft_log_folder + '/' + read_name +'_rRNA_index_log'
    rRNA_index_folder = tmp_file_path+'/'+ribosome_name+'_index'
    if not os.path.exists(rRNA_index_folder):
         os.makedirs(rRNA_index_folder)

    os.system('bowtie-build --threads {} {} {}/{} > {}'.format(thread, ribosome_fasta, rRNA_index_folder, ribosome_name, rRNA_index_log))

    print(get_time(), 'make transcript index...')
    trans_index_log = soft_log_folder + '/' + read_name +'_transcprit_index_log'
    trans_index_folder = tmp_file_path+'/'+transcript_name+'_index'
    if not os.path.exists(trans_index_folder):
         os.makedirs(trans_index_folder)

    os.system('bowtie-build --threads {} {} {}/{} > {}'.format(thread, transcript_fasta, trans_index_folder, transcript_name, trans_index_log))

def deal_raw_reads(raw_reads, tmp_file_path, thread, transcript_fasta, ribosome_fasta, trimmomatic, ribotype, soft_log_folder):

    read_name = str(raw_reads).split('/')[-1].split('.')[0]

    if ribotype == 'sra':
        fastq_reads = tmp_file_path + '/' + read_name + '.fastq'
    else:
        fastq_reads = raw_reads

    # cut predicted adapters
    if len(adapter_full_dic) != 0:
        for item in adapter_full_seq_list:
            print(get_time(), 'cut predicted adapter {}...'.format(item))
            cut_adapter_window = 5
            adapter_command = ''
            for i in range(0,len(item)-cut_adapter_window+1):
                adapter_command += ' -a ' + item[i:i+cut_adapter_window]
            fastq_reads_cutadapt = tmp_file_path + '/' + read_name + '_cutadapt_adapter'
            cutadapt_log = soft_log_folder + '/' + read_name + '_cutadapt_log_adapter'

            ############## --quiet 
            os.system('cutadapt -j {}{} -m 16 -e 0.2 --trim-n -o {} {} > {}'.format(thread, adapter_command, fastq_reads_cutadapt+ '.fastq', fastq_reads, cutadapt_log))
                
    elif len(adapter_full_dic) == 0:
        # define double for no adapter trim
        fastq_reads_cutadapt = fastq_reads.split('.fastq')[0]

    # Filter out low quality reads by Trimmomatic
    print(get_time(), 'filter out low quality reads...')
    fastq_reads_cutadapt_clean = tmp_file_path + '/' + fastq_reads_cutadapt.split('/')[-1] + '_clean'
    trim_log = soft_log_folder + '/' + read_name + '_trim_log'
    os.system('java -jar {} SE -threads {} {} {} LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:16 > {}'.format(trimmomatic, thread, fastq_reads_cutadapt+'.fastq', fastq_reads_cutadapt_clean+'.fastq', trim_log))

    # Map clean reads to transcript sequence by bowtie
    print(get_time(), 'clean reads from transcript...')
    unmapped_transcript_reads = fastq_reads_cutadapt_clean + '_untranscipt'
    transcript_name = str(transcript_fasta).split('/')[-1].split('.')[0]
    transcript_index = tmp_file_path + '/' + transcript_name + '_index/' + transcript_name
    log_trans = soft_log_folder + '/' + transcript_name + '.map_to_transcript.sam'
    os.system('bowtie -p {} --un {} --norc {} {} > {}'.format(thread, unmapped_transcript_reads+'.fastq', transcript_index, fastq_reads_cutadapt_clean+'.fastq', log_trans))

    # Map clean reads to ribosome sequence by bowtie
    print(get_time(), 'clean reads from rRNA...')
    unmapped_rRNA_reads = tmp_file_path + '/' + read_name + '_final_clean'
    ribosome_name = str(ribosome_fasta).split('/')[-1].split('.')[0]
    rRNA_index = tmp_file_path + '/' + ribosome_name + '_index/' + ribosome_name
    log_out = soft_log_folder + '/' + ribosome_name + '.map_to_rRNA.sam'
    os.system('bowtie -p {} --un {} --norc {} {} > {}'.format(thread, unmapped_rRNA_reads+'.fastq', rRNA_index, unmapped_transcript_reads+'.fastq', log_out))
        
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
    #trimmomatic = fileload['trimmomatic_jar']
    thread = fileload['thread']
    ribotype = fileload['ribotype']
    raw_circrnas = fileload['circrnas']
    circrna_gtf = fileload['circrna_gtf']

    trimmomatic = './requiredSoft/trimmomatic-0.38.jar'

    if not os.path.exists(tmp_file_path):
         os.makedirs(tmp_file_path)     

    soft_log_folder = tmp_file_path + '/' + 'log'

    if not os.path.exists(soft_log_folder):
         os.makedirs(soft_log_folder)

    for item in raw_reads:
        raw_reads = item
        predict_adapter(raw_reads, ribotype, tmp_file_path)
        make_mrna_rrna_index(transcript_fasta, ribosome_fasta, thread, tmp_file_path, raw_reads, soft_log_folder)
        deal_raw_reads(raw_reads, tmp_file_path, thread, transcript_fasta, ribosome_fasta, trimmomatic, ribotype, soft_log_folder)

    print(get_time(), 'deal_raw_reads finished!!!')

if __name__ == '__main__':
    main()

