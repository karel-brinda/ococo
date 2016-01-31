#ifndef _OCOCO_IMPL_H_
#define _OCOCO_IMPL_H_

#pragma once

/****************************
 *** Consensus parameters ***
 ****************************/

ococo::consensus_params_t::consensus_params_t():
mode(BATCH),
strategy(MAJORITY),
min_mapq(1),
min_baseq(0),
init_ref_weight(2),
vcf_fo(nullptr),
fasta_cons_fo(nullptr)
{}


/***********************
 *** Main statistics ***
 ***********************/

template<typename T, int counter_size, int refbase_size>
ococo::stats_t<T,counter_size,refbase_size>::stats_t(ococo::consensus_params_t parameters,bam_hdr_t &h):
n_seqs(h.n_targets),
seq_active(new (std::nothrow) bool[n_seqs]()),
seq_len(new (std::nothrow) int64_t[n_seqs]()),
seq_name(new (std::nothrow) std::string[n_seqs]()),
seq_comment(new (std::nothrow) std::string[n_seqs]()),
seq_stats(new (std::nothrow) T*[n_seqs]()),
params(parameters)
{
    for (int seqid=0;seqid<n_seqs;seqid++){
        seq_len[seqid]=h.target_len[seqid];
        seq_active[seqid]=true;
        seq_name[seqid]=std::string(h.target_name[seqid]);

        seq_stats[seqid]=new (std::nothrow) T[seq_len[seqid]]();
    }
    
}


template<typename T, int counter_size, int refbase_size>
ococo::stats_t<T,counter_size,refbase_size>::~stats_t(){
    for (int32_t seqid=0;seqid<n_seqs;seqid++){
        delete[] seq_stats[seqid];
    }
    delete[] seq_active;
    delete[] seq_len;
    delete[] seq_name;
    delete[] seq_comment;
    delete[] seq_stats;
}


template<typename T, int counter_size, int refbase_size>
int ococo::stats_t<T,counter_size,refbase_size>::
load_fasta(const std::string &fasta_fn) {
    gzFile fp;
    kseq_t *seq;
    int l;
    fp = gzopen(fasta_fn.c_str(), "r");
    seq = kseq_init(fp);
    
    constexpr int32_t max_counter_value = ococo::right_full_mask<T, counter_size>();
    
    for(int seqid=0;(l = kseq_read(seq)) >= 0;seqid++) {
        
        assert(seq_name[seqid].compare(seq->name.s)==0);
        assert(static_cast<int64_t>(seq->seq.l) == seq_len[seqid]);
        
        if (seq->comment.l && seq_comment[seqid].empty()){
            seq_comment[seqid]=std::string(seq->comment.s);
        }
        
        for(int64_t pos=0; pos < static_cast<int64_t>(seq->seq.l); pos++){
            assert(seq_stats[seqid][pos]==0);
            
            pos_stats_uncompr_t psu = {0,{0,0,0,0},0};
            psu.nt16=nt256_nt16[static_cast<int32_t>(seq->seq.s[pos])];
            
            if (psu.nt16!=nt256_nt16[static_cast<int32_t>('N')]){
                for(int32_t i=0;i<4;i++){
                    psu.counters[i]= ( (0x1<<i) & psu.nt16) ? std::min(params.init_ref_weight,max_counter_value) : 0;
                }
            }
            
            psu.sum=psu.counters[0]+psu.counters[1]+psu.counters[2]+psu.counters[3];
            seq_stats[seqid][pos]=compress_position_stats(psu);
        }
    }
    kseq_destroy(seq); // STEP 5: destroy seq
    gzclose(fp); // STEP 6: close the file handler
    return 0;
}


template<typename T, int counter_size, int refbase_size>
int ococo::stats_t<T,counter_size,refbase_size>::save_fasta() const {
    assert(check_allocation());
    assert(params.fasta_cons_fo!=nullptr);
    
    char fasta_buffer[fasta_line_l];
    for(int s=0;s<n_seqs;s++){
        //printf("%s\n",seq_name[s]);
        if(!seq_comment[s].empty()){
            fprintf(params.fasta_cons_fo,">%s %s\n",seq_name[s].c_str(),seq_comment[s].c_str());
        }
        else{
            fprintf(params.fasta_cons_fo,">%s\n",seq_name[s].c_str());
        }
        
        for (int64_t i=0,j=0;i<seq_len[s];i++,j++){
            get_nucl_nt256(s,i,fasta_buffer[j]);
            
            if(j==fasta_line_l-1 || i==seq_len[s]-1){
                fwrite(fasta_buffer,1,j+1,params.fasta_cons_fo);
                fwrite("\n",1,1,params.fasta_cons_fo);
                j=-1;
            }
        }
    }
    
    return 0;
}


template<typename T, int counter_size, int refbase_size>
bool ococo::stats_t<T,counter_size,refbase_size>::check_allocation() const {

    if(seq_active==nullptr || seq_len==nullptr || seq_stats==nullptr || seq_name==nullptr || seq_comment==nullptr){
        return false;
    }
    
    for(int i=0;i<n_seqs;i++){
        if(seq_stats[i]==nullptr){
            return false;
        }
    }
    
    return true;
}


template<typename T, int counter_size, int refbase_size>
bool ococo::stats_t<T,counter_size,refbase_size>::check_headers_bam_hdr(const bam_hdr_t &h) const {
    if (!check_allocation()) return false;
    
    for(int32_t seqid=0;seqid<n_seqs;seqid++){
        if(seq_len[seqid]!=static_cast<int64_t>(h.target_len[seqid])) {
            return false;
        }
        if(seq_name[seqid].compare(h.target_name[seqid])!=0) {
            return false;
        }
    }
    
    return true;
}


template<typename T, int counter_size, int refbase_size>
int ococo::stats_t<T,counter_size,refbase_size>::import_stats(const std::string &stats_fn){
    assert(check_allocation());
    
    FILE *fo=fopen(stats_fn.c_str(),"r");
    
    /* number of seqs */
    int32_t n_seqs_loaded;
    fread(&n_seqs_loaded,sizeof(int32_t),1,fo);
    assert(n_seqs_loaded = n_seqs);
    
    for(int seqid=0;seqid<n_seqs;seqid++){
        /* sequence */
        
        single_seq_serial_t seq_ser;
        fread(&seq_ser,sizeof(single_seq_serial_t),1,fo);
        assert(seq_ser.seq_active == seq_active[seqid]);
        assert(seq_ser.seq_len == seq_len[seqid]);
        assert(seq_name[seqid].compare(seq_ser.seq_name)==0);
        
        fread(seq_stats[seqid],sizeof(T),seq_len[seqid],fo);
    }
    
    fclose(fo);
    return 0;
}

template<typename T, int counter_size, int refbase_size>
int ococo::stats_t<T,counter_size,refbase_size>::export_stats(const std::string &stats_fn) const {
    assert(check_allocation());
    
    FILE *fo=fopen(stats_fn.c_str(),"w+");
    
    /* number of seqs */
    fwrite(&n_seqs,sizeof(int32_t),1,fo);
    
    for(int seqid=0;seqid<n_seqs;seqid++){
        /* sequence */
        single_seq_serial_t seq_ser={0};
        seq_ser.seq_active=seq_active[seqid];
        seq_ser.seq_len=seq_len[seqid];
        strncpy(seq_ser.seq_name,seq_name[seqid].c_str(),999);
        strncpy(seq_ser.seq_name,seq_name[seqid].c_str(),999);
        fwrite(&seq_ser,sizeof(single_seq_serial_t),1,fo);
        fwrite(seq_stats[seqid],sizeof(T),seq_len[seqid],fo);
    }
    
    fclose(fo);
    return 0;
}


template<typename T, int counter_size, int refbase_size>
int ococo::stats_t<T,counter_size,refbase_size>::call_consensus() {
    assert(check_allocation());
    
    for(int32_t seqid=0;seqid<n_seqs;seqid++){
        for (int64_t pos=0;pos<seq_len[seqid];pos++){
            call_consensus_position(seqid, pos);
        }
    }
    
    return 0;
}

template<typename T, int counter_size, int refbase_size>
int ococo::stats_t<T,counter_size,refbase_size>::call_consensus_position(int32_t seqid, int64_t pos) {
    pos_stats_uncompr_t psu;
    decompress_position_stats(seq_stats[seqid][pos], psu);
    
    char old_base_nt256;
    get_nucl_nt256(seqid,pos,old_base_nt256);
    const char new_base_nt256=rand_nucl(psu);
    
    if(old_base_nt256!=new_base_nt256){
        if(params.vcf_fo!=nullptr){
            print_vcf_substitution(seqid,pos,old_base_nt256,new_base_nt256,psu);
        }
        
        set_nucl_nt256(seqid,pos,new_base_nt256);
    }
    
    return 0;
}


template<typename T, int counter_size, int refbase_size>
T ococo::stats_t<T,counter_size,refbase_size>::compress_position_stats(const pos_stats_uncompr_t &psu) {
    
    T psc=0;
    
    for(int32_t i=0;i<4;i++){
        psc <<= counter_size;
        psc |= psu.counters[i] & right_full_mask<T,counter_size>();
    }
    
    psc <<= refbase_size;
    psc |= psu.nt16;
    
    return psc;
}


template<typename T, int counter_size, int refbase_size>
void ococo::stats_t<T,counter_size,refbase_size>::decompress_position_stats(T psc, pos_stats_uncompr_t &psu) {
    
    psu.nt16= psc & right_full_mask<T,refbase_size>();
    psc >>= refbase_size;
    
    psu.sum=0;
    for(int32_t i=3;i>=0;i--){
        psu.counters[i] = psc & right_full_mask<T,counter_size>();
        psu.sum += psu.counters[i];
        psc >>= counter_size;
    }
}

template<typename T, int counter_size, int refbase_size>
int ococo::stats_t<T,counter_size,refbase_size>::print_vcf_header(std::string cmd,std::string fasta) const {
    assert(check_allocation());
    assert(params.vcf_fo!=nullptr);
    
    std::time_t tt = std::time(nullptr);
    tm *tm = localtime(&tt);
    
    
    fprintf(params.vcf_fo,
            "##fileformat=VCFv4.3\n"
            "##fileDate=%04d%02d%02d\n"
            "##source=Ococo\n",
            tm->tm_year + 1900,
            tm->tm_mon + 1,
            tm->tm_mday
            );
    
    if(!cmd.empty()){
        fprintf(params.vcf_fo, "##ococo_command=%s\n", cmd.c_str());
    }
    fprintf(params.vcf_fo,"##ococo_stats_datatype_size=%zubits\n",8*sizeof(T));
    fprintf(params.vcf_fo,"##ococo_counter_size=%dbits\n", counter_size);
    
    if(!fasta.empty()){
        fprintf(params.vcf_fo, "##reference=%s\n", fasta.c_str());
    }
    
    for (int seqid=0;seqid<n_seqs;seqid++){
        fprintf(params.vcf_fo,"##contig=<ID=%s,length=%" PRId64 ">\n",seq_name[seqid].c_str(),seq_len[seqid]);
    }
    
    fprintf(params.vcf_fo,"##INFO=<ID=CS,Number=4,Type=Integer,Description=\"Values of A,C,G,T counters.\">\n");
    fprintf(params.vcf_fo,"##INFO=<ID=SUM,Number=1,Type=Integer,Description=\"Sum of A,C,G,T counters.\">\n");
    fprintf(params.vcf_fo,"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
    
    return 0;
}

template<typename T, int counter_size, int refbase_size>
int ococo::stats_t<T,counter_size,refbase_size>::print_vcf_substitution(int32_t seqid, int64_t pos, char old_base, char new_base, const pos_stats_uncompr_t &psu) const {
    assert(check_allocation());
    assert(params.vcf_fo!=nullptr);
    
    fprintf(params.vcf_fo,"%s\t%" PRId64 "\t.\t%c\t%c\t100\tPASS\tCS=%" PRId32 ",%" PRId32 ",%" PRId32 ",%" PRId32 ";SUM=%" PRId32 "\n",
            seq_name[seqid].c_str(),
            pos+1,
            old_base,
            new_base,
            psu.counters[0],
            psu.counters[1],
            psu.counters[2],
            psu.counters[3],
            psu.sum
            );
    
    return 0;
}

template<typename T, int counter_size, int refbase_size>
std::string ococo::stats_t<T,counter_size,refbase_size>::debug_str_counters(int32_t seqid, int64_t pos) const {
    pos_stats_uncompr_t psu;
    decompress_position_stats(seq_stats[seqid][pos], psu);
    std::stringstream ss;
    ss << "(" << psu.counters[0] << "," << psu.counters[1] << "," << psu.counters[2] << "," << psu.counters[3] << ")";
    return ss.str();
}

template<typename T, int counter_size, int refbase_size>
void ococo::stats_t<T,counter_size,refbase_size>::debug_print_counters() const {
    for(int seqid=0;seqid<n_seqs;seqid++){
        fprintf(stderr,"%s\n",seq_name[seqid]);
        for (int64_t pos=0;pos<seq_len[seqid];pos++){
            fprintf(stderr,"%8" PRId64 " %04x \n",pos,seq_stats[seqid][pos]);
        }
    }
}

template<typename T, int counter_size, int refbase_size>
inline int ococo::stats_t<T,counter_size,refbase_size>::set_nucl_nt256(int32_t seqid, int64_t pos, const char &nt256) {
    nt16_t nt16 =nt256_nt16[static_cast<int>(nt256)];
    T n_psc=seq_stats[seqid][pos];
    n_psc>>=refbase_size;
    n_psc<<=refbase_size;
    n_psc |= nt16 & right_full_mask<T,refbase_size>();
    seq_stats[seqid][pos]=n_psc;
    return 0;
}

template<typename T, int counter_size, int refbase_size>
inline int ococo::stats_t<T,counter_size,refbase_size>::get_nucl_nt256(int32_t seqid, int64_t pos, char &nt256) const {
    nt256=nt16_nt256[
                     seq_stats[seqid][pos] & right_full_mask<T,refbase_size>()
                     ];
    return 0;
}

template<typename T, int counter_size, int refbase_size>
T ococo::stats_t<T,counter_size,refbase_size>::increment(T psc, nt4_t nt4){
    assert(0 <= nt4 && nt4 < 4);
    
    pos_stats_uncompr_t psu;
    decompress_position_stats(psc, psu);
    
    if(psu.counters[nt4]==right_full_mask<uint16_t, counter_size>()){
        psu.counters[0]>>=1;
        psu.counters[1]>>=1;
        psu.counters[2]>>=1;
        psu.counters[3]>>=1;
    }
    
    psu.counters[nt4]++;
    
    return compress_position_stats(psu);
}

#endif