#ifndef _OCOCO_IMPL_H_
#define _OCOCO_IMPL_H_

#pragma once


/*************************
 *** Consensus calling ***
 *************************/

char ococo::cons_call_stoch(const pos_stats_uncompr_t &psu){
    if(psu.sum==0){
        return 'N';
    }
    
    const int32_t prefsum[]={
        psu.counters[0],
        psu.counters[0]+psu.counters[1],
        psu.counters[0]+psu.counters[1]+psu.counters[2],
        psu.counters[0]+psu.counters[1]+psu.counters[2]+psu.counters[3]
    };
    
    assert(prefsum[3]==psu.sum);
    
    const int32_t rn=rand() % psu.sum;
    for(int32_t i=0;i<4;i++){
        if (rn < prefsum[i]) {
            return nt4_nt256[i];
        }
    }
    
    return 'n';
}

char ococo::cons_call_stoch_amb(const pos_stats_uncompr_t &psu){
    if(psu.sum==0){
        return 'N';
    }
    
    nt16_t nucl_nt16=0;
    
    for(int32_t i=0;i<4;i++){
        const int32_t rn=rand() % psu.sum;
        
        if(rn < psu.counters[i]){
            nucl_nt16|=nt4_nt16[i];
        }
        
    }
    return 0;
}

char ococo::cons_call_maj(const pos_stats_uncompr_t &psu){
    char nucl_nt256='N';
    int32_t max=0;
    for(int32_t i=0;i<4;i++){
        if (psu.counters[i]>max){
            max=psu.counters[i];
            nucl_nt256=nt4_nt256[i];
        }
    }
    return nucl_nt256;
}

char ococo::cons_call_maj_amb(const pos_stats_uncompr_t &psu){
    nt16_t nucl_nt16=0;
    int32_t max=0;
    for(int32_t i=0;i<4;i++){
        if (psu.counters[i]>max){
            max=psu.counters[i];
            nucl_nt16=nt4_nt16[i];
        }
        else if(psu.counters[i]>max){
            nucl_nt16|=nt4_nt16[i];
        }
    }
    return nt16_nt256[nucl_nt16];
}



/****************************
 *** Consensus parameters ***
 ****************************/

ococo::consensus_params_t::consensus_params_t():
mode(BATCH),
strategy(MAJORITY),
min_mapq(1),
min_baseq(0),
init_ref_weight(2)
{
    cons_alg[strategy_t::STOCHASTIC]=&cons_call_stoch;
    cons_alg[strategy_t::STOCHASTIC_AMB]=&cons_call_stoch_amb;
    cons_alg[strategy_t::MAJORITY]=&cons_call_maj;
    cons_alg[strategy_t::MAJORITY_AMB]=&cons_call_maj_amb;
}


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
    if(seq_stats!=nullptr){
        for (int32_t seqid=0;seqid<n_seqs;seqid++){
            delete[] seq_stats[seqid];
        }
    }
    delete[] seq_active;
    delete[] seq_len;
    delete[] seq_name;
    delete[] seq_comment;
    delete[] seq_stats;
}


template<typename T, int counter_size, int refbase_size>
int ococo::stats_t<T,counter_size,refbase_size>::load_fasta(const std::string &fasta_fn) {
    gzFile fp;
    kseq_t *seq;
    int l;
    fp = gzopen(fasta_fn.c_str(), "r");
    seq = kseq_init(fp);
    
    constexpr int32_t max_counter_value = ococo::right_full_mask<T, counter_size>();
    
    for(int seqid=0;(l = kseq_read(seq)) >= 0;seqid++) {
        
        if(seq_name[seqid].compare(seq->name.s)!=0){
            error("Sequence names in BAM/SAM and in FASTA do not correspond ('%s'!='%s').\n",seq_name[seqid].c_str(),seq->name.s);
            return -1;
        }
        
        if(seq_len[seqid] != static_cast<int64_t>(seq->seq.l)){
            error("Sequence lengths in BAM/SAM and in FASTA do not correspond (%" PRId64 "!=%" PRId64 ").\n",static_cast<int64_t>(seq->seq.l), static_cast<int64_t>(seq_len[seqid]));
            return -1;
        }
        
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
int ococo::stats_t<T,counter_size,refbase_size>::save_fasta(const std::string &fasta_fn) const {
    assert(check_allocation());

    FILE *fasta_file=nullptr;
    fasta_file=fopen(fasta_fn.c_str(),"w+");
    
    char fasta_buffer[fasta_line_l];
    for(int s=0;s<n_seqs;s++){
        //printf("%s\n",seq_name[s]);
        if(!seq_comment[s].empty()){
            fprintf(fasta_file,">%s %s\n",seq_name[s].c_str(),seq_comment[s].c_str());
        }
        else{
            fprintf(fasta_file,">%s\n",seq_name[s].c_str());
        }
        
        for (int64_t i=0,j=0;i<seq_len[s];i++,j++){
            get_nucl_nt256(s,i,fasta_buffer[j]);
            
            if(j==fasta_line_l-1 || i==seq_len[s]-1){
                fwrite(fasta_buffer,1,j+1,fasta_file);
                fwrite("\n",1,1,fasta_file);
                j=-1;
            }
        }
    }
    
    
    fclose(fasta_file);
    
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
    
    if(n_seqs_loaded != n_seqs){
        error("Numbers of sequences in stats and SAM/BAM do not correspond %" PRId32 "!=%" PRId32 ").\n",n_seqs_loaded, n_seqs);
        return -1;
    }
    
    for(int seqid=0;seqid<n_seqs;seqid++){
        /* sequence */
        
        single_seq_serial_t seq_ser;
        fread(&seq_ser,sizeof(single_seq_serial_t),1,fo);
        
        if(seq_ser.seq_active != seq_active[seqid]){
            error("Active sequences in stats and SAM/BAM do not correspond (seqid %" PRId32 ").\n", seqid);
            return -1;
        }
        
        if(seq_ser.seq_len != seq_len[seqid]){
            error("Sequence lengths in stats and SAM/BAM do not correspond (seqid %" PRId32 ", %" PRId64 "!=%" PRId64 ").\n",seqid,seq_ser.seq_len,seq_len[seqid]);
            return -1;
        }
        
        if(seq_name[seqid].compare(seq_ser.seq_name)!=0){
            error("Sequence names in stats and SAM/BAM do not correspond (seqid %" PRId32 ", '%s'!='%s').\n", seqid, seq_ser.seq_name,seq_name[seqid].c_str());
            return -1;
        }
        
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
int ococo::stats_t<T,counter_size,refbase_size>::call_consensus(FILE *vcf_file) {
    assert(check_allocation());
    
    for(int32_t seqid=0;seqid<n_seqs;seqid++){
        for (int64_t pos=0;pos<seq_len[seqid];pos++){
            call_consensus_position(vcf_file,seqid, pos);
        }
    }
    
    return 0;
}

template<typename T, int counter_size, int refbase_size>
int ococo::stats_t<T,counter_size,refbase_size>::call_consensus_position(FILE *vcf_file, int32_t seqid, int64_t pos) {
    pos_stats_uncompr_t psu;
    decompress_position_stats(seq_stats[seqid][pos], psu);
    
    char old_base_nt256;
    get_nucl_nt256(seqid,pos,old_base_nt256);
    //const char new_base_nt256=cons_call_maj(psu);
    const char new_base_nt256=(params.cons_alg[params.strategy])(psu);
    
    if(old_base_nt256!=new_base_nt256){
        if(vcf_file!=nullptr){
            print_vcf_substitution(vcf_file,seqid,pos,old_base_nt256,new_base_nt256,psu);
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
int ococo::stats_t<T,counter_size,refbase_size>::print_vcf_header(FILE *vcf_file, std::string cmd,std::string fasta) const {
    assert(check_allocation());
    assert(vcf_file!=nullptr);
    
    std::time_t tt = std::time(nullptr);
    tm *tm = localtime(&tt);
    
    
    fprintf(vcf_file,
            "##fileformat=VCFv4.3\n"
            "##fileDate=%04d%02d%02d\n"
            "##source=Ococo\n",
            tm->tm_year + 1900,
            tm->tm_mon + 1,
            tm->tm_mday
            );
    
    if(!cmd.empty()){
        fprintf(vcf_file, "##ococo_command=%s\n", cmd.c_str());
    }
    fprintf(vcf_file,"##ococo_stats_datatype_size=%zubits\n",8*sizeof(T));
    fprintf(vcf_file,"##ococo_counter_size=%dbits\n", counter_size);
    
    if(!fasta.empty()){
        fprintf(vcf_file, "##reference=%s\n", fasta.c_str());
    }
    
    for (int seqid=0;seqid<n_seqs;seqid++){
        fprintf(vcf_file,"##contig=<ID=%s,length=%" PRId64 ">\n",seq_name[seqid].c_str(),seq_len[seqid]);
    }
    
    fprintf(vcf_file,"##INFO=<ID=CS,Number=4,Type=Integer,Description=\"Values of A,C,G,T counters.\">\n");
    fprintf(vcf_file,"##INFO=<ID=SUM,Number=1,Type=Integer,Description=\"Sum of A,C,G,T counters.\">\n");
    fprintf(vcf_file,"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
    
    return 0;
}

template<typename T, int counter_size, int refbase_size>
int ococo::stats_t<T,counter_size,refbase_size>::print_vcf_substitution(FILE *vcf_file, int32_t seqid, int64_t pos, char old_base, char new_base, const pos_stats_uncompr_t &psu) const {
    assert(check_allocation());
    assert(vcf_file!=nullptr);
    
    fprintf(vcf_file,"%s\t%" PRId64 "\t.\t%c\t%c\t100\tPASS\tCS=%" PRId32 ",%" PRId32 ",%" PRId32 ",%" PRId32 ";SUM=%" PRId32 "\n",
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