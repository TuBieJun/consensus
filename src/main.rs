extern crate rust_htslib;
extern crate clap;




use std::env;
use std::collections::HashSet;
use std::collections::HashMap;
use rust_htslib::bam;
use rust_htslib::prelude::*;
use rust_htslib::bam::Read;
use std::str;
use rust_htslib::bam::Record;
use std::io;
use std::io::prelude::*;
use std::io::Write;
use std::fs::File;
use std::path::{Path, PathBuf};
use clap::{Arg, App, SubCommand};

fn reverse_complement(v:&Vec<u8>) -> Vec<u8> {
    let mut new_v:Vec<u8> = Vec::new();
    for c in v.iter().rev() {
       let c_r = match c {
            &65u8 => 84u8, // A -> T
            &84u8 => 65u8, // T -> A
            &67u8 => 71u8, // C -> G
            &71u8 => 67u8, // G -> C
            _  => 78u8,

       };
       new_v.push(c_r);
    }
    new_v
}

fn reverse_only(v:&[u8]) -> Vec<u8> {
    let mut new_v:Vec<u8> = Vec::new();
    for c in v.iter().rev() {
        new_v.push(*c);
    }
    new_v
}

fn check_flag(flag:&u16) -> bool {
    // check the record of flag 
    let flag_l:Vec<u16> = vec![181, 117, 137, 133, 
                               73, 89, 69, 153,
                               97, 81, 161, 145, 
                               129, 65, 177, 113,
                               99, 83, 163, 147,
                               0, 16];
    let flag_set:HashSet<&u16> = flag_l.iter().collect();
    if flag_set.contains(flag) {
        true
    } else {
        false
    }
}

fn bamLineParse(record:&Record) -> (String, String, String, i32, i32){

    let qname_str = String::from_utf8_lossy(record.qname()).to_string();
    let temp_v:Vec<&str> = qname_str.split('#').collect();
    let index = temp_v[1].to_string();
    let insert_size = record.insert_size();
    let mut pair;
    let mut strand;
    let mut aln_pos;
    if record.is_first_in_template() {
        pair = "R1";
        if record.is_reverse() {
            strand = "-+";
            aln_pos = record.mpos(); 
        } else {
            strand = "+-";
            aln_pos = record.pos(); 
        }
    } else {
        pair = "R2";
        if record.is_reverse() {
            strand = "+-";
            aln_pos = record.mpos();
        } else {
            strand = "-+";
            aln_pos = record.pos();
        }
    }

    (index, pair.to_string(), strand.to_string(), insert_size.abs(), aln_pos)
}

fn dup2consensus_work(group_k:&str, pair_dup:& Vec<(Vec<u8>, Vec<u8>, Vec<u8>)>, pair_flag: &str, o_buff: & mut File, min_t_s: &u32, min_t_p: &f32)  {
    let consensus_base:Vec<u8> = Vec::new();
    let index_identity = -1;
    let mut base_position:HashMap<usize, HashMap<&u8, u32>> = HashMap::new();
    let mut num_position:HashMap<usize, u32> = HashMap::new();
    let mut qual_position:HashMap<usize, u32> = HashMap::new();
    let mut consensus_source_reads = String::new();
    let mut seq_consensus:Vec<u8>; 
    let mut qual_consnesus:Vec<u8>;
    let mut baseEach_num_record = String::new();
    let mut baseEach_percent_record = String::new();

    
    if pair_dup.len() == 1 {
        seq_consensus = Vec::new();
        seq_consensus.extend(pair_dup[0].0.iter().cloned());
        let real_qual:Vec<u8> = pair_dup[0].1.iter().map(|x| x + 33).collect();
        qual_consnesus = real_qual;
        baseEach_num_record =  format!("1,{}.", vec![".,";seq_consensus.len()-2].join(""));
        baseEach_percent_record = format!("{}.", vec![".,";seq_consensus.len()-1].join(""));
    } else {
        for tup in pair_dup {
            for (pos, base)  in tup.0.iter().enumerate() {
                //println!("{}\t{}", pos, base);
                let n_p = num_position.entry(pos).or_insert(0);
                let q_p = qual_position.entry(pos).or_insert(0);
                *n_p += 1;
                *q_p += tup.1[pos] as u32;

                if !base_position.contains_key(&pos) {
                let mut base_count:HashMap<&u8, u32> = HashMap::new();
                *base_count.entry(base).or_insert(0) += 1;
                base_position.insert(pos, base_count);
                } else {
                *base_position.get_mut(&pos).unwrap().entry(base).or_insert(0) += 1;
                }
            }
        }
        seq_consensus = vec![78;base_position.len()];
        qual_consnesus =  vec![78;base_position.len()];
        let mut top_base_percent = 0.0f32;
        for i in 0..base_position.len() {
            let p_n = *num_position.get(&i).unwrap();
            if num_position.get(&i).unwrap() >= min_t_s {
                let top_base = base_position.get(&i).unwrap().iter()
                                        .max_by(|a, b| a.1.cmp(b.1))
                                        .unwrap().0;
                let top_base_n = *base_position.get(&i).unwrap().get(top_base).unwrap() as f32;
                let p_n_f = p_n as f32;
                top_base_percent = top_base_n / p_n_f;
                if top_base_percent >= *min_t_p {
                    seq_consensus[i] = **top_base;
                }
            }
            let sum_qual = *qual_position.get(&i).unwrap();
            qual_consnesus[i] = (sum_qual / p_n + 33) as u8 ;

            if (num_position[&(i as usize)] == num_position[&(0 as usize)] && baseEach_num_record.len() > 0) {
                baseEach_num_record.push_str(",.");
            } else {
                baseEach_num_record.push_str(&format!("{}", num_position[&(0 as usize)]));
            }

            let mut str_percent:String;
            if top_base_percent < 1.0f32 {
                str_percent = format!("{}", top_base_percent);
            } else {
                str_percent = String::from(".");
            }

            if (baseEach_percent_record.len() > 0) {
                baseEach_percent_record.push_str(&format!(",{}", str_percent));
            } else {
                baseEach_percent_record.push_str(&str_percent);
            
            }
            
        }
    }
    
    let temp_v:Vec<&str> = group_k.split('_').collect();
    let chrom = temp_v[1];
    let aln_pos = temp_v[2];
    let index = temp_v[0];
    let tlen_flag = temp_v[3];
    let fq_id = format!("@{}:{}|{}|{}#{} {}:N:0:{}:{}", chrom, aln_pos, tlen_flag, index, index_identity, pair_flag, baseEach_num_record, baseEach_percent_record);
    o_buff.write(&fq_id.as_bytes()).unwrap();
    for member_record in pair_dup {
        o_buff.write(b"\t").unwrap();
        o_buff.write(&member_record.2);
        o_buff.write(b"|").unwrap();
        o_buff.write(&member_record.0);
        o_buff.write(b"|").unwrap();
        let real_qual_temp:Vec<u8> = member_record.1.iter().map(|x| x + 33).collect();
        o_buff.write(&real_qual_temp);
    }
    o_buff.write(b"\n");

    o_buff.write(&seq_consensus).unwrap();
    o_buff.write(b"\n");
    o_buff.write(b"+\n").unwrap();
    o_buff.write(&qual_consnesus).unwrap(); 
    o_buff.write(b"\n");
}

fn  block_consensus(pd:& mut HashMap<String, HashMap<String, Vec<(Vec<u8>, Vec<u8>, Vec<u8>)>>>, chrom: i32, use_check:bool, min_t_s: &u32, min_t_p: &f32, o_r1: & mut File, o_r2: & mut File) {
    let mut can_consensus_v:Vec<String> = Vec::new();
    
    for (key_dupGroup, pair_dup) in pd.iter() {
        if pair_dup.len() == 2 {
            if (use_check && pair_dup.get("R1").unwrap().len() == pair_dup.get("R2").unwrap().len()) || (!use_check){
                can_consensus_v.push(key_dupGroup.to_string());
            }
        }
    }

    for k in can_consensus_v {
        dup2consensus_work(&k, pd.get(&k).unwrap().get("R1").unwrap(),"1", o_r1, min_t_s, min_t_p);
        dup2consensus_work(&k, pd.get(&k).unwrap().get("R2").unwrap(),"2", o_r2, min_t_s, min_t_p);
        pd.remove(&k);
    }     
}

fn main() {

    // Parse the argument line
    let matches = App::new("ctdna consensus program")
                          .version("1.0")
                          .author("Teng Li. <707704459@qq.com>")
                          .about("conensus the bam file  to fastq by index")
                          .arg(Arg::with_name("bam").short("b").long("bam").value_name("BAM")
                               .required(true)
                               .help("the input bam path "))
                          .arg(Arg::with_name("outdir").short("o").long("outdir")
                               .value_name("OUTDIR")
                               .required(true)
                               .help("the output dir path"))
                          .arg(Arg::with_name("prefix").long("prefix")
                               .value_name("PREFIX")
                               .required(true)
                               .help("the output fastq file prefix"))
                          .arg(Arg::with_name("min_group").short("n")
                               .value_name("MIN_GROUP")
                               .default_value("2")
                               .help("the number of min top base support"))
                          .arg(Arg::with_name("min_percent").short("p")
                               .value_name("MIN_PERCENT")
                               .default_value("0.8")
                               .help("the number of min top base percent"))
                          .get_matches();

    let bam_path = matches.value_of("bam").unwrap();
    let outdir = matches.value_of("outdir").unwrap();
    let prefix = matches.value_of("prefix").unwrap();
    let min_t_s:u32 = matches.value_of("min_group").unwrap().trim().parse()
                     .expect("the -n must be a int number");
    let min_t_p:f32 = matches.value_of("min_percent").unwrap().trim().parse()
                     .expect("the -p must be a float number");

    let f_r1_path = Path::new(outdir).join(format!("{}_consensus_R1_fastq", prefix));
    let f_r2_path = Path::new(outdir).join(format!("{}_consensus_R2_fastq", prefix));

    
    let mut f_r1 = File::create(f_r1_path).unwrap();
    let mut f_r2 = File::create(f_r2_path).unwrap();


    let mut bam = bam::Reader::from_path(&bam_path).unwrap();
    let mut iter_pos = 0;
    let mut chrom = 0;
    
    let mut position_dupGroup: HashMap<String, HashMap<String, Vec<(Vec<u8>, Vec<u8>, Vec<u8>)>>> = HashMap::new();
    for r in bam.records() {
        let record = r.unwrap();
        if record.pos() != iter_pos {
            block_consensus(&mut position_dupGroup, chrom, true, & min_t_s, & min_t_p,  & mut f_r1, & mut f_r2);
        }
        iter_pos = record.pos();
        if check_flag(&record.flags()) {
            
            let (index, pair, strand, insert_size, aln_pos) = bamLineParse(&record);
            let chrom = record.tid();
            let raw_seq_b = record.seq().as_bytes();
            let raw_qual_b = record.qual();
            let mut new_seq_b:Vec<u8> = Vec::new();
            let mut new_qual_b:Vec<u8> = Vec::new();
            let mut qname:Vec<u8> = Vec::new();

            if record.is_reverse() {
                new_seq_b = reverse_complement(&raw_seq_b);
                new_qual_b = reverse_only(raw_qual_b);
            } else {
                new_seq_b = raw_seq_b;
                new_qual_b.extend_from_slice(raw_qual_b);
            }
            qname.extend_from_slice(record.qname());
            let key_dupGroup = format!("{}_{}_{}_{}{}", index, chrom, aln_pos, strand, insert_size);
            if position_dupGroup.contains_key(&key_dupGroup) {
                    position_dupGroup.get_mut(&key_dupGroup).unwrap().entry(pair).or_insert(Vec::new()).push((new_seq_b, new_qual_b, qname));

            } else {
                let mut pair_dup = HashMap::new();
                pair_dup.insert(pair,vec![(new_seq_b, new_qual_b, qname)]);
                position_dupGroup.insert(key_dupGroup, pair_dup);
            }
        }
        
    }
    block_consensus(& mut position_dupGroup, chrom, false, &min_t_s, &min_t_p, & mut f_r1, & mut f_r2);
}
