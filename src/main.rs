use rust_htslib::bcf::{Reader as bcf_Reader, Read as bcf_Read, Writer as bcf_Writer};
use rust_htslib::tbx::{Reader as tbx_Reader, Read as tbx_Read};
use csv::{Reader as csv_Reader};
use std::collections::HashSet;
use clap::Parser;
use std::path::PathBuf;
use std::error::Error;
use std::fs::File;

#[derive(Parser,Debug)]
struct Cli {
    /// Path to indexed VCF/BCF
    #[arg(short = 'i', long = "invcf")]
    bcfpath: PathBuf,
    /// Path to a bed file with cCREs coordinates obtained from screenDB.
    #[arg(short = 'c', long = "ccres")]
    ccre_bed: PathBuf,
    /// Path to gene links to ccre files formatted as bed.
    #[arg(short = 'g', long = "gene_link",value_delimiter = ',', num_args = 1..)]
    ccre_to_gene: Vec<PathBuf>,
    /// Path to biodata bed files per condition (cond1:path1,cond2,path2,...) 
    #[arg(short = 'b', long = "biodata",value_delimiter = ',', num_args = 1.., value_parser = parse_key_val::<String, PathBuf,f32>)]
    biodata: Option<Vec<(String,PathBuf,f32)>>,
    /// Output vcf
    #[arg(short = 'b', long = "biodata",default_value= "screen_annotated.vcf.gz")]
    outfile: PathBuf,
    /// gene annotation files in tsv format (First column should hold the gene name information)
    #[arg(short = 'g', long = "gene_link",value_delimiter = ',', num_args = 1..)]
    gene_annot_files: Option<Vec<PathBuf>>,
}

fn parse_key_val<T, U, O>(s: &str) -> Result<(T, U, O), Box<dyn Error + Send + Sync + 'static>>
    where
        T: std::str::FromStr,
        T::Err: Error + Send + Sync + 'static,
        U: std::str::FromStr,
        U::Err: Error + Send + Sync + 'static,
        O: std::str::FromStr,
        O::Err: Error + Send + Sync + 'static,
    {
        let args=s.split(':').map(|e| e.to_string()).collect::<Vec<String>>();
        let result = match args.len(){
            3 => Ok((args[0].parse()?, args[1].parse()?, args[2].parse()?)),
            2 => Ok((args[0].parse()?, args[1].parse()?, "2.0".to_string().parse()?)),
            _ => panic!("Cannot parse biodata argument")
        };
        result
    }


fn gene_annot(gene_name: &str, gene_annot: &mut csv_Reader<File>) -> Vec<String>{
    for rec in gene_annot.records(){
        let record_raw=rec.expect("Could not parse record");
        let record:Vec<&str>=record_raw.iter().collect();
        if record[0].trim().to_uppercase()==gene_name.trim().to_uppercase(){
            return record[1..record.len()].iter().map(|e| e.to_string()).collect::<Vec<String>>()
        }
    }
    vec![String::new()]
}

fn fetch_ccres(contig: &u64, start: &u64, end: &u64, ccres_bed:&mut tbx_Reader ) -> Vec<[String;3]> {
    // Retrieves ccreID of elements overlapping variant in ccres_bed if any
    ccres_bed.fetch(*contig,*start,*end).expect("Could not fecth region");
    let mut ccre_list:Vec<[String;3]>= Vec::new();
    for rec in ccres_bed.records(){
        let record= rec.expect("could not parse record");
        let ccre_name:String=std::str::from_utf8(&record)
            .expect("Fail to parse UTF8 string")
            .split("\t")
            .collect::<Vec<&str>>()[4]
            .to_string();
        let ccre_name2:String=std::str::from_utf8(&record)
            .expect("Fail to parse UTF8 string")
            .split("\t")
            .collect::<Vec<&str>>()[3]
            .to_string();
        let ccre_type_all:String=std::str::from_utf8(&record)
            .expect("Fail to parse UTF8 string")
            .split("\t")
            .collect::<Vec<&str>>()[5]
            .to_string();
        let ccre_type_list:Vec<&str>= ccre_type_all.split(",")
            .collect::<Vec<&str>>();
        for ccre_type in ccre_type_list{
            ccre_list.push([ccre_name.clone(),ccre_name2.clone(),ccre_type.to_string()]);
        }
        
    }
    ccre_list

}

fn fetch_genes(contig:u64, start:u64, end:u64, ccre_name: &str, ccre_to_gene_files: &mut Vec<tbx_Reader>, gene_annot_files: &mut Vec<csv_Reader<File>>)-> Vec<String>{
    // Get the genes linked to the ccre
    let mut gene_list:Vec<String>=Vec::new();
    for ccre_to_gene in ccre_to_gene_files {
        ccre_to_gene.fetch(contig,start,end).expect("cCould not fetch region.");
        
        for rec in ccre_to_gene.records(){
            let raw_record = rec.expect("Could not parse record");
            let record = std::str::from_utf8(&raw_record)
                .expect("Fail to parse UTF8 string")
                .split("\t")
                .collect::<Vec<&str>>();
            if &record[3]==&ccre_name{    
                let gene_rec=&record[5].trim().to_string();
                let mut annots:String=String::new();
                for gene_annot_file in &mut *gene_annot_files{
                    annots=annots+&gene_annot(gene_rec, gene_annot_file).join(":")+"|";
                }
                let gene_rec=gene_rec.to_owned()+"|"+&record[&record.len()-5]+"|"+&record[&record.len()-3].replace(",","_")+"|"+&record[&record.len()-1]+"|"+&annots;
                if !(gene_list.contains(&gene_rec)){
                    gene_list.push(gene_rec);
                }
            }
        } 
    }
    gene_list
}

fn fetch_zscores(contig:u64, start:u64, end:u64, ccre_id:&str, zscore_files:&mut Vec<(String,tbx_Reader,f32)>) -> Vec<String>{
    let mut cell_types:Vec<String>=Vec::new();
    for (condition,zscores,zscore_threshold) in zscore_files{
        zscores.fetch(contig,start,end).expect("Could not fetch region.");

        let raw_header_zscore=zscores.header()
            .clone()[0]
            .to_string();
        let header_zscore= raw_header_zscore.split("\t")
            .collect::<Vec<&str>>();
        for rec in zscores.records(){
            let raw_record = rec.expect("Could not parse record");
            let record = std::str::from_utf8(&raw_record)
                .expect("Fail to parse UTF8 string")
                .split("\t")
                .collect::<Vec<&str>>();

            if &record[3]==&ccre_id{
                for (i,val) in record.iter().enumerate().skip(4){
                    if val.parse::<f32>().expect("Could not parse value") > *zscore_threshold {
                        cell_types.push(condition.to_owned()+"|"+header_zscore[i])
                    }
                }
                break;  
            }
        }
    }
    cell_types
}


fn annot_var(variant_record: &mut rust_htslib::bcf::record::Record, mut ccres_bed:&mut tbx_Reader, mut ccre_to_gene: &mut Vec<tbx_Reader>, mut zscores: &mut Vec<(String,tbx_Reader,f32)>, gene_annot_files: &mut Vec<csv_Reader<File>>) {
    let contig:u64 = variant_record.rid().unwrap() as u64;
    let start:u64 = variant_record.pos() as u64;
    let end:u64 = variant_record.end() as u64;

    let ccre_details=fetch_ccres(&contig,&start,&end,&mut ccres_bed);
    let mut rec_type:Vec<Vec<u8>>= Vec::new();
    let mut gene_annot:Vec<Vec<u8>>= Vec::new();
    let mut biodata_annot:Vec<Vec<u8>>= Vec::new();
    for ccre_rec in &ccre_details{
        let ccre_annot:Vec<u8>=(ccre_rec[0].clone()+"|"+&ccre_rec[2]).to_string().into_bytes();
        if !(rec_type.contains(&ccre_annot)){
            rec_type.push(ccre_annot)
        }
    }
    let unique_ccre_names:&Vec<String>=&ccre_details.iter().map(|r| r[0].clone()).collect();
    let unique_ccre_names:HashSet<&String>=HashSet::from_iter(unique_ccre_names.iter());
    for ccre_name in &unique_ccre_names{
        let genes=fetch_genes(contig,start,end,ccre_name,&mut ccre_to_gene, gene_annot_files);
        for gene in genes {
            let ann_rec=(ccre_name.to_owned().clone()+"|"+&gene).into_bytes();
            gene_annot.push(ann_rec);
        }
    }
    let unique_ccre_ids:&Vec<String>=&ccre_details.iter().map(|r| r[1].clone()).collect();
    let unique_ccre_ids:HashSet<&String>=HashSet::from_iter(unique_ccre_ids.iter());

    for ccre_id in &unique_ccre_ids{
        let biodata_scores=fetch_zscores(contig,start,end,ccre_id,&mut zscores);
        for biodata in biodata_scores {
            let ann_rec=biodata.replace(",","_").into_bytes();
            if !(biodata_annot.contains(&ann_rec)){
                biodata_annot.push(ann_rec);
            }
        }
    } 
    let rec_type:Vec<&[u8]> = rec_type.iter().map(|a| a.as_slice()).collect::<Vec<&[u8]>>();
    let gene_annot:Vec<&[u8]> = gene_annot.iter().map(|a| a.as_slice()).collect::<Vec<&[u8]>>();
    let biodata_annot:Vec<&[u8]> = biodata_annot.iter().map(|a| a.as_slice()).collect::<Vec<&[u8]>>();
    variant_record.push_info_string(b"ScreenElements",&rec_type).expect("Failed to add annot to variant");
    variant_record.push_info_string(b"ScreenGenes",&gene_annot).expect("Failed to add annot to variant");
    variant_record.push_info_string(b"ScreenBiodata",&biodata_annot).expect("Failed to add annot to variant");
}


fn main() {
    let args=Cli::parse();
    println!("{:?}", args);

    let mut invcf=bcf_Reader::from_path(args.bcfpath).expect("File not found.");
    let mut ccres_bed=tbx_Reader::from_path(args.ccre_bed).expect("File not found");
    let mut biodata_files:Vec<(String,tbx_Reader,f32)>= match args.biodata {
        Some(v) => v.iter()
        .map(|(cond,path,thr)|(cond.clone(),tbx_Reader::from_path(path).expect("File not found"),*thr))
        .collect(),
        None => Vec::new()
    };
    let mut ccre_to_gene_files:Vec<tbx_Reader>=args.ccre_to_gene
        .iter()
        .map(|path| tbx_Reader::from_path(path).expect("File not found"))
        .collect();
    let mut gene_annot_files:Vec<csv_Reader<File>>= match args.gene_annot_files{
        Some(v) => v.iter().map(|f| csv::ReaderBuilder::new()
                .has_headers(true)
                .delimiter(b'\t')
                .from_path(f).expect("File not found"))
            .collect(),
        None => Vec::new()
    };
        

    let mut out_header=rust_htslib::bcf::Header::from_template(invcf.header());
    out_header.push_record(b"##INFO=<ID=ScreenElements, Type=String, Number=., Description=\"Regulatory elements from Screen\">");
    out_header.push_record(b"##INFO=<ID=ScreenGenes, Type=String, Number=., Description=\"Regulatory elements from Screen\">");
    out_header.push_record(b"##INFO=<ID=ScreenBiodata, Type=String, Number=., Description=\"Regulatory elements from Screen\">");
    let mut outvcf=bcf_Writer::from_path(args.outfile, &out_header, false,rust_htslib::bcf::Format::Vcf ).expect("Could not open file.");


    
    for rec in invcf.records(){
        let mut record=rec.expect("Could not parse record");
        outvcf.translate(&mut record);
        annot_var(&mut record, &mut ccres_bed, &mut ccre_to_gene_files, &mut biodata_files, &mut gene_annot_files);
        outvcf.write(&record).expect("Unable to write record");
    }
}
