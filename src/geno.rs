use ndarray::prelude::*;
use std::fs::File;
use std::io::{self, prelude::*, BufReader, Error, ErrorKind, SeekFrom};
use std::str;
use std::sync::{Arc, Mutex};

use crate::helpers::*;
use crate::structs_and_traits::*;

impl Parse<LocusFrequencies> for String {
    fn lparse(&self) -> io::Result<Box<LocusFrequencies>> {
        // Ignore commented-out lines (i.e. '#' => 35)
        if self.as_bytes()[0] == 35_u8 {
            return Err(Error::new(ErrorKind::Other, "Commented out line"));
        }
        let vec_line: Vec<&str> = self.split('\t').collect();
        let vec_line: Vec<&str> = if vec_line.len() == 1 {
            self.split(',').collect()
        } else {
            vec_line
        };
        let vec_line: Vec<&str> = if vec_line.len() == 1 {
            self.split(';').collect()
        } else {
            vec_line
        };
        // println!("vec_line={:?}", vec_line);
        let l = vec_line.len();
        let n = l - 3;
        let chromosome: String = vec_line[0].to_owned();
        let position = match vec_line[1].parse::<u64>() {
            Ok(x) => x,
            Err(_) => {
                return Err(Error::new(
                    ErrorKind::Other,
                    "Please check format of the file: position is not and integer.",
                ))
            }
        };
        let alleles_vector: Vec<String> = vec![vec_line[2].to_owned()];
        let matrix: Array2<f64> = Array2::from_shape_vec(
            (n, 1),
            vec_line[3..l]
                .iter()
                .map(|x| {
                    match x.parse::<f64>() {
                        Ok(x) => x,
                        Err(_) => f64::NAN,
                    }
                })
                .collect::<Vec<f64>>(),
        )
        .expect("Error parsing the allele frequency table text file within the lparse() method for parsing String into LocusFrequencies struct.");
        let freq_line = LocusFrequencies {
            chromosome,
            position,
            alleles_vector,
            matrix,
        };
        // println!("freq_line={:?}", freq_line);
        Ok(Box::new(freq_line))
    }
}

impl LoadAll for FileGeno {
    fn per_chunk_load(
        &self,
        start: &u64,
        end: &u64,
        _filter_stats: &FilterStats,
        _keep_p_minus_1: bool,
    ) -> io::Result<(Vec<LocusFrequencies>, Vec<LocusCounts>)> {
        // Input syn file
        let fname = self.filename.clone();
        // Prepare output vectors
        let mut freq: Vec<LocusFrequencies> = Vec::new();
        let cnts: Vec<LocusCounts> = Vec::new(); // Empty and will remain empty as each line corresponds to just an allele of a locus
                                                 // Input file chunk
        let file = File::open(fname.clone()).expect("Error opening the allele frequency table text file within the per_chunk_load() method for FileGeno struct.");
        let mut reader = BufReader::new(file);
        // Navigate to the start of the chunk
        let mut i: u64 = *start;
        reader.seek(SeekFrom::Start(*start)).expect("Error navigating across the allele frequency table text file within the per_chunk_load() method for FileGeno struct.");
        // Read and parse until the end of the chunk
        while i < *end {
            // Instantiate the line
            let mut line = String::new();
            // Read the line which automatically moves the cursor position to the next line
            let _ = reader.read_line(&mut line).expect("Error reading the allele frequency table text file within the per_chunk_load() method for FileGeno struct.");
            // Find the new cursor position
            i = reader.stream_position().expect("Error navigating across the allele frequency table text file within the per_chunk_load() method for FileGeno struct.");
            // Remove trailing newline character in Unix-like (\n) and Windows (\r)
            if line.ends_with('\n') {
                line.pop();
                if line.ends_with('\r') {
                    line.pop();
                }
            }
            // Parse the geno line
            let allele_freqs: LocusFrequencies = match line.lparse() {
                Ok(x) => *x,
                Err(x) => match x.kind() {
                    ErrorKind::Other => continue,
                    _ => {
                        return Err(Error::new(
                            ErrorKind::Other,
                            "T_T Input sync file error, i.e. '".to_owned()
                                + &fname
                                + "' at line with the first 20 characters as: "
                                + &line[0..20]
                                + ".",
                        ))
                    }
                },
            };
            freq.push(allele_freqs);
        }
        Ok((freq, cnts))
    }

    fn load(
        &self,
        filter_stats: &FilterStats,
        keep_p_minus_1: bool,
        n_threads: &usize,
    ) -> io::Result<(Vec<LocusFrequencies>, Vec<LocusCounts>)> {
        let fname = self.filename.clone();
        // Find the positions whereto split the file into n_threads pieces
        let chunks = find_file_splits(&fname, n_threads).expect("Error splitting the allele frequency table file format given the number of threads suppplied within load() method for FileGeno struct.");
        let n_threads = chunks.len() - 1;
        println!("Chunks: {:?}", chunks);
        // Tuple arguments of pileup2sync_chunks
        // Instantiate thread object for parallel execution
        let mut thread_objects = Vec::new();
        // Vector holding all returns from pileup2sync_chunk()
        let thread_ouputs_freq: Arc<Mutex<Vec<LocusFrequencies>>> =
            Arc::new(Mutex::new(Vec::new())); // Mutated within each thread worker
        let thread_ouputs_cnts: Arc<Mutex<Vec<LocusCounts>>> = Arc::new(Mutex::new(Vec::new())); // Mutated within each thread worker
                                                                                                 // Making four separate threads calling the `search_for_word` function
        for i in 0..n_threads {
            // Clone pileup2sync_chunk parameters
            let self_clone = self.clone();
            let start = chunks[i];
            let end = chunks[i + 1];
            let filter_stats = filter_stats.clone();
            let thread_ouputs_freq_clone = thread_ouputs_freq.clone(); // Mutated within the current thread worker
            let thread_ouputs_cnts_clone = thread_ouputs_cnts.clone(); // Mutated within the current thread worker
            let thread = std::thread::spawn(move || {
                let (mut freq, mut cnts) = self_clone
                    .per_chunk_load(&start, &end, &filter_stats, keep_p_minus_1)
                    .expect(
                        "Error calling per_chunk_load() within load() method for FileGeno struct.",
                    );
                thread_ouputs_freq_clone
                    .lock()
                    .expect("Thread error within load() method for FileGeno struct.")
                    .append(&mut freq);
                thread_ouputs_cnts_clone
                    .lock()
                    .expect("Thread error within load() method for FileGeno struct.")
                    .append(&mut cnts);
            });
            thread_objects.push(thread);
        }
        // Waiting for all threads to finish
        for thread in thread_objects {
            thread.join().expect("Unknown thread error occured.");
        }
        // Extract output filenames from each thread into a vector and sort them
        let mut freq: Vec<LocusFrequencies> = Vec::new();
        let cnts: Vec<LocusCounts> = Vec::new(); // Empty and will remain empty as each line corresponds to just an allele of a locus
        for x in thread_ouputs_freq.lock().expect("Error unlocking the threads after multi-threaded execution of per_chunk_load() within load() method for FileGeno struct.").iter() {
            freq.push(x.clone());
        }
        freq.sort_by(|a, b| {
            a.chromosome
                .cmp(&b.chromosome)
                .then(a.position.cmp(&b.position))
        });
        Ok((freq, cnts))
    }

    fn convert_into_genotypes_and_phenotypes(
        &self,
        filter_stats: &FilterStats,
        keep_p_minus_1: bool,
        n_threads: &usize,
    ) -> io::Result<GenotypesAndPhenotypes> {
        // No filtering! Just loading the allele frequency data
        // Extract pool names
        let file: File = File::open(self.filename.clone()).expect("Error opening the allele frequency table text file within the convert_into_genotypes_and_phenotypes() method for FileGeno struct.");
        let reader = io::BufReader::new(file);
        let mut header: String = match reader.lines().next() {
            Some(x) => x.expect("Error reading the allele frequency table text file within the convert_into_genotypes_and_phenotypes() method for FileGeno struct."),
            None => return Err(Error::new(ErrorKind::Other, "No header line found.")),
        };
        if header.ends_with('\n') {
            header.pop();
            if header.ends_with('\r') {
                header.pop();
            }
        }
        let vec_header: Vec<&str> = header.split('\t').collect();
        let vec_header: Vec<&str> = if vec_header.len() == 1 {
            header.split(',').collect()
        } else {
            vec_header
        };
        let vec_header: Vec<&str> = if vec_header.len() == 1 {
            header.split(';').collect()
        } else {
            vec_header
        };

        let pool_names: Vec<String> = vec_header[3..vec_header.len()]
            .iter()
            .map(|&x| x.to_owned())
            .collect();
        // Load allele frequencies
        let (freqs, _cnts) = self.load(filter_stats, keep_p_minus_1, n_threads).expect("Error calling load() within the convert_into_genotypes_and_phenotypes() method for FileGeno struct.");
        let n = freqs[0].matrix.nrows();
        assert_eq!(
            n,
            pool_names.len(),
            "Header names and allele frequency data does not have the same number of samples."
        );
        let mut p = freqs.len(); // total number of alleles across all loci
        p += 1; // include the intercept
        let mut chromosome: Vec<String> = Vec::with_capacity(p);
        chromosome.push("intercept".to_owned());
        let mut position: Vec<u64> = Vec::with_capacity(p);
        position.push(0);
        let mut allele: Vec<String> = Vec::with_capacity(p);
        allele.push("intercept".to_owned());
        // Populate the GenotypeAndPhenotypes struct fields
        let mut l: usize = 0; // loci counter
        let mut mat: Array2<f64> = Array2::from_elem((n, p), 1.0);
        // Count the number of loci (Note: assumes the loci are sorted) and extract the loci coordinates
        let mut loci_idx: Vec<usize> = vec![];
        let mut loci_chr: Vec<String> = vec![];
        let mut loci_pos: Vec<u64> = vec![];
        for j in 1..p {
            let f = &freqs[j - 1];
            chromosome.push(f.chromosome.clone());
            position.push(f.position);
            allele.push(f.alleles_vector[0].clone());
            if (chromosome[j - 1] != chromosome[j]) || (position[j - 1] != position[j]) {
                l += 1;
                loci_idx.push(j);
                loci_chr.push(chromosome[j].to_owned());
                loci_pos.push(position[j]);
            }
            for i in 0..n {
                mat[(i, j)] = f.matrix[(i, 0)];
            }
        }
        // Add the last allele of the last locus
        loci_idx.push(p);
        loci_chr.push(chromosome.last().expect("Error push chromosome within the convert_into_genotypes_and_phenotypes() method for FileGeno struct.").to_owned());
        loci_pos.push(position.last().expect("Error push position within the convert_into_genotypes_and_phenotypes() method for FileGeno struct.").to_owned());
        // Add alternative alleles if the allele frequencies per locus do not add up to 1.00 (~or if only one allele per locus is present~)
        // Count how many allele we have to add
        for j in 0..l {
            let idx_ini = loci_idx[j];
            let idx_fin = loci_idx[j + 1];
            let _n_alleles = idx_fin - idx_ini;
            let mut freq_sum_less_than_one = false;
            for i in 0..n {
                if mat.slice(s![i, idx_ini..idx_fin]).sum() < 1.0 {
                    freq_sum_less_than_one = true;
                    break;
                }
            }
            // if (n_alleles == 1) || freq_sum_less_than_one {
            if freq_sum_less_than_one {
                p += 1;
            }
        }
        let mut chromosome_new: Vec<String> = Vec::with_capacity(p);
        chromosome_new.push("intercept".to_owned());
        let mut position_new: Vec<u64> = Vec::with_capacity(p);
        position_new.push(0);
        let mut allele_new: Vec<String> = Vec::with_capacity(p);
        allele_new.push("intercept".to_owned());
        let mut mat_new: Array2<f64> = Array2::from_elem((n, p), 1.0);
        let mut j = 1;
        for j_orig in 0..l {
            let idx_ini = loci_idx[j_orig];
            let idx_fin = loci_idx[j_orig + 1];
            let n_alleles = idx_fin - idx_ini;
            // println!("j={}; j_orig={}; idx_ini={}; idx_fin={}; n_alleles={}", j, j_orig, idx_ini, idx_fin, n_alleles);
            for a in idx_ini..idx_fin {
                chromosome_new.push(chromosome[a].to_owned());
                position_new.push(position[a]);
                allele_new.push(allele[a].to_owned());
                for i in 0..n {
                    mat_new[(i, j)] = mat[(i, a)];
                }
                j += 1;
            }
            let mut freq_sum_less_than_one = false;
            for i in 0..n {
                if mat.slice(s![i, idx_ini..idx_fin]).sum() < 1.0 {
                    freq_sum_less_than_one = true;
                    break;
                }
            }
            if (n_alleles == 1) || freq_sum_less_than_one {
                chromosome_new.push(chromosome[idx_ini].to_owned());
                position_new.push(position[idx_ini]);
                allele_new.push("U".to_owned()); // unknown alternative allele
                for i in 0..n {
                    let alt = 1.00 - mat.slice(s![i, idx_ini..idx_fin]).sum();
                    mat_new[(i, j)] = alt;
                }
                j += 1;
            }
        }
        Ok(GenotypesAndPhenotypes {
            chromosome: chromosome_new,
            position: position_new,
            allele: allele_new,
            intercept_and_allele_frequencies: mat_new,
            phenotypes: Array2::from_shape_vec((n, 1), vec![f64::NAN; n]).expect("Error generating dummy phenotype data within the convert_into_genotypes_and_phenotypes() method for FileGeno struct."),
            pool_names,
            coverages: Array2::from_elem((n, l), 1_000.0),
        })
    }
}

pub fn load_geno<'a, 'b>(
    fname: &'a str,
    filter_stats: &'b mut FilterStats,
    _fname_out_prefix: &'a str,
    _rand_id: &'a str,
    n_threads: &'a usize,
) -> io::Result<(GenotypesAndPhenotypes, &'b FilterStats)> {
    // Extract pool names from the txt file
    let file: File = File::open(fname).expect("Error reading the allele frequency table file.");
    let reader = io::BufReader::new(file);
    let mut header: String = reader
        .lines()
        .next()
        .expect("Error reading the allele frequency table file.")
        .expect("Please check the format of the allele frequency table text file.");
    if header.ends_with('\n') {
        header.pop();
        if header.ends_with('\r') {
            header.pop();
        }
    }
    let vec_header: Vec<&str> = header.split('\t').collect();
    let vec_header: Vec<&str> = if vec_header.len() == 1 {
        header.split(',').collect()
    } else {
        vec_header
    };
    let vec_header: Vec<&str> = if vec_header.len() == 1 {
        header.split(';').collect()
    } else {
        vec_header
    };
    let pool_names: Vec<String> = vec_header[3..vec_header.len()]
        .iter()
        .map(|&x| x.to_owned())
        .collect();
    let n = pool_names.len();
    // If a single pool size was supplied then we are assuming the same sizes across all pools
    if filter_stats.pool_sizes.len() == 1 {
        filter_stats.pool_sizes = vec![filter_stats.pool_sizes[0]; n];
    }
    assert_eq!(
        filter_stats.pool_sizes.len(),
        n,
        "Error: the number of pools and the pool sizes do not match."
    );
    let file_geno = FileGeno {
        filename: fname.to_owned(),
    };
    Ok((file_geno
        .convert_into_genotypes_and_phenotypes(filter_stats, false, n_threads)
        .expect("Error parsing the genotype data (extracted from allele frequency table text file) via convert_into_genotypes_and_phenotypes() method within impute()."), 
        filter_stats))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_geno() {
        let line: String = "chr_1\t12345\tA\t0.25\t0.1\t0.9\t0.45\t0.9\t0.85".to_owned();
        let _genotypes_and_phenotypes: LocusFrequencies = *line.lparse().unwrap();

        let file_geno = FileGeno {
            filename: "tests/test.txt".to_owned(),
        };
        let filter_stats = FilterStats {
            remove_ns: true,
            min_quality: 0.005,
            min_coverage: 1,
            min_allele_frequency: 0.005,
            max_missingness_rate: 0.0,
            pool_sizes: vec![20., 20., 20., 20., 20.],
        };
        let n_threads = 8;

        let genotypes_and_phenotype = file_geno
            .convert_into_genotypes_and_phenotypes(&filter_stats, false, &n_threads)
            .unwrap();
        // println!("genotypes_and_phenotype={:?}", genotypes_and_phenotype);
        assert!(
            genotypes_and_phenotype.pool_names
                == vec![
                    "G1".to_owned(),
                    "G2".to_owned(),
                    "G3".to_owned(),
                    "G5".to_owned()
                ]
        );
        assert!(genotypes_and_phenotype.chromosome[0] == "intercept".to_owned());
        assert!(genotypes_and_phenotype.chromosome[1] == "Chromosome1".to_owned());
        assert!(genotypes_and_phenotype.chromosome[2] == "Chromosome1".to_owned());
        assert!(genotypes_and_phenotype.chromosome[3] == "Chromosome1".to_owned());
        assert!(genotypes_and_phenotype.chromosome[4] == "Chromosome1".to_owned());
        assert!(genotypes_and_phenotype.position[0] == 0);
        assert!(genotypes_and_phenotype.position[1] == 456527);
        assert!(genotypes_and_phenotype.position[2] == 456527);
        assert!(genotypes_and_phenotype.position[3] == 1133215);
        assert!(genotypes_and_phenotype.position[4] == 1133215);
        assert!(genotypes_and_phenotype.allele[0] == "intercept".to_owned());
        assert!(genotypes_and_phenotype.allele[1] == "T".to_owned());
        assert!(genotypes_and_phenotype.allele[2] == "C".to_owned());
        assert!(genotypes_and_phenotype.allele[3] == "A".to_owned());
        assert!(genotypes_and_phenotype.allele[4] == "C".to_owned());
        assert!(genotypes_and_phenotype.intercept_and_allele_frequencies[(0, 0)] == 1.00);
        assert!(genotypes_and_phenotype.intercept_and_allele_frequencies[(0, 1)] == 0.238095);
        assert!(genotypes_and_phenotype.intercept_and_allele_frequencies[(0, 2)] == 0.761905);
        assert!(genotypes_and_phenotype.intercept_and_allele_frequencies[(1, 3)] == 0.0);
        assert!(genotypes_and_phenotype.intercept_and_allele_frequencies[(2, 4)] == 1.0);
        // assert!(0 == 1);
    }
}
