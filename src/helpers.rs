//! Helper functions including parallel selective matrix multiplication and pseudo-inverse via singular value decomposition

use ndarray::prelude::*;

use statrs::distribution::{ContinuousCDF, StudentsT};

use std::fs::File;
use std::io::{self, prelude::*, BufReader, SeekFrom};
use std::io::{Error, ErrorKind};

// use crate::structs_and_traits::*;

/// Find the start position of the next line given the current position,`pos`.
/// Positions are coded as the nth UTF8 character count in the file counting the newline characters at the end of each line.
/// This is used in file splitting to allocate a chunk of the file to a single thread for parallel processing.
fn find_start_of_next_line(fname: &str, pos: u64) -> u64 {
    let mut out = pos;
    if out > 0 {
        let mut file = File::open(fname).expect("Error opening file.");
        let _ = file.seek(SeekFrom::Start(out));
        let mut reader = BufReader::new(file);
        let mut line = String::new();
        let _ = reader.read_line(&mut line).expect("Error reading file.");
        out = reader.stream_position().expect("Error navigating file.");
    }
    out
}

// Define chunks which respect loci groupings
pub fn define_chunks(
    loci_idx: &[usize],
    n_threads: &usize,
) -> io::Result<(Vec<usize>, Vec<usize>)> {
    let mut n_chunks = *n_threads; // can be more than the number of threads due to unequal division of the loci
    let l = loci_idx.len();
    let chunk_size = (l as f64 / n_chunks as f64).floor() as usize;
    // Define the indices of the indices of loci
    let mut vec_idx_all: Vec<usize> = if chunk_size < l {
        (0..l).step_by(chunk_size).collect()
    } else {
        vec![0, l]
    };
    // println!("l={:?}", l);
    // println!("n_chunks={:?}", n_chunks);
    // println!("chunk_size={:?}", chunk_size);
    // println!("vec_idx_all[vec_idx_all.len()-1]={:?}", vec_idx_all[vec_idx_all.len()-1]);
    if vec_idx_all.len() < (n_chunks + 1) {
        vec_idx_all.push(l - 1);
    } else {
        vec_idx_all.pop();
        vec_idx_all.push(l - 1);
    }
    // println!("vec_idx_all[vec_idx_all.len()-1]={:?}", vec_idx_all[vec_idx_all.len()-1]);
    // println!("vec_idx_all={:?}", vec_idx_all);
    n_chunks = vec_idx_all.len();
    let vec_idx_loci_idx_ini: Vec<usize> = vec_idx_all[0..(n_chunks - 1)].to_owned();
    let vec_idx_loci_idx_fin: Vec<usize> = vec_idx_all[1..n_chunks].to_owned();
    // println!("vec_idx_loci_idx_ini={:?}", vec_idx_loci_idx_ini);
    // println!("vec_idx_loci_idx_fin={:?}", vec_idx_loci_idx_fin);
    assert_eq!(vec_idx_loci_idx_ini.len(), vec_idx_loci_idx_fin.len());
    Ok((vec_idx_loci_idx_ini, vec_idx_loci_idx_fin))
}

/// Detect the cursor positions across the input file corresponding to the splits for parallel computation
pub fn find_file_splits(fname: &str, n_threads: &usize) -> io::Result<Vec<u64>> {
    let mut file = match File::open(fname) {
        Ok(x) => x,
        Err(_) => return Err(Error::new(ErrorKind::Other, "The input file: ".to_owned() + fname + " does not exist. Please make sure you are entering the correct filename and/or the correct path.")),
    };
    let _ = file.seek(SeekFrom::End(0));
    let mut reader = BufReader::new(file);
    let end = reader.stream_position().expect("Error navigating file.");
    let mut out = (0..end)
        .step_by((end as usize) / n_threads)
        .collect::<Vec<u64>>();
    out.push(end);
    for x in &mut out {
        *x = find_start_of_next_line(fname, x.to_owned());
    }
    out.dedup();
    Ok(out)
}

/// Round-up an `f64` to `n_digits` decimal points
pub fn sensible_round(x: f64, n_digits: usize) -> f64 {
    let factor = ("1e".to_owned() + &n_digits.to_string())
        .parse::<f64>()
        .expect("Error parsing String into f64.");
    (x * factor).round() / factor
}

/// Round-up an `f64` to `n_digits` decimal points and cast into a `String`
pub fn parse_f64_roundup_and_own(x: f64, n_digits: usize) -> String {
    let s = x.to_string();
    if s.len() < n_digits {
        return s;
    }
    sensible_round(x, n_digits).to_string()
}

/// Calculate the mean of a 1D array ignoring NaN
pub fn mean_array1_ignore_nan(
    x: &ArrayBase<ndarray::ViewRepr<&f64>, Dim<[usize; 1]>>,
) -> io::Result<f64> {
    // let sum = x.fold(0.0, |sum, &a| if a.is_nan() { sum } else { sum + a });
    // let counts = x.iter().filter(|&a| !a.is_nan()).count() as f64;
    // Ok(sum / counts)
    let mut sum = 0.0;
    let mut n = 0.0;
    for &x in x.iter() {
        if !x.is_nan() {
            sum += x;
            n += 1.0;
        }
    }
    Ok(sum / n)
}

/// Calculate the axis-wise means of an array while ignoring NaN
pub fn mean_axis_ignore_nan<D>(
    a: &Array<f64, D>,
    axis: usize,
) -> io::Result<Array<f64, <D>::Smaller>>
where
    D: ndarray::Dimension + ndarray::RemoveAxis,
{
    let sum: Array<f64, <D>::Smaller> =
        a.fold_axis(
            Axis(axis),
            0.0,
            |&sum, &x| {
                if x.is_nan() {
                    sum
                } else {
                    sum + x
                }
            },
        );
    let n: Array<f64, <D>::Smaller> = a.map_axis(Axis(axis), |x| {
        x.iter().filter(|&&y| !y.is_nan()).count() as f64
    });
    let out: Array<f64, <D>::Smaller> = sum / n;
    Ok(out)
}

/// Pearson's product moment correlation using only complete data across pair of vectors
pub fn pearsons_correlation_pairwise_complete(
    x: &ArrayBase<ndarray::ViewRepr<&f64>, Dim<[usize; 1]>>,
    y: &ArrayBase<ndarray::ViewRepr<&f64>, Dim<[usize; 1]>>,
) -> io::Result<(f64, f64)> {
    let n = x.len();
    if n != y.len() {
        return Err(Error::new(
            ErrorKind::Other,
            "Input vectors are not the same size.",
        ));
    }
    // Using pairs of values with non-missing data across the pair of vectors
    // Note that this may result in unreasonable correlations is used for a matrix, i.e. column vectors may be comparing different sets of rows
    let filtered_vectors: (Vec<f64>, Vec<f64>) = x
        .iter()
        .zip(y.iter())
        .filter(|&(&x, &y)| (!x.is_nan()) && (!y.is_nan()))
        .unzip();
    // println!("q.0={:?}; q.1={:?}", q.0, q.1);
    let x = Array1::from_vec(filtered_vectors.0);
    let y = Array1::from_vec(filtered_vectors.1);
    // Make sure we are handling NAN properly
    let mu_x = mean_array1_ignore_nan(&x.view())
        .expect("Error calculating the mean of x while ignoring NANs.");
    let mu_y = mean_array1_ignore_nan(&y.view())
        .expect("Error calculating the mean of y while ignoring NANs.");
    let x_less_mu_x = x
        .iter()
        .filter(|&x| !x.is_nan())
        .map(|x| x - mu_x)
        .collect::<Array1<f64>>();
    let y_less_mu_y = y
        .iter()
        .filter(|&y| !y.is_nan())
        .map(|y| y - mu_y)
        .collect::<Array1<f64>>();
    let x_less_mu_x_squared = x_less_mu_x.map(|x| x.powf(2.0));
    let y_less_mu_y_squared = y_less_mu_y.map(|y| y.powf(2.0));
    let numerator = (x_less_mu_x * y_less_mu_y).sum();
    let denominator = x_less_mu_x_squared.sum().sqrt() * y_less_mu_y_squared.sum().sqrt();
    let r_tmp = numerator / denominator;
    let r = match r_tmp.is_nan() {
        true => {
            if numerator == 0.0 {
                0.0
            } else {
                return Ok((f64::NAN, f64::NAN));
            }
        }
        false => r_tmp,
    };
    let sigma_r_denominator = (1.0 - r.powf(2.0)) / (n as f64 - 2.0);
    if sigma_r_denominator <= 0.0 {
        // Essentially no variance in r2, hence very significant
        return Ok((r, f64::EPSILON));
    }
    let sigma_r = sigma_r_denominator.sqrt();
    let t = r / sigma_r;
    let pval = if n > 2 {
        let d = StudentsT::new(0.0, 1.0, n as f64 - 2.0)
            .expect("Error defining Student's t-distribution.");
        2.00 * (1.00 - d.cdf(t.abs()))
    } else {
        f64::NAN
    };
    Ok((sensible_round(r, 7), pval))
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;
    #[test]
    fn test_helpers() {
        if cfg!(windows) {
            assert_eq!(
                find_start_of_next_line(&"./tests/test.sync".to_owned(), 77),
                79
            ); // line 1 has a total of 78 characters in windows and only 77 in unix, so the next line should start at the 79th character position
        } else {
            assert_eq!(
                find_start_of_next_line(&"./tests/test.sync".to_owned(), 77),
                78
            ); // 77th character is the newline character as line 1 has a total of 77 characters, so the next line should start at the 78th character position
        }
        assert_eq!(
            find_file_splits(&"./tests/test.sync".to_owned(), &2)
                .unwrap()
                .len(),
            3
        );
        assert_eq!(
            define_chunks(&vec![0, 1, 2, 3, 4, 5, 6], &2).unwrap(),
            (vec![0, 3], vec![3, 6])
        );
        assert_eq!(
            define_chunks(&vec![0, 1, 2, 3, 4, 5, 6], &3).unwrap(),
            (vec![0, 2, 4], vec![2, 4, 6])
        );
        assert_eq!(
            define_chunks(&vec![0, 1, 2, 3, 4, 5, 6], &5).unwrap(),
            (vec![0, 1, 2, 3, 4, 5], vec![1, 2, 3, 4, 5, 6])
        );
        assert_eq!(sensible_round(0.420000012435, 4), 0.42);
        assert_eq!(
            parse_f64_roundup_and_own(0.690000012435, 4),
            "0.69".to_owned()
        );
        let _a: Array2<f64> =
            Array2::from_shape_vec((5, 3), (0..15).map(|x| x as f64).collect::<Vec<f64>>())
                .unwrap();
        let _b: Array2<f64> = Array2::from_shape_vec(
            (5, 3),
            (0..15).map(|x| x as f64 / 2.0).collect::<Vec<f64>>(),
        )
        .unwrap();
        let _idx_w3: Vec<usize> = vec![1, 3, 4];
        let _idx_x2: Vec<usize> = vec![0, 2];
        let _idx_y2: Vec<usize> = vec![1, 3];
        let _idx_z2: Vec<usize> = vec![0, 1];

        let array1d: Array1<f64> = Array1::from_vec(vec![0.1, 0.2, 0.3, f64::NAN, 0.5]);
        assert_eq!(mean_array1_ignore_nan(&array1d.view()).unwrap(), 0.275);
        let mut array2d: Array2<f64> =
            Array2::from_shape_vec((2, 5), (0..10).map(|x| x as f64).collect::<Vec<f64>>())
                .unwrap();
        array2d[(0, 0)] = f64::NAN;
        println!("array2d={:?}", array2d);
        assert_eq!(
            Array1::from_shape_vec(5, vec![5.0, 3.5, 4.5, 5.5, 6.5]).unwrap(),
            mean_axis_ignore_nan(&array2d, 0).unwrap()
        );
        assert_eq!(
            Array1::from_shape_vec(2, vec![2.5, 7.0]).unwrap(),
            mean_axis_ignore_nan(&array2d, 1).unwrap()
        );
        let _fname = "./tests/test.tsv".to_owned();
        let _delimiter = ",".to_owned();
        let _chr_col = 0;
        let _pos_start_col = 1;
        let _pos_end_col = 1;
        let _data_start_col = 2;
        let _data_end_col = 6;

        let x = Array1::from_vec(vec![0.1, 0.2, 0.3, 0.4, 0.5]);
        let y = Array1::from_vec(vec![0.2, 0.4, 0.6, 0.8, 1.0]);
        let z = Array1::from_vec(vec![0.0, 0.1, 0.0, 0.1, 0.0]);
        let corr1 = pearsons_correlation_pairwise_complete(&x.view(), &y.view()).unwrap();
        let corr2 = pearsons_correlation_pairwise_complete(&x.view(), &z.view()).unwrap();
        println!("corr1={:?}", corr1);
        println!("corr2={:?}", corr2);
        assert!(corr1.0 == 1.00);
        assert!(corr1.1 < 0.0001);
        assert!(corr2.0 == 0.00);
        assert!(corr2.1 == 1.00);
    }
}
