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
fn find_start_of_next_line(fname: &String, pos: u64) -> u64 {
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

/// Detect the cursor positions across the input file corresponding to the splits for parallel computation
pub fn find_file_splits(fname: &String, n_threads: &usize) -> io::Result<Vec<u64>> {
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
    for i in 0..out.len() {
        out[i] = find_start_of_next_line(fname, out[i]);
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

// /// Multi-threaded selective matrix multiplication variant: $AB$
// /// Matrix multiplication of two 2-dimensional arrays where only the specified rows and columns are used
// /// This is an attempt at multi-threadded selective matrix multiplication while minimising memory allocation by using pointers to arrays and a subset of its data, instead of coping data.
// pub fn multiply_views_xx(
//     a: &Array2<f64>,
//     b: &Array2<f64>,
//     a_rows: &Vec<usize>,
//     a_cols: &Vec<usize>,
//     b_rows: &Vec<usize>,
//     b_cols: &Vec<usize>,
// ) -> io::Result<Array2<f64>> {
//     let n = a_rows.len();
//     let m = b_cols.len();
//     if a_cols.len() != b_rows.len() {
//         return Err(Error::new(
//             ErrorKind::Other,
//             "The two matrices are incompatible.",
//         ));
//     }
//     let mut out: Array2<f64> = Array2::zeros((n, m));
//     let a_rows_mat = Array2::from_shape_vec((m, n), a_rows.repeat(m))
//         .unwrap()
//         .reversed_axes();
//     let b_cols_mat = Array2::from_shape_vec((n, m), b_cols.repeat(n)).unwrap();
//     Zip::from(&mut out)
//         .and(&a_rows_mat)
//         .and(&b_cols_mat)
//         .par_for_each(|x, &a_i, &b_j| {
//             for k in 0..a_cols.len() {
//                 let a_j = a_cols[k];
//                 let b_i = b_rows[k];
//                 *x += a[(a_i, a_j)] * b[(b_i, b_j)];
//             }
//         });
//     Ok(out)
// }

// /// Multi-threaded selective matrix multiplication variant: $A^{T}B$
// pub fn multiply_views_xtx(
//     a: &Array2<f64>,
//     b: &Array2<f64>,
//     a_rows: &Vec<usize>,
//     a_cols: &Vec<usize>,
//     b_rows: &Vec<usize>,
//     b_cols: &Vec<usize>,
// ) -> io::Result<Array2<f64>> {
//     let n = a_cols.len(); // reversed a
//     let m = b_cols.len();
//     if a_rows.len() != b_rows.len() {
//         // reversed a
//         return Err(Error::new(
//             ErrorKind::Other,
//             "The two matrices are incompatible.",
//         ));
//     }
//     let mut out: Array2<f64> = Array2::zeros((n, m));
//     let a_cols_mat = Array2::from_shape_vec((m, n), a_cols.repeat(m))
//         .unwrap()
//         .reversed_axes();
//     let b_cols_mat = Array2::from_shape_vec((n, m), b_cols.repeat(n)).unwrap();
//     Zip::from(&mut out)
//         .and(&a_cols_mat)
//         .and(&b_cols_mat)
//         .par_for_each(|x, &a_j, &b_j| {
//             for k in 0..a_rows.len() {
//                 let a_i = a_rows[k];
//                 let b_i = b_rows[k];
//                 *x += a[(a_i, a_j)] * b[(b_i, b_j)];
//             }
//         });
//     Ok(out)
// }

// /// Multi-threaded selective matrix multiplication variant: $AB^{T}$
// pub fn multiply_views_xxt(
//     a: &Array2<f64>,
//     b: &Array2<f64>,
//     a_rows: &Vec<usize>,
//     a_cols: &Vec<usize>,
//     b_rows: &Vec<usize>,
//     b_cols: &Vec<usize>,
// ) -> io::Result<Array2<f64>> {
//     let n = a_rows.len();
//     let m = b_rows.len(); // reversed b
//     if a_cols.len() != b_cols.len() {
//         // reversed b
//         return Err(Error::new(
//             ErrorKind::Other,
//             "The two matrices are incompatible.",
//         ));
//     }
//     let mut out: Array2<f64> = Array2::zeros((n, m));
//     let a_rows_mat = Array2::from_shape_vec((m, n), a_rows.repeat(m))
//         .unwrap()
//         .reversed_axes();
//     let b_rows_mat = Array2::from_shape_vec((n, m), b_rows.repeat(n)).unwrap();
//     Zip::from(&mut out)
//         .and(&a_rows_mat)
//         .and(&b_rows_mat)
//         .par_for_each(|x, &a_i, &b_i| {
//             for k in 0..a_cols.len() {
//                 let a_j = a_cols[k];
//                 let b_j = b_cols[k];
//                 *x += a[(a_i, a_j)] * b[(b_i, b_j)];
//             }
//         });
//     Ok(out)
// }

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

/// Extract the coordinates of each sliding window (can accommodate redundant and non-redundant loci)
pub fn define_sliding_windows(
    loci_chr: &Vec<String>,
    loci_pos: &Vec<u64>,
    window_size_bp: &u64,
    window_slide_size_bp: &u64,
    min_loci_per_window: &u64,
) -> io::Result<(Vec<usize>, Vec<usize>)> {
    assert_eq!(loci_chr.len(), loci_pos.len());
    let l = loci_chr.len();
    // Indices, chromosome names, and positions of the start and end of the window, respectively (will be filtered to remove redundant tails to remove windows which are complete subsets of a bigger window near the end of chromosomes or scaffolds)
    let mut idx_head: Vec<usize> = vec![0];
    let mut idx_tail = idx_head.clone();
    let mut chr_head: Vec<String> = vec![loci_chr[0].to_owned()];
    let mut chr_tail = chr_head.clone();
    let mut pos_head: Vec<u64> = vec![loci_pos[0]];
    let mut pos_tail = pos_head.clone();
    // Number of genotyped loci included per window
    let mut cov: Vec<u64> = vec![1];
    // A boolean to mark whether the start of the next window has been found before the end of the current window has been reached
    let mut marker_next_window_head: bool = false;
    // The index of the sart of the next window
    let mut idx_next_head: usize = 0;
    // Index of the current locus
    let mut i: usize = 1;
    // We iterate across the genome until we reach the last locus (may not be consecutive as the start of the next window may be within the body of the current window)
    while i < l {
        let chr = loci_chr[i].to_owned();
        let pos = loci_pos[i];
        // println!("i={:?}", i);
        // println!("idx_head={:?}", idx_head);
        // println!("idx_tail={:?}", idx_tail);
        // Did we reach the end of the chromosome or the end of the window according to window size?
        if (&chr != chr_head.last().expect("Error extracting the last character of chr_head.")) | (pos > (pos_head.last().expect("Error extracting the last character of pos_head.") + window_size_bp))
        {
            // If we found the start of the next window in body of the current (ending) window then,
            //  we use the next window head as the start of the next slide not the end of the window,
            //  otherwise we use the end of the window.
            i = if marker_next_window_head {
                idx_next_head
            } else {
                i
            };
            let chr = loci_chr[i].to_owned();
            let pos = loci_pos[i];
            // Do we have the minimum number of required loci in the current window?
            if cov.last().expect("Error extracting the last value of cov.") >= min_loci_per_window {
                // If we have enough loci covered in the current (ending) window:
                // We also add the details of the start of the next window
                idx_head.push(i);
                idx_tail.push(i);
                chr_head.push(chr.to_owned());
                chr_tail.push(chr.to_owned());
                pos_head.push(pos);
                pos_tail.push(pos);
                cov.push(1);
            } else {
                // If we did no have enough loci covered in the current (ending) window:
                // We ditch the current (ending) window and replace it with the start of the next window
                let i_ = idx_head.len() - 1;
                idx_head[i_] = i;
                chr_head[i_] = chr;
                pos_head[i_] = pos;
                cov[i_] = 1;
            }
            // Reset the marker for the start of the next window
            marker_next_window_head = false;
        } else {
            // If we have yet to reach the end of the current window or the end of the chromosome,
            // then we just replace the tail of the current window with the current locus
            // and add the another locus to the coverage counter.
            let i_ = idx_tail.len() - 1;
            idx_tail[i_] = i;
            chr_tail[i_] = chr;
            pos_tail[i_] = pos;
            cov[i_] += 1;
            // We also check if we have reached the start of the next window and note the index if we have
            if !marker_next_window_head & (pos >= (pos_head.last().expect("Error extracting the last character of pos_head.") + window_slide_size_bp))
            {
                marker_next_window_head = true;
                idx_next_head = i;
            }
        }
        // Move to the next locus
        i += 1;
    }
    // Remove redundant tails
    assert_eq!(idx_head.len(), idx_tail.len());
    let n = idx_head.len();
    let mut out_idx_head: Vec<usize> = vec![idx_head[0]];
    let mut out_idx_tail: Vec<usize> = vec![idx_tail[0]];
    for i in 1..n {
        // println!("out_idx_tail={:?}", out_idx_tail);
        if &idx_tail[i] != out_idx_tail.last().expect("Error extracting the last value of out_idx_tail.") {
            out_idx_head.push(idx_head[i]);
            out_idx_tail.push(idx_tail[i]);
        }
    }
    // println!("#################################");
    // println!("loci_chr={:?}", loci_chr);
    // println!("loci_pos={:?}", loci_pos);
    // println!("window_size_bp={:?}", window_size_bp);
    // println!("min_loci_per_window={:?}", min_loci_per_window);
    // println!("cov={:?}", cov);
    // println!("idx_head={:?}", idx_head);
    // println!("idx_tail={:?}", idx_tail);
    // println!("out_idx_head={:?}", out_idx_head);
    // println!("out_idx_tail={:?}", out_idx_tail);
    Ok((out_idx_head, out_idx_tail))
}

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
        .filter(|&(&x, &y)| (!x.is_nan()) & (!y.is_nan()))
        .unzip();
    // println!("q.0={:?}; q.1={:?}", q.0, q.1);
    let x = Array1::from_vec(filtered_vectors.0);
    let y = Array1::from_vec(filtered_vectors.1);
    // Make sure we are handling NAN properly
    let mu_x = mean_array1_ignore_nan(&x.view()).expect("Error calculating the mean of x while ignoring NANs.");
    let mu_y = mean_array1_ignore_nan(&y.view()).expect("Error calculating the mean of y while ignoring NANs.");
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
        true => return Ok((f64::NAN, f64::NAN)),
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
        let d = StudentsT::new(0.0, 1.0, n as f64 - 2.0).expect("Error defining Student's t-distribution.");
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
                find_start_of_next_line(&"./tests/test.pileup".to_owned(), 56),
                58
            ); // line 1 has a total of 57 characters in windows and only 56 in unix, so the next line should start at the 58th character position
        } else {
            assert_eq!(
                find_start_of_next_line(&"./tests/test.pileup".to_owned(), 56),
                57
            ); // 56th character is the newline character as line 1 has a total of 56 characters, so the next line should start at the 57th character position
        }
        assert_eq!(
            find_file_splits(&"./tests/test.pileup".to_owned(), &2)
                .unwrap()
                .len(),
            3
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
        let _fname = "./tests/test.csv".to_owned();
        let _delimiter = ",".to_owned();
        let _chr_col = 0;
        let _pos_start_col = 1;
        let _pos_end_col = 1;
        let _data_start_col = 2;
        let _data_end_col = 6;

        // Define sliding windows (non-redundant loci, i.e. per locus list with alleles ID removed)
        let loci_chr: Vec<String> = [
            "chr1", "chr1", "chr1", "chr1", "chr2", "chr2", "chr2", "chr2", "chr3",
        ]
        .iter()
        .map(|&x| x.to_owned())
        .collect();
        let loci_pos: Vec<u64> = vec![123, 174, 220, 254, 55, 56, 100, 500, 765];
        let window_size_bp: u64 = 100;
        let window_slide_size_bp: u64 = 50;
        let min_loci_per_window: u64 = 1;
        let (windows_idx_head, windows_idx_tail) = define_sliding_windows(
            &loci_chr,
            &loci_pos,
            &window_size_bp,
            &window_slide_size_bp,
            &min_loci_per_window,
        )
        .unwrap();
        println!("windows_idx_head={:?}", windows_idx_head);
        println!("windows_idx_tail={:?}", windows_idx_tail);
        assert_eq!(windows_idx_head, vec![0, 1, 4, 7, 8]); // filtered out window start:2-end:3 which is a complete subset of window start:1-end:3
        assert_eq!(windows_idx_tail, vec![2, 3, 6, 7, 8]); // filtered out window start:2-end:3 which is a complete subset of window start:1-end:3
                                                           // Define sliding windows (redundant loci, i.e. per allele per locus)
        let loci_chr: Vec<String> = ["X", "X", "X", "Y", "Y"]
            .iter()
            .map(|&x| x.to_owned())
            .collect();
        let loci_pos: Vec<u64> = vec![123, 123, 123, 456, 456];
        let window_size_bp: u64 = 100;
        let window_slide_size_bp: u64 = 50;
        let min_loci_per_window: u64 = 1;
        let (windows_idx_head, windows_idx_tail) = define_sliding_windows(
            &loci_chr,
            &loci_pos,
            &window_size_bp,
            &window_slide_size_bp,
            &min_loci_per_window,
        )
        .unwrap();
        println!("windows_idx_head={:?}", windows_idx_head);
        println!("windows_idx_tail={:?}", windows_idx_tail);
        assert_eq!(windows_idx_head, vec![0, 3]); // filtered out window start:2-end:3 which is a complete subset of window start:1-end:3
        assert_eq!(windows_idx_tail, vec![2, 4]); // filtered out window start:2-end:3 which is a complete subset of window start:1-end:3

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
