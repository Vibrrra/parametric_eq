use std::f32::consts::PI;

use slint::SharedString;

/// Function to generate logarithmically spaced frequency points
fn log_space(start: f32, end: f32, num_points: usize) -> Vec<f32> {
    let log_start = start.ln();
    let log_end = end.ln();
    let step = (log_end - log_start) / (num_points as f32 - 1.0);

    (0..num_points).map(|i| (log_start + step * i as f32).exp()).collect()
}

fn log10_space(start: f32, end: f32, num_points: usize) -> Vec<f32> {
    let log_start = start.log10();
    let log_end = end.log10();
    let step = (log_end - log_start) / (num_points as f32 - 1.0);

    (0..num_points).map(|i| 10f32.powf(log_start + step * i as f32)).collect()
}


/// similar function to the one above, but it takes a specific frequency 
/// that calculates the x value of the frequency in the log space
/// and returns the index of the frequency in the log space
/// # Arguments
/// * `freq` - The frequency to find in the log space
/// * `start` - The start of the log space
///     
/// 
/// # Returns
///    

pub fn f_to_lin(y0: f32, y1: f32, x0: f32, x1: f32, y: f32) -> f32 {
    (y.log10() - y0.log10()) * (x1 - x0) / (y1.log10() - y0.log10()) + x0
}

// inverse function from above
pub fn lin_to_f(y0: f32, y1: f32, x0: f32, x1: f32, x: f32) -> f32 {
    10.0f32.powf( (x - x0) * (y1.log10() - y0.log10()) / (x1 - x0) + y0.log10() )
}

fn find_freq_index(freq: f32, start: f32, num_points: usize) -> f32 {
    let log_start = start.ln();
    let log_freq = freq.ln();
    let step = (log_freq - log_start) / (num_points as f32 - 1.0);
    let index = ((log_freq - log_start) / step) ;
    index
}

fn find_freq_index2(freq: f32, start: f32, num_points: usize) -> f32 {
    let log_start = start.log10();
    let log_freq = freq.log10();
    let step = (log_freq - log_start) / (num_points as f32 - 1.0);
    let index = ((log_freq - log_start) / step) ;
    index
}


// A function that computes a logarithmically spaced frequency grid with skewing
// The skew factor of the logarithmic scale is calculated by the formula:
// - basic transform really is min + ( max - min ) * pow( position, 1 / skewFactor )
// 

fn log_space_skewed(start: f32, end: f32, num_points: usize, skew_f: f32) -> Vec<f32> {
    let mut frequencies = Vec::with_capacity(num_points);
    for i in 0..num_points {
        let position = i as f32 / (num_points as f32 - 1.0);
        let skewed_position = position.powf(1.0/skew_f);
        let frequency = start + (end - start) * skewed_position;
        frequencies.push(frequency);
    }
    frequencies
}

// similar to the log_space_skewed function, but it takes a single input
pub fn lin_to_f_skew(start: f32, end: f32, num_points: usize, skew_f: f32, x: f32) -> f32 {
    let position = x / (num_points as f32 - 1.0);
    let skewed_position = position.powf(1.0/skew_f);
    start + (end - start) * skewed_position
}

// inverse of the above function
pub fn f_to_lin_skew(start: f32, end: f32, num_points: usize, skew_f: f32, y: f32) -> f32 {
    let position = (y - start) / (end - start);
    let skewed_position = position.powf(skew_f);
    (num_points as f32 - 1.0) * skewed_position
}

/// Function to compute the magnitude response of a biquad filter
/// # Arguments
/// * `b` - A tuple with the biquad numerator coefficients (b0, b1, b2)
/// * `a` - A tuple with the biquad denominator coefficients (a0, a1, a2)
/// * `sample_rate` - The sample rate of the system (e.g., 44100 Hz)
/// * `num_points` - The number of frequency points to calculate
/// 
/// # Returns
/// A tuple of two Vecs: 
/// - A Vec containing the magnitude response in dB
/// - A Vec containing the logarithmically spaced frequency points
pub fn calculate_biquad_response(
    b: (f32, f32, f32),
    a: (f32, f32, f32),
    sample_rate: f32,
    num_points: usize,
) -> (Vec<f32>, Vec<f32>) {
    let (b0, b1, b2) = b;
    let (a0, a1, a2) = a;

    // Logarithmically spaced frequency points between 20 Hz and 20 kHz
    // let frequencies = log10_space(20.0, 20000.0, num_points);
    let frequencies = log_space_skewed(20.0, 20000.0, num_points, 0.3);

    // Vector to hold the magnitude response
    let mut magnitude_response: Vec<f32> = Vec::with_capacity(num_points);

    for &f in &frequencies {
        // Normalized frequency: f_normalized = 2 * pi * f / sample_rate
        let omega = 2.0 * PI * f / sample_rate;

        // Compute the real and imaginary parts of the numerator (b)
        let num_real = b0 + b1 * (-omega).cos() + b2 * (-2.0 * omega).cos();
        let num_imag = b1 * (-omega).sin() + b2 * (-2.0 * omega).sin();

        // Compute the real and imaginary parts of the denominator (a)
        let denom_real = a0 + a1 * (-omega).cos() + a2 * (-2.0 * omega).cos();
        let denom_imag = a1 * (-omega).sin() + a2 * (-2.0 * omega).sin();

        // Compute magnitude of H(f) = |H(f)| = sqrt(real^2 + imag^2)
        let numerator_magnitude = (num_real.powi(2) + num_imag.powi(2)).sqrt();
        let denominator_magnitude = (denom_real.powi(2) + denom_imag.powi(2)).sqrt();

        let h_f = numerator_magnitude / denominator_magnitude;

        // Convert to dB: 20 * log10(|H(f)|)
        let h_db = 20.0 * h_f.log10();

        // Store the result in magnitude_response
        magnitude_response.push(h_db);
    }

    (magnitude_response, frequencies)
}
pub fn calculate_biquad_response2(
    biqad: &crate::BiquadFilter<f32>,
    sample_rate: f32,
    num_points: usize,
) -> (Vec<f32>, Vec<f32>) {
    let (b0, b1, b2) = (biqad.biquad_coeffs.b0, biqad.biquad_coeffs.b1, biqad.biquad_coeffs.b2);
    let (a0, a1, a2) = (1.0, biqad.biquad_coeffs.a1, biqad.biquad_coeffs.a2);

    // Logarithmically spaced frequency points between 20 Hz and 20 kHz
    // let mut frequencies = log_space(20.0, 20000.0, num_points);
    let mut frequencies = log_space_skewed(20.0, 20000.0, num_points, 0.3);
    frequencies.push(biqad.freq);
    frequencies.sort_by(|a, b| a.partial_cmp(b).unwrap());
    // Vector to hold the magnitude response
    let mut magnitude_response: Vec<f32> = Vec::with_capacity(num_points);

    for &f in &frequencies {
        // Normalized frequency: f_normalized = 2 * pi * f / sample_rate
        let omega = 2.0 * PI * f / sample_rate;

        // Compute the real and imaginary parts of the numerator (b)
        let num_real = b0 + b1 * (-omega).cos() + b2 * (-2.0 * omega).cos();
        let num_imag = b1 * (-omega).sin() + b2 * (-2.0 * omega).sin();

        // Compute the real and imaginary parts of the denominator (a)
        let denom_real = a0 + a1 * (-omega).cos() + a2 * (-2.0 * omega).cos();
        let denom_imag = a1 * (-omega).sin() + a2 * (-2.0 * omega).sin();

        // Compute magnitude of H(f) = |H(f)| = sqrt(real^2 + imag^2)
        let numerator_magnitude = (num_real.powi(2) + num_imag.powi(2)).sqrt();
        let denominator_magnitude = (denom_real.powi(2) + denom_imag.powi(2)).sqrt();

        let h_f = numerator_magnitude / denominator_magnitude;

        // Convert to dB: 20 * log10(|H(f)|)
        let h_db = 20.0 * h_f.log10();

        // Store the result in magnitude_response
        magnitude_response.push(h_db);
    }

    (magnitude_response, frequencies)
}


pub fn points_to_smooth_svg_path(points: &[(f32, f32)], smoothness: f32) -> SharedString {
    let mut path = SharedString::new();
    let mut prev = points[0];
    path.push_str(&format!("M {} {}", -2.0, 201.0));
    path.push_str(&format!("L {} {}", prev.0-1.0, prev.1));
    path.push_str(&format!("L {} {}", prev.0, prev.1));
    for i in 1..points.len() {
        let curr = points[i];
        let mid = ((prev.0 + curr.0) / 2.0, (prev.1 + curr.1) / 2.0);
        let cp1 = ((mid.0 + prev.0) / 2.0, (mid.1 + prev.1) / 2.0);
        let cp2 = ((mid.0 + curr.0) / 2.0, (mid.1 + curr.1) / 2.0);
        path.push_str(&format!(" C {} {} {} {} {} {}", cp1.0, cp1.1, cp2.0, cp2.1, curr.0, curr.1));
        prev = curr;
    }
    let last_point = points.last().unwrap();
    path.push_str(&format!("L {} {}", last_point.0+2.0, last_point.1));
    path.push_str(&format!("L {} {}", last_point.0+2.0, 201.0));
    path
}

pub fn create_grid_svg(freq_points: Vec<f32>, gain_points: Vec<f32>) -> SharedString {
    let mut path = SharedString::new();
    for i in 0..freq_points.len() {
        let point = freq_points[i];
        path.push_str(&format!("M {} {}", point, -1000.0));
        path.push_str(&format!("L {} {}", point, 1000.0));
    }   for i in 0..gain_points.len() {
        let point = gain_points[i];
        path.push_str(&format!("M {} {}", -1000,point));
        path.push_str(&format!("L {} {}", 1000, point));
    }

    path
}

mod tests {
#[cfg(test)]
#[test]
fn test() {
    // Example: biquad coefficients for a peak filter
    let b = (0.2929, -0.5858, 0.2929); // Example numerator coefficients
    let a = (1.0, -0.0000, 0.1716); // Example denominator coefficients

    // Sampling rate
    let sample_rate = 44100.0;

    // Number of points for the frequency response
    let num_points = 10;

    // Calculate the frequency response
    let (magnitude_response, frequencies) = super::calculate_biquad_response(b, a, sample_rate, num_points);

    // Print the results
    for (freq, mag) in frequencies.iter().zip(magnitude_response.iter()) {
        dbg!("Frequency: {:.2} Hz, Magnitude: {:.2} dB", freq, mag);
    }
    0;
}

#[cfg(test)]
#[test]
fn test_f_to_lin() {
    let y0 = 20.0;
    let y1 = 20000.0;
    let x0 = 0.0;
    let x1 = 225.0;
    let y = 1000.0;
    let x = super::f_to_lin(y0, y1, x0, x1, y);
    dbg!(x);
}

#[cfg(test)]
#[test]
fn test_lin_to_f() {
    let y0 = 20.0;
    let y1 = 24000.0;
    let x0 = 0.0;
    let x1 = 250.0;
    let x = 120.0;
    let y = super::lin_to_f(y0, y1, x0, x1, x);
    dbg!(y);
}

#[cfg(test)]
#[test]
fn test_log_space_skewed() {
    let start = 20.0;
    let end = 20000.0;
    let num_points = 10;
    let skew_f = 0.3; // Example skew factor
    let frequencies = super::log_space_skewed(start, end, num_points, skew_f);
    for freq in frequencies.iter().enumerate() {
        println!(r" {}:, {}", freq.0,freq.1);
    }
}

#[test]
fn test_lin_to_f_skew() {
    let start = 20.0;
    let end = 20000.0;
    let num_points = 500;
    let skew_f = 0.3; // Example skew factor
    let x = 250.0; // Example x value
    let frequency = super::lin_to_f_skew(start, end, num_points, skew_f, x);
    dbg!(frequency);
}

#[test]
fn test_f_to_lin_skew() {
    let start = 20.0;
    let end = 20000.0;
    let num_points = 500;
    let skew_f = 0.3; // Example skew factor
    let y = 1000.0; // Example y value
    let x = super::f_to_lin_skew(start, end, num_points, skew_f, y);
    dbg!(x); }
}